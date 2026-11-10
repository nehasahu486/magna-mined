"""
Bulk Paper Screening & Multi-Gene Literature Ingestion Engine for Magna-Mined.
Scans newly published M. oryzae papers directly rather than querying gene by gene:
1. Streams candidate M. oryzae experimental papers from Europe PMC.
2. Screens abstracts to skip reviews and identifies all mentioned genes and aliases.
3. Retrieves full text (Europe PMC XML, PDF, or abstract fallback) once per paper.
4. Extracts quote-verified phenotype and knockout evidence without duplicate work.
5. Tracks evaluated papers in SQLite (processed_papers) so re-runs resume cleanly.
"""

from typing import Dict, Any, List, Optional, Iterator
from rich.console import Console
from rich.table import Table

from .database import DatabaseManager
from .literature import EuropePMCClient
from .extractor import EvidenceExtractor
from .schema import GeneEvidenceRecord, EvidenceSource, EvidenceConfidence

console = Console()

DEFAULT_SEARCH_QUERY = (
    '("Magnaporthe oryzae" OR "Pyricularia oryzae" OR "rice blast") '
    'AND (mutant OR knockout OR deletion OR virulence OR appressorium OR conidiation OR pathogenesis)'
)


class PaperMiner:
    """
    Orchestrates automated paper mining, multi-gene extraction, and processed-paper tracking.
    """

    def __init__(
        self,
        db_path: str = "magnaporthe.db",
        api_key: Optional[str] = None,
        model_name: str = "gemini-2.0-flash",
        enable_pdf: bool = True,
        unclassified_only: bool = True,
        check_additional_functions: bool = False
    ):
        self.db = DatabaseManager(db_path)
        self.enable_pdf = enable_pdf
        self.unclassified_only = unclassified_only
        self.check_additional_functions = check_additional_functions
        self.lit_client = EuropePMCClient()
        self.extractor = EvidenceExtractor(api_key=api_key, model_name=model_name)
        self.linker = self.db.linker

    def process_paper(
        self,
        article: Dict[str, Any],
        verbose: bool = True,
        pdf_sections_override: Optional[Dict[str, str]] = None,
        force: bool = False
    ) -> Dict[str, Any]:
        """
        Process a single literature article end-to-end:
        - Check sync ledger to prevent duplicate work
        - Triage abstract for functional experimental genetics
        - Discover all M. oryzae genes mentioned in title, abstract, or full-text
        - Extract grounded phenotype evidence for each detected gene
        - Record provenance in processed_papers and extractions tables

        pdf_sections_override: pre-parsed PDF sections (e.g. a local PDF supplied by a
            curator) used when JATS XML is unavailable.
        force: re-process even if already in the ledger, and skip abstract triage
            (the curator has explicitly chosen this paper). Prior extractions for the
            same paper are replaced so upgraded evidence doesn't duplicate.
        """
        pmid = article.get("pmid")
        pmcid = article.get("pmcid")
        doi = article.get("doi")
        title = article.get("title", "")
        abstract = article.get("abstract", "")
        pub_year = article.get("pub_year")
        is_oa = article.get("is_open_access", False)

        # 1. Check sync ledger (idempotent / resumable)
        if force:
            self.db.delete_extractions_for_paper(pmid=pmid, pmcid=pmcid, doi=doi)
        elif self.db.is_paper_processed(
            pmid=pmid,
            pmcid=pmcid,
            require_additional_check=self.check_additional_functions
        ):
            return {
                "status": "already_processed",
                "pmid": pmid,
                "pmcid": pmcid,
                "genes_extracted": 0,
                "extractions": []
            }

        # 2. Fast abstract triage (skipped when a curator forces a paper in)
        triage = self.extractor.triage_abstract(pmid or "N/A", title, abstract)
        if force:
            triage.is_experimental = True
            triage.rationale = f"Curator-supplied document. {triage.rationale}"
        if not triage.is_experimental:
            self.db.record_processed_paper(
                pmid=pmid,
                pmcid=pmcid,
                doi=doi,
                title=title,
                pub_year=pub_year,
                is_open_access=is_oa,
                is_experimental=False,
                triage_rationale=triage.rationale,
                genes_found=0
            )
            return {
                "status": "non_experimental",
                "pmid": pmid,
                "pmcid": pmcid,
                "genes_extracted": 0,
                "extractions": []
            }

        # 3. Discover paper-local inline aliases (e.g. 'MoAcb1 (MGG_06177)') & title/abstract mentions
        inline_aliases = self.linker.extract_inline_aliases(f"{title}\n{abstract}")
        candidate_mentions = self.linker.find_mentions(f"{title} {abstract}", extra_aliases=inline_aliases)
        candidate_genes: Dict[str, Optional[str]] = {
            m.canonical_id: m.canonical_name for m in candidate_mentions
        }

        # 4. Ingest full text using the source cascade:
        #    JATS XML (Europe PMC) -> layout-aware PDF -> abstract-only fallback
        full_text_body = ""
        all_section_text = ""
        evidence_source = EvidenceSource.ABSTRACT_ONLY

        if pmcid:
            xml_sections = self.lit_client.fetch_fulltext_xml(pmcid)
            if xml_sections:
                full_text_body = f"{xml_sections.get('results', '')}\n{xml_sections.get('methods', '')}"
                all_section_text = "\n".join(xml_sections.values())
                if full_text_body.strip():
                    evidence_source = EvidenceSource.JATS_XML

        if evidence_source == EvidenceSource.ABSTRACT_ONLY and pdf_sections_override is not None:
            pdf_sections = pdf_sections_override
        elif evidence_source == EvidenceSource.ABSTRACT_ONLY and self.enable_pdf and (article.get("pdf_url") or doi):
            pdf_sections = self.lit_client.fetch_fulltext_pdf(
                pmcid=pmcid, pdf_url=article.get("pdf_url"), doi=doi
            )
        else:
            pdf_sections = None

        if pdf_sections:
            body = f"{pdf_sections.get('results', '')}\n{pdf_sections.get('methods', '')}"
            # Heading detection can fail on PDFs; fall back to the whole document
            full_text_body = body if body.strip() else pdf_sections.get("other", "")
            all_section_text = "\n".join(pdf_sections.values())
            if full_text_body.strip():
                evidence_source = EvidenceSource.PDF

        # Learn inline aliases from full text (including Introduction / 'other')
        if all_section_text.strip():
            for mgg_id, syms in self.linker.extract_inline_aliases(all_section_text).items():
                bucket = inline_aliases.setdefault(mgg_id, [])
                for s in syms:
                    if s not in bucket:
                        bucket.append(s)

        if full_text_body.strip():
            # Scan full text for additional genes not mentioned in abstract
            for bm in self.linker.find_mentions(full_text_body, extra_aliases=inline_aliases):
                if bm.canonical_id not in candidate_genes:
                    candidate_genes[bm.canonical_id] = bm.canonical_name
        else:
            full_text_body = ""
            evidence_source = EvidenceSource.ABSTRACT_ONLY

        # Filter candidate genes:
        # - Main mode (unclassified_only=True, check_additional_functions=False): only unclassified genes
        # - Double-check mode (check_additional_functions=True): inspect already-classified genes for NEW additional roles
        skipped_baseline = 0
        unclassified_ids = self.db.get_unclassified_gene_ids()
        if self.unclassified_only and not self.check_additional_functions and candidate_genes:
            filtered_genes = {
                gid: gname for gid, gname in candidate_genes.items()
                if gid.upper() in unclassified_ids
            }
            skipped_baseline = len(candidate_genes) - len(filtered_genes)
            candidate_genes = filtered_genes

        # If no target genes remain, record paper in ledger with 0 mapped genes
        if not candidate_genes:
            rationale_suffix = (
                f" ({skipped_baseline} baseline-classified genes skipped)"
                if skipped_baseline else " (No MGG genes mapped)"
            )
            self.db.record_processed_paper(
                pmid=pmid,
                pmcid=pmcid,
                doi=doi,
                title=title,
                pub_year=pub_year,
                is_open_access=is_oa,
                is_experimental=True,
                triage_rationale=f"{triage.rationale}{rationale_suffix}",
                genes_found=0,
                evidence_source=evidence_source.value,
                checked_additional=self.check_additional_functions
            )
            return {
                "status": "no_genes_mapped",
                "pmid": pmid,
                "pmcid": pmcid,
                "genes_extracted": 0,
                "skipped_baseline_genes": skipped_baseline,
                "additional_functions_found": 0,
                "extractions": []
            }

        # 5. Extract grounded evidence for each discovered gene
        saved_extractions: List[GeneEvidenceRecord] = []
        additional_hits = 0
        for gene_id, gene_name in candidate_genes.items():
            is_unclassified = gene_id.upper() in unclassified_ids
            if self.check_additional_functions and self.unclassified_only and is_unclassified:
                # Double-check mode can run alongside unclassified curation seamlessly
                pass

            local_aliases = inline_aliases.get(gene_id, [])
            aliases = list(dict.fromkeys(self.linker.get_aliases_for_gene(gene_id) + local_aliases))
            effective_name = gene_name or (local_aliases[0] if local_aliases else None)
            snippets = []
            record_source = evidence_source
            if full_text_body:
                snippets = self.lit_client.extract_targeted_snippets(
                    full_text_body,
                    gene_id=gene_id,
                    aliases=aliases
                )

            if not snippets:
                # No full-text mention of this gene: ground in the abstract only
                snippets = [abstract] if abstract else [title]
                record_source = EvidenceSource.ABSTRACT_ONLY

            record = self.extractor.extract_evidence_from_snippets(
                gene_id=gene_id,
                gene_name=effective_name,
                title=title,
                abstract=abstract,
                snippets=snippets,
                pmid=pmid,
                pmcid=pmcid,
                doi=doi,
                pub_year=pub_year
            )

            if record:
                record.evidence_source = record_source
                if record_source == EvidenceSource.ABSTRACT_ONLY:
                    # Abstract-only evidence can never exceed the lowest-confidence tier
                    # so curators know full text was unavailable / unverified.
                    record.confidence_tier = EvidenceConfidence.TIER_4

                # For baseline-classified genes under double-check mode, only save if a NEW
                # additional function (or knockout report) is discovered!
                if not is_unclassified and self.check_additional_functions:
                    gene_row = self.db.get_gene(gene_id) or {}
                    base_set = self.db.normalize_process_set(gene_row.get("baseline_process"))
                    rec_procs = [record.process_affected.value] + [p.value for p in record.secondary_processes]
                    new_procs = [
                        p for p in rec_procs
                        if p.lower() not in ("not reported", "dispensable / wild-type-like")
                        and p.lower() not in base_set
                    ]
                    needs_ko = (
                        not gene_row.get("knockout_reported")
                        and record.experimental_method.value in self.db._KNOCKOUT_METHODS
                    )
                    if not new_procs and not needs_ko:
                        continue
                    if new_procs:
                        additional_hits += 1

                try:
                    self.db.save_evidence(record, preserve_baseline_primary=(not is_unclassified))
                    saved_extractions.append(record)
                except Exception as e:
                    if verbose:
                        console.print(f"[yellow]Warning: Could not save evidence for {gene_id}: {e}[/yellow]")

        # 6. Record completion in ledger
        self.db.record_processed_paper(
            pmid=pmid,
            pmcid=pmcid,
            doi=doi,
            title=title,
            pub_year=pub_year,
            is_open_access=is_oa,
            is_experimental=True,
            triage_rationale=triage.rationale,
            genes_found=len(saved_extractions),
            evidence_source=evidence_source.value,
            checked_additional=self.check_additional_functions
        )

        return {
            "status": "curated",
            "pmid": pmid,
            "pmcid": pmcid,
            "genes_extracted": len(saved_extractions),
            "skipped_baseline_genes": skipped_baseline,
            "additional_functions_found": additional_hits,
            "extractions": saved_extractions
        }

    def mine_papers(
        self,
        max_papers: int = 50,
        page_size: int = 25,
        open_access_only: bool = False,
        query: Optional[str] = None,
        year_from: Optional[int] = None,
        year_to: Optional[int] = None,
        days: Optional[int] = None,
        since: Optional[str] = None,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Stream and process M. oryzae papers in batches.
        - `year_from` / `year_to`: iterates year-by-year (up to `max_papers`/year).
        - `days` (e.g. 90 for quarterly, 180 for semi-annual) or `since` (YYYY-MM-DD):
          restricts to papers indexed/published in that recent window.
        """
        from datetime import date, timedelta

        search_query = query or DEFAULT_SEARCH_QUERY
        if year_from is not None and year_to is not None and year_from > year_to:
            raise ValueError(f"year_from ({year_from}) cannot exceed year_to ({year_to})")

        date_from: Optional[str] = since
        date_to: Optional[str] = None
        if days is not None:
            today = date.today()
            date_from = (today - timedelta(days=days)).isoformat()
            date_to = today.isoformat()

        year_slices: List[Optional[int]]
        if not date_from and year_from is not None and year_to is not None:
            year_slices = list(range(year_from, year_to + 1))
        else:
            year_slices = [None]

        if verbose:
            console.print(f"\n[bold blue]=== Magna-Mined Paper Mining Starting ===[/bold blue]")
            if self.check_additional_functions and self.unclassified_only:
                mode_label = "Unclassified MagnaGenes + Double-Check Additional Functions"
            elif self.check_additional_functions:
                mode_label = "Double-Check Additional Functions"
            elif self.unclassified_only:
                mode_label = "Unclassified MagnaGenes Only (default)"
            else:
                mode_label = "All Genes (including baseline-classified)"
            console.print(f"Curation Scope: [cyan]{mode_label}[/cyan]")
            if date_from:
                console.print(f"Date Window: [cyan]{date_from} -> {date_to or 'today'}[/cyan] (max [cyan]{max_papers}[/cyan] papers)")
            elif year_from is not None and year_to is not None:
                console.print(f"Publication Years: [cyan]{year_from} -> {year_to}[/cyan] (up to [cyan]{max_papers}[/cyan] papers/year)")
            elif year_from is not None or year_to is not None:
                console.print(f"Publication Year Filter: [cyan]{year_from or '..'} to {year_to or '..'}[/cyan] (max [cyan]{max_papers}[/cyan] papers)")
            else:
                console.print(f"Target Max Papers: [cyan]{max_papers}[/cyan]")
            console.print(f"Open Access Only: [cyan]{open_access_only}[/cyan]")

        stats = {
            "streamed": 0,
            "already_processed": 0,
            "non_experimental": 0,
            "curated": 0,
            "extractions_created": 0,
            "additional_functions_discovered": 0,
            "skipped_baseline_genes": 0,
            "genes_affected": set()
        }

        for yr in year_slices:
            if verbose and yr is not None:
                console.print(f"\n[bold magenta]--- Year {yr} ---[/bold magenta]")

            stream_kwargs: Dict[str, Any] = {
                "query": search_query,
                "max_papers": max_papers,
                "page_size": page_size,
                "open_access_only": open_access_only,
            }
            if date_from:
                stream_kwargs["date_from"] = date_from
                stream_kwargs["date_to"] = date_to
            elif yr is not None:
                stream_kwargs["pub_year"] = yr
            elif year_from is not None or year_to is not None:
                stream_kwargs["year_from"] = year_from
                stream_kwargs["year_to"] = year_to

            stream = self.lit_client.stream_papers(**stream_kwargs)

            for article in stream:
                stats["streamed"] += 1
                pmid = article.get("pmid")
                art_year = article.get("pub_year") or yr or "?"
                title = article.get("title", "")[:70]

                result = self.process_paper(article, verbose=verbose)
                status = result["status"]

                if status == "already_processed":
                    stats["already_processed"] += 1
                    if verbose:
                        console.print(f"[dim]• [{art_year} | PMID:{pmid}] Already processed - skipping[/dim]")
                elif status == "non_experimental":
                    stats["non_experimental"] += 1
                    if verbose:
                        console.print(f"[yellow]• [{art_year} | PMID:{pmid}] Filtered out (non-experimental):[/yellow] {title}...")
                elif status in ("curated", "no_genes_mapped"):
                    stats["curated"] += 1
                    num_genes = result["genes_extracted"]
                    skipped_b = result.get("skipped_baseline_genes", 0)
                    add_f = result.get("additional_functions_found", 0)
                    stats["skipped_baseline_genes"] += skipped_b
                    stats["additional_functions_discovered"] += add_f
                    stats["extractions_created"] += num_genes
                    for rec in result.get("extractions", []):
                        stats["genes_affected"].add(rec.gene_id)
                    if verbose:
                        skip_note = f" (skipped {skipped_b} baseline-classified)" if skipped_b else ""
                        console.print(f"[bold green]✓ [{art_year} | PMID:{pmid}] Curated {num_genes} genes{skip_note}:[/bold green] {title}...")
                        for rec in result.get("extractions", []):
                            console.print(f"    - [cyan]{rec.gene_id}[/cyan] ({rec.gene_name or ''}): {rec.process_affected.value} [{rec.confidence_tier.value}, {rec.evidence_source.value}]")

        stats["unique_genes_curated"] = len(stats["genes_affected"])
        stats.pop("genes_affected", None)

        if verbose:
            console.print(f"\n[bold green]=== Paper Mining Complete ===[/bold green]")
            table = Table(title="Paper Mining Summary")
            table.add_column("Metric", style="cyan")
            table.add_column("Count", style="green")
            for k, v in stats.items():
                table.add_row(k.replace("_", " ").title(), str(v))
            console.print(table)

        return stats

    def ingest_local_pdf(
        self,
        pdf_path: str,
        pmid: Optional[str] = None,
        doi: Optional[str] = None,
        title: Optional[str] = None,
        pub_year: Optional[int] = None,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Ingest a PDF already on disk (e.g. a paywalled paper the lab has access to).
        If a PMID is given, title/abstract/year are pulled from Europe PMC metadata.
        """
        from pathlib import Path
        path = Path(pdf_path)
        if not path.exists():
            raise FileNotFoundError(pdf_path)

        sections = self.lit_client.pdf_to_sections(path.read_bytes())
        if not sections:
            raise ValueError(f"Could not convert PDF to text: {pdf_path} (is pymupdf4llm installed?)")

        article: Dict[str, Any] = {
            "pmid": pmid, "pmcid": None, "doi": doi,
            "title": title or "", "abstract": sections.get("abstract", ""),
            "pub_year": pub_year, "is_open_access": False,
        }
        if pmid:
            meta = self.lit_client.search_literature(f"EXT_ID:{pmid} AND SRC:MED", max_results=1)
            if meta:
                m = meta[0]
                article.update({
                    "pmcid": m.get("pmcid"),
                    "doi": doi or m.get("doi"),
                    "title": title or m.get("title", ""),
                    "abstract": m.get("abstract") or article["abstract"],
                    "pub_year": pub_year or m.get("pub_year"),
                    "is_open_access": m.get("is_open_access", False),
                })
        if not article["pmid"]:
            # Ledger key: fall back to a stable file-derived id when no PMID is known
            article["pmid"] = f"local:{doi or path.name}"
        if not article["title"]:
            article["title"] = path.stem

        result = self.process_paper(article, verbose=verbose, pdf_sections_override=sections, force=True)
        if verbose:
            console.print(f"[bold green]✓ Ingested {path.name}:[/bold green] {result['genes_extracted']} gene records")
            for rec in result.get("extractions", []):
                console.print(
                    f"    - [cyan]{rec.gene_id}[/cyan] ({rec.gene_name or ''}): "
                    f"{rec.process_affected.value} [{rec.confidence_tier.value}, {rec.evidence_source.value}]"
                )
        return result
