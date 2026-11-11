"""
End-to-end Curation Pipeline & CLI for Magna-Mined.
Demonstrates literature search, triage, JATS XML extraction, and grounded evidence logging.
"""

import sys
from pathlib import Path
from typing import Optional, List, Dict, Any
import click
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

from .database import DatabaseManager
from .literature import EuropePMCClient
from .extractor import EvidenceExtractor
from .schema import GeneEvidenceRecord, PhenotypeCategory

console = Console()


class CurationPipeline:
    """Orchestrates literature discovery, XML parsing, and grounded extraction."""

    def __init__(self, db_path: str = "magnaporthe.db", api_key: Optional[str] = None):
        self.db = DatabaseManager(db_path)
        self.lit_client = EuropePMCClient()
        self.extractor = EvidenceExtractor(api_key=api_key)

    def curate_single_gene(
        self,
        gene_id_or_name: str,
        max_papers: int = 5,
        verbose: bool = True
    ) -> List[GeneEvidenceRecord]:
        """
        Execute full curation lifecycle for a target gene:
        1. Resolve gene identifier from catalog
        2. Query consolidated literature (PubMed / Europe PMC)
        3. Triage abstracts for experimental relevance
        4. Ingest clean JATS XML (or abstract fallback)
        5. Extract grounded evidence and persist to SQLite
        """
        gene_record = self.db.get_gene(gene_id_or_name)
        if not gene_record:
            if verbose:
                console.print(f"[bold red]Gene '{gene_id_or_name}' not found in database catalog.[/bold red]")
            return []

        gene_id = gene_record["gene"]
        gene_name = gene_record["gene_name"]

        if verbose:
            console.print(f"\n[bold cyan]Curating Gene:[/bold cyan] {gene_id} ({gene_name or 'unnamed'})")

        # 1. Literature Search (expanded with all known aliases)
        aliases = self.db.linker.get_aliases_for_gene(gene_id)
        articles = self.lit_client.search_gene(
            gene_id,
            gene_name=gene_name,
            aliases=aliases,
            max_results=max_papers
        )
        if not articles:
            if verbose:
                console.print(f"[yellow]No literature found for {gene_id}.[/yellow]")
            return []

        if verbose:
            console.print(f"Found [green]{len(articles)}[/green] candidate articles.")

        extractions: List[GeneEvidenceRecord] = []

        # 2. Process articles
        for art in articles:
            pmid = art.get("pmid")
            pmcid = art.get("pmcid")
            title = art.get("title", "")
            abstract = art.get("abstract", "")

            # Triage
            triage = self.extractor.triage_abstract(pmid or "N/A", title, abstract)
            if not triage.is_experimental:
                continue

            # Ingest full text if open access XML is available
            snippets = []
            if pmcid:
                xml_sections = self.lit_client.fetch_fulltext_xml(pmcid)
                if xml_sections:
                    combined_body = f"{xml_sections.get('results', '')}\n{xml_sections.get('methods', '')}"
                    snippets = self.lit_client.extract_targeted_snippets(
                        combined_body,
                        gene_id=gene_id,
                        aliases=aliases
                    )

            if not snippets:
                snippets = [abstract]

            evidence = self.extractor.extract_evidence_from_snippets(
                gene_id=gene_id,
                gene_name=gene_name,
                title=title,
                abstract=abstract,
                snippets=snippets,
                pmid=pmid,
                pmcid=pmcid,
                doi=art.get("doi"),
                pub_year=art.get("pub_year")
            )

            if evidence:
                try:
                    self.db.save_evidence(evidence)
                    extractions.append(evidence)
                    if verbose:
                        console.print(f"  [green]✓ Extracted:[/green] {evidence.process_affected.value} ({evidence.confidence_tier.value}) from PMID:{pmid}")
                except Exception as e:
                    if verbose:
                        console.print(f"  [red]Failed to save evidence:[/red] {e}")

        return extractions


@click.group()
def cli():
    """Magna-Mined: Automated Magnaporthe Literature Miner."""
    pass


@cli.command()
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--excel", default="data/Magnagenes.xlsx", help="Path to Magnagenes.xlsx")
@click.option("--csv", default="data/MagnaGenes.csv", help="Path to MagnaGenes.csv fallback")
def init_db(db, excel, csv):
    """Initialize database and seed baseline genes from Excel or CSV."""
    manager = DatabaseManager(db)
    excel_file = Path(excel)
    csv_file = Path(csv)

    if excel_file.exists():
        console.print(f"Seeding genes from Excel: [cyan]{excel}[/cyan]...")
        count = manager.seed_from_excel(str(excel_file))
        console.print(f"[bold green]Successfully seeded {count} genes into {db} from {excel}![/bold green]")
    elif csv_file.exists():
        console.print(f"Seeding genes from CSV: [cyan]{csv}[/cyan]...")
        count = manager.seed_from_csv(str(csv_file))
        console.print(f"[bold green]Successfully seeded {count} genes into {db} from {csv}![/bold green]")
    else:
        console.print(f"[yellow]Schema initialized at {db}. Neither {excel} nor {csv} found.[/yellow]")


@cli.command(name="sync-aliases")
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
def sync_aliases(db):
    """Fetch M. oryzae gene names & synonyms from UniProt and index them in `gene_aliases`."""
    from .entity_linker import sync_uniprot_aliases

    manager = DatabaseManager(db)
    console.print("Fetching M. oryzae gene names & synonyms from UniProt...")
    res = sync_uniprot_aliases(manager)
    console.print(
        f"[bold green]✓ Alias sync complete:[/bold green] "
        f"scanned {res['uniprot_rows_scanned']} UniProt entries, "
        f"added {res['aliases_inserted']} new gene aliases."
    )


@cli.command()
@click.argument("gene")
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--max-papers", default=5, help="Max papers to inspect")
def curate(gene, db, max_papers):
    """Curate a specific gene (e.g. MGG_01001 or PMK1)."""
    pipeline = CurationPipeline(db_path=db)
    pipeline.curate_single_gene(gene, max_papers=max_papers, verbose=True)


@cli.command()
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
def benchmark(db):
    """Run extraction benchmark on 4 canonical genes (ENG1, PMK1, MAC1, control)."""
    pipeline = CurationPipeline(db_path=db)
    test_genes = [
        ("MGG_01001", "ENG1"),
        ("MGG_09565", "PMK1"),
        ("MGG_09898", "MAC1"),
        ("MGG_16142", "Uncharacterized")
    ]

    console.print(Panel.fit("[bold blue]Running Magna-Mined Canonical Benchmark[/bold blue]"))

    table = Table(title="Benchmark Results")
    table.add_column("Gene ID", style="cyan")
    table.add_column("Gene Name", style="magenta")
    table.add_column("Extractions Found", justify="right")
    table.add_column("Top Process", style="green")

    for gene_id, name in test_genes:
        extractions = pipeline.curate_single_gene(gene_id, max_papers=3, verbose=False)
        top_process = extractions[0].process_affected.value if extractions else "No literature"
        table.add_row(gene_id, name, str(len(extractions)), top_process)

    console.print(table)


@cli.command(name="mine-papers")
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--max-papers", default=25, help="Max papers to fetch (per year when --year-from/--year-to are set; 0 = unlimited)")
@click.option("--page-size", default=25, help="Batch size per Europe PMC query")
@click.option("--open-access-only/--all", default=False, help="Restrict to open access papers with full text XML")
@click.option("--query", default=None, help="Custom Europe PMC search query")
@click.option("--year-from", type=int, default=None, help="Start publication year (e.g. 2022)")
@click.option("--year-to", type=int, default=None, help="End publication year (e.g. 2023)")
@click.option("--days", type=int, default=None, help="Only fetch papers from the last N days (e.g. 90 for quarterly, 180 for 6-month sync)")
@click.option("--since", default=None, help="Only fetch papers published on/after YYYY-MM-DD")
@click.option("--pdf/--no-pdf", default=True, help="Fall back to PDF text extraction when full-text XML is unavailable")
@click.option("--unclassified-only/--include-classified", default=True, help="Curate only genes unclassified ('Not reported') in baseline MagnaGenes (default)")
@click.option("--double-check", is_flag=True, help="Also inspect already-classified baseline genes for newly discovered additional functions")
def mine_papers(db, max_papers, page_size, open_access_only, query, year_from, year_to, days, since, pdf, unclassified_only, double_check):
    """Scan newly published M. oryzae papers and curate unclassified MagnaGenes by default."""
    from .paper_miner import PaperMiner
    engine = PaperMiner(
        db_path=db,
        enable_pdf=pdf,
        unclassified_only=unclassified_only,
        check_additional_functions=double_check,
    )
    engine.mine_papers(
        max_papers=max_papers,
        page_size=page_size,
        open_access_only=open_access_only,
        query=query,
        year_from=year_from,
        year_to=year_to,
        days=days,
        since=since,
        verbose=True
    )


@cli.command(name="double-check")
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--max-papers", default=25, help="Max papers to inspect (per year when --year-from/--year-to are set)")
@click.option("--page-size", default=25, help="Batch size per Europe PMC query")
@click.option("--open-access-only/--all", default=False, help="Restrict to open access papers with full text XML")
@click.option("--query", default=None, help="Custom Europe PMC search query")
@click.option("--year-from", type=int, default=None, help="Start publication year (e.g. 2022)")
@click.option("--year-to", type=int, default=None, help="End publication year (e.g. 2023)")
@click.option("--days", type=int, default=None, help="Only check papers from the last N days (e.g. 90 or 180)")
@click.option("--since", default=None, help="Only check papers published on/after YYYY-MM-DD")
@click.option("--pdf/--no-pdf", default=True, help="Fall back to PDF text extraction when full-text XML is unavailable")
def double_check(db, max_papers, page_size, open_access_only, query, year_from, year_to, days, since, pdf):
    """Double-check classified baseline genes for newly reported additional functions & knockouts."""
    from .paper_miner import PaperMiner
    engine = PaperMiner(
        db_path=db,
        enable_pdf=pdf,
        unclassified_only=True,
        check_additional_functions=True,
    )
    engine.mine_papers(
        max_papers=max_papers,
        page_size=page_size,
        open_access_only=open_access_only,
        query=query,
        year_from=year_from,
        year_to=year_to,
        days=days,
        since=since,
        verbose=True
    )


@cli.command()
@click.argument("pdf_path", type=click.Path(exists=True, dir_okay=False))
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--pmid", default=None, help="PubMed ID (metadata is fetched from Europe PMC)")
@click.option("--doi", default=None, help="DOI of the paper")
@click.option("--title", default=None, help="Paper title (if no PMID)")
@click.option("--year", type=int, default=None, help="Publication year (if no PMID)")
@click.option("--unclassified-only/--include-classified", default=True, help="Curate only genes unclassified in baseline MagnaGenes (default)")
def ingest_pdf(pdf_path, db, pmid, doi, title, year, unclassified_only):
    """Ingest a local PDF (such as a paywalled paper stored locally)."""
    from .paper_miner import PaperMiner
    engine = PaperMiner(db_path=db, unclassified_only=unclassified_only)
    engine.ingest_local_pdf(pdf_path, pmid=pmid, doi=doi, title=title, pub_year=year, verbose=True)


@cli.command()
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--review-id", type=int, default=None, help="Mark a specific extraction ID as human-verified")
@click.option("--notes", default="", help="Curator review notes (used with --review-id)")
def reconcile(db, review_id, notes):
    """Reconcile competing evidence extractions and inspect the human review queue."""
    from .reconciler import EvidenceReconciler

    rec = EvidenceReconciler(db_path=db)
    if review_id is not None:
        res = rec.mark_reviewed(review_id, notes=notes)
        console.print(
            f"[bold green]✓ Marked extraction #{review_id} human-verified:[/bold green] "
            f"{res['gene_id']} -> {res['consensus_process']} ({res['conflict_status']})"
        )
        return

    summary = rec.reconcile_all()
    table = Table(title="Evidence Reconciliation Summary")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="green")
    table.add_row("Genes Reconciled", str(summary["genes_reconciled"]))
    table.add_row("Conflicts Detected", str(summary["conflicts_detected"]))
    table.add_row("Resolved By Human", str(summary["resolved_by_human"]))
    for tier, count in sorted(summary["tier_breakdown"].items()):
        table.add_row(f"Consensus Tier ({tier})", str(count))
    console.print(table)

    queue = rec.get_review_queue()
    if queue:
        console.print(f"\n[bold yellow]Human Review Queue ({len(queue)} conflicts):[/bold yellow]")
        for item in queue[:10]:
            console.print(f"  [cyan]{item['gene_id']}[/cyan] ({item['gene_name'] or 'unnamed'}) -> current: {item['current_consensus']}")
            for ext in item["competing_extractions"]:
                console.print(
                    f"    • [id={ext['id']}] PMID:{ext['pmid']} | {ext['process_affected']} "
                    f"[{ext['confidence_tier']}, {ext['evidence_source']}]"
                )


@cli.command(name="eval-benchmark")
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--mode", type=click.Choice(["db", "live"]), default="db", help="'db' evaluates stored extractions; 'live' runs blind held-out Zenodo sample")
@click.option("--sample-size", default=8, help="Number of canonical Zenodo genes for live benchmark")
def eval_benchmark(db, mode, sample_size):
    """Evaluate extraction accuracy against Zenodo 4647766 ground truth."""
    from .benchmark import ZenodoBenchmark

    bench = ZenodoBenchmark(db_path=db)
    if mode == "live":
        console.print(Panel.fit(f"[bold blue]Running Blind Live Zenodo Benchmark (sample={sample_size})[/bold blue]"))
        metrics = bench.run_live_zenodo_benchmark(sample_size=sample_size)
    else:
        console.print(Panel.fit("[bold blue]Evaluating Current Database vs Zenodo Ground Truth[/bold blue]"))
        metrics = bench.evaluate_current_db()

    summary = Table(title="Zenodo Ground-Truth Benchmark Metrics")
    summary.add_column("Metric", style="cyan")
    summary.add_column("Score / Count", style="green")
    for key in (
        "evaluated_genes", "genes_with_extractions", "target_genes",
        "genes_with_predictions", "retrieval_coverage", "quote_grounding_rate",
        "micro_precision", "micro_recall", "micro_f1",
        "any_overlap_rate", "exact_match_rate", "novel_discoveries_count",
    ):
        if key in metrics:
            summary.add_row(key.replace("_", " ").title(), str(metrics[key]))
    console.print(summary)

    if metrics.get("per_category"):
        cat_table = Table(title="Per-Category Breakdown")
        cat_table.add_column("Category", style="cyan")
        cat_table.add_column("Precision", justify="right")
        cat_table.add_column("Recall", justify="right")
        cat_table.add_column("F1", justify="right", style="green")
        cat_table.add_column("TP / FP / FN", justify="right")
        for cat, vals in metrics["per_category"].items():
            cat_table.add_row(
                cat,
                f"{vals['precision']:.2f}",
                f"{vals['recall']:.2f}",
                f"{vals['f1']:.2f}",
                f"{vals['tp']} / {vals['fp']} / {vals['fn']}",
            )
        console.print(cat_table)


@cli.command()
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--query", "-q", default=None, help="Free-text filter (MGG ID, gene name, Pfam, Blast2GO)")
@click.option("--process", default=None, help="Filter by affected biological process (e.g. virulence, conidiation)")
@click.option("--confidence", default=None, help="Filter by confidence tier (e.g. tier_1_high, tier_2_high)")
@click.option("--source", default=None, help="Filter by evidence source (jats_xml, pdf, abstract_only)")
@click.option("--has-extractions", is_flag=True, help="Only show genes with automated extraction records")
@click.option("--newly-classified", is_flag=True, help="Only show genes unclassified in baseline MagnaGenes that now have extractions")
@click.option("--knockout-only", is_flag=True, help="Only show genes with a reported knockout/CRISPR mutant")
@click.option("--with-additional", is_flag=True, help="Only show baseline-classified genes with newly discovered additional functions")
@click.option("--conflicts", is_flag=True, help="Only show genes flagged with evidence conflicts")
@click.option("--limit", default=25, help="Maximum rows to display")
def search(db, query, process, confidence, source, has_extractions, newly_classified, knockout_only, with_additional, conflicts, limit):
    """Search and filter the curated gene catalog."""
    from .explorer import CatalogExplorer

    explorer = CatalogExplorer(db_path=db)
    rows = explorer.search_genes(
        query=query,
        process=process,
        confidence=confidence,
        source=source,
        has_extractions=has_extractions,
        newly_classified=newly_classified,
        knockout_only=knockout_only,
        has_additional_functions=with_additional,
        conflicts_only=conflicts,
        limit=limit,
    )

    table = Table(title=f"Magna-Mined Gene Search ({len(rows)} results)")
    table.add_column("Gene ID", style="cyan")
    table.add_column("Name", style="magenta")
    table.add_column("Baseline", style="dim")
    table.add_column("Curated Process", style="green")
    table.add_column("Additional Funcs", style="yellow")
    table.add_column("KO?")
    table.add_column("Tier")
    table.add_column("Conflict")

    for r in rows:
        table.add_row(
            r["gene"],
            r.get("gene_name") or "-",
            r.get("baseline_process") or "Not reported",
            r.get("process_affected") or "not reported",
            r.get("additional_functions") or "-",
            "Yes" if r.get("knockout_reported") else "No",
            r.get("consensus_tier") or "-",
            r.get("conflict_status") or "none",
        )
    console.print(table)


@cli.command()
@click.argument("gene")
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
def show(gene, db):
    """Inspect a gene's full catalog profile and grounded extraction quotes."""
    from .explorer import CatalogExplorer

    explorer = CatalogExplorer(db_path=db)
    dossier = explorer.get_gene_dossier(gene)
    if not dossier:
        console.print(f"[bold red]Gene '{gene}' not found in catalog.[/bold red]")
        return

    g = dossier["gene"]
    ko_str = f"Yes: {g.get('knockout_phenotype')}" if g.get("knockout_reported") else "No"
    add_str = g.get("additional_functions") or "-"
    add_ref = g.get("additional_function_references") or "-"
    header = (
        f"[bold cyan]{g['gene']}[/bold cyan] ({g.get('gene_name') or 'unnamed'})\n"
        f"Baseline: {g.get('baseline_process') or 'Not reported'} | "
        f"Curated Process: [green]{g.get('process_affected')}[/green]\n"
        f"Knockout Reported: {ko_str}\n"
        f"Additional Functions: [yellow]{add_str}[/yellow] (Refs: {add_ref})\n"
        f"Consensus Tier: {g.get('consensus_tier') or 'baseline'} | "
        f"Conflict: {g.get('conflict_status') or 'none'}\n"
        f"Aliases: {', '.join(dossier['aliases'])}\n"
        f"Blast2GO: {g.get('blast2go_annotation') or '-'}\n"
        f"Pfam: {g.get('pfam') or '-'}"
    )
    console.print(Panel(header, title="Gene Dossier"))

    exts = dossier["extractions"]
    if not exts:
        console.print("[dim]No automated literature extractions logged yet for this gene.[/dim]")
        return

    for idx, e in enumerate(exts, 1):
        rev = " [HUMAN VERIFIED]" if e.get("reviewed_by_human") else ""
        console.print(
            f"\n[bold]#{idx}[/bold] [cyan]PMID:{e.get('pmid')}[/cyan] ({e.get('publication_year') or '?'}) "
            f"| {e['process_affected']} | [{e['confidence_tier']}, {e.get('evidence_source', 'jats_xml')}]{rev}"
        )
        console.print(f"   Method: {e['experimental_method']} (Complementation: {bool(e['complementation_tested'])})")
        console.print(f"   Quote: [italic]\"{e['verbatim_quote']}\"[/italic]")


@cli.command()
@click.argument("output_path")
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--format", "fmt", type=click.Choice(["csv", "tsv", "xlsx", "json", "html"]), default=None, help="Output format (inferred from extension if omitted)")
@click.option("--query", "-q", default=None, help="Optional text filter")
@click.option("--process", default=None, help="Filter by affected biological process")
@click.option("--confidence", default=None, help="Filter by confidence tier")
@click.option("--source", default=None, help="Filter by evidence source")
@click.option("--has-extractions", is_flag=True, help="Export only genes with automated extractions")
@click.option("--newly-classified", is_flag=True, help="Export only newly classified genes (previously 'Not reported' in MagnaGenes)")
@click.option("--knockout-only", is_flag=True, help="Export only genes with a reported knockout/CRISPR mutant")
@click.option("--with-additional", is_flag=True, help="Export only baseline-classified genes with newly discovered additional functions")
@click.option("--conflicts", is_flag=True, help="Export only genes flagged with conflicts")
def export(output_path, db, fmt, query, process, confidence, source, has_extractions, newly_classified, knockout_only, with_additional, conflicts):
    """Export curated gene catalog and extractions to CSV, TSV, Excel, JSON, or interactive HTML."""
    from .explorer import CatalogExplorer

    explorer = CatalogExplorer(db_path=db)
    count = explorer.export_catalog(
        output_path=output_path,
        fmt=fmt,
        query=query,
        process=process,
        confidence=confidence,
        source=source,
        has_extractions=has_extractions,
        newly_classified=newly_classified,
        knockout_only=knockout_only,
        has_additional_functions=with_additional,
        conflicts_only=conflicts,
    )
    console.print(f"[bold green]✓ Exported {count} genes to {output_path}[/bold green]")


@cli.command(name="export-magnagenes")
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
@click.option("--output-dir", default="exports/magnagenes", help="Directory where fin_data_YYYY-MM-DD.tsv is written (e.g. magnagenes_v2/data-raw/source/annotation)")
@click.option("--base-fin-data", default=None, type=click.Path(exists=True, dir_okay=False), help="Optional existing fin_data.tsv to merge updates into")
@click.option("--merge-existing/--no-merge-existing", default=True, help="Automatically merge into the latest fin_data*.tsv in --output-dir if present")
@click.option("--scope", type=click.Choice(["curated", "full"]), default="curated", help="'curated' exports updated/new genes; 'full' exports entire catalog")
@click.option("--date", "export_date", default=None, help="Override date tag (YYYY-MM-DD, defaults to today)")
@click.option("--archive/--no-archive", default=True, help="Move older fin_data*.tsv files into archive/ with SHA-256 deduplication")
def export_magnagenes(db, output_dir, base_fin_data, merge_existing, scope, export_date, archive):
    """Export dated fin_data_YYYY-MM-DD.tsv for MagnaGenes v2.0 and deduplicate archive/."""
    from .magnagenes_sync import MagnaGenesSyncExporter

    exporter = MagnaGenesSyncExporter(db_path=db)
    res = exporter.export_dated_fin_data(
        output_dir=output_dir,
        export_date=export_date,
        base_fin_data=base_fin_data,
        merge_existing=merge_existing,
        scope=scope,
        archive_existing=archive,
    )

    table = Table(title="MagnaGenes v2.0 Export Summary")
    table.add_column("Field", style="cyan")
    table.add_column("Value", style="green")
    table.add_row("Status", res["status"])
    table.add_row("Active File", res["active_file"])
    table.add_row("Genes Included", str(res["genes_included"]))
    table.add_row("Long-Format Rows", str(res["rows_written"]))
    table.add_row("SHA-256", res["sha256"][:16] + "...")
    if res.get("archived_files"):
        table.add_row("Archived Older Files", ", ".join(res["archived_files"]))
    if res.get("duplicates_removed"):
        table.add_row("Duplicates Skipped", ", ".join(res["duplicates_removed"]))
    console.print(table)


@cli.command()
@click.option("--db", default="magnaporthe.db", help="Path to SQLite database")
def stats(db):
    """Display database statistics."""
    manager = DatabaseManager(db)
    stats_data = manager.get_stats()
    table = Table(title="Database Statistics")
    table.add_column("Metric", style="cyan")
    table.add_column("Count", style="green")

    for k, v in stats_data.items():
        table.add_row(k.replace("_", " ").title(), str(v))

    console.print(table)


if __name__ == "__main__":
    cli()
