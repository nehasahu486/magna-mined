"""
SQLite Database Layer for Magna-Mined.
Maintains the complete 12,827 M. oryzae proteome catalog and the automated evidence ledger.
"""

import sqlite3
import csv
from pathlib import Path
from typing import Dict, Any, Optional, List
from datetime import datetime

from .schema import GeneEvidenceRecord, PhenotypeCategory


class DatabaseManager:
    """Manages SQLite storage for genes, clusters, and automated extractions."""

    def __init__(self, db_path: str = "magnaporthe.db"):
        self.db_path = Path(db_path)
        self.init_database()

    def get_connection(self) -> sqlite3.Connection:
        """Create and return a database connection with foreign key support."""
        conn = sqlite3.connect(self.db_path)
        conn.execute("PRAGMA foreign_keys = ON")
        conn.row_factory = sqlite3.Row
        return conn

    def init_database(self) -> None:
        """Initialize the SQLite schema."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        with self.get_connection() as conn:
            cursor = conn.cursor()

            # Master genes table (matches 12,827 genome catalog)
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS genes (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    gene TEXT NOT NULL UNIQUE,
                    phylostratum TEXT,
                    class TEXT,
                    size_aa INTEGER,
                    ensembl_link TEXT,
                    uniprot TEXT,
                    gene_name TEXT,
                    process_affected TEXT DEFAULT 'not reported',
                    baseline_process TEXT DEFAULT 'not reported',
                    knockout_reported BOOLEAN DEFAULT 0,
                    knockout_phenotype TEXT,
                    additional_functions TEXT,
                    additional_function_references TEXT,
                    detailed_mutant_phenotype TEXT,
                    pmid_or_link_to_publication TEXT,
                    pmid_or_link_to_publication_1 TEXT,
                    blast2go_annotation TEXT,
                    pfam TEXT
                )
            """)

            # Functional clusters
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS clusters (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT NOT NULL UNIQUE,
                    description TEXT
                )
            """)

            # Many-to-many cluster assignments
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS cluster_assignments (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    gene_id INTEGER NOT NULL,
                    cluster_id INTEGER NOT NULL,
                    FOREIGN KEY (gene_id) REFERENCES genes(id) ON DELETE CASCADE,
                    FOREIGN KEY (cluster_id) REFERENCES clusters(id) ON DELETE CASCADE,
                    UNIQUE(gene_id, cluster_id)
                )
            """)

            # Automated extractions ledger (provenance & audit log)
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS extractions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    gene_id INTEGER NOT NULL,
                    pmid TEXT,
                    pmcid TEXT,
                    doi TEXT,
                    process_affected TEXT NOT NULL,
                    detailed_phenotype TEXT NOT NULL,
                    experimental_method TEXT NOT NULL,
                    complementation_tested BOOLEAN DEFAULT 0,
                    confidence_tier TEXT NOT NULL,
                    verbatim_quote TEXT NOT NULL,
                    paper_title TEXT,
                    publication_year INTEGER,
                    species_side_note TEXT,
                    extracted_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (gene_id) REFERENCES genes(id) ON DELETE CASCADE
                )
            """)

            # External & baseline gene aliases (UniProt / Ensembl / literature synonyms)
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS gene_aliases (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    gene_id INTEGER NOT NULL,
                    alias TEXT NOT NULL,
                    source TEXT DEFAULT 'uniprot',
                    FOREIGN KEY (gene_id) REFERENCES genes(id) ON DELETE CASCADE,
                    UNIQUE(gene_id, alias)
                )
            """)

            # Processed literature tracking ledger
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS processed_papers (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pmid TEXT UNIQUE,
                    pmcid TEXT,
                    doi TEXT,
                    title TEXT,
                    pub_year INTEGER,
                    is_open_access BOOLEAN DEFAULT 0,
                    is_experimental BOOLEAN DEFAULT 0,
                    triage_rationale TEXT,
                    genes_found INTEGER DEFAULT 0,
                    checked_additional BOOLEAN DEFAULT 0,
                    processed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)

            # Lightweight migrations for databases created before Phase 5 / Phase 6
            for table in ("extractions", "processed_papers"):
                cursor.execute(f"PRAGMA table_info({table})")
                existing_cols = {r[1] for r in cursor.fetchall()}
                if "evidence_source" not in existing_cols:
                    default = "'jats_xml'" if table == "extractions" else "NULL"
                    cursor.execute(f"ALTER TABLE {table} ADD COLUMN evidence_source TEXT DEFAULT {default}")
                if table == "extractions":
                    if "reviewed_by_human" not in existing_cols:
                        cursor.execute("ALTER TABLE extractions ADD COLUMN reviewed_by_human BOOLEAN DEFAULT 0")
                    if "review_notes" not in existing_cols:
                        cursor.execute("ALTER TABLE extractions ADD COLUMN review_notes TEXT")
                    if "species_side_note" not in existing_cols:
                        cursor.execute("ALTER TABLE extractions ADD COLUMN species_side_note TEXT")
                if table == "processed_papers":
                    if "checked_additional" not in existing_cols:
                        cursor.execute("ALTER TABLE processed_papers ADD COLUMN checked_additional BOOLEAN DEFAULT 0")

            cursor.execute("PRAGMA table_info(genes)")
            gene_cols = {r[1] for r in cursor.fetchall()}
            if "conflict_status" not in gene_cols:
                cursor.execute("ALTER TABLE genes ADD COLUMN conflict_status TEXT DEFAULT 'none'")
            if "consensus_tier" not in gene_cols:
                cursor.execute("ALTER TABLE genes ADD COLUMN consensus_tier TEXT")
            if "baseline_process" not in gene_cols:
                cursor.execute("ALTER TABLE genes ADD COLUMN baseline_process TEXT DEFAULT 'not reported'")
                cursor.execute("UPDATE genes SET baseline_process = process_affected WHERE baseline_process IS NULL")
            if "knockout_reported" not in gene_cols:
                cursor.execute("ALTER TABLE genes ADD COLUMN knockout_reported BOOLEAN DEFAULT 0")
            if "knockout_phenotype" not in gene_cols:
                cursor.execute("ALTER TABLE genes ADD COLUMN knockout_phenotype TEXT")
            if "additional_functions" not in gene_cols:
                cursor.execute("ALTER TABLE genes ADD COLUMN additional_functions TEXT")
            if "additional_function_references" not in gene_cols:
                cursor.execute("ALTER TABLE genes ADD COLUMN additional_function_references TEXT")
            if "other_fungi_notes" not in gene_cols:
                cursor.execute("ALTER TABLE genes ADD COLUMN other_fungi_notes TEXT")

            # Indexes for high-throughput queries
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_genes_gene ON genes(gene)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene_aliases_alias ON gene_aliases(alias)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_genes_name ON genes(gene_name)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_genes_process ON genes(process_affected)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_extractions_gene ON extractions(gene_id)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_extractions_pmid ON extractions(pmid)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_processed_pmid ON processed_papers(pmid)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_processed_pmcid ON processed_papers(pmcid)")

            conn.commit()

    def seed_from_excel(self, excel_path: str) -> int:
        """Seed genes table from Magnagenes.xlsx (Zenodo record 4647766)."""
        import openpyxl

        path = Path(excel_path)
        if not path.exists():
            raise FileNotFoundError(f"Excel file not found: {excel_path}")

        wb = openpyxl.load_workbook(str(path), read_only=True)
        sheet = wb[wb.sheetnames[0]]

        inserted_count = 0
        rows = sheet.iter_rows(values_only=True)
        header = next(rows)  # skip header

        with self.get_connection() as conn:
            cursor = conn.cursor()
            for row in rows:
                if not row or not row[0]:
                    continue
                gene_id = str(row[0]).strip()
                if not gene_id.startswith("MGG_"):
                    continue

                clean_str = lambda val: str(val).strip() if val is not None and str(val).strip() else None
                size_aa = int(row[1]) if len(row) > 1 and str(row[1]).isdigit() else None
                ensembl = clean_str(row[2]) if len(row) > 2 else None
                uniprot = clean_str(row[3]) if len(row) > 3 else None
                gene_name = clean_str(row[4]) if len(row) > 4 else None
                process = clean_str(row[5]) if len(row) > 5 else "not reported"
                process = process or "not reported"
                detailed = clean_str(row[6]) if len(row) > 6 else None
                pmid_1 = clean_str(row[7]) if len(row) > 7 else None
                pmid_2 = clean_str(row[8]) if len(row) > 8 else None
                b2g = clean_str(row[9]) if len(row) > 9 else None
                pfam = clean_str(row[10]) if len(row) > 10 else None

                cursor.execute("""
                    INSERT OR IGNORE INTO genes (
                        gene, size_aa, ensembl_link, uniprot,
                        gene_name, process_affected, baseline_process, detailed_mutant_phenotype,
                        pmid_or_link_to_publication, pmid_or_link_to_publication_1,
                        blast2go_annotation, pfam
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    gene_id, size_aa, ensembl, uniprot,
                    gene_name, process, process, detailed,
                    pmid_1, pmid_2, b2g, pfam
                ))
                inserted_count += 1
            conn.commit()
        wb.close()
        return inserted_count

    def seed_from_csv(self, magna_genes_csv: str) -> int:
        """Seed genes table from MagnaGenes.csv."""
        csv_path = Path(magna_genes_csv)
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {magna_genes_csv}")

        clean_str = lambda val: str(val).strip() if val is not None and str(val).strip() else None
        inserted_count = 0
        with open(csv_path, mode="r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            with self.get_connection() as conn:
                cursor = conn.cursor()
                for row in reader:
                    gene_id = clean_str(row.get("Gene") or row.get("gene"))
                    if not gene_id:
                        continue

                    size_raw = clean_str(row.get("size (aa)"))
                    csv_process = clean_str(row.get("process affected")) or "not reported"
                    cursor.execute("""
                        INSERT OR IGNORE INTO genes (
                            gene, phylostratum, class, size_aa, ensembl_link, uniprot,
                            gene_name, process_affected, baseline_process, detailed_mutant_phenotype,
                            pmid_or_link_to_publication, pmid_or_link_to_publication_1,
                            blast2go_annotation, pfam
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """, (
                        gene_id,
                        clean_str(row.get("Phylostratum")),
                        clean_str(row.get("Class")),
                        int(size_raw) if size_raw and size_raw.isdigit() else None,
                        clean_str(row.get("Ensembl link")),
                        clean_str(row.get("Uniprot")),
                        clean_str(row.get("gene name")),
                        csv_process,
                        csv_process,
                        clean_str(row.get("detailed mutant phenotype")),
                        clean_str(row.get("PMID or link to publication")),
                        clean_str(row.get("PMID or link to publication.1")),
                        clean_str(row.get("BLAST2GO annotation")),
                        clean_str(row.get("Pfam")),
                    ))
                    inserted_count += 1
                conn.commit()
        return inserted_count

    @property
    def linker(self):
        """Lazy-loaded EntityLinker for alias and mutant disambiguation."""
        if not hasattr(self, "_linker") or self._linker is None:
            from .entity_linker import EntityLinker
            self._linker = EntityLinker(db=self)
        return self._linker

    def get_gene(self, gene_query: str) -> Optional[Dict[str, Any]]:
        """Look up gene by systematic ID (MGG_#####), mutant, alias, or standard name."""
        clean_query = gene_query.strip()
        resolved_mgg = self.linker.resolve(clean_query)
        target = resolved_mgg or clean_query

        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT * FROM genes 
                WHERE UPPER(gene) = UPPER(?) 
                   OR UPPER(gene_name) = UPPER(?)
                   OR UPPER(gene_name) = UPPER(?)
                LIMIT 1
            """, (target, clean_query, clean_query.replace("MO", "")))
            row = cursor.fetchone()
            return dict(row) if row else None

    _TIER_RANK = {
        "tier_1_high": 4,
        "tier_2_high": 3,
        "tier_3_medium": 2,
        "tier_4_low": 1,
    }

    _KNOCKOUT_METHODS = {"knockout/deletion", "CRISPR/Cas9", "complementation"}

    _PROCESS_SYNONYMS = {
        "pathogenicity": "virulence",
        "virulenece": "virulence",
        "effector": "virulence",
        "elicitor": "virulence",
        "appresorium formation": "appressorium formation",
        "penetration": "appressorium formation",
        "turgor": "appressorium formation",
        "adhesion": "appressorium formation",
        "condiation": "conidiation",
        "pigmentation": "stress response",
        "pimentation": "stress response",
        "growth": "growth/development",
        "germination": "growth/development",
        "sexual development": "growth/development",
        "putative essential gene": "growth/development",
        "none": "dispensable / wild-type-like",
    }

    @classmethod
    def normalize_process_set(cls, proc_str: Optional[str]) -> set:
        """Normalize comma-separated biological processes for overlap/novelty checks."""
        if not proc_str:
            return set()
        out = set()
        for tok in str(proc_str).split(","):
            norm = tok.strip().lower()
            if not norm or norm in ("not reported", "mixed outcomes"):
                continue
            out.add(cls._PROCESS_SYNONYMS.get(norm, norm))
        return out

    def save_evidence(self, evidence: GeneEvidenceRecord, preserve_baseline_primary: bool = False) -> int:
        """
        Save a validated evidence record and update master gene record:
        - Updates `knockout_reported` / `knockout_phenotype` when knockout/CRISPR evidence exists.
        - For already-classified baseline genes, captures newly discovered biological roles in
          `additional_functions` and `additional_function_references`.
        - Updates `process_affected` using tier precedence (unless `preserve_baseline_primary=True`).
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                SELECT id, process_affected, baseline_process, consensus_tier,
                       knockout_reported, knockout_phenotype,
                       additional_functions, additional_function_references,
                       other_fungi_notes
                FROM genes WHERE UPPER(gene) = UPPER(?)
                """,
                (evidence.gene_id,)
            )
            row = cursor.fetchone()
            if not row:
                raise ValueError(f"Gene {evidence.gene_id} not found in database catalog.")

            db_gene_id = row["id"]

            # Build combined process string for pleiotropic phenotypes
            all_procs = [evidence.process_affected.value]
            if hasattr(evidence, "secondary_processes") and evidence.secondary_processes:
                for p in evidence.secondary_processes:
                    val = p.value if hasattr(p, "value") else str(p)
                    if val not in all_procs and val != "not reported":
                        all_procs.append(val)
            combined_process_str = ", ".join(all_procs)

            # Avoid inserting duplicate extraction rows for the same gene, paper, and quote
            existing = cursor.execute(
                """
                SELECT id FROM extractions
                WHERE gene_id = ? AND COALESCE(pmid, '') = COALESCE(?, '')
                  AND COALESCE(pmcid, '') = COALESCE(?, '')
                  AND verbatim_quote = ?
                LIMIT 1
                """,
                (db_gene_id, evidence.pmid, evidence.pmcid, evidence.verbatim_quote),
            ).fetchone()
            if existing:
                return existing["id"]

            # Insert extraction record
            cursor.execute("""
                INSERT INTO extractions (
                    gene_id, pmid, pmcid, doi, process_affected, detailed_phenotype,
                    experimental_method, complementation_tested, confidence_tier,
                    verbatim_quote, paper_title, publication_year, evidence_source,
                    species_side_note
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                db_gene_id,
                evidence.pmid,
                evidence.pmcid,
                evidence.doi,
                combined_process_str,
                evidence.detailed_phenotype,
                evidence.experimental_method.value,
                1 if evidence.complementation_tested else 0,
                evidence.confidence_tier.value,
                evidence.verbatim_quote,
                evidence.paper_title,
                evidence.publication_year,
                evidence.evidence_source.value,
                getattr(evidence, "species_side_note", None),
            ))
            extraction_id = cursor.lastrowid

            disp_tag = "[dispensable] " if getattr(evidence, "is_dispensable", False) else ""

            # 1. Track knockout status & knockout phenotype explicitly
            ko_reported = int(row["knockout_reported"] or 0)
            ko_phenotype = row["knockout_phenotype"]
            if evidence.experimental_method.value in self._KNOCKOUT_METHODS:
                ko_reported = 1
                ref_label = f"PMID:{evidence.pmid}" if evidence.pmid else (evidence.doi or "N/A")
                new_ko = f"[{evidence.confidence_tier.value}] {combined_process_str}: {disp_tag}{evidence.detailed_phenotype} ({ref_label})"
                if not ko_phenotype:
                    ko_phenotype = new_ko
                elif new_ko not in ko_phenotype:
                    ko_phenotype = f"{ko_phenotype}; {new_ko}"

            # 2. Track additional functions discovered for already-classified baseline genes
            base_proc = (row["baseline_process"] or "not reported").strip()
            base_is_classified = base_proc.lower() not in ("not reported", "")
            add_funcs_str = row["additional_functions"]
            add_refs_str = row["additional_function_references"]

            if base_is_classified:
                base_set = self.normalize_process_set(base_proc)
                existing_add_list = [p.strip() for p in (add_funcs_str or "").split(",") if p.strip()]
                existing_add_set = {p.lower() for p in existing_add_list}

                new_procs = [
                    p for p in all_procs
                    if p.lower() not in ("not reported", "dispensable / wild-type-like")
                    and p.lower() not in base_set
                    and p.lower() not in existing_add_set
                ]
                if new_procs:
                    updated_add_list = existing_add_list + new_procs
                    add_funcs_str = ", ".join(updated_add_list)
                    ref_label = f"PMID:{evidence.pmid}" if evidence.pmid else (evidence.doi or "N/A")
                    yr_label = f" ({evidence.publication_year})" if evidence.publication_year else ""
                    new_ref = f"[{', '.join(new_procs)}] {ref_label}{yr_label}: {evidence.detailed_phenotype}"
                    if not add_refs_str:
                        add_refs_str = new_ref
                    elif new_ref not in add_refs_str:
                        add_refs_str = f"{add_refs_str}; {new_ref}"

            # 2b. Track cross-species pathogenicity side notes (e.g. Colletotrichum)
            other_notes_str = row["other_fungi_notes"]
            side_note = (getattr(evidence, "species_side_note", None) or "").strip()
            if side_note:
                ref_label = f"PMID:{evidence.pmid}" if evidence.pmid else (evidence.doi or "N/A")
                formatted_note = f"{side_note} ({ref_label})"
                if not other_notes_str:
                    other_notes_str = formatted_note
                elif formatted_note not in other_notes_str and side_note not in other_notes_str:
                    other_notes_str = f"{other_notes_str}; {formatted_note}"

            # 3. Update primary process_affected using tier precedence (unless preserving baseline primary)
            current_process = (row["process_affected"] or "not reported").strip().lower()
            existing_rank = self._TIER_RANK.get(row["consensus_tier"] or "", 0)
            new_rank = self._TIER_RANK.get(evidence.confidence_tier.value, 0)
            should_promote = (
                not preserve_baseline_primary
                and (
                    (current_process in ("not reported", "none", "") and evidence.process_affected != PhenotypeCategory.NOT_REPORTED)
                    or (evidence.process_affected != PhenotypeCategory.NOT_REPORTED and new_rank >= existing_rank)
                )
            )

            if should_promote:
                cursor.execute("""
                    UPDATE genes SET
                        process_affected = ?,
                        consensus_tier = ?,
                        knockout_reported = ?,
                        knockout_phenotype = ?,
                        additional_functions = ?,
                        additional_function_references = ?,
                        other_fungi_notes = ?,
                        detailed_mutant_phenotype = COALESCE(detailed_mutant_phenotype || '; ', '') || ?,
                        pmid_or_link_to_publication = COALESCE(pmid_or_link_to_publication, ?)
                    WHERE id = ?
                """, (
                    combined_process_str,
                    evidence.confidence_tier.value,
                    ko_reported,
                    ko_phenotype,
                    add_funcs_str,
                    add_refs_str,
                    other_notes_str,
                    f"[{evidence.confidence_tier.value}] {disp_tag}{evidence.detailed_phenotype}",
                    evidence.pmid or evidence.doi,
                    db_gene_id
                ))
            else:
                cursor.execute("""
                    UPDATE genes SET
                        knockout_reported = ?,
                        knockout_phenotype = ?,
                        additional_functions = ?,
                        additional_function_references = ?,
                        other_fungi_notes = ?
                    WHERE id = ?
                """, (
                    ko_reported,
                    ko_phenotype,
                    add_funcs_str,
                    add_refs_str,
                    other_notes_str,
                    db_gene_id
                ))

            conn.commit()
            return extraction_id

    def add_gene_aliases_bulk(self, pairs: List[tuple], source: str = "uniprot") -> int:
        """
        Insert `(canonical_mgg_id, alias)` mappings into `gene_aliases`.
        Returns the number of new alias rows inserted.
        """
        if not pairs:
            return 0
        with self.get_connection() as conn:
            mgg_to_pk = {
                r["gene"].strip().upper(): r["id"]
                for r in conn.execute("SELECT id, gene FROM genes").fetchall()
            }
            before = conn.total_changes
            to_insert = []
            for mgg_id, alias in pairs:
                pk = mgg_to_pk.get(str(mgg_id).strip().upper())
                al = str(alias).strip()
                if pk and al and len(al) >= 2:
                    to_insert.append((pk, al, source))
            conn.executemany(
                "INSERT OR IGNORE INTO gene_aliases (gene_id, alias, source) VALUES (?, ?, ?)",
                to_insert,
            )
            conn.commit()
            inserted = conn.total_changes - before
        self._linker = None
        return inserted

    def get_all_gene_aliases(self) -> Dict[str, List[str]]:
        """Return `{canonical_mgg_upper: [alias, ...]}` from the `gene_aliases` table."""
        with self.get_connection() as conn:
            rows = conn.execute(
                "SELECT g.gene, a.alias FROM gene_aliases a JOIN genes g ON g.id = a.gene_id"
            ).fetchall()
        out: Dict[str, List[str]] = {}
        for r in rows:
            out.setdefault(r["gene"].strip().upper(), []).append(r["alias"].strip())
        return out

    def delete_extractions_for_paper(
        self, pmid: Optional[str] = None, pmcid: Optional[str] = None, doi: Optional[str] = None
    ) -> int:
        """Remove prior extraction rows for a paper (used when re-ingesting with better full text)."""
        clauses, params = [], []
        for col, val in (("pmid", pmid), ("pmcid", pmcid), ("doi", doi)):
            if val:
                clauses.append(f"{col} = ?")
                params.append(val)
        if not clauses:
            return 0
        with self.get_connection() as conn:
            cur = conn.execute(f"DELETE FROM extractions WHERE {' OR '.join(clauses)}", params)
            conn.commit()
            return cur.rowcount

    def is_paper_processed(
        self,
        pmid: Optional[str] = None,
        pmcid: Optional[str] = None,
        require_additional_check: bool = False
    ) -> bool:
        """Check if a paper has already been triaged/extracted (and optionally double-checked)."""
        if not pmid and not pmcid:
            return False

        extra_clause = " AND checked_additional = 1" if require_additional_check else ""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            if pmid and pmcid:
                cursor.execute(
                    f"SELECT 1 FROM processed_papers WHERE (pmid = ? OR pmcid = ?){extra_clause} LIMIT 1",
                    (pmid, pmcid)
                )
            elif pmid:
                cursor.execute(f"SELECT 1 FROM processed_papers WHERE pmid = ?{extra_clause} LIMIT 1", (pmid,))
            else:
                cursor.execute(f"SELECT 1 FROM processed_papers WHERE pmcid = ?{extra_clause} LIMIT 1", (pmcid,))
            return cursor.fetchone() is not None

    def record_processed_paper(
        self,
        pmid: Optional[str],
        pmcid: Optional[str],
        doi: Optional[str],
        title: str,
        pub_year: Optional[int],
        is_open_access: bool,
        is_experimental: bool,
        triage_rationale: str = "",
        genes_found: int = 0,
        evidence_source: Optional[str] = None,
        checked_additional: bool = False
    ) -> int:
        """Record or update a paper in the processed-paper tracking ledger."""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO processed_papers (
                    pmid, pmcid, doi, title, pub_year,
                    is_open_access, is_experimental, triage_rationale, genes_found,
                    evidence_source, checked_additional
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(pmid) DO UPDATE SET
                    pmcid = COALESCE(excluded.pmcid, processed_papers.pmcid),
                    doi = COALESCE(excluded.doi, processed_papers.doi),
                    is_experimental = excluded.is_experimental,
                    triage_rationale = excluded.triage_rationale,
                    genes_found = MAX(processed_papers.genes_found, excluded.genes_found),
                    evidence_source = COALESCE(excluded.evidence_source, processed_papers.evidence_source),
                    checked_additional = MAX(processed_papers.checked_additional, excluded.checked_additional),
                    processed_at = CURRENT_TIMESTAMP
            """, (
                pmid,
                pmcid,
                doi,
                title,
                pub_year,
                1 if is_open_access else 0,
                1 if is_experimental else 0,
                triage_rationale,
                genes_found,
                evidence_source,
                1 if checked_additional else 0
            ))
            conn.commit()
            return cursor.lastrowid

    def get_unclassified_gene_ids(self) -> set:
        """Return uppercase MGG IDs whose baseline status in MagnaGenes was unclassified ('Not reported')."""
        with self.get_connection() as conn:
            rows = conn.execute(
                """
                SELECT UPPER(gene) AS gene
                FROM genes
                WHERE LOWER(TRIM(COALESCE(baseline_process, 'not reported'))) IN ('not reported', '')
                """
            ).fetchall()
            return {r["gene"] for r in rows}

    def is_gene_unclassified_in_baseline(self, gene_id: str) -> bool:
        """Return True if gene was unclassified ('Not reported') in the baseline MagnaGenes catalog."""
        with self.get_connection() as conn:
            row = conn.execute(
                """
                SELECT COALESCE(baseline_process, 'not reported') AS base_proc
                FROM genes
                WHERE UPPER(gene) = UPPER(?)
                LIMIT 1
                """,
                (gene_id.strip(),)
            ).fetchone()
            if not row:
                return False
            return row["base_proc"].strip().lower() in ("not reported", "")

    def get_stats(self) -> Dict[str, Any]:
        """Return summary statistics of the database."""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) as total FROM genes")
            total_genes = cursor.fetchone()["total"]

            cursor.execute("SELECT COUNT(*) as with_phenotype FROM genes WHERE LOWER(TRIM(process_affected)) NOT IN ('not reported', 'none', '')")
            with_phenotype = cursor.fetchone()["with_phenotype"]

            cursor.execute("SELECT COUNT(*) as base_unclassified FROM genes WHERE LOWER(TRIM(COALESCE(baseline_process, 'not reported'))) IN ('not reported', '')")
            base_unclassified = cursor.fetchone()["base_unclassified"]

            cursor.execute("""
                SELECT COUNT(*) as newly_classified FROM genes
                WHERE LOWER(TRIM(COALESCE(baseline_process, 'not reported'))) IN ('not reported', '')
                  AND EXISTS (SELECT 1 FROM extractions e WHERE e.gene_id = genes.id)
            """)
            newly_classified = cursor.fetchone()["newly_classified"]

            cursor.execute("SELECT COUNT(*) as with_lit FROM genes WHERE pmid_or_link_to_publication IS NOT NULL")
            with_lit = cursor.fetchone()["with_lit"]

            cursor.execute("SELECT COUNT(*) as total_extractions FROM extractions")
            total_extractions = cursor.fetchone()["total_extractions"]

            cursor.execute("SELECT COUNT(*) as papers_triaged FROM processed_papers")
            papers_triaged = cursor.fetchone()["papers_triaged"]

            cursor.execute("SELECT COUNT(*) as exp_papers FROM processed_papers WHERE is_experimental = 1")
            exp_papers = cursor.fetchone()["exp_papers"]

            cursor.execute("SELECT COUNT(*) as conflicts FROM genes WHERE conflict_status = 'conflict'")
            conflicts = cursor.fetchone()["conflicts"]

            cursor.execute("SELECT COUNT(*) as ko_genes FROM genes WHERE knockout_reported = 1")
            ko_genes = cursor.fetchone()["ko_genes"]

            cursor.execute("SELECT COUNT(*) as add_funcs FROM genes WHERE additional_functions IS NOT NULL AND TRIM(additional_functions) != ''")
            add_funcs = cursor.fetchone()["add_funcs"]

            return {
                "total_genes": total_genes,
                "baseline_unclassified": base_unclassified,
                "newly_classified_from_unclassified": newly_classified,
                "knockout_reported_genes": ko_genes,
                "genes_with_additional_functions": add_funcs,
                "genes_with_phenotype": with_phenotype,
                "genes_with_literature": with_lit,
                "total_evidence_extractions": total_extractions,
                "papers_triaged": papers_triaged,
                "experimental_papers": exp_papers,
                "conflicts_flagged": conflicts
            }
