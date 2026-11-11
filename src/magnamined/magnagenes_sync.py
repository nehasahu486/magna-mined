"""
MagnaGenes v2.0 Integration & Dated Export/Archive Manager for Magna-Mined.

Produces long-format `(Gene, variable, value)` TSV files matching
`magnagenes_v2/data-raw/source/annotation/fin_data.tsv`, using dated filenames
(`fin_data_YYYY-MM-DD.tsv`) and maintaining a SHA-256 deduplicated `archive/`
folder for previously used files.
"""

import csv
import hashlib
import re
import shutil
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .database import DatabaseManager

FIN_DATA_PATTERN = re.compile(r"^fin_data(?:_(\d{4}-\d{2}-\d{2}))?(?:_[0-9a-f]{8})?\.tsv$")


def sha256_file(path: Path) -> str:
    """Compute SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _format_pubmed_url(raw_ref: Optional[str]) -> Optional[str]:
    if not raw_ref:
        return None
    s = str(raw_ref).strip()
    if not s:
        return None
    if s.startswith("http://") or s.startswith("https://"):
        return s
    if s.isdigit():
        return f"https://pubmed.ncbi.nlm.nih.gov/{s}/"
    if s.startswith("10."):
        return f"https://doi.org/{s}"
    return s


def _format_quote(e: Dict[str, Any]) -> Optional[str]:
    """Format one extraction as "[tier] PMID:x (year): quote" on a single line.
    Tabs, newlines and pipes are replaced so the TSV stays one row per value."""
    q = re.sub(r"\s+", " ", str(e.get("verbatim_quote") or "")).replace("|", "/").strip()
    if not q:
        return None
    src = f"PMID:{e['pmid']}" if e.get("pmid") else "no PMID"
    if e.get("publication_year"):
        src += f" ({e['publication_year']})"
    tier = e.get("confidence_tier") or "untiered"
    return f"[{tier}] {src}: {q}"


def find_latest_fin_data(directory: Path) -> Optional[Path]:
    """
    Return the latest dated `fin_data_YYYY-MM-DD.tsv` (or fallback `fin_data.tsv`)
    in `directory`, ignoring subdirectories such as `archive/`.
    """
    if not directory.exists():
        return None
    candidates: List[Tuple[str, float, Path]] = []
    for p in directory.iterdir():
        if not p.is_file():
            continue
        m = FIN_DATA_PATTERN.match(p.name)
        if m:
            date_key = m.group(1) or "0000-00-00"
            candidates.append((date_key, p.stat().st_mtime, p))
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return candidates[0][2]


def archive_old_fin_data(
    output_dir: Path,
    active_file: Path,
    archive_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """
    Move any older `fin_data*.tsv` files in `output_dir` (except `active_file`)
    into `archive_dir` (`output_dir/archive/` by default), deduplicating by SHA-256
    so `archive/` never stores duplicate copies.
    """
    arch = archive_dir or (output_dir / "archive")
    arch.mkdir(parents=True, exist_ok=True)

    # Index SHA-256 hashes already present in archive/
    existing_hashes: Dict[str, Path] = {}
    for p in arch.iterdir():
        if p.is_file() and p.suffix == ".tsv":
            existing_hashes[sha256_file(p)] = p

    moved: List[str] = []
    deduplicated: List[str] = []

    active_resolved = active_file.resolve()
    for p in sorted(output_dir.iterdir()):
        if not p.is_file() or p.resolve() == active_resolved:
            continue
        if not FIN_DATA_PATTERN.match(p.name):
            continue

        digest = sha256_file(p)
        if digest in existing_hashes:
            # Identical content already archived -> remove duplicate instead of storing twice
            p.unlink()
            deduplicated.append(p.name)
            continue

        dest = arch / p.name
        if dest.exists():
            # Same filename but different content -> preserve with short hash suffix
            dest = arch / f"{p.stem}_{digest[:8]}{p.suffix}"
            if dest.exists():
                p.unlink()
                deduplicated.append(p.name)
                continue

        shutil.move(str(p), str(dest))
        existing_hashes[digest] = dest
        moved.append(dest.name)

    return {
        "archive_dir": str(arch),
        "archived_files": moved,
        "duplicates_removed": deduplicated,
    }


class MagnaGenesSyncExporter:
    """Exports Magna-Mined curations into MagnaGenes v2.0 `fin_data_YYYY-MM-DD.tsv` format."""

    def __init__(self, db: Optional[DatabaseManager] = None, db_path: str = "magnaporthe.db"):
        self.db = db or DatabaseManager(db_path)

    def _get_gene_updates(self, scope: str = "curated") -> List[Dict[str, Any]]:
        """Fetch genes with automated curations/knockouts/additional functions (or full catalog)."""
        with self.db.get_connection() as conn:
            if scope == "full":
                sql = """
                    SELECT g.*,
                           (SELECT GROUP_CONCAT(DISTINCT e.pmid) FROM extractions e WHERE e.gene_id = g.id AND e.pmid IS NOT NULL) AS ext_pmids
                    FROM genes g
                    ORDER BY g.gene ASC
                """
            else:
                sql = """
                    SELECT g.*,
                           (SELECT GROUP_CONCAT(DISTINCT e.pmid) FROM extractions e WHERE e.gene_id = g.id AND e.pmid IS NOT NULL) AS ext_pmids
                    FROM genes g
                    WHERE EXISTS (SELECT 1 FROM extractions e WHERE e.gene_id = g.id)
                       OR g.knockout_reported = 1
                       OR (g.additional_functions IS NOT NULL AND TRIM(g.additional_functions) != '')
                    ORDER BY g.gene ASC
                """
            rows = [dict(r) for r in conn.execute(sql).fetchall()]
            quotes: Dict[int, List[Dict[str, Any]]] = {}
            for e in conn.execute(
                """SELECT gene_id, pmid, publication_year, confidence_tier, verbatim_quote
                   FROM extractions
                   WHERE verbatim_quote IS NOT NULL AND TRIM(verbatim_quote) != ''
                   ORDER BY gene_id, confidence_tier, publication_year DESC, pmid"""
            ).fetchall():
                quotes.setdefault(e["gene_id"], []).append(dict(e))
            for g in rows:
                g["evidence_quotes"] = quotes.get(g["id"], [])
            return rows

    def _gene_to_long_rows(self, g: Dict[str, Any], include_baseline_meta: bool = False) -> List[Tuple[str, str, str]]:
        """Convert a gene record into `(Gene, variable, value)` tuples compatible with `prepare_data.R`."""
        gid = g["gene"]
        rows: List[Tuple[str, str, str]] = []

        def _add(var: str, val: Optional[str]) -> None:
            if val is not None and str(val).strip():
                rows.append((gid, var, str(val).strip()))

        if g.get("gene_name"):
            _add("Gene.Name", g["gene_name"])

        proc = (g.get("process_affected") or "").strip()
        if proc and proc.lower() != "not reported":
            _add("Process.Affected", proc)

        if g.get("detailed_mutant_phenotype"):
            _add("Mutant.Phenotype", g["detailed_mutant_phenotype"])

        primary_ref = _format_pubmed_url(g.get("pmid_or_link_to_publication"))
        if primary_ref:
            _add("PMID.or.link", primary_ref)
        elif g.get("ext_pmids"):
            for pm in str(g["ext_pmids"]).split(","):
                url = _format_pubmed_url(pm)
                if url:
                    _add("PMID.or.link", url)

        if g.get("pmid_or_link_to_publication_1"):
            _add("Other.links.to.publication", _format_pubmed_url(g["pmid_or_link_to_publication_1"]))

        if g.get("knockout_reported"):
            _add("Knockout.Reported", "Yes")
        if g.get("knockout_phenotype"):
            _add("Knockout.Phenotype", g["knockout_phenotype"])
        if g.get("additional_functions"):
            _add("Additional.Functions", g["additional_functions"])
        if g.get("additional_function_references"):
            _add("Other.links.to.publication", g["additional_function_references"])
        if g.get("other_fungi_notes"):
            _add("Other.Fungi.Notes", g["other_fungi_notes"])
        # One row per supporting excerpt: "[tier] PMID:x (year): quote"
        for e in g.get("evidence_quotes") or []:
            _add("Evidence.Quotes", _format_quote(e))

        if include_baseline_meta:
            if g.get("size_aa"):
                _add("Size.(aa)", str(g["size_aa"]))
            _add("BLAST2GO.Annotation", g.get("blast2go_annotation"))
            _add("Pfam", g.get("pfam"))
            _add("Uniprot", g.get("uniprot"))
            _add("Ensembl.Link", g.get("ensembl_link"))

        return rows

    def build_fin_data_rows(
        self,
        scope: str = "curated",
        base_fin_data: Optional[Path] = None,
    ) -> Tuple[List[Tuple[str, str, str]], int]:
        """
        Build `(Gene, variable, value)` rows.
        If `base_fin_data` is given, merges updates cleanly into the existing
        `magnagenes_v2` `fin_data.tsv` (replacing `'Not reported'` values and
        preserving all existing columns like `InterPro`, `KEGG.Pathway`, `AlphaFoldDB`).
        Returns `(rows, updated_gene_count)`.
        """
        updates = self._get_gene_updates(scope=scope)
        update_by_gene = {g["gene"].upper(): g for g in updates}

        if base_fin_data is None:
            all_rows: List[Tuple[str, str, str]] = []
            for g in updates:
                all_rows.extend(self._gene_to_long_rows(g, include_baseline_meta=(scope == "full")))
            return all_rows, len(updates)

        # Merge mode: load existing fin_data.tsv and patch updated genes in place
        if not base_fin_data.exists():
            raise FileNotFoundError(f"Base fin_data file not found: {base_fin_data}")

        merged: List[Tuple[str, str, str]] = []
        seen_triples = set()
        replace_vars_for_newly_classified = {"Process.Affected", "Mutant.Phenotype", "PMID.or.link"}

        newly_classified_set = {
            gid for gid, g in update_by_gene.items()
            if (g.get("baseline_process") or "not reported").strip().lower() in ("not reported", "")
        }

        with open(base_fin_data, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                gene = (r.get("Gene") or "").strip()
                var = (r.get("variable") or "").strip()
                val = (r.get("value") or "").strip()
                if not gene or not var:
                    continue
                g_up = gene.upper()
                if g_up in newly_classified_set and var in replace_vars_for_newly_classified:
                    if val.lower() in ("not reported", ""):
                        continue
                triple = (gene, var, val)
                if triple not in seen_triples:
                    seen_triples.add(triple)
                    merged.append(triple)

        for g in updates:
            for triple in self._gene_to_long_rows(g, include_baseline_meta=False):
                if triple not in seen_triples:
                    seen_triples.add(triple)
                    merged.append(triple)

        return merged, len(updates)

    def export_dated_fin_data(
        self,
        output_dir: str = "exports/magnagenes",
        export_date: Optional[str] = None,
        base_fin_data: Optional[str] = None,
        merge_existing: bool = False,
        scope: str = "curated",
        archive_existing: bool = True,
    ) -> Dict[str, Any]:
        """
        Write `fin_data_YYYY-MM-DD.tsv` into `output_dir`, archive older `fin_data*.tsv`
        files into `output_dir/archive/`, and deduplicate identical files by SHA-256.
        If `merge_existing=True` and `base_fin_data` is omitted, automatically merges
        with the latest existing `fin_data*.tsv` found in `output_dir`.
        """
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        date_tag = export_date or date.today().isoformat()
        target_path = out_dir / f"fin_data_{date_tag}.tsv"

        latest_existing = find_latest_fin_data(out_dir)
        if base_fin_data:
            base_path: Optional[Path] = Path(base_fin_data)
        elif merge_existing and latest_existing and latest_existing.resolve() != target_path.resolve():
            base_path = latest_existing
        else:
            base_path = None

        rows, gene_count = self.build_fin_data_rows(scope=scope, base_fin_data=base_path)

        # Write atomically via temp file
        tmp_path = out_dir / f".tmp_fin_data_{date_tag}.tsv"
        with open(tmp_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["Gene", "variable", "value"])
            writer.writerows(rows)

        new_hash = sha256_file(tmp_path)

        # Check if latest existing file in output_dir is already byte-identical
        if latest_existing and latest_existing.resolve() != target_path.resolve():
            if sha256_file(latest_existing) == new_hash:
                # Content is identical to current active file; skip churning archive
                tmp_path.unlink()
                return {
                    "status": "unchanged",
                    "active_file": str(latest_existing),
                    "genes_included": gene_count,
                    "rows_written": len(rows),
                    "sha256": new_hash,
                    "archived_files": [],
                    "duplicates_removed": [],
                }

        shutil.move(str(tmp_path), str(target_path))

        archive_res = {"archived_files": [], "duplicates_removed": [], "archive_dir": str(out_dir / "archive")}
        if archive_existing:
            archive_res = archive_old_fin_data(out_dir, active_file=target_path)

        return {
            "status": "exported",
            "active_file": str(target_path),
            "genes_included": gene_count,
            "rows_written": len(rows),
            "sha256": new_hash,
            **archive_res,
        }
