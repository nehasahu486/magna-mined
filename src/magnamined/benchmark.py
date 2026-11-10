"""
Zenodo Ground-Truth Validation Benchmark for Magna-Mined (Phase 6).
Compares automated literature extractions against the Zenodo Record 4647766 baseline
(12,522 genes; ~1,274 manually curated phenotypes) to compute Precision, Recall, F1,
novel discovery counts, and discrepancy audits.
"""

import csv
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from .database import DatabaseManager
from .reconciler import EvidenceReconciler
from .schema import PhenotypeCategory

# Faithful mapping from Zenodo Magnagenes.xlsx vocabulary (including original typos)
ZENODO_PROCESS_MAP: Dict[str, PhenotypeCategory] = {
    "virulence": PhenotypeCategory.VIRULENCE,
    "virulenece": PhenotypeCategory.VIRULENCE,
    "effector": PhenotypeCategory.VIRULENCE,
    "elicitor": PhenotypeCategory.VIRULENCE,
    "appressorium formation": PhenotypeCategory.APPRESSORIUM,
    "appresorium formation": PhenotypeCategory.APPRESSORIUM,
    "penetration": PhenotypeCategory.APPRESSORIUM,
    "turgor": PhenotypeCategory.APPRESSORIUM,
    "adhesion": PhenotypeCategory.APPRESSORIUM,
    "conidiation": PhenotypeCategory.CONIDIATION,
    "condiation": PhenotypeCategory.CONIDIATION,
    "stress response": PhenotypeCategory.STRESS_RESPONSE,
    "pigmentation": PhenotypeCategory.STRESS_RESPONSE,
    "pimentation": PhenotypeCategory.STRESS_RESPONSE,
    "growth": PhenotypeCategory.GROWTH_DEVELOPMENT,
    "growth/development": PhenotypeCategory.GROWTH_DEVELOPMENT,
    "germination": PhenotypeCategory.GROWTH_DEVELOPMENT,
    "sexual development": PhenotypeCategory.GROWTH_DEVELOPMENT,
    "putative essential gene": PhenotypeCategory.GROWTH_DEVELOPMENT,
    "none": PhenotypeCategory.DISPENSABLE,
    "dispensable / wild-type-like": PhenotypeCategory.DISPENSABLE,
    "not reported": PhenotypeCategory.NOT_REPORTED,
    "mixed outcomes": PhenotypeCategory.NOT_REPORTED,
}


def parse_zenodo_processes(raw: Optional[str]) -> Set[PhenotypeCategory]:
    """Map comma-separated Zenodo process labels to controlled PhenotypeCategory set."""
    if not raw:
        return {PhenotypeCategory.NOT_REPORTED}
    cats: Set[PhenotypeCategory] = set()
    for tok in str(raw).split(","):
        norm = tok.strip().lower()
        if not norm:
            continue
        mapped = ZENODO_PROCESS_MAP.get(norm)
        if mapped and mapped != PhenotypeCategory.NOT_REPORTED:
            cats.add(mapped)
    return cats or {PhenotypeCategory.NOT_REPORTED}


def extract_pmid(raw_link: Optional[str]) -> Optional[str]:
    """Extract numeric PMID from a PubMed URL or raw PMID string."""
    if not raw_link:
        return None
    m = re.search(r"(\d{6,9})", str(raw_link))
    return m.group(1) if m else None


def load_zenodo_ground_truth(
    excel_path: str = "data/Magnagenes.xlsx",
    csv_path: str = "data/MagnaGenes.csv",
) -> Dict[str, Dict[str, Any]]:
    """Load baseline gene labels and PMIDs directly from the Zenodo dataset."""
    truth: Dict[str, Dict[str, Any]] = {}
    ep = Path(excel_path)
    cp = Path(csv_path)

    if ep.exists():
        import openpyxl
        wb = openpyxl.load_workbook(str(ep), read_only=True)
        ws = wb[wb.sheetnames[0]]
        rows = ws.iter_rows(values_only=True)
        next(rows, None)
        for row in rows:
            if not row or not row[0]:
                continue
            gid = str(row[0]).strip().upper()
            if not gid.startswith("MGG_"):
                continue
            raw_proc = str(row[5]).strip() if len(row) > 5 and row[5] else "Not reported"
            truth[gid] = {
                "gene_id": gid,
                "gene_name": str(row[4]).strip() if len(row) > 4 and row[4] and str(row[4]).strip() else None,
                "raw_process": raw_proc,
                "processes": parse_zenodo_processes(raw_proc),
                "pmid": extract_pmid(row[7] if len(row) > 7 else None),
            }
        wb.close()
        return truth

    if cp.exists():
        with open(cp, mode="r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                gid = (row.get("Gene") or row.get("gene") or "").strip().upper()
                if not gid.startswith("MGG_"):
                    continue
                raw_proc = (row.get("process affected") or "Not reported").strip()
                gname = (row.get("gene name") or "").strip() or None
                truth[gid] = {
                    "gene_id": gid,
                    "gene_name": gname,
                    "raw_process": raw_proc,
                    "processes": parse_zenodo_processes(raw_proc),
                    "pmid": extract_pmid(row.get("PMID or link to publication")),
                }
        return truth

    raise FileNotFoundError(f"Neither {excel_path} nor {csv_path} found.")


def _compute_multilabel_metrics(
    gold_map: Dict[str, Set[PhenotypeCategory]],
    pred_map: Dict[str, Set[PhenotypeCategory]],
) -> Dict[str, Any]:
    """Compute micro Precision/Recall/F1, exact match, any-overlap match, and per-category scores."""
    eval_cats = [c for c in PhenotypeCategory if c != PhenotypeCategory.NOT_REPORTED]
    per_cat_counts = {c: {"tp": 0, "fp": 0, "fn": 0} for c in eval_cats}

    exact_matches = 0
    any_overlap_matches = 0
    discrepancies: List[Dict[str, Any]] = []

    overlap_genes = sorted(set(gold_map.keys()) & set(pred_map.keys()))
    for gid in overlap_genes:
        g_set = {c for c in gold_map[gid] if c != PhenotypeCategory.NOT_REPORTED}
        p_set = {c for c in pred_map[gid] if c != PhenotypeCategory.NOT_REPORTED}

        if g_set == p_set:
            exact_matches += 1
        if g_set & p_set:
            any_overlap_matches += 1
        else:
            discrepancies.append({
                "gene_id": gid,
                "zenodo_expected": sorted(c.value for c in g_set),
                "magna_mined_predicted": sorted(c.value for c in p_set),
            })

        for c in eval_cats:
            in_g = c in g_set
            in_p = c in p_set
            if in_g and in_p:
                per_cat_counts[c]["tp"] += 1
            elif in_p and not in_g:
                per_cat_counts[c]["fp"] += 1
            elif in_g and not in_p:
                per_cat_counts[c]["fn"] += 1

    total_tp = sum(v["tp"] for v in per_cat_counts.values())
    total_fp = sum(v["fp"] for v in per_cat_counts.values())
    total_fn = sum(v["fn"] for v in per_cat_counts.values())

    def _prf(tp: int, fp: int, fn: int) -> Dict[str, float]:
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) > 0 else 0.0
        return {"precision": round(prec, 4), "recall": round(rec, 4), "f1": round(f1, 4)}

    n = len(overlap_genes)
    per_category = {}
    for c, counts in per_cat_counts.items():
        if counts["tp"] + counts["fp"] + counts["fn"] > 0:
            per_category[c.value] = {**_prf(**counts), **counts}

    micro = _prf(total_tp, total_fp, total_fn)
    return {
        "evaluated_genes": n,
        "exact_match_rate": round(exact_matches / n, 4) if n else 0.0,
        "any_overlap_rate": round(any_overlap_matches / n, 4) if n else 0.0,
        "micro_precision": micro["precision"],
        "micro_recall": micro["recall"],
        "micro_f1": micro["f1"],
        "tp": total_tp,
        "fp": total_fp,
        "fn": total_fn,
        "per_category": per_category,
        "discrepancies": discrepancies,
    }


class ZenodoBenchmark:
    """Evaluates Magna-Mined curation outputs against Zenodo Record 4647766."""

    def __init__(
        self,
        db: Optional[DatabaseManager] = None,
        db_path: str = "magnaporthe.db",
        excel_path: str = "data/Magnagenes.xlsx",
        csv_path: str = "data/MagnaGenes.csv",
    ):
        self.db = db or DatabaseManager(db_path)
        self.ground_truth = load_zenodo_ground_truth(excel_path=excel_path, csv_path=csv_path)

    def evaluate_current_db(self) -> Dict[str, Any]:
        """
        Evaluate all automated `extractions` currently stored in SQLite against
        the Zenodo baseline ground truth.
        """
        reconciler = EvidenceReconciler(db=self.db)
        with self.db.get_connection() as conn:
            rows = conn.execute(
                """
                SELECT g.gene, e.*
                FROM extractions e
                JOIN genes g ON e.gene_id = g.id
                ORDER BY g.gene, e.id
                """
            ).fetchall()

        by_gene: Dict[str, List[Dict[str, Any]]] = {}
        for r in rows:
            by_gene.setdefault(r["gene"].upper(), []).append(dict(r))

        gold_overlap: Dict[str, Set[PhenotypeCategory]] = {}
        pred_overlap: Dict[str, Set[PhenotypeCategory]] = {}
        novel_genes: List[Dict[str, Any]] = []

        for gid, ext_rows in by_gene.items():
            ranked = sorted(ext_rows, key=reconciler.rank_extraction, reverse=True)
            meaningful = [
                r for r in ranked
                if any(p != PhenotypeCategory.NOT_REPORTED.value for p in reconciler._parse_processes(r.get("process_affected")))
            ]
            if not meaningful:
                continue

            winner = meaningful[0]
            pred_cats = parse_zenodo_processes(winner.get("process_affected"))
            # Combine same-tier corroborations
            for other in meaningful[1:]:
                if other.get("confidence_tier") == winner.get("confidence_tier"):
                    other_cats = parse_zenodo_processes(other.get("process_affected"))
                    if PhenotypeCategory.DISPENSABLE not in pred_cats and PhenotypeCategory.DISPENSABLE not in other_cats:
                        pred_cats |= other_cats
            pred_cats.discard(PhenotypeCategory.NOT_REPORTED)
            if not pred_cats:
                continue

            gt = self.ground_truth.get(gid)
            if not gt:
                continue

            gt_cats = gt["processes"]
            if gt_cats == {PhenotypeCategory.NOT_REPORTED}:
                novel_genes.append({
                    "gene_id": gid,
                    "gene_name": gt.get("gene_name"),
                    "predicted": sorted(c.value for c in pred_cats),
                    "confidence_tier": winner.get("confidence_tier"),
                    "pmid": winner.get("pmid"),
                })
            else:
                gold_overlap[gid] = gt_cats
                pred_overlap[gid] = pred_cats

        metrics = _compute_multilabel_metrics(gold_overlap, pred_overlap)
        metrics["genes_with_extractions"] = len(by_gene)
        metrics["novel_discoveries_count"] = len(novel_genes)
        metrics["novel_discoveries"] = novel_genes[:15]
        return metrics

    def run_live_zenodo_benchmark(
        self,
        gene_ids: Optional[List[str]] = None,
        sample_size: int = 8,
        api_key: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Blind held-out benchmark:
        Seeds an isolated temp DB with curated Zenodo genes (with phenotypes blanked out),
        fetches their cited papers from Europe PMC, runs extraction, and scores
        Precision/Recall/F1 + quote grounding against the hidden Zenodo labels.
        """
        from .paper_miner import PaperMiner

        default_canonical = [
            "MGG_09565",  # PMK1
            "MGG_09898",  # MAC1
            "MGG_12958",  # MST12
            "MGG_05871",  # PTH11
            "MGG_06368",  # CPKA
            "MGG_02252",  # BUF1
            "MGG_07219",  # ALB1
            "MGG_03977",  # COS1
            "MGG_00184",  # HTF1
            "MGG_05059",  # RSY1
        ]
        targets = gene_ids or [g for g in default_canonical if g in self.ground_truth][:sample_size]

        with tempfile.NamedTemporaryFile(suffix=".db") as tmp:
            tmp_db = DatabaseManager(tmp.name)
            with tmp_db.get_connection() as conn:
                for gid in targets:
                    gt = self.ground_truth.get(gid)
                    if not gt:
                        continue
                    conn.execute(
                        "INSERT INTO genes (gene, gene_name, process_affected) VALUES (?, ?, 'Not reported')",
                        (gid, gt.get("gene_name")),
                    )
                conn.commit()

            engine = PaperMiner(db_path=tmp.name, api_key=api_key)
            gold_map: Dict[str, Set[PhenotypeCategory]] = {}
            pred_map: Dict[str, Set[PhenotypeCategory]] = {}
            quotes_total = 0
            quotes_verified = 0

            for gid in targets:
                gt = self.ground_truth.get(gid)
                if not gt:
                    continue
                gold_map[gid] = gt["processes"]
                pmid = gt.get("pmid")
                articles = []
                if pmid:
                    articles = engine.lit_client.search_literature(f"EXT_ID:{pmid} AND SRC:MED", max_results=1)
                if not articles:
                    articles = engine.lit_client.search_gene(gid, gene_name=gt.get("gene_name"), max_results=2)

                for art in articles:
                    res = engine.process_paper(art, verbose=False, force=True)
                    for rec in res.get("extractions", []):
                        if rec.gene_id.upper() == gid:
                            quotes_total += 1
                            if rec.verbatim_quote:
                                quotes_verified += 1
                            if rec.process_affected != PhenotypeCategory.NOT_REPORTED:
                                pset = pred_map.setdefault(gid, set())
                                pset.add(rec.process_affected)
                                for sp in rec.secondary_processes:
                                    if sp != PhenotypeCategory.NOT_REPORTED:
                                        pset.add(sp)

            metrics = _compute_multilabel_metrics(gold_map, pred_map)
            metrics["target_genes"] = len(targets)
            metrics["genes_with_predictions"] = len(pred_map)
            metrics["retrieval_coverage"] = round(len(pred_map) / len(targets), 4) if targets else 0.0
            metrics["quotes_extracted"] = quotes_total
            metrics["quote_grounding_rate"] = round(quotes_verified / quotes_total, 4) if quotes_total else 0.0
            return metrics
