"""
Evidence Conflict Resolution & Audit Ledger for Magna-Mined (Phase 6).
Ranks competing literature extractions by experimental rigor, source provenance,
and recency; detects polarity/phenotype conflicts; updates consensus gene records;
and maintains a human curator review queue.
"""

from typing import Any, Dict, List, Optional, Set, Tuple

from .database import DatabaseManager
from .schema import PhenotypeCategory

TIER_RANK: Dict[str, int] = {
    "tier_1_high": 4,
    "tier_2_high": 3,
    "tier_3_medium": 2,
    "tier_4_low": 1,
}

SOURCE_RANK: Dict[str, int] = {
    "jats_xml": 3,
    "pdf": 2,
    "abstract_only": 1,
}

METHOD_RANK: Dict[str, int] = {
    "CRISPR/Cas9": 3,
    "knockout/deletion": 3,
    "complementation": 3,
    "overexpression": 2,
    "RNAi silencing": 2,
    "biochemical assay": 1,
    "computational only": 1,
}


class EvidenceReconciler:
    """Deterministic conflict detector and consensus resolver over SQLite extractions."""

    def __init__(self, db: Optional[DatabaseManager] = None, db_path: str = "magnaporthe.db"):
        self.db = db or DatabaseManager(db_path)

    @staticmethod
    def rank_extraction(rec: Dict[str, Any]) -> Tuple[int, int, int, int, int, int]:
        """
        Compute deterministic precedence tuple (higher wins):
        (human_reviewed, confidence_tier, evidence_source, method, complementation, pub_year)
        """
        return (
            1 if rec.get("reviewed_by_human") else 0,
            TIER_RANK.get(rec.get("confidence_tier") or "", 0),
            SOURCE_RANK.get(rec.get("evidence_source") or "jats_xml", 1),
            METHOD_RANK.get(rec.get("experimental_method") or "", 0),
            1 if rec.get("complementation_tested") else 0,
            int(rec.get("publication_year") or 0),
        )

    @staticmethod
    def _parse_processes(proc_str: Optional[str]) -> List[str]:
        if not proc_str:
            return []
        return [p.strip().lower() for p in proc_str.split(",") if p.strip()]

    @classmethod
    def detect_conflicts(cls, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Inspect extraction records for a single gene and detect biological conflicts:
        1. Polarity conflict: one paper claims dispensable/wild-type-like while another
           claims a mutant defect (virulence, appressorium, conidiation, etc.).
        2. High-tier disjoint conflict: >=2 Tier 1/2 records assert completely disjoint
           non-dispensable phenotypes across distinct papers.
        """
        if len(records) < 2:
            return {"has_conflict": False, "conflict_type": "none", "details": ""}

        if any(r.get("reviewed_by_human") for r in records):
            return {"has_conflict": False, "conflict_type": "resolved", "details": "Resolved by human curator."}

        has_dispensable = False
        has_defect = False
        high_tier_sets: List[Tuple[str, Set[str]]] = []

        for r in records:
            procs = [p for p in cls._parse_processes(r.get("process_affected")) if p != PhenotypeCategory.NOT_REPORTED.value]
            if not procs:
                continue
            is_disp = any(p == PhenotypeCategory.DISPENSABLE.value or p == "none" for p in procs)
            if is_disp:
                has_dispensable = True
            else:
                has_defect = True
                if TIER_RANK.get(r.get("confidence_tier") or "", 0) >= 3:
                    high_tier_sets.append((str(r.get("pmid") or r.get("id")), set(procs)))

        if has_dispensable and has_defect:
            return {
                "has_conflict": True,
                "conflict_type": "polarity_conflict",
                "details": "Dispensable/wild-type result conflicts with reported mutant phenotype defect.",
            }

        if len(high_tier_sets) >= 2:
            # Check if every high-tier record shares at least one process with the top-ranked set
            anchor_procs = high_tier_sets[0][1]
            disjoint = [pmid for pmid, pset in high_tier_sets[1:] if not (anchor_procs & pset)]
            if disjoint:
                return {
                    "has_conflict": True,
                    "conflict_type": "phenotype_discrepancy",
                    "details": f"High-confidence papers report disjoint phenotypes ({sorted(anchor_procs)} vs other PMIDs {disjoint}).",
                }

        return {"has_conflict": False, "conflict_type": "none", "details": ""}

    def reconcile_gene(self, gene_id: str) -> Dict[str, Any]:
        """Reconcile all extractions for one gene and update the master `genes` row."""
        with self.db.get_connection() as conn:
            gene_row = conn.execute(
                "SELECT id, gene, gene_name, process_affected FROM genes WHERE UPPER(gene) = UPPER(?)",
                (gene_id.strip(),)
            ).fetchone()
            if not gene_row:
                raise ValueError(f"Gene {gene_id} not found in database.")

            rows = conn.execute(
                "SELECT * FROM extractions WHERE gene_id = ? ORDER BY id ASC",
                (gene_row["id"],)
            ).fetchall()
            records = [dict(r) for r in rows]

            if not records:
                return {
                    "gene_id": gene_row["gene"],
                    "status": "no_extractions",
                    "consensus_process": gene_row["process_affected"],
                    "consensus_tier": None,
                    "conflict_status": "none",
                }

            # Rank records (highest precedence first)
            ranked = sorted(records, key=self.rank_extraction, reverse=True)
            conflict_info = self.detect_conflicts(ranked)

            # Filter meaningful records (ignore pure 'not reported' if anything stronger exists)
            meaningful = [
                r for r in ranked
                if any(p != PhenotypeCategory.NOT_REPORTED.value for p in self._parse_processes(r.get("process_affected")))
            ]
            winner = meaningful[0] if meaningful else ranked[0]

            # Build consensus processes: winner's processes + non-conflicting same-tier corroborations
            winner_procs = self._parse_processes(winner.get("process_affected"))
            consensus_procs = list(winner_procs)
            winner_tier = winner.get("confidence_tier")
            winner_is_disp = PhenotypeCategory.DISPENSABLE.value in winner_procs

            if not winner.get("reviewed_by_human") and not winner_is_disp:
                for other in meaningful[1:]:
                    if other.get("confidence_tier") == winner_tier:
                        other_procs = self._parse_processes(other.get("process_affected"))
                        if PhenotypeCategory.DISPENSABLE.value not in other_procs:
                            for p in other_procs:
                                if p not in consensus_procs and p != PhenotypeCategory.NOT_REPORTED.value:
                                    consensus_procs.append(p)

            consensus_str = ", ".join(consensus_procs) if consensus_procs else PhenotypeCategory.NOT_REPORTED.value

            if any(r.get("reviewed_by_human") for r in records):
                conflict_status = "resolved"
            elif conflict_info["has_conflict"]:
                conflict_status = "conflict"
            else:
                conflict_status = "none"

            conn.execute(
                """
                UPDATE genes SET
                    process_affected = ?,
                    consensus_tier = ?,
                    conflict_status = ?,
                    pmid_or_link_to_publication = COALESCE(pmid_or_link_to_publication, ?)
                WHERE id = ?
                """,
                (
                    consensus_str,
                    winner_tier,
                    conflict_status,
                    winner.get("pmid") or winner.get("doi"),
                    gene_row["id"],
                ),
            )
            conn.commit()

            return {
                "gene_id": gene_row["gene"],
                "gene_name": gene_row["gene_name"],
                "status": "reconciled",
                "extraction_count": len(records),
                "winning_extraction_id": winner["id"],
                "winning_pmid": winner.get("pmid"),
                "consensus_process": consensus_str,
                "consensus_tier": winner_tier,
                "conflict_status": conflict_status,
                "conflict_details": conflict_info["details"],
            }

    def reconcile_all(self) -> Dict[str, Any]:
        """Reconcile all genes that have automated extraction records."""
        with self.db.get_connection() as conn:
            gene_rows = conn.execute(
                """
                SELECT DISTINCT g.gene
                FROM genes g
                JOIN extractions e ON e.gene_id = g.id
                ORDER BY g.gene
                """
            ).fetchall()

        results = [self.reconcile_gene(r["gene"]) for r in gene_rows]
        tier_counts: Dict[str, int] = {}
        conflicts: List[Dict[str, Any]] = []
        resolved = 0

        for res in results:
            t = res.get("consensus_tier") or "unknown"
            tier_counts[t] = tier_counts.get(t, 0) + 1
            if res["conflict_status"] == "conflict":
                conflicts.append(res)
            elif res["conflict_status"] == "resolved":
                resolved += 1

        return {
            "genes_reconciled": len(results),
            "conflicts_detected": len(conflicts),
            "resolved_by_human": resolved,
            "tier_breakdown": tier_counts,
            "conflicts": conflicts,
        }

    def get_review_queue(self) -> List[Dict[str, Any]]:
        """Return all genes flagged for human review with their competing evidence records."""
        with self.db.get_connection() as conn:
            flagged = conn.execute(
                "SELECT id, gene, gene_name, process_affected, consensus_tier FROM genes WHERE conflict_status = 'conflict'"
            ).fetchall()

            queue: List[Dict[str, Any]] = []
            for g in flagged:
                exts = conn.execute(
                    """
                    SELECT id, pmid, process_affected, experimental_method, confidence_tier,
                           evidence_source, publication_year, verbatim_quote, reviewed_by_human
                    FROM extractions
                    WHERE gene_id = ?
                    ORDER BY id ASC
                    """,
                    (g["id"],)
                ).fetchall()
                queue.append({
                    "gene_id": g["gene"],
                    "gene_name": g["gene_name"],
                    "current_consensus": g["process_affected"],
                    "consensus_tier": g["consensus_tier"],
                    "competing_extractions": [dict(e) for e in exts],
                })
            return queue

    def mark_reviewed(self, extraction_id: int, notes: str = "") -> Dict[str, Any]:
        """Mark a specific extraction as human-verified and re-reconcile its gene."""
        with self.db.get_connection() as conn:
            row = conn.execute(
                """
                SELECT e.id, g.gene
                FROM extractions e
                JOIN genes g ON e.gene_id = g.id
                WHERE e.id = ?
                """,
                (extraction_id,)
            ).fetchone()
            if not row:
                raise ValueError(f"Extraction {extraction_id} not found.")

            conn.execute(
                "UPDATE extractions SET reviewed_by_human = 1, review_notes = ? WHERE id = ?",
                (notes, extraction_id)
            )
            conn.commit()
            target_gene = row["gene"]

        return self.reconcile_gene(target_gene)
