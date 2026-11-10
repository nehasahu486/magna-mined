"""
Unit tests for Phase 6: Evidence Conflict Resolution & Zenodo Ground-Truth Benchmark.
"""

import tempfile
import pytest

from magnamined.benchmark import ZenodoBenchmark, parse_zenodo_processes
from magnamined.database import DatabaseManager
from magnamined.reconciler import EvidenceReconciler
from magnamined.schema import (
    EvidenceConfidence,
    EvidenceSource,
    ExperimentalMethod,
    GeneEvidenceRecord,
    PhenotypeCategory,
)


@pytest.fixture
def rec_db():
    with tempfile.NamedTemporaryFile(suffix=".db") as tmp:
        db = DatabaseManager(tmp.name)
        with db.get_connection() as conn:
            conn.execute("INSERT INTO genes (gene, gene_name, process_affected) VALUES ('MGG_09565', 'PMK1', 'Not reported')")
            conn.execute("INSERT INTO genes (gene, gene_name, process_affected) VALUES ('MGG_01001', 'ENG1', 'Not reported')")
            conn.commit()
        yield db


def _make_record(
    gene_id="MGG_09565",
    process=PhenotypeCategory.VIRULENCE,
    tier=EvidenceConfidence.TIER_2,
    source=EvidenceSource.JATS_XML,
    method=ExperimentalMethod.KNOCKOUT_DELETION,
    comp=False,
    dispensable=False,
    pmid="100",
    year=2022,
):
    return GeneEvidenceRecord(
        gene_id=gene_id,
        process_affected=process,
        detailed_phenotype=f"Observed {process.value}",
        is_dispensable=dispensable,
        experimental_method=method,
        complementation_tested=comp,
        confidence_tier=tier,
        verbatim_quote="Sample verbatim quote from paper.",
        evidence_source=source,
        pmid=pmid,
        publication_year=year,
    )


def test_save_evidence_respects_tier_precedence(rec_db):
    # Save high-tier record first, then low-tier record
    rec_db.save_evidence(_make_record(process=PhenotypeCategory.VIRULENCE, tier=EvidenceConfidence.TIER_1, comp=True, pmid="1"))
    rec_db.save_evidence(_make_record(process=PhenotypeCategory.CONIDIATION, tier=EvidenceConfidence.TIER_4, source=EvidenceSource.ABSTRACT_ONLY, pmid="2"))

    gene = rec_db.get_gene("MGG_09565")
    assert gene["process_affected"] == "virulence"
    assert gene["consensus_tier"] == "tier_1_high"


def test_reconciler_hierarchy_and_conflict_review(rec_db):
    # MGG_09565: Tier 4 says dispensable, Tier 1 (with complementation) says appressorium formation
    id_low = rec_db.save_evidence(_make_record(
        gene_id="MGG_09565",
        process=PhenotypeCategory.DISPENSABLE,
        dispensable=True,
        tier=EvidenceConfidence.TIER_4,
        source=EvidenceSource.ABSTRACT_ONLY,
        method=ExperimentalMethod.BIOCHEMICAL_ASSAY,
        pmid="111",
        year=2015,
    ))
    id_high = rec_db.save_evidence(_make_record(
        gene_id="MGG_09565",
        process=PhenotypeCategory.APPRESSORIUM,
        tier=EvidenceConfidence.TIER_1,
        source=EvidenceSource.JATS_XML,
        method=ExperimentalMethod.KNOCKOUT_DELETION,
        comp=True,
        pmid="222",
        year=2023,
    ))

    reconciler = EvidenceReconciler(db=rec_db)
    summary = reconciler.reconcile_all()

    assert summary["genes_reconciled"] == 1
    assert summary["conflicts_detected"] == 1

    # Consensus winner must be Tier 1 appressorium formation, flagged as conflict
    gene = rec_db.get_gene("MGG_09565")
    assert gene["process_affected"] == "appressorium formation"
    assert gene["consensus_tier"] == "tier_1_high"
    assert gene["conflict_status"] == "conflict"

    queue = reconciler.get_review_queue()
    assert len(queue) == 1
    assert queue[0]["gene_id"] == "MGG_09565"
    assert len(queue[0]["competing_extractions"]) == 2

    # Human curator reviews and overrides with id_low
    resolved = reconciler.mark_reviewed(id_low, notes="Verified negative result in clean background")
    assert resolved["conflict_status"] == "resolved"
    assert resolved["consensus_process"] == "dispensable / wild-type-like"
    assert reconciler.get_review_queue() == []


def test_parse_zenodo_processes_handles_typos_and_synonyms():
    assert parse_zenodo_processes("virulenece, penetration") == {
        PhenotypeCategory.VIRULENCE,
        PhenotypeCategory.APPRESSORIUM,
    }
    assert parse_zenodo_processes("appresorium formation, condiation, pimentation") == {
        PhenotypeCategory.APPRESSORIUM,
        PhenotypeCategory.CONIDIATION,
        PhenotypeCategory.STRESS_RESPONSE,
    }
    assert parse_zenodo_processes("none") == {PhenotypeCategory.DISPENSABLE}
    assert parse_zenodo_processes("Not reported") == {PhenotypeCategory.NOT_REPORTED}


def test_evaluate_current_db_benchmark(tmp_path):
    csv_file = tmp_path / "zenodo.csv"
    csv_file.write_text(
        "Gene,size (aa),Ensembl link,Uniprot,gene name,process affected,detailed mutant phenotype,"
        "PMID or link to publication,PMID or link to publication.1,BLAST2GO annotation,Pfam\n"
        "MGG_09565,356,e,u,PMK1,\"appressorium formation, virulence\",obs,8946911,,,\n"
        "MGG_01001,400,e,u,ENG1,conidiation,obs,1234567,,,\n"
        "MGG_09999,200,e,u,NEW1,Not reported,,,,,\n",
        encoding="utf-8",
    )

    db = DatabaseManager(str(tmp_path / "bench.db"))
    db.seed_from_csv(str(csv_file))

    # Exact match for MGG_01001, partial match for MGG_09565 (predicts virulence only), novel discovery for MGG_09999
    db.save_evidence(_make_record(gene_id="MGG_09565", process=PhenotypeCategory.VIRULENCE, tier=EvidenceConfidence.TIER_1, comp=True))
    db.save_evidence(_make_record(gene_id="MGG_01001", process=PhenotypeCategory.CONIDIATION, tier=EvidenceConfidence.TIER_2))
    db.save_evidence(_make_record(gene_id="MGG_09999", process=PhenotypeCategory.STRESS_RESPONSE, tier=EvidenceConfidence.TIER_2))

    bench = ZenodoBenchmark(db=db, excel_path=str(tmp_path / "missing.xlsx"), csv_path=str(csv_file))
    metrics = bench.evaluate_current_db()

    # Overlap = 2 genes (MGG_09565, MGG_01001); Gold = 3 labels, Pred = 2 labels -> TP=2, FP=0, FN=1
    assert metrics["evaluated_genes"] == 2
    assert metrics["micro_precision"] == 1.0
    assert metrics["micro_recall"] == round(2 / 3, 4)
    assert metrics["any_overlap_rate"] == 1.0
    assert metrics["exact_match_rate"] == 0.5
    assert metrics["novel_discoveries_count"] == 1
    assert metrics["novel_discoveries"][0]["gene_id"] == "MGG_09999"
