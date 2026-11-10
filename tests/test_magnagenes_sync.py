"""
Unit tests for MagnaGenes v2.0 integration: dated fin_data_YYYY-MM-DD.tsv export,
base fin_data merging, and SHA-256 deduplicated archiving.
"""

import csv
import shutil
from pathlib import Path

import pytest

from magnamined.database import DatabaseManager
from magnamined.magnagenes_sync import (
    MagnaGenesSyncExporter,
    archive_old_fin_data,
    find_latest_fin_data,
    sha256_file,
)
from magnamined.schema import (
    EvidenceConfidence,
    EvidenceSource,
    ExperimentalMethod,
    GeneEvidenceRecord,
    PhenotypeCategory,
)


@pytest.fixture
def sync_db(tmp_path):
    db_file = tmp_path / "sync.db"
    db = DatabaseManager(str(db_file))
    with db.get_connection() as conn:
        # Unclassified gene in baseline
        conn.execute(
            "INSERT INTO genes (gene, gene_name, process_affected, baseline_process) "
            "VALUES ('MGG_01311', NULL, 'Not reported', 'Not reported')"
        )
        # Already classified gene in baseline
        conn.execute(
            "INSERT INTO genes (gene, gene_name, process_affected, baseline_process) "
            "VALUES ('MGG_05059', 'RSY1', 'virulence', 'virulence')"
        )
        conn.commit()

    # Newly classify MGG_01311
    db.save_evidence(GeneEvidenceRecord(
        gene_id="MGG_01311",
        gene_name="MoPah1",
        process_affected=PhenotypeCategory.VIRULENCE,
        detailed_phenotype="Deletion of MoPah1 attenuated virulence.",
        experimental_method=ExperimentalMethod.KNOCKOUT_DELETION,
        complementation_tested=True,
        confidence_tier=EvidenceConfidence.TIER_1,
        verbatim_quote="Deletion of MoPah1 attenuated virulence on rice.",
        pmid="35191164",
        publication_year=2022,
    ))

    # Double-check MGG_05059 -> adds additional function conidiation while preserving virulence
    db.save_evidence(GeneEvidenceRecord(
        gene_id="MGG_05059",
        gene_name="RSY1",
        process_affected=PhenotypeCategory.CONIDIATION,
        detailed_phenotype="RSY1 deletion reduced conidiation.",
        experimental_method=ExperimentalMethod.KNOCKOUT_DELETION,
        complementation_tested=True,
        confidence_tier=EvidenceConfidence.TIER_1,
        verbatim_quote="Deletion of RSY1 reduced conidiation.",
        pmid="35295315",
        publication_year=2022,
    ), preserve_baseline_primary=True)

    return db


def test_export_dated_fin_data_and_merge_with_base(sync_db, tmp_path):
    out_dir = tmp_path / "annotation"
    out_dir.mkdir()

    # Seed an initial base fin_data.tsv like magnagenes_v2 has
    base_file = out_dir / "fin_data.tsv"
    base_file.write_text(
        "Gene\tvariable\tvalue\n"
        "MGG_01311\tProcess.Affected\tNot reported\n"
        "MGG_01311\tInterPro\tIPR000001\n"
        "MGG_05059\tProcess.Affected\tvirulence\n"
        "MGG_05059\tKEGG.Pathway\tmap00900\n",
        encoding="utf-8",
    )

    exporter = MagnaGenesSyncExporter(db=sync_db)
    res = exporter.export_dated_fin_data(
        output_dir=str(out_dir),
        export_date="2026-11-10",
        merge_existing=True,
    )

    assert res["status"] == "exported"
    active = Path(res["active_file"])
    assert active.name == "fin_data_2026-11-10.tsv"
    assert active.exists()

    # Old fin_data.tsv was archived cleanly
    assert not base_file.exists()
    assert (out_dir / "archive" / "fin_data.tsv").exists()

    # Read exported rows and verify magnagenes_v2 compatibility
    with open(active, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    pairs = {(r["Gene"], r["variable"]): r["value"] for r in rows}
    all_tuples = {(r["Gene"], r["variable"], r["value"]) for r in rows}

    # Preserved InterPro / KEGG from base
    assert ("MGG_01311", "InterPro", "IPR000001") in all_tuples
    assert ("MGG_05059", "KEGG.Pathway", "map00900") in all_tuples

    # Replaced "Not reported" for newly classified MGG_01311
    assert ("MGG_01311", "Process.Affected", "Not reported") not in all_tuples
    assert ("MGG_01311", "Process.Affected", "virulence") in all_tuples
    assert ("MGG_01311", "PMID.or.link", "https://pubmed.ncbi.nlm.nih.gov/35191164/") in all_tuples
    assert ("MGG_01311", "Knockout.Reported", "Yes") in all_tuples
    assert "Knockout.Phenotype" in {r["variable"] for r in rows if r["Gene"] == "MGG_01311"}

    # Preserved virulence on MGG_05059 and added Additional.Functions + Other.links.to.publication
    assert ("MGG_05059", "Process.Affected", "virulence") in all_tuples
    assert pairs[("MGG_05059", "Additional.Functions")] == "conidiation"
    assert "PMID:35295315" in pairs[("MGG_05059", "Other.links.to.publication")]


def test_find_latest_and_deduplicated_archive(sync_db, tmp_path):
    out_dir = tmp_path / "ann"
    exporter = MagnaGenesSyncExporter(db=sync_db)

    # 1. Export on 2026-08-01
    r1 = exporter.export_dated_fin_data(output_dir=str(out_dir), export_date="2026-08-01")
    assert Path(r1["active_file"]).name == "fin_data_2026-08-01.tsv"

    # 2. Export on 2026-09-01 with identical content -> unchanged, no archive churn
    r2 = exporter.export_dated_fin_data(output_dir=str(out_dir), export_date="2026-09-01")
    assert r2["status"] == "unchanged"
    assert not (out_dir / "archive").exists() or list((out_dir / "archive").iterdir()) == []

    # 3. Add a new gene curation so content changes, then export on 2026-11-10
    with sync_db.get_connection() as conn:
        conn.execute(
            "INSERT INTO genes (gene, gene_name, process_affected, baseline_process) "
            "VALUES ('MGG_09999', 'NEW1', 'Not reported', 'Not reported')"
        )
        conn.commit()
    sync_db.save_evidence(GeneEvidenceRecord(
        gene_id="MGG_09999",
        process_affected=PhenotypeCategory.APPRESSORIUM,
        detailed_phenotype="Appressorium defect.",
        experimental_method=ExperimentalMethod.KNOCKOUT_DELETION,
        confidence_tier=EvidenceConfidence.TIER_2,
        verbatim_quote="Deletion of NEW1 caused appressorium defect.",
        pmid="445566",
    ))

    r3 = exporter.export_dated_fin_data(output_dir=str(out_dir), export_date="2026-11-10")
    assert r3["status"] == "exported"
    assert r3["archived_files"] == ["fin_data_2026-08-01.tsv"]
    assert find_latest_fin_data(out_dir).name == "fin_data_2026-11-10.tsv"

    # 4. Drop a byte-identical copy of the 2026-08-01 file under another date name -> deduplicated!
    dup = out_dir / "fin_data_2026-08-02.tsv"
    shutil.copy(out_dir / "archive" / "fin_data_2026-08-01.tsv", dup)
    arch_res = archive_old_fin_data(out_dir, active_file=Path(r3["active_file"]))
    assert arch_res["duplicates_removed"] == ["fin_data_2026-08-02.tsv"]
    assert arch_res["archived_files"] == []

    # Verify archive/ contains strictly unique SHA-256 hashes
    arch_files = list((out_dir / "archive").iterdir())
    hashes = [sha256_file(p) for p in arch_files]
    assert len(hashes) == len(set(hashes)) == 1
