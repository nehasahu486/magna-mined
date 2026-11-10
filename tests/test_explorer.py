"""
Unit tests for Phase 7: CatalogExplorer querying, gene dossiers, and multi-format exports.
"""

import csv
import json
import openpyxl
import pytest

from magnamined.database import DatabaseManager
from magnamined.explorer import CatalogExplorer
from magnamined.schema import (
    EvidenceConfidence,
    EvidenceSource,
    ExperimentalMethod,
    GeneEvidenceRecord,
    PhenotypeCategory,
)


@pytest.fixture
def populated_db(tmp_path):
    db_file = tmp_path / "exp.db"
    db = DatabaseManager(str(db_file))
    with db.get_connection() as conn:
        conn.execute(
            "INSERT INTO genes (gene, gene_name, size_aa, process_affected, blast2go_annotation, pfam) "
            "VALUES ('MGG_09565', 'PMK1', 356, 'Not reported', 'mitogen-activated protein kinase', 'PF00069')"
        )
        conn.execute(
            "INSERT INTO genes (gene, gene_name, size_aa, process_affected, blast2go_annotation, pfam) "
            "VALUES ('MGG_01001', 'ENG1', 400, 'Not reported', 'endo-glucanase', 'PF00722')"
        )
        conn.commit()

    db.save_evidence(GeneEvidenceRecord(
        gene_id="MGG_09565",
        gene_name="PMK1",
        process_affected=PhenotypeCategory.APPRESSORIUM,
        secondary_processes=[PhenotypeCategory.VIRULENCE],
        detailed_phenotype="Abolished appressorium formation and virulence.",
        experimental_method=ExperimentalMethod.KNOCKOUT_DELETION,
        complementation_tested=True,
        confidence_tier=EvidenceConfidence.TIER_1,
        verbatim_quote="Deletion of MoPMK1 abolished appressorium formation.",
        evidence_source=EvidenceSource.JATS_XML,
        pmid="8946911",
        publication_year=2023,
    ))
    return db


def test_search_genes_filters(populated_db):
    explorer = CatalogExplorer(db=populated_db)

    all_rows = explorer.search_genes()
    assert len(all_rows) == 2

    by_proc = explorer.search_genes(process="appressorium")
    assert [r["gene"] for r in by_proc] == ["MGG_09565"]

    by_tier = explorer.search_genes(confidence="tier_1_high")
    assert [r["gene"] for r in by_tier] == ["MGG_09565"]

    by_src = explorer.search_genes(source="jats_xml")
    assert [r["gene"] for r in by_src] == ["MGG_09565"]

    by_ext = explorer.search_genes(has_extractions=True)
    assert [r["gene"] for r in by_ext] == ["MGG_09565"]

    by_kw = explorer.search_genes(query="endo-glucanase")
    assert [r["gene"] for r in by_kw] == ["MGG_01001"]

    # Newly classified filter excludes genes already classified in MagnaGenes baseline
    with populated_db.get_connection() as conn:
        conn.execute("UPDATE genes SET baseline_process = 'virulence' WHERE gene = 'MGG_09565'")
        conn.commit()
    populated_db.save_evidence(GeneEvidenceRecord(
        gene_id="MGG_01001",
        gene_name="ENG1",
        process_affected=PhenotypeCategory.CONIDIATION,
        detailed_phenotype="Reduced conidiation.",
        experimental_method=ExperimentalMethod.KNOCKOUT_DELETION,
        confidence_tier=EvidenceConfidence.TIER_2,
        verbatim_quote="Deletion of ENG1 reduced conidiation.",
        pmid="99999",
    ))
    newly = explorer.search_genes(newly_classified=True)
    assert [r["gene"] for r in newly] == ["MGG_01001"]


def test_get_gene_dossier_with_alias(populated_db):
    explorer = CatalogExplorer(db=populated_db)
    dossier = explorer.get_gene_dossier("MoPMK1")
    assert dossier is not None
    assert dossier["gene"]["gene"] == "MGG_09565"
    assert "PMK1" in dossier["aliases"]
    assert len(dossier["extractions"]) == 1
    assert dossier["extractions"][0]["pmid"] == "8946911"


def test_export_all_formats(populated_db, tmp_path):
    explorer = CatalogExplorer(db=populated_db)

    # 1. CSV
    csv_out = tmp_path / "out.csv"
    assert explorer.export_catalog(str(csv_out), has_extractions=True) == 1
    with open(csv_out, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert rows[0]["gene"] == "MGG_09565"
    assert rows[0]["knockout_reported"] == "1"
    assert "appressorium formation" in rows[0]["knockout_phenotype"]
    assert "additional_functions" in rows[0]
    assert "additional_function_references" in rows[0]

    # 2. TSV
    tsv_out = tmp_path / "out.tsv"
    assert explorer.export_catalog(str(tsv_out)) == 2
    with open(tsv_out, encoding="utf-8") as f:
        trows = list(csv.DictReader(f, delimiter="\t"))
    assert len(trows) == 2

    # 3. JSON (includes nested extractions)
    json_out = tmp_path / "out.json"
    assert explorer.export_catalog(str(json_out), has_extractions=True) == 1
    data = json.loads(json_out.read_text(encoding="utf-8"))
    assert len(data) == 1
    assert data[0]["extractions"][0]["verbatim_quote"] == "Deletion of MoPMK1 abolished appressorium formation."

    # 4. Excel (.xlsx with Genes + Extractions sheets)
    xlsx_out = tmp_path / "out.xlsx"
    assert explorer.export_catalog(str(xlsx_out), has_extractions=True) == 1
    wb = openpyxl.load_workbook(str(xlsx_out))
    assert wb.sheetnames == ["Genes", "Extractions"]
    assert wb["Genes"].max_row == 2
    assert wb["Extractions"].max_row == 2
    wb.close()

    # 5. Self-contained HTML Explorer
    html_out = tmp_path / "report.html"
    assert explorer.export_catalog(str(html_out), has_extractions=True) == 1
    html_text = html_out.read_text(encoding="utf-8")
    assert "Magna-Mined Gene Curation Explorer" in html_text
    assert "Prior to Search" in html_text
    assert "chartCoverage" in html_text
    assert "chartProcesses" in html_text
    assert "chartTiers" in html_text
    assert "chartYears" in html_text
    assert "MGG_09565" in html_text
    assert "Deletion of MoPMK1 abolished appressorium formation." in html_text
