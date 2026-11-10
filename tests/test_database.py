"""
Unit tests for DatabaseManager operations.
"""

import tempfile
from pathlib import Path
import pytest

from magnamined.database import DatabaseManager
from magnamined.schema import (
    GeneEvidenceRecord,
    PhenotypeCategory,
    ExperimentalMethod,
    EvidenceConfidence,
)


@pytest.fixture
def temp_db():
    with tempfile.NamedTemporaryFile(suffix=".db") as tmp:
        db = DatabaseManager(tmp.name)
        # Insert a sample gene for testing
        with db.get_connection() as conn:
            conn.execute("""
                INSERT INTO genes (gene, gene_name, process_affected)
                VALUES ('MGG_01001', 'ENG1', 'not reported')
            """)
            conn.commit()
        yield db


def test_database_initialization(temp_db):
    stats = temp_db.get_stats()
    assert stats["total_genes"] == 1
    assert stats["genes_with_phenotype"] == 0


def test_get_gene(temp_db):
    by_id = temp_db.get_gene("MGG_01001")
    assert by_id is not None
    assert by_id["gene"] == "MGG_01001"

    by_name = temp_db.get_gene("ENG1")
    assert by_name is not None
    assert by_name["gene"] == "MGG_01001"


def test_save_evidence(temp_db):
    evidence = GeneEvidenceRecord(
        gene_id="MGG_01001",
        gene_name="ENG1",
        process_affected=PhenotypeCategory.VIRULENCE,
        detailed_phenotype="Significantly attenuated pathogenicity",
        experimental_method=ExperimentalMethod.KNOCKOUT_DELETION,
        complementation_tested=True,
        confidence_tier=EvidenceConfidence.TIER_1,
        verbatim_quote="The eng1 deletion mutant failed to cause full disease lesions.",
        pmid="98765432"
    )

    extraction_id = temp_db.save_evidence(evidence)
    assert extraction_id > 0

    stats = temp_db.get_stats()
    assert stats["genes_with_phenotype"] == 1
    assert stats["total_evidence_extractions"] == 1

    updated_gene = temp_db.get_gene("MGG_01001")
    assert updated_gene["process_affected"] == "virulence"


def test_seed_from_csv_nulls_empty_fields(tmp_path):
    csv_file = tmp_path / "mini.csv"
    csv_file.write_text(
        "Gene,size (aa),Ensembl link,Uniprot,gene name,process affected,detailed mutant phenotype,"
        "PMID or link to publication,PMID or link to publication.1,BLAST2GO annotation,Pfam\n"
        "MGG_00001,100,https://e,https://u,,Not reported,,,,hypothetical,\n"
        "MGG_00002,200,https://e,https://u,PMK1,virulence,lost virulence,12345,,kinase,\n",
        encoding="utf-8",
    )
    db = DatabaseManager(str(tmp_path / "csv.db"))
    assert db.seed_from_csv(str(csv_file)) == 2
    stats = db.get_stats()
    assert stats["total_genes"] == 2
    assert stats["genes_with_phenotype"] == 1
    assert stats["genes_with_literature"] == 1
