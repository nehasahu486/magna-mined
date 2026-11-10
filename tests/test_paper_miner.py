"""
Unit tests for PaperMiner and paper-first bulk literature mining.
"""

from unittest.mock import Mock, patch
import pytest
import tempfile
from pathlib import Path

from magnamined.database import DatabaseManager
from magnamined.paper_miner import PaperMiner


@pytest.fixture
def temp_engine():
    with tempfile.NamedTemporaryFile(suffix=".db") as tmp:
        # Seed test db with test genes
        db = DatabaseManager(tmp.name)
        with db.get_connection() as conn:
            conn.execute("INSERT INTO genes (gene, gene_name, process_affected) VALUES ('MGG_09565', 'PMK1', 'not reported')")
            conn.execute("INSERT INTO genes (gene, gene_name, process_affected) VALUES ('MGG_01001', 'ENG1', 'not reported')")
            conn.commit()

        engine = PaperMiner(db_path=tmp.name)
        # Keep tests hermetic: no network full-text retrieval by default
        engine.lit_client.fetch_fulltext_xml = Mock(return_value=None)
        engine.lit_client.fetch_fulltext_pdf = Mock(return_value=None)
        yield engine


def test_process_paper_multi_gene(temp_engine):
    # A single paper characterizing two genes (MoPMK1 and MoENG1)
    article = {
        "pmid": "99887766",
        "pmcid": "PMC998877",
        "doi": "10.1000/test",
        "title": "Characterization of MoPMK1 and MoENG1 in Magnaporthe oryzae Pathogenicity",
        "abstract": (
            "We investigated Magnaporthe oryzae infection mechanisms. "
            "Deletion of MoPMK1 caused complete loss of virulence on rice leaves. "
            "In addition, the Δeng1 mutant exhibited altered conidiation."
        ),
        "pub_year": 2024,
        "is_open_access": True
    }

    result = temp_engine.process_paper(article, verbose=False)

    assert result["status"] == "curated"
    assert result["genes_extracted"] >= 1
    assert temp_engine.db.is_paper_processed(pmid="99887766") is True

    # Check that database records both extractions and updated master genes
    stats = temp_engine.db.get_stats()
    assert stats["papers_triaged"] == 1
    assert stats["experimental_papers"] == 1
    assert stats["total_evidence_extractions"] >= 1


def test_process_paper_non_experimental(temp_engine):
    article = {
        "pmid": "11223344",
        "pmcid": None,
        "doi": None,
        "title": "A computational survey of fungal secondary metabolites",
        "abstract": "We mined public databases for secondary metabolite gene clusters across 50 ascomycetes.",
        "pub_year": 2023,
        "is_open_access": False
    }

    result = temp_engine.process_paper(article, verbose=False)

    assert result["status"] == "non_experimental"
    assert result["genes_extracted"] == 0
    assert temp_engine.db.is_paper_processed(pmid="11223344") is True

    stats = temp_engine.db.get_stats()
    assert stats["papers_triaged"] == 1
    assert stats["experimental_papers"] == 0
    assert stats["total_evidence_extractions"] == 0


def test_process_paper_idempotency(temp_engine):
    article = {
        "pmid": "55443322",
        "title": "Role of MoPMK1 in appressorium formation in Magnaporthe oryzae",
        "abstract": "Targeted deletion of MoPMK1 resulted in defective appressorium formation and lost virulence.",
        "pub_year": 2023,
        "is_open_access": False
    }

    first_run = temp_engine.process_paper(article, verbose=False)
    assert first_run["status"] == "curated"

    # Second run should skip immediately
    second_run = temp_engine.process_paper(article, verbose=False)
    assert second_run["status"] == "already_processed"
    assert second_run["genes_extracted"] == 0


def test_mine_papers_mocked(temp_engine):
    mock_articles = [
        {
            "pmid": "1001",
            "title": "MoPMK1 deletion mutant in Magnaporthe oryzae",
            "abstract": "The Δmopmk1 deletion mutant lost virulence on rice.",
            "pub_year": 2023,
            "is_open_access": False
        },
        {
            "pmid": "1002",
            "title": "Review of fungal plant pathogens",
            "abstract": "This article reviews the epidemiology of cereal fungi without experimental assays.",
            "pub_year": 2022,
            "is_open_access": False
        }
    ]

    with patch.object(temp_engine.lit_client, "stream_papers", return_value=iter(mock_articles)):
        stats = temp_engine.mine_papers(max_papers=2, verbose=False)

        assert stats["streamed"] == 2
        assert stats["curated"] == 1
        assert stats["non_experimental"] == 1
        assert stats["extractions_created"] >= 1


def test_mine_papers_year_by_year(temp_engine):
    calls = []

    def fake_stream(**kwargs):
        yr = kwargs.get("pub_year")
        calls.append(yr)
        return iter([
            {
                "pmid": f"yr-{yr}",
                "title": f"MoPMK1 study in {yr}",
                "abstract": "Deletion of MoPMK1 in Magnaporthe oryzae abolished virulence.",
                "pub_year": yr,
                "is_open_access": False,
            }
        ])

    with patch.object(temp_engine.lit_client, "stream_papers", side_effect=fake_stream):
        stats = temp_engine.mine_papers(max_papers=5, year_from=2022, year_to=2023, verbose=False)

    assert calls == [2022, 2023]
    assert stats["streamed"] == 2
    assert stats["curated"] == 2
    assert stats["extractions_created"] == 2


def test_process_paper_inline_alias_extraction(temp_engine):
    # Seed an unnamed gene (gene_name=None) like MoAcb1 in 2023 literature
    with temp_engine.db.get_connection() as conn:
        conn.execute("INSERT INTO genes (gene, gene_name, process_affected) VALUES ('MGG_06177', NULL, 'Not reported')")
        conn.commit()
    temp_engine.db._linker = None  # refresh index
    temp_engine.linker = temp_engine.db.linker

    temp_engine.lit_client.fetch_fulltext_xml = Mock(return_value={
        "other": "Introduction: We identified MoAcb1 (MGG_06177) in Magnaporthe oryzae.",
        "results": "Deletion of ∆Moacb1 caused severe conidiation defects and reduced virulence on rice.",
        "methods": "Strains were cultured on CM."
    })

    article = {
        "pmid": "37187543",
        "pmcid": "PMC999999",
        "title": "Acyl-coenzyme A binding protein MoAcb1 regulates conidiation and pathogenicity",
        "abstract": "We characterized MoAcb1 in Magnaporthe oryzae.",
        "pub_year": 2023,
        "is_open_access": True,
    }
    res = temp_engine.process_paper(article, verbose=False)
    assert res["status"] == "curated"
    assert res["genes_extracted"] == 1
    rec = res["extractions"][0]
    assert rec.gene_id == "MGG_06177"
    assert rec.gene_name == "MoAcb1"
    assert rec.confidence_tier.value in ("tier_1_high", "tier_2_high")


def test_unclassified_only_is_default_main_mode(tmp_path):
    db_file = tmp_path / "filter.db"
    db = DatabaseManager(str(db_file))
    with db.get_connection() as conn:
        # Already classified in MagnaGenes baseline
        conn.execute(
            "INSERT INTO genes (gene, gene_name, process_affected, baseline_process) "
            "VALUES ('MGG_09565', 'PMK1', 'virulence', 'virulence')"
        )
        # Unclassified in MagnaGenes baseline
        conn.execute(
            "INSERT INTO genes (gene, gene_name, process_affected, baseline_process) "
            "VALUES ('MGG_01001', 'ENG1', 'Not reported', 'Not reported')"
        )
        conn.commit()

    engine = PaperMiner(db_path=str(db_file))  # default unclassified_only=True
    engine.lit_client.fetch_fulltext_xml = Mock(return_value=None)
    engine.lit_client.fetch_fulltext_pdf = Mock(return_value=None)

    article = {
        "pmid": "777777",
        "title": "Study of MoPMK1 and Δeng1 in Magnaporthe oryzae",
        "abstract": "Deletion of MoPMK1 abolished virulence. Targeted deletion of Δeng1 caused conidiation defects.",
        "pub_year": 2023,
        "is_open_access": False,
    }
    res = engine.process_paper(article, verbose=False)
    assert res["status"] == "curated"
    assert res["skipped_baseline_genes"] == 1
    assert [r.gene_id for r in res["extractions"]] == ["MGG_01001"]


def test_double_check_additional_functions_and_knockout_columns(tmp_path):
    db_file = tmp_path / "dc.db"
    db = DatabaseManager(str(db_file))
    with db.get_connection() as conn:
        conn.execute(
            "INSERT INTO genes (gene, gene_name, process_affected, baseline_process) "
            "VALUES ('MGG_09565', 'PMK1', 'virulence', 'virulence')"
        )
        conn.commit()

    article = {
        "pmid": "991122",
        "title": "New role of MoPMK1 in conidiation in Magnaporthe oryzae",
        "abstract": "Deletion of MoPMK1 caused severe conidiation defects in Magnaporthe oryzae.",
        "pub_year": 2025,
        "is_open_access": False,
    }

    # 1. Default unclassified-only pass skips baseline-classified PMK1
    main_engine = PaperMiner(db_path=str(db_file))
    main_engine.lit_client.fetch_fulltext_xml = Mock(return_value=None)
    main_engine.lit_client.fetch_fulltext_pdf = Mock(return_value=None)
    r1 = main_engine.process_paper(article, verbose=False)
    assert r1["skipped_baseline_genes"] == 1
    assert r1["genes_extracted"] == 0

    # 2. Double-check pass inspects same paper for additional functions & knockouts
    dc_engine = PaperMiner(db_path=str(db_file), check_additional_functions=True)
    dc_engine.lit_client.fetch_fulltext_xml = Mock(return_value=None)
    dc_engine.lit_client.fetch_fulltext_pdf = Mock(return_value=None)
    r2 = dc_engine.process_paper(article, verbose=False)
    assert r2["status"] == "curated"
    assert r2["additional_functions_found"] == 1

    gene = db.get_gene("MGG_09565")
    assert gene["process_affected"] == "virulence"  # baseline primary preserved
    assert "conidiation" in gene["additional_functions"]
    assert "PMID:991122" in gene["additional_function_references"]
    assert gene["knockout_reported"] == 1
    assert "conidiation" in gene["knockout_phenotype"]

    # 3. Re-running double-check skips already double-checked paper
    r3 = dc_engine.process_paper(article, verbose=False)
    assert r3["status"] == "already_processed"


def test_mine_papers_days_window(temp_engine):
    captured = {}

    def fake_stream(**kwargs):
        captured.update(kwargs)
        return iter([])

    with patch.object(temp_engine.lit_client, "stream_papers", side_effect=fake_stream):
        temp_engine.mine_papers(max_papers=5, days=90, verbose=False)

    assert "date_from" in captured and captured["date_from"] is not None
    assert "date_to" in captured and captured["date_to"] is not None
