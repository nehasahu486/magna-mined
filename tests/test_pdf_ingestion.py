"""
Unit tests for Phase 5: PDF ingestion, source cascade, and abstract-only fallback.
"""

import tempfile
from unittest.mock import Mock

import pytest

from magnamined.paper_miner import PaperMiner
from magnamined.database import DatabaseManager
from magnamined.literature import EuropePMCClient
from magnamined.schema import EvidenceConfidence, EvidenceSource


@pytest.fixture
def engine():
    with tempfile.NamedTemporaryFile(suffix=".db") as tmp:
        db = DatabaseManager(tmp.name)
        with db.get_connection() as conn:
            conn.execute("INSERT INTO genes (gene, gene_name, process_affected) VALUES ('MGG_09565', 'PMK1', 'not reported')")
            conn.commit()
        eng = PaperMiner(db_path=tmp.name)
        eng.lit_client.fetch_fulltext_xml = Mock(return_value=None)
        eng.lit_client.fetch_fulltext_pdf = Mock(return_value=None)
        yield eng


def _article(pmid, pmcid="PMC1"):
    return {
        "pmid": pmid,
        "pmcid": pmcid,
        "doi": f"10.1000/{pmid}",
        "title": "Role of MoPMK1 in Magnaporthe oryzae",
        "abstract": "Deletion of MoPMK1 in Magnaporthe oryzae caused reduced virulence on rice.",
        "pub_year": 2024,
        "is_open_access": True,
    }


def test_split_markdown_sections():
    md = (
        "# Abstract\nSummary text.\n\n"
        "## Results\nThe MoPMK1 mutant was avirulent.\n\n"
        "## Materials and Methods\nStrains were grown on oatmeal agar.\n\n"
        "## References\n1. Someone et al."
    )
    s = EuropePMCClient.split_markdown_sections(md)
    assert "Summary text." in s["abstract"]
    assert "avirulent" in s["results"]
    assert "oatmeal agar" in s["methods"]
    assert "Someone" not in s["results"] + s["methods"]
    assert s["other"] == md  # full document always preserved


def test_split_markdown_sections_pymupdf4llm_bold_headings():
    # Matches actual heading style emitted by pymupdf4llm on PLOS PDFs
    md = (
        "## **Abstract** \nAbstract body.\n\n"
        "#### **Results** \nResults body.\n\n"
        "#### **Materials and methods** \nMethods body.\n\n"
        "#### **Discussion** \nDiscussion body.\n"
    )
    s = EuropePMCClient.split_markdown_sections(md)
    assert "Abstract body." in s["abstract"]
    assert "Results body." in s["results"]
    assert "Methods body." in s["methods"]
    assert "Discussion body." in s["discussion"]


def test_pdf_bytes_to_markdown_roundtrip():
    pymupdf = pytest.importorskip("pymupdf")
    pytest.importorskip("pymupdf4llm")
    doc = pymupdf.open()
    page = doc.new_page()
    page.insert_text((72, 72), "The MoPMK1 deletion mutant was avirulent on rice.")
    pdf_bytes = doc.tobytes()
    doc.close()

    md = EuropePMCClient.pdf_bytes_to_markdown(pdf_bytes)
    assert md is not None
    assert "MoPMK1" in md


def test_pick_pdf_url():
    item = {"fullTextUrlList": {"fullTextUrl": [
        {"documentStyle": "html", "availability": "Open access", "url": "http://x/html"},
        {"documentStyle": "pdf", "availability": "Subscription required", "url": "http://x/paywalled"},
        {"documentStyle": "pdf", "availability": "Open access", "url": "https://europepmc.org/articles/PMC1?pdf=render"},
        {"documentStyle": "pdf", "availability": "Open access", "url": "http://x/open.pdf"},
    ]}}
    # Europe PMC render links are skipped (they 403 for scripted clients)
    assert EuropePMCClient._pick_pdf_url(item) == "http://x/open.pdf"
    assert EuropePMCClient._pick_pdf_url({}) is None


def test_unpaywall_requires_email(monkeypatch):
    monkeypatch.delenv("MAGNAMINED_CONTACT_EMAIL", raising=False)
    client = EuropePMCClient()
    client.session.get = Mock()
    assert client.find_unpaywall_pdf("10.1/x") is None
    client.session.get.assert_not_called()


def test_unpaywall_lookup(monkeypatch):
    monkeypatch.setenv("MAGNAMINED_CONTACT_EMAIL", "curator@example.org")
    client = EuropePMCClient()
    resp = Mock(status_code=200)
    resp.json.return_value = {"best_oa_location": {"url_for_pdf": "https://pub.example/paper.pdf"}}
    client.session.get = Mock(return_value=resp)
    assert client.find_unpaywall_pdf("10.1/x") == "https://pub.example/paper.pdf"


def test_ingest_local_pdf(engine, tmp_path):
    pymupdf = pytest.importorskip("pymupdf")
    pytest.importorskip("pymupdf4llm")
    doc = pymupdf.open()
    page = doc.new_page()
    page.insert_text((72, 72), "Results")
    page.insert_text((72, 100), "The MoPMK1 deletion mutant showed reduced virulence on rice.")
    pdf = tmp_path / "paper.pdf"
    doc.save(str(pdf))
    doc.close()

    res = engine.ingest_local_pdf(str(pdf), title="Local MoPMK1 study", verbose=False)
    assert res["status"] == "curated"
    assert res["extractions"]
    assert all(r.evidence_source == EvidenceSource.PDF for r in res["extractions"])

    # Re-ingesting replaces rather than duplicates extractions
    engine.ingest_local_pdf(str(pdf), title="Local MoPMK1 study", verbose=False)
    with engine.db.get_connection() as conn:
        n = conn.execute("SELECT COUNT(*) FROM extractions WHERE pmid = 'local:paper.pdf'").fetchone()[0]
    assert n == len(res["extractions"])


def test_cascade_prefers_xml(engine):
    engine.lit_client.fetch_fulltext_xml = Mock(return_value={
        "results": "Deletion of MoPMK1 caused reduced virulence. Complementation rescued the phenotype.",
        "methods": "Strains were grown.",
    })
    res = engine.process_paper(_article("1"), verbose=False)
    assert res["status"] == "curated"
    assert all(r.evidence_source == EvidenceSource.JATS_XML for r in res["extractions"])
    engine.lit_client.fetch_fulltext_pdf.assert_not_called()


def test_cascade_falls_back_to_pdf(engine):
    engine.lit_client.fetch_fulltext_pdf = Mock(return_value={
        "results": "", "methods": "",
        "other": "The MoPMK1 deletion mutant showed reduced virulence on rice seedlings.",
    })
    res = engine.process_paper(_article("2"), verbose=False)
    assert res["status"] == "curated"
    assert res["extractions"], "expected at least one extraction"
    assert all(r.evidence_source == EvidenceSource.PDF for r in res["extractions"])


def test_abstract_only_fallback_caps_tier(engine):
    res = engine.process_paper(_article("3"), verbose=False)
    assert res["status"] == "curated"
    assert res["extractions"]
    for r in res["extractions"]:
        assert r.evidence_source == EvidenceSource.ABSTRACT_ONLY
        assert r.confidence_tier == EvidenceConfidence.TIER_4

    # Persisted to the ledger and extractions table
    with engine.db.get_connection() as conn:
        row = conn.execute("SELECT evidence_source FROM processed_papers WHERE pmid = '3'").fetchone()
        assert row["evidence_source"] == "abstract_only"
        src = conn.execute("SELECT DISTINCT evidence_source FROM extractions").fetchall()
        assert [s["evidence_source"] for s in src] == ["abstract_only"]


def test_pdf_disabled_skips_pdf(engine):
    engine.enable_pdf = False
    engine.process_paper(_article("4"), verbose=False)
    engine.lit_client.fetch_fulltext_pdf.assert_not_called()


def test_migration_adds_evidence_source_column():
    with tempfile.NamedTemporaryFile(suffix=".db") as tmp:
        import sqlite3
        conn = sqlite3.connect(tmp.name)
        # Simulate a pre-Phase-5 schema
        conn.execute("CREATE TABLE extractions (id INTEGER PRIMARY KEY, gene_id INTEGER, process_affected TEXT, detailed_phenotype TEXT, experimental_method TEXT, confidence_tier TEXT, verbatim_quote TEXT, pmid TEXT, pmcid TEXT, doi TEXT, complementation_tested BOOLEAN, paper_title TEXT, publication_year INTEGER, extracted_at TIMESTAMP)")
        conn.commit()
        conn.close()
        db = DatabaseManager(tmp.name)
        with db.get_connection() as c:
            cols = {r[1] for r in c.execute("PRAGMA table_info(extractions)")}
            assert "evidence_source" in cols
