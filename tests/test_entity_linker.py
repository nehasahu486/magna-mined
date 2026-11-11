"""
Unit tests for EntityLinker gene disambiguation and alias resolution.
"""

import pytest
from magnamined.entity_linker import EntityLinker, GeneMention, BIOLOGICAL_STOPWORDS
from magnamined.database import DatabaseManager


@pytest.fixture(scope="module")
def linker():
    db = DatabaseManager("magnaporthe.db")
    return EntityLinker(db=db)


def test_normalize_token(linker):
    # Systematic ID normalization
    assert linker.normalize_token("MGG_01001") == "MGG_01001"
    assert linker.normalize_token("mgg_01001T0") == "MGG_01001"
    assert linker.normalize_token("MGG-01001") == "MGG_01001"
    assert linker.normalize_token("MGG01001") == "MGG_01001"
    assert linker.normalize_token("ΔMGG_01001") == "MGG_01001"

    # Fungal mutant and organism prefix normalization
    assert linker.normalize_token("MoENG1") == "ENG1"
    assert linker.normalize_token("PoENG1") == "ENG1"
    assert linker.normalize_token("Δeng1") == "eng1"
    assert linker.normalize_token("delta-eng1") == "eng1"
    assert linker.normalize_token("eng1Δ") == "eng1"
    assert linker.normalize_token("pmk1-null") == "pmk1"
    assert linker.normalize_token("Δmopmk1") == "pmk1"

    # Ensure gene symbols ending with T and digits are not mangled
    assert linker.normalize_token("MST12") == "MST12"
    assert linker.normalize_token("MoMST12") == "MST12"
    assert linker.normalize_token("Δmst12") == "mst12"


def test_resolve_canonical_genes(linker):
    assert linker.resolve("MGG_01001") == "MGG_01001"
    assert linker.resolve("ENG1") == "MGG_01001"
    assert linker.resolve("MoENG1") == "MGG_01001"
    assert linker.resolve("Δeng1") == "MGG_01001"

    assert linker.resolve("PMK1") == "MGG_09565"
    assert linker.resolve("MoPMK1") == "MGG_09565"
    assert linker.resolve("Δmopmk1") == "MGG_09565"
    assert linker.resolve("pmk1Δ") == "MGG_09565"

    assert linker.resolve("MAC1") == "MGG_09898"
    assert linker.resolve("MoMAC1") == "MGG_09898"
    assert linker.resolve("MST12") == "MGG_12958"
    assert linker.resolve("PTH11") == "MGG_05871"
    assert linker.resolve("CPKA") == "MGG_06368"


def test_resolve_legacy_broad_ids(linker):
    # Older papers citing Broad Institute v6 IDs
    assert linker.resolve("MGG_00692") == "MGG_09565"  # PMK1
    assert linker.resolve("MGG_01130") == "MGG_09898"  # MAC1
    assert linker.resolve("MGG_02846") == "MGG_12958"  # MST12
    assert linker.resolve("MGG_00659") == "MGG_06368"  # CPKA


def test_find_mentions_in_text(linker):
    text = (
        "In Magnaporthe oryzae, the MAP kinase MoPMK1 (MGG_09565) regulates appressorium turgor. "
        "Deletion mutant Δmopmk1 is completely avirulent. "
        "Conversely, Δeng1 shows reduced conidiation. "
        "We also evaluated MAC1 and downstream target MoCPKA. "
        "Controls were transformed with empty GFP vector."
    )

    mentions = linker.find_mentions(text)
    detected_ids = [m.canonical_id for m in mentions]

    assert "MGG_09565" in detected_ids  # PMK1
    assert "MGG_01001" in detected_ids  # ENG1
    assert "MGG_09898" in detected_ids  # MAC1
    assert "MGG_06368" in detected_ids  # CPKA

    # Check mutant identification
    pmk1_mut = next(m for m in mentions if m.raw_text == "Δmopmk1")
    assert pmk1_mut.is_mutant is True
    assert pmk1_mut.canonical_id == "MGG_09565"

    eng1_mut = next(m for m in mentions if m.raw_text == "Δeng1")
    assert eng1_mut.is_mutant is True
    assert eng1_mut.canonical_id == "MGG_01001"

    # Stopwords like GFP should not be extracted as genes
    assert not any(m.raw_text == "GFP" for m in mentions)


def test_get_aliases_for_gene(linker):
    pmk1_aliases = linker.get_aliases_for_gene("MGG_09565")
    assert "MGG_09565" in pmk1_aliases
    assert "MGG_00692" in pmk1_aliases  # legacy
    assert "PMK1" in pmk1_aliases
    assert "MoPMK1" in pmk1_aliases


def test_database_get_gene_resolution(linker):
    # Verify DatabaseManager.get_gene leverages linker
    db = linker.db
    by_mutant = db.get_gene("Δmopmk1")
    assert by_mutant is not None
    assert by_mutant["gene"] == "MGG_09565"

    by_legacy = db.get_gene("MGG_00692")
    assert by_legacy is not None
    assert by_legacy["gene"] == "MGG_09565"


def test_extract_inline_aliases_for_unnamed_genes(linker):
    # MGG_06177 (MoAcb1) and MGG_01003 (MoPct1) have no gene_name in baseline DB
    text = (
        "We identified an Acb1 homolog MoAcb1 (MGG_06177) and studied MGG_01003 (MoPct1). "
        "Note that the hypothetical protein (MGG_16142) should not bind the word 'protein'. "
        "Deletion of ∆Moacb1 reduced conidiation, and ΔMopct1 impaired vegetative growth."
    )
    inline = linker.extract_inline_aliases(text)
    assert "MGG_06177" in inline
    assert "MoAcb1" in inline["MGG_06177"]
    assert "MGG_01003" in inline
    assert "MoPct1" in inline["MGG_01003"]
    assert "MGG_16142" not in inline  # 'protein' rejected as generic word

    mentions = linker.find_mentions(text, extra_aliases=inline)
    by_id = {m.canonical_id for m in mentions}
    assert "MGG_06177" in by_id
    assert "MGG_01003" in by_id
    assert any(m.raw_text == "∆Moacb1" and m.canonical_id == "MGG_06177" and m.is_mutant for m in mentions)
    assert any(m.raw_text == "ΔMopct1" and m.canonical_id == "MGG_01003" and m.is_mutant for m in mentions)


def test_gene_aliases_table_indexed_for_unnamed_genes(tmp_path):
    db = DatabaseManager(str(tmp_path / "alias_test.db"))
    with db.get_connection() as conn:
        conn.execute("INSERT INTO genes (gene, gene_name, process_affected) VALUES ('MGG_06177', NULL, 'Not reported')")
        conn.commit()
    added = db.add_gene_aliases_bulk([("MGG_06177", "MoPcs60")], source="uniprot")
    assert added == 1

    mentions = db.linker.find_mentions("Deletion of ΔMoPcs60 caused conidiation defects in Magnaporthe oryzae.")
    assert [m.canonical_id for m in mentions] == ["MGG_06177"]
    assert mentions[0].canonical_name == "MoPcs60"
