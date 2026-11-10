"""
Unit tests for EvidenceExtractor triage, grounding, and quote verification.
"""

import pytest
from magnamined.extractor import EvidenceExtractor
from magnamined.schema import PhenotypeCategory, EvidenceConfidence


def test_verify_verbatim_quote():
    source = "Deletion of MGG_01001 resulted in loss of appressorium formation and avirulence on susceptible rice cultivars."
    
    # Exact match
    assert EvidenceExtractor.verify_verbatim_quote(
        "Deletion of MGG_01001 resulted in loss of appressorium formation",
        source
    ) is True

    # Whitespace variance
    assert EvidenceExtractor.verify_verbatim_quote(
        "Deletion  of  MGG_01001 resulted in loss\nof appressorium formation",
        source
    ) is True

    # Fabricated / hallucinated quote
    assert EvidenceExtractor.verify_verbatim_quote(
        "MGG_01001 was shown to regulate cell cycle progression during mitosis",
        source
    ) is False


def test_triage_abstract():
    extractor = EvidenceExtractor()
    title = "Functional characterization of ENG1 in Magnaporthe oryzae pathogenicity"
    abstract = "Targeted deletion of the MGG_01001 gene abolished virulence and appressorium penetration on host rice."

    triage = extractor.triage_abstract("12345", title, abstract)
    assert triage.is_experimental is True
    assert triage.target_organism_confirmed is True
    assert PhenotypeCategory.VIRULENCE in triage.suggested_phenotypes
    assert "MGG_01001" in triage.mentioned_genes


def test_extract_evidence_heuristic():
    extractor = EvidenceExtractor()
    title = "Role of ENG1 in blast disease"
    abstract = "We characterized MGG_01001 in Magnaporthe oryzae."
    snippets = [
        "Results: Targeted deletion of MGG_01001 caused reduced virulence on rice seedlings. Genetic complementation rescued the wild type."
    ]

    record = extractor.extract_evidence_from_snippets(
        gene_id="MGG_01001",
        gene_name="ENG1",
        title=title,
        abstract=abstract,
        snippets=snippets,
        pmid="123456"
    )

    assert record is not None
    assert record.gene_id == "MGG_01001"
    assert record.process_affected == PhenotypeCategory.VIRULENCE
    assert record.complementation_tested is True
    assert record.confidence_tier == EvidenceConfidence.TIER_1
    assert "MGG_01001" in record.verbatim_quote


def test_verify_verbatim_quote_fuzzy():
    # Handles HTML/XML tags and unicode hyphens
    source = "The <i>MoPMK1</i>-mediated signaling pathway was essential for plant infection."
    quote = "The MoPMK1-mediated signaling pathway was essential for plant infection."
    assert EvidenceExtractor.verify_verbatim_quote(quote, source) is True

    # Handles soft hyphens and line-break whitespace
    source_wrapped = "Deletion of MoENG1 resulted in reduced conidia-\ntion on oatmeal agar."
    quote_clean = "Deletion of MoENG1 resulted in reduced conidiation on oatmeal agar."
    assert EvidenceExtractor.verify_verbatim_quote(quote_clean, source_wrapped) is True


def test_extract_evidence_negative_result():
    extractor = EvidenceExtractor()
    title = "Analysis of putative transcription factor MGG_05000"
    abstract = "We deleted MGG_05000 in Magnaporthe oryzae."
    snippets = [
        "Results: The MGG_05000 deletion mutant grew normally and was dispensable for virulence on barley seedlings."
    ]

    record = extractor.extract_evidence_from_snippets(
        gene_id="MGG_05000",
        gene_name="TF1",
        title=title,
        abstract=abstract,
        snippets=snippets,
        pmid="778899"
    )

    assert record is not None
    assert record.is_dispensable is True
    assert record.process_affected == PhenotypeCategory.DISPENSABLE
    assert "dispensable for virulence" in record.verbatim_quote


def test_extract_evidence_pleiotropy():
    extractor = EvidenceExtractor()
    title = "Dual roles of MoPMK1 in infection"
    abstract = "MoPMK1 is a key MAP kinase in Magnaporthe oryzae."
    snippets = [
        "Results: Deletion of MoPMK1 resulted in defective appressorium formation and altered conidiation in culture."
    ]

    record = extractor.extract_evidence_from_snippets(
        gene_id="MGG_09565",
        gene_name="PMK1",
        title=title,
        abstract=abstract,
        snippets=snippets,
        pmid="889900"
    )

    assert record is not None
    assert record.process_affected == PhenotypeCategory.APPRESSORIUM
    assert PhenotypeCategory.CONIDIATION in record.secondary_processes
    assert len(record.secondary_processes) >= 1


def test_llm_extraction_preserves_metadata_and_rejects_hallucination():
    import json
    from unittest.mock import MagicMock, patch

    extractor = EvidenceExtractor(api_key="dummy-key")
    snippets = ["Deletion of MoPMK1 (MGG_09565) abolished appressorium formation and virulence."]

    # 1. Valid grounded LLM JSON (omitting publication metadata)
    valid_payload = json.dumps({
        "gene_id": "MGG_09565",
        "gene_name": "PMK1",
        "process_affected": "appressorium formation",
        "secondary_processes": ["virulence"],
        "detailed_phenotype": "Abolished appressorium formation and virulence.",
        "is_dispensable": False,
        "experimental_method": "knockout/deletion",
        "complementation_tested": False,
        "confidence_tier": "tier_2_high",
        "verbatim_quote": "Deletion of MoPMK1 (MGG_09565) abolished appressorium formation and virulence."
    })

    mock_client = MagicMock()
    mock_client.models.generate_content.return_value = MagicMock(text=valid_payload)

    with patch("google.genai.Client", return_value=mock_client):
        rec = extractor.extract_evidence_from_snippets(
            gene_id="MGG_09565", gene_name="PMK1", title="PMK1 Study",
            abstract="Abstract.", snippets=snippets, pmid="42381406",
            pmcid="PMC999", doi="10.1/pmk1", pub_year=2025
        )
        assert rec is not None
        assert rec.detailed_phenotype == "Abolished appressorium formation and virulence."
        assert rec.pmid == "42381406"
        assert rec.pmcid == "PMC999"
        assert rec.doi == "10.1/pmk1"
        assert rec.paper_title == "PMK1 Study"
        assert rec.publication_year == 2025

    # 2. Hallucinated quote is rejected and warns before falling back
    hallucinated_payload = json.dumps({
        "gene_id": "MGG_09565",
        "process_affected": "virulence",
        "detailed_phenotype": "Hallucinated claim.",
        "experimental_method": "knockout/deletion",
        "confidence_tier": "tier_1_high",
        "verbatim_quote": "This completely fabricated sentence does not appear anywhere in the paper."
    })
    mock_client.models.generate_content.return_value = MagicMock(text=hallucinated_payload)

    with patch("google.genai.Client", return_value=mock_client):
        with pytest.warns(RuntimeWarning, match="grounding verification"):
            rec2 = extractor.extract_evidence_from_snippets(
                gene_id="MGG_09565", gene_name="PMK1", title="PMK1 Study",
                abstract="Abstract.", snippets=snippets, pmid="42381406"
            )
        assert rec2 is not None
        assert "Abolished" not in rec2.detailed_phenotype
        assert extractor.verify_verbatim_quote(rec2.verbatim_quote, snippets[0])


def test_heuristic_rejects_bystanders_and_leakage():
    extractor = EvidenceExtractor()

    # 1. qPCR expression target + localization marker in another gene's deletion paper
    snippet_bystanders = [
        "Deletion of MoSUN4 caused reduced conidiation and attenuated virulence. "
        "We further showed that expression levels of BUF1 and ALB1 were downregulated by qRT-PCR in Δmosun4. "
        "Accordingly, fluorescence was co-localized with the well-established apoplastic effector Bas4."
    ]
    assert extractor.extract_evidence_from_snippets(
        gene_id="MGG_02252", gene_name="BUF1", title="MoSun4 study",
        abstract="MoSun4 regulates blast infection.", snippets=snippet_bystanders, pmid="42646132"
    ) is None
    assert extractor.extract_evidence_from_snippets(
        gene_id="MGG_10914", gene_name="BAS4", title="MoSun4 study",
        abstract="MoSun4 regulates blast infection.", snippets=snippet_bystanders, pmid="42646132"
    ) is None

    # 2. Multi-gene figure / heatmap list
    snippet_list = [
        "The right side label reads, Glucanase genes, followed by gene identifiers listed top to bottom: "
        "MGG_10083, MGG_01001, MGG_10400, MGG_04689."
    ]
    assert extractor.extract_evidence_from_snippets(
        gene_id="MGG_01001", gene_name="ENG1", title="Kmt1 study",
        abstract="MoKmt1 regulates gene expression.", snippets=snippet_list, pmid="42727022"
    ) is None

    # 3. Distant 'dispensable' keyword does not leak onto another gene, and OCR comments are stripped
    snippet_ocr = [
        "Deletion of MGG_03580 (MrATG4) caused severe conidiation defects.\n"
        "<!-- Start of picture text -->\n"
        "In contrast, ATG5 is discussed here only as a conserved autophagy factor. "
        "Earlier work showed another kinase was dispensable for virulence."
    ]
    atg4_rec = extractor.extract_evidence_from_snippets(
        gene_id="MGG_03580", gene_name="ATG4", title="MrATG4 study",
        abstract="MrATG4 regulates conidiation.", snippets=snippet_ocr, pmid="41557747"
    )
    assert atg4_rec is not None
    assert "<!--" not in atg4_rec.verbatim_quote
    assert atg4_rec.is_dispensable is False
    assert atg4_rec.process_affected == PhenotypeCategory.CONIDIATION

    atg5_rec = extractor.extract_evidence_from_snippets(
        gene_id="MGG_09262", gene_name="ATG5", title="MrATG4 study",
        abstract="MrATG4 regulates conidiation.", snippets=snippet_ocr, pmid="41557747"
    )
    assert atg5_rec is None or atg5_rec.is_dispensable is False


def test_cyp51b_prior_study_unicode_hyphen_and_no_effect_negations():
    extractor = EvidenceExtractor()

    # 1. Unicode hyphen RT‐qPCR DEG list must not be extracted as pathogenicity defect
    deg_snippet = [
        "Analysis by reverse transcription quantitative PCR (RT\u2010qPCR) confirmed the reduced expression "
        "in the Mopah1 mutant of DEGs involved in fungal growth (FLB3, HOX1, HOX4, CRZ1), "
        "conidiation (TDG1, TDG7, FAR2, CYP51B), and pathogenicity (RSY1, Dac, CYP51B, PTH11) (Figure 6c)."
    ]
    assert extractor.extract_evidence_from_snippets(
        gene_id="MGG_04432", gene_name="CYP51b", title="MoPah1 study",
        abstract="MoPah1 regulates fungal development.", snippets=deg_snippet, pmid="35191164"
    ) is None

    # 2. Prior-study background citation must be filtered out
    prior_snippet = [
        "A previous study indicated that MoCYP51A is a core determinant of DMI sensitivity in M. oryzae, "
        "while deletion of MoCYP51B has little phenotypic effect [20]."
    ]
    assert extractor.extract_evidence_from_snippets(
        gene_id="MGG_04432", gene_name="CYP51b", title="DMI study",
        abstract="Azole sensitivity in Magnaporthe oryzae.", snippets=prior_snippet, pmid="42750243"
    ) is None

    # 3. Direct experimental 'little phenotypic effect' or 'no effect on pathogenicity' -> DISPENSABLE, intact sentence
    direct_neg_snippet = [
        "The pathogenicity of these mutant strains was analyzed, and results showed that the deletion of "
        "MoCYP51B has no effect on the pathogenicity of M. oryzae."
    ]
    neg_rec = extractor.extract_evidence_from_snippets(
        gene_id="MGG_04432", gene_name="CYP51b", title="CYP51B deletion study",
        abstract="Targeted deletion of MoCYP51B.", snippets=direct_neg_snippet, pmid="999001"
    )
    assert neg_rec is not None
    assert neg_rec.is_dispensable is True
    assert neg_rec.process_affected == PhenotypeCategory.DISPENSABLE
    assert "M. oryzae." in neg_rec.verbatim_quote

    # 4. Pull-down / co-IP mention in another mutant must not assign phenotype to Mst7
    pulldown_snippet = [
        "However, the PARylation mutant GRF1m3-5 did not pull-down either Mst7 or Pmk1 (Fig. 3)."
    ]
    assert extractor.extract_evidence_from_snippets(
        gene_id="MGG_00800", gene_name="MST7", title="GRF1 study",
        abstract="GRF1 interacts with MAPK components.", snippets=pulldown_snippet, pmid="39277621"
    ) is None
