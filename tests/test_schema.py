"""
Unit tests for data schemas and Pydantic validation.
"""

import pytest
from magnamined.schema import (
    PhenotypeCategory,
    ExperimentalMethod,
    EvidenceConfidence,
    GeneEvidenceRecord,
    AbstractTriageResult
)


def test_schema_valid_record():
    record = GeneEvidenceRecord(
        gene_id="MGG_01001",
        gene_name="ENG1",
        process_affected=PhenotypeCategory.VIRULENCE,
        detailed_phenotype="Δeng1 showed reduced appressorium turgor and impaired virulence",
        experimental_method=ExperimentalMethod.KNOCKOUT_DELETION,
        complementation_tested=True,
        confidence_tier=EvidenceConfidence.TIER_1,
        verbatim_quote="Deletion of ENG1 caused reduced lesion development on rice seedlings.",
        pmid="12345678",
        publication_year=2023
    )
    assert record.gene_id == "MGG_01001"
    assert record.process_affected == PhenotypeCategory.VIRULENCE
    assert record.complementation_tested is True
    assert record.confidence_tier == EvidenceConfidence.TIER_1


def test_schema_triage_result():
    triage = AbstractTriageResult(
        pmid="12345678",
        is_experimental=True,
        target_organism_confirmed=True,
        mentioned_genes=["MGG_01001", "ENG1"],
        suggested_phenotypes=[PhenotypeCategory.VIRULENCE, PhenotypeCategory.APPRESSORIUM],
        relevance_score=0.9,
        rationale="Deletion of ENG1 mutant tested in blast fungus."
    )
    assert triage.is_experimental is True
    assert len(triage.suggested_phenotypes) == 2
