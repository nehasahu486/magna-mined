"""
Data schemas and Pydantic models for Magna-Mined.
"""

from enum import Enum
from typing import List, Optional
from pydantic import BaseModel, Field


class PhenotypeCategory(str, Enum):
    """Controlled vocabulary for primary biological processes affected."""
    VIRULENCE = "virulence"
    APPRESSORIUM = "appressorium formation"
    CONIDIATION = "conidiation"
    STRESS_RESPONSE = "stress response"
    GROWTH_DEVELOPMENT = "growth/development"
    DISPENSABLE = "dispensable / wild-type-like"
    NOT_REPORTED = "not reported"


class ExperimentalMethod(str, Enum):
    """Controlled vocabulary for experimental evidence types."""
    KNOCKOUT_DELETION = "knockout/deletion"
    COMPLEMENTATION = "complementation"
    OVEREXPRESSION = "overexpression"
    CRISPR_CAS9 = "CRISPR/Cas9"
    RNAI_SILENCING = "RNAi silencing"
    BIOCHEMICAL_ASSAY = "biochemical assay"
    COMPUTATIONAL_ONLY = "computational only"


class EvidenceConfidence(str, Enum):
    """Evidence quality tiers based on experimental rigor."""
    TIER_1 = "tier_1_high"       # Deletion + genetic complementation
    TIER_2 = "tier_2_high"       # Targeted deletion with validated mutant
    TIER_3 = "tier_3_medium"     # Overexpression, silencing, or chemical assay
    TIER_4 = "tier_4_low"        # Expression correlation or computational only


class EvidenceSource(str, Enum):
    """Provenance of the text an extraction was grounded in."""
    JATS_XML = "jats_xml"          # Structured full text from Europe PMC
    PDF = "pdf"                    # Layout-aware PDF -> Markdown conversion
    ABSTRACT_ONLY = "abstract_only"  # Full text unavailable (e.g. paywalled)


class AbstractTriageResult(BaseModel):
    """Result of first-pass abstract screening."""
    pmid: str = Field(description="PubMed ID")
    is_experimental: bool = Field(description="True if paper describes wet-lab functional genetics in M. oryzae")
    target_organism_confirmed: bool = Field(description="True if experiment performed in Magnaporthe oryzae")
    mentioned_genes: List[str] = Field(default_factory=list, description="Gene IDs or standard names detected in abstract")
    suggested_phenotypes: List[PhenotypeCategory] = Field(default_factory=list, description="Preliminary phenotype categories observed")
    relevance_score: float = Field(ge=0.0, le=1.0, description="Confidence score for experimental relevance")
    rationale: str = Field(description="One-sentence rationale explaining the decision")


class GeneEvidenceRecord(BaseModel):
    """Structured evidence record extracted from literature for a specific gene."""
    gene_id: str = Field(description="Systematic gene ID in MGG_##### format")
    gene_name: Optional[str] = Field(default=None, description="Common gene name (e.g. ENG1, PMK1)")
    process_affected: PhenotypeCategory = Field(description="Primary biological process affected")
    secondary_processes: List[PhenotypeCategory] = Field(default_factory=list, description="Additional affected processes for pleiotropic genes")
    detailed_phenotype: str = Field(description="Specific experimental observation (e.g. reduced appressoria, lost virulence on rice)")
    is_dispensable: bool = Field(default=False, description="True if knockout/mutant shows wild-type phenotype without defect (negative result)")
    experimental_method: ExperimentalMethod = Field(description="Genetic or experimental method used")
    complementation_tested: bool = Field(default=False, description="Whether wild-type phenotype was restored by complementation")
    confidence_tier: EvidenceConfidence = Field(description="Confidence tier of the evidence")
    verbatim_quote: str = Field(description="Exact verbatim excerpt from the paper text supporting the claim (mandatory for grounding)")
    evidence_source: EvidenceSource = Field(default=EvidenceSource.JATS_XML, description="Document type the evidence was extracted from (jats_xml, pdf, abstract_only)")
    pmid: Optional[str] = Field(default=None, description="PubMed ID")
    pmcid: Optional[str] = Field(default=None, description="PubMed Central ID if open access")
    doi: Optional[str] = Field(default=None, description="Digital Object Identifier")
    paper_title: Optional[str] = Field(default=None, description="Title of the paper")
    publication_year: Optional[int] = Field(default=None, description="Year of publication")
