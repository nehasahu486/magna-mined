"""
Magna-Mined: Scalable Literature Curation Engine for Magnaporthe oryzae
"""

__version__ = "0.1.0"

from .schema import (
    PhenotypeCategory,
    ExperimentalMethod,
    EvidenceConfidence,
    GeneEvidenceRecord,
    AbstractTriageResult,
)
from .database import DatabaseManager
from .literature import EuropePMCClient
from .extractor import EvidenceExtractor

__all__ = [
    "PhenotypeCategory",
    "ExperimentalMethod",
    "EvidenceConfidence",
    "GeneEvidenceRecord",
    "AbstractTriageResult",
    "DatabaseManager",
    "EuropePMCClient",
    "EvidenceExtractor",
]
