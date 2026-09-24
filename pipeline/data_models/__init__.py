"""Canonical data contracts shared across independently replaceable stages.

The package deliberately contains no solver imports.  It describes values
crossing component boundaries while each scientific implementation retains
ownership of its internal data structures.
"""

from pipeline.data_models.core import (
    CandidateId,
    ConvergenceStatus,
    EvidenceLevel,
    FidelityLevel,
    JsonValue,
)

__all__ = [
    "CandidateId",
    "ConvergenceStatus",
    "EvidenceLevel",
    "FidelityLevel",
    "JsonValue",
]
