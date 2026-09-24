"""Small dependency-free types used at scientific component boundaries."""

from __future__ import annotations

from enum import Enum
from typing import NewType, TypeAlias


CandidateId = NewType("CandidateId", str)

JsonScalar: TypeAlias = str | int | float | bool | None
JsonValue: TypeAlias = JsonScalar | list["JsonValue"] | dict[str, "JsonValue"]


class StringEnum(str, Enum):
    """Enum whose values serialize to the established JSON strings."""

    def __str__(self) -> str:
        return self.value


class EvidenceLevel(StringEnum):
    """Highest evidence actually supporting a reported scientific value."""

    ENUMERATED = "enumerated"
    SURROGATE = "surrogate"
    DFT = "dft"
    MULTIPHYSICS = "multiphysics"
    EXPERIMENTAL = "experimental"


class FidelityLevel(StringEnum):
    """Computational fidelity used to produce a result."""

    SCREENING = "screening"
    REDUCED_PHYSICS = "reduced_physics"
    FULL_PHYSICS = "full_physics"
    EXPERIMENTAL = "experimental"


class ConvergenceStatus(StringEnum):
    """Explicit convergence state; absence is never interpreted as success."""

    NOT_RUN = "not_run"
    PENDING = "pending"
    CONVERGED = "converged"
    NOT_CONVERGED = "not_converged"
    FAILED = "failed"
