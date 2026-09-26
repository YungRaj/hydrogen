"""Typed records for comparing solver outputs with published references."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping


SolverFamily = Literal["quantum_espresso", "cudaq_vqe"]
ReferenceKind = Literal["experiment", "published_computation"]


@dataclass(frozen=True)
class ReferenceObservation:
    """A literature value and the exact protocol under which it is comparable."""

    reference_id: str
    material_id: str
    observable: str
    value: float
    unit: str
    tolerance: float
    reference_kind: ReferenceKind
    solver_family: SolverFamily
    citation: str
    source_url: str
    conditions: Mapping[str, object]
    protocol: str
    tolerance_basis: str


@dataclass(frozen=True)
class ComputedObservation:
    """A physical observable derived from a completed solver calculation."""

    reference_id: str
    material_id: str
    observable: str
    value: float
    unit: str
    solver_family: SolverFamily
    protocol: str
    conditions: Mapping[str, object]
    calculation_conditions: Mapping[str, object]
    converged: bool
    mock: bool
    candidate_specific: bool
    benchmarked: bool
    artifact_sha256: str


@dataclass(frozen=True)
class ReferenceComparison:
    """Result of a like-for-like solver-to-literature comparison."""

    reference_id: str
    absolute_error: float
    tolerance: float
    passed: bool
