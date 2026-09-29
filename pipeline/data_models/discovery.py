"""Named records exchanged by discovery selection and reactor evaluation."""

from __future__ import annotations

from typing import Literal, TypedDict


class _AdmissibilityOptional(TypedDict, total=False):
    reason: str
    admissible_count: int
    dropped_count: int


class AdmissibilitySummary(_AdmissibilityOptional):
    """Selection-only phase-stability filter outcome from ADR 0001."""

    filter: Literal['phase_stable_at_application_T'] | None


class EquilibriumSummary(TypedDict):
    """Compact comparison between reduced kinetics and equilibrium bounds."""

    within_tolerance: bool | None
    worst_abs_error: float | None
