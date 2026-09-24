"""Readable data shapes shared by mode-specific reactor implementations."""

from __future__ import annotations

from typing import Any, Literal, TypedDict

from pipeline.reactors.modes import ReactorTypeName


ReactorStatus = Literal[
    'complete', 'partial', 'failed', 'validation_required', 'not_applicable']


class ReactorResult(TypedDict, total=False):
    """Common superset carried by one reactor condition result.

    Mode-specific solvers may add fields, but shared consumers use only this
    documented set. Unit-bearing names prevent ambiguous numerical handoffs.
    """

    status: ReactorStatus
    valid: bool
    reactor_type: ReactorTypeName
    pathway_mode: str
    catalyst_name: str
    candidate_id: str
    material_class: str
    T_K: float
    P_Pa: float
    CH4_conversion: float
    single_pass_CH4_conversion: float
    H2_selectivity: float
    solid_C_selectivity: float | None
    residence_time_s: float
    specific_energy_kWh_kg_H2: float
    faradaic_efficiency_H2: float
    current_density_A_cm2: float
    cell_voltage_V: float
    electrical_power_density_W_cm2: float
    conversion_basis: str
    thermal_mode: str
    reactor_evidence_tier: str
    converged: bool
    mock: bool
    can_exclude_candidate: bool
    limitations: list[str]
    multiphysics_evidence: dict[str, Any]


class ReactorSweepSummary(TypedDict):
    """Evidence-aware summary of all conditions for one candidate."""

    best_condition: ReactorResult
    sweep_status: ReactorStatus
    completed_conditions: int
    usable_conditions: int
    failed_conditions: int
    pending_conditions: int
    not_applicable_conditions: int
    unresolved_conditions: int
    can_exclude_candidate: bool


class CandidateReactorResult(ReactorSweepSummary):
    """Candidate identity, mechanism provenance, and its complete sweep."""

    catalyst: str
    candidate_id: str
    material_class: str | None
    pathway_mode: str
    E_act: float
    mechanism_file: str | None
    sweep: list[ReactorResult]
