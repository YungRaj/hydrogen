"""Readable data shapes shared by mode-specific reactor implementations."""

from __future__ import annotations

from typing import Literal, TypedDict

from pipeline.reactors.modes import ReactorTypeName
from pipeline.data_models.artifacts import ArtifactValidationResult


ReactorStatus = Literal[
    'complete', 'partial', 'failed', 'validation_required', 'not_applicable'
]


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
    multiphysics_evidence: ArtifactValidationResult | None


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


SolidsMetricRow = TypedDict(
    'SolidsMetricRow',
    {
        'catalyst_name': str | None,
        'reactor_type': ReactorTypeName | None,
        'T_K': float | None,
        'single_pass_CH4_conversion': float,
        'active_sv_1_m': float | None,
        'WHSV_h-1': float | None,
        'ergun_delta_p_Pa': float | None,
        'ergun_delta_p_bar': float | None,
        'ergun_ok': bool | None,
        'catalyst_E_act_eV': float | None,
        'catalyst_dE_H_eV': float | None,
        'h_parked': bool,
        'catalyst_particle_mm': float | None,
        'metal_loading': float | None,
        'metal_dispersion': float | None,
        'exceeds_equilibrium': bool,
        'carbon_balance_ok': bool | None,
    },
    total=False,
)
"""Unit-bearing reactor fields used to rank solid catalysts."""


class SolidsScorecard(TypedDict):
    """Fail-closed comparison of PFR/fluidized solids and MMBCR output."""

    judge_catalyst: str | None
    headline_catalyst: str | None
    judge_catalyst_requested: str | None
    judge_reason: str
    headline_t_min: float
    headline_t_max: float | None
    headline_band_note: str
    headline: dict[str, SolidsMetricRow]
    headline_solids_conversion: float | None
    solids_max_excluding_h_parked: SolidsMetricRow | None
    h_parked_excluded: list[SolidsMetricRow]
    mmbcr_max_conversion: float | None
    mmbcr_note: str
    n_solids_records: int
    n_mmbcr_records: int
