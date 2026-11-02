"""Readable data shapes shared by mode-specific reactor implementations."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Protocol, TypedDict

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt

from pipeline.reactors.modes import ReactorTypeName
from pipeline.data_models.artifacts import ArtifactValidationResult


ReactorStatus = Literal[
    'complete', 'partial', 'failed', 'validation_required', 'not_applicable'
]


class CanteraGasPhase(Protocol):
    """Typed protocol for a Cantera gas-phase thermo/kinetics object."""

    T: float
    P: float
    TP: tuple[float, float]
    TPX: object
    X: npt.NDArray[np.float64]
    species_names: list[str]

    def species_index(self, name: str) -> int:
        """Return the zero-based index of a gas species by name.

        Args:
            name: Species identifier in the Cantera gas phase.

        Returns:
            Zero-based species index in `species_names` and `X`.
        """
        ...


class CanteraSurfacePhase(Protocol):
    """Typed protocol for a Cantera catalytic interface phase."""

    coverages: npt.NDArray[np.float64]
    species_names: list[str]
    n_species: int

    def species_index(self, name: str) -> int:
        """Return the zero-based index of a surface species by name.

        Args:
            name: Species identifier on the Cantera catalytic interface.

        Returns:
            Zero-based species index in `species_names` and `coverages`.
        """
        ...


class LoadStatusFields(TypedDict):
    """Phase-loading status flags returned by load_status_fields."""

    gas_only: bool
    surface_loaded: bool
    graphite_loaded: bool
    surface_name: str


class KineticsProvenanceFields(TypedDict):
    """Candidate kinetic descriptor fields recorded on reactor outputs."""

    catalyst_E_act_eV: float
    catalyst_dE_H_eV: float


class KineticsEvidenceFields(TypedDict):
    """Kinetic provenance and exclusion-authority fields for a reactor run."""

    kinetics_status: str
    carbon_phase_model: str
    reactor_evidence_tier: str
    can_exclude_candidate: bool
    reactor_evidence_limitations: list[str]


class SolidsInventoryFields(TypedDict):
    """Catalyst surface-to-volume and Ergun pressure-drop summary fields."""

    geometric_sv_1_m: float
    active_sv_1_m: float
    active_area_multiplier: float
    ergun_delta_p_Pa: float
    ergun_delta_p_bar: float
    ergun_ok: bool


class PBRPelletEffectivenessMetrics(TypedDict):
    """Dimensionless and unit-bearing pellet mass-transfer diagnostics for PBR."""

    thiele_modulus: float
    internal_effectiveness_factor: float
    overall_effectiveness_factor: float
    external_damkohler_number: float
    external_mass_transfer_coeff_m_s: float
    effective_diffusivity_m2_s: float
    knudsen_diffusivity_m2_s: float
    molecular_diffusivity_m2_s: float
    sherwood_number: float
    particle_reynolds_number: float


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
    exit_P_Pa: float
    coupled_ergun_delta_p_Pa: float
    coupled_ergun_delta_p_bar: float
    catalyst_weight_kg: float
    stage_catalyst_weight_kg: float
    WHSV_mass_kg_feed_per_kg_cat_h: float
    thiele_modulus: float
    internal_effectiveness_factor: float
    overall_effectiveness_factor: float
    mean_overall_effectiveness_factor: float
    external_damkohler_number: float
    external_mass_transfer_coeff_m_s: float
    effective_diffusivity_m2_s: float
    knudsen_diffusivity_m2_s: float
    molecular_diffusivity_m2_s: float
    sherwood_number: float
    particle_reynolds_number: float
    w_positions_kg: list[float]
    pressure_profile_Pa: list[float]
    bed_voidage_profile: list[float]
    exit_bed_porosity: float
    thiele_modulus_profile: list[float]
    effectiveness_factor_profile: list[float]
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


class CarbonPolicyConfig(TypedDict, total=False):
    """Explicit carbon-removal and regeneration policy applied to a reactor sweep."""

    co2_permitted: bool
    fluidized_mode: Literal['batch_regen', 'circulating']
    max_regen_cycles: int
    regen_mechanism: Literal['mechanical', 'consumable', 'oxidative']


class AnalyticalClosureEvidence(TypedDict):
    """Provenance record for an analytical hydrodynamic closure fallback."""

    source: str
    candidate_exclusion_authorized: bool
    basis: str
    multiphysics_evidence: ArtifactValidationResult | dict[str, object]


class ResolvedKineticsBarriers(TypedDict):
    """Elementary activation barriers in eV consumed by CandidateKinetics."""

    methane_activation_eV: float
    ch3_dehydrogenation_eV: float
    ch2_dehydrogenation_eV: float
    ch_dehydrogenation_eV: float
    h2_desorption_eV: float
    carbon_transfer_eV: float
    carbon_encapsulation_eV: float


class MlipNebValidationRecord(TypedDict):
    """Elementary activation barriers and convergence count from 5-step MLIP CI-NEB."""

    E_act_neb: float
    ch3_dehydrogenation_eV: float
    ch2_dehydrogenation_eV: float
    ch_dehydrogenation_eV: float
    h2_desorption_eV: float
    carbon_transfer_eV: float
    carbon_encapsulation_eV: float
    mlip_neb_steps_converged: int


class _PyrolysisValidationOptional(TypedDict, total=False):
    evidence_level: str
    mlip_neb_steps_converged: int


class PyrolysisValidationArtifact(_PyrolysisValidationOptional):
    """Persisted pyrolysis_validation.json consumed by Phase 2 Cantera mechanisms."""

    candidate_id: str
    complete: bool
    resolved_kinetics_eV: ResolvedKineticsBarriers


class _SolidsScorecardOptional(TypedDict, total=False):
    candidate_champion_conversion: float | None
    candidate_champion_catalyst: str | None
    candidate_champion_reactor: str | None
    candidate_champion_T_K: float | None
    candidate_champion_kinetics_status: str | None


class SolidsScorecard(_SolidsScorecardOptional):
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
    ranking_authority: bool
    solids_max_excluding_h_parked: float | None
    h_parked_excluded: list[SolidsMetricRow]
    mmbcr_max_conversion: float | None
    mmbcr_note: str
    ntec_max_conversion: float | None
    electrochemical_max_conversion: float | None
    n_solids_records: int
    n_mmbcr_records: int
    n_ntec_records: int
    n_electrochemical_records: int


class B5ClosureGateResult(TypedDict):
    """Pass/fail evidence summary for the B5 solids ranking authority gate."""

    b5_accepted: bool
    judge_ok: bool
    pfr_roi_conversion: float
    pbr_roi_conversion: float
    fluidized_roi_conversion: float
    pfr_conversion_ok: bool
    pbr_conversion_ok: bool
    fluidized_conversion_ok: bool
    b66_criteria_ok: bool
    b67_joint_roi_ok: bool

