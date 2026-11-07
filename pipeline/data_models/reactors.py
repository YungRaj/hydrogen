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

FeedstockName = Literal[
    'CH4',
    'ASSOC_GAS',
    'C2H6',
    'C3H8',
    'C4H10',
    'C5H12',
    'C6H14',
    'C6H12',
    'C8H18',
    'C12H26',
    'C16H34',
    'C10H20_PE',
    'BIOMETHANE',
    'BIOGAS',
    'C18H34O2_LIPID',
    'C3H8O3_GLYCEROL',
    'C6H10O2_BIOOIL',
    'C6H10O4_BIOMASS',
]

FeedstockCategory = Literal[
    'natural_gas',
    'light_alkane',
    'liquid_hydrocarbon',
    'cycloalkane_olefin',
    'waste_plastic',
    'biomass_biogenic',
]


class H2C2SelectivityRateMetrics(TypedDict):
    """Explicit selectivity, yield, C2 mole-fraction, and H2 rate metrics.

    Note:
        Both feedstock-general keys (``exit_x_feedstock``,
        ``feedstock_conversion_rate_mol_m3_s``) and legacy CH4-default keys
        (``exit_x_CH4``, ``ch4_conversion_rate_mol_m3_s``) are populated so
        callers can inspect either the general feedstock name or the historical
        methane-default schema key.
    """

    H2_selectivity: float
    H2_atom_balance: float
    H2_yield: float
    solid_C_selectivity: float
    C2_carbon_selectivity: float
    exit_x_feedstock: float
    exit_x_CH4: float
    exit_x_C2H6: float
    exit_x_C2H4: float
    exit_x_C2H2: float
    feedstock_conversion_rate_mol_m3_s: float
    ch4_conversion_rate_mol_m3_s: float
    h2_production_rate_mol_m3_s: float
    h2_production_rate_kg_m3_h: float
    h2_turnover_frequency_1_s: float | None


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


class CanteraSpeciesDescriptor(Protocol):
    """Typed protocol for a Cantera species descriptor returned by surf.species(k)."""

    size: float
    composition: dict[str, float]


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

    def species(self, index: int | str) -> CanteraSpeciesDescriptor:
        """Return the species descriptor for a zero-based index or species name.

        Args:
            index: Zero-based species index or species identifier.

        Returns:
            Cantera species descriptor exposing `size` and `composition`.
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
    When ``feedstock`` is a non-methane hydrocarbon or biomass feedstock,
    feedstock-general fields (``feedstock_conversion``,
    ``single_pass_feedstock_conversion``, ``cycle_averaged_feedstock_conversion``,
    ``exit_x_feedstock``, ``feedstock_conversion_rate_mol_m3_s``,
    ``endothermic_heat_duty_kJ_per_mol_feed_fed``,
    ``endothermic_heat_duty_T_kJ_per_mol_feed_fed``) carry the active feedstock
    values, and legacy ``CH4_*`` keys are synchronized to those values for
    backward compatibility with CH4-default consumers.
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
    # Active feedstock conversion (synchronized with CH4_* keys for backward compatibility)
    feedstock_conversion: float
    single_pass_feedstock_conversion: float
    cycle_averaged_feedstock_conversion: float
    CH4_conversion: float
    single_pass_CH4_conversion: float
    cycle_averaged_CH4_conversion: float
    X_eq_local_tp: float
    X_eq_table: float
    gas_holdup_fraction: float
    bed_porosity: float
    active_sv_1_m: float
    interfacial_sv_ratio_1_m: float
    site_density_mol_cm2: float
    catalyst_E_act_eV: float
    catalyst_dE_H_eV: float
    H2_selectivity: float
    H2_atom_balance: float
    H2_yield: float
    solid_C_selectivity: float | None
    C2_carbon_selectivity: float
    exit_x_feedstock: float
    exit_x_CH4: float
    exit_x_H2: float
    exit_x_C2H6: float | None
    exit_x_C2H4: float | None
    exit_x_C2H2: float | None
    feedstock_conversion_rate_mol_m3_s: float
    ch4_conversion_rate_mol_m3_s: float
    h2_production_rate_mol_m3_s: float
    h2_production_rate_kg_m3_h: float
    h2_turnover_frequency_1_s: float | None
    feedstock: FeedstockName
    feedstock_name: str
    feedstock_category: FeedstockCategory
    feedstock_n_C: int
    feedstock_n_H: int
    feedstock_n_O: int
    feedstock_h_to_c_ratio: float
    feedstock_nu_H2: float
    feedstock_nu_C: float
    feedstock_nu_CO: float
    feedstock_kg_per_kg_h2: float
    feedstock_kg_c_per_kg_h2: float
    exit_x_CO: float
    catalyst_E_act_CH4_reference_eV: float
    endothermic_heat_duty_kJ_per_mol_feed_fed: float
    endothermic_heat_duty_T_kJ_per_mol_feed_fed: float
    endothermic_heat_duty_kJ_per_mol_CH4_fed: float
    endothermic_heat_duty_T_kJ_per_mol_CH4_fed: float
    c2_tracked: bool
    residence_time_s: float
    specific_energy_kWh_kg_H2: float
    faradaic_efficiency_H2: float
    current_density_A_cm2: float
    cell_voltage_V: float
    electrical_power_density_W_cm2: float
    conversion_basis: str
    thermal_mode: str
    reactor_evidence_tier: str
    reactor_evidence_limitations: list[str]
    converged: bool
    mock: bool
    can_exclude_candidate: bool
    limitations: list[str]
    multiphysics_evidence: ArtifactValidationResult | None
    reactor_closure_evidence: AnalyticalClosureEvidence | None


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
        'literature_catalyst_name': str | None,
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
    feedstock: str
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


class PBRErgunDiagnosticsSummary(TypedDict):
    """Summary of PBR intraparticle effectiveness, Ergun ΔP, pilot residence-time, and regenerated-cycle conversion deltas."""

    n_completed: int
    n_pbr_records: int
    mean_pbr_overall_effectiveness_factor: float | None
    mean_pbr_thiele_modulus: float | None
    n_ergun_checked: int
    n_ergun_ok: int
    n_pilot_residence_checked: int
    n_within_pilot_residence_band: int
    n_plugging_or_overloaded_flagged: int
    mean_regenerated_cycle_minus_single_pass_X: float | None


class InventoryGridCell(TypedDict):
    """One quantized (particle diameter, metal loading, metal dispersion) solids inventory cell."""

    catalyst_particle_mm: float
    metal_loading: float
    metal_dispersion: float


class MMBCRHydrodynamicsResult(TypedDict):
    """Mendelson bubble rise, gas holdup, residence time, and melt property diagnostics for MMBCR."""

    bubble_rise_velocity_m_s: float
    bubble_rise_basis: str
    gas_holdup_fraction: float
    gas_holdup_basis: str | None
    residence_time_s: float
    residence_time_basis: str
    interfacial_sv_ratio_1_m: float
    melt_surface_tension_N_m: float
    melt_density_kg_m3: float
    melt_property_basis: str


class KuniiLevenspielExchangeMetrics(TypedDict):
    """Kunii-Levenspiel (1991) bubble-to-emulsion mass-transfer metrics for fluidized beds."""

    umf_wen_yu_m_s: float
    bubble_diameter_kl_m: float
    bubble_rise_velocity_kl_m_s: float
    k_be_1_s: float
    interphase_transfer_units_N_be: float
    bubble_exchange_efficiency: float
    effective_bypass_fraction: float


class SupportedNiFilamentClosureGateResult(TypedDict):
    """Pass/fail evidence summary for the supported-Ni(111)/SiO2 filament closure & solids ranking authority gate."""

    solids_ranking_authority_accepted: bool
    judge_ok: bool
    pfr_roi_conversion: float
    pbr_roi_conversion: float
    fluidized_roi_conversion: float
    pfr_conversion_ok: bool
    pbr_conversion_ok: bool
    fluidized_conversion_ok: bool
    baker_helveg_carbon_criteria_ok: bool
    ermakova_takenaka_joint_roi_ok: bool


class _SolidsScorecardOptional(TypedDict, total=False):
    judge_catalyst_literature_name: str | None
    headline_catalyst_literature_name: str | None
    candidate_champion_conversion: float | None
    candidate_champion_catalyst: str | None
    candidate_champion_literature_name: str | None
    candidate_champion_reactor: str | None
    candidate_champion_T_K: float | None
    candidate_champion_kinetics_status: str | None
    pbr_ergun_residence_diagnostics: PBRErgunDiagnosticsSummary
    supported_ni_filament_closure_gate: SupportedNiFilamentClosureGateResult | None
    validation_gate_registry: dict[str, dict[str, object]]


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


class _SupportedNiFilamentClosureSummaryOptional(TypedDict, total=False):
    content_sha256: str


class SupportedNiFilamentClosureSummary(_SupportedNiFilamentClosureSummaryOptional):
    """Full Alves-Bengaard-Ermakova-Takenaka supported-Ni(111)/SiO2 filament closure artifact."""

    schema_version: str
    milestone: str
    judge_catalyst: str
    judge_catalyst_literature_name: str
    roi_temperature_band_K: list[float]
    supported_ni_filament_closure_gate: SupportedNiFilamentClosureGateResult
    ranking_authority: bool
    solids_max_excluding_h_parked: float | None
    baker_helveg_criteria_summary: dict[str, object]
    ermakova_takenaka_joint_summary: dict[str, object]
    calibrated_solids_scorecard: SolidsScorecard
    validation_gate_registry: dict[str, dict[str, object]]



