"""Named evidence, selection, campaign-status, and lineage documents."""

from __future__ import annotations

from typing import Literal, TypedDict

from pipeline.data_models.core import JsonValue
from pipeline.data_models.reactors import (
    FeedstockName,
    PBRErgunDiagnosticsSummary,
    SupportedNiFilamentClosureSummary,
)


CandidateDisposition = Literal[
    'hard_excluded', 'validation_required', 'quantitative_screening'
]


class CampaignStatus(TypedDict):
    """Fail-closed summary of campaign acceptance criteria."""

    ready: bool
    criteria: dict[str, bool]
    missing: list[str]
    evidence_manifest_valid: bool
    evidence_errors: list[str]


class LedgerEvent(TypedDict):
    """One hash-linked, JSON-compatible multi-fidelity campaign event."""

    sequence: int
    event_type: str
    previous_event_sha256: str | None
    payload: dict[str, JsonValue]
    event_sha256: str


class LedgerDocument(TypedDict):
    """Versioned persisted campaign event chain."""

    schema_version: int
    events: list[LedgerEvent]


class TEAAssumptionsDict(TypedDict):
    """Serialized techno-economic assumption fields."""

    natural_gas_usd_mmbtu: float
    electricity_usd_kwh: float
    base_energy_kwh_kg_h2: float
    capex_usd_kg_h2: float
    carbon_value_usd_kg_h2: float
    source_id: str


class _TEAScenarioEstimateOptional(TypedDict, total=False):
    feedstock: FeedstockName
    feedstock_name: str
    kg_feed_per_kg_h2: float
    kg_c_per_kg_h2: float
    net_feedstock_utilization: float
    psa_feedstock_recovery: float


class TEAScenarioEstimate(_TEAScenarioEstimateOptional):
    """Single-scenario levelized hydrogen cost estimate with PSA hydrocarbon recycle.

    Note:
        ``net_feedstock_utilization`` and ``psa_feedstock_recovery`` carry the
        feedstock-general utilization and PSA recovery fraction; legacy keys
        ``net_ch4_utilization`` and ``psa_ch4_recovery`` are populated with the
        same values for backward compatibility with CH4-default callers.
    """

    scenario: str
    h2_cost_usd_kg: float
    conversion: float
    recycle_ratio: float
    net_ch4_utilization: float
    psa_ch4_recovery: float
    solid_c_selectivity: float
    carbon_credit_usd_kg_h2: float
    energy_input_kwh_kg_h2: float
    effective_capex_usd_kg_h2: float
    assumptions: TEAAssumptionsDict
    evidence_level: Literal['screening_scenario_not_measured_tea']


class TEAScenarioRangeResult(TypedDict):
    """Sensitivity range across optimistic, base, and conservative TEA scenarios."""

    estimates: dict[str, TEAScenarioEstimate]
    min_usd_kg: float
    max_usd_kg: float
    evidence_level: Literal['screening_sensitivity_range']


class _DiscoveryChampionRowOptional(TypedDict, total=False):
    qe_dft_converged: bool
    qe_dft_energy_Ry: float | None
    qe_dft_energy_eV: float | None
    qe_max_force_Ry_bohr: float | None
    qe_neb_converged: bool
    qe_neb_barrier_eV: float


class DiscoveryChampionRow(_DiscoveryChampionRowOptional):
    """Phase 1 turquoise-hydrogen discovery champion annotated with interfacial and QE DFT evidence."""

    rank: int
    candidate_id: str
    literature_catalyst_name: str
    material_class: str
    genome: str
    E_act_eV: float
    calibrated_activation_eV: float | None
    interfacial_carbon_displacement_A: float | None
    interfacial_pyrolysis_score_eV: float | None
    calibrated_interfacial_score_eV: float | None
    coking_index_eV: float | None


class PathwayLeaderboardRow(TypedDict, total=False):
    """Champion reactor and TEA summary for one pathway mode."""

    status: str
    n_usable_records: int
    champion_catalyst: str
    champion_catalyst_literature_name: str
    reactor_type: str
    T_K: float
    single_pass_CH4_conversion: float
    cycle_averaged_CH4_conversion: float | None
    overall_effectiveness_factor: float | None
    thiele_modulus: float | None
    ergun_delta_p_bar: float | None
    H2_selectivity: float
    H2_yield: float | None
    solid_C_selectivity: float | None
    C2_carbon_selectivity: float | None
    feedstock_conversion_rate_mol_m3_s: float | None
    ch4_conversion_rate_mol_m3_s: float | None
    h2_production_rate_mol_m3_s: float | None
    h2_production_rate_kg_m3_h: float | None
    h2_turnover_frequency_1_s: float | None
    reactor_evidence_tier: str
    kinetics_status: str
    tea_base_h2_usd_kg: float | None
    tea_optimistic_h2_usd_kg: float | None
    tea_conservative_h2_usd_kg: float | None
    energy_input_kwh_kg_h2_base: float | None


class DFTCalibratedReactorMetrics(TypedDict):
    """DFT/NEB-calibrated reactor conversion and effectiveness metrics."""

    material_class: str
    dft_barrier_shift_eV: float
    arrhenius_rate_scaling: float
    dft_calibrated_single_pass_CH4_conversion: float
    _x_cal_exact: float
    dft_calibrated_overall_effectiveness_factor: float | None
    dft_calibrated_thiele_modulus: float | None
    kinetics_status: str


class ExtractedSelectivityRateMetrics(TypedDict):
    """Selectivity, H2 yield, and volumetric/TOF rate metrics extracted for scorecard rows."""

    H2_selectivity: float
    H2_yield: float
    solid_C_selectivity: float
    C2_carbon_selectivity: float
    feedstock_conversion_rate_mol_m3_s: float
    ch4_conversion_rate_mol_m3_s: float
    h2_production_rate_mol_m3_s: float
    h2_production_rate_kg_m3_h: float
    h2_turnover_frequency_1_s: float | None


class DFTCalibratedPathwayRow(TypedDict, total=False):
    """QE-DFT/NEB-calibrated champion reactor and TEA summary for one pathway mode."""

    status: str
    n_usable_records: int
    champion_catalyst: str
    champion_catalyst_literature_name: str
    material_class: str
    reactor_type: str
    T_K: float
    uncalibrated_single_pass_CH4_conversion: float
    dft_barrier_shift_eV: float
    arrhenius_rate_scaling: float
    dft_calibrated_single_pass_CH4_conversion: float
    dft_calibrated_overall_effectiveness_factor: float | None
    dft_calibrated_thiele_modulus: float | None
    ergun_delta_p_bar: float | None
    H2_selectivity: float
    H2_yield: float | None
    solid_C_selectivity: float | None
    C2_carbon_selectivity: float | None
    dft_calibrated_feedstock_conversion_rate_mol_m3_s: float | None
    dft_calibrated_ch4_conversion_rate_mol_m3_s: float | None
    dft_calibrated_h2_production_rate_mol_m3_s: float | None
    dft_calibrated_h2_production_rate_kg_m3_h: float | None
    dft_calibrated_h2_turnover_frequency_1_s: float | None
    kinetics_status: str
    dft_calibrated_tea_base_h2_usd_kg: float | None
    dft_calibrated_tea_optimistic_h2_usd_kg: float | None
    dft_calibrated_tea_conservative_h2_usd_kg: float | None
    dft_calibrated_energy_input_kwh_kg_h2_base: float | None


class PEMFCCategoryChampionSummary(TypedDict):
    """Category-level (PGM-free or low-PGM) acid-stable PEMFC cathode champion."""

    cathode_catalyst: str
    cathode_catalyst_literature_name: str
    material_class: str
    genome: str
    orr_overpotential_gc_V: float
    pourbaix_dissolution_margin_V: float
    peak_power_W_cm2: float
    cell_efficiency_at_peak: float


class PEMFCStackSummaryRow(TypedDict, total=False):
    """Screening or QE-GC-confirmed PEMFC cathode and 400-cell stack summary."""

    status: str
    n_single_cell_simulations: int
    n_qe_converged_acid_stable_cathodes: int
    champion_cathode: str
    champion_cathode_literature_name: str
    champion_material_class: str
    champion_genome: str
    champion_membrane: str
    orr_overpotential_V: float
    orr_overpotential_che_V: float
    orr_overpotential_gc_V: float
    orr_tafel_slope_mV_dec: float
    pourbaix_dissolution_margin_V: float
    qe_dft_energy_Ry: float | None
    peak_power_W_cm2: float
    rated_power_at_0p6V_W_cm2: float | None
    cell_efficiency_at_peak: float
    stack_cells: int
    stack_net_power_kW: float
    stack_system_efficiency: float
    stack_gravimetric_W_kg: float
    stack_volumetric_W_L: float
    stack_cost_usd_per_kW: float
    pgm_free_champion: PEMFCCategoryChampionSummary | None
    low_pgm_champion: PEMFCCategoryChampionSummary | None


class CoupledPathwayTEARow(TypedDict, total=False):
    """Coupled turquoise-H2 to PEMFC electricity fuel-cost summary for one pathway."""

    h2_cost_base_usd_kg: float
    h2_cost_optimistic_usd_kg: float
    pemfc_fuel_cost_base_usd_kwh_e: float
    pemfc_fuel_cost_optimistic_usd_kwh_e: float
    dft_calibrated_h2_cost_base_usd_kg: float
    dft_calibrated_h2_cost_optimistic_usd_kg: float
    dft_calibrated_h2_cost_conservative_usd_kg: float | None
    dft_calibrated_pemfc_fuel_cost_base_usd_kwh_e: float
    dft_calibrated_pemfc_fuel_cost_optimistic_usd_kwh_e: float


class MultipathwayChampionScorecard(TypedDict):
    """End-to-end 7-pathway discovery, reactor, DFT, PEMFC, and TEA champion scorecard."""

    discovery_champions: list[DiscoveryChampionRow]
    pathway_leaderboard: dict[str, PathwayLeaderboardRow]
    dft_calibrated_pathway_leaderboard: dict[str, DFTCalibratedPathwayRow]
    pbr_ergun_residence_diagnostics: PBRErgunDiagnosticsSummary
    supported_ni_filament_closure_summary: SupportedNiFilamentClosureSummary | None
    validation_gate_registry: dict[str, dict[str, object]]
    literature_validation_gates: dict[str, dict[str, object]]
    named_reference_catalysts: dict[str, dict[str, str]] | dict[str, dict[str, object]]
    pemfc_stack_summary: PEMFCStackSummaryRow
    dft_calibrated_pemfc_stack_summary: PEMFCStackSummaryRow
    coupled_turquoise_pemfc_tea: dict[str, CoupledPathwayTEARow]
    dft_calibrated_coupled_turquoise_pemfc_tea: dict[str, CoupledPathwayTEARow]
    evidence_level: str


class CalibrationBucketStats(TypedDict):
    """Empirical MAE, bias, productivity, and best observed value for a region or class."""

    n: int
    mae: float
    bias: float
    productivity: float
    best_observed: float


class CalibrationObservation(TypedDict):
    """Paired surrogate prediction and DFT observation recorded for active-learning calibration."""

    application: str
    candidate_id: str
    material_class: str
    region: str
    genome: str
    predicted_eV: float
    observed_eV: float
    residual_eV: float
    fidelity: str
    productive: bool
    provenance: dict[str, float | str | bool]


class ApplicationCalibrationSummary(TypedDict):
    """Empirical-Bayes class and regional bias maps for one scientific application."""

    n_validation_samples: int
    class_bias_eV: dict[str, float]
    class_bias_delta_eV: dict[str, float]
    region_bias_eV: dict[str, float]
    class_calibration_stats: dict[str, CalibrationBucketStats]
    regional_calibration_stats: dict[str, CalibrationBucketStats]


class ActiveLearningCalibrationArtifact(TypedDict):
    """Persisted active-learning QE DFT calibration document with SHA-256 provenance."""

    schema_version: int
    n_validation_samples: int
    prior_rows: float
    class_bias_eV: dict[str, float]
    class_bias_delta_eV: dict[str, float]
    region_bias_eV: dict[str, float]
    applications: dict[str, ApplicationCalibrationSummary]
    observations: list[CalibrationObservation]
    calibration_sha256: str


class Round2BatchReport(TypedDict):
    """Per-batch hit counts and stratified random baseline for Round-2 evaluation."""

    batch_id: str
    n_eligible: int
    n_valid: int
    hit_count: int
    policy_hits: dict[str, int]
    random_mean_hits: float


class Round2ApplicationReport(TypedDict, total=False):
    """Aggregated Round-2 prospective policy evaluation for one application."""

    n_batches: int
    total_budget_selected: int
    policy_hits_gc_acid_stable_top20pct: dict[str, int]
    policy_hits_interfacial_calibrated_top20pct: dict[str, int]
    doe_035V_gc_acid_stable_hits: dict[str, int]
    dft_calibrated_hit_lift_vs_uncalibrated: int
    dft_calibrated_doe_hit_lift_vs_uncalibrated: int
    policy_matched_random_mean_hits: float
    dft_calibrated_enrichment_vs_random: float
    dft_calibrated_p_value: float
    uncalibrated_p_value: float
    acceptance_passed: bool
    batches: list[Round2BatchReport]


class Round2CampaignReport(TypedDict):
    """Full Round-2 DFT-calibrated prospective evaluation report with SHA-256 digest."""

    schema_version: int
    calibration_sha256: str | None
    fuel_cell_orr: Round2ApplicationReport
    turquoise_hydrogen: Round2ApplicationReport
    report_sha256: str


class _ORRSlateCandidateRowBase(TypedDict):
    """Base fields for a candidate in the 50-catalyst ORR confirmatory slate."""

    batch_id: str
    candidate_id: str
    genome: str
    material_class: str
    predicted_overpotential_V: float
    dG_OH_eV: float
    dG_O_eV: float
    dG_OOH_eV: float
    model_confidence: float
    policy_membership: list[str]


class ORRSlateSelectedRow(_ORRSlateCandidateRowBase, total=False):
    """Selected ORR slate row with assigned selection role and optional rank."""

    selection_role: str
    rank: int


class ORRValidatedSlateRecord(ORRSlateSelectedRow):
    """QE-DFT and Grand-Canonical Pourbaix validated ORR slate record."""

    orr_overpotential_gc_V: float
    gc_field_shift_OH_eV: float
    pourbaix_dissolution_onset_V: float
    pourbaix_dissolution_margin_V: float
    pourbaix_acid_stable: bool
    qe_converged: bool
    qe_dft_energy_Ry: float | None
    qe_dft_energy_eV: float | None
    qe_max_force_Ry_bohr: float | None
    doe_target_hit: bool
    doe_target_hit_gc: bool


class ORRRoleSummary(TypedDict):
    """Per-role convergence, Pourbaix stability, and DOE target hit summary."""

    count: int
    n_qe_converged: int
    n_pourbaix_acid_stable: int
    mean_overpotential_V: float
    min_overpotential_V: float
    mean_overpotential_gc_V: float
    min_overpotential_gc_V: float
    doe_target_hits_le_0_35V: int
    doe_target_hits_gc_le_0_35V: int


class ORRPolicyEnrichmentSummary(TypedDict):
    """Guided vs. control DOE target hit-rate enrichment summary."""

    guided_and_champion_count: int
    control_count: int
    guided_doe_hit_rate: float
    control_doe_hit_rate: float
    enrichment_ratio: float | None


class ORRConfirmatoryDFTSummary(TypedDict):
    """Persisted 50-candidate ORR confirmatory DFT summary artifact."""

    schema_version: int
    application: str
    slate_sha256: str
    n_slate_candidates: int
    n_qe_converged: int
    qe_convergence_rate: float
    n_pourbaix_acid_stable: int
    role_summary: dict[str, ORRRoleSummary]
    policy_enrichment: ORRPolicyEnrichmentSummary
    records: list[ORRValidatedSlateRecord]
    summary_sha256: str

