#!/usr/bin/env python3
# Turquoise-hydrogen reactor process models.
"""
Cantera Reactor-Scale Simulation Models for Methane Pyrolysis.

Three Cantera-backed reactor archetypes with distinct carbon-handling physics:
  A. MMBCR — melt ODE to tabulated X_eq with bubble-area flotation. Residence
     time is column height over the Mendelson bubble rise velocity; gas holdup
     is derived, not an input. Carbon leaves the bubble; no site lattice.
  B. PFR — one shared Langmuir surface marched through stages (θ_C is
     time-on-stream, not axial); optional discrete non-oxidative regen.
  C. Fluidized — emulsion pass with in-step C_s removal (circulating) or
     batch regen, then bubble bypass mixed on molar flows via the Ar tracer.

NTEC and electrochemical pathways consume validated external multiphysics
artifacts. Cantera supplies chemistry within those coupled calculations, but
does not itself solve their mechanical/electrical/charge-transport physics.

Solid carbon is never a gas-phase species. Surface C_s blocks sites on solid
paths until removed by a named policy or a gated off-site channel (B6).
Oxidative regen requires co2_permitted. Γ is a monolayer and is locked (B1).

Reactor routing is owned by ``pipeline.process.pathway_modes``: a reactor
type must belong to the selected pathway mode, and MMBCR is applicable only
to MoltenMetal candidates.
"""

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

try:
    import cantera as ct
    HAS_CANTERA = True
except ImportError:
    HAS_CANTERA = False

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.common.utils import (
    REACTOR_DIR,
    setup_logger, print_banner, save_json,
)
from pipeline.process.pathway_modes import (
    DEFAULT_MODE, REACTOR_MODELS, reactor_applicability,
    reactor_types_for_mode, validate_mode_reactors)
from pipeline.process.reactor_mechanisms import MONOLAYER_SITE_DENSITY_MOL_CM2

logger = setup_logger('reactor_models', 'reactor/reactor_simulation.log')

REGEN_MECHANICAL = 'mechanical'
REGEN_CONSUMABLE = 'consumable'
REGEN_OXIDATIVE = 'oxidative'
ALLOWED_REGEN = frozenset({REGEN_MECHANICAL, REGEN_CONSUMABLE, REGEN_OXIDATIVE})
FLUIDIZED_BATCH = 'batch_regen'
FLUIDIZED_CIRCULATING = 'circulating'
ALLOWED_FLUIDIZED = frozenset({FLUIDIZED_BATCH, FLUIDIZED_CIRCULATING})

# Packed-bed ΔP flag (B4). Cells above this are still run and marked.
ERGUN_DP_LIMIT_PA = 101325.0

# Melt-side interfacial prefactor (B3). Order-of-magnitude bubble-interface
# velocity, not a DFT barrier and not a Upham-fitted mass-transfer k.
# Upham 2017; Chen 2023; Abdollahi 2024. The 1 m/s cap is a guardrail
# against inventing Da, not a physical bound like loading ≤ 1.
MMBCR_INTERFACIAL_K0_DEFAULT = 0.01
MMBCR_INTERFACIAL_K0_MAX = 1.0
MMBCR_INTERFACIAL_K0_MAX_KIND = 'guardrail_not_physical_bound'
# None = unconstrained flotation (carbon leaves as produced). 0 = fouled
# interface. Finite = k_float / (k_float + k_if * a) on the ODE.
MMBCR_FLOTATION_UNCONSTRAINED = None
# Melt properties for the Mendelson bubble rise velocity. Ni-Bi order of
# magnitude (σ ~ 0.4 N/m, ρ ~ 9000 kg/m³); not fitted to any column.
MMBCR_MELT_SURFACE_TENSION_N_M = 0.4
MMBCR_MELT_DENSITY_KG_M3 = 9000.0
# Above this derived holdup the column is churn-turbulent and the bubbly-flow
# CSTR-cascade / single-bubble ODE picture no longer holds. Fail closed.
MMBCR_MAX_GAS_HOLDUP = 0.3
MMBCR_BUBBLE_RISE_BASIS = 'mendelson_sqrt(2*sigma/(rho*d_b) + g*d_b/2)'
MMBCR_RESIDENCE_BASIS = 'column_height_over_bubble_rise_velocity'
MMBCR_ARTIFACT_RESIDENCE_BASIS = 'gas_holdup_times_column_height_over_superficial_velocity'
G_M_S2 = 9.80665
R_J_MOL_K = 8.314462618  # ct.gas_constant is per kmol; do not mix
FLUIDIZED_REMOVAL_SUBSTEPS = 20
# Emulsion voidage at minimum fluidization (Geldart B order of magnitude).
# Emulsion area is per emulsion volume; the reacting gas volume in 1 m³ of
# emulsion is ε_mf. Same area basis as PFR (area per bed volume × bed volume).
FLUIDIZED_EMULSION_VOIDAGE = 0.45

# Closure tiers for the two Cantera reactors whose hydrodynamics upstream
# expects from an external solver. A validated OpenFOAM artifact or a
# calibrated surrogate always wins; absent both, these reactors run on an
# explicitly labelled analytical bubbly-flow closure that can never exclude
# a candidate. NTEC / Electrochemical have no analytical closure.
ANALYTICAL_CLOSURE_SOURCE = 'analytical_hydrodynamic_closure'
ANALYTICAL_CLOSURE_REACTORS = frozenset({'Fluidized', 'MMBCR'})
ANALYTICAL_CLOSURE_BASIS = {
    'Fluidized': 'bubble_fraction = clip((u0-umf)/u0, 0.01, 0.5); emulsion at umf',
    'MMBCR': MMBCR_BUBBLE_RISE_BASIS + '; eps_g = u_sup/u_b; tau = H/u_b',
}
EXTERNAL_CLOSURE_SOURCES = frozenset({
    'validated_full_physics', 'calibrated_transport_surrogate'})

# Single-reactor pathway mode for callers that sweep a reactor list one
# reactor at a time (yaml_sweep, inventory_sweep, eact_sensitivity).
# Upstream rejects a sweep that mixes solids and melt modes.
SINGLE_REACTOR_MODE = {
    'PFR': 'thermocatalytic_pfr',
    'Fluidized': 'thermocatalytic_fluidized',
    'MMBCR': 'mmbcr',
}
# Material class for synthetic template-kinetics diagnostics (E_act sweeps,
# ablations) that have no screening row. Labelled diagnostic, not a candidate.
DIAGNOSTIC_MATERIAL_CLASS = {
    'PFR': 'SolidCatalyst',
    'Fluidized': 'SolidCatalyst',
    'MMBCR': 'MoltenMetal',
}

# Production solids particle size (B1-3). ROI map: last Ergun-legal
# envelope cell with margin is 0.10 mm (0.67 bar); 0.08 mm fails.
# 0.13 mm is in-band (~0.40 bar) and ~15× geometric a vs the old 2 mm.
DEFAULT_SOLIDS_PARTICLE_MM = 0.13

# Production metal inventory (B1-4). Area-fraction levers, not wt% / BET.
# 0.5 × 0.3 is the TCD-like ROI cell: supported Ni, not a bulk-metal
# pellet. Alves 2021 / Sánchez-Bastardo 2021: TCD Ni is supported at
# tens of wt%; Gili 2024: accessible metal dies to encapsulation.
# Product a = 0.15 × a_geom. Neither factor may exceed 1.
DEFAULT_METAL_LOADING = 0.5
DEFAULT_METAL_DISPERSION = 0.3

# Quantized B1-2 coarse grid (archive only). Production defaults are
# DEFAULT_SOLIDS_PARTICLE_MM × DEFAULT_METAL_LOADING × DEFAULT_METAL_DISPERSION.
INVENTORY_PARTICLE_MM = (2.0, 0.5, 0.2, 0.1)
INVENTORY_METAL_LOADING = (1.0, 0.5, 0.2)
INVENTORY_METAL_DISPERSION = (1.0, 0.3, 0.1)

# B1-2 refine ROI from the coarse grid: X only became material at
# d_p <= 0.2 mm; Ergun at 0.1 mm was 0.67 bar so ~0.08 mm is the 1 bar wall.
# Drop loading=0.2 / disp=0.1 (they only recreate the 2 mm cell).
INVENTORY_ROI_PARTICLE_MM = (0.25, 0.20, 0.16, 0.13, 0.10, 0.08)
INVENTORY_ROI_METAL_LOADING = (1.0, 0.7, 0.5)
INVENTORY_ROI_METAL_DISPERSION = (1.0, 0.5, 0.3)
INVENTORY_ROI_REASON = (
    'coarse grid: X proportional to a; gain starts at d_p<=0.5 mm and is '
    'material at <=0.2 mm; Ergun wall ~0.08 mm on this 0.5 m / 0.05 m/s bed'
)


# ═══════════════════════════════════════════════════════════════════════════════
# REACTOR CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ReactorConfig:
    """Configuration for reactor simulation."""
    # Operating conditions
    T_inlet_K: float = 1000.0        # Inlet temperature
    P_inlet_Pa: float = 101325.0     # Inlet pressure (1 atm)
    inlet_composition: str = 'CH4:0.95, Ar:0.05'  # Feed composition (Ar tracer)

    # Bubble column specific. Column height is the residence-time design
    # lever: τ = H / u_b with u_b from Mendelson (melt σ, ρ, d_b). Gas holdup
    # is derived as u_sup / u_b and must not be set by the caller.
    column_height_m: float = 1.5     # Molten metal column height
    column_diameter_m: float = 0.10  # Internal column diameter
    bubble_diameter_mm: float = 5.0  # Average bubble diameter
    gas_velocity_m_s: float = 0.05   # Superficial gas velocity
    gas_holdup_fraction: Optional[float] = None  # DERIVED (u_sup / u_b); None only
    n_cstr_stages: int = 20          # Number of melt ODE stages
    melt_surface_tension_N_m: float = MMBCR_MELT_SURFACE_TENSION_N_M
    melt_density_kg_m3: float = MMBCR_MELT_DENSITY_KG_M3

    # Packed bed specific
    bed_length_m: float = 0.5        # Catalyst bed length
    bed_diameter_m: float = 0.05     # Bed diameter (lab-scale tube reactor)
    catalyst_particle_mm: float = DEFAULT_SOLIDS_PARTICLE_MM
    bed_porosity: float = 0.4        # Void fraction
    # Γ is a monolayer. Do not raise to force Da (B1).
    site_density_mol_cm2: float = MONOLAYER_SITE_DENSITY_MOL_CM2
    # Fraction of geometric pellet surface that is metal, and of that metal
    # that is surface-available. Defaults are the supported-TCD proxy
    # (B1-4; Alves 2021; Sánchez-Bastardo 2021; Gili 2024). Neither may
    # exceed 1 — extra area is not a BET/Γ invention (B1).
    metal_loading: float = DEFAULT_METAL_LOADING
    metal_dispersion: float = DEFAULT_METAL_DISPERSION

    # Fluidized bed specific
    u_mf_m_s: float = 0.02          # Minimum fluidization velocity
    bed_height_m: float = 0.8       # Static bed height
    catalyst_density_kg_m3: float = 2500.0  # Catalyst particle density
    fluidized_bubble_fraction: Optional[float] = None  # None = (u0-umf)/u0 clipped

    # General
    reactor_type: str = 'PFR'
    mechanism_file: str = ''
    catalyst_name: str = 'test'
    material_class: Optional[str] = None
    pathway_mode: str = DEFAULT_MODE
    candidate_id: str = 'unknown'
    multiphysics_results_dir: Optional[str] = None
    multiphysics_artifact: Optional[dict] = None
    closure_features: Optional[dict] = None
    reactor_closure_evidence: Optional[dict] = None
    max_residence_time_s: float = 60.0
    catalyst_E_act_eV: float = 0.8   # Catalyst activation barrier (mock + MMBCR k_if)
    catalyst_dE_H_eV: float = 0.0

    # --- Carbon handling (reactor-specific; not one shared "decoke" flag) ---
    # MMBCR: bubble S/V + flotation. Interfacial k0 [m/s] is a melt-side
    # prefactor (not DFT). Carbon does not occupy a solid site lattice.
    # Removal rate is flotation frequency [1/s]; None = unconstrained.
    mmbcr_carbon_removal_rate_1_s: Optional[float] = MMBCR_FLOTATION_UNCONSTRAINED
    mmbcr_interfacial_k0_m_s: float = MMBCR_INTERFACIAL_K0_DEFAULT
    # PFR / batch fluidized: produce → mechanical outfeed/clear → return.
    regen_coverage_threshold: float = 0.8
    regen_mechanism: str = REGEN_MECHANICAL
    max_regen_cycles: int = 3
    # Fluidized: must be chosen explicitly.
    fluidized_mode: str = FLUIDIZED_CIRCULATING
    # Circulating fluidized / continuous removal rate [1/s].
    circulating_carbon_removal_rate_1_s: float = 0.5
    # Oxidative regen locked unless explicitly enabled for testing.
    co2_permitted: bool = False
    # Fail-closed: surface and graphite must load unless the caller asked
    # for a gas-only run. A name mismatch must not become X ≈ 0.
    gas_only: bool = False


def _validate_carbon_policy(config: ReactorConfig) -> None:
    if config.regen_mechanism not in ALLOWED_REGEN:
        raise ValueError(f'Unknown regen_mechanism={config.regen_mechanism!r}')
    if config.fluidized_mode not in ALLOWED_FLUIDIZED:
        raise ValueError(f'Unknown fluidized_mode={config.fluidized_mode!r}')
    if config.regen_mechanism == REGEN_OXIDATIVE and not config.co2_permitted:
        raise RuntimeError(
            'oxidative regen requires co2_permitted=True (default False; '
            'turquoise-compliant runs must not burn carbon to CO2)')
    if abs(config.site_density_mol_cm2 - MONOLAYER_SITE_DENSITY_MOL_CM2) > 1e-15:
        raise ValueError(
            f'site_density_mol_cm2={config.site_density_mol_cm2}; B1 locks Γ at '
            f'{MONOLAYER_SITE_DENSITY_MOL_CM2} mol/cm^2')
    if not (0.0 < config.metal_loading <= 1.0):
        raise ValueError(
            f'metal_loading={config.metal_loading} must be in (0, 1]; '
            'do not invent area above geometric')
    if not (0.0 < config.metal_dispersion <= 1.0):
        raise ValueError(
            f'metal_dispersion={config.metal_dispersion} must be in (0, 1]; '
            'do not invent area above geometric')
    if not (0.0 < config.mmbcr_interfacial_k0_m_s <= MMBCR_INTERFACIAL_K0_MAX):
        raise ValueError(
            f'mmbcr_interfacial_k0_m_s={config.mmbcr_interfacial_k0_m_s} '
            f'must be in (0, {MMBCR_INTERFACIAL_K0_MAX}]; '
            f'{MMBCR_INTERFACIAL_K0_MAX} m/s is a guardrail against inventing '
            'Da, not a physical bound like loading <= 1 (B3)')
    if (config.mmbcr_carbon_removal_rate_1_s is not None
            and config.mmbcr_carbon_removal_rate_1_s < 0.0):
        raise ValueError(
            'mmbcr_carbon_removal_rate_1_s must be None (unconstrained) or >= 0')
    if config.circulating_carbon_removal_rate_1_s < 0.0:
        raise ValueError('circulating_carbon_removal_rate_1_s must be >= 0')


def _validate_reactor_config(config: ReactorConfig) -> None:
    """Reject geometries that do not represent the selected physical bed."""
    if config.T_inlet_K <= 0 or config.P_inlet_Pa <= 0:
        raise ValueError('reactor temperature and pressure must be positive')
    if config.reactor_type == 'PFR':
        if config.bed_length_m <= 0 or config.bed_diameter_m <= 0 or \
                config.catalyst_particle_mm <= 0 or config.gas_velocity_m_s <= 0:
            raise ValueError('packed-bed dimensions, particle size, and flow must be positive')
        if not 0 < config.bed_porosity < 1:
            raise ValueError('packed-bed porosity must lie strictly between zero and one')
    elif config.reactor_type == 'Fluidized':
        if config.bed_height_m <= 0 or config.catalyst_particle_mm <= 0 or \
                config.u_mf_m_s <= 0:
            raise ValueError('fluidized-bed dimensions, particle size, and umf must be positive')
        if config.gas_velocity_m_s <= config.u_mf_m_s:
            raise ValueError('gas velocity must exceed minimum fluidization velocity')
    elif config.reactor_type == 'MMBCR':
        if config.column_height_m <= 0 or config.column_diameter_m <= 0 or \
                config.bubble_diameter_mm <= 0 or config.gas_velocity_m_s <= 0 or \
                config.n_cstr_stages < 1:
            raise ValueError('MMBCR dimensions, flow, bubbles, and stage count must be positive')
        if config.melt_surface_tension_N_m <= 0 or config.melt_density_kg_m3 <= 0:
            raise ValueError('MMBCR melt surface tension and density must be positive')
        _mmbcr_hydrodynamics(config)  # rejects caller-set holdup and churn-turbulent flow


def _mechanism_metadata(config: ReactorConfig) -> dict:
    path = Path(config.mechanism_file).with_suffix('.kinetics.json') if config.mechanism_file else None
    if path is None or not path.exists():
        return {'inputs': {'quantitative_status': 'missing_provenance'},
                'carbon_phase_model': 'unknown'}
    try:
        return json.loads(path.read_text(encoding='utf-8'))
    except Exception:
        return {'inputs': {'quantitative_status': 'invalid_provenance'},
                'carbon_phase_model': 'unknown'}


def _kinetics_evidence(config: ReactorConfig) -> dict:
    """Declare whether reactor output may make a candidate-level decision."""
    metadata = _mechanism_metadata(config)
    status = metadata.get('inputs', {}).get(
        'quantitative_status', 'missing_provenance')
    carbon_model = metadata.get('carbon_phase_model', 'unknown')
    limitations = []
    if status != 'candidate_specific':
        limitations.append('incomplete_candidate_kinetics')
    if carbon_model == 'legacy_gas_tracer':
        limitations.append('legacy_gas_carbon_tracer')
    elif carbon_model == 'unknown':
        limitations.append('unknown_carbon_phase_model')
    complete = not limitations
    return {
        'kinetics_status': status,
        'carbon_phase_model': carbon_model,
        'reactor_evidence_tier': (
            'candidate_specific_kinetics' if complete else
            'diagnostic_screening_template'),
        # Incomplete template kinetics can guide sensitivity/validation but may
        # never eliminate a candidate or count as reactor validation.
        'can_exclude_candidate': bool(complete),
        'reactor_evidence_limitations': limitations,
    }


def _kinetics_metadata(config: ReactorConfig) -> dict:
    return _mechanism_metadata(config).get('inputs', {})


def _policy_metadata(config: ReactorConfig) -> Dict:
    return {
        'co2_permitted': bool(config.co2_permitted),
        'regen_mechanism': config.regen_mechanism,
        'max_regen_cycles': int(config.max_regen_cycles),
        'regen_coverage_threshold': float(config.regen_coverage_threshold),
        'mmbcr_carbon_removal_rate_1_s': config.mmbcr_carbon_removal_rate_1_s,
        'mmbcr_carbon_removal_role': 'interface_flotation_frequency_1_s',
        'mmbcr_flotation': _mmbcr_flotation_label(config),
        'mmbcr_carbon_removal_wired': True,
        'fluidized_mode': config.fluidized_mode,
        'circulating_carbon_removal_rate_1_s': float(
            config.circulating_carbon_removal_rate_1_s),
        'circulating_carbon_removal_when': 'during_integrate_substeps',
        'fluidized_bypass_basis': 'ar_tracer_molar_mix',
        'carbon_phase_model': 'condensed_graphite_plus_surface_C_s',
        'mmbcr_rate_model': 'bubble_area_flotation',
        'mmbcr_interfacial_k0_m_s': float(config.mmbcr_interfacial_k0_m_s),
        'mmbcr_interfacial_k0_basis': (
            'calibrated melt-side prefactor; not DFT or Upham-fitted k; '
            'B3; Upham 2017; Chen 2023; Abdollahi 2024'),
        'mmbcr_interfacial_k0_max_m_s': MMBCR_INTERFACIAL_K0_MAX,
        'mmbcr_interfacial_k0_max_kind': MMBCR_INTERFACIAL_K0_MAX_KIND,
        'mmbcr_flotation_default': 'unconstrained (eta=1); knob wired but off',
        'mmbcr_bubble_rise_basis': MMBCR_BUBBLE_RISE_BASIS,
        'mmbcr_residence_time_basis': MMBCR_RESIDENCE_BASIS,
        'mmbcr_gas_holdup_basis': 'derived u_sup / u_b; not an input',
        'mmbcr_max_gas_holdup': MMBCR_MAX_GAS_HOLDUP,
        'melt_surface_tension_N_m': float(config.melt_surface_tension_N_m),
        'melt_density_kg_m3': float(config.melt_density_kg_m3),
        'melt_property_basis': 'Ni-Bi order of magnitude; not fitted',
        'solids_thermal_mode': 'isothermal_energy_disabled',
        'solids_conversion_basis': 'argon_tracer',
        'h2_metric_note': (
            'H2_atom_balance is not branching selectivity; lumped mechanism '
            'has no C2 competition branch for true H2 selectivity'),
        'site_density_mol_cm2': float(config.site_density_mol_cm2),
        'site_density_basis': 'monolayer_2.5e-9_mol_cm2',
        'metal_loading': float(config.metal_loading),
        'metal_dispersion': float(config.metal_dispersion),
        'solids_area_model': 'a_geom * loading * dispersion (both <= 1)',
        'solids_loading_basis': (
            'supported_tcd_area_fraction; Alves 2021; '
            'Sanchez-Bastardo 2021; Gili 2024'),
    }


def geometric_sv_pfr(config: ReactorConfig) -> float:
    """External pellet area per bed volume: 6(1−ε)/d_p."""
    d_p = config.catalyst_particle_mm * 1e-3
    if d_p <= 0:
        raise ValueError('catalyst_particle_mm must be positive')
    return 6.0 * (1.0 - config.bed_porosity) / d_p


def geometric_sv_fluidized(config: ReactorConfig) -> float:
    """Emulsion solids area per emulsion volume: 6·(1−ε_mf)/d_p."""
    d_p = config.catalyst_particle_mm * 1e-3
    if d_p <= 0:
        raise ValueError('catalyst_particle_mm must be positive')
    return 6.0 * (1.0 - FLUIDIZED_EMULSION_VOIDAGE) / d_p


def active_area_multiplier(config: ReactorConfig) -> float:
    return float(config.metal_loading) * float(config.metal_dispersion)


def active_sv(geometric_sv: float, config: ReactorConfig) -> float:
    return float(geometric_sv) * active_area_multiplier(config)


def ch4_feed_density_kg_m3(T_K: float, P_Pa: float) -> float:
    """Ideal-gas density for CH4:0.95 / Ar:0.05."""
    M = 0.01604 * 0.95 + 0.03995 * 0.05
    return P_Pa * M / (8.314462618 * T_K)


def ch4_viscosity_pa_s(T_K: float) -> float:
    """Sutherland estimate for CH4 (μ0=1.03e-5 Pa·s at 273.15 K, S=164 K)."""
    T0, mu0, S = 273.15, 1.03e-5, 164.0
    return mu0 * (T_K / T0) ** 1.5 * (T0 + S) / (T_K + S)


def ergun_delta_p_pa(config: ReactorConfig, T_K: float = None,
                     P_Pa: float = None) -> float:
    """Packed-bed Ergun ΔP over bed_length_m (B4). Not used for MMBCR."""
    T = config.T_inlet_K if T_K is None else T_K
    P = config.P_inlet_Pa if P_Pa is None else P_Pa
    d_p = config.catalyst_particle_mm * 1e-3
    eps = config.bed_porosity
    u = config.gas_velocity_m_s if config.gas_velocity_m_s > 0 else 0.1
    mu = ch4_viscosity_pa_s(T)
    rho = ch4_feed_density_kg_m3(T, P)
    viscous = 150.0 * mu * (1.0 - eps) ** 2 / (eps ** 3 * d_p ** 2) * u
    inertial = 1.75 * rho * (1.0 - eps) / (eps ** 3 * d_p) * u ** 2
    return float(config.bed_length_m * (viscous + inertial))


def reciprocal_residence_h(tau_s: float) -> float:
    """Space velocity as 1/τ in h⁻¹. Field name WHSV is historical."""
    return 3600.0 / float(tau_s) if tau_s and float(tau_s) > 0 else 0.0


def kinetics_fields(config: ReactorConfig) -> Dict:
    return {
        'catalyst_E_act_eV': float(config.catalyst_E_act_eV),
        'catalyst_dE_H_eV': float(config.catalyst_dE_H_eV),
    }


def solids_inventory_fields(config: ReactorConfig, geometric_sv: float) -> Dict:
    a = active_sv(geometric_sv, config)
    dp = ergun_delta_p_pa(config)
    return {
        'geometric_sv_1_m': float(geometric_sv),
        'active_sv_1_m': float(a),
        'active_area_multiplier': active_area_multiplier(config),
        'ergun_delta_p_Pa': dp,
        'ergun_delta_p_bar': dp / 1e5,
        'ergun_ok': dp <= ERGUN_DP_LIMIT_PA,
    }


def inventory_grid_cells(particle_mm=None, loadings=None, dispersions=None):
    """Quantized (d_p, loading, dispersion) grid. Defaults to the coarse B1-2 set."""
    d_levels = INVENTORY_PARTICLE_MM if particle_mm is None else particle_mm
    w_levels = INVENTORY_METAL_LOADING if loadings is None else loadings
    s_levels = INVENTORY_METAL_DISPERSION if dispersions is None else dispersions
    cells = []
    for d_p in d_levels:
        for loading in w_levels:
            for dispersion in s_levels:
                cells.append({
                    'catalyst_particle_mm': d_p,
                    'metal_loading': loading,
                    'metal_dispersion': dispersion,
                })
    return cells


def inventory_roi_grid_cells():
    return inventory_grid_cells(
        INVENTORY_ROI_PARTICLE_MM,
        INVENTORY_ROI_METAL_LOADING,
        INVENTORY_ROI_METAL_DISPERSION,
    )


def _load_candidate_phases(config: ReactorConfig):
    """Upstream name: gas plus the required candidate-specific surface."""
    gas, _graphite, surf = _load_gas_and_surface(config)
    if surf is None:
        raise RuntimeError(
            f'candidate surface phase failed to load: '
            f'{config.catalyst_name}_surface')
    return gas, surf


def _load_gas_and_surface(config: ReactorConfig):
    gas = ct.Solution(config.mechanism_file, 'gas')
    if 'C_graphite' in gas.species_names:
        raise RuntimeError(
            'mechanism still contains gas-phase C_graphite; regenerate YAML')
    graphite = None
    try:
        graphite = ct.Solution(config.mechanism_file, 'graphite')
    except Exception as exc:
        if not config.gas_only:
            raise RuntimeError(
                f'graphite phase missing from {config.mechanism_file}; '
                'set gas_only=True only for an explicit gas-only run'
            ) from exc
    surf = None
    surf_name = f'{config.catalyst_name}_surface'
    try:
        adjacent = [gas]
        if graphite is not None:
            adjacent.append(graphite)
        try:
            surf = ct.Interface(config.mechanism_file, surf_name, adjacent)
        except Exception:
            surf = ct.Interface(config.mechanism_file, surf_name, [gas])
    except Exception as exc:
        if not config.gas_only:
            raise RuntimeError(
                f'surface {surf_name!r} missing from {config.mechanism_file} '
                f'(catalyst_name={config.catalyst_name!r} must match the '
                'YAML written by write_full_mechanism); set gas_only=True '
                'only for an explicit gas-only run'
            ) from exc
    return gas, graphite, surf


def _load_status_fields(config: ReactorConfig, graphite, surf) -> dict:
    return {
        'gas_only': bool(config.gas_only),
        'surface_loaded': surf is not None,
        'graphite_loaded': graphite is not None,
        'surface_name': f'{config.catalyst_name}_surface',
    }


def _species_x(gas, name: str) -> float:
    if name not in gas.species_names:
        return 0.0
    return float(gas.X[gas.species_index(name)])


def _coverage(surf, name: str) -> float:
    if surf is None or name not in surf.species_names:
        return 0.0
    return float(surf.coverages[surf.species_index(name)])


def _apply_continuous_carbon_removal(surf, rate_1_s: float, dt: float) -> float:
    """
    Transport-style continuous removal of C_s → free sites.

    Returns approximate coverage of C removed (not moles). This is a mass-transport
    lump wearing kinetics clothing — not Arrhenius chemistry.
    """
    if surf is None or rate_1_s <= 0 or dt <= 0:
        return 0.0
    if 'C_s' not in surf.species_names or 'site' not in surf.species_names:
        return 0.0
    cov = np.array(surf.coverages, dtype=float)
    i_c = surf.species_index('C_s')
    i_site = surf.species_index('site')
    c_before = cov[i_c]
    removed = c_before * (1.0 - np.exp(-rate_1_s * dt))
    cov[i_c] = c_before - removed
    cov[i_site] += removed
    # Renormalize site-occupying coverages only (exclude sites==0 species if any).
    site_mask = np.array([surf.species(n).size > 0 for n in range(surf.n_species)])
    s = cov[site_mask].sum()
    if s > 0:
        cov[site_mask] /= s
    surf.coverages = cov
    return float(removed)


def _reset_surface_carbon(surf) -> None:
    """Mechanical / consumable regen: clear adsorbates and restore free sites.

    Encapsulating carbon (``C_encap_s``, B6 Cδ) is time-on-stream death and
    is never cleared: mechanical outfeed moves the particle, it does not
    strip graphene off the face. Only oxidative burn-off would, and that is
    not modelled (non-turquoise).
    """
    if surf is None or 'C_s' not in surf.species_names:
        return
    encap = _coverage(surf, 'C_encap_s') if 'C_encap_s' in surf.species_names else 0.0
    cov = np.zeros(surf.n_species)
    if 'C_encap_s' in surf.species_names:
        cov[surf.species_index('C_encap_s')] = encap
    if 'site' in surf.species_names:
        cov[surf.species_index('site')] = max(0.0, 1.0 - encap)
    surf.coverages = cov


METAL_MOLAR_MASS_G_MOL = {'Ni': 58.6934, 'Fe': 55.845, 'Co': 58.9332}
# Ermakova 2000 / Takenaka: 40-384 gC/gNi over 4-50 h on high-Ni/SiO2 -> ~8-10 gC/(gNi h).
NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H = (8.0, 10.0)
# Same TOS literature: the catalyst dies (encapsulated) after 4-50 h.
NI_TOS_LIFETIME_BAND_H = (4.0, 50.0)
# B6-6 encapsulation-onset markers. Primary: θ_encap at the bed exit has
# reached ENCAP_ONSET_THETA. Secondary: Cγ/Cδ below ENCAP_ONSET_RATIO.
# Lifetime: hours for θ_encap to reach ENCAP_DEAD_THETA at the pass-
# averaged Cδ rate (linear extrapolation; Cδ ∝ θ_C² makes this a lower
# bound while θ_C is still climbing).
ENCAP_ONSET_THETA = 0.1
ENCAP_ONSET_RATIO = 10.0
ENCAP_DEAD_THETA = 0.5


def _surface_carbon_coverage(surf, cov: np.ndarray) -> float:
    """Total θ of carbon-bearing site-occupying species (C_s, CH_x_s, C_encap_s)."""
    total = 0.0
    for k in range(surf.n_species):
        sp = surf.species(k)
        if sp.size > 0 and sp.composition.get('C', 0) > 0:
            total += float(cov[k])
    return total


def _off_site_carbon_metrics(config: ReactorConfig, gas, surf, *,
                             ch4_initial: float, ar_initial: float,
                             pass_conversion: float, n_sites_mol: float,
                             n_ch4_fed_mol: float, n_parcel_in_mol: float,
                             pass_time_s: float,
                             cov_start: Optional[np.ndarray],
                             removed_surface_carbon_mol: float = 0.0,
                             basis: str = 'Gamma*A_stage / (c_CH4*eps*V_stage), PFR parcel'
                             ) -> Dict:
    """Per-pass carbon accounting for the B5/B6 closure criterion.

    ``site_inventory_bound_X`` is the conversion a stoichiometric monolayer
    can deliver (Γ·a / (ε·c_CH4) on the PFR parcel basis). Turnovers above
    1 need Cγ returning sites. Cγ is obtained by carbon balance (solid carbon
    from the Ar tracer minus carbon still on the surface minus carbon taken
    off by a circulating outfeed, ``removed_surface_carbon_mol``); Cδ is
    the change in ``C_encap_s`` coverage. The yield rate is pass-averaged.
    """
    x_bound = n_sites_mol / n_ch4_fed_mol if n_ch4_fed_mol > 0 else None
    x_eq = _tabulated_x_eq(config.T_inlet_K)
    out = {
        'site_inventory_bound_X': x_bound,
        'site_inventory_bound_basis': basis,
        # Cγ is irreversible into a graphite sink; the surface mechanism
        # does not enforce CH4 <=> C(gr) + 2 H2 from the solid side, so a
        # fast Cγ can overshoot equilibrium. Flagged, never clipped.
        'X_eq_table': x_eq,
        'exceeds_equilibrium': bool(pass_conversion > x_eq + 1e-6),
    }
    if surf is None:
        return out
    cov_end = np.array(surf.coverages, dtype=float)
    if cov_start is None:
        cov_start = np.zeros_like(cov_end)
    # Solid carbon this pass: CH4 converted minus carbon left in C2 gas (Ar basis).
    x_ar = _species_x(gas, 'Ar') or ar_initial
    n_ar = ar_initial * n_parcel_in_mol
    n_c2 = sum(
        2.0 * (_species_x(gas, sp) or 0.0) / max(x_ar, 1e-30) * n_ar
        for sp in ('C2H2', 'C2H4', 'C2H6'))
    n_solid = max(0.0, pass_conversion * n_ch4_fed_mol - n_c2)
    d_surface_c = (_surface_carbon_coverage(surf, cov_end)
                   - _surface_carbon_coverage(surf, cov_start)) * n_sites_mol
    n_removed_raw = max(0.0, float(removed_surface_carbon_mol))
    # Exported carbon cannot exceed CH4 converted onto the solids this pass
    # minus carbon that remained on the surface. A raw outfeed larger than
    # that is a state-accounting leak (B2 circulating Interface vs
    # ReactorSurface), not extra Cγ. Conserve before scoring filaments.
    available_to_export = max(0.0, n_solid - d_surface_c)
    n_removed = min(n_removed_raw, available_to_export)
    balance_residual = n_solid - d_surface_c - n_removed_raw
    n_gamma = max(0.0, n_solid - d_surface_c - n_removed)
    has_encap = 'C_encap_s' in surf.species_names
    n_delta = 0.0
    if has_encap:
        i = surf.species_index('C_encap_s')
        n_delta = max(0.0, float(cov_end[i] - cov_start[i])) * n_sites_mol
    turnovers = n_solid / n_sites_mol if n_sites_mol > 0 else None
    ratio = (n_gamma / n_delta) if (has_encap and n_delta > 0) else None
    theta_encap = (
        float(cov_end[surf.species_index('C_encap_s')]) if has_encap else None)
    # Lifetime: θ_encap → ENCAP_DEAD_THETA at the pass-averaged Cδ rate.
    lifetime_h = None
    if has_encap and n_sites_mol > 0 and pass_time_s > 0:
        d_theta_dt = n_delta / n_sites_mol / pass_time_s
        # None when no Cδ happened this pass (no finite lifetime to report).
        lifetime_h = (ENCAP_DEAD_THETA / d_theta_dt / 3600.0
                      if d_theta_dt > 0 else None)
    out.update({
        'carbon_turnovers_per_site': turnovers,
        'turnover_factor_vs_bound': (
            pass_conversion / x_bound if x_bound else None),
        'solid_carbon_mol_per_pass': n_solid,
        'outfeed_carbon_mol_per_pass': n_removed,
        'outfeed_carbon_mol_raw': n_removed_raw,
        'carbon_balance_residual_mol': balance_residual,
        'carbon_balance_ok': bool(
            balance_residual >= -1e-3 * max(n_sites_mol, n_solid, 1e-12)),
        'c_gamma_mol_per_pass': n_gamma if has_encap else None,
        'c_delta_mol_per_pass': n_delta if has_encap else None,
        'c_gamma_to_c_delta_ratio': ratio,
        'exit_theta_C_encap': theta_encap,
        'off_site_carbon_active': has_encap,
    })
    if has_encap:
        out.update({
            'encapsulation_onset': bool(theta_encap >= ENCAP_ONSET_THETA),
            'encapsulation_onset_basis': f'exit theta_encap >= {ENCAP_ONSET_THETA}',
            'c_delta_competitive': (
                None if ratio is None else bool(ratio < ENCAP_ONSET_RATIO)),
            'c_delta_competitive_basis': f'C_gamma/C_delta < {ENCAP_ONSET_RATIO}',
            'encapsulation_lifetime_h': lifetime_h,
            'encapsulation_lifetime_basis': (
                f'theta_encap -> {ENCAP_DEAD_THETA} at the pass-averaged '
                'C_delta rate; linear extrapolation'),
            'tos_lifetime_literature_band_h': list(NI_TOS_LIFETIME_BAND_H),
            'lifetime_within_tos_band': (
                None if lifetime_h is None else bool(
                    NI_TOS_LIFETIME_BAND_H[0] <= lifetime_h
                    <= NI_TOS_LIFETIME_BAND_H[1])),
        })
    if has_encap:
        metadata = _mechanism_metadata(config)
        genome = metadata.get('inputs', {}).get('genome')
        metal = None
        try:
            from pipeline.process.reactor_mechanisms import (
                OFF_SITE_CARBON_METALS, _nanoparticle_metals, parse_catalyst_genome)
            metals = _nanoparticle_metals(
                parse_catalyst_genome(genome),
                metadata.get('inputs', {}).get('material_class')) & OFF_SITE_CARBON_METALS
            if len(metals) == 1:
                metal = next(iter(metals))
        except Exception:
            metal = None
        dispersion = float(config.metal_dispersion)
        yield_rate = None
        if metal and dispersion > 0 and pass_time_s > 0 and n_sites_mol > 0:
            n_metal = n_sites_mol / dispersion
            yield_rate = (n_gamma / pass_time_s) * 3600.0 * 12.011 / (
                n_metal * METAL_MOLAR_MASS_G_MOL[metal])
        out.update({
            'filament_metal': metal,
            'filament_yield_gC_per_gMetal_h': yield_rate,
            'filament_yield_basis': (
                'pass-averaged C_gamma rate; metal moles = sites / dispersion'),
            'filament_yield_literature_band_gC_per_gNi_h': list(
                NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H),
            'filament_yield_within_band': (
                None if yield_rate is None else bool(
                    NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H[0] <= yield_rate
                    <= NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H[1])),
        })
    return out


def _h2_atom_balance_metric(ch4_initial: float, final_conv: float, x_h2: float) -> float:
    """H-atom balance metric — not true branching H2 selectivity."""
    if final_conv <= 0.01:
        return 0.0
    h_in_ch4 = 4.0 * ch4_initial
    h_in_h2 = 2.0 * x_h2
    return float(np.clip(h_in_h2 / max(h_in_ch4 * final_conv, 1e-10), 0, 1))


def _solid_c_from_balance(ch4_initial: float, final_conv: float,
                          x_c2h2: float, x_c2h4: float, x_c2h6: float) -> float:
    c_in_c2 = 2.0 * (x_c2h2 + x_c2h4 + x_c2h6)
    c_to_solid = final_conv * ch4_initial - c_in_c2
    if final_conv <= 0.01:
        return 0.0
    return float(np.clip(c_to_solid / max(final_conv * ch4_initial, 1e-10), 0, 1))


def _tabulated_x_eq(T_K: float) -> float:
    from pipeline.process.equilibrium_check import TABULATED_X_CH4_1BAR
    if T_K in TABULATED_X_CH4_1BAR:
        return float(TABULATED_X_CH4_1BAR[T_K])
    nearest = min(TABULATED_X_CH4_1BAR, key=lambda t: abs(t - T_K))
    return float(TABULATED_X_CH4_1BAR[nearest])


def _mmbcr_interfacial_k_m_s(E_act_eV: float, T_K: float, k0_m_s: float) -> float:
    k_B_eV = 8.617333262e-5
    return float(k0_m_s * np.exp(-E_act_eV / max(k_B_eV * T_K, 1e-12)))


def _mmbcr_flotation_label(config: ReactorConfig) -> str:
    rate = config.mmbcr_carbon_removal_rate_1_s
    if rate is None:
        return 'unconstrained'
    if rate <= 0.0:
        return 'blocked'
    return 'finite'


def _mmbcr_flotation_eta(k_if_m_s: float, sv_ratio_1_m: float,
                         k_float_1_s: Optional[float]) -> float:
    """Available-interface factor. None = unconstrained; 0 = fouled."""
    if k_float_1_s is None:
        return 1.0
    if k_float_1_s <= 0.0:
        return 0.0
    denom = k_float_1_s + max(k_if_m_s, 0.0) * max(sv_ratio_1_m, 0.0)
    if denom <= 0:
        return 0.0
    return float(k_float_1_s / denom)


def mendelson_bubble_rise_velocity_m_s(d_b_m: float, sigma_N_m: float,
                                       rho_kg_m3: float) -> float:
    """Mendelson (1967) terminal rise velocity: sqrt(2σ/(ρ d_b) + g d_b/2)."""
    if d_b_m <= 0 or sigma_N_m <= 0 or rho_kg_m3 <= 0:
        raise ValueError('bubble diameter, melt surface tension, and density must be positive')
    return float(np.sqrt(2.0 * sigma_N_m / (rho_kg_m3 * d_b_m) + G_M_S2 * d_b_m / 2.0))


def _mmbcr_hydrodynamics(config: ReactorConfig) -> Dict:
    """Derive u_b, τ, and gas holdup from geometry; fail closed outside bubbly flow.

    Column height is the residence-time lever (τ = H / u_b). Gas holdup is
    ε_g = u_sup / u_b and is not an input; above MMBCR_MAX_GAS_HOLDUP the
    column is churn-turbulent and the single-bubble ODE no longer applies.
    """
    d_b = config.bubble_diameter_mm * 1e-3
    evidence = config.reactor_closure_evidence or {}
    external = evidence.get('source') in EXTERNAL_CLOSURE_SOURCES
    if config.gas_holdup_fraction is not None and not external:
        raise ValueError(
            f'gas_holdup_fraction={config.gas_holdup_fraction} was set explicitly '
            'without a validated artifact or calibrated surrogate; MMBCR holdup '
            'is derived as gas_velocity_m_s / u_b (Mendelson) and must be left '
            'as None')
    if external and config.gas_holdup_fraction is not None:
        # Artifact / surrogate closure wins: eps_g is an external output,
        # u_b follows from it, tau = eps_g * H / u_sup (gas residence).
        eps_g = float(config.gas_holdup_fraction)
        if not 0.0 < eps_g < 1.0:
            raise ValueError(
                f'artifact gas_holdup_fraction={eps_g} must lie in (0, 1)')
        u_b = config.gas_velocity_m_s / eps_g
        tau_total = config.column_height_m * eps_g / config.gas_velocity_m_s
        rise_basis = f'u_sup / eps_g from {evidence.get("source")}'
        holdup_basis = evidence.get('source')
        residence_basis = MMBCR_ARTIFACT_RESIDENCE_BASIS
    else:
        u_b = mendelson_bubble_rise_velocity_m_s(
            d_b, config.melt_surface_tension_N_m, config.melt_density_kg_m3)
        eps_g = config.gas_velocity_m_s / u_b
        if eps_g > MMBCR_MAX_GAS_HOLDUP:
            raise ValueError(
                f'derived gas holdup eps_g={eps_g:.3f} exceeds {MMBCR_MAX_GAS_HOLDUP} '
                f'(churn-turbulent; bubbly-flow melt ODE invalid). u_b={u_b:.3f} m/s '
                f'from Mendelson at d_b={config.bubble_diameter_mm} mm. Reduce '
                f'gas_velocity_m_s={config.gas_velocity_m_s} or raise column_height_m'
                f'={config.column_height_m} for residence time instead')
        tau_total = config.column_height_m / u_b
        rise_basis = MMBCR_BUBBLE_RISE_BASIS
        holdup_basis = 'derived u_sup / u_b'
        residence_basis = MMBCR_RESIDENCE_BASIS
    return {
        'bubble_rise_velocity_m_s': float(u_b),
        'bubble_rise_basis': rise_basis,
        'gas_holdup_fraction': float(eps_g),
        'gas_holdup_basis': holdup_basis,
        'residence_time_s': float(tau_total),
        'residence_time_basis': residence_basis,
        'interfacial_sv_ratio_1_m': 6.0 / d_b,
        'melt_surface_tension_N_m': float(config.melt_surface_tension_N_m),
        'melt_density_kg_m3': float(config.melt_density_kg_m3),
        'melt_property_basis': 'Ni-Bi order of magnitude; not fitted',
    }


def _ch4_extent(gas, x_ch4_feed: float, x_ar_feed: float) -> float:
    from pipeline.process.equilibrium_check import ch4_conversion_from_argon_tracer
    return ch4_conversion_from_argon_tracer(
        _species_x(gas, 'CH4'), _species_x(gas, 'Ar'), x_ch4_feed, x_ar_feed)


def _require_ar_tracer(x_ch4_feed: float, x_ar_feed: float) -> None:
    if x_ar_feed <= 0 or x_ch4_feed <= 0:
        raise ValueError(
            'solids CH4 conversion uses an Ar mole-fraction tracer; '
            'inlet_composition must include Ar (default CH4:0.95, Ar:0.05)')


def _disable_reactor_energy(reactor) -> str:
    """Hold solids stages at inlet T so they match isothermal MMBCR."""
    if hasattr(reactor, 'energy_enabled'):
        reactor.energy_enabled = False
        return 'isothermal_energy_disabled'
    raise RuntimeError(
        'Cantera reactor has no energy_enabled; cannot force isothermal solids')


def _mole_fraction_drop(gas, x_ch4_feed: float) -> float:
    if x_ch4_feed <= 0:
        return 0.0
    return float(max(0.0, 1.0 - _species_x(gas, 'CH4') / x_ch4_feed))


def _set_gas_from_ch4_conversion(gas, T_K: float, P_Pa: float,
                                 x_ch4_feed: float, x_ar_feed: float, X: float):
    """CH4 → C(s) + 2 H2; C leaves the bubble by flotation (not in the gas)."""
    X = float(np.clip(X, 0.0, 1.0))
    n_ch4 = x_ch4_feed * (1.0 - X)
    n_h2 = 2.0 * x_ch4_feed * X
    n_ar = x_ar_feed
    n_tot = n_ch4 + n_h2 + n_ar
    if n_tot <= 0:
        return
    gas.TPX = T_K, P_Pa, f'CH4:{n_ch4 / n_tot}, H2:{n_h2 / n_tot}, Ar:{n_ar / n_tot}'


def _mix_bubble_bypass(gas, inlet_x: np.ndarray, delta: float,
                       x_ar_feed: float) -> None:
    """Mix a bypass stream (fraction δ of feed, unreacted) with the emulsion outlet
    on a molar basis. Ar is conserved, so the emulsion stream's molar expansion
    per mole of feed is x_Ar,feed / x_Ar,emulsion."""
    x_em = np.array(gas.X, dtype=float)
    x_ar_em = _species_x(gas, 'Ar')
    if x_ar_em <= 0:
        raise RuntimeError('Ar tracer vanished from emulsion outlet; cannot mix bypass')
    expansion = x_ar_feed / x_ar_em
    moles = delta * np.asarray(inlet_x, dtype=float) + (1.0 - delta) * expansion * x_em
    total = float(moles.sum())
    if total <= 0:
        raise RuntimeError('bubble bypass mixing produced no moles')
    gas.TPX = gas.T, gas.P, moles / total


# ═══════════════════════════════════════════════════════════════════════════════
# A. MMBCR — bubble area + carbon flotation (no solid site lattice)
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_mmbcr(config: ReactorConfig) -> Dict:
    """Melt ODE to tabulated X_eq with bubble-area flotation (B3).

    dX/dz-style first-order approach: per stage
    ``X ← X_eq − (X_eq − X)·exp(−Da_stage)`` with
    ``Da = k_if(E_act, T)·(6/d_b)·τ·η``. τ = H / u_b (Mendelson). Column
    height is the design lever; gas holdup is derived. Carbon floats out of
    the bubble and is reconstructed from the CH4/H2/Ar balance.

    Args:
        config: Configuration controlling this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    _validate_carbon_policy(config)
    hydro = _mmbcr_hydrodynamics(config)
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'MMBCR')

    logger.info(f"Simulating MMBCR: {config.catalyst_name} at {config.T_inlet_K} K")
    gas, graphite, surf = _load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    x_ch4_feed = _species_x(gas, 'CH4') or 0.95
    x_ar_feed = _species_x(gas, 'Ar')
    if x_ar_feed <= 0:
        x_ar_feed = max(0.0, 1.0 - x_ch4_feed)

    tau_total = hydro['residence_time_s']
    tau_stage = tau_total / config.n_cstr_stages
    sv_ratio = hydro['interfacial_sv_ratio_1_m']  # bubble S/V (m² interface / m³ gas)
    x_eq = _tabulated_x_eq(config.T_inlet_K)
    k_if = _mmbcr_interfacial_k_m_s(
        config.catalyst_E_act_eV, config.T_inlet_K, config.mmbcr_interfacial_k0_m_s)
    eta_float = _mmbcr_flotation_eta(
        k_if, sv_ratio, config.mmbcr_carbon_removal_rate_1_s)
    da_stage = k_if * sv_ratio * tau_stage * eta_float

    z_positions = np.linspace(0, config.column_height_m, config.n_cstr_stages + 1)
    conversion_profile = [0.0]
    temperature_profile = [config.T_inlet_K]
    species_profiles = {sp: [_species_x(gas, sp)] for sp in
                        ['CH4', 'H2', 'C2H2', 'C2H4', 'C2H6']}
    conv = 0.0
    carbon_removed_coverage = 0.0

    for _stage in range(config.n_cstr_stages):
        # First-order approach to melt/gas equilibrium; C floats out of the bubble.
        conv = x_eq - (x_eq - conv) * np.exp(-da_stage)
        carbon_removed_coverage += max(0.0, conv - conversion_profile[-1]) * x_ch4_feed
        _set_gas_from_ch4_conversion(
            gas, config.T_inlet_K, config.P_inlet_Pa, x_ch4_feed, x_ar_feed, conv)
        conversion_profile.append(float(conv))
        temperature_profile.append(config.T_inlet_K)
        for sp in species_profiles:
            species_profiles[sp].append(_species_x(gas, sp))

    final_conv = conversion_profile[-1]
    x_h2 = _species_x(gas, 'H2')

    result = {
        'reactor_type': 'MMBCR',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'P_Pa': config.P_inlet_Pa,
        'column_height_m': config.column_height_m,
        'column_diameter_m': config.column_diameter_m,
        'gas_velocity_m_s': config.gas_velocity_m_s,
        'bubble_diameter_mm': config.bubble_diameter_mm,
        **hydro,
        'CH4_conversion': float(final_conv),
        'H2_atom_balance': _h2_atom_balance_metric(x_ch4_feed, final_conv, x_h2),
        # Backward-compatible alias; not true selectivity.
        'H2_selectivity': _h2_atom_balance_metric(x_ch4_feed, final_conv, x_h2),
        'solid_C_selectivity': None,
        'solid_C_selectivity_note': (
            'melt reconstructs CH4/H2/Ar only; C2s are not in the bubble gas'),
        'c2_tracked': False,
        'exit_x_H2': x_h2,
        'exit_x_CH4': _species_x(gas, 'CH4'),
        'exit_x_C2H2': None,
        'exit_x_C2H4': None,
        'exit_x_C2H6': None,
        'exit_T_K': float(gas.T),
        'exit_theta_C': 0.0,
        'carbon_removed_coverage_proxy': float(carbon_removed_coverage),
        'X_eq_table': float(x_eq),
        'mmbcr_k_if_m_s': float(k_if),
        'mmbcr_flotation_eta': float(eta_float),
        'mmbcr_Da': float(k_if * sv_ratio * tau_total * eta_float),
        'conversion_basis': 'melt_ode_to_Xeq',
        'thermal_mode': 'isothermal_by_construction',
        **kinetics_fields(config),
        'z_positions': z_positions.tolist(),
        'conversion_profile': conversion_profile,
        'temperature_profile': temperature_profile,
        'theta_C_profile': [0.0] * len(conversion_profile),
        **_kinetics_evidence(config),
        **_policy_metadata(config),
        **_load_status_fields(config, graphite, surf),
    }
    logger.info(
        f"  MMBCR result: conversion={final_conv:.2%}, "
        f"H2_atom_balance={result['H2_atom_balance']:.2%}, "
        f"τ={tau_total:.1f}s (u_b={hydro['bubble_rise_velocity_m_s']:.3f} m/s, "
        f"eps_g={hydro['gas_holdup_fraction']:.3f})"
    )
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# B. PFR (one shared surface marched through stages = time-on-stream)
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_pfr(config: ReactorConfig) -> Dict:
    """Packed bed as a staged Lagrangian PFR with one shared surface.

    Surface area per stage is ``sv·V_bed_stage`` where ``sv = 6(1−ε)/d_p``
    is area per *bed* volume (not ε-scaled), times loading × dispersion.
    The gas parcel volume is ``ε·V_bed_stage``. Stages are isothermal at
    ``T_inlet``; CH4 conversion is the Ar-tracer extent.

    Args:
        config: Configuration controlling this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    _validate_carbon_policy(config)
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'PFR')

    logger.info(f"Simulating PFR: {config.catalyst_name} at {config.T_inlet_K} K")
    gas, graphite, surf = _load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    eps = config.bed_porosity
    sv_ratio = active_sv(geometric_sv_pfr(config), config)  # m² active / m³ bed
    ch4_initial = _species_x(gas, 'CH4') or 1.0
    ar_initial = _species_x(gas, 'Ar')
    _require_ar_tracer(ch4_initial, ar_initial)

    n_stages = 50
    bed_cross_area = np.pi * (config.bed_diameter_m / 2) ** 2
    stage_length = config.bed_length_m / n_stages
    stage_bed_volume = bed_cross_area * stage_length
    stage_gas_volume = stage_bed_volume * eps
    stage_surface_area = sv_ratio * stage_bed_volume
    u_sup = config.gas_velocity_m_s if config.gas_velocity_m_s > 0 else 0.1
    tau_total = config.bed_length_m * eps / u_sup
    tau_stage = tau_total / n_stages

    z_positions = np.linspace(0, config.bed_length_m, n_stages + 1)
    conversion_profile = [0.0]
    temperature_profile = [config.T_inlet_K]
    theta_C_tos = [_coverage(surf, 'C_s')]

    cycles_completed = 0
    per_cycle_conversion: List[float] = []
    produce_time_s = 0.0
    pass_start_cov: List[Optional[np.ndarray]] = [None]

    def _advance_bed():
        nonlocal produce_time_s
        conversion_profile.clear()
        conversion_profile.append(0.0)
        temperature_profile.clear()
        temperature_profile.append(config.T_inlet_K)
        theta_C_tos.clear()
        theta_C_tos.append(_coverage(surf, 'C_s'))
        pass_start_cov[0] = (np.array(surf.coverages, dtype=float)
                             if surf is not None else None)
        for _ in range(n_stages):
            reactor = ct.IdealGasReactor(gas)
            reactor.volume = stage_gas_volume
            _disable_reactor_energy(reactor)
            if surf is not None:
                ct.ReactorSurface(surf, reactor, A=stage_surface_area)
            net = ct.ReactorNet([reactor])
            net.advance(tau_stage)
            gas.TPX = reactor.thermo.T, reactor.thermo.P, reactor.thermo.X
            produce_time_s += tau_stage
            conversion_profile.append(_ch4_extent(gas, ch4_initial, ar_initial))
            temperature_profile.append(float(gas.T))
            theta_C_tos.append(_coverage(surf, 'C_s'))

    # One produce pass (always). Optional discrete regen cycles if configured.
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    if surf is not None:
        _reset_surface_carbon(surf)
    _advance_bed()
    per_cycle_conversion.append(conversion_profile[-1])

    while (config.max_regen_cycles > 0
           and cycles_completed < config.max_regen_cycles
           and theta_C_tos and max(theta_C_tos) >= config.regen_coverage_threshold):
        if config.regen_mechanism == REGEN_OXIDATIVE and not config.co2_permitted:
            raise RuntimeError('oxidative regen blocked (co2_permitted=False)')
        if config.regen_mechanism == REGEN_OXIDATIVE:
            logger.warning('Oxidative PFR regen enabled via co2_permitted=True (test only)')
        # Mechanical / consumable: free-site reset without CO2 chemistry in-model.
        _reset_surface_carbon(surf)
        gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
        cycles_completed += 1
        _advance_bed()
        per_cycle_conversion.append(conversion_profile[-1])

    final_conv = conversion_profile[-1]
    x_h2 = _species_x(gas, 'H2')

    # Carbon accounting on the parcel basis (last pass).
    gamma_mol_m2 = float(config.site_density_mol_cm2) * 1e4
    n_sites_mol = gamma_mol_m2 * stage_surface_area
    c_total = config.P_inlet_Pa / (R_J_MOL_K * config.T_inlet_K)  # mol/m³
    n_parcel_in_mol = c_total * stage_gas_volume
    n_ch4_fed_mol = ch4_initial * n_parcel_in_mol
    off_site = _off_site_carbon_metrics(
        config, gas, surf,
        ch4_initial=ch4_initial, ar_initial=ar_initial,
        pass_conversion=float(final_conv), n_sites_mol=n_sites_mol,
        n_ch4_fed_mol=n_ch4_fed_mol, n_parcel_in_mol=n_parcel_in_mol,
        pass_time_s=tau_total, cov_start=pass_start_cov[0])

    result = {
        'reactor_type': 'PFR',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'bed_length_m': config.bed_length_m,
        'bed_diameter_m': config.bed_diameter_m,
        'catalyst_particle_mm': config.catalyst_particle_mm,
        'residence_time_s': tau_total,
        'produce_time_s': produce_time_s,
        'WHSV_h-1': reciprocal_residence_h(tau_total),
        'CH4_conversion': float(final_conv),
        'single_pass_CH4_conversion': float(
            per_cycle_conversion[0] if per_cycle_conversion else final_conv),
        'CH4_mole_fraction_drop': _mole_fraction_drop(gas, ch4_initial),
        'conversion_basis': 'argon_tracer',
        'thermal_mode': 'isothermal_energy_disabled',
        'surface_area_basis': 'sv_per_bed_volume_times_bed_volume',
        'stage_surface_area_m2': float(stage_surface_area),
        'stage_gas_volume_m3': float(stage_gas_volume),
        'per_cycle_CH4_conversion': per_cycle_conversion,
        'regen_cycles_completed': cycles_completed,
        'exit_x_H2': x_h2,
        'exit_T_K': float(gas.T),
        'z_positions': z_positions.tolist(),
        'conversion_profile': conversion_profile,
        'temperature_profile': temperature_profile,
        'theta_C_time_on_stream': theta_C_tos,
        'theta_C_profile_basis': 'single_shared_surface_cumulative_time_on_stream',
        'theta_C_axial': theta_C_tos,
        'inlet_theta_C': float(theta_C_tos[1] if len(theta_C_tos) > 1 else 0.0),
        'max_theta_C': float(max(theta_C_tos) if theta_C_tos else 0.0),
        'exit_theta_C': float(theta_C_tos[-1]) if theta_C_tos else None,
        **off_site,
        **kinetics_fields(config),
        **solids_inventory_fields(config, geometric_sv_pfr(config)),
        **_kinetics_evidence(config),
        **_policy_metadata(config),
        **_load_status_fields(config, graphite, surf),
    }
    logger.info(f"  PFR result: conversion={final_conv:.2%}, τ={tau_total:.1f}s, "
                f"regen_cycles={cycles_completed}")
    return result


def _integrate_fluidized_pass(gas, surf, tau: float, sv_ratio: float,
                              removal_rate_1_s: float) -> float:
    """Advance emulsion residence with C_s removal *during* integrate (B2).

    Basis: 1 m³ of emulsion. Gas volume is ε_mf; solids area is
    ``sv_ratio`` (area per emulsion volume). Pairing per-emulsion-volume
    area with 1 m³ of gas would undercount area by 1/ε_mf.

    Each substep builds a new ``ReactorNet`` from the current Interface
    state. Mutating ``Interface.coverages`` on a live ``ReactorSurface``
    does not change the integrated surface; the next ``advance`` would
    restore the pre-removal coverages and report outfeed that never left
    the reactor.
    """
    carbon_removed = 0.0
    n = max(1, int(FLUIDIZED_REMOVAL_SUBSTEPS))
    dt = tau / n
    for _ in range(n):
        reactor_em = ct.IdealGasReactor(gas)
        reactor_em.volume = FLUIDIZED_EMULSION_VOIDAGE
        _disable_reactor_energy(reactor_em)
        if surf is not None:
            ct.ReactorSurface(surf, reactor_em, A=sv_ratio)
        net = ct.ReactorNet([reactor_em])
        net.advance(dt)
        gas.TPX = reactor_em.thermo.T, reactor_em.thermo.P, reactor_em.thermo.X
        if removal_rate_1_s > 0 and surf is not None:
            carbon_removed += _apply_continuous_carbon_removal(
                surf, removal_rate_1_s, dt)
    return carbon_removed


# ═══════════════════════════════════════════════════════════════════════════════
# C. Fluidized bed
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_fluidized_bed(config: ReactorConfig) -> Dict:
    """Two-phase fluidized bed: reacting emulsion plus bubble bypass.

    The emulsion phase reacts at the minimum-fluidization residence time with
    in-step C_s removal (circulating) or batch regen. A bubble fraction δ of
    the feed bypasses unreacted and is mixed with the emulsion outlet on molar
    flows via the Ar tracer, so ``CH4_conversion = (1−δ)·X_emulsion``.
    Interphase mass transfer is not resolved; this is a screening approximation.

    Args:
        config: Configuration controlling this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    _validate_carbon_policy(config)
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'Fluidized')

    logger.info(
        f"Simulating Fluidized ({config.fluidized_mode}): "
        f"{config.catalyst_name} at {config.T_inlet_K} K"
    )
    gas, graphite, surf = _load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    inlet_x = np.array(gas.X, dtype=float)
    ch4_initial = _species_x(gas, 'CH4') or 1.0
    ar_initial = _species_x(gas, 'Ar')
    _require_ar_tracer(ch4_initial, ar_initial)

    u0 = max(config.gas_velocity_m_s, 0.05)
    umf = config.u_mf_m_s
    delta = (float(config.fluidized_bubble_fraction)
             if config.fluidized_bubble_fraction is not None
             else min(0.5, max(0.01, (u0 - umf) / u0)))
    if not 0 <= delta < 1:
        raise ValueError('fluidized-bed bubble fraction must lie in [0, 1)')
    tau_emulsion = config.bed_height_m * (1 - delta) / umf

    sv_ratio = active_sv(geometric_sv_fluidized(config), config)
    circulating = config.fluidized_mode == FLUIDIZED_CIRCULATING
    removal_rate = config.circulating_carbon_removal_rate_1_s if circulating else 0.0
    if surf is not None:
        _reset_surface_carbon(surf)
    pass_start_cov = (np.array(surf.coverages, dtype=float)
                      if surf is not None else None)
    carbon_removed = _integrate_fluidized_pass(
        gas, surf, tau_emulsion, sv_ratio, removal_rate)
    last_pass_removed = carbon_removed

    regen_cycles = 0
    per_cycle = []

    if not circulating:
        theta = _coverage(surf, 'C_s')
        per_cycle.append(_ch4_extent(gas, ch4_initial, ar_initial))
        while (config.max_regen_cycles > 0
               and regen_cycles < config.max_regen_cycles
               and theta >= config.regen_coverage_threshold):
            if config.regen_mechanism == REGEN_OXIDATIVE and not config.co2_permitted:
                raise RuntimeError('oxidative regen blocked (co2_permitted=False)')
            _reset_surface_carbon(surf)
            gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
            pass_start_cov = np.array(surf.coverages, dtype=float)
            last_pass_removed = _integrate_fluidized_pass(
                gas, surf, tau_emulsion, sv_ratio, 0.0)
            carbon_removed += last_pass_removed
            theta = _coverage(surf, 'C_s')
            regen_cycles += 1
            per_cycle.append(_ch4_extent(gas, ch4_initial, ar_initial))

    emulsion_conv = _ch4_extent(gas, ch4_initial, ar_initial)
    exit_theta_c = _coverage(surf, 'C_s')

    # B5/B6 carbon accounting on the emulsion basis (1 m³ emulsion: gas
    # volume ε_mf, area sv_ratio), before the bubble bypass is mixed in.
    # Carbon taken off by the circulating outfeed is B2 solids product,
    # not Cγ, and is subtracted from the filament balance.
    gamma_mol_m2 = float(config.site_density_mol_cm2) * 1e4
    n_sites_mol = gamma_mol_m2 * sv_ratio
    c_total = config.P_inlet_Pa / (R_J_MOL_K * config.T_inlet_K)
    n_parcel_in_mol = c_total * FLUIDIZED_EMULSION_VOIDAGE
    n_ch4_fed_mol = ch4_initial * n_parcel_in_mol
    off_site = _off_site_carbon_metrics(
        config, gas, surf,
        ch4_initial=ch4_initial, ar_initial=ar_initial,
        pass_conversion=float(emulsion_conv), n_sites_mol=n_sites_mol,
        n_ch4_fed_mol=n_ch4_fed_mol, n_parcel_in_mol=n_parcel_in_mol,
        pass_time_s=tau_emulsion, cov_start=pass_start_cov,
        removed_surface_carbon_mol=last_pass_removed * n_sites_mol,
        basis='Gamma*sv_emulsion / (c_CH4*eps_mf), emulsion parcel before bypass')

    # Bubble bypass: δ of the feed passes unreacted; mix on molar flows (Ar tracer).
    _mix_bubble_bypass(gas, inlet_x, delta, ar_initial)
    final_conv = _ch4_extent(gas, ch4_initial, ar_initial)
    expected = (1.0 - delta) * emulsion_conv
    if abs(final_conv - expected) > 1e-6 + 1e-6 * abs(expected):
        raise RuntimeError(
            f'bypass mixing inconsistent: X_mixed={final_conv:.6f} vs '
            f'(1-delta)*X_em={expected:.6f}')
    x_h2 = _species_x(gas, 'H2')

    result = {
        'reactor_type': 'Fluidized',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'catalyst_particle_mm': config.catalyst_particle_mm,
        'bed_height_m': config.bed_height_m,
        'u0_m_s': u0,
        'umf_m_s': umf,
        'bubble_fraction': delta,
        'bubble_fraction_basis': (
            'override' if config.fluidized_bubble_fraction is not None
            else 'clip((u0-umf)/u0, 0.01, 0.5)'),
        'residence_time_s': tau_emulsion,
        'WHSV_h-1': reciprocal_residence_h(tau_emulsion),
        'emulsion_voidage': FLUIDIZED_EMULSION_VOIDAGE,
        'surface_area_basis': 'sv_per_emulsion_volume_times_emulsion_volume',
        'emulsion_CH4_conversion': float(emulsion_conv),
        'CH4_conversion': float(final_conv),
        'single_pass_CH4_conversion': float(final_conv),
        'bypass_basis': 'ar_tracer_molar_mix',
        'CH4_mole_fraction_drop': _mole_fraction_drop(gas, ch4_initial),
        'conversion_basis': 'argon_tracer',
        'thermal_mode': 'isothermal_energy_disabled',
        'exit_x_H2': float(x_h2),
        'exit_T_K': float(gas.T),
        'exit_theta_C': exit_theta_c,
        'carbon_removed_coverage_proxy': float(carbon_removed),
        'regen_cycles_completed': regen_cycles,
        'per_cycle_CH4_conversion': per_cycle,
        'per_cycle_basis': 'emulsion_only_before_bypass_mix',
        **off_site,
        **kinetics_fields(config),
        **solids_inventory_fields(config, geometric_sv_fluidized(config)),
        **_kinetics_evidence(config),
        **_policy_metadata(config),
        **_load_status_fields(config, graphite, surf),
    }
    logger.info(
        f"  Fluidized result: conversion={final_conv:.2%} "
        f"(emulsion {emulsion_conv:.2%}, delta={delta:.2f}) mode={config.fluidized_mode}")
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# MOCK RESULTS (for testing without Cantera)
# ═══════════════════════════════════════════════════════════════════════════════

def _mock_reactor_result(config: ReactorConfig, reactor_type: str) -> Dict:
    logger.warning(f"Cantera not available. Generating mock {reactor_type} results.")
    _validate_carbon_policy(config)
    E_act = config.catalyst_E_act_eV
    k_B_eV = 8.617e-5
    k = 1e13 * np.exp(-E_act / (k_B_eV * config.T_inlet_K))
    extra: Dict = {}
    if reactor_type == 'MMBCR':
        hydro = _mmbcr_hydrodynamics(config)
        tau = hydro['residence_time_s']
        extra = hydro
    elif reactor_type == 'PFR':
        tau = config.bed_length_m * config.bed_porosity / max(config.gas_velocity_m_s, 0.1)
    else:
        tau = config.bed_height_m / max(config.gas_velocity_m_s, 0.05)
    conversion = float(np.clip(1.0 - np.exp(-k * tau * 1e-12), 0.01, 0.99))
    return {
        'reactor_type': reactor_type,
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'catalyst_E_act_eV': E_act,
        'residence_time_s': tau,
        **extra,
        'WHSV_h-1': reciprocal_residence_h(tau),
        'CH4_conversion': conversion,
        'single_pass_CH4_conversion': conversion,
        'conversion_basis': (
            'melt_ode_to_Xeq' if reactor_type == 'MMBCR' else 'argon_tracer'),
        'thermal_mode': (
            'isothermal_by_construction' if reactor_type == 'MMBCR'
            else 'isothermal_energy_disabled'),
        'exit_T_K': config.T_inlet_K,
        **kinetics_fields(config),
        'H2_atom_balance': 0.95,
        'H2_selectivity': 0.95,
        'solid_C_selectivity': 0.90,
        'exit_x_H2': conversion * 0.95 * 2.0 / (1.0 + conversion * 0.95),
        'mock': True,
        **_kinetics_evidence(config),
        **_policy_metadata(config),
        **_load_status_fields(config, None, None),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED SIMULATION INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_ntec_pathway(config: ReactorConfig) -> Dict:
    """Describe NTEC readiness without substituting an unrelated reactor model.

    Args:
        config: Configuration controlling this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    from pipeline.process.ntec_model import (
        conditions_from_environment, ntec_assistance)

    assistance = ntec_assistance(conditions_from_environment())
    if config.multiphysics_artifact:
        evidence = config.multiphysics_artifact
        artifact = evidence['artifact']
        outputs = artifact['outputs']
        return {
            'status': 'complete', 'valid': True,
            'reactor_type': 'NTEC', 'pathway_mode': 'ntec',
            'catalyst_name': config.catalyst_name,
            'candidate_id': config.candidate_id,
            'T_K': config.T_inlet_K,
            'CH4_conversion': float(outputs['CH4_conversion']),
            'H2_selectivity': float(outputs['H2_selectivity']),
            'solid_C_selectivity': float(outputs['solid_C_selectivity']),
            'specific_energy_kWh_kg_H2': float(
                outputs['specific_energy_kWh_kg_H2']),
            'ntec_assistance': assistance,
            'multiphysics_evidence': evidence,
            'reactor_evidence_tier': 'calibrated_multiphysics_screening',
            'can_exclude_candidate': False,
        }
    return {
        'status': 'validation_required',
        'valid': False,
        'reactor_type': 'NTEC',
        'pathway_mode': 'ntec',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'ntec_assistance': assistance,
        'reactor_evidence_tier': 'pathway_model_pending',
        'can_exclude_candidate': False,
        'limitations': [
            'candidate_specific_ntec_pathway_kinetics_required',
            'liquid_solid_hydrodynamic_model_required',
        ],
    }


def simulate_electrochemical_pathway(config: ReactorConfig) -> Dict:
    """Report electrochemical evidence without inventing a Cantera conversion.

    Args:
        config: Configuration controlling this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    from pipeline.process.electrochemical_model import (
        conditions_from_environment, electrochemical_evidence)

    evidence = electrochemical_evidence(conditions_from_environment())
    phase = evidence['conditions'].get('electrolyte_phase')
    if config.multiphysics_artifact:
        solver_evidence = config.multiphysics_artifact
        artifact = solver_evidence['artifact']
        outputs = artifact['outputs']
        phase = artifact.get('electrolyte_phase', phase)
        return {
            'status': 'complete', 'valid': True,
            'reactor_type': 'Electrochemical',
            'pathway_mode': 'electrochemical',
            'electrolyte_phase': phase,
            'catalyst_name': config.catalyst_name,
            'candidate_id': config.candidate_id,
            'T_K': config.T_inlet_K,
            **{name: float(outputs[name]) for name in (
                'CH4_conversion', 'H2_selectivity',
                'faradaic_efficiency_H2', 'current_density_A_cm2',
                'cell_voltage_V', 'electrical_power_density_W_cm2')},
            'electrochemical_evidence': evidence,
            'multiphysics_evidence': solver_evidence,
            'reactor_evidence_tier': 'mechanistic_multiphysics_screening',
            'can_exclude_candidate': False,
        }
    return {
        'status': 'validation_required',
        'valid': False,
        'reactor_type': 'Electrochemical',
        'pathway_mode': 'electrochemical',
        'electrolyte_phase': phase,
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'electrochemical_evidence': evidence,
        'reactor_evidence_tier': 'pathway_model_pending',
        'can_exclude_candidate': False,
        'limitations': [
            'candidate_specific_electrochemical_kinetics_required',
            f'{phase or "unspecified"}_electrolyte_transport_model_required',
        ],
    }


def simulate_reactor(config: ReactorConfig, coupling_services=None) -> Dict:
    """Run the appropriate reactor simulation based on config.reactor_type.

    Args:
        config: Configuration controlling this operation.
        coupling_services: Coupling services used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    simulators = {
        'MMBCR': simulate_mmbcr,
        'PFR': simulate_pfr,
        'Fluidized': simulate_fluidized_bed,
        'NTEC': simulate_ntec_pathway,
        'Electrochemical': simulate_electrochemical_pathway,
    }
    if config.reactor_type not in simulators:
        raise ValueError(f'unsupported reactor type: {config.reactor_type}')
    selected_mode = reactor_types_for_mode(config.pathway_mode)
    if config.reactor_type not in selected_mode:
        raise ValueError(
            f'reactor {config.reactor_type!r} is not valid for pathway mode '
            f'{config.pathway_mode!r}; expected {selected_mode}')
    applicable, reason = reactor_applicability(
        config.reactor_type, config.material_class)
    spec = REACTOR_MODELS[config.reactor_type]
    if not applicable:
        result = {
            'status': 'not_applicable', 'valid': False,
            'reactor_type': config.reactor_type,
            'pathway_mode': config.pathway_mode,
            'catalyst_name': config.catalyst_name,
            'T_K': config.T_inlet_K,
            'material_class': config.material_class,
            'reason': reason, 'can_exclude_candidate': False,
            'reactor_evidence_tier': 'not_applicable_non_excluding',
            'bed_or_interface': spec.bed_or_interface,
            'cantera_reactor_model': spec.cantera_model,
            'reaction_domain': spec.reaction_domain,
            'reactor_model_fidelity': spec.fidelity,
        }
    else:
        from pipeline.process.multiphysics_contract import EXTERNAL_SOLVERS
        if config.reactor_type in EXTERNAL_SOLVERS:
            from pipeline.process.reactor_coupling import (
                default_reactor_coupling_services)
            coupling_services = (coupling_services or
                                 default_reactor_coupling_services())
            loaded = coupling_services.load_artifact(
                config.multiphysics_results_dir, config.candidate_id,
                config.pathway_mode, config.reactor_type, config.T_inlet_K)
            loaded = coupling_services.validate_compatibility(config, loaded)
            if not loaded['valid']:
                from pipeline.process.closure_provider import (
                    resolve_reactor_closure)
                surrogate = (coupling_services.load_surrogate(
                    config.pathway_mode, config.reactor_type)
                    if coupling_services.load_surrogate else None)
                closure = resolve_reactor_closure(
                    full_physics=loaded, surrogate=surrogate,
                    features=config.closure_features,
                    pathway_mode=config.pathway_mode,
                    reactor_type=config.reactor_type,
                    temperature_K=config.T_inlet_K)
                if closure['available']:
                    from pipeline.process.reactor_coupling import (
                        couple_surrogate_closure)
                    couple_surrogate_closure(config, closure)
                else:
                    loaded = {**loaded, 'closure_resolution': closure}
            if (not loaded['valid'] and not config.reactor_closure_evidence
                    and config.reactor_type in ANALYTICAL_CLOSURE_REACTORS):
                # No validated artifact or calibrated surrogate. Fluidized and
                # MMBCR carry an explicit analytical bubbly-flow closure
                # (Mendelson u_b / derived eps_g; clipped (u0-umf)/u0) so the
                # screening run can proceed with its provenance labelled.
                # It can never exclude a candidate. NTEC / Electrochemical
                # have no such closure and stay validation_required.
                config.reactor_closure_evidence = {
                    'source': ANALYTICAL_CLOSURE_SOURCE,
                    'candidate_exclusion_authorized': False,
                    'basis': ANALYTICAL_CLOSURE_BASIS[config.reactor_type],
                    'multiphysics_evidence': loaded,
                }
            if not loaded['valid'] and not config.reactor_closure_evidence:
                result = {
                    'status': 'validation_required', 'valid': False,
                    'reactor_type': config.reactor_type,
                    'pathway_mode': config.pathway_mode,
                    'catalyst_name': config.catalyst_name,
                    'T_K': config.T_inlet_K,
                    'candidate_id': config.candidate_id,
                    'material_class': config.material_class,
                    'multiphysics_evidence': loaded,
                    'can_exclude_candidate': False,
                    'reactor_evidence_tier': 'required_solver_evidence_missing',
                    'bed_or_interface': spec.bed_or_interface,
                    'cantera_reactor_model': spec.cantera_model,
                    'reaction_domain': spec.reaction_domain,
                    'reactor_model_fidelity': spec.fidelity,
                }
                REACTOR_DIR.mkdir(parents=True, exist_ok=True)
                fname = (f"{config.reactor_type}_{config.catalyst_name}_"
                         f"{int(config.T_inlet_K)}K.json")
                save_json(result, fname, subdir='reactor')
                return result
            if loaded['valid']:
                coupling_services.couple_evidence(config, loaded)
        _validate_reactor_config(config)
        _validate_carbon_policy(config)
        result = simulators[config.reactor_type](config)
        result.setdefault('status', 'complete')
        result.update({
            'pathway_mode': config.pathway_mode,
            'material_class': config.material_class,
            'bed_or_interface': spec.bed_or_interface,
            'cantera_reactor_model': spec.cantera_model,
            'reaction_domain': spec.reaction_domain,
            'reactor_model_fidelity': spec.fidelity,
            'multiphysics_evidence': config.multiphysics_artifact,
            'reactor_closure_evidence': config.reactor_closure_evidence,
        })
        if spec.fidelity != 'validated_predictive':
            limitations = list(result.get('reactor_evidence_limitations', []))
            limitation = f'{config.reactor_type.lower()}_screening_model_not_validated'
            if limitation not in limitations:
                limitations.append(limitation)
            result['reactor_evidence_limitations'] = limitations
            result['can_exclude_candidate'] = False

    # Save result
    REACTOR_DIR.mkdir(parents=True, exist_ok=True)
    fname = f"{config.reactor_type}_{config.catalyst_name}_{int(config.T_inlet_K)}K.json"
    save_json(result, fname, subdir='reactor')

    return result


def run_reactor_sweep(catalyst_name: str, mechanism_file: str,
                      temperatures: List[float] = None,
                      reactor_types: List[str] = None,
                      catalyst_E_act_eV: float = 0.8,
                      pathway_mode: str = DEFAULT_MODE,
                      material_class: Optional[str] = None,
                      candidate_id: str = 'unknown',
                      multiphysics_results_dir: Optional[str] = None,
                      coupling_services=None,
                      catalyst_dE_H_eV: float = 0.0,
                      reactor_config_kwargs: Optional[Dict] = None) -> List[Dict]:
    """Sweep operating conditions across the reactors owned by a pathway.

    Args:
        catalyst_name: Human-readable catalyst identifier used in output names.
        mechanism_file: Candidate-specific Cantera mechanism path.
        temperatures: Absolute operating temperatures in kelvin; defaults to the
            standard four-point screening sweep.
        reactor_types: Ordered reactor implementations; defaults to those owned
            by ``pathway_mode``.
        catalyst_E_act_eV: Activation barrier in eV. Used by the MMBCR melt
            k_if and by the explicitly labeled mock path.
        pathway_mode: Methane-conversion pathway controlling valid reactors.
        material_class: Catalyst class used for physical-bed compatibility.
        candidate_id: Stable identity required by multiphysics evidence.
        multiphysics_results_dir: Root containing validated external artifacts.
        coupling_services: Optional injected evidence-loading/coupling adapters.
        catalyst_dE_H_eV: Screening H adsorption energy, carried into results.
        reactor_config_kwargs: Extra ``ReactorConfig`` fields (carbon policy,
            solids inventory, melt properties) applied to every condition.

    Returns:
        One result record per requested reactor and temperature. Individual
        condition failures remain non-excluding records rather than aborting the
        complete sweep.
    """
    if temperatures is None:
        temperatures = [773.15, 900.0, 1100.0, 1300.0]
    if reactor_types is None:
        # The default family is conventional thermocatalysis. MMBCR, NTEC, and
        # electrochemical physics must be selected explicitly.
        reactor_types = list(reactor_types_for_mode(pathway_mode))
    validate_mode_reactors(pathway_mode, reactor_types)
    extra = dict(reactor_config_kwargs or {})

    results = []
    for rt in reactor_types:
        for T in temperatures:
            config = ReactorConfig(
                T_inlet_K=T,
                reactor_type=rt,
                mechanism_file=str(mechanism_file),
                catalyst_name=catalyst_name,
                catalyst_E_act_eV=catalyst_E_act_eV,
                catalyst_dE_H_eV=catalyst_dE_H_eV,
                pathway_mode=pathway_mode,
                material_class=material_class,
                candidate_id=candidate_id,
                multiphysics_results_dir=multiphysics_results_dir,
                **extra,
            )
            try:
                result = (simulate_reactor(config) if coupling_services is None
                          else simulate_reactor(
                              config, coupling_services=coupling_services))
            except Exception as exc:
                # A stiff condition must not discard the other temperatures or
                # reactor types. Preserve it as non-excluding failed evidence.
                result = {
                    'status': 'failed',
                    'valid': False,
                    'reactor_type': rt,
                    'pathway_mode': pathway_mode,
                    'material_class': material_class,
                    'temperature_K': float(T),
                    'catalyst': catalyst_name,
                    'error_type': type(exc).__name__,
                    'error': str(exc),
                    'can_exclude_candidate': False,
                    'reactor_evidence_tier': 'failed_simulation',
                }
                REACTOR_DIR.mkdir(parents=True, exist_ok=True)
                fname = f"{rt}_{catalyst_name}_{int(T)}K.json"
                save_json(result, fname, subdir='reactor')
                logger.error(
                    'Reactor condition failed without excluding candidate: '
                    f'{rt} {catalyst_name} {T} K: {exc}')
            results.append(result)

    return results


if __name__ == '__main__':
    print_banner("REACTOR SIMULATION TEST")

    # Test with mock data (no Cantera needed)
    config = ReactorConfig(
        T_inlet_K=1000.0,
        reactor_type='MMBCR',
        catalyst_name='NiBi_test',
        material_class='MoltenMetal',
        pathway_mode='mmbcr',
    )
    result = simulate_reactor(config)
    print(json.dumps(result, indent=2))
