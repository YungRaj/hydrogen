"""Shared configuration and physical constants for reactor implementations.

This module contains the vocabulary common to every reactor family.  It does
not select a reactor or execute a solver; those responsibilities remain with
the family implementations and the compatibility router in ``models``.
"""

from dataclasses import dataclass
from pathlib import Path
import json
from typing import Dict, Optional

import numpy as np

try:
    import cantera as ct
    HAS_CANTERA = True
except ImportError:
    HAS_CANTERA = False

from pipeline.reactors.mechanisms import MONOLAYER_SITE_DENSITY_MOL_CM2
from pipeline.reactors.modes import DEFAULT_MODE
from pipeline.utils import setup_logger

logger = setup_logger('reactor_core', 'reactor/reactor_simulation.log')


REGEN_MECHANICAL = 'mechanical'
REGEN_CONSUMABLE = 'consumable'
REGEN_OXIDATIVE = 'oxidative'
ALLOWED_REGEN = frozenset({REGEN_MECHANICAL, REGEN_CONSUMABLE, REGEN_OXIDATIVE})
FLUIDIZED_BATCH = 'batch_regen'
FLUIDIZED_CIRCULATING = 'circulating'
ALLOWED_FLUIDIZED = frozenset({FLUIDIZED_BATCH, FLUIDIZED_CIRCULATING})

ERGUN_DP_LIMIT_PA = 101325.0
MMBCR_INTERFACIAL_K0_DEFAULT = 0.01
MMBCR_INTERFACIAL_K0_MAX = 1.0
MMBCR_INTERFACIAL_K0_MAX_KIND = 'guardrail_not_physical_bound'
MMBCR_FLOTATION_UNCONSTRAINED = None
MMBCR_MELT_SURFACE_TENSION_N_M = 0.4
MMBCR_MELT_DENSITY_KG_M3 = 9000.0
MMBCR_MAX_GAS_HOLDUP = 0.3
MMBCR_BUBBLE_RISE_BASIS = 'mendelson_sqrt(2*sigma/(rho*d_b) + g*d_b/2)'
MMBCR_RESIDENCE_BASIS = 'column_height_over_bubble_rise_velocity'
MMBCR_ARTIFACT_RESIDENCE_BASIS = (
    'gas_holdup_times_column_height_over_superficial_velocity')
G_M_S2 = 9.80665
R_J_MOL_K = 8.314462618
FLUIDIZED_REMOVAL_SUBSTEPS = 20
FLUIDIZED_EMULSION_VOIDAGE = 0.45

ANALYTICAL_CLOSURE_SOURCE = 'analytical_hydrodynamic_closure'
ANALYTICAL_CLOSURE_REACTORS = frozenset({'Fluidized', 'MMBCR'})
ANALYTICAL_CLOSURE_BASIS = {
    'Fluidized': 'bubble_fraction = clip((u0-umf)/u0, 0.01, 0.5); emulsion at umf',
    'MMBCR': MMBCR_BUBBLE_RISE_BASIS + '; eps_g = u_sup/u_b; tau = H/u_b',
}
EXTERNAL_CLOSURE_SOURCES = frozenset({
    'validated_full_physics', 'calibrated_transport_surrogate'})

SINGLE_REACTOR_MODE = {
    'PFR': 'thermocatalytic_pfr',
    'Fluidized': 'thermocatalytic_fluidized',
    'MMBCR': 'mmbcr',
}
DIAGNOSTIC_MATERIAL_CLASS = {
    'PFR': 'SolidCatalyst',
    'Fluidized': 'SolidCatalyst',
    'MMBCR': 'MoltenMetal',
}

DEFAULT_SOLIDS_PARTICLE_MM = 0.13
DEFAULT_METAL_LOADING = 0.5
DEFAULT_METAL_DISPERSION = 0.3

INVENTORY_PARTICLE_MM = (2.0, 0.5, 0.2, 0.1)
INVENTORY_METAL_LOADING = (1.0, 0.5, 0.2)
INVENTORY_METAL_DISPERSION = (1.0, 0.3, 0.1)
INVENTORY_ROI_PARTICLE_MM = (0.25, 0.20, 0.16, 0.13, 0.10, 0.08)
INVENTORY_ROI_METAL_LOADING = (1.0, 0.7, 0.5)
INVENTORY_ROI_METAL_DISPERSION = (1.0, 0.5, 0.3)
INVENTORY_ROI_REASON = (
    'coarse grid: X proportional to a; gain starts at d_p<=0.5 mm and is '
    'material at <=0.2 mm; Ergun wall ~0.08 mm on this 0.5 m / 0.05 m/s bed')


@dataclass
class ReactorConfig:
    """Configuration shared by routed reactor simulations."""

    T_inlet_K: float = 1000.0
    P_inlet_Pa: float = 101325.0
    inlet_composition: str = 'CH4:0.95, Ar:0.05'

    column_height_m: float = 1.5
    column_diameter_m: float = 0.10
    bubble_diameter_mm: float = 5.0
    gas_velocity_m_s: float = 0.05
    gas_holdup_fraction: Optional[float] = None
    n_cstr_stages: int = 20
    melt_surface_tension_N_m: float = MMBCR_MELT_SURFACE_TENSION_N_M
    melt_density_kg_m3: float = MMBCR_MELT_DENSITY_KG_M3

    bed_length_m: float = 0.5
    bed_diameter_m: float = 0.05
    catalyst_particle_mm: float = DEFAULT_SOLIDS_PARTICLE_MM
    bed_porosity: float = 0.4
    site_density_mol_cm2: float = MONOLAYER_SITE_DENSITY_MOL_CM2
    metal_loading: float = DEFAULT_METAL_LOADING
    metal_dispersion: float = DEFAULT_METAL_DISPERSION

    u_mf_m_s: float = 0.02
    bed_height_m: float = 0.8
    catalyst_density_kg_m3: float = 2500.0
    fluidized_bubble_fraction: Optional[float] = None

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
    catalyst_E_act_eV: float = 0.8
    catalyst_dE_H_eV: float = 0.0

    mmbcr_carbon_removal_rate_1_s: Optional[float] = (
        MMBCR_FLOTATION_UNCONSTRAINED)
    mmbcr_interfacial_k0_m_s: float = MMBCR_INTERFACIAL_K0_DEFAULT
    regen_coverage_threshold: float = 0.8
    regen_mechanism: str = REGEN_MECHANICAL
    max_regen_cycles: int = 3
    fluidized_mode: str = FLUIDIZED_CIRCULATING
    circulating_carbon_removal_rate_1_s: float = 0.5
    co2_permitted: bool = False
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
    """External pellet area per bed volume: 6(1−ε)/d_p.

    Args:
        config: Input controlling config.

    Returns:
        Validated float output for this operation.
    """
    d_p = config.catalyst_particle_mm * 1e-3
    if d_p <= 0:
        raise ValueError('catalyst_particle_mm must be positive')
    return 6.0 * (1.0 - config.bed_porosity) / d_p


def geometric_sv_fluidized(config: ReactorConfig) -> float:
    """Emulsion solids area per emulsion volume: 6·(1−ε_mf)/d_p.

    Args:
        config: Input controlling config.

    Returns:
        Validated float output for this operation.
    """
    d_p = config.catalyst_particle_mm * 1e-3
    if d_p <= 0:
        raise ValueError('catalyst_particle_mm must be positive')
    return 6.0 * (1.0 - FLUIDIZED_EMULSION_VOIDAGE) / d_p


def active_area_multiplier(config: ReactorConfig) -> float:
    """Return the loading-times-dispersion active-area fraction.

    Args:
        config: Input controlling config.

    Returns:
        Validated float output for this operation.
    """
    return float(config.metal_loading) * float(config.metal_dispersion)


def active_sv(geometric_sv: float, config: ReactorConfig) -> float:
    """Apply catalyst utilization to a geometric surface-to-volume ratio.

    Args:
        geometric_sv: Input controlling geometric sv.
        config: Input controlling config.

    Returns:
        Validated float output for this operation.
    """
    return float(geometric_sv) * active_area_multiplier(config)


def ch4_feed_density_kg_m3(T_K: float, P_Pa: float) -> float:
    """Ideal-gas density for CH4:0.95 / Ar:0.05.

    Args:
        T_K: Input controlling T K.
        P_Pa: Input controlling P Pa.

    Returns:
        Validated float output for this operation.
    """
    M = 0.01604 * 0.95 + 0.03995 * 0.05
    return P_Pa * M / (8.314462618 * T_K)


def ch4_viscosity_pa_s(T_K: float) -> float:
    """Sutherland estimate for CH4 (μ0=1.03e-5 Pa·s at 273.15 K, S=164 K).

    Args:
        T_K: Input controlling T K.

    Returns:
        Validated float output for this operation.
    """
    T0, mu0, S = 273.15, 1.03e-5, 164.0
    return mu0 * (T_K / T0) ** 1.5 * (T0 + S) / (T_K + S)


def ergun_delta_p_pa(config: ReactorConfig, T_K: float = None,
                     P_Pa: float = None) -> float:
    """Packed-bed Ergun ΔP over bed_length_m (B4). Not used for MMBCR.

    Args:
        config: Input controlling config.
        T_K: Input controlling T K.
        P_Pa: Input controlling P Pa.

    Returns:
        Validated float output for this operation.
    """
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
    """Space velocity as 1/τ in h⁻¹. Field name WHSV is historical.

    Args:
        tau_s: Input controlling tau s.

    Returns:
        Validated float output for this operation.
    """
    return 3600.0 / float(tau_s) if tau_s and float(tau_s) > 0 else 0.0


def kinetics_fields(config: ReactorConfig) -> Dict:
    """Read candidate kinetic provenance fields from mechanism metadata.

    Args:
        config: Input controlling config.

    Returns:
        Validated Dict output for this operation.
    """
    return {
        'catalyst_E_act_eV': float(config.catalyst_E_act_eV),
        'catalyst_dE_H_eV': float(config.catalyst_dE_H_eV),
    }


def solids_inventory_fields(config: ReactorConfig, geometric_sv: float) -> Dict:
    """Calculate catalyst inventory, area, WHSV, and pressure-drop fields.

    Args:
        config: Input controlling config.
        geometric_sv: Input controlling geometric sv.

    Returns:
        Validated Dict output for this operation.
    """
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
    """Quantized (d_p, loading, dispersion) grid. Defaults to the coarse B1-2 set.

    Args:
        particle_mm: Input controlling particle mm.
        loadings: Input controlling loadings.
        dispersions: Input controlling dispersions.
    """
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
    """Return the predefined region-of-interest inventory grid.
    """
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
            from pipeline.reactors.mechanisms import (
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
    from pipeline.reactors.equilibrium import TABULATED_X_CH4_1BAR
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
    """Mendelson (1967) terminal rise velocity: sqrt(2σ/(ρ d_b) + g d_b/2).

    Args:
        d_b_m: Input controlling d b m.
        sigma_N_m: Input controlling sigma N m.
        rho_kg_m3: Input controlling rho kg m3.

    Returns:
        Validated float output for this operation.
    """
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
    from pipeline.reactors.equilibrium import ch4_conversion_from_argon_tracer
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


# Export private compatibility helpers deliberately: models.py remains the
# stable import and monkeypatch boundary while implementations are decomposed.
_INFRASTRUCTURE_EXPORTS = {
    'Dict', 'Optional', 'Path', 'ct', 'dataclass', 'json', 'logger', 'np',
    'setup_logger', '_INFRASTRUCTURE_EXPORTS',
}
__all__ = [
    name for name in globals()
    if not name.startswith('__') and name not in _INFRASTRUCTURE_EXPORTS]
