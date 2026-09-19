#!/usr/bin/env python3
# Turquoise-hydrogen reactor process models.
"""
Cantera Reactor-Scale Simulation Models for Methane Pyrolysis.

Three reactor archetypes with distinct carbon-handling physics:
  A. MMBCR — continuous buoyant/transport carbon removal (steady; no site lattice claim)
  B. PFR — one shared surface marched through stages (time-on-stream, not axial);
     optional discrete non-oxidative regen
  C. Fluidized — explicit batch_regen vs circulating mode

Solid carbon is never a gas-phase species. Surface C_s blocks sites on solid
paths until removed by a named policy. Oxidative regen requires co2_permitted.
"""

import json
import sys
from dataclasses import dataclass, field
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
FLUIDIZED_REMOVAL_SUBSTEPS = 20

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


@dataclass
class ReactorConfig:
    """Configuration for reactor simulation."""
    T_inlet_K: float = 1000.0
    P_inlet_Pa: float = 101325.0
    inlet_composition: str = 'CH4:0.95, Ar:0.05'

    column_height_m: float = 1.5
    bubble_diameter_mm: float = 5.0
    gas_velocity_m_s: float = 0.05
    n_cstr_stages: int = 20

    bed_length_m: float = 0.5
    bed_diameter_m: float = 0.05
    catalyst_particle_mm: float = DEFAULT_SOLIDS_PARTICLE_MM
    bed_porosity: float = 0.4
    # Γ is a monolayer. Do not raise to force Da (B1).
    site_density_mol_cm2: float = MONOLAYER_SITE_DENSITY_MOL_CM2
    # Fraction of geometric pellet surface that is metal, and of that metal
    # that is surface-available. Defaults are the supported-TCD proxy
    # (B1-4; Alves 2021; Sánchez-Bastardo 2021; Gili 2024). Neither may
    # exceed 1 — extra area is not a BET/Γ invention (B1).
    metal_loading: float = DEFAULT_METAL_LOADING
    metal_dispersion: float = DEFAULT_METAL_DISPERSION

    u_mf_m_s: float = 0.02
    bed_height_m: float = 0.8
    catalyst_density_kg_m3: float = 2500.0

    reactor_type: str = 'MMBCR'
    mechanism_file: str = ''
    catalyst_name: str = 'test'
    max_residence_time_s: float = 60.0
    catalyst_E_act_eV: float = 0.8
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
        'carbon_phase_model': 'condensed_graphite_plus_surface_C_s',
        'mmbcr_rate_model': 'bubble_area_flotation',
        'mmbcr_interfacial_k0_m_s': float(config.mmbcr_interfacial_k0_m_s),
        'mmbcr_interfacial_k0_basis': (
            'calibrated melt-side prefactor; not DFT or Upham-fitted k; '
            'B3; Upham 2017; Chen 2023; Abdollahi 2024'),
        'mmbcr_interfacial_k0_max_m_s': MMBCR_INTERFACIAL_K0_MAX,
        'mmbcr_interfacial_k0_max_kind': MMBCR_INTERFACIAL_K0_MAX_KIND,
        'mmbcr_flotation_default': 'unconstrained (eta=1); knob wired but off',
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
    """Emulsion solids area per emulsion volume: 6×0.55/d_p."""
    d_p = config.catalyst_particle_mm * 1e-3
    if d_p <= 0:
        raise ValueError('catalyst_particle_mm must be positive')
    return 6.0 * 0.55 / d_p


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
    """Mechanical / consumable regen: clear C_s and restore free sites."""
    if surf is None or 'C_s' not in surf.species_names:
        return
    cov = np.zeros(surf.n_species)
    if 'site' in surf.species_names:
        cov[surf.species_index('site')] = 1.0
    surf.coverages = cov


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


# ═══════════════════════════════════════════════════════════════════════════════
# A. MMBCR — bubble area + carbon flotation (no solid site lattice)
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_mmbcr(config: ReactorConfig) -> Dict:
    _validate_carbon_policy(config)
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'MMBCR')

    logger.info(f"Simulating MMBCR: {config.catalyst_name} at {config.T_inlet_K} K")
    gas, graphite, surf = _load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    x_ch4_feed = _species_x(gas, 'CH4') or 0.95
    x_ar_feed = _species_x(gas, 'Ar')
    if x_ar_feed <= 0:
        x_ar_feed = max(0.0, 1.0 - x_ch4_feed)

    tau_total = config.column_height_m / config.gas_velocity_m_s
    tau_stage = tau_total / config.n_cstr_stages
    d_b = config.bubble_diameter_mm * 1e-3
    sv_ratio = 6.0 / d_b  # bubble S/V (m² interface / m³ bubble)
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
        'gas_velocity_m_s': config.gas_velocity_m_s,
        'bubble_diameter_mm': config.bubble_diameter_mm,
        'interfacial_sv_ratio_1_m': sv_ratio,
        'residence_time_s': tau_total,
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
        f"H2_atom_balance={result['H2_atom_balance']:.2%}, τ={tau_total:.1f}s"
    )
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# B. PFR (one shared surface marched through stages = time-on-stream)
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_pfr(config: ReactorConfig) -> Dict:
    _validate_carbon_policy(config)
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'PFR')

    logger.info(f"Simulating PFR: {config.catalyst_name} at {config.T_inlet_K} K")
    gas, graphite, surf = _load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    eps = config.bed_porosity
    sv_ratio = active_sv(geometric_sv_pfr(config), config)
    ch4_initial = _species_x(gas, 'CH4') or 1.0
    ar_initial = _species_x(gas, 'Ar')
    _require_ar_tracer(ch4_initial, ar_initial)

    n_stages = 50
    bed_cross_area = np.pi * (config.bed_diameter_m / 2) ** 2
    stage_length = config.bed_length_m / n_stages
    stage_volume = bed_cross_area * stage_length * eps
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

    def _advance_bed():
        nonlocal produce_time_s
        conversion_profile.clear()
        conversion_profile.append(0.0)
        temperature_profile.clear()
        temperature_profile.append(config.T_inlet_K)
        theta_C_tos.clear()
        theta_C_tos.append(_coverage(surf, 'C_s'))
        for _ in range(n_stages):
            reactor = ct.IdealGasReactor(gas)
            reactor.volume = stage_volume
            _disable_reactor_energy(reactor)
            if surf is not None:
                ct.ReactorSurface(surf, reactor, A=sv_ratio * stage_volume)
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
    """Advance emulsion residence with C_s removal *during* integrate (B2)."""
    reactor_em = ct.IdealGasReactor(gas)
    reactor_em.volume = 1.0
    _disable_reactor_energy(reactor_em)
    if surf is not None:
        ct.ReactorSurface(surf, reactor_em, A=sv_ratio)
    net = ct.ReactorNet([reactor_em])
    carbon_removed = 0.0
    n = max(1, int(FLUIDIZED_REMOVAL_SUBSTEPS))
    dt = tau / n
    t = 0.0
    for _ in range(n):
        t += dt
        net.advance(t)
        if removal_rate_1_s > 0 and surf is not None:
            carbon_removed += _apply_continuous_carbon_removal(
                surf, removal_rate_1_s, dt)
    gas.TPX = reactor_em.thermo.T, reactor_em.thermo.P, reactor_em.thermo.X
    return carbon_removed


# ═══════════════════════════════════════════════════════════════════════════════
# C. Fluidized bed
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_fluidized_bed(config: ReactorConfig) -> Dict:
    _validate_carbon_policy(config)
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'Fluidized')

    logger.info(
        f"Simulating Fluidized ({config.fluidized_mode}): "
        f"{config.catalyst_name} at {config.T_inlet_K} K"
    )
    gas, graphite, surf = _load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    ch4_initial = _species_x(gas, 'CH4') or 1.0
    ar_initial = _species_x(gas, 'Ar')
    _require_ar_tracer(ch4_initial, ar_initial)

    u0 = max(config.gas_velocity_m_s, 0.05)
    umf = config.u_mf_m_s
    delta = min(0.5, max(0.01, (u0 - umf) / u0))
    tau_emulsion = config.bed_height_m * (1 - delta) / umf

    sv_ratio = active_sv(geometric_sv_fluidized(config), config)
    circulating = config.fluidized_mode == FLUIDIZED_CIRCULATING
    removal_rate = config.circulating_carbon_removal_rate_1_s if circulating else 0.0
    carbon_removed = _integrate_fluidized_pass(
        gas, surf, tau_emulsion, sv_ratio, removal_rate)

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
            carbon_removed += _integrate_fluidized_pass(
                gas, surf, tau_emulsion, sv_ratio, 0.0)
            theta = _coverage(surf, 'C_s')
            regen_cycles += 1
            per_cycle.append(_ch4_extent(gas, ch4_initial, ar_initial))

    final_conv = _ch4_extent(gas, ch4_initial, ar_initial)
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
        'residence_time_s': tau_emulsion,
        'WHSV_h-1': reciprocal_residence_h(tau_emulsion),
        'CH4_conversion': float(final_conv),
        'single_pass_CH4_conversion': float(final_conv),
        'CH4_mole_fraction_drop': _mole_fraction_drop(gas, ch4_initial),
        'conversion_basis': 'argon_tracer',
        'thermal_mode': 'isothermal_energy_disabled',
        'exit_x_H2': float(x_h2),
        'exit_T_K': float(gas.T),
        'exit_theta_C': _coverage(surf, 'C_s'),
        'carbon_removed_coverage_proxy': float(carbon_removed),
        'regen_cycles_completed': regen_cycles,
        'per_cycle_CH4_conversion': per_cycle,
        **kinetics_fields(config),
        **solids_inventory_fields(config, geometric_sv_fluidized(config)),
        **_kinetics_evidence(config),
        **_policy_metadata(config),
        **_load_status_fields(config, graphite, surf),
    }
    logger.info(f"  Fluidized result: conversion={final_conv:.2%} mode={config.fluidized_mode}")
    return result


def _mock_reactor_result(config: ReactorConfig, reactor_type: str) -> Dict:
    logger.warning(f"Cantera not available. Generating mock {reactor_type} results.")
    _validate_carbon_policy(config)
    E_act = config.catalyst_E_act_eV
    k_B_eV = 8.617e-5
    k = 1e13 * np.exp(-E_act / (k_B_eV * config.T_inlet_K))
    if reactor_type == 'MMBCR':
        tau = config.column_height_m / config.gas_velocity_m_s
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


def simulate_reactor(config: ReactorConfig) -> Dict:
    simulators = {
        'MMBCR': simulate_mmbcr,
        'PFR': simulate_pfr,
        'Fluidized': simulate_fluidized_bed,
    }
    sim = simulators.get(config.reactor_type, simulate_mmbcr)
    result = sim(config)
    REACTOR_DIR.mkdir(parents=True, exist_ok=True)
    fname = f"{config.reactor_type}_{config.catalyst_name}_{int(config.T_inlet_K)}K.json"
    save_json(result, fname, subdir='reactor')
    return result


def run_reactor_sweep(catalyst_name: str, mechanism_file: str,
                      temperatures: List[float] = None,
                      reactor_types: List[str] = None,
                      catalyst_E_act_eV: float = 0.8,
                      catalyst_dE_H_eV: float = 0.0,
                      reactor_config_kwargs: Optional[Dict] = None) -> List[Dict]:
    if temperatures is None:
        temperatures = [773.15, 900.0, 1100.0, 1300.0]
    if reactor_types is None:
        reactor_types = ['MMBCR', 'PFR', 'Fluidized']
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
                **extra,
            )
            results.append(simulate_reactor(config))
    return results


if __name__ == '__main__':
    print_banner("REACTOR SIMULATION TEST")
    config = ReactorConfig(
        T_inlet_K=1000.0,
        reactor_type='MMBCR',
        catalyst_name='NiBi_test',
    )
    result = simulate_reactor(config)
    print(json.dumps(result, indent=2))
