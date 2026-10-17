"""B5/B6 off-site carbon accounting and Cantera surface state helpers."""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

try:
    import cantera as ct
except ImportError:
    ct = None

METAL_MOLAR_MASS_G_MOL = {'Ni': 58.6934, 'Fe': 55.845, 'Co': 58.9332}
# Ermakova 2000 / Takenaka: 40-384 gC/gNi over 4-50 h on high-Ni/SiO2 -> ~8-10 gC/(gNi h).
NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H = (8.0, 10.0)
# Same TOS literature: the catalyst dies (encapsulated) after 4-50 h.
NI_TOS_LIFETIME_BAND_H = (4.0, 50.0)
# B6-6 encapsulation-onset markers.
ENCAP_ONSET_THETA = 0.1
ENCAP_ONSET_RATIO = 10.0
ENCAP_DEAD_THETA = 0.5
PFR_PILOT_TAU_BAND_S = (2.0, 12.0)
CH4_PYROLYSIS_DELTA_H_KJ_MOL = 74.85


def _load_gas_and_surface(config):
    gas = ct.Solution(config.mechanism_file, 'gas')
    if 'C_graphite' in gas.species_names:
        raise RuntimeError(
            'mechanism still contains gas-phase C_graphite; regenerate YAML'
        )
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


def _load_status_fields(config, graphite, surf) -> dict:
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
    """Transport-style continuous removal of C_s → free sites.

    Returns approximate coverage of C removed (not moles).
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
    site_mask = np.array([surf.species(n).size > 0 for n in range(surf.n_species)])
    s = cov[site_mask].sum()
    if s > 0:
        cov[site_mask] /= s
    surf.coverages = cov
    return float(removed)


def _reset_surface_carbon(surf) -> None:
    """Mechanical / consumable regen: clear adsorbates and restore free sites.

    Encapsulating carbon (``C_encap_s``, B6 Cδ) is time-on-stream death and
    is never cleared.
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


def _surface_carbon_coverage(surf, cov: np.ndarray) -> float:
    """Total θ of carbon-bearing site-occupying species (C_s, CH_x_s, C_encap_s)."""
    total = 0.0
    for k in range(surf.n_species):
        sp = surf.species(k)
        if sp.size > 0 and sp.composition.get('C', 0) > 0:
            total += float(cov[k])
    return total


def _tabulated_x_eq(T_K: float) -> float:
    from pipeline.reactors.equilibrium import interpolate_x_eq

    return interpolate_x_eq(T_K)


def _off_site_carbon_metrics(
    config,
    gas,
    surf,
    *,
    ch4_initial: float,
    ar_initial: float,
    pass_conversion: float,
    n_sites_mol: float,
    n_ch4_fed_mol: float,
    n_parcel_in_mol: float,
    pass_time_s: float,
    cov_start: Optional[np.ndarray],
    removed_surface_carbon_mol: float = 0.0,
    basis: str = 'Gamma*A_stage / (c_CH4*eps*V_stage), PFR parcel',
    mechanism_metadata_fn=None,
) -> Dict:
    """Per-pass carbon accounting for the B5/B6 closure criterion."""
    x_bound = n_sites_mol / n_ch4_fed_mol if n_ch4_fed_mol > 0 else None
    x_eq = _tabulated_x_eq(config.T_inlet_K)
    out = {
        'site_inventory_bound_X': x_bound,
        'site_inventory_bound_basis': basis,
        'X_eq_table': x_eq,
        'exceeds_equilibrium': bool(pass_conversion > x_eq + 1e-6),
    }
    if surf is None:
        return out
    cov_end = np.array(surf.coverages, dtype=float)
    if cov_start is None:
        cov_start = np.zeros_like(cov_end)
    x_ar = _species_x(gas, 'Ar') or ar_initial
    n_ar = ar_initial * n_parcel_in_mol
    n_c2 = sum(
        2.0 * (_species_x(gas, sp) or 0.0) / max(x_ar, 1e-30) * n_ar
        for sp in ('C2H2', 'C2H4', 'C2H6')
    )
    n_solid = max(0.0, pass_conversion * n_ch4_fed_mol - n_c2)
    d_surface_c = (
        _surface_carbon_coverage(surf, cov_end)
        - _surface_carbon_coverage(surf, cov_start)
    ) * n_sites_mol
    n_removed_raw = max(0.0, float(removed_surface_carbon_mol))
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
        float(cov_end[surf.species_index('C_encap_s')]) if has_encap else None
    )
    lifetime_h = None
    if has_encap and n_sites_mol > 0 and pass_time_s > 0:
        d_theta_dt = n_delta / n_sites_mol / pass_time_s
        lifetime_h = (
            ENCAP_DEAD_THETA / d_theta_dt / 3600.0 if d_theta_dt > 0 else None
        )
    out.update(
        {
            'carbon_turnovers_per_site': turnovers,
            'turnover_factor_vs_bound': (
                pass_conversion / x_bound if x_bound else None
            ),
            'solid_carbon_mol_per_pass': n_solid,
            'outfeed_carbon_mol_per_pass': n_removed,
            'outfeed_carbon_mol_raw': n_removed_raw,
            'carbon_balance_residual_mol': balance_residual,
            'carbon_balance_ok': bool(
                balance_residual >= -1e-3 * max(n_sites_mol, n_solid, 1e-12)
            ),
            'c_gamma_mol_per_pass': n_gamma if has_encap else None,
            'c_delta_mol_per_pass': n_delta if has_encap else None,
            'c_gamma_to_c_delta_ratio': ratio,
            'exit_theta_C_encap': theta_encap,
            'off_site_carbon_active': has_encap,
        }
    )
    if has_encap and theta_encap is not None:
        out.update(
            {
                'encapsulation_onset': bool(theta_encap >= ENCAP_ONSET_THETA),
                'encapsulation_onset_basis': f'exit theta_encap >= {ENCAP_ONSET_THETA}',
                'c_delta_competitive': (
                    None if ratio is None else bool(ratio < ENCAP_ONSET_RATIO)
                ),
                'c_delta_competitive_basis': f'C_gamma/C_delta < {ENCAP_ONSET_RATIO}',
                'encapsulation_lifetime_h': lifetime_h,
                'encapsulation_lifetime_basis': (
                    f'theta_encap -> {ENCAP_DEAD_THETA} at the pass-averaged '
                    'C_delta rate; linear extrapolation'
                ),
                'tos_lifetime_literature_band_h': list(NI_TOS_LIFETIME_BAND_H),
                'lifetime_within_tos_band': (
                    None
                    if lifetime_h is None
                    else bool(
                        NI_TOS_LIFETIME_BAND_H[0]
                        <= lifetime_h
                        <= NI_TOS_LIFETIME_BAND_H[1]
                    )
                ),
            }
        )
    if has_encap:
        metadata = (
            mechanism_metadata_fn(config) if mechanism_metadata_fn is not None else {}
        )
        genome = metadata.get('inputs', {}).get('genome')
        metal = None
        try:
            from pipeline.reactors.mechanisms import (
                OFF_SITE_CARBON_METALS,
                _nanoparticle_metals,
                parse_catalyst_genome,
            )

            metals = (
                _nanoparticle_metals(
                    parse_catalyst_genome(genome),
                    metadata.get('inputs', {}).get('material_class'),
                )
                & OFF_SITE_CARBON_METALS
            )
            if len(metals) == 1:
                metal = next(iter(metals))
        except Exception:
            metal = None
        dispersion = float(config.metal_dispersion)
        yield_rate = None
        if metal and dispersion > 0 and pass_time_s > 0 and n_sites_mol > 0:
            n_metal = n_sites_mol / dispersion
            yield_rate = (
                (n_gamma / pass_time_s)
                * 3600.0
                * 12.011
                / (n_metal * METAL_MOLAR_MASS_G_MOL[metal])
            )
        out.update(
            {
                'filament_metal': metal,
                'filament_yield_gC_per_gMetal_h': yield_rate,
                'filament_yield_basis': (
                    'pass-averaged C_gamma rate; metal moles = sites / dispersion'
                ),
                'filament_yield_literature_band_gC_per_gNi_h': list(
                    NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H
                ),
                'filament_yield_within_band': (
                    None
                    if yield_rate is None
                    else bool(
                        NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H[0]
                        <= yield_rate
                        <= NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H[1]
                    )
                ),
            }
        )
    return out


def _h2_atom_balance_metric(
    ch4_initial: float, final_conv: float, x_h2: float
) -> float:
    """H-atom balance metric — not true branching H2 selectivity."""
    if final_conv <= 0.01:
        return 0.0
    h_in_ch4 = 4.0 * ch4_initial
    h_in_h2 = 2.0 * x_h2
    return float(np.clip(h_in_h2 / max(h_in_ch4 * final_conv, 1e-10), 0, 1))


def _ch4_extent(gas, x_ch4_feed: float, x_ar_feed: float) -> float:
    from pipeline.reactors.equilibrium import ch4_conversion_from_argon_tracer

    return ch4_conversion_from_argon_tracer(
        _species_x(gas, 'CH4'), _species_x(gas, 'Ar'), x_ch4_feed, x_ar_feed
    )


def _require_ar_tracer(x_ch4_feed: float, x_ar_feed: float) -> None:
    if x_ar_feed <= 0 or x_ch4_feed <= 0:
        raise ValueError(
            'solids CH4 conversion uses an Ar mole-fraction tracer; '
            'inlet_composition must include Ar (default CH4:0.95, Ar:0.05)'
        )


def _disable_reactor_energy(reactor) -> str:
    """Hold solids stages at inlet T so they match isothermal MMBCR."""
    if hasattr(reactor, 'energy_enabled'):
        reactor.energy_enabled = False
        return 'isothermal_energy_disabled'
    raise RuntimeError(
        'Cantera reactor has no energy_enabled; cannot force isothermal solids'
    )


def _advance_reactor_stage(
    gas, surf, volume_m3: float, surface_area_m2: float, dt: float
) -> None:
    """Advance an isothermal gas+surface stage, falling back to substepped tight tolerances on stiff CVODES failure."""
    gas_t, gas_p = float(gas.T), float(gas.P)
    gas_x = np.array(gas.X, dtype=float)
    surf_cov = np.array(surf.coverages, dtype=float) if surf is not None else None

    def _step_stage(step_dt: float, *, tight: bool = False) -> None:
        try:
            reactor = ct.IdealGasReactor(gas, clone=False)
        except TypeError:
            reactor = ct.IdealGasReactor(gas)
        reactor.volume = volume_m3
        _disable_reactor_energy(reactor)
        if surf is not None:
            try:
                ct.ReactorSurface(surf, reactor, A=surface_area_m2, clone=False)
            except TypeError:
                ct.ReactorSurface(surf, reactor, A=surface_area_m2)
        net = ct.ReactorNet([reactor])
        if tight:
            net.rtol = 1e-11
            net.atol = 1e-18
            net.max_time_step = max(step_dt / 10.0, 1e-7)
        net.advance(step_dt)
        phase = getattr(reactor, 'phase', None) or reactor.thermo
        gas.TPX = phase.T, phase.P, phase.X

    try:
        _step_stage(dt, tight=False)
        return
    except Exception:
        pass

    gas.TPX = gas_t, gas_p, gas_x
    if surf is not None and surf_cov is not None:
        surf.coverages = surf_cov
    n_sub = 25
    sub_dt = dt / n_sub
    for _ in range(n_sub):
        _step_stage(sub_dt, tight=True)


def _mole_fraction_drop(gas, x_ch4_feed: float) -> float:
    if x_ch4_feed <= 0:
        return 0.0
    return float(max(0.0, 1.0 - _species_x(gas, 'CH4') / x_ch4_feed))


def _set_gas_from_ch4_conversion(
    gas, T_K: float, P_Pa: float, x_ch4_feed: float, x_ar_feed: float, X: float
):
    """CH4 → C(s) + 2 H2; C leaves the bubble by flotation (not in the gas)."""
    X = float(np.clip(X, 0.0, 1.0))
    n_ch4 = x_ch4_feed * (1.0 - X)
    n_h2 = 2.0 * x_ch4_feed * X
    n_ar = x_ar_feed
    n_tot = n_ch4 + n_h2 + n_ar
    if n_tot <= 0:
        return
    gas.TPX = (
        T_K,
        P_Pa,
        f'CH4:{n_ch4 / n_tot}, H2:{n_h2 / n_tot}, Ar:{n_ar / n_tot}',
    )


def _mix_bubble_bypass(
    gas, inlet_x: np.ndarray, delta: float, x_ar_feed: float
) -> None:
    """Mix a bypass stream (fraction δ of feed, unreacted) with the emulsion outlet on a molar basis."""
    x_em = np.array(gas.X, dtype=float)
    x_ar_em = _species_x(gas, 'Ar')
    if x_ar_em <= 0:
        raise RuntimeError(
            'Ar tracer vanished from emulsion outlet; cannot mix bypass'
        )
    expansion = x_ar_feed / x_ar_em
    moles = (
        delta * np.asarray(inlet_x, dtype=float)
        + (1.0 - delta) * expansion * x_em
    )
    total = float(moles.sum())
    if total <= 0:
        raise RuntimeError('bubble bypass mixing produced no moles')
    gas.TPX = gas.T, gas.P, moles / total
