"""Thermocatalytic two-phase fluidized-bed reactor model."""

import numpy as np

try:
    import cantera as ct
except ImportError:  # The facade selects the explicit mock path first.
    ct = None

from pipeline.data_models.reactors import ReactorResult
from pipeline.reactors.reactor_core import *  # noqa: F403
from pipeline.reactors.reactor_core import logger


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

def simulate_fluidized_bed(config: ReactorConfig, *, cantera_available: bool,
                             mock_result) -> ReactorResult:
    """Two-phase fluidized bed: reacting emulsion plus bubble bypass.

    The emulsion phase reacts at the minimum-fluidization residence time with
    in-step C_s removal (circulating) or batch regen. A bubble fraction δ of
    the feed bypasses unreacted and is mixed with the emulsion outlet on molar
    flows via the Ar tracer, so ``CH4_conversion = (1−δ)·X_emulsion``.
    Interphase mass transfer is not resolved; this is a screening approximation.

    Args:
        config: Configuration controlling this operation.
        cantera_available: Whether the configured Cantera backend is usable.
        mock_result: Explicit fallback used when Cantera is unavailable.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    _validate_carbon_policy(config)
    if not cantera_available:
        return mock_result(config, 'Fluidized')

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
