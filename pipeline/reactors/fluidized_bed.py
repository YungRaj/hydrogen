"""Thermocatalytic two-phase fluidized-bed reactor model."""

from __future__ import annotations

from collections.abc import Callable
from typing import cast

import numpy as np
import numpy.typing as npt

from pipeline.data_models.reactors import (
    CanteraGasPhase,
    CanteraSurfacePhase,
    H2C2SelectivityRateMetrics,
    ReactorResult,
)
from pipeline.reactors.carbon_accounting import (
    advance_reactor_stage,
    apply_continuous_carbon_removal,
    clamp_gas_to_equilibrium,
    compute_h2_c2_selectivity_and_rates,
    feedstock_extent,
    load_gas_and_surface,
    mix_bubble_bypass,
    mole_fraction_drop,
    require_ar_tracer,
    reset_surface_carbon,
    species_x,
    surface_coverage,
)
from pipeline.reactors.reactor_core import (
    FLUIDIZED_CIRCULATING,
    FLUIDIZED_EMULSION_VOIDAGE,
    FLUIDIZED_REMOVAL_SUBSTEPS,
    REGEN_OXIDATIVE,
    R_J_MOL_K,
    ReactorConfig,
    active_sv,
    common_reactor_result_fields,
    geometric_sv_fluidized,
    kunii_levenspiel_bubble_exchange,
    logger,
    off_site_carbon_metrics,
    solids_inventory_fields,
    validate_carbon_policy,
)


def _integrate_fluidized_pass(
    gas: CanteraGasPhase,
    surf: CanteraSurfacePhase | None,
    tau: float,
    sv_ratio: float,
    removal_rate_1_s: float,
    feed_initial: float = 0.95,
    ar_initial: float = 0.05,
) -> float:
    """Advance emulsion residence with C_s removal *during* integrate (Muradov 2005; Parkinson 2019).

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
        advance_reactor_stage(
            gas, surf, FLUIDIZED_EMULSION_VOIDAGE, sv_ratio, dt
        )
        clamp_gas_to_equilibrium(gas, feed_initial, ar_initial)
        if removal_rate_1_s > 0 and surf is not None:
            carbon_removed += apply_continuous_carbon_removal(
                surf, removal_rate_1_s, dt
            )
    return carbon_removed


# ═══════════════════════════════════════════════════════════════════════════════
# C. Fluidized bed
# ═══════════════════════════════════════════════════════════════════════════════


def simulate_fluidized_bed(
    config: ReactorConfig,
    *,
    cantera_available: bool,
    mock_result: Callable[[ReactorConfig, str], ReactorResult],
) -> ReactorResult:
    """Two-phase fluidized bed: reacting emulsion plus bubble bypass.

    The emulsion phase reacts at the minimum-fluidization residence time with
    in-step C_s removal (circulating) or batch regen. A bubble fraction δ of
    the feed bypasses unreacted and is mixed with the emulsion outlet on molar
    flows via the Ar tracer, so ``feedstock_conversion = (1−δ)·X_emulsion``.
    Interphase mass transfer is not resolved; this is a screening approximation.

    Args:
        config: Configuration controlling this operation.
        cantera_available: Whether the configured Cantera backend is usable.
        mock_result: Explicit fallback used when Cantera is unavailable.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    validate_carbon_policy(config)
    if not cantera_available:
        return mock_result(config, 'Fluidized')

    logger.info(
        f"Simulating Fluidized ({config.fluidized_mode}): "
        f"{config.catalyst_name} at {config.T_inlet_K} K"
    )
    gas, graphite, surf = load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    inlet_x = np.array(gas.X, dtype=float)
    feed_initial = species_x(gas, 'CH4') or 1.0
    ar_initial = species_x(gas, 'Ar')
    require_ar_tracer(feed_initial, ar_initial)

    u0 = max(config.gas_velocity_m_s, 0.05)  # Superficial gas velocity floor (0.05 m/s)
    umf = config.u_mf_m_s
    # Kunii-Levenspiel two-phase theory bubble volume fraction (Kunii & Levenspiel,
    # Fluidization Engineering, 1991): delta ~ (u_0 - u_mf) / u_0, clamped to
    # [0.01, 0.50] for the bubbling fluidization regime prior to slugging/turbulent transition.
    delta = (
        float(config.fluidized_bubble_fraction)
        if config.fluidized_bubble_fraction is not None
        else min(0.5, max(0.01, (u0 - umf) / u0))
    )
    if not 0 <= delta < 1:
        raise ValueError('fluidized-bed bubble fraction must lie in [0, 1)')
    # Interstitial gas residence in the emulsion: the emulsion occupies
    # H·(1−δ) of the bed, its gas fraction is ε_mf and the gas moves through it
    # at superficial u_mf, i.e. interstitial u_mf/ε_mf. This matches the PFR
    # convention (L·ε/u) and the ε_mf gas volume integrated per m³ emulsion.
    tau_emulsion = (
        config.bed_height_m * (1 - delta) * FLUIDIZED_EMULSION_VOIDAGE / umf
    )

    sv_ratio = active_sv(geometric_sv_fluidized(config), config)
    circulating = config.fluidized_mode == FLUIDIZED_CIRCULATING
    removal_rate = config.circulating_carbon_removal_rate_1_s if circulating else 0.0
    if surf is not None:
        reset_surface_carbon(surf)
    pass_start_cov: npt.NDArray[np.float64] | None = (
        np.array(surf.coverages, dtype=np.float64) if surf is not None else None
    )
    carbon_removed = _integrate_fluidized_pass(
        gas, surf, tau_emulsion, sv_ratio, removal_rate, feed_initial, ar_initial
    )
    last_pass_removed = carbon_removed

    regen_cycles = 0
    per_cycle: list[float] = []

    if not circulating:
        theta = surface_coverage(surf, 'C_s')
        per_cycle.append(feedstock_extent(gas, feed_initial, ar_initial))
        while (
            config.max_regen_cycles > 0
            and regen_cycles < config.max_regen_cycles
            and theta >= config.regen_coverage_threshold
        ):
            if config.regen_mechanism == REGEN_OXIDATIVE and not config.co2_permitted:
                raise RuntimeError('oxidative regen blocked (co2_permitted=False)')
            reset_surface_carbon(surf)
            gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
            pass_start_cov = (
                np.array(surf.coverages, dtype=np.float64)
                if surf is not None
                else None
            )
            last_pass_removed = _integrate_fluidized_pass(
                gas, surf, tau_emulsion, sv_ratio, 0.0, feed_initial, ar_initial
            )
            carbon_removed += last_pass_removed
            theta = surface_coverage(surf, 'C_s')
            regen_cycles += 1
            per_cycle.append(feedstock_extent(gas, feed_initial, ar_initial))

    emulsion_conv = feedstock_extent(gas, feed_initial, ar_initial)
    exit_theta_c = surface_coverage(surf, 'C_s')

    # Baker-Helveg-Abild-Pedersen off-site carbon accounting on the emulsion basis
    # (1 m³ emulsion: gas volume ε_mf, area sv_ratio), before the bubble bypass is mixed in.
    # Carbon taken off by the circulating outfeed is Muradov-Parkinson solids product,
    # not Cγ, and is subtracted from the filament balance.
    gamma_mol_m2 = float(config.site_density_mol_cm2) * 1e4
    n_sites_mol = gamma_mol_m2 * sv_ratio
    c_total = config.P_inlet_Pa / (R_J_MOL_K * config.T_inlet_K)
    n_parcel_in_mol = c_total * FLUIDIZED_EMULSION_VOIDAGE
    n_feed_fed_mol = feed_initial * n_parcel_in_mol
    off_site = off_site_carbon_metrics(
        config,
        gas,
        surf,
        feed_initial=feed_initial,
        ar_initial=ar_initial,
        pass_conversion=float(emulsion_conv),
        n_sites_mol=n_sites_mol,
        n_feed_fed_mol=n_feed_fed_mol,
        n_parcel_in_mol=n_parcel_in_mol,
        pass_time_s=tau_emulsion,
        cov_start=pass_start_cov,
        removed_surface_carbon_mol=last_pass_removed * n_sites_mol,
        basis='Gamma*sv_emulsion / (c_CH4*eps_mf), emulsion parcel before bypass',
    )

    # Bubble bypass: δ of the feed passes unreacted; mix on molar flows (Ar tracer).
    mix_bubble_bypass(gas, inlet_x, delta, ar_initial)
    final_conv = feedstock_extent(gas, feed_initial, ar_initial)
    expected = (1.0 - delta) * emulsion_conv
    if abs(final_conv - expected) > 1e-6 + 1e-6 * abs(expected):
        raise RuntimeError(
            f'bypass mixing inconsistent: X_mixed={final_conv:.6f} vs '
            f'(1-delta)*X_em={expected:.6f}'
        )
    x_h2 = species_x(gas, 'H2')
    kl_metrics = kunii_levenspiel_bubble_exchange(config, delta, u0, umf)

    circulating_outfeed_active = bool(circulating and removal_rate > 0.0)
    carbon_outfeed_policy_summary = (
        f"fluidized_mode={config.fluidized_mode};"
        f"circulating_outfeed={circulating_outfeed_active};"
        f"removal_rate_1_s={removal_rate:.3f};"
        f"exit_theta_C={float(exit_theta_c or 0.0):.3f}"
    )

    h2_c2_metrics: H2C2SelectivityRateMetrics = compute_h2_c2_selectivity_and_rates(
        gas,
        feed_initial=feed_initial,
        ar_initial=ar_initial,
        final_conv=float(final_conv),
        tau_s=tau_emulsion,
        P_inlet_Pa=config.P_inlet_Pa,
        T_inlet_K=config.T_inlet_K,
        active_sv_1_m=sv_ratio * (1.0 - delta),
        site_density_mol_cm2=float(config.site_density_mol_cm2),
        bed_porosity_or_holdup=FLUIDIZED_EMULSION_VOIDAGE * (1.0 - delta) + delta,
    )

    single_pass_conv: float = float(
        (1.0 - delta) * per_cycle[0] if per_cycle else final_conv
    )
    cycle_avg_conv: float = float(
        (1.0 - delta) * np.mean(per_cycle) if per_cycle else final_conv
    )
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
            'override'
            if config.fluidized_bubble_fraction is not None
            else 'clip((u0-umf)/u0, 0.01, 0.5)'
        ),
        **kl_metrics,
        'residence_time_s': tau_emulsion,
        'emulsion_voidage': FLUIDIZED_EMULSION_VOIDAGE,
        'surface_area_basis': 'sv_per_emulsion_volume_times_emulsion_volume',
        'emulsion_CH4_conversion': float(emulsion_conv),
        'feedstock_conversion': float(final_conv),
        'single_pass_feedstock_conversion': single_pass_conv,
        'cycle_averaged_feedstock_conversion': cycle_avg_conv,
        'CH4_conversion': float(final_conv),
        'kunii_levenspiel_coupled_CH4_conversion': float(
            (1.0 - kl_metrics['effective_bypass_fraction']) * emulsion_conv
        ),
        'single_pass_CH4_conversion': single_pass_conv,
        'cycle_averaged_CH4_conversion': cycle_avg_conv,
        'circulating_outfeed_active': circulating_outfeed_active,
        'carbon_outfeed_policy_summary': carbon_outfeed_policy_summary,
        'bypass_basis': 'ar_tracer_molar_mix',
        'CH4_mole_fraction_drop': mole_fraction_drop(gas, feed_initial),
        'conversion_basis': 'argon_tracer',
        'thermal_mode': 'isothermal_energy_disabled',
        'exit_x_H2': float(x_h2),
        'exit_T_K': float(gas.T),
        'exit_theta_C': exit_theta_c,
        **h2_c2_metrics,
        'carbon_removed_coverage_proxy': float(carbon_removed),
        'regen_cycles_completed': regen_cycles,
        'per_cycle_CH4_conversion': per_cycle,
        'per_cycle_basis': 'emulsion_only_before_bypass_mix',
        **off_site,
        **solids_inventory_fields(config, geometric_sv_fluidized(config)),
        **common_reactor_result_fields(
            config, final_conv, tau_emulsion, feed_initial, graphite, surf
        ),
    }
    logger.info(
        f"  Fluidized result: conversion={final_conv:.2%} "
        f"(emulsion {emulsion_conv:.2%}, delta={delta:.2f}) mode={config.fluidized_mode}"
    )
    return cast(ReactorResult, result)
