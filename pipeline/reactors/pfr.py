"""Thermocatalytic packed-bed plug-flow reactor model."""

from __future__ import annotations

from collections.abc import Callable
from typing import cast

import numpy as np
import numpy.typing as npt

from pipeline.data_models.reactors import (
    CanteraGasPhase,
    CanteraSurfacePhase,
    ReactorResult,
)
from pipeline.reactors.reactor_core import (
    PFR_PILOT_TAU_BAND_S,
    REGEN_OXIDATIVE,
    R_J_MOL_K,
    ReactorConfig,
    active_sv,
    advance_reactor_stage,
    ch4_extent,
    clamp_gas_to_equilibrium,
    common_reactor_result_fields,
    compute_h2_c2_selectivity_and_rates,
    geometric_sv_pfr,
    load_gas_and_surface,
    logger,
    mole_fraction_drop,
    off_site_carbon_metrics,
    require_ar_tracer,
    reset_surface_carbon,
    solids_inventory_fields,
    species_x,
    surface_coverage,
    validate_carbon_policy,
)


def _run_pfr_pass(
    gas: CanteraGasPhase,
    surf: CanteraSurfacePhase | None,
    *,
    n_stages: int,
    t_inlet_k: float,
    ch4_initial: float,
    ar_initial: float,
    stage_gas_volume: float,
    stage_surface_area: float,
    tau_stage: float,
) -> tuple[list[float], list[float], list[float], npt.NDArray[np.float64] | None, float]:
    """Advance one full Lagrangian PFR pass across all axial stages without nonlocal state."""
    conversion_profile: list[float] = [0.0]
    temperature_profile: list[float] = [float(t_inlet_k)]
    theta_c_tos: list[float] = [surface_coverage(surf, 'C_s')]
    pass_start_cov: npt.NDArray[np.float64] | None = (
        np.array(surf.coverages, dtype=np.float64)
        if surf is not None
        else None
    )
    pass_produce_time_s = 0.0
    for _ in range(n_stages):
        advance_reactor_stage(
            gas, surf, stage_gas_volume, stage_surface_area, tau_stage
        )
        clamp_gas_to_equilibrium(gas, ch4_initial, ar_initial)
        pass_produce_time_s += tau_stage
        conversion_profile.append(ch4_extent(gas, ch4_initial, ar_initial))
        temperature_profile.append(float(gas.T))
        theta_c_tos.append(surface_coverage(surf, 'C_s'))
    return (
        conversion_profile,
        temperature_profile,
        theta_c_tos,
        pass_start_cov,
        pass_produce_time_s,
    )


def simulate_pfr(
    config: ReactorConfig,
    *,
    cantera_available: bool,
    mock_result: Callable[[ReactorConfig, str], ReactorResult],
) -> ReactorResult:
    """Packed bed as a staged Lagrangian PFR with one shared surface.

    Surface area per stage is ``sv·V_bed_stage`` where ``sv = 6(1−ε)/d_p``
    is area per *bed* volume (not ε-scaled), times loading × dispersion.
    The gas parcel volume is ``ε·V_bed_stage``. Stages are isothermal at
    ``T_inlet``; CH4 conversion is the Ar-tracer extent.

    Args:
        config: Configuration controlling this operation.
        cantera_available: Whether the configured Cantera backend is usable.
        mock_result: Explicit fallback used when Cantera is unavailable.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    validate_carbon_policy(config)
    if not cantera_available:
        return mock_result(config, 'PFR')

    logger.info(f"Simulating PFR: {config.catalyst_name} at {config.T_inlet_K} K")
    gas, graphite, surf = load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    eps = config.bed_porosity
    sv_ratio = active_sv(geometric_sv_pfr(config), config)  # m² active / m³ bed
    ch4_initial = species_x(gas, 'CH4') or 1.0
    ar_initial = species_x(gas, 'Ar')
    require_ar_tracer(ch4_initial, ar_initial)

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
    cycles_completed = 0
    per_cycle_conversion: list[float] = []
    produce_time_s = 0.0

    # One produce pass (always). Optional discrete regen cycles if configured.
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    if surf is not None:
        reset_surface_carbon(surf)
    (
        conversion_profile,
        temperature_profile,
        theta_C_tos,
        pass_start_cov,
        pass_dt,
    ) = _run_pfr_pass(
        gas,
        surf,
        n_stages=n_stages,
        t_inlet_k=config.T_inlet_K,
        ch4_initial=ch4_initial,
        ar_initial=ar_initial,
        stage_gas_volume=stage_gas_volume,
        stage_surface_area=stage_surface_area,
        tau_stage=tau_stage,
    )
    produce_time_s += pass_dt
    per_cycle_conversion.append(conversion_profile[-1])

    while (
        config.max_regen_cycles > 0
        and cycles_completed < config.max_regen_cycles
        and theta_C_tos
        and max(float(v or 0.0) for v in theta_C_tos) >= config.regen_coverage_threshold
    ):
        if config.regen_mechanism == REGEN_OXIDATIVE and not config.co2_permitted:
            raise RuntimeError('oxidative regen blocked (co2_permitted=False)')
        if config.regen_mechanism == REGEN_OXIDATIVE:
            logger.warning(
                'Oxidative PFR regen enabled via co2_permitted=True (test only)'
            )
        # Mechanical / consumable: free-site reset without CO2 chemistry in-model.
        reset_surface_carbon(surf)
        gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
        cycles_completed += 1
        (
            conversion_profile,
            temperature_profile,
            theta_C_tos,
            pass_start_cov,
            pass_dt,
        ) = _run_pfr_pass(
            gas,
            surf,
            n_stages=n_stages,
            t_inlet_k=config.T_inlet_K,
            ch4_initial=ch4_initial,
            ar_initial=ar_initial,
            stage_gas_volume=stage_gas_volume,
            stage_surface_area=stage_surface_area,
            tau_stage=tau_stage,
        )
        produce_time_s += pass_dt
        per_cycle_conversion.append(conversion_profile[-1])

    final_conv = conversion_profile[-1]
    x_h2 = species_x(gas, 'H2')

    # Carbon accounting on the parcel basis (last pass).
    gamma_mol_m2 = float(config.site_density_mol_cm2) * 1e4
    n_sites_mol = gamma_mol_m2 * stage_surface_area
    c_total = config.P_inlet_Pa / (R_J_MOL_K * config.T_inlet_K)  # mol/m³
    n_parcel_in_mol = c_total * stage_gas_volume
    n_ch4_fed_mol = ch4_initial * n_parcel_in_mol
    off_site = off_site_carbon_metrics(
        config,
        gas,
        surf,
        ch4_initial=ch4_initial,
        ar_initial=ar_initial,
        pass_conversion=float(final_conv),
        n_sites_mol=n_sites_mol,
        n_ch4_fed_mol=n_ch4_fed_mol,
        n_parcel_in_mol=n_parcel_in_mol,
        pass_time_s=tau_total,
        cov_start=pass_start_cov,
    )

    inv_fields = solids_inventory_fields(config, geometric_sv_pfr(config))
    tau_in_band = bool(
        PFR_PILOT_TAU_BAND_S[0] <= tau_total <= PFR_PILOT_TAU_BAND_S[1]
    )
    ergun_ok = bool(inv_fields.get('ergun_ok', True))
    max_theta_c = float(max(theta_C_tos) if theta_C_tos else 0.0)
    overloaded_bed_warning = bool(not tau_in_band or not ergun_ok)
    plugging_risk = bool(
        max_theta_c >= config.regen_coverage_threshold or not ergun_ok
    )
    carbon_outfeed_policy_summary = (
        f"regen_mechanism={config.regen_mechanism};"
        f"regen_cycles={cycles_completed};"
        f"max_theta_C={max_theta_c:.3f};"
        f"plugging_risk={plugging_risk}"
    )

    h2_c2_metrics = compute_h2_c2_selectivity_and_rates(
        gas,
        ch4_initial=ch4_initial,
        ar_initial=ar_initial,
        final_conv=float(final_conv),
        tau_s=tau_total,
        P_inlet_Pa=config.P_inlet_Pa,
        T_inlet_K=config.T_inlet_K,
        active_sv_1_m=sv_ratio,
        site_density_mol_cm2=float(config.site_density_mol_cm2),
        bed_porosity_or_holdup=eps,
    )

    result = {
        'reactor_type': 'PFR',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'bed_length_m': config.bed_length_m,
        'bed_diameter_m': config.bed_diameter_m,
        'catalyst_particle_mm': config.catalyst_particle_mm,
        'residence_time_s': tau_total,
        'residence_time_within_pilot_band': tau_in_band,
        'pilot_residence_time_band_s': list(PFR_PILOT_TAU_BAND_S),
        'overloaded_bed_warning': overloaded_bed_warning,
        'plugging_risk': plugging_risk,
        'carbon_outfeed_policy_summary': carbon_outfeed_policy_summary,
        'produce_time_s': produce_time_s,
        'CH4_conversion': float(final_conv),
        'single_pass_CH4_conversion': float(
            per_cycle_conversion[0] if per_cycle_conversion else final_conv
        ),
        'cycle_averaged_CH4_conversion': float(
            np.mean(per_cycle_conversion) if per_cycle_conversion else final_conv
        ),
        'CH4_mole_fraction_drop': mole_fraction_drop(gas, ch4_initial),
        'conversion_basis': 'argon_tracer',
        'thermal_mode': 'isothermal_energy_disabled',
        'surface_area_basis': 'sv_per_bed_volume_times_bed_volume',
        'stage_surface_area_m2': float(stage_surface_area),
        'stage_gas_volume_m3': float(stage_gas_volume),
        'per_cycle_CH4_conversion': per_cycle_conversion,
        'regen_cycles_completed': cycles_completed,
        'exit_x_H2': x_h2,
        'exit_T_K': float(gas.T),
        **h2_c2_metrics,
        'z_positions': z_positions.tolist(),
        'conversion_profile': conversion_profile,
        'temperature_profile': temperature_profile,
        'theta_C_time_on_stream': theta_C_tos,
        'theta_C_profile_basis': 'single_shared_surface_cumulative_time_on_stream',
        'theta_C_axial': theta_C_tos,
        'inlet_theta_C': float(theta_C_tos[1] if len(theta_C_tos) > 1 else 0.0),
        'max_theta_C': max_theta_c,
        'exit_theta_C': float(theta_C_tos[-1]) if theta_C_tos else None,
        **off_site,
        **inv_fields,
        **common_reactor_result_fields(
            config, final_conv, tau_total, ch4_initial, graphite, surf
        ),
    }
    logger.info(
        f"  PFR result: conversion={final_conv:.2%}, τ={tau_total:.1f}s, "
        f"regen_cycles={cycles_completed}"
    )
    return cast(ReactorResult, result)
