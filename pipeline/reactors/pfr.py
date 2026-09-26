"""Thermocatalytic packed-bed plug-flow reactor model."""

from typing import List, Optional

import numpy as np

try:
    import cantera as ct
except ImportError:  # The facade selects the explicit mock path first.
    ct = None

from pipeline.data_models.reactors import ReactorResult
from pipeline.reactors.reactor_core import *  # noqa: F403
from pipeline.reactors.reactor_core import logger


def simulate_pfr(config: ReactorConfig, *, cantera_available: bool,
             mock_result) -> ReactorResult:
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
    _validate_carbon_policy(config)
    if not cantera_available:
        return mock_result(config, 'PFR')

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
