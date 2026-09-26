"""Molten-metal bubble-column methane-conversion model."""

import numpy as np

from pipeline.data_models.reactors import ReactorResult
from pipeline.reactors.reactor_core import *  # noqa: F403
from pipeline.reactors.reactor_core import logger


def simulate_mmbcr(config: ReactorConfig, *, cantera_available: bool,
             mock_result) -> ReactorResult:
    """Melt ODE to tabulated X_eq with bubble-area flotation (B3).

    dX/dz-style first-order approach: per stage
    ``X ← X_eq − (X_eq − X)·exp(−Da_stage)`` with
    ``Da = k_if(E_act, T)·(6/d_b)·τ·η``. τ = H / u_b (Mendelson). Column
    height is the design lever; gas holdup is derived. Carbon floats out of
    the bubble and is reconstructed from the CH4/H2/Ar balance.

    Args:
        config: Configuration controlling this operation.
        cantera_available: Whether the configured Cantera backend is usable.
        mock_result: Explicit fallback used when Cantera is unavailable.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    _validate_carbon_policy(config)
    hydro = _mmbcr_hydrodynamics(config)
    if not cantera_available:
        return mock_result(config, 'MMBCR')

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
