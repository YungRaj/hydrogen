"""Molten-metal bubble-column methane-conversion model."""

from __future__ import annotations

from collections.abc import Callable
from typing import cast

import numpy as np

from pipeline.data_models.reactors import ReactorResult
from pipeline.reactors.reactor_core import (
    G_M_S2,
    R_J_MOL_K,
    ReactorConfig,
    common_reactor_result_fields,
    compute_h2_c2_selectivity_and_rates,
    load_gas_and_surface,
    local_tp_x_eq,
    logger,
    mmbcr_flotation_eta,
    mmbcr_hydrodynamics,
    mmbcr_interfacial_k_m_s,
    set_gas_from_ch4_conversion,
    species_x,
    tabulated_x_eq,
    validate_carbon_policy,
)


def simulate_mmbcr(
    config: ReactorConfig,
    *,
    cantera_available: bool,
    mock_result: Callable[[ReactorConfig, str], ReactorResult],
) -> ReactorResult:
    """Melt ODE to tabulated X_eq with bubble-area flotation (Upham 2017 / Chen 2023).

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
    validate_carbon_policy(config)
    hydro = mmbcr_hydrodynamics(config)
    if not cantera_available:
        return mock_result(config, 'MMBCR')

    logger.info(f"Simulating MMBCR: {config.catalyst_name} at {config.T_inlet_K} K")
    gas, graphite, surf = load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    x_ch4_feed = species_x(gas, 'CH4') or 0.95
    x_ar_feed = species_x(gas, 'Ar')
    if x_ar_feed <= 0:
        x_ar_feed = max(0.0, 1.0 - x_ch4_feed)

    tau_total = hydro['residence_time_s']
    tau_stage = tau_total / config.n_cstr_stages
    sv_ratio = hydro['interfacial_sv_ratio_1_m']  # bubble S/V (m² interface / m³ gas)
    x_eq = tabulated_x_eq(config.T_inlet_K)
    k_if = mmbcr_interfacial_k_m_s(
        config.catalyst_E_act_eV, config.T_inlet_K, config.mmbcr_interfacial_k0_m_s
    )
    eta_float = mmbcr_flotation_eta(
        k_if, sv_ratio, config.mmbcr_carbon_removal_rate_1_s
    )
    # Dimensionless interfacial Damkohler number per stage:
    # Da_stage = k_if [m/s] * (6 / d_b) [m^-1] * tau_stage [s] * eta_float
    da_stage = k_if * sv_ratio * tau_stage * eta_float

    rho_l = float(hydro['melt_density_kg_m3'])
    eps_g = float(hydro['gas_holdup_fraction'])
    # Sparger-level pressure including hydrostatic head of aerated liquid column:
    # P_bottom = P_top + rho_L * g * (1 - eps_g) * H
    inlet_bottom_p_pa = float(
        config.P_inlet_Pa
        + rho_l * G_M_S2 * max(0.0, 1.0 - eps_g) * float(config.column_height_m)
    )

    z_positions = np.linspace(0, config.column_height_m, config.n_cstr_stages + 1)
    conversion_profile = [0.0]
    temperature_profile = [config.T_inlet_K]
    species_profiles = {
        sp: [species_x(gas, sp)] for sp in ['CH4', 'H2', 'C2H2', 'C2H4', 'C2H6']
    }
    conv = 0.0
    carbon_removed_coverage = 0.0

    for _stage in range(config.n_cstr_stages):
        # Local hydrostatic pressure P(z) = P_top + rho_L * g * (1 - eps_g) * (H - z)
        # and pressure-adjusted equilibrium ceiling
        z_mid = float(
            0.5 * (z_positions[_stage] + z_positions[_stage + 1])
        )
        p_stage = float(
            config.P_inlet_Pa
            + rho_l
            * G_M_S2
            * max(0.0, 1.0 - eps_g)
            * max(0.0, float(config.column_height_m) - z_mid)
        )
        x_eq_local = min(
            x_eq,
            local_tp_x_eq(config.T_inlet_K, p_stage, x_ch4_feed),
        )
        # Account for gas molar expansion (1 + y_CH4,0 * X) from CH4 -> 2 H2 and
        # hydrostatic decompression (P_bottom / P(z), guarded by a 10 kPa = 0.1 bar
        # floor), which accelerates bubbles and reduces local stage residence time.
        expansion_factor = (1.0 + x_ch4_feed * max(0.0, conv)) * (
            inlet_bottom_p_pa / max(p_stage, 10000.0)
        )
        da_stage_local = da_stage / max(expansion_factor, 1.0)
        # First-order approach to local melt/gas equilibrium; C floats out of the bubble.
        if conv < x_eq_local:
            conv = x_eq_local - (x_eq_local - conv) * float(np.exp(-da_stage_local))
        conv = float(min(conv, x_eq))
        carbon_removed_coverage += max(0.0, conv - conversion_profile[-1]) * x_ch4_feed
        set_gas_from_ch4_conversion(
            gas, config.T_inlet_K, p_stage, x_ch4_feed, x_ar_feed, conv
        )
        conversion_profile.append(float(conv))
        temperature_profile.append(config.T_inlet_K)
        for sp in species_profiles:
            species_profiles[sp].append(species_x(gas, sp))

    # Top of column is at P_inlet_Pa
    set_gas_from_ch4_conversion(
        gas, config.T_inlet_K, config.P_inlet_Pa, x_ch4_feed, x_ar_feed, conv
    )
    final_conv = conversion_profile[-1]

    # Competitive homogeneous bubble-interior C2 pyrolysis (2 CH4 -> C2H6 + H2 -> C2H4 -> C2H2,
    # Ea = 356 kJ/mol, A_1st = 4.5e13 1/s, Dean J. Phys. Chem. 1990 / Palmer et al. ACS Catal. 2020)
    # vs. heterogeneous melt-bubble interfacial dehydrogenation (k_if * sv_ratio):
    k_hom_1_s = float(
        4.5e13
        * np.exp(-356000.0 / (R_J_MOL_K * max(float(config.T_inlet_K), 300.0)))
    )
    k_het_1_s = float(k_if * sv_ratio * eta_float)
    if final_conv > 1e-6 and (k_hom_1_s + k_het_1_s) > 0.0:
        # C2 carbon branching fraction f_C2 bounded by thermodynamic C2 ceiling (~8% at 1300 K)
        f_c2_branch = float(
            np.clip(0.08 * k_hom_1_s / (k_het_1_s + k_hom_1_s), 0.0, 0.08)
        )
    else:
        f_c2_branch = 0.0

    # Partition C2 carbon among C2H4 (60%), C2H6 (25%), and C2H2 (15%) at bubble exit
    dn_ch4 = x_ch4_feed * final_conv
    n_c2_total = 0.5 * dn_ch4 * f_c2_branch
    n_c2h6 = 0.25 * n_c2_total
    n_c2h4 = 0.60 * n_c2_total
    n_c2h2 = 0.15 * n_c2_total
    n_ch4_out = x_ch4_feed * (1.0 - final_conv)
    n_h2_out = max(
        0.0,
        2.0 * dn_ch4 - (3.0 * n_c2h6 + 2.0 * n_c2h4 + 1.0 * n_c2h2),
    )
    n_ar_out = x_ar_feed
    n_tot_out = max(
        n_ch4_out + n_h2_out + n_ar_out + n_c2h6 + n_c2h4 + n_c2h2, 1e-12
    )
    if all(
        sp in gas.species_names for sp in ('CH4', 'H2', 'Ar', 'C2H6', 'C2H4', 'C2H2')
    ):
        gas.TPX = (
            config.T_inlet_K,
            config.P_inlet_Pa,
            f"CH4:{n_ch4_out / n_tot_out}, H2:{n_h2_out / n_tot_out}, "
            f"C2H6:{n_c2h6 / n_tot_out}, C2H4:{n_c2h4 / n_tot_out}, "
            f"C2H2:{n_c2h2 / n_tot_out}, Ar:{n_ar_out / n_tot_out}",
        )
    x_h2 = species_x(gas, 'H2')

    h2_c2_metrics = compute_h2_c2_selectivity_and_rates(
        gas,
        ch4_initial=x_ch4_feed,
        ar_initial=x_ar_feed,
        final_conv=float(final_conv),
        tau_s=tau_total,
        P_inlet_Pa=config.P_inlet_Pa,
        T_inlet_K=config.T_inlet_K,
        active_sv_1_m=sv_ratio * eps_g,
        site_density_mol_cm2=float(config.site_density_mol_cm2),
        bed_porosity_or_holdup=eps_g,
    )

    result = {
        'reactor_type': 'MMBCR',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'P_Pa': config.P_inlet_Pa,
        'inlet_bottom_P_Pa': inlet_bottom_p_pa,
        'column_height_m': config.column_height_m,
        'column_diameter_m': config.column_diameter_m,
        'gas_velocity_m_s': config.gas_velocity_m_s,
        'bubble_diameter_mm': config.bubble_diameter_mm,
        **hydro,
        'CH4_conversion': float(final_conv),
        **h2_c2_metrics,
        'solid_C_selectivity_note': (
            'competitive interfacial melt dehydrogenation vs bubble-interior C2 pyrolysis'
        ),
        'c2_tracked': True,
        'exit_x_H2': x_h2,
        'exit_T_K': float(gas.T),
        'exit_theta_C': 0.0,
        'carbon_removed_coverage_proxy': float(carbon_removed_coverage),
        'X_eq_table': float(x_eq),
        'X_eq_local_tp': float(
            local_tp_x_eq(config.T_inlet_K, config.P_inlet_Pa, x_ch4_feed)
        ),
        'mmbcr_k_if_m_s': float(k_if),
        'mmbcr_flotation_eta': float(eta_float),
        'mmbcr_Da': float(k_if * sv_ratio * tau_total * eta_float),
        'conversion_basis': 'melt_ode_to_Xeq',
        'thermal_mode': 'isothermal_by_construction',
        'z_positions': z_positions.tolist(),
        'conversion_profile': conversion_profile,
        'temperature_profile': temperature_profile,
        'theta_C_profile': [0.0] * len(conversion_profile),
        **common_reactor_result_fields(
            config, final_conv, tau_total, x_ch4_feed, graphite, surf
        ),
    }
    logger.info(
        f"  MMBCR result: conversion={final_conv:.2%}, "
        f"H2_selectivity={result['H2_selectivity']:.2%}, "
        f"H2_yield={result['H2_yield']:.2%}, "
        f"τ={tau_total:.1f}s (u_b={hydro['bubble_rise_velocity_m_s']:.3f} m/s, "
        f"eps_g={hydro['gas_holdup_fraction']:.3f})"
    )
    return cast(ReactorResult, result)
