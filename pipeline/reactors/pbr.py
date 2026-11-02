"""Thermocatalytic packed-bed reactor (PBR) model with differential Ergun pressure drop and Thiele pellet diffusion."""

from __future__ import annotations

from typing import Callable, cast

import numpy as np
import numpy.typing as npt

from pipeline.data_models.reactors import (
    PBRPelletEffectivenessMetrics,
    ReactorResult,
)
from pipeline.reactors.reactor_core import (
    CH4_PYROLYSIS_DELTA_H_KJ_MOL,
    ERGUN_DP_LIMIT_PA,
    PFR_PILOT_TAU_BAND_S,
    REGEN_OXIDATIVE,
    R_J_MOL_K,
    ReactorConfig,
    active_area_multiplier,
    active_sv,
    advance_reactor_stage,
    ch4_extent,
    ch4_feed_density_kg_m3,
    ch4_viscosity_pa_s,
    geometric_sv_pfr,
    kinetics_evidence,
    kinetics_fields,
    load_gas_and_surface,
    load_status_fields,
    logger,
    mole_fraction_drop,
    off_site_carbon_metrics,
    policy_metadata,
    reciprocal_residence_h,
    require_ar_tracer,
    reset_surface_carbon,
    solids_inventory_fields,
    species_x,
    surface_coverage,
    validate_carbon_policy,
)

# Solid carbon skeletal density [kg/m³] for dynamic packed-bed voidage reduction.
SOLID_CARBON_DENSITY_KG_M3: float = 2000.0
MIN_PACKED_BED_POROSITY: float = 0.10
MIN_STAGE_PRESSURE_PA: float = 10000.0
CH4_MOLAR_MASS_KG_MOL: float = 0.016043
DEFAULT_STICKING_COEFF: float = 0.05


def pbr_pellet_effectiveness(
    config: ReactorConfig,
    T_K: float,
    P_Pa: float,
    u_sup_m_s: float,
    eps_bed: float,
) -> PBRPelletEffectivenessMetrics:
    """Compute spherical-pellet Thiele modulus, internal η, and Wakao-Funazkri overall Ω.

    Uses Bosanquet combination of molecular (Fuller) and Knudsen pore
    diffusivities inside the catalyst pellet and the Wakao-Funazkri (1979)
    packed-bed correlation ``Sh = 2.0 + 1.1 Re_p^0.6 Sc^(1/3)`` for external
    film mass transfer.

    Args:
        config: Reactor configuration supplying pellet and kinetic parameters.
        T_K: Local stage temperature in kelvin.
        P_Pa: Local stage pressure in pascals.
        u_sup_m_s: Local superficial gas velocity in m/s.
        eps_bed: Local packed-bed void fraction in (0, 1).

    Returns:
        Typed record of dimensionless and unit-bearing mass-transfer metrics.
    """
    _ = eps_bed
    d_p = max(float(config.catalyst_particle_mm) * 1e-3, 1e-6)
    eps_p = float(config.pellet_porosity)
    tau_p = float(config.pellet_tortuosity)
    r_pore_m = float(config.pellet_pore_radius_nm) * 1e-9

    # Fuller-type binary molecular diffusivity CH4 in H2/Ar mixture [m²/s]
    d_mol = (
        2.2e-5
        * ((T_K / 298.15) ** 1.75)
        * (101325.0 / max(P_Pa, MIN_STAGE_PRESSURE_PA))
    )
    # Knudsen diffusivity in mesopores [m²/s]
    v_thermal = float(
        np.sqrt(8.0 * R_J_MOL_K * T_K / (np.pi * CH4_MOLAR_MASS_KG_MOL))
    )
    d_knudsen = (2.0 / 3.0) * r_pore_m * v_thermal
    # Bosanquet effective intraparticle diffusivity [m²/s]
    d_bosanquet = 1.0 / (1.0 / max(d_mol, 1e-12) + 1.0 / max(d_knudsen, 1e-12))
    d_eff = max((eps_p / max(tau_p, 1e-6)) * d_bosanquet, 1e-12)

    # Intrinsic first-order surface reaction velocity k_s'' [m/s] from kinetic theory
    k_b_ev = 8.617333262e-5
    e_act = max(float(config.catalyst_E_act_eV), 0.0)
    k_surf_m_s = (
        0.25
        * DEFAULT_STICKING_COEFF
        * v_thermal
        * float(np.exp(-e_act / max(k_b_ev * T_K, 1e-9)))
    )

    # Active area per unit pellet volume [m² active / m³ pellet]
    av_pellet = (6.0 / d_p) * active_area_multiplier(config)
    k_vol_1_s = av_pellet * k_surf_m_s

    # Spherical pellet characteristic length L_c = d_p / 6
    l_char = d_p / 6.0
    thiele_phi = float(l_char * np.sqrt(max(k_vol_1_s, 0.0) / d_eff))
    if thiele_phi < 1e-2:
        eta_internal = float(
            1.0 - (thiele_phi ** 2) / 15.0 + (2.0 * thiele_phi ** 4) / 315.0
        )
    elif thiele_phi > 200.0:
        eta_internal = float(3.0 / thiele_phi)
    else:
        eta_internal = float(
            (3.0 / (thiele_phi ** 2)) * (thiele_phi / np.tanh(thiele_phi) - 1.0)
        )
    eta_internal = min(max(eta_internal, 1e-6), 1.0)

    # External film mass transfer via Wakao-Funazkri (1979) packed-bed correlation
    rho_g = ch4_feed_density_kg_m3(T_K, P_Pa)
    mu_g = ch4_viscosity_pa_s(T_K)
    re_p = max(rho_g * max(u_sup_m_s, 1e-6) * d_p / max(mu_g, 1e-9), 1e-6)
    sc = max(mu_g / max(rho_g * d_mol, 1e-12), 1e-6)
    sherwood = float(2.0 + 1.1 * (re_p ** 0.6) * (sc ** (1.0 / 3.0)))
    k_c_m_s = max(sherwood * d_mol / d_p, 1e-9)

    da_external = float((eta_internal * k_surf_m_s) / k_c_m_s)
    omega_overall = min(max(eta_internal / (1.0 + da_external), 1e-6), 1.0)

    return {
        'thiele_modulus': thiele_phi,
        'internal_effectiveness_factor': eta_internal,
        'overall_effectiveness_factor': omega_overall,
        'external_damkohler_number': da_external,
        'external_mass_transfer_coeff_m_s': float(k_c_m_s),
        'effective_diffusivity_m2_s': float(d_eff),
        'knudsen_diffusivity_m2_s': float(d_knudsen),
        'molecular_diffusivity_m2_s': float(d_mol),
        'sherwood_number': sherwood,
        'particle_reynolds_number': float(re_p),
    }


def _stage_ergun_gradient_pa_m(
    T_K: float,
    P_Pa: float,
    u_sup_m_s: float,
    eps_bed: float,
    d_p_m: float,
) -> float:
    """Differential Ergun pressure gradient -dP/dz [Pa/m] at local bed state."""
    eps = min(max(float(eps_bed), MIN_PACKED_BED_POROSITY), 0.99)
    mu = ch4_viscosity_pa_s(T_K)
    rho = ch4_feed_density_kg_m3(T_K, P_Pa)
    u = max(float(u_sup_m_s), 1e-6)
    viscous = 150.0 * mu * ((1.0 - eps) ** 2) / ((eps ** 3) * (d_p_m ** 2)) * u
    inertial = 1.75 * rho * (1.0 - eps) / ((eps ** 3) * d_p_m) * (u ** 2)
    return float(viscous + inertial)


def simulate_pbr(
    config: ReactorConfig,
    *,
    cantera_available: bool,
    mock_result: Callable[[ReactorConfig, str], ReactorResult],
) -> ReactorResult:
    """Catalytic Packed-Bed Reactor (PBR) integrated on a catalyst-weight basis.

    Distinguishes a physical packed bed from an idealized plug-flow tube by:
      1. Formulating the bed on a differential catalyst-weight basis
         ``dW = (1 - ε_0)·ρ_cat·A_cross·dz`` with explicit WHSV [kg_feed/(kg_cat·h)].
      2. Coupling the differential Ergun equation ``-dP/dz`` stage-by-stage
         with gas-phase molar expansion ``(1 + y_CH4,0·X)``, decompression,
         and dynamic bed voidage ``ε(z)`` as solid carbon accumulates.
      3. Resolving intraparticle pore diffusion (Bosanquet molecular + Knudsen
         diffusivity, spherical-pellet Thiele modulus ``φ``, internal
         effectiveness factor ``η_int``) and external Wakao-Funazkri film mass
         transfer (overall effectiveness ``Ω <= 1``).

    Args:
        config: Configuration controlling this operation.
        cantera_available: Whether the configured Cantera backend is usable.
        mock_result: Explicit fallback used when Cantera is unavailable.

    Returns:
        Typed ReactorResult record with PBR axial profiles and diagnostics.
    """
    validate_carbon_policy(config)
    if not cantera_available:
        return mock_result(config, 'PBR')

    logger.info(f"Simulating PBR: {config.catalyst_name} at {config.T_inlet_K} K")
    gas, graphite, surf = load_gas_and_surface(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    eps_0 = float(config.bed_porosity)
    d_p_m = max(float(config.catalyst_particle_mm) * 1e-3, 1e-6)
    ch4_initial = float(species_x(gas, 'CH4') or 1.0)
    ar_initial = float(species_x(gas, 'Ar'))
    require_ar_tracer(ch4_initial, ar_initial)

    n_stages = 50
    bed_cross_area = float(np.pi * (config.bed_diameter_m / 2.0) ** 2)
    total_bed_volume = bed_cross_area * float(config.bed_length_m)
    total_catalyst_weight_kg = (
        (1.0 - eps_0) * float(config.catalyst_density_kg_m3) * total_bed_volume
    )
    stage_length = float(config.bed_length_m) / n_stages
    stage_bed_volume = bed_cross_area * stage_length
    stage_catalyst_weight_kg = total_catalyst_weight_kg / n_stages

    u_sup_inlet = float(config.gas_velocity_m_s if config.gas_velocity_m_s > 0 else 0.1)
    rho_inlet = ch4_feed_density_kg_m3(config.T_inlet_K, config.P_inlet_Pa)
    mass_flow_kg_s = rho_inlet * u_sup_inlet * bed_cross_area
    whsv_mass_h_inv = (
        (mass_flow_kg_s * 3600.0) / total_catalyst_weight_kg
        if total_catalyst_weight_kg > 0
        else 0.0
    )

    z_positions = np.linspace(0.0, config.bed_length_m, n_stages + 1)
    w_positions = np.linspace(0.0, total_catalyst_weight_kg, n_stages + 1)
    conversion_profile: list[float] = [0.0]
    temperature_profile: list[float] = [float(config.T_inlet_K)]
    pressure_profile: list[float] = [float(config.P_inlet_Pa)]
    bed_voidage_profile: list[float] = [eps_0]
    thiele_profile: list[float] = []
    effectiveness_profile: list[float] = []
    theta_C_tos: list[float] = [float(surface_coverage(surf, 'C_s') or 0.0)]

    cycles_completed = 0
    per_cycle_conversion: list[float] = []
    produce_time_s = 0.0
    tau_total_pass_s = 0.0
    pass_start_cov: list[npt.NDArray[np.float64] | None] = [None]
    inlet_pellet_metrics = pbr_pellet_effectiveness(
        config, config.T_inlet_K, config.P_inlet_Pa, u_sup_inlet, eps_0
    )

    def _advance_bed() -> None:
        nonlocal produce_time_s, tau_total_pass_s
        conversion_profile.clear()
        conversion_profile.append(0.0)
        temperature_profile.clear()
        temperature_profile.append(float(config.T_inlet_K))
        pressure_profile.clear()
        pressure_profile.append(float(config.P_inlet_Pa))
        bed_voidage_profile.clear()
        bed_voidage_profile.append(eps_0)
        thiele_profile.clear()
        effectiveness_profile.clear()
        theta_C_tos.clear()
        theta_C_tos.append(float(surface_coverage(surf, 'C_s') or 0.0))
        pass_start_cov[0] = (
            np.array(surf.coverages, dtype=np.float64)
            if surf is not None
            else None
        )

        p_local = float(config.P_inlet_Pa)
        tau_total_pass_s = 0.0

        for _ in range(n_stages):
            x_prev = conversion_profile[-1]
            t_local = float(gas.T)
            theta_c_curr = float(surface_coverage(surf, 'C_s') or 0.0)

            # Dynamic voidage reduction from surface/deposited carbon holdup
            c_holdup_vol_frac = min(
                0.5 * eps_0,
                0.08 * max(theta_c_curr, 0.0) + 0.02 * max(x_prev, 0.0),
            )
            eps_local = max(MIN_PACKED_BED_POROSITY, eps_0 - c_holdup_vol_frac)

            # Gas velocity expansion from CH4 -> C(s) + 2 H2 and decompression
            molar_expansion = 1.0 + ch4_initial * max(x_prev, 0.0)
            u_local = (
                u_sup_inlet
                * (config.P_inlet_Pa / max(p_local, MIN_STAGE_PRESSURE_PA))
                * (t_local / config.T_inlet_K)
                * molar_expansion
            )

            # Pellet Thiele modulus and overall effectiveness factor Ω_i
            pellet = pbr_pellet_effectiveness(
                config, t_local, p_local, u_local, eps_local
            )
            omega_i = pellet['overall_effectiveness_factor']
            thiele_profile.append(pellet['thiele_modulus'])
            effectiveness_profile.append(omega_i)

            # Stage geometric and effective active catalyst surface area
            sv_geom_local = 6.0 * (1.0 - eps_local) / d_p_m
            sv_active_local = active_sv(sv_geom_local, config)
            stage_gas_vol = stage_bed_volume * eps_local
            stage_eff_area = omega_i * sv_active_local * stage_bed_volume
            tau_stage = stage_length * eps_local / max(u_local, 1e-6)

            gas.TP = t_local, p_local
            advance_reactor_stage(
                gas, surf, stage_gas_vol, stage_eff_area, tau_stage
            )

            # Differential Ergun pressure drop across stage dz
            dp_dz = _stage_ergun_gradient_pa_m(
                t_local, p_local, u_local, eps_local, d_p_m
            )
            p_local = max(MIN_STAGE_PRESSURE_PA, p_local - dp_dz * stage_length)
            gas.TP = float(gas.T), p_local

            produce_time_s += tau_stage
            tau_total_pass_s += tau_stage
            conversion_profile.append(float(ch4_extent(gas, ch4_initial, ar_initial)))
            temperature_profile.append(float(gas.T))
            pressure_profile.append(p_local)
            bed_voidage_profile.append(eps_local)
            theta_C_tos.append(float(surface_coverage(surf, 'C_s') or 0.0))

    # Initial produce pass, followed by optional discrete non-oxidative regen cycles
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    if surf is not None:
        reset_surface_carbon(surf)
    _advance_bed()
    per_cycle_conversion.append(conversion_profile[-1])

    while (
        config.max_regen_cycles > 0
        and cycles_completed < config.max_regen_cycles
        and theta_C_tos
        and max(theta_C_tos) >= config.regen_coverage_threshold
    ):
        if config.regen_mechanism == REGEN_OXIDATIVE and not config.co2_permitted:
            raise RuntimeError('oxidative regen blocked (co2_permitted=False)')
        if config.regen_mechanism == REGEN_OXIDATIVE:
            logger.warning(
                'Oxidative PBR regen enabled via co2_permitted=True (test only)'
            )
        reset_surface_carbon(surf)
        gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
        cycles_completed += 1
        _advance_bed()
        per_cycle_conversion.append(conversion_profile[-1])

    final_conv = conversion_profile[-1]
    x_h2 = float(species_x(gas, 'H2'))

    # Carbon accounting on the stage parcel basis (last pass)
    sv_ratio_0 = active_sv(geometric_sv_pfr(config), config)
    stage_surface_area_0 = sv_ratio_0 * stage_bed_volume
    stage_gas_volume_0 = stage_bed_volume * eps_0
    gamma_mol_m2 = float(config.site_density_mol_cm2) * 1e4
    n_sites_mol = gamma_mol_m2 * stage_surface_area_0
    c_total = config.P_inlet_Pa / (R_J_MOL_K * config.T_inlet_K)
    n_parcel_in_mol = c_total * stage_gas_volume_0
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
        pass_time_s=tau_total_pass_s,
        cov_start=pass_start_cov[0],
        basis='Gamma*A_stage / (c_CH4*eps*V_stage), PBR differential catalyst-weight parcel',
    )

    inv_fields = solids_inventory_fields(config, geometric_sv_pfr(config))
    coupled_dp_pa = float(config.P_inlet_Pa - pressure_profile[-1])
    coupled_ergun_ok = bool(coupled_dp_pa <= ERGUN_DP_LIMIT_PA)
    inv_fields['ergun_delta_p_Pa'] = coupled_dp_pa
    inv_fields['ergun_delta_p_bar'] = coupled_dp_pa / 1e5
    inv_fields['ergun_ok'] = coupled_ergun_ok

    tau_in_band = bool(
        PFR_PILOT_TAU_BAND_S[0] <= tau_total_pass_s <= PFR_PILOT_TAU_BAND_S[1]
    )
    max_theta_c = float(max(theta_C_tos) if theta_C_tos else 0.0)
    overloaded_bed_warning = bool(not tau_in_band or not coupled_ergun_ok)
    plugging_risk = bool(
        max_theta_c >= config.regen_coverage_threshold or not coupled_ergun_ok
    )
    b2_carbon_policy_summary = (
        f"regen_mechanism={config.regen_mechanism};"
        f"regen_cycles={cycles_completed};"
        f"max_theta_C={max_theta_c:.3f};"
        f"plugging_risk={plugging_risk};"
        f"coupled_dp_bar={coupled_dp_pa / 1e5:.4f}"
    )

    z_list = [float(v) for v in z_positions.tolist()]
    w_list = [float(v) for v in w_positions.tolist()]
    kinetics_ev = kinetics_evidence(config)
    kinetics_ev_limitations = list(kinetics_ev['reactor_evidence_limitations'])
    if 'pbr_screening_model_not_validated' not in kinetics_ev_limitations:
        kinetics_ev_limitations.append('pbr_screening_model_not_validated')

    result = cast(
        ReactorResult,
        {
            'reactor_type': 'PBR',
            'bed_or_interface': 'catalytic_packed_bed_pellet',
            'cantera_reactor_model': 'staged_lagrangian_packed_bed_reactors',
            'catalyst_name': config.catalyst_name,
            'T_K': config.T_inlet_K,
            'P_Pa': config.P_inlet_Pa,
            'exit_P_Pa': float(pressure_profile[-1]),
            'coupled_ergun_delta_p_Pa': coupled_dp_pa,
            'coupled_ergun_delta_p_bar': coupled_dp_pa / 1e5,
            'bed_length_m': config.bed_length_m,
            'bed_diameter_m': config.bed_diameter_m,
            'catalyst_particle_mm': config.catalyst_particle_mm,
            'catalyst_weight_kg': float(total_catalyst_weight_kg),
            'stage_catalyst_weight_kg': float(stage_catalyst_weight_kg),
            'WHSV_mass_kg_feed_per_kg_cat_h': float(whsv_mass_h_inv),
            'residence_time_s': float(tau_total_pass_s),
            'residence_time_within_pilot_band': tau_in_band,
            'pilot_residence_time_band_s': list(PFR_PILOT_TAU_BAND_S),
            'overloaded_bed_warning': overloaded_bed_warning,
            'plugging_risk': plugging_risk,
            'b2_carbon_policy_summary': b2_carbon_policy_summary,
            'produce_time_s': float(produce_time_s),
            'WHSV_h-1': reciprocal_residence_h(tau_total_pass_s),
            'CH4_conversion': float(final_conv),
            'single_pass_CH4_conversion': float(
                per_cycle_conversion[0] if per_cycle_conversion else final_conv
            ),
            'cycle_averaged_CH4_conversion': float(
                np.mean(per_cycle_conversion) if per_cycle_conversion else final_conv
            ),
            'endothermic_heat_duty_kJ_per_mol_CH4_fed': float(
                CH4_PYROLYSIS_DELTA_H_KJ_MOL * final_conv
            ),
            'CH4_mole_fraction_drop': float(mole_fraction_drop(gas, ch4_initial)),
            'conversion_basis': 'argon_tracer',
            'thermal_mode': 'isothermal_energy_disabled',
            'surface_area_basis': (
                'overall_pellet_effectiveness_times_sv_per_bed_volume_times_bed_volume'
            ),
            'stage_surface_area_m2': float(
                inlet_pellet_metrics['overall_effectiveness_factor']
                * stage_surface_area_0
            ),
            'stage_gas_volume_m3': float(stage_gas_volume_0),
            'per_cycle_CH4_conversion': per_cycle_conversion,
            'regen_cycles_completed': cycles_completed,
            'exit_x_H2': x_h2,
            'exit_T_K': float(gas.T),
            'z_positions': z_list,
            'z_positions_m': z_list,
            'w_positions_kg': w_list,
            'conversion_profile': conversion_profile,
            'temperature_profile': temperature_profile,
            'pressure_profile_Pa': pressure_profile,
            'bed_voidage_profile': bed_voidage_profile,
            'exit_bed_porosity': float(bed_voidage_profile[-1]),
            'thiele_modulus_profile': thiele_profile,
            'effectiveness_factor_profile': effectiveness_profile,
            'mean_overall_effectiveness_factor': float(
                np.mean(effectiveness_profile) if effectiveness_profile else 1.0
            ),
            **inlet_pellet_metrics,
            'theta_C_time_on_stream': theta_C_tos,
            'theta_C_profile_basis': 'single_shared_surface_cumulative_time_on_stream',
            'theta_C_axial': theta_C_tos,
            'inlet_theta_C': float(theta_C_tos[1] if len(theta_C_tos) > 1 else 0.0),
            'max_theta_C': max_theta_c,
            'exit_theta_C': float(theta_C_tos[-1]) if theta_C_tos else None,
            **off_site,
            **kinetics_fields(config),
            **inv_fields,
            **kinetics_ev,
            'reactor_evidence_limitations': kinetics_ev_limitations,
            **policy_metadata(config),
            **load_status_fields(config, graphite, surf),
        },
    )
    logger.info(
        f"  PBR result: conversion={final_conv:.2%}, τ={tau_total_pass_s:.1f}s, "
        f"Ω={inlet_pellet_metrics['overall_effectiveness_factor']:.3f}, "
        f"ΔP={coupled_dp_pa / 1e5:.4f} bar, regen_cycles={cycles_completed}"
    )
    return result
