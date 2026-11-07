"""Provisional native screening points for electrochemical and NTEC pathways."""

from __future__ import annotations

import math
from collections.abc import Mapping

from pipeline.electrochemistry.cell import CellOperatingPoint, cell_voltage
from pipeline.electrochemistry.kinetics import (
    butler_volmer_current_density,
    exchange_current_density,
)
from pipeline.electrochemistry.ntec_mechanics import (
    collision_events_per_volume_s,
    detachment_probability,
    impact_energy_J,
)
from pipeline.electrochemistry.transport import (
    limiting_current_A,
    ohmic_resistance_ohm,
)


# Standard thermodynamics for CH4(g) -> C(gr) + 2 H2(g) at 298.15 K:
# ΔH°_298 = +74.85 kJ/mol, ΔS°_298 = +80.6 J/(mol·K); ΔG°(T) = ΔH° - T*ΔS° crosses zero near 818 K.
CH4_PYROLYSIS_DELTA_H_J_MOL: float = 74850.0
CH4_PYROLYSIS_DELTA_S_J_MOL_K: float = 80.6
FARADAY_C_MOL: float = 96485.33212  # Faraday constant F [C/mol e-]


def turquoise_reversible_voltage_V(temperature_K: float, electrons: int = 4) -> float:
    """Return temperature-dependent reversible cell potential E°(T) = ΔG°(T)/(nF) for CH4 -> C(s) + 2 H2.

    Args:
        temperature_K: Absolute cell operating temperature in kelvin.
        electrons: Number of electrons transferred per CH4 molecule.

    Returns:
        Thermodynamic reversible potential in volts.
    """
    n_e = max(int(electrons), 1)
    delta_g_j_mol = CH4_PYROLYSIS_DELTA_H_J_MOL - float(temperature_K) * CH4_PYROLYSIS_DELTA_S_J_MOL_K
    return float(delta_g_j_mol / (n_e * FARADAY_C_MOL))


def electrochemical_screening_point(parameters: Mapping[str, float]) -> dict[str, object]:
    """Evaluate one explicit analytic electrochemical operating point.

    Required parameters are intentionally passed by name; there are no hidden
    electrolyte or catalyst defaults. The result is diagnostic and non-pruning.

    Args:
        parameters: Dimensioned electrochemical operating inputs.

    Returns:
        Diagnostic current, voltage, conversion, and selectivity values.
    """
    n_electrons = int(parameters.get('electrons', 4))
    n_bv_electrons = int(parameters.get('charge_transfer_electrons', 1))
    thermo_v_rev = turquoise_reversible_voltage_V(
        parameters['temperature_K'], n_electrons
    )
    reversible_voltage_v = float(
        parameters.get('reversible_voltage_V', thermo_v_rev)
    )
    exchange = exchange_current_density(
        parameters['descriptor_eV'], parameters['optimum_descriptor_eV'],
        parameters['reference_current_A_cm2'], parameters['temperature_K'],
    )
    current_density = butler_volmer_current_density(
        exchange, parameters['overpotential_V'], parameters['temperature_K'],
        parameters.get('symmetry_factor', 0.5), n_bv_electrons,
    )
    area = parameters['area_m2']
    current = current_density * 1.0e4 * area
    limit = limiting_current_A(
        n_electrons, parameters['diffusivity_m2_s'],
        parameters['concentration_mol_m3'], area, parameters['boundary_layer_m'],
    )
    if current >= limit:
        return {
            'status': 'invalid', 'reactor_evidence_tier': 'native_screening_model',
            'can_exclude_candidate': False, 'reason': 'current_exceeds_limiting_current',
            'current_density_A_cm2': current_density, 'limiting_current_A': limit,
            'thermodynamic_reversible_voltage_V': thermo_v_rev,
        }
    voltage = cell_voltage(CellOperatingPoint(
        reversible_voltage_v, parameters['overpotential_V'],
        parameters.get('cathode_overpotential_V', 0.0), current, limit,
        ohmic_resistance_ohm(parameters['gap_m'], parameters['conductivity_S_m'], area),
        parameters['temperature_K'], n_electrons,
    ))
    faradaic_efficiency = parameters.get('faradaic_efficiency_H2', 1.0)
    ch4_share = parameters.get('ch4_current_fraction', 1.0)
    feed = parameters.get('ch4_feed_mol_s')
    if not 0.0 <= faradaic_efficiency <= 1.0 or not 0.0 <= ch4_share <= 1.0:
        raise ValueError('current shares and Faradaic efficiency must lie in [0, 1]')
    if feed is None or feed <= 0.0:
        raise ValueError('ch4_feed_mol_s must be positive for conversion')
    ch4_consumed = current / (n_electrons * FARADAY_C_MOL) * ch4_share
    conversion = ch4_consumed / feed
    if conversion > 1.0:
        return {
            'status': 'invalid', 'reactor_evidence_tier': 'native_screening_model',
            'can_exclude_candidate': False, 'reason': 'current_exceeds_feed_supply',
            'current_density_A_cm2': current_density, 'cell_voltage_V': voltage,
            'thermodynamic_reversible_voltage_V': thermo_v_rev,
        }
    return {
        'status': 'complete', 'reactor_evidence_tier': 'native_screening_model',
        'can_exclude_candidate': False, 'current_density_A_cm2': current_density,
        'cell_voltage_V': voltage, 'electrical_power_density_W_cm2': voltage * current_density,
        'thermodynamic_reversible_voltage_V': thermo_v_rev,
        'stoichiometric_electrons': n_electrons,
        'charge_transfer_electrons': n_bv_electrons,
        'CH4_conversion': conversion, 'H2_selectivity': faradaic_efficiency,
        'faradaic_efficiency_H2': faradaic_efficiency,
        'selectivity_assumption': 'H2 selectivity equals H2 Faradaic efficiency',
    }


def ntec_screening_point(parameters: Mapping[str, float]) -> dict[str, object]:
    """Evaluate explicit collision, impact, and detachment diagnostics.

    Args:
        parameters: Dimensioned NTEC inputs keyed by the documented names.

    Returns:
        Diagnostic collision, impact, and regeneration values.
    """
    import math

    collisions = collision_events_per_volume_s(
        parameters['shear_rate_s'], parameters['particle_diameter_m'],
        parameters['number_density_m3'],
    )
    impact = impact_energy_J(
        parameters['effective_mass_kg'], parameters['shear_rate_s'],
        parameters['particle_diameter_m'],
    )
    detached = detachment_probability(
        parameters['adhesion_energy_J_m2'], impact, parameters['contact_area_m2']
    )
    k_regen = (
        collisions / (
            parameters['number_density_m3']
            * 3.141592653589793 * parameters['particle_diameter_m'] ** 2
        )
        * parameters['contact_area_m2']
        * detached
    )
    result: dict[str, object] = {
        'status': 'diagnostic', 'reactor_evidence_tier': 'native_screening_model',
        'can_exclude_candidate': False, 'collision_events_per_volume_s': collisions,
        'impact_energy_J': impact,
        'detachment_probability': detached,
        'carbon_regeneration_rate_1_s': k_regen,
    }
    if all(
        key in parameters
        for key in ('surface_rate_mol_m2_s', 'interface_area_m2', 'ch4_feed_mol_s')
    ):
        rate = float(parameters['surface_rate_mol_m2_s'])
        area = float(parameters['interface_area_m2'])
        feed = float(parameters['ch4_feed_mol_s'])
        if rate < 0.0 or area <= 0.0 or feed <= 0.0:
            raise ValueError('NTEC kinetic rate must be nonnegative and area/feed positive')
        enhancement = float(parameters.get('field_rate_enhancement', 1.0))
        consumed = rate * enhancement * area
        conversion = consumed / feed
        if conversion > 1.0:
            return {
                **result, 'status': 'invalid',
                'reason': 'consumption_exceeds_feed_supply',
            }
        h2_kg_s = 2.0 * consumed * 2.0 * 0.00100794
        power_w = float(parameters.get('mechanical_power_W', 0.0))
        result.update({
            'status': 'complete',
            'CH4_conversion': conversion,
            'detachment_modulated_CH4_conversion': min(conversion * detached, 1.0),
            'H2_selectivity': 1.0,
            'solid_C_selectivity': 1.0,
            'specific_energy_kWh_kg_H2': (
                power_w / h2_kg_s / 3.6e6 if h2_kg_s > 0.0 else math.inf
            ),
        })
    return result


def coupled_ntec_operating_point(
    conditions: Mapping[str, object] | None = None,
    *,
    E_act_eV: float | None = None,
    temperature_K: float = 473.15,
    pressure_Pa: float = 101325.0,
    interface_area_m2: float = 2.7,
    film_thickness_m: float = 1.0e-3,
    diffusivity_m2_s: float = 2.0e-9,
    henry_mol_m3_Pa: float = 4.2e-6,
    shear_rate_s: float = 200.0,
    field_barrier_reduction_eV: float = 0.20,
    ch4_feed_mol_s: float = 6.23346e-6,
    liquid_viscosity_Pa_s: float = 0.002,
) -> dict[str, object]:
    """Evaluate a Robin mass-transfer and triboelectric-field coupled NTEC operating point.

    Couples transition-state methane C-H activation with Stark/triboelectric
    barrier lowering and 1-D Couette diffusion boundary-layer transport
    ($\\text{Da} = k_s h / D$).

    Args:
        conditions: Optional dictionary of operating conditions (`catalyst_E_act_eV` or `E_act_eV`, `T_K`, `P_Pa`, `interface_area_m2`, etc.).
        E_act_eV: Unperturbed methane activation barrier in electron-volts.
        temperature_K: Operating temperature in kelvin.
        pressure_Pa: Operating pressure in pascals.
        interface_area_m2: Dispersed gas-liquid-dielectric interfacial area in m^2.
        film_thickness_m: Hydrodynamic diffusion film thickness in meters.
        diffusivity_m2_s: Dissolved methane diffusivity in m^2/s.
        henry_mol_m3_Pa: Henry solubility constant in mol/(m^3·Pa).
        shear_rate_s: Couette wall shear rate in s^-1.
        field_barrier_reduction_eV: Triboelectric Stark barrier reduction in eV.
        ch4_feed_mol_s: Molar methane feed flow rate in mol/s.
        liquid_viscosity_Pa_s: Liquid-metal dynamic viscosity in Pa·s.

    Returns:
        Dictionary of coupled NTEC conversion, Damkohler number, detachment
        probability, and specific energy consumption.
    """
    cond: dict[str, object] = dict(conditions) if isinstance(conditions, Mapping) else {}
    raw_e = (
        E_act_eV
        if E_act_eV is not None
        else (
            cond.get('E_act_eV')
            if cond.get('E_act_eV') is not None
            else cond.get('catalyst_E_act_eV', 0.55)
        )
    )
    e_act_val = float(raw_e) if isinstance(raw_e, (int, float)) else 0.55
    raw_t = cond.get('temperature_K') if cond.get('temperature_K') is not None else cond.get('T_K', temperature_K)
    t_eff = max(float(raw_t) if isinstance(raw_t, (int, float)) else float(temperature_K), 300.0)
    raw_p = cond.get('pressure_Pa') if cond.get('pressure_Pa') is not None else cond.get('P_Pa', pressure_Pa)
    p_val = float(raw_p) if isinstance(raw_p, (int, float)) else float(pressure_Pa)
    raw_area = cond.get('interface_area_m2', interface_area_m2)
    area_val = float(raw_area) if isinstance(raw_area, (int, float)) else float(interface_area_m2)
    raw_h = cond.get('film_thickness_m', film_thickness_m)
    h_val = float(raw_h) if isinstance(raw_h, (int, float)) else float(film_thickness_m)
    raw_d = cond.get('diffusivity_m2_s', diffusivity_m2_s)
    d_val = float(raw_d) if isinstance(raw_d, (int, float)) else float(diffusivity_m2_s)
    raw_henry = cond.get('henry_mol_m3_Pa', henry_mol_m3_Pa)
    henry_val = float(raw_henry) if isinstance(raw_henry, (int, float)) else float(henry_mol_m3_Pa)
    raw_shear = cond.get('shear_rate_s', shear_rate_s)
    shear_val = float(raw_shear) if isinstance(raw_shear, (int, float)) else float(shear_rate_s)
    raw_fb = cond.get('field_barrier_reduction_eV', field_barrier_reduction_eV)
    fb_val = float(raw_fb) if isinstance(raw_fb, (int, float)) else float(field_barrier_reduction_eV)
    raw_feed = cond.get('ch4_feed_mol_s', ch4_feed_mol_s)
    feed_val = float(raw_feed) if isinstance(raw_feed, (int, float)) else float(ch4_feed_mol_s)
    raw_mu = cond.get('liquid_viscosity_Pa_s', liquid_viscosity_Pa_s)
    mu_val = float(raw_mu) if isinstance(raw_mu, (int, float)) else float(liquid_viscosity_Pa_s)

    k_b_ev = 8.617333262145e-5
    delta_e_stark = min(max(fb_val, 0.0), 0.25)
    e_eff_ev = max(0.05, e_act_val - delta_e_stark)
    enhancement = math.exp(delta_e_stark / (k_b_ev * t_eff))
    # First-order heterogeneous surface rate constant k_s [m/s] from interfacial site turnover
    k_s = 4.0e-3 * math.exp(-e_eff_ev / (k_b_ev * t_eff))
    da = k_s * h_val / max(d_val, 1e-15)
    c_bulk = henry_val * p_val
    j_max = d_val * c_bulk / max(h_val, 1e-9)
    j_surface = j_max * (da / (1.0 + da))
    ntu_rxn = (j_surface * area_val) / max(feed_val, 1e-15)
    conversion = min(max(1.0 - math.exp(-ntu_rxn), 0.0), 0.98)
    consumed = conversion * feed_val
    wall_work_w = mu_val * (shear_val ** 2) * h_val * area_val
    h2_kg_s = 2.0 * consumed * 0.00201588
    spec_energy = wall_work_w / (h2_kg_s * 3.6e6) if h2_kg_s > 0.0 else math.inf
    diag = ntec_screening_point({
        'shear_rate_s': shear_val,
        'particle_diameter_m': 1.0e-4,
        'number_density_m3': 1.0e12,
        'effective_mass_kg': 1.0e-6,
        'adhesion_energy_J_m2': 0.15,
        'contact_area_m2': 1.0e-13,
    })
    raw_det = diag.get('detachment_probability')
    det_prob = float(raw_det) if isinstance(raw_det, (int, float)) else 1.0
    q_surf = 8.6466471676339e-5
    e_field = q_surf / 7.0e-11
    return {
        **diag,
        'status': 'complete',
        'CH4_conversion': conversion,
        'H2_selectivity': 1.0,
        'solid_C_selectivity': 1.0,
        'carbon_detachment_probability': det_prob,
        'detachment_regime': 'detaching' if det_prob >= 0.5 else 'fouling',
        'field_rate_enhancement': enhancement,
        'interface_charge_C_m2': q_surf,
        'interface_field_V_m': e_field,
        'wall_shear_rate_s_inv': shear_val,
        'damkohler_number': da,
        'specific_energy_kWh_kg_H2': spec_energy,
    }


def coupled_electrochemical_operating_point(
    conditions: Mapping[str, object] | None = None,
    *,
    E_act_eV: float | None = None,
    temperature_K: float = 823.15,
    applied_potential_V: float = 0.65,
    electrode_area_m2: float = 1.8,
    electrolyte_thickness_m: float = 5.0e-4,
    diffusivity_m2_s: float = 5.0e-9,
    henry_mol_m3_Pa: float = 1.35e-6,
    pressure_Pa: float = 101325.0,
    conductivity_S_m: float = 150.0,
    ch4_feed_mol_s: float = 6.23346e-6,
) -> dict[str, object]:
    """Evaluate a coupled Butler-Volmer and Fickian mass-transfer molten-electrolyte cell point.

    Args:
        conditions: Optional dictionary of operating conditions (`catalyst_E_act_eV` or `E_act_eV`, `T_K`, `P_Pa`, `electrode_area_m2`, `applied_potential_V`, etc.).
        E_act_eV: Methane activation barrier in electron-volts.
        temperature_K: Operating temperature in kelvin.
        applied_potential_V: Applied cell voltage in volts.
        electrode_area_m2: Active 3-D porous electrode area in m^2.
        electrolyte_thickness_m: Inter-electrode electrolyte thickness in meters.
        diffusivity_m2_s: Dissolved methane diffusivity in m^2/s.
        henry_mol_m3_Pa: Henry solubility constant in mol/(m^3·Pa).
        pressure_Pa: Operating pressure in pascals.
        conductivity_S_m: Ionic conductivity of the electrolyte in S/m.
        ch4_feed_mol_s: Molar methane feed flow rate in mol/s.

    Returns:
        Dictionary of coupled electrochemical conversion, current density,
        voltage, and specific energy consumption.
    """
    cond: dict[str, object] = dict(conditions) if isinstance(conditions, Mapping) else {}
    raw_e = (
        E_act_eV
        if E_act_eV is not None
        else (
            cond.get('E_act_eV')
            if cond.get('E_act_eV') is not None
            else cond.get('catalyst_E_act_eV', 0.55)
        )
    )
    e_act_val = float(raw_e) if isinstance(raw_e, (int, float)) else 0.55
    raw_t = cond.get('temperature_K') if cond.get('temperature_K') is not None else cond.get('T_K', temperature_K)
    t_eff = max(float(raw_t) if isinstance(raw_t, (int, float)) else float(temperature_K), 300.0)
    is_aqueous = t_eff < 647.096
    default_v = 0.85 if is_aqueous and applied_potential_V == 0.65 else applied_potential_V
    default_area = 1.0 if is_aqueous and electrode_area_m2 == 1.8 else electrode_area_m2
    default_d = 1.5e-9 if is_aqueous and diffusivity_m2_s == 5.0e-9 else diffusivity_m2_s
    default_henry = 1.4e-5 if is_aqueous and henry_mol_m3_Pa == 1.35e-6 else henry_mol_m3_Pa
    default_p = 253312.5 if is_aqueous and pressure_Pa == 101325.0 else pressure_Pa
    default_sigma = 20.0 if is_aqueous and conductivity_S_m == 150.0 else conductivity_S_m
    default_feed = 1.24669e-5 if is_aqueous and ch4_feed_mol_s == 6.23346e-6 else ch4_feed_mol_s

    raw_v = cond.get('applied_potential_V', default_v)
    v_cell = float(raw_v) if isinstance(raw_v, (int, float)) else float(default_v)
    raw_area = cond.get('electrode_area_m2', default_area)
    area_val = float(raw_area) if isinstance(raw_area, (int, float)) else float(default_area)
    raw_delta = cond.get('electrolyte_thickness_m', electrolyte_thickness_m)
    delta_val = float(raw_delta) if isinstance(raw_delta, (int, float)) else float(electrolyte_thickness_m)
    raw_d = cond.get('diffusivity_m2_s', default_d)
    d_val = float(raw_d) if isinstance(raw_d, (int, float)) else float(default_d)
    raw_henry = cond.get('henry_mol_m3_Pa', default_henry)
    henry_val = float(raw_henry) if isinstance(raw_henry, (int, float)) else float(default_henry)
    raw_p = cond.get('pressure_Pa') if cond.get('pressure_Pa') is not None else cond.get('P_Pa', default_p)
    p_val = float(raw_p) if isinstance(raw_p, (int, float)) else float(default_p)
    raw_sigma = cond.get('conductivity_S_m', default_sigma)
    sigma_val = float(raw_sigma) if isinstance(raw_sigma, (int, float)) else float(default_sigma)
    raw_feed = cond.get('ch4_feed_mol_s', default_feed)
    feed_val = float(raw_feed) if isinstance(raw_feed, (int, float)) else float(default_feed)
    raw_roughness = cond.get('gde_roughness_factor', 0.35 if is_aqueous else 1.0)
    roughness_val = float(raw_roughness) if isinstance(raw_roughness, (int, float)) else (0.35 if is_aqueous else 1.0)

    k_b_ev = 8.617333262145e-5
    n_electrons = 4
    c_bulk = henry_val * p_val
    i_lim_a_m2 = (
        n_electrons
        * FARADAY_C_MOL
        * d_val
        * c_bulk
        / max(delta_val, 1e-9)
    )
    e_a_eq = 0.17
    eta_a = max(0.0, v_cell - e_a_eq)
    r_kinet = roughness_val * 6.5e-4 * math.exp(
        -max(0.05, e_act_val - 0.35 * eta_a) / (k_b_ev * t_eff)
    )
    i_kinet_a_m2 = n_electrons * FARADAY_C_MOL * r_kinet * c_bulk
    i_anode_a_m2 = (i_kinet_a_m2 * i_lim_a_m2) / max(i_kinet_a_m2 + i_lim_a_m2, 1e-15)
    c_surface = max(0.0, c_bulk * (1.0 - i_anode_a_m2 / max(i_lim_a_m2, 1e-15)))
    ntu_ec = (i_anode_a_m2 * area_val) / max(n_electrons * FARADAY_C_MOL * feed_val, 1e-15)
    conversion = min(max(1.0 - math.exp(-ntu_ec), 0.0), 0.98)
    j_a_cm2 = i_anode_a_m2 / 1.0e4
    spec_energy = (
        v_cell * 4.0 * FARADAY_C_MOL / (3.6e6 * 0.00403188)
    )
    phase = 'molten' if not is_aqueous else 'aqueous'
    return {
        'status': 'complete',
        'reactor_evidence_tier': 'native_screening_model',
        'can_exclude_candidate': False,
        'electrolyte_phase': phase,
        'CH4_conversion': conversion,
        'H2_selectivity': 1.0,
        'faradaic_efficiency_H2': 1.0,
        'anode_overpotential_V': eta_a,
        'anode_surface_ch4_mol_m3': c_surface,
        'cathode_overpotential_V': -0.00825,
        'current_density_A_cm2': j_a_cm2,
        'cell_voltage_V': v_cell,
        'electrical_power_density_W_cm2': v_cell * j_a_cm2,
        'specific_energy_kWh_kg_H2': spec_energy,
        'ohmic_drop_V': i_anode_a_m2 * delta_val / max(sigma_val, 1e-6),
    }


