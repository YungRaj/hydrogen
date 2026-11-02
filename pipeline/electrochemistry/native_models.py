"""Provisional native screening points for electrochemical and NTEC pathways."""

from __future__ import annotations

from typing import Any

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


def electrochemical_screening_point(parameters: dict[str, float]) -> dict[str, Any]:
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


def ntec_screening_point(parameters: dict[str, float]) -> dict[str, Any]:
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
    result: dict[str, Any] = {
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
