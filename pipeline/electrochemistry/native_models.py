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


def electrochemical_screening_point(parameters: dict[str, float]) -> dict[str, Any]:
    """Evaluate one explicit analytic electrochemical operating point.

    Required parameters are intentionally passed by name; there are no hidden
    electrolyte or catalyst defaults. The result is diagnostic and non-pruning.
    """
    exchange = exchange_current_density(
        parameters['descriptor_eV'], parameters['optimum_descriptor_eV'],
        parameters['reference_current_A_cm2'], parameters['temperature_K'],
    )
    current_density = butler_volmer_current_density(
        exchange, parameters['overpotential_V'], parameters['temperature_K'],
        parameters.get('symmetry_factor', 0.5), int(parameters.get('electrons', 2)),
    )
    area = parameters['area_m2']
    current = current_density * 1.0e4 * area
    limit = limiting_current_A(
        int(parameters.get('electrons', 2)), parameters['diffusivity_m2_s'],
        parameters['concentration_mol_m3'], area, parameters['boundary_layer_m'],
    )
    if current >= limit:
        return {
            'status': 'invalid', 'reactor_evidence_tier': 'native_screening_model',
            'can_exclude_candidate': False, 'reason': 'current_exceeds_limiting_current',
            'current_density_A_cm2': current_density, 'limiting_current_A': limit,
        }
    voltage = cell_voltage(CellOperatingPoint(
        parameters['reversible_voltage_V'], parameters['overpotential_V'],
        parameters.get('cathode_overpotential_V', 0.0), current, limit,
        ohmic_resistance_ohm(parameters['gap_m'], parameters['conductivity_S_m'], area),
        parameters['temperature_K'], int(parameters.get('electrons', 2)),
    ))
    return {
        'status': 'diagnostic', 'reactor_evidence_tier': 'native_screening_model',
        'can_exclude_candidate': False, 'current_density_A_cm2': current_density,
        'cell_voltage_V': voltage, 'electrical_power_density_W_cm2': voltage * current_density,
    }


def ntec_screening_point(parameters: dict[str, float]) -> dict[str, Any]:
    """Evaluate explicit collision, impact, and detachment diagnostics."""
    collisions = collision_events_per_volume_s(
        parameters['shear_rate_s'], parameters['particle_diameter_m'],
        parameters['number_density_m3'],
    )
    impact = impact_energy_J(
        parameters['effective_mass_kg'], parameters['shear_rate_s'],
        parameters['particle_diameter_m'],
    )
    return {
        'status': 'diagnostic', 'reactor_evidence_tier': 'native_screening_model',
        'can_exclude_candidate': False, 'collision_events_per_volume_s': collisions,
        'impact_energy_J': impact,
        'detachment_probability': detachment_probability(parameters['adhesion_energy_J'], impact),
    }
