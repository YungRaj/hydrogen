"""Validated, unit-bearing physical inputs for multiphysics reactor cases."""

from __future__ import annotations

import json
import math
from pathlib import Path


CASE_SCHEMA_VERSION = 1

_REQUIRED = {
    'Fluidized': {
        'geometry': ('column_diameter_m', 'bed_height_m'),
        'operating': ('temperature_K', 'pressure_Pa', 'methane_mass_flow_kg_s'),
        'properties': ('particle_diameter_m', 'particle_density_kg_m3',
                       'gas_viscosity_Pa_s'),
        'models': ('drag_model', 'heat_transfer_model'),
    },
    'MMBCR': {
        'geometry': ('column_diameter_m', 'liquid_height_m',
                     'sparger_orifice_diameter_m'),
        'operating': ('temperature_K', 'pressure_Pa', 'methane_mass_flow_kg_s'),
        'properties': ('liquid_density_kg_m3', 'liquid_viscosity_Pa_s',
                       'surface_tension_N_m'),
        'models': ('bubble_breakup_model', 'bubble_coalescence_model'),
    },
    'NTEC': {
        'geometry': ('reactor_volume_m3', 'interface_area_m2'),
        'operating': ('temperature_K', 'pressure_Pa', 'methane_mass_flow_kg_s',
                      'shear_rate_s-1', 'mechanical_power_W_kg'),
        'properties': ('liquid_viscosity_Pa_s', 'permittivity_F_m',
                       'ionic_conductivity_S_m'),
        'models': ('contact_electrification_model', 'species_transport_model'),
    },
    'Electrochemical': {
        'geometry': ('electrode_area_m2', 'electrolyte_thickness_m'),
        'operating': ('temperature_K', 'pressure_Pa', 'methane_mass_flow_kg_s',
                      'applied_potential_V'),
        'properties': ('ionic_conductivity_S_m', 'electronic_conductivity_S_m',
                       'methane_diffusivity_m2_s'),
        'models': ('charge_transfer_model', 'species_transport_model'),
    },
}


def _positive_number(value) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool) and
            math.isfinite(float(value)) and float(value) > 0)


def load_physical_case(path: str | Path, *, candidate_id: str, mode: str,
                       reactor_type: str, temperature_K: float) -> dict:
    """Load a physical case and fail on missing physics or provenance."""
    source = Path(path)
    try:
        case = json.loads(source.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f'physical case is unreadable: {exc}') from exc
    if reactor_type not in _REQUIRED:
        raise ValueError(f'no physical-case contract for {reactor_type}')
    identity = {
        'schema_version': case.get('schema_version') == CASE_SCHEMA_VERSION,
        'candidate_id': case.get('candidate_id') == candidate_id,
        'pathway_mode': case.get('pathway_mode') == mode,
        'reactor_type': case.get('reactor_type') == reactor_type,
    }
    operating = case.get('operating', {})
    try:
        identity['temperature_K'] = abs(
            float(operating.get('temperature_K')) - float(temperature_K)) < 1e-6
    except (TypeError, ValueError):
        identity['temperature_K'] = False
    failures = [name for name, valid in identity.items() if not valid]
    parameter_sources = case.get('parameter_sources', {})
    for section, names in _REQUIRED[reactor_type].items():
        values = case.get(section, {})
        if not isinstance(values, dict):
            failures.append(section)
            continue
        for name in names:
            value = values.get(name)
            if section == 'models':
                if not isinstance(value, str) or not value.strip():
                    failures.append(f'{section}.{name}')
            elif not _positive_number(value):
                failures.append(f'{section}.{name}')
            elif not parameter_sources.get(f'{section}.{name}'):
                failures.append(f'parameter_sources.{section}.{name}')
    feed = case.get('feed', {})
    if not isinstance(feed, dict) or not feed.get('composition') or not feed.get('source'):
        failures.append('feed.composition_and_source')
    kinetics = case.get('kinetics', {})
    if not isinstance(kinetics, dict) or not kinetics.get('source'):
        failures.append('kinetics.source')
    calibration = case.get('calibration', {})
    training = calibration.get('training_ids', []) if isinstance(calibration, dict) else []
    validation = calibration.get('validation_ids', []) if isinstance(calibration, dict) else []
    if not training or not validation or set(training).intersection(validation):
        failures.append('calibration.disjoint_training_and_validation_ids')
    if not isinstance(calibration, dict) or not calibration.get('source'):
        failures.append('calibration.source')
    if reactor_type == 'NTEC' and calibration.get('paired_control') is not True:
        failures.append('calibration.paired_control')
    if reactor_type == 'Electrochemical' and case.get(
            'electrolyte_phase') not in {'aqueous', 'molten'}:
        failures.append('electrolyte_phase')
    if failures:
        raise ValueError('physical case failed: ' + ', '.join(sorted(set(failures))))
    return case


def case_summary(case: dict) -> dict:
    """Return traceable metadata without duplicating the full solver input."""
    calibration = case['calibration']
    return {
        'schema_version': case['schema_version'],
        'kinetics_source': case['kinetics']['source'],
        'feed_source': case['feed']['source'],
        'calibration_source': calibration['source'],
        'calibration_count': len(calibration['training_ids']),
        'holdout_validation_count': len(calibration['validation_ids']),
        'disjoint_holdout': True,
    }
