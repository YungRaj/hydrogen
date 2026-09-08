"""Validated, unit-bearing physical inputs for multiphysics reactor cases."""

from __future__ import annotations

import json
import math
import argparse
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
                      'shear_rate_s_inv', 'mechanical_power_W_kg'),
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


def _declared(value) -> bool:
    return (isinstance(value, str) and bool(value.strip()) and
            'REPLACE_' not in value.upper())


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
    from pipeline.process.pathway_modes import reactor_types_for_mode
    try:
        routed = reactor_type in reactor_types_for_mode(mode)
    except ValueError:
        routed = False
    identity = {
        'schema_version': case.get('schema_version') == CASE_SCHEMA_VERSION,
        'candidate_id': case.get('candidate_id') == candidate_id,
        'pathway_mode': case.get('pathway_mode') == mode,
        'reactor_type': case.get('reactor_type') == reactor_type,
        'mode_reactor_routing': routed,
    }
    operating = case.get('operating', {})
    try:
        identity['temperature_K'] = abs(
            float(operating.get('temperature_K')) - float(temperature_K)) < 1e-6
    except (TypeError, ValueError):
        identity['temperature_K'] = False
    failures = [name for name, valid in identity.items() if not valid]
    if case.get('template') is True:
        failures.append('template_must_be_completed')
    parameter_sources = case.get('parameter_sources', {})
    for section, names in _REQUIRED[reactor_type].items():
        values = case.get(section, {})
        if not isinstance(values, dict):
            failures.append(section)
            continue
        for name in names:
            value = values.get(name)
            if section == 'models':
                if not _declared(value):
                    failures.append(f'{section}.{name}')
            elif not _positive_number(value):
                failures.append(f'{section}.{name}')
            elif not _declared(parameter_sources.get(f'{section}.{name}')):
                failures.append(f'parameter_sources.{section}.{name}')
    feed = case.get('feed', {})
    composition = feed.get('composition', {}) if isinstance(feed, dict) else {}
    if (not composition or not _declared(feed.get('source')) or
            any(not isinstance(value, (int, float)) or value < 0
                for value in composition.values()) or
            sum(composition.values()) <= 0):
        failures.append('feed.composition_and_source')
    kinetics = case.get('kinetics', {})
    if not isinstance(kinetics, dict) or not _declared(kinetics.get('source')):
        failures.append('kinetics.source')
    calibration = case.get('calibration', {})
    training = calibration.get('training_ids', []) if isinstance(calibration, dict) else []
    validation = calibration.get('validation_ids', []) if isinstance(calibration, dict) else []
    if (not training or not validation or set(training).intersection(validation)
            or any(not _declared(value) for value in training + validation)):
        failures.append('calibration.disjoint_training_and_validation_ids')
    if not isinstance(calibration, dict) or not _declared(calibration.get('source')):
        failures.append('calibration.source')
    try:
        threshold = float(calibration.get('acceptance_threshold'))
        if (calibration.get('metric') not in {'rmse', 'mae', 'relative_rmse'} or
                not math.isfinite(threshold) or threshold < 0):
            failures.append('calibration.metric_and_acceptance_threshold')
    except (TypeError, ValueError):
        failures.append('calibration.metric_and_acceptance_threshold')
    if reactor_type == 'NTEC' and calibration.get('paired_control') is not True:
        failures.append('calibration.paired_control')
    if reactor_type == 'Electrochemical' and case.get(
            'electrolyte_phase') not in {'aqueous', 'molten'}:
        failures.append('electrolyte_phase')
    if reactor_type == 'Electrochemical':
        electrolyte = case.get('electrolyte', {})
        if (not _declared(electrolyte.get('identity')) or
                not _declared(electrolyte.get('source'))):
            failures.append('electrolyte.identity_and_source')
        if (case.get('electrolyte_phase') == 'aqueous' and
                _positive_number(operating.get('temperature_K')) and
                float(operating['temperature_K']) >= 647.096):
            failures.append('aqueous_temperature_below_water_critical_point')
    geometry = case.get('geometry', {})
    properties = case.get('properties', {})
    if reactor_type == 'Fluidized' and all(_positive_number(value) for value in (
            geometry.get('column_diameter_m'), properties.get('particle_diameter_m'))):
        if properties['particle_diameter_m'] >= geometry['column_diameter_m']:
            failures.append('particle_smaller_than_column')
    if reactor_type == 'MMBCR' and all(_positive_number(value) for value in (
            geometry.get('column_diameter_m'), geometry.get('sparger_orifice_diameter_m'))):
        if geometry['sparger_orifice_diameter_m'] >= geometry['column_diameter_m']:
            failures.append('sparger_orifice_smaller_than_column')
    if failures:
        raise ValueError('physical case failed: ' + ', '.join(sorted(set(failures))))
    return case


def case_template(*, candidate_id: str, mode: str, reactor_type: str,
                  temperature_K: float, electrolyte_phase: str | None = None) -> dict:
    """Create a deliberately non-runnable template with all required keys."""
    from pipeline.process.pathway_modes import reactor_types_for_mode
    if reactor_type not in reactor_types_for_mode(mode):
        raise ValueError(f'{reactor_type} is not routed by mode {mode}')
    sections = {}
    defaults = {'temperature_K': temperature_K, 'pressure_Pa': 101325.0}
    for section, fields in _REQUIRED[reactor_type].items():
        sections[section] = {
            name: ('REPLACE_WITH_MODEL' if section == 'models' else
                   defaults.get(name, 1.0)) for name in fields}
    sources = {
        f'{section}.{name}': 'REPLACE_WITH_SOURCE'
        for section in ('geometry', 'operating', 'properties')
        for name in sections[section]}
    result = {
        'schema_version': CASE_SCHEMA_VERSION, 'template': True,
        'candidate_id': candidate_id, 'pathway_mode': mode,
        'reactor_type': reactor_type, **sections,
        'feed': {'composition': {'CH4': 1.0}, 'source': 'REPLACE_WITH_SOURCE'},
        'kinetics': {'source': 'REPLACE_WITH_SOURCE'},
        'parameter_sources': sources,
        'calibration': {
            'training_ids': ['REPLACE_TRAIN_ID'],
            'validation_ids': ['REPLACE_HOLDOUT_ID'],
            'source': 'REPLACE_WITH_SOURCE', 'metric': 'relative_rmse',
            'acceptance_threshold': 0.1},
    }
    if reactor_type == 'NTEC':
        result['calibration']['paired_control'] = True
    if reactor_type == 'Electrochemical':
        result['electrolyte_phase'] = electrolyte_phase or 'aqueous'
        result['electrolyte'] = {
            'identity': 'REPLACE_WITH_ELECTROLYTE',
            'source': 'REPLACE_WITH_SOURCE'}
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    init = sub.add_parser('init')
    for command in (init,):
        command.add_argument('--candidate-id', required=True)
        command.add_argument('--mode', required=True)
        command.add_argument('--reactor-type', required=True,
                             choices=sorted(_REQUIRED))
        command.add_argument('--temperature-K', required=True, type=float)
    init.add_argument('--output', required=True)
    init.add_argument('--electrolyte-phase', choices=('aqueous', 'molten'))
    validate = sub.add_parser('validate')
    validate.add_argument('path')
    validate.add_argument('--candidate-id', required=True)
    validate.add_argument('--mode', required=True)
    validate.add_argument('--reactor-type', required=True)
    validate.add_argument('--temperature-K', required=True, type=float)
    args = parser.parse_args()
    if args.command == 'init':
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        target = output / 'hydrogen_case.json'
        if target.exists():
            raise SystemExit(f'refusing to overwrite {target}')
        target.write_text(json.dumps(case_template(
            candidate_id=args.candidate_id, mode=args.mode,
            reactor_type=args.reactor_type, temperature_K=args.temperature_K,
            electrolyte_phase=args.electrolyte_phase), indent=2) + '\n')
        print(target)
    else:
        try:
            load_physical_case(
                args.path, candidate_id=args.candidate_id, mode=args.mode,
                reactor_type=args.reactor_type, temperature_K=args.temperature_K)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        print('physical case: ready')


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


if __name__ == '__main__':
    main()
