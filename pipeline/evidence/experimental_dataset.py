"""Validate checksum-bound experimental calibration and holdout datasets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


COMMON = (
    'measurement_id', 'experimental_unit_id', 'candidate_id', 'application',
    'measurement_type', 'split', 'replicate_id', 'apparatus_id',
    'measured_at_utc', 'split_assigned_at_utc', 'split_assignment_source',
    'raw_source_path', 'raw_source_sha256', 'conditions', 'observations',
    'uncertainties', 'uncertainty_type',
)

REQUIRED = {
    'pyrolysis_reactor': {
        'application': 'turquoise_hydrogen',
        'conditions': ('temperature_K', 'pressure_Pa', 'inlet_CH4_mol_s',
                       'duration_h', 'feed_mole_fraction'),
        'observations': ('CH4_conversion_fraction', 'H2_selectivity_fraction',
                         'solid_carbon_yield_fraction',
                         'carbon_balance_closure_fraction',
                         'hydrogen_balance_closure_fraction',
                         'net_energy_kWh_kg_H2',
                         'deactivation_fraction_per_h'),
    },
    'ntec_pair': {
        'application': 'turquoise_hydrogen',
        'conditions': ('temperature_K', 'pressure_Pa', 'inlet_CH4_mol_s',
                       'duration_h', 'feed_mole_fraction', 'treatment',
                       'shear_rate_s_inv', 'mechanical_power_W'),
        'observations': ('CH4_conversion_fraction', 'H2_selectivity_fraction',
                         'solid_carbon_yield_fraction',
                         'carbon_balance_closure_fraction',
                         'hydrogen_balance_closure_fraction',
                         'net_energy_kWh_kg_H2',
                         'deactivation_fraction_per_h'),
    },
    'mea_performance': {
        'application': 'fuel_cell',
        'conditions': ('temperature_K', 'pressure_Pa', 'relative_humidity_anode',
                       'relative_humidity_cathode', 'duration_h'),
        'observations': ('peak_power_W_cm2', 'system_efficiency_fraction',
                         'orr_overpotential_V'),
    },
    'durability': {
        'application': 'fuel_cell',
        'conditions': ('temperature_K', 'pressure_Pa', 'duration_h'),
        'observations': ('voltage_degradation_uV_h',
                         'power_retention_fraction'),
    },
    'hydrogen_impurity': {
        'application': 'fuel_cell',
        'conditions': ('temperature_K', 'pressure_Pa', 'duration_h',
                       'hydrogen_impurities_mole_fraction'),
        'observations': ('peak_power_W_cm2', 'power_retention_fraction'),
    },
}

FRACTIONS = {
    'CH4_conversion_fraction', 'H2_selectivity_fraction',
    'solid_carbon_yield_fraction', 'system_efficiency_fraction',
    'carbon_balance_closure_fraction', 'hydrogen_balance_closure_fraction',
    'power_retention_fraction', 'relative_humidity_anode',
    'relative_humidity_cathode',
}


def _digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def _number(value, label: str, *, minimum: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{label} must be numeric') from exc
    if not math.isfinite(number) or number < minimum:
        raise ValueError(f'{label} must be finite and >= {minimum}')
    return number


def _utc(value: str, label: str) -> datetime:
    try:
        parsed = datetime.fromisoformat(str(value).replace('Z', '+00:00'))
    except ValueError as exc:
        raise ValueError(f'{label} must be ISO-8601') from exc
    if parsed.tzinfo is None:
        raise ValueError(f'{label} must include a timezone')
    return parsed


def _validate_values(record: dict, index: int) -> None:
    label = f'records[{index}]'
    kind = record['measurement_type']
    if kind not in REQUIRED:
        raise ValueError(f'{label}.measurement_type is unsupported')
    spec = REQUIRED[kind]
    if record['application'] != spec['application']:
        raise ValueError(f'{label}.application does not match measurement_type')
    for section in ('conditions', 'observations'):
        values = record[section]
        if not isinstance(values, dict):
            raise ValueError(f'{label}.{section} must be an object')
        missing = [key for key in spec[section] if key not in values]
        if missing:
            raise ValueError(f'{label}.{section} missing {missing}')
        for key in spec[section]:
            if key in ('treatment', 'hydrogen_impurities_mole_fraction',
                       'feed_mole_fraction'):
                continue
            number = _number(values[key], f'{label}.{section}.{key}')
            if key in FRACTIONS and number > 1:
                raise ValueError(f'{label}.{section}.{key} must be <= 1')
    if _number(record['conditions']['temperature_K'], label + '.temperature_K') <= 0:
        raise ValueError(f'{label}.temperature_K must be positive')
    if _number(record['conditions']['pressure_Pa'], label + '.pressure_Pa') <= 0:
        raise ValueError(f'{label}.pressure_Pa must be positive')
    if _number(record['conditions']['duration_h'], label + '.duration_h') <= 0:
        raise ValueError(f'{label}.duration_h must be positive')
    uncertainties = record['uncertainties']
    if set(uncertainties) != set(spec['observations']):
        raise ValueError(f'{label}.uncertainties must cover every observation exactly')
    for key, value in uncertainties.items():
        _number(value, f'{label}.uncertainties.{key}')
    if kind == 'ntec_pair' and record['conditions']['treatment'] not in ('ntec', 'control'):
        raise ValueError(f'{label}.conditions.treatment must be ntec or control')
    if kind in ('pyrolysis_reactor', 'ntec_pair'):
        feed = record['conditions']['feed_mole_fraction']
        if not isinstance(feed, dict) or 'CH4' not in feed:
            raise ValueError(f'{label} feed must include CH4 mole fraction')
        total = sum(_number(value, f'{label}.feed.{key}')
                    for key, value in feed.items())
        if not math.isclose(total, 1.0, rel_tol=0, abs_tol=1e-6):
            raise ValueError(f'{label} feed mole fractions must sum to 1')
    if kind == 'ntec_pair':
        shear = float(record['conditions']['shear_rate_s_inv'])
        power = float(record['conditions']['mechanical_power_W'])
        treatment = record['conditions']['treatment']
        if treatment == 'control' and (shear != 0 or power != 0):
            raise ValueError(f'{label} NTEC control intervention must be zero')
        if treatment == 'ntec' and (shear <= 0 or power <= 0):
            raise ValueError(f'{label} NTEC intervention must be positive')
    if kind in ('pyrolysis_reactor', 'ntec_pair') and float(
            record['observations']['solid_carbon_yield_fraction']) > float(
                record['observations']['CH4_conversion_fraction']):
        raise ValueError(f'{label} carbon yield cannot exceed methane conversion')
    if kind == 'hydrogen_impurity':
        composition = record['conditions']['hydrogen_impurities_mole_fraction']
        if not isinstance(composition, dict) or not composition:
            raise ValueError(f'{label} requires a nonempty impurity composition')
        total = sum(_number(value, f'{label}.impurity.{key}')
                    for key, value in composition.items())
        if total >= 1:
            raise ValueError(f'{label} impurity mole fractions must sum below 1')


def validate_experimental_dataset(path: str | Path) -> dict:
    """Validate raw provenance, physical fields, split isolation, and NTEC pairs."""
    dataset_path = Path(path).expanduser().resolve()
    payload = json.loads(dataset_path.read_text())
    for key in ('schema_version', 'dataset_id', 'protocol_id', 'records',
                'run_log_path', 'run_log_sha256', 'analysis_plan'):
        if payload.get(key) in (None, '', []):
            raise ValueError(f'dataset missing {key}')
    if payload['schema_version'] != 1 or not isinstance(payload['records'], list):
        raise ValueError('experimental dataset schema mismatch')
    run_log = Path(payload['run_log_path']).expanduser()
    if not run_log.is_absolute():
        run_log = dataset_path.parent / run_log
    if not run_log.is_file() or _digest(run_log) != str(
            payload['run_log_sha256']).lower():
        raise ValueError('run log is missing or checksum mismatches')
    with run_log.open(newline='') as handle:
        attempts = list(csv.DictReader(handle))
    required_attempt = {'measurement_id', 'experimental_unit_id', 'split',
                        'status', 'failure_reason'}
    if not attempts or not required_attempt.issubset(attempts[0]):
        raise ValueError('run log schema is incomplete')
    attempt_ids = [row['measurement_id'] for row in attempts]
    if any(not value for value in attempt_ids) or len(attempt_ids) != len(set(attempt_ids)):
        raise ValueError('run log measurement IDs must be nonempty and unique')
    for row in attempts:
        if row['split'] not in ('calibration', 'holdout'):
            raise ValueError('run log split is invalid')
        if row['status'] not in ('completed', 'failed', 'excluded'):
            raise ValueError('run log status is invalid')
        if row['status'] != 'completed' and not row['failure_reason'].strip():
            raise ValueError('failed/excluded run log entries require a reason')
    plan = payload['analysis_plan']
    required_plan = ('plan_id', 'locked_at_utc', 'source_path', 'source_sha256',
                     'minimum_calibration_units', 'minimum_holdout_units')
    missing_plan = [key for key in required_plan if plan.get(key) in (None, '')]
    if missing_plan:
        raise ValueError(f'analysis_plan missing {missing_plan}')
    plan_source = Path(plan['source_path']).expanduser()
    if not plan_source.is_absolute():
        plan_source = dataset_path.parent / plan_source
    if not plan_source.is_file() or _digest(plan_source) != str(
            plan['source_sha256']).lower():
        raise ValueError('analysis plan is missing or checksum mismatches')
    plan_locked = _utc(plan['locked_at_utc'], 'analysis_plan.locked_at_utc')
    minimum_calibration_raw = _number(
        plan['minimum_calibration_units'], 'minimum_calibration_units', minimum=1)
    minimum_holdout_raw = _number(
        plan['minimum_holdout_units'], 'minimum_holdout_units', minimum=1)
    if not minimum_calibration_raw.is_integer() or not minimum_holdout_raw.is_integer():
        raise ValueError('minimum independent-unit counts must be integers')
    minimum_calibration, minimum_holdout = (
        int(minimum_calibration_raw), int(minimum_holdout_raw))
    identities, unit_splits, pairs = set(), defaultdict(set), defaultdict(list)
    counts = Counter()
    for index, record in enumerate(payload['records']):
        label = f'records[{index}]'
        missing = [key for key in COMMON if record.get(key) in (None, '')]
        if missing:
            raise ValueError(f'{label} missing {missing}')
        identity = str(record['measurement_id'])
        if identity in identities:
            raise ValueError(f'{label} duplicates measurement_id')
        identities.add(identity)
        split = str(record['split'])
        if split not in ('calibration', 'holdout'):
            raise ValueError(f'{label}.split must be calibration or holdout')
        assigned = _utc(record['split_assigned_at_utc'], label + '.split_assigned_at_utc')
        measured = _utc(record['measured_at_utc'], label + '.measured_at_utc')
        if split == 'holdout' and assigned > measured:
            raise ValueError(f'{label} holdout was assigned after measurement')
        if split == 'holdout' and plan_locked > measured:
            raise ValueError(f'{label} analysis plan was locked after measurement')
        if split == 'holdout' and record.get('blinded') is not True:
            raise ValueError(f'{label} holdout must be explicitly blinded')
        unit_splits[str(record['experimental_unit_id'])].add(split)
        source = Path(record['raw_source_path']).expanduser()
        if not source.is_absolute():
            source = dataset_path.parent / source
        if not source.is_file() or _digest(source) != str(
                record['raw_source_sha256']).lower():
            raise ValueError(f'{label} raw source is missing or checksum mismatches')
        _validate_values(record, index)
        if record['uncertainty_type'] not in (
                'standard_deviation', 'standard_error', 'expanded_uncertainty'):
            raise ValueError(f'{label}.uncertainty_type is unsupported')
        counts[(record['candidate_id'], record['measurement_type'], split)] += 1
        if record['measurement_type'] == 'ntec_pair':
            if not record.get('control_pair_id'):
                raise ValueError(f'{label} ntec_pair requires control_pair_id')
            pairs[str(record['control_pair_id'])].append(record)
    leaking = sorted(unit for unit, splits in unit_splits.items() if len(splits) > 1)
    if leaking:
        raise ValueError(f'experimental units leak across splits: {leaking}')
    completed_ids = {row['measurement_id'] for row in attempts
                     if row['status'] == 'completed'}
    if completed_ids != identities:
        raise ValueError('completed run log IDs must exactly match dataset records')
    for row in attempts:
        unit_splits[row['experimental_unit_id']].add(row['split'])
    leaking = sorted(unit for unit, splits in unit_splits.items() if len(splits) > 1)
    if leaking:
        raise ValueError(f'run-log experimental units leak across splits: {leaking}')
    for pair_id, records in pairs.items():
        if len(records) != 2 or {r['conditions']['treatment'] for r in records} != {
                'ntec', 'control'}:
            raise ValueError(f'NTEC pair {pair_id} must contain one ntec and one control')
        invariant = ('candidate_id', 'apparatus_id', 'split')
        if any(len({str(r[key]) for r in records}) != 1 for key in invariant):
            raise ValueError(f'NTEC pair {pair_id} identity or split mismatch')
        if len({str(r['experimental_unit_id']) for r in records}) != 2:
            raise ValueError(f'NTEC pair {pair_id} requires distinct experimental units')
        baseline = ('temperature_K', 'pressure_Pa', 'inlet_CH4_mol_s', 'duration_h')
        if any(len({float(r['conditions'][key]) for r in records}) != 1
               for key in baseline):
            raise ValueError(f'NTEC pair {pair_id} baseline conditions mismatch')
        if len({json.dumps(r['conditions']['feed_mole_fraction'], sort_keys=True)
                for r in records}) != 1:
            raise ValueError(f'NTEC pair {pair_id} feed composition mismatch')
    groups = []
    for candidate, kind in sorted({(c, k) for c, k, _ in counts}):
        relevant = [record for record in payload['records']
                    if record['candidate_id'] == candidate and
                    record['measurement_type'] == kind]
        calibration = len({record['experimental_unit_id'] for record in relevant
                           if record['split'] == 'calibration'})
        holdout = len({record['experimental_unit_id'] for record in relevant
                       if record['split'] == 'holdout'})
        groups.append({'candidate_id': candidate, 'measurement_type': kind,
                       'independent_calibration_units': calibration,
                       'independent_holdout_units': holdout,
                       'split_complete': calibration >= minimum_calibration and
                                         holdout >= minimum_holdout})
    return {
        'valid': True, 'dataset_id': payload['dataset_id'],
        'protocol_id': payload['protocol_id'], 'record_count': len(identities),
        'ntec_pair_count': len(pairs), 'groups': groups,
        'split_complete_groups': sum(group['split_complete'] for group in groups),
        'minimum_calibration_units': minimum_calibration,
        'minimum_holdout_units': minimum_holdout,
        'run_log_sha256': payload['run_log_sha256'],
        'failed_or_excluded_attempts': sum(
            row['status'] != 'completed' for row in attempts),
        'analysis_plan_sha256': plan['source_sha256'],
        'dataset_sha256': _digest(dataset_path),
    }


def experimental_evidence_records(path: str | Path) -> dict:
    """Create evidence-manifest fragments only for split-complete groups."""
    report = validate_experimental_dataset(path)
    mapping = {
        'pyrolysis_reactor': ('measured_reactor', 'measured_deactivation'),
        'ntec_pair': ('ntec_control_pair', 'measured_reactor',
                      'measured_deactivation'),
        'mea_performance': ('measured_mea',),
        'durability': ('measured_durability',),
        'hydrogen_impurity': ('hydrogen_impurity_test',),
    }
    records = defaultdict(list)
    for group in report['groups']:
        if not group['split_complete']:
            continue
        for evidence_type in mapping[group['measurement_type']]:
            records[evidence_type].append({
                'candidate_id': group['candidate_id'],
                'source_path': str(Path(path).expanduser().resolve()),
                'sha256': report['dataset_sha256'],
                'protocol_id': report['protocol_id'],
                'status': 'measured',
                'dataset_id': report['dataset_id'],
                'measurement_type': group['measurement_type'],
            })
    return dict(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dataset')
    parser.add_argument('--output')
    parser.add_argument('--evidence-output')
    args = parser.parse_args()
    report = validate_experimental_dataset(args.dataset)
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        Path(args.output).write_text(rendered)
    if args.evidence_output:
        Path(args.evidence_output).write_text(json.dumps(
            experimental_evidence_records(args.dataset), indent=2,
            sort_keys=True) + '\n')
    print(rendered, end='')


if __name__ == '__main__':
    main()
