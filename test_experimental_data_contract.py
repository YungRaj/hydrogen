#!/usr/bin/env python3
"""Regression tests for experimental evidence and split integrity."""

import hashlib
import json
import tempfile
from copy import deepcopy
from pathlib import Path

from pipeline.evidence.experimental_dataset import (
    experimental_evidence_records, validate_experimental_dataset)


def _record(root, measurement_id, unit_id, split='calibration', kind='pyrolysis_reactor'):
    raw = root / f'{measurement_id}.csv'
    raw.write_text('time_s,ch4_mol_s,h2_mol_s\n0,0.001,0\n')
    record = {
        'measurement_id': measurement_id, 'experimental_unit_id': unit_id,
        'candidate_id': 'candidate-1', 'application': 'turquoise_hydrogen',
        'measurement_type': kind, 'split': split, 'replicate_id': 'rep-1',
        'apparatus_id': 'rig-1', 'measured_at_utc': '2026-01-02T00:00:00Z',
        'split_assigned_at_utc': '2026-01-01T00:00:00Z',
        'split_assignment_source': 'randomization-v1',
        'raw_source_path': raw.name,
        'raw_source_sha256': hashlib.sha256(raw.read_bytes()).hexdigest(),
        'uncertainty_type': 'standard_deviation',
        'conditions': {'temperature_K': 1000, 'pressure_Pa': 101325,
                       'inlet_CH4_mol_s': 0.001, 'duration_h': 10,
                       'feed_mole_fraction': {'CH4': 0.95, 'Ar': 0.05}},
        'observations': {'CH4_conversion_fraction': 0.75,
                         'H2_selectivity_fraction': 0.97,
                         'solid_carbon_yield_fraction': 0.70,
                         'carbon_balance_closure_fraction': 0.99,
                         'hydrogen_balance_closure_fraction': 0.98,
                         'net_energy_kWh_kg_H2': 14,
                         'deactivation_fraction_per_h': 0.005},
        'uncertainties': {'CH4_conversion_fraction': 0.02,
                          'H2_selectivity_fraction': 0.01,
                          'solid_carbon_yield_fraction': 0.03,
                          'carbon_balance_closure_fraction': 0.01,
                          'hydrogen_balance_closure_fraction': 0.01,
                          'net_energy_kWh_kg_H2': 0.5,
                          'deactivation_fraction_per_h': 0.001},
    }
    if split == 'holdout':
        record['blinded'] = True
    return record


def _write(root, records):
    run_log = root / 'run_log.csv'
    run_log.write_text(
        'measurement_id,experimental_unit_id,split,status,failure_reason\n' +
        ''.join(f"{record['measurement_id']},{record['experimental_unit_id']},"
                f"{record['split']},completed,\n" for record in records))
    plan = root / 'analysis_plan.txt'
    plan.write_text('Locked analysis plan fixture\n')
    path = root / 'dataset.json'
    path.write_text(json.dumps({'schema_version': 1, 'dataset_id': 'dataset-1',
                                'protocol_id': 'protocol-1',
                                'run_log_path': run_log.name,
                                'run_log_sha256': hashlib.sha256(
                                    run_log.read_bytes()).hexdigest(),
                                'analysis_plan': {
                                    'plan_id': 'plan-1',
                                    'locked_at_utc': '2025-12-31T00:00:00Z',
                                    'source_path': plan.name,
                                    'source_sha256': hashlib.sha256(
                                        plan.read_bytes()).hexdigest(),
                                    'minimum_calibration_units': 1,
                                    'minimum_holdout_units': 1},
                                'records': records}))
    return path


def _reject(path, phrase):
    try:
        validate_experimental_dataset(path)
    except ValueError as exc:
        assert phrase in str(exc), str(exc)
    else:
        raise AssertionError(f'expected rejection containing {phrase}')


def test_calibration_and_blinded_holdout_are_ready():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        path = _write(root, [_record(root, 'cal-1', 'batch-1'),
                             _record(root, 'hold-1', 'batch-2', 'holdout')])
        report = validate_experimental_dataset(path)
        assert report['valid'] and report['split_complete_groups'] == 1
        assert report['groups'][0]['independent_calibration_units'] == 1
        assert report['groups'][0]['independent_holdout_units'] == 1
        evidence = experimental_evidence_records(path)
        assert len(evidence['measured_reactor']) == 1
        assert len(evidence['measured_deactivation']) == 1


def test_experimental_unit_cannot_leak_across_split():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        path = _write(root, [_record(root, 'cal-1', 'same-batch'),
                             _record(root, 'hold-1', 'same-batch', 'holdout')])
        _reject(path, 'leak across splits')


def test_holdout_must_be_preassigned_and_blinded():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        record = _record(root, 'hold-1', 'batch-1', 'holdout')
        record['split_assigned_at_utc'] = '2026-01-03T00:00:00Z'
        _reject(_write(root, [record]), 'assigned after measurement')
        record['split_assigned_at_utc'] = '2026-01-01T00:00:00Z'
        record['blinded'] = False
        _reject(_write(root, [record]), 'explicitly blinded')


def test_ntec_requires_matched_treatment_and_control():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        control = _record(root, 'control-1', 'pair-unit', kind='ntec_pair')
        control['control_pair_id'] = 'pair-1'
        control['conditions']['treatment'] = 'control'
        control['conditions']['shear_rate_s_inv'] = 0
        control['conditions']['mechanical_power_W'] = 0
        treatment = deepcopy(control)
        treatment['measurement_id'] = 'ntec-1'
        treatment['experimental_unit_id'] = 'pair-unit-treatment'
        raw = root / 'ntec-1.csv'
        raw.write_text('time_s,ch4_mol_s,h2_mol_s\n0,0.001,0\n')
        treatment['raw_source_path'] = raw.name
        treatment['raw_source_sha256'] = hashlib.sha256(raw.read_bytes()).hexdigest()
        treatment['conditions']['treatment'] = 'ntec'
        treatment['conditions']['shear_rate_s_inv'] = 100
        treatment['conditions']['mechanical_power_W'] = 10
        report = validate_experimental_dataset(_write(root, [control, treatment]))
        assert report['ntec_pair_count'] == 1
        broken = deepcopy(treatment)
        broken['conditions']['temperature_K'] = 999
        _reject(_write(root, [control, broken]), 'baseline conditions mismatch')


def test_checksum_uncertainty_and_carbon_bounds_fail_closed():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        record = _record(root, 'cal-1', 'batch-1')
        record['raw_source_sha256'] = '0' * 64
        _reject(_write(root, [record]), 'checksum mismatches')
        record = _record(root, 'cal-2', 'batch-2')
        record['uncertainties'].pop('CH4_conversion_fraction')
        _reject(_write(root, [record]), 'cover every observation exactly')
        record = _record(root, 'cal-3', 'batch-3')
        record['observations']['solid_carbon_yield_fraction'] = 0.8
        _reject(_write(root, [record]), 'cannot exceed methane conversion')


def test_run_log_preserves_failures_and_blocks_omitted_successes():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        record = _record(root, 'cal-1', 'batch-1')
        path = _write(root, [record])
        run_log = root / 'run_log.csv'
        with run_log.open('a') as handle:
            handle.write('failed-1,batch-failed,calibration,failed,leak detected\n')
        payload = json.loads(path.read_text())
        payload['run_log_sha256'] = hashlib.sha256(run_log.read_bytes()).hexdigest()
        path.write_text(json.dumps(payload))
        assert validate_experimental_dataset(path)['failed_or_excluded_attempts'] == 1
        with run_log.open('a') as handle:
            handle.write('omitted-1,batch-2,calibration,completed,\n')
        payload['run_log_sha256'] = hashlib.sha256(run_log.read_bytes()).hexdigest()
        path.write_text(json.dumps(payload))
        _reject(path, 'exactly match dataset records')


def main():
    tests = [value for name, value in sorted(globals().items())
             if name.startswith('test_') and callable(value)]
    for test in tests:
        test()
        print('PASS', test.__name__)
    print(f'{len(tests)}/{len(tests)} experimental-data contracts passed')


if __name__ == '__main__':
    main()
