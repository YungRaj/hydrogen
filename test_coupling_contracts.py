#!/usr/bin/env python3
"""Contracts for external multiphysics solver handoffs and provenance."""

import hashlib
import json
import tempfile
from pathlib import Path

from pipeline.process.coupling_contract import (
    require_pristine_case, validate_hydrodynamic_handoff,
    validate_solver_coupling)
from pipeline.process.multiphysics_runner import _tree_digest


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _handoff(case):
    field = case / 'hydrodynamics.xdmf'
    field.write_text('<Xdmf/>\n')
    path = case / 'hydrogen_hydrodynamics.json'
    path.write_text(json.dumps({
        'schema_version': 1, 'candidate_id': 'candidate-1',
        'pathway_mode': 'ntec', 'reactor_type': 'NTEC',
        'temperature_K': 300.0,
        'fields': {'velocity_m_s': 0.1, 'pressure_Pa': 101325,
                   'temperature_K': 300.0, 'liquid_volume_fraction': 0.9,
                   'shear_rate_s_inv': 100.0},
        'field_artifact': {'path': field.name, 'sha256': _sha(field),
                           'format': 'XDMF', 'mesh_id': 'mesh-1',
                           'coordinate_system': 'cartesian-m'}}))
    return path


def _coupling(case):
    values = {
        'mechanism': ('mechanism.yaml', 'phases: []\n'),
        'cantera_log': ('cantera.log', 'Cantera 3 execution\n'),
        'rate_exchange': ('rates.json', json.dumps({
            'schema_version': 1, 'candidate_id': 'candidate-1',
            'temperature_K': 300.0,
            'reaction_rates_mol_m3_s': {'CH4': -0.1, 'H2': 0.2}})),
        'coupling_history': ('history.json', json.dumps({'iterations': [
            {'iteration': 1, 'residual_relative': 0.02},
            {'iteration': 2, 'residual_relative': 1e-5}]})),
    }
    proof = {
        'schema_version': 1, 'cantera_used': True,
        'coupling_method': 'iterative_two_way', 'coupling_iterations': 2,
        'coupling_residual_relative': 1e-5,
        'coupling_tolerance_relative': 1e-4,
        'exchanged_fields': ['species', 'temperature', 'reaction_heat',
                             'reaction_rates', 'momentum', 'charge',
                             'potential']}
    for stem, (name, content) in values.items():
        path = case / name
        path.write_text(content)
        proof[f'{stem}_path'] = name
        proof[f'{stem}_sha256'] = _sha(path)
    return proof


def _reject(call, phrase):
    try:
        call()
    except ValueError as exc:
        assert phrase in str(exc), str(exc)
    else:
        raise AssertionError(f'expected rejection containing {phrase}')


def test_hydrodynamic_handoff_checks_identity_fields_and_hash():
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        path = _handoff(case)
        valid = validate_hydrodynamic_handoff(
            case, candidate_id='candidate-1', mode='ntec',
            reactor_type='NTEC', temperature_K=300.0)
        assert valid['mesh_id'] == 'mesh-1'
        value = json.loads(path.read_text())
        value['candidate_id'] = 'wrong'
        path.write_text(json.dumps(value))
        _reject(lambda: validate_hydrodynamic_handoff(
            case, candidate_id='candidate-1', mode='ntec',
            reactor_type='NTEC', temperature_K=300.0), 'candidate_id')


def test_cantera_boolean_is_not_coupling_proof():
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        _reject(lambda: validate_solver_coupling(
            case, {'cantera_used': True}, candidate_id='candidate-1',
            temperature_K=300.0, reactor_type='NTEC'), 'proof missing')
        proof = _coupling(case)
        result = validate_solver_coupling(
            case, proof, candidate_id='candidate-1', temperature_K=300.0,
            reactor_type='NTEC')
        assert result['converged'] and result['iterations'] == 2
        proof['coupling_residual_relative'] = 0.1
        _reject(lambda: validate_solver_coupling(
            case, proof, candidate_id='candidate-1', temperature_K=300.0,
            reactor_type='NTEC'), 'not converged')


def test_pristine_case_rejects_outputs_and_openfoam_time_directories():
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        (case / '0').mkdir()
        require_pristine_case(case)
        (case / 'hydrogen_outputs.json').write_text('{}')
        _reject(lambda: require_pristine_case(case), 'stale generated outputs')
        (case / 'hydrogen_outputs.json').unlink()
        (case / '0.1').mkdir()
        _reject(lambda: require_pristine_case(case), '0.1')


def test_external_fenics_script_changes_input_digest():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        case = root / 'case'
        case.mkdir()
        (case / 'input').write_text('case input\n')
        script = root / 'model.py'
        script.write_text('value = 1\n')
        first = _tree_digest(case, (script,))
        script.write_text('value = 2\n')
        assert _tree_digest(case, (script,)) != first


def main():
    tests = [value for name, value in sorted(globals().items())
             if name.startswith('test_') and callable(value)]
    for test in tests:
        test()
        print('PASS', test.__name__)
    print(f'{len(tests)}/{len(tests)} coupling contracts passed')


if __name__ == '__main__':
    main()
