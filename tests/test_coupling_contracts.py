#!/usr/bin/env python3
"""Contracts for external multiphysics solver handoffs and provenance."""

import hashlib
import json
import tempfile
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.simulation.solver_handoff import (
    require_pristine_case,
    validate_hydrodynamic_handoff,
    validate_solver_coupling,
)
from pipeline.simulation.external_runner import _tree_digest


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _handoff(case):
    field = case / 'hydrodynamics.xdmf'
    field.write_text('<Xdmf/>\n')
    path = case / 'hydrogen_hydrodynamics.json'
    path.write_text(
        json.dumps(
            {
                'schema_version': 1,
                'candidate_id': 'candidate-1',
                'pathway_mode': 'ntec',
                'reactor_type': 'NTEC',
                'temperature_K': 300.0,
                'fields': {
                    'velocity_m_s': 0.1,
                    'pressure_Pa': 101325,
                    'temperature_K': 300.0,
                    'liquid_volume_fraction': 0.9,
                    'shear_rate_s_inv': 100.0,
                },
                'field_artifact': {
                    'path': field.name,
                    'sha256': _sha(field),
                    'format': 'XDMF',
                    'mesh_id': 'mesh-1',
                    'coordinate_system': 'cartesian-m',
                },
            }
        )
    )
    return path


def _coupling(case):
    values = {
        'mechanism': ('mechanism.yaml', 'phases: []\n'),
        'cantera_log': ('cantera.log', 'Cantera 3 execution\n'),
        'rate_exchange': (
            'rates.json',
            json.dumps(
                {
                    'schema_version': 1,
                    'candidate_id': 'candidate-1',
                    'temperature_K': 300.0,
                    'reaction_rates_mol_m3_s': {'CH4': -0.1, 'H2': 0.2},
                }
            ),
        ),
        'coupling_history': (
            'history.json',
            json.dumps(
                {
                    'iterations': [
                        {'iteration': 1, 'residual_relative': 0.02},
                        {'iteration': 2, 'residual_relative': 1e-5},
                    ]
                }
            ),
        ),
    }
    proof = {
        'schema_version': 1,
        'cantera_used': True,
        'coupling_method': 'iterative_two_way',
        'coupling_iterations': 2,
        'coupling_residual_relative': 1e-5,
        'coupling_tolerance_relative': 1e-4,
        'exchanged_fields': [
            'species',
            'temperature',
            'reaction_heat',
            'reaction_rates',
            'momentum',
            'charge',
            'potential',
        ],
    }
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
            case,
            candidate_id='candidate-1',
            mode='ntec',
            reactor_type='NTEC',
            temperature_K=300.0,
        )
        assert valid['mesh_id'] == 'mesh-1'
        value = json.loads(path.read_text())
        value['candidate_id'] = 'wrong'
        path.write_text(json.dumps(value))
        _reject(
            lambda: validate_hydrodynamic_handoff(
                case,
                candidate_id='candidate-1',
                mode='ntec',
                reactor_type='NTEC',
                temperature_K=300.0,
            ),
            'candidate_id',
        )


def test_cantera_boolean_is_not_coupling_proof():
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        _reject(
            lambda: validate_solver_coupling(
                case,
                {'cantera_used': True},
                candidate_id='candidate-1',
                temperature_K=300.0,
                reactor_type='NTEC',
            ),
            'proof missing',
        )
        proof = _coupling(case)
        result = validate_solver_coupling(
            case,
            proof,
            candidate_id='candidate-1',
            temperature_K=300.0,
            reactor_type='NTEC',
        )
        assert result['converged'] and result['iterations'] == 2
        proof['coupling_residual_relative'] = 0.1
        _reject(
            lambda: validate_solver_coupling(
                case,
                proof,
                candidate_id='candidate-1',
                temperature_K=300.0,
                reactor_type='NTEC',
            ),
            'not converged',
        )


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


def test_external_moose_input_changes_input_digest():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        case = root / 'case'
        case.mkdir()
        (case / 'input').write_text('case input\n')
        model = root / 'model.i'
        model.write_text('value = 1\n')
        first = _tree_digest(case, (model,))
        model.write_text('value = 2\n')
        assert _tree_digest(case, (model,)) != first


def test_moose_command_passes_case_values_as_overrides_and_records_them():
    from pipeline.simulation.moose_backend import build_moose_command, moose_case_parameters

    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        physical = {
            'geometry': {'electrode_area_m2': 0.01},
            'operating': {'temperature_K': 353.15, 'pressure_Pa': 101325},
            'properties': {'ionic_conductivity_S_m': 20.0},
            'model_parameters': {'electrons_per_ch4': 4, 'label_only': 'ignored'},
        }
        parameters = moose_case_parameters(physical, case, 1)
        assert parameters['coupling_iteration'] == 1 and 'label_only' not in parameters
        command = build_moose_command('/bin/moose-opt', Path('/m/cell.i'), case, parameters)
        assert command[:3] == ['/bin/moose-opt', '-i', '/m/cell.i']
        assert '--allow-unused' in command
        assert f'Outputs/file_base={case / "hydrogen_moose"}' in command
        assert 'temperature_K=353.15' in command and 'electrons_per_ch4=4.0' in command
        record = json.loads((case / 'hydrogen_moose_parameters.json').read_text())
        assert record['command'] == command


def test_moose_case_parameters_reject_ambiguous_or_missing_exchange():
    from pipeline.simulation.moose_backend import moose_case_parameters, moose_overrides

    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        physical = {'geometry': {'x_m': 1.0}, 'properties': {'x_m': 2.0}}
        _reject(lambda: moose_case_parameters(physical, case, 1), 'declared twice')
        single = {'geometry': {'x_m': 1.0}}
        try:
            moose_case_parameters(single, case, 2)
        except RuntimeError as exc:
            assert 'coupling feedback missing' in str(exc)
        else:
            raise AssertionError('iteration 2 ran without Cantera feedback')
        (case / 'hydrogen_interface_rates.json').write_text(
            json.dumps({'moose_parameters': {'x_m': 3.0}}))
        _reject(lambda: moose_case_parameters(single, case, 2), 'overrides case field')
        (case / 'hydrogen_interface_rates.json').write_text(
            json.dumps({'moose_parameters': {'k_m_s': 3.0}}))
        assert moose_case_parameters(single, case, 2)['k_m_s'] == 3.0
    _reject(lambda: moose_overrides({'bad name': 1.0}), 'invalid MOOSE parameter name')
    _reject(lambda: moose_overrides({'x': float('nan')}), 'finite')


def test_malformed_multiphysics_artifact_is_rejected_not_raised():
    from pipeline.simulation.result_contract import artifact_path, load_validated_artifact

    with tempfile.TemporaryDirectory() as tmp:
        path = artifact_path(tmp, 'c', 'ntec', 'NTEC', 300.0)
        path.parent.mkdir(parents=True)
        for document in ([], 'text', {'provenance': 'not-an-object'},
                         {'convergence': [], 'provenance': {}}, {'backend_solvers': 3}):
            path.write_text(json.dumps(document))
            result = load_validated_artifact(tmp, 'c', 'ntec', 'NTEC', 300.0)
            assert result['valid'] is False, document
            assert result['reason'] == 'multiphysics_artifact_invalid', result


def main():
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith('test_') and callable(value)
    ]
    for test in tests:
        test()
        print('PASS', test.__name__)
    print(f'{len(tests)}/{len(tests)} coupling contracts passed')


if __name__ == '__main__':
    main()
