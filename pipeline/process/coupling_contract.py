"""Verifiable file contracts between OpenFOAM, FEniCSx, and Cantera."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


GENERATED_FILES = {
    'hydrogen_outputs.json', 'hydrogen_convergence.json',
    'hydrogen_metadata.json', 'hydrogen_hydrodynamics.json',
    'hydrogen_coupling_state.json', 'hydrogen_feedback.json',
    'hydrogen_coupling_request.json',
    'hydrogen_openfoam.stdout.log', 'hydrogen_openfoam.stderr.log',
    'hydrogen_fenicsx.stdout.log', 'hydrogen_fenicsx.stderr.log',
}


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def case_file(case: Path, relative: str, label: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ValueError(f'{label} must be a relative case path')
    root = case.resolve()
    target = (root / relative).resolve()
    if not target.is_relative_to(root) or not target.is_file():
        raise ValueError(f'{label} is missing or escapes the case directory')
    return target


def require_pristine_case(case: str | Path) -> None:
    """Reject outputs from an earlier run before calculating the input hash."""
    root = Path(case)
    stale = sorted(name for name in GENERATED_FILES if (root / name).exists())
    stale.extend(path.name for pattern in ('hydrogen_*.stdout.log',
                                           'hydrogen_*.stderr.log')
                 for path in root.glob(pattern))
    numeric_times = []
    for child in root.iterdir():
        if not child.is_dir() or child.name == '0':
            continue
        try:
            float(child.name)
        except ValueError:
            continue
        numeric_times.append(child.name)
    if stale or numeric_times:
        raise ValueError(
            f'case contains stale generated outputs: {stale + sorted(numeric_times)}')


def validate_hydrodynamic_handoff(case: str | Path, *, candidate_id: str,
                                  mode: str, reactor_type: str,
                                  temperature_K: float,
                                  iteration: int | None = None) -> dict:
    root = Path(case).resolve()
    path = root / 'hydrogen_hydrodynamics.json'
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f'hydrodynamic handoff unreadable: {exc}') from exc
    checks = {
        'schema_version': value.get('schema_version') == 1,
        'candidate_id': value.get('candidate_id') == candidate_id,
        'pathway_mode': value.get('pathway_mode') == mode,
        'reactor_type': value.get('reactor_type') == reactor_type,
    }
    if iteration is not None:
        checks['iteration'] = value.get('iteration') == iteration
    try:
        checks['temperature_K'] = abs(
            float(value['temperature_K']) - float(temperature_K)) < 1e-6
        fields = value['fields']
        for name in ('velocity_m_s', 'pressure_Pa', 'temperature_K',
                     'liquid_volume_fraction', 'shear_rate_s_inv'):
            number = float(fields[name])
            checks[name] = math.isfinite(number)
        checks['pressure_positive'] = float(fields['pressure_Pa']) > 0
        checks['temperature_positive'] = float(fields['temperature_K']) > 0
        checks['field_temperature_matches'] = abs(
            float(fields['temperature_K']) - float(temperature_K)) < 1e-6
        checks['liquid_fraction'] = 0 <= float(
            fields['liquid_volume_fraction']) <= 1
        checks['shear_nonnegative'] = float(fields['shear_rate_s_inv']) >= 0
    except (KeyError, TypeError, ValueError):
        checks['field_schema'] = False
    artifact = value.get('field_artifact', {})
    try:
        field_path = case_file(root, artifact['path'], 'hydrodynamic field artifact')
        checks['field_format'] = artifact.get('format') in {'XDMF', 'HDF5', 'VTU'}
        checks['field_sha256'] = sha256(field_path) == str(
            artifact['sha256']).lower()
        checks['mesh_id'] = bool(artifact.get('mesh_id'))
        checks['coordinate_system'] = bool(artifact.get('coordinate_system'))
    except (KeyError, TypeError, ValueError):
        checks['field_artifact'] = False
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError('invalid hydrodynamic handoff: ' + ', '.join(failed))
    return {'sha256': sha256(path), 'field_sha256': artifact['sha256'],
            'mesh_id': artifact['mesh_id'], 'fields': sorted(fields)}


def validate_coupling_state(case: str | Path, *, candidate_id: str,
                            reactor_type: str, temperature_K: float,
                            iteration: int) -> dict:
    """Validate one feedback state observed by the runner's outer loop."""
    root = Path(case).resolve()
    path = root / 'hydrogen_coupling_state.json'
    try:
        value = json.loads(path.read_text())
        residual = float(value['residual_relative'])
        tolerance = float(value['tolerance_relative'])
        feedback = value['feedback_artifact']
        feedback_path = case_file(root, feedback['path'], 'coupling feedback')
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise ValueError(f'coupling state is unreadable or incomplete: {exc}') from exc
    checks = {
        'schema_version': value.get('schema_version') == 1,
        'candidate_id': value.get('candidate_id') == candidate_id,
        'reactor_type': value.get('reactor_type') == reactor_type,
        'temperature_K': abs(float(value.get('temperature_K', -1)) -
                             temperature_K) < 1e-6,
        'iteration': value.get('iteration') == iteration,
        'residual': math.isfinite(residual) and residual >= 0,
        'tolerance': math.isfinite(tolerance) and tolerance >= 0,
        'feedback_sha256': sha256(feedback_path) == str(
            feedback.get('sha256', '')).lower(),
        'converged_consistent': value.get('converged') is (residual <= tolerance),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError('invalid coupling state: ' + ', '.join(failed))
    return {'iteration': iteration, 'residual_relative': residual,
            'tolerance_relative': tolerance,
            'converged': value['converged'],
            'feedback_sha256': feedback['sha256']}


def validate_solver_coupling(case: str | Path, coupling: dict, *,
                             candidate_id: str, temperature_K: float,
                             reactor_type: str) -> dict:
    """Validate Cantera execution and exchanged rates, not a boolean claim."""
    root = Path(case).resolve()
    required = ('schema_version', 'cantera_used', 'mechanism_path',
                'mechanism_sha256', 'cantera_log_path', 'cantera_log_sha256',
                'rate_exchange_path', 'rate_exchange_sha256',
                'coupling_history_path', 'coupling_history_sha256',
                'coupling_method', 'coupling_iterations',
                'coupling_residual_relative', 'coupling_tolerance_relative',
                'exchanged_fields')
    missing = [key for key in required if coupling.get(key) in (None, '', [])]
    if missing:
        raise ValueError(f'solver coupling proof missing: {missing}')
    if coupling['schema_version'] != 1 or coupling['cantera_used'] is not True:
        raise ValueError('invalid solver coupling identity')
    files = {}
    for stem in ('mechanism', 'cantera_log', 'rate_exchange', 'coupling_history'):
        target = case_file(root, coupling[f'{stem}_path'], stem)
        if target.stat().st_size == 0:
            raise ValueError(f'{stem} is empty')
        if sha256(target) != str(coupling[f'{stem}_sha256']).lower():
            raise ValueError(f'{stem} checksum mismatch')
        files[stem] = target
    try:
        exchange = json.loads(files['rate_exchange'].read_text())
        rates = exchange['reaction_rates_mol_m3_s']
        identity = (exchange.get('schema_version') == 1 and
                    exchange.get('candidate_id') == candidate_id and
                    abs(float(exchange['temperature_K']) - temperature_K) < 1e-6)
        rates_ok = isinstance(rates, dict) and bool(rates) and all(
            math.isfinite(float(value)) for value in rates.values())
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError):
        identity = rates_ok = False
    if not identity or not rates_ok:
        raise ValueError('invalid candidate-specific Cantera rate exchange')
    method = coupling['coupling_method']
    if method not in {'one_way', 'iterative_two_way'}:
        raise ValueError('unsupported coupling method')
    iterations = int(coupling['coupling_iterations'])
    residual = float(coupling['coupling_residual_relative'])
    tolerance = float(coupling['coupling_tolerance_relative'])
    if not all(math.isfinite(x) for x in (residual, tolerance)) or tolerance < 0:
        raise ValueError('invalid coupling residual or tolerance')
    if method == 'iterative_two_way' and (
            iterations < 2 or residual > tolerance):
        raise ValueError('two-way solver coupling is not converged')
    try:
        history = json.loads(files['coupling_history'].read_text())['iterations']
        history_residuals = [float(row['residual_relative']) for row in history]
        history_ids = [int(row['iteration']) for row in history]
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError):
        raise ValueError('coupling history is invalid')
    if (len(history) != iterations or history_ids != list(range(1, iterations + 1))
            or any(not math.isfinite(value) or value < 0
                   for value in history_residuals)
            or not math.isclose(history_residuals[-1], residual,
                                rel_tol=1e-9, abs_tol=1e-12)):
        raise ValueError('coupling history does not match declared convergence')
    required_fields = {'species', 'temperature', 'reaction_heat',
                       'reaction_rates'}
    if reactor_type == 'NTEC':
        required_fields |= {'momentum', 'charge', 'potential'}
    else:
        required_fields |= {'charge', 'potential'}
    if not required_fields.issubset(set(coupling['exchanged_fields'])):
        raise ValueError('solver coupling omits required exchanged fields')
    return {
        'method': method, 'iterations': iterations,
        'residual_relative': residual, 'tolerance_relative': tolerance,
        'converged': method == 'iterative_two_way' and residual <= tolerance,
        'mechanism_sha256': coupling['mechanism_sha256'],
        'cantera_log_sha256': coupling['cantera_log_sha256'],
        'rate_exchange_sha256': coupling['rate_exchange_sha256'],
        'coupling_history_sha256': coupling['coupling_history_sha256'],
        'exchanged_fields': sorted(set(coupling['exchanged_fields'])),
    }
