"""Fail-closed contracts for mode-owned external multiphysics calculations."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import math
from functools import lru_cache
from pathlib import Path


SCHEMA_VERSION = 1
MODE_SOLVERS = {
    'thermocatalytic': {'cantera', 'openfoam'},
    'thermocatalytic_pfr': {'cantera'},
    'thermocatalytic_fluidized': {'cantera', 'openfoam'},
    'mmbcr': {'cantera', 'openfoam'},
    'ntec': {'cantera', 'openfoam', 'fenicsx'},
    'electrochemical': {'cantera', 'fenicsx'},
}
EXTERNAL_SOLVERS = {
    'Fluidized': {'openfoam'},
    'MMBCR': {'openfoam'},
    'NTEC': {'openfoam', 'fenicsx', 'cantera'},
    'Electrochemical': {'fenicsx', 'cantera'},
}

REQUIRED_OUTPUTS = {
    'Fluidized': {'gas_velocity_m_s', 'u_mf_m_s', 'bubble_fraction'},
    'MMBCR': {'gas_velocity_m_s', 'gas_holdup_fraction',
              'bubble_diameter_mm'},
    'NTEC': {'CH4_conversion', 'H2_selectivity',
             'solid_C_selectivity', 'specific_energy_kWh_kg_H2'},
    'Electrochemical': {'CH4_conversion', 'H2_selectivity',
                        'faradaic_efficiency_H2', 'current_density_A_cm2',
                        'cell_voltage_V', 'electrical_power_density_W_cm2'},
}
REQUIRED_BALANCES = {
    'Fluidized': {'mass', 'carbon', 'hydrogen', 'energy'},
    'MMBCR': {'mass', 'carbon', 'hydrogen', 'energy'},
    'NTEC': {'mass', 'carbon', 'hydrogen', 'energy', 'charge'},
    'Electrochemical': {'mass', 'carbon', 'hydrogen', 'charge'},
}


def verify_numerics(convergence: dict, reactor_type: str | None = None) -> dict:
    """Independently compute mesh stability and conservation residuals."""
    series = convergence.get('mesh_series', [])
    tolerance = convergence.get('mesh_tolerance_relative')
    mesh_ok = False
    mesh_change = math.inf
    try:
        cells = [int(row['cells']) for row in series]
        values = [float(row['observable']) for row in series]
        tolerance = float(tolerance)
        if (len(series) >= 3 and all(b > a > 0 for a, b in zip(cells, cells[1:]))
                and all(math.isfinite(value) for value in values)
                and math.isfinite(tolerance) and tolerance >= 0):
            mesh_change = abs(values[-1] - values[-2]) / max(
                abs(values[-1]), 1e-12)
            mesh_ok = mesh_change <= tolerance
    except (KeyError, TypeError, ValueError):
        pass
    budgets = convergence.get('conservation_budgets', {})
    residuals = {}
    if isinstance(budgets, dict):
        for name, budget in budgets.items():
            try:
                inlet, outlet = float(budget['inlet']), float(budget['outlet'])
                if math.isfinite(inlet) and math.isfinite(outlet):
                    residuals[name] = abs(outlet - inlet) / max(abs(inlet), 1e-12)
            except (KeyError, TypeError, ValueError):
                continue
    required_balances = REQUIRED_BALANCES.get(reactor_type, set())
    conservation_ok = (bool(residuals) and len(residuals) == len(budgets) and
                       required_balances.issubset(residuals) and all(
                           value <= 1e-5 for value in residuals.values()))
    return {
        'mesh_independent': mesh_ok,
        'mesh_relative_change': mesh_change,
        'conservation_relative_residuals': residuals,
        'conservation_satisfied': conservation_ok,
    }


def verify_physical_outputs(reactor_type: str, outputs: dict,
                            physical_case: dict | None = None) -> list[str]:
    """Return violated mode-specific identities and physical bounds."""
    failed = []
    try:
        if reactor_type == 'Fluidized':
            if not float(outputs['gas_velocity_m_s']) > float(outputs['u_mf_m_s']) > 0:
                failed.append('fluidization_velocity')
            if not 0 < float(outputs['bubble_fraction']) < 1:
                failed.append('bubble_fraction')
        elif reactor_type == 'MMBCR':
            if float(outputs['gas_velocity_m_s']) <= 0:
                failed.append('gas_velocity')
            if not 0 < float(outputs['gas_holdup_fraction']) < 1:
                failed.append('gas_holdup_fraction')
            if float(outputs['bubble_diameter_mm']) <= 0:
                failed.append('bubble_diameter')
            if physical_case and float(outputs['bubble_diameter_mm']) / 1000 >= float(
                    physical_case['geometry']['column_diameter_m']):
                failed.append('bubble_smaller_than_column')
        elif reactor_type == 'NTEC':
            if float(outputs['specific_energy_kWh_kg_H2']) <= 0:
                failed.append('specific_energy')
        elif reactor_type == 'Electrochemical':
            current = float(outputs['current_density_A_cm2'])
            voltage = float(outputs['cell_voltage_V'])
            if current < 0 or voltage <= 0:
                failed.append('electrochemical_current_voltage')
            expected = current * voltage
            actual = float(outputs['electrical_power_density_W_cm2'])
            if actual < 0:
                failed.append('electrical_power_nonnegative')
            if abs(actual - expected) > 1e-6 * max(abs(expected), 1.0):
                failed.append('electrical_power_identity')
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        failed.append('physical_output_types')
    return failed


def _conda_module_available(environment: str, module: str) -> bool:
    conda = shutil.which('conda')
    if not conda:
        return False
    probe = subprocess.run(
        [conda, 'run', '-n', environment, 'python', '-c',
         f'import {module}'], capture_output=True, timeout=30)
    return probe.returncode == 0


@lru_cache(maxsize=1)
def solver_preflight() -> dict:
    """Discover supported solvers without assuming installation paths."""
    from pipeline.common.executables import resolve_executable

    cantera = (importlib.util.find_spec('cantera') is not None or
               _conda_module_available('cp2k-env', 'cantera') or
               _conda_module_available('fenicsx-env', 'cantera'))
    fenicsx = importlib.util.find_spec('dolfinx') is not None or \
        _conda_module_available('fenicsx-env', 'dolfinx')
    openfoam = resolve_executable(
        'multiphaseEulerFoam', env_var='OPENFOAM_SOLVER',
        conda_env='openfoam-env', required=False)
    return {
        'cantera': {'available': cantera},
        'openfoam': {'available': openfoam is not None, 'executable': openfoam},
        'fenicsx': {'available': fenicsx},
    }


def mode_preflight(mode: str) -> dict:
    available = solver_preflight()
    required = sorted(MODE_SOLVERS[mode])
    missing = [name for name in required if not available[name]['available']]
    return {'mode': mode, 'required': required, 'missing': missing,
            'ready': not missing, 'solvers': available}


def artifact_path(root: str | Path, candidate_id: str, pathway_mode: str,
                  reactor_type: str, temperature_K: float) -> Path:
    safe = Path(candidate_id).name
    if safe != candidate_id:
        raise ValueError('candidate_id is not a safe path component')
    return (Path(root) / safe / pathway_mode /
            f'{reactor_type}_{int(temperature_K)}K.json')


def load_validated_artifact(root: str | Path | None, candidate_id: str,
                            pathway_mode: str, reactor_type: str,
                            temperature_K: float) -> dict:
    """Load solver evidence and reject incomplete, mismatched, or unconverged data."""
    if not root:
        return {'valid': False, 'reason': 'multiphysics_results_dir_not_configured'}
    path = artifact_path(root, candidate_id, pathway_mode, reactor_type,
                         temperature_K)
    if not path.is_file():
        return {'valid': False, 'reason': 'multiphysics_artifact_missing',
                'path': str(path)}
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return {'valid': False, 'reason': 'multiphysics_artifact_unreadable',
                'error': str(exc), 'path': str(path)}
    try:
        temperature_matches = abs(float(value.get('temperature_K', -1)) -
                                  float(temperature_K)) < 1e-6
    except (TypeError, ValueError):
        temperature_matches = False
    checks = {
        'schema_version': value.get('schema_version') == SCHEMA_VERSION,
        'candidate_id': str(value.get('candidate_id')) == candidate_id,
        'pathway_mode': value.get('pathway_mode') == pathway_mode,
        'reactor_type': value.get('reactor_type') == reactor_type,
        'temperature_K': temperature_matches,
        'complete': value.get('complete') is True,
        'converged': value.get('convergence', {}).get('converged') is True,
    }
    declared = set(value.get('backend_solvers', {}))
    required = EXTERNAL_SOLVERS.get(reactor_type, set())
    checks['required_solvers'] = required.issubset(declared)
    checks['outputs'] = isinstance(value.get('outputs'), dict)
    outputs = value.get('outputs', {})
    checks['required_outputs'] = REQUIRED_OUTPUTS[reactor_type].issubset(outputs)
    convergence = value.get('convergence', {})
    numerical = verify_numerics(convergence, reactor_type)
    checks['mesh_independent'] = numerical['mesh_independent']
    checks['conservation_residuals'] = numerical['conservation_satisfied']
    checks['backend_versions'] = all(
        isinstance(value.get('backend_solvers', {}).get(name), str) and
        bool(value['backend_solvers'][name].strip()) for name in required)
    provenance = value.get('provenance', {})
    digest = provenance.get('input_sha256', '')
    checks['input_digest'] = (isinstance(digest, str) and len(digest) == 64 and
                              all(char in '0123456789abcdef' for char in digest.lower()))
    checks['model_source'] = bool(provenance.get('model_source'))
    if 'fenicsx' in required:
        script_digest = provenance.get('fenics_model_sha256', '')
        checks['fenics_model_digest'] = (
            isinstance(script_digest, str) and len(script_digest) == 64 and
            all(char in '0123456789abcdef' for char in script_digest.lower()))
    physical_case = value.get('physical_case', {})
    try:
        calibration_count = int(physical_case.get('calibration_count', 0))
        holdout_count = int(physical_case.get('holdout_validation_count', 0))
    except (TypeError, ValueError):
        calibration_count = holdout_count = 0
    checks['physical_case'] = (
        physical_case.get('schema_version') == 1 and
        physical_case.get('disjoint_holdout') is True and
        calibration_count > 0 and holdout_count > 0 and
        bool(physical_case.get('kinetics_source')) and
        bool(physical_case.get('feed_source')) and
        bool(physical_case.get('calibration_source')))
    validation = value.get('model_validation', {})
    try:
        holdout_error = float(validation.get('holdout_error'))
        acceptance_threshold = float(validation.get('acceptance_threshold'))
    except (TypeError, ValueError):
        holdout_error = acceptance_threshold = math.nan
    checks['holdout_validation'] = (
        validation.get('metric') in {'rmse', 'mae', 'relative_rmse'} and
        math.isfinite(holdout_error) and holdout_error >= 0 and
        math.isfinite(acceptance_threshold) and acceptance_threshold >= 0 and
        holdout_error <= acceptance_threshold and
        validation.get('passed') is True and
        bool(validation.get('record_source')))
    if reactor_type == 'NTEC':
        calibration = value.get('calibration', {})
        checks['paired_control_calibration'] = (
            calibration.get('paired_control') is True and
            bool(calibration.get('paired_control_source')))
        handoff = provenance.get('hydrodynamic_handoff', {})
        checks['hydrodynamic_handoff'] = (
            isinstance(handoff.get('sha256'), str) and
            len(handoff['sha256']) == 64 and bool(handoff.get('mesh_id')) and
            isinstance(handoff.get('field_sha256'), str) and
            len(handoff['field_sha256']) == 64)
    if reactor_type == 'Electrochemical':
        mechanism = value.get('mechanism', {})
        checks['electrochemical_mechanism'] = (
            mechanism.get('complete') is True and
            bool(mechanism.get('source')))
        checks['electrolyte_phase'] = value.get(
            'electrolyte_phase') in {'aqueous', 'molten'}
    if reactor_type in {'NTEC', 'Electrochemical'}:
        coupling = value.get('solver_coupling', {})
        try:
            iterations = int(coupling.get('iterations', 0))
            coupling_residual = float(coupling.get(
                'residual_relative', math.inf))
            coupling_tolerance = float(coupling.get(
                'tolerance_relative', -1))
        except (TypeError, ValueError):
            iterations = 0
            coupling_residual, coupling_tolerance = math.inf, -1
        checks['iterative_solver_coupling'] = (
            coupling.get('method') == 'iterative_two_way' and
            coupling.get('converged') is True and
            iterations >= 2 and coupling_residual <= coupling_tolerance and
            all(isinstance(coupling.get(name), str) and
                len(coupling[name]) == 64 for name in (
                    'mechanism_sha256', 'cantera_log_sha256',
                    'rate_exchange_sha256', 'coupling_history_sha256')))
        observed = provenance.get('observed_coupling_iterations', [])
        try:
            observed_ok = (
                isinstance(observed, list) and len(observed) == iterations and
                iterations >= 2 and observed[-1].get('converged') is True and
                int(observed[-1].get('iteration', 0)) == iterations and
                math.isclose(float(observed[-1].get(
                    'residual_relative', math.inf)), coupling_residual,
                    rel_tol=1e-9, abs_tol=1e-12))
        except (AttributeError, TypeError, ValueError):
            observed_ok = False
        checks['runner_observed_coupling'] = observed_ok
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        return {'valid': False, 'reason': 'multiphysics_artifact_invalid',
                'failed_checks': failed, 'path': str(path)}
    bounds = {
        name: outputs.get(name) for name in (
            'CH4_conversion', 'H2_selectivity', 'solid_C_selectivity',
            'faradaic_efficiency_H2') if name in outputs}
    if any(not isinstance(number, (int, float)) or
           not math.isfinite(float(number)) or not 0 <= float(number) <= 1
           for number in bounds.values()):
        return {'valid': False, 'reason': 'multiphysics_artifact_invalid',
                'failed_checks': ['bounded_outputs'], 'path': str(path)}
    physical_failures = verify_physical_outputs(reactor_type, outputs)
    if physical_failures:
        return {'valid': False, 'reason': 'multiphysics_artifact_invalid',
                'failed_checks': physical_failures, 'path': str(path)}
    return {'valid': True, 'artifact': value, 'path': str(path)}
