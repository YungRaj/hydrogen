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
        'mesh_independent': value.get('convergence', {}).get(
            'mesh_independent') is True,
    }
    declared = set(value.get('backend_solvers', {}))
    required = EXTERNAL_SOLVERS.get(reactor_type, set())
    checks['required_solvers'] = required.issubset(declared)
    checks['outputs'] = isinstance(value.get('outputs'), dict)
    outputs = value.get('outputs', {})
    checks['required_outputs'] = REQUIRED_OUTPUTS[reactor_type].issubset(outputs)
    convergence = value.get('convergence', {})
    residuals = convergence.get('conservation_relative_residuals', {})
    checks['conservation_residuals'] = bool(residuals) and all(
        isinstance(number, (int, float)) and math.isfinite(float(number)) and
        abs(float(number)) <= 1e-5 for number in residuals.values())
    checks['backend_versions'] = all(
        isinstance(value.get('backend_solvers', {}).get(name), str) and
        bool(value['backend_solvers'][name].strip()) for name in required)
    provenance = value.get('provenance', {})
    digest = provenance.get('input_sha256', '')
    checks['input_digest'] = (isinstance(digest, str) and len(digest) == 64 and
                              all(char in '0123456789abcdef' for char in digest.lower()))
    checks['model_source'] = bool(provenance.get('model_source'))
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
        validation.get('passed') is True)
    if reactor_type == 'NTEC':
        calibration = value.get('calibration', {})
        checks['paired_control_calibration'] = (
            calibration.get('paired_control') is True and
            bool(calibration.get('paired_control_source')))
    if reactor_type == 'Electrochemical':
        mechanism = value.get('mechanism', {})
        checks['electrochemical_mechanism'] = (
            mechanism.get('complete') is True and
            bool(mechanism.get('source')))
        checks['electrolyte_phase'] = value.get(
            'electrolyte_phase') in {'aqueous', 'molten'}
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
    return {'valid': True, 'artifact': value, 'path': str(path)}
