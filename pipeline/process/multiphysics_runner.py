"""Portable execution of required OpenFOAM and FEniCSx pathway backends.

Solver cases/models own their discretization and write ``hydrogen_outputs.json``
plus ``hydrogen_convergence.json``. This runner executes them without shell
interpolation and packages their outputs into the pipeline's strict artifact.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

from pipeline.common.executables import resolve_executable
from pipeline.process.multiphysics_contract import (
    EXTERNAL_SOLVERS, SCHEMA_VERSION, artifact_path, load_validated_artifact,
    mode_preflight)
from pipeline.process.pathway_modes import MODE_CHOICES, reactor_types_for_mode
from pipeline.process.physical_case import case_summary, load_physical_case


def _tree_digest(path: Path) -> str:
    digest = hashlib.sha256()
    for item in sorted(value for value in path.rglob('*') if value.is_file()):
        if item.is_symlink():
            continue
        digest.update(str(item.relative_to(path)).encode())
        digest.update(item.read_bytes())
    return digest.hexdigest()


def _read_json(path: Path) -> dict:
    if not path.is_file():
        raise RuntimeError(f'required solver output is missing: {path}')
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise RuntimeError(f'solver output must be a JSON object: {path}')
    return value


def _run(command: list[str], cwd: Path, timeout_s: int) -> None:
    completed = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, timeout=timeout_s)
    (cwd / 'hydrogen_solver.stdout.log').write_text(completed.stdout)
    (cwd / 'hydrogen_solver.stderr.log').write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError(
            f'solver exited {completed.returncode}; see logs in {cwd}')


def _fenics_command(model_script: Path) -> tuple[list[str], str]:
    if importlib.util.find_spec('dolfinx') is not None:
        import dolfinx
        return [sys.executable, str(model_script)], str(dolfinx.__version__)
    conda = shutil.which('conda')
    if not conda:
        raise RuntimeError('FEniCSx is required; install dolfinx or fenicsx-env')
    probe = subprocess.run(
        [conda, 'run', '-n', 'fenicsx-env', 'python', '-c',
         'import dolfinx; print(dolfinx.__version__)'],
        capture_output=True, text=True, timeout=30)
    if probe.returncode:
        raise RuntimeError('FEniCSx is required in the documented fenicsx-env')
    version = probe.stdout.strip().splitlines()[-1]
    return [conda, 'run', '-n', 'fenicsx-env', 'python', str(model_script)], version


def _cantera_version() -> str:
    if importlib.util.find_spec('cantera') is not None:
        import cantera
        return str(cantera.__version__)
    conda = shutil.which('conda')
    if conda:
        for environment in ('cp2k-env', 'fenicsx-env'):
            probe = subprocess.run(
                [conda, 'run', '-n', environment, 'python', '-c',
                 'import cantera; print(cantera.__version__)'],
                capture_output=True, text=True, timeout=30)
            if probe.returncode == 0 and probe.stdout.strip():
                return probe.stdout.strip().splitlines()[-1]
    raise RuntimeError('Cantera is required in the current, cp2k-env, or fenicsx-env environment')


def run_backend(*, mode: str, reactor_type: str, candidate_id: str,
                temperature_K: float, case_dir: str | Path,
                results_dir: str | Path, model_source: str,
                fenics_model: str | Path | None = None,
                timeout_s: int = 86400) -> Path:
    """Execute the external case backends and emit one validated artifact.

    Thermal Fluidized/MMBCR artifacts supply OpenFOAM hydrodynamics which the
    reactor layer subsequently couples to Cantera kinetics. NTEC and
    electrochemical cases must perform and declare their Cantera coupling in
    the external model because those modes do not use the thermal reactor.
    """
    if reactor_type not in reactor_types_for_mode(mode):
        raise ValueError(f'{reactor_type} is not assigned to mode {mode}')
    required = EXTERNAL_SOLVERS.get(reactor_type, set())
    if not required:
        raise ValueError(f'{reactor_type} does not require an external backend')
    preflight = mode_preflight(mode)
    missing = sorted(required.intersection(preflight['missing']))
    if missing:
        raise RuntimeError(f'required solver(s) unavailable: {missing}')

    case = Path(case_dir).expanduser().resolve()
    if not case.is_dir():
        raise ValueError(f'case directory does not exist: {case}')
    physical_case = load_physical_case(
        case / 'hydrogen_case.json', candidate_id=candidate_id, mode=mode,
        reactor_type=reactor_type, temperature_K=temperature_K)
    # Hash the pristine model/case, including hydrogen_input.json, before a
    # backend creates time directories, logs, or result JSON.
    input_digest = _tree_digest(case)
    versions = {}
    if 'openfoam' in required:
        executable = preflight['solvers']['openfoam']['executable']
        _run([executable, '-case', str(case)], case, timeout_s)
        versions['openfoam'] = Path(executable).name
    if 'fenicsx' in required:
        if fenics_model is None:
            raise ValueError('a FEniCSx model script is required for this mode')
        script = Path(fenics_model).expanduser().resolve()
        if not script.is_file():
            raise ValueError(f'FEniCSx model script does not exist: {script}')
        command, version = _fenics_command(script)
        _run(command, case, timeout_s)
        versions['fenicsx'] = version
    outputs = _read_json(case / 'hydrogen_outputs.json')
    convergence = _read_json(case / 'hydrogen_convergence.json')
    metadata = _read_json(case / 'hydrogen_metadata.json') \
        if (case / 'hydrogen_metadata.json').is_file() else {}
    if 'cantera' in required:
        if metadata.get('solver_coupling', {}).get('cantera_used') is not True:
            raise RuntimeError(
                'backend must declare solver_coupling.cantera_used=true')
        versions['cantera'] = _cantera_version()
    artifact = {
        'schema_version': SCHEMA_VERSION,
        'candidate_id': candidate_id,
        'pathway_mode': mode,
        'reactor_type': reactor_type,
        'temperature_K': float(temperature_K),
        'complete': True,
        'backend_solvers': versions,
        'convergence': convergence,
        'outputs': outputs,
        'provenance': {
            'input_sha256': input_digest,
            'model_source': model_source,
        },
        'physical_case': case_summary(physical_case),
    }
    for key in ('calibration', 'mechanism', 'electrolyte_phase',
                'solver_coupling', 'model_validation'):
        if key in metadata:
            artifact[key] = metadata[key]
    target = artifact_path(
        results_dir, candidate_id, mode, reactor_type, temperature_K)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(artifact, indent=2, sort_keys=True) + '\n')
    validated = load_validated_artifact(
        results_dir, candidate_id, mode, reactor_type, temperature_K)
    if not validated['valid']:
        target.unlink(missing_ok=True)
        failed = ', '.join(validated.get('failed_checks', ()))
        raise RuntimeError(
            'solver output did not satisfy the artifact contract: ' +
            (failed or validated['reason']))
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', required=True, choices=MODE_CHOICES)
    parser.add_argument('--reactor-type', required=True)
    parser.add_argument('--candidate-id', required=True)
    parser.add_argument('--temperature-K', required=True, type=float)
    parser.add_argument('--case-dir', required=True)
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--model-source', required=True)
    parser.add_argument('--fenics-model')
    parser.add_argument('--timeout-s', type=int, default=86400)
    args = parser.parse_args()
    print(run_backend(
        mode=args.mode, reactor_type=args.reactor_type,
        candidate_id=args.candidate_id, temperature_K=args.temperature_K,
        case_dir=args.case_dir, results_dir=args.results_dir,
        model_source=args.model_source, fenics_model=args.fenics_model,
        timeout_s=args.timeout_s))


if __name__ == '__main__':
    main()
