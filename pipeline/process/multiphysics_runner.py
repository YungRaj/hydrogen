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
import math
import shutil
import subprocess
import sys
from pathlib import Path

from pipeline.common.executables import resolve_executable
from pipeline.process.multiphysics_contract import (
    EXTERNAL_SOLVERS, SCHEMA_VERSION, artifact_path, load_validated_artifact,
    mode_preflight, verify_numerics, verify_physical_outputs)
from pipeline.process.pathway_modes import MODE_CHOICES, reactor_types_for_mode
from pipeline.process.physical_case import case_summary, load_physical_case
from pipeline.process.model_validation import score_holdout
from pipeline.process.coupling_contract import (
    require_pristine_case, sha256, validate_hydrodynamic_handoff,
    validate_coupling_state, validate_solver_coupling)


def _tree_digest(path: Path, extra_files: tuple[Path, ...] = ()) -> str:
    digest = hashlib.sha256()
    for item in sorted(value for value in path.rglob('*') if value.is_file()):
        if item.is_symlink():
            continue
        digest.update(str(item.relative_to(path)).encode())
        digest.update(item.read_bytes())
    for item in sorted(set(file.resolve() for file in extra_files)):
        if not item.is_file() or item.is_relative_to(path.resolve()):
            continue
        digest.update(b'external-input\0')
        digest.update(item.name.encode())
        digest.update(item.read_bytes())
    return digest.hexdigest()


def _read_json(path: Path) -> dict:
    if not path.is_file():
        raise RuntimeError(f'required solver output is missing: {path}')
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise RuntimeError(f'solver output must be a JSON object: {path}')
    return value


def _run(command: list[str], cwd: Path, timeout_s: int,
         backend: str = 'solver') -> None:
    completed = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, timeout=timeout_s)
    (cwd / f'hydrogen_{backend}.stdout.log').write_text(completed.stdout)
    (cwd / f'hydrogen_{backend}.stderr.log').write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError(
            f'solver exited {completed.returncode}; see logs in {cwd}')


def _openfoam_version(executable: str) -> str:
    probe = subprocess.run(
        [executable, '-help'], capture_output=True, text=True, timeout=30)
    combined = probe.stdout + probe.stderr
    for line in combined.splitlines():
        if 'OpenFOAM-v' in line:
            return line.strip()
    return Path(executable).name


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
                timeout_s: int = 86400,
                max_coupling_iterations: int = 20) -> Path:
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
    require_pristine_case(case)
    physical_case = load_physical_case(
        case / 'hydrogen_case.json', candidate_id=candidate_id, mode=mode,
        reactor_type=reactor_type, temperature_K=temperature_K)
    # Hash the pristine model/case, including hydrogen_input.json, before a
    # backend creates time directories, logs, or result JSON.
    script = None
    if 'fenicsx' in required:
        if fenics_model is None:
            raise ValueError('a FEniCSx model script is required for this mode')
        script = Path(fenics_model).expanduser().resolve()
        if not script.is_file():
            raise ValueError(f'FEniCSx model script does not exist: {script}')
    input_digest = _tree_digest(case, (script,) if script else ())
    versions, handoff, observed_coupling = {}, None, []
    executable = preflight['solvers']['openfoam']['executable'] \
        if 'openfoam' in required else None
    if executable:
        versions['openfoam'] = _openfoam_version(executable)
    fenics_command = None
    if 'fenicsx' in required:
        fenics_command, version = _fenics_command(script)
        versions['fenicsx'] = version
    specialized = reactor_type in {'NTEC', 'Electrochemical'}
    iterations = range(1, max_coupling_iterations + 1) if specialized else (1,)
    if specialized and max_coupling_iterations < 2:
        raise ValueError('max_coupling_iterations must be at least 2')
    for iteration in iterations:
        if specialized:
            (case / 'hydrogen_coupling_request.json').write_text(json.dumps({
                'schema_version': 1, 'candidate_id': candidate_id,
                'reactor_type': reactor_type, 'temperature_K': temperature_K,
                'iteration': iteration}, sort_keys=True) + '\n')
        if executable:
            _run([executable, '-case', str(case)], case, timeout_s,
                 f'openfoam.iteration-{iteration}')
            if reactor_type == 'NTEC':
                try:
                    handoff = validate_hydrodynamic_handoff(
                        case, candidate_id=candidate_id, mode=mode,
                        reactor_type=reactor_type, temperature_K=temperature_K,
                        iteration=iteration)
                except ValueError as exc:
                    raise RuntimeError(str(exc)) from exc
        if fenics_command:
            _run(fenics_command, case, timeout_s,
                 f'fenicsx.iteration-{iteration}')
        if specialized:
            try:
                state = validate_coupling_state(
                    case, candidate_id=candidate_id,
                    reactor_type=reactor_type, temperature_K=temperature_K,
                    iteration=iteration)
            except ValueError as exc:
                raise RuntimeError(str(exc)) from exc
            observed_coupling.append(state)
            if state['converged'] and iteration >= 2:
                break
    if specialized and (not observed_coupling or
                        not observed_coupling[-1]['converged']):
        raise RuntimeError('iterative solver coupling did not converge')
    outputs = _read_json(case / 'hydrogen_outputs.json')
    convergence = _read_json(case / 'hydrogen_convergence.json')
    numerical = verify_numerics(convergence, reactor_type)
    convergence.update(numerical)
    physical_failures = verify_physical_outputs(
        reactor_type, outputs, physical_case)
    if physical_failures:
        raise RuntimeError(
            'solver outputs violate physical checks: ' +
            ', '.join(physical_failures))
    metadata = _read_json(case / 'hydrogen_metadata.json') \
        if (case / 'hydrogen_metadata.json').is_file() else {}
    metadata['model_validation'] = score_holdout(
        case / 'hydrogen_validation_records.json',
        physical_case['calibration'])
    if 'cantera' in required:
        try:
            coupling = validate_solver_coupling(
                case, metadata.get('solver_coupling', {}),
                candidate_id=candidate_id, temperature_K=temperature_K,
                reactor_type=reactor_type)
        except ValueError as exc:
            raise RuntimeError(str(exc)) from exc
        if coupling['method'] != 'iterative_two_way' or not coupling['converged']:
            raise RuntimeError(
                'NTEC/electrochemical production artifacts require converged '
                'iterative two-way Cantera coupling')
        if (coupling['iterations'] != len(observed_coupling) or
                not math.isclose(coupling['residual_relative'],
                                 observed_coupling[-1]['residual_relative'],
                                 rel_tol=1e-9, abs_tol=1e-12) or
                not math.isclose(coupling['tolerance_relative'],
                                 observed_coupling[-1]['tolerance_relative'],
                                 rel_tol=1e-9, abs_tol=1e-12)):
            raise RuntimeError('coupling proof does not match runner-observed iterations')
        metadata['solver_coupling'] = coupling
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
            'fenics_model_sha256': sha256(script) if script else None,
            'hydrodynamic_handoff': handoff,
            'observed_coupling_iterations': observed_coupling,
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
    parser.add_argument('--max-coupling-iterations', type=int, default=20)
    args = parser.parse_args()
    print(run_backend(
        mode=args.mode, reactor_type=args.reactor_type,
        candidate_id=args.candidate_id, temperature_K=args.temperature_K,
        case_dir=args.case_dir, results_dir=args.results_dir,
        model_source=args.model_source, fenics_model=args.fenics_model,
        timeout_s=args.timeout_s,
        max_coupling_iterations=args.max_coupling_iterations))


if __name__ == '__main__':
    main()
