"""Portable execution of required OpenFOAM and MOOSE pathway backends.

OpenFOAM cases and repository MOOSE models own their discretization and write
``hydrogen_outputs.json`` plus ``hydrogen_convergence.json``. This runner executes them without shell
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

from pipeline.simulation.result_contract import (
    EXTERNAL_SOLVERS, SCHEMA_VERSION, artifact_path, load_validated_artifact,
    mode_preflight, verify_numerics, verify_physical_outputs)
from pipeline.reactors.modes import MODE_CHOICES, reactor_types_for_mode
from pipeline.simulation.physical_case import (
    case_summary, load_physical_case, surrogate_inputs)
from pipeline.simulation.model_validation import score_holdout
from pipeline.simulation.solver_handoff import (
    require_pristine_case, sha256, validate_hydrodynamic_handoff,
    validate_coupling_state, validate_solver_coupling)
from pipeline.simulation.solver_execution import SolverExecutionServices
from pipeline.simulation.artifact_store import persist_validated_artifact
from pipeline.simulation.moose_backend import (
    build_moose_command, default_moose_input, moose_case_parameters,
    moose_version, resolve_moose_executable)
from pipeline.data_models.artifacts import ExternalSolverArtifact


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


def _moose_command(input_file: Path, case: Path,
                   parameters: dict) -> tuple[list[str], str]:
    executable = resolve_moose_executable(required=True)
    return (build_moose_command(executable, input_file, case, parameters),
            moose_version(executable))


def _coupling_command(case: Path, reactor_type: str) -> list[str]:
    """Run the repository-owned Cantera interface step for one iteration."""
    return [sys.executable, '-m', 'pipeline.electrochemistry.cantera_coupling',
            '--case-dir', str(case), '--reactor-type', reactor_type]


def _cantera_version() -> str:
    if importlib.util.find_spec('cantera') is not None:
        import cantera
        return str(cantera.__version__)
    conda = shutil.which('conda')
    if conda:
        for environment in ('cp2k-env', 'moose-env'):
            probe = subprocess.run(
                [conda, 'run', '-n', environment, 'python', '-c',
                 'import cantera; print(cantera.__version__)'],
                capture_output=True, text=True, timeout=30)
            if probe.returncode == 0 and probe.stdout.strip():
                return probe.stdout.strip().splitlines()[-1]
    raise RuntimeError('Cantera is required in the current, cp2k-env, or moose-env environment')


def default_solver_execution_services() -> SolverExecutionServices:
    """Bind production solver operations at call time for patchability/tests.

    Returns:
        Computed `SolverExecutionServices` result.
    """
    return SolverExecutionServices(
        preflight=mode_preflight, execute=_run,
        openfoam_version=_openfoam_version,
        moose_command=_moose_command, coupling_command=_coupling_command,
        cantera_version=_cantera_version)


def run_backend(*, mode: str, reactor_type: str, candidate_id: str,
                temperature_K: float, case_dir: str | Path,
                results_dir: str | Path, model_source: str,
                moose_input: str | Path | None = None,
                timeout_s: int = 86400,
                max_coupling_iterations: int = 20,
                execution: SolverExecutionServices | None = None) -> Path:
    """Execute the external case backends and emit one validated artifact.

        Thermal Fluidized/MMBCR artifacts supply OpenFOAM hydrodynamics which the
        reactor layer subsequently couples to Cantera kinetics. NTEC and
        electrochemical cases run MOOSE for the spatial fields and the
        repository's Cantera interface step for kinetics, iterated two-way
        under this runner's control.

    Args:
        mode: Configured methane-conversion pathway.
        reactor_type: Physical reactor implementation identifier.
        candidate_id: Stable candidate identifier.
        temperature_K: Absolute temperature in kelvin.
        case_dir: Filesystem location used for case dir.
        results_dir: Directory containing or receiving calculation results.
        model_source: Model source used by this operation.
        moose_input: MOOSE input file; defaults to the repository template.
        timeout_s: Maximum allowed wall time in seconds.
        max_coupling_iterations: Bound controlling max coupling iterations.
        execution: Execution used by this operation.

    Returns:
        Filesystem path produced or resolved by the operation.
    """
    execution = execution or default_solver_execution_services()
    if reactor_type not in reactor_types_for_mode(mode):
        raise ValueError(f'{reactor_type} is not assigned to mode {mode}')
    required = EXTERNAL_SOLVERS.get(reactor_type, set())
    if not required:
        raise ValueError(f'{reactor_type} does not require an external backend')
    preflight = execution.preflight(mode)
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
    input_file = None
    if 'moose' in required:
        input_file = Path(moose_input).expanduser().resolve() \
            if moose_input is not None else default_moose_input(reactor_type)
        if not input_file.is_file():
            raise ValueError(f'MOOSE input file does not exist: {input_file}')
    input_digest = _tree_digest(case, (input_file,) if input_file else ())
    versions, handoff, observed_coupling = {}, None, []
    executable = preflight['solvers']['openfoam']['executable'] \
        if 'openfoam' in required else None
    if executable:
        versions['openfoam'] = execution.openfoam_version(executable)
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
            execution.execute(
                [executable, '-case', str(case)], case, timeout_s,
                f'openfoam.iteration-{iteration}')
        if input_file is not None:
            command, versions['moose'] = execution.moose_command(
                input_file, case,
                moose_case_parameters(physical_case, case, iteration))
            execution.execute(command, case, timeout_s,
                              f'moose.iteration-{iteration}')
        if specialized:
            execution.execute(
                execution.coupling_command(case, reactor_type), case,
                timeout_s, f'cantera.iteration-{iteration}')
            if reactor_type == 'NTEC':
                # The Cantera step packages MOOSE's flow field as the handoff.
                try:
                    handoff = validate_hydrodynamic_handoff(
                        case, candidate_id=candidate_id, mode=mode,
                        reactor_type=reactor_type, temperature_K=temperature_K,
                        iteration=iteration)
                except ValueError as exc:
                    raise RuntimeError(str(exc)) from exc
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
        versions['cantera'] = execution.cantera_version()
    artifact: ExternalSolverArtifact = {
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
            'moose_input_sha256': sha256(input_file) if input_file else None,
            'hydrodynamic_handoff': handoff,
            'observed_coupling_iterations': observed_coupling,
        },
        'physical_case': case_summary(physical_case),
        # A compact numeric snapshot permits independently trained transport
        # surrogates without divorcing a row from this validated artifact.
        'surrogate_inputs': surrogate_inputs(physical_case),
    }
    for key in ('calibration', 'mechanism', 'electrolyte_phase',
                'solver_coupling', 'model_validation'):
        if key in metadata:
            artifact[key] = metadata[key]
    return persist_validated_artifact(
        artifact, results_dir=results_dir, candidate_id=candidate_id,
        mode=mode, reactor_type=reactor_type, temperature_K=temperature_K,
        path_builder=artifact_path, validator=load_validated_artifact)


def main() -> None:
    """Run the module command-line entry point.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', required=True, choices=MODE_CHOICES)
    parser.add_argument('--reactor-type', required=True)
    parser.add_argument('--candidate-id', required=True)
    parser.add_argument('--temperature-K', required=True, type=float)
    parser.add_argument('--case-dir', required=True)
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--model-source', required=True)
    parser.add_argument('--moose-input')
    parser.add_argument('--timeout-s', type=int, default=86400)
    parser.add_argument('--max-coupling-iterations', type=int, default=20)
    args = parser.parse_args()
    print(run_backend(
        mode=args.mode, reactor_type=args.reactor_type,
        candidate_id=args.candidate_id, temperature_K=args.temperature_K,
        case_dir=args.case_dir, results_dir=args.results_dir,
        model_source=args.model_source, moose_input=args.moose_input,
        timeout_s=args.timeout_s,
        max_coupling_iterations=args.max_coupling_iterations))


if __name__ == '__main__':
    main()
