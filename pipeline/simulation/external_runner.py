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

import pipeline.simulation.result_contract as _result_contract
from pipeline.simulation.result_contract import (
    EXTERNAL_SOLVERS, REQUIRED_OUTPUTS, SCHEMA_VERSION, artifact_path,
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
    moose_environment,
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
    environment = (
        moose_environment(resolve_moose_executable(required=False))
        if backend.startswith('moose') else None
    )
    completed = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, timeout=timeout_s,
        env=environment)
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
    # Run by path: the working directory is the case, not the repository.
    script = Path(__file__).resolve().parents[1] / 'electrochemistry' / 'cantera_coupling.py'
    return [sys.executable, str(script), '--case-dir', str(case),
            '--reactor-type', reactor_type]


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


SCREENING_TIER = 'uncalibrated_multiphysics_screening'
SCREENING_RESULT = 'hydrogen_screening_result.json'
SCREENING_REACTORS = frozenset({'NTEC', 'Electrochemical'})


def _execute_coupled_case(*, mode: str, reactor_type: str, candidate_id: str,
                          temperature_K: float, case_dir: str | Path,
                          moose_input: str | Path | None, timeout_s: int,
                          max_coupling_iterations: int,
                          execution: SolverExecutionServices,
                          require_calibration: bool) -> dict:
    """Run the solver loop shared by production artifacts and screening.

    Performs routing/preflight checks, loads the case, hashes the pristine
    inputs, runs the (MOOSE <-> Cantera iterative or OpenFOAM) backends and
    applies the numeric and physical output checks. Holdout scoring, coupling
    proof validation and persistence stay with the callers.
    """
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
        reactor_type=reactor_type, temperature_K=temperature_K,
        require_calibration=require_calibration)
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
    return {
        'case': case, 'physical_case': physical_case, 'required': required,
        'input_file': input_file, 'input_digest': input_digest,
        'versions': versions, 'handoff': handoff,
        'observed_coupling': observed_coupling, 'outputs': outputs,
        'convergence': convergence, 'metadata': metadata,
    }


def _validate_cantera_coupling(run: dict, *, candidate_id: str,
                               temperature_K: float, reactor_type: str,
                               execution: SolverExecutionServices) -> None:
    """Check the declared coupling proof against runner-observed iterations."""
    case, metadata = run['case'], run['metadata']
    observed_coupling = run['observed_coupling']
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
    run['versions']['cantera'] = execution.cantera_version()


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
    run = _execute_coupled_case(
        mode=mode, reactor_type=reactor_type, candidate_id=candidate_id,
        temperature_K=temperature_K, case_dir=case_dir,
        moose_input=moose_input, timeout_s=timeout_s,
        max_coupling_iterations=max_coupling_iterations,
        execution=execution, require_calibration=True)
    case, physical_case = run['case'], run['physical_case']
    metadata, versions = run['metadata'], run['versions']
    input_file = run['input_file']
    metadata['model_validation'] = score_holdout(
        case / 'hydrogen_validation_records.json',
        physical_case['calibration'])
    if 'cantera' in run['required']:
        _validate_cantera_coupling(
            run, candidate_id=candidate_id, temperature_K=temperature_K,
            reactor_type=reactor_type, execution=execution)
    artifact: ExternalSolverArtifact = {
        'schema_version': SCHEMA_VERSION,
        'candidate_id': candidate_id,
        'pathway_mode': mode,
        'reactor_type': reactor_type,
        'temperature_K': float(temperature_K),
        'complete': True,
        'backend_solvers': versions,
        'convergence': run['convergence'],
        'outputs': run['outputs'],
        'provenance': {
            'input_sha256': run['input_digest'],
            'model_source': model_source,
            'moose_input_sha256': sha256(input_file) if input_file else None,
            'hydrodynamic_handoff': run['handoff'],
            'observed_coupling_iterations': run['observed_coupling'],
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
        path_builder=artifact_path,
        validator=_result_contract.load_validated_artifact)


def screening_case_digest(case_dir: str | Path, reactor_type: str) -> str:
    """Return the input digest ``run_screening_case`` records for a pristine case.

    Args:
        case_dir: Pristine screening case directory (inputs only).
        reactor_type: ``NTEC`` or ``Electrochemical``; selects the MOOSE input.

    Returns:
        Hex SHA-256 over the case tree and the repository MOOSE input template.
    """
    case = Path(case_dir).expanduser().resolve()
    return _tree_digest(case, (default_moose_input(reactor_type),))


def run_screening_case(*, mode: str, reactor_type: str, candidate_id: str,
                       temperature_K: float, case_dir: str | Path,
                       timeout_s: int = 3600,
                       max_coupling_iterations: int = 20,
                       execution: SolverExecutionServices | None = None) -> dict:
    """Run an UNCALIBRATED MOOSE <-> Cantera screening case for one candidate.

    Uses the same iterative two-way coupling, numeric (mesh/conservation)
    checks, physical output checks and coupling-proof validation as
    ``run_backend``, but the case carries assumed representative operating
    conditions and no calibration block. Therefore no holdout score is
    computed, no production artifact is persisted under the results tree, and
    the result is labelled ``uncalibrated_multiphysics_screening`` and never
    authorizes candidate exclusion. The result is written to
    ``hydrogen_screening_result.json`` in ``case_dir``.

    Args:
        mode: Pathway mode; must route ``reactor_type``.
        reactor_type: ``NTEC`` or ``Electrochemical`` only.
        candidate_id: Stable candidate identifier recorded in the case.
        temperature_K: Absolute operating temperature in kelvin.
        case_dir: Pristine case directory built by
            ``pipeline.simulation.screening_cases.build_screening_case``.
        timeout_s: Wall-time bound for each solver process in seconds.
        max_coupling_iterations: Maximum MOOSE <-> Cantera iterations (>= 2).
        execution: Injectable solver process boundary; defaults to production.

    Returns:
        Screening record with ``valid`` True, ``tier``, ``outputs``,
        ``convergence``, ``backend_solvers``, ``solver_coupling``,
        ``case_dir``, ``input_sha256`` and ``candidate_exclusion_authorized``
        False.

    Raises:
        ValueError: For a non-screening reactor type or invalid case.
        RuntimeError: On solver, convergence, numeric or physical failure.
    """
    if reactor_type not in SCREENING_REACTORS:
        raise ValueError(f'{reactor_type} has no uncalibrated screening tier')
    execution = execution or default_solver_execution_services()
    run = _execute_coupled_case(
        mode=mode, reactor_type=reactor_type, candidate_id=candidate_id,
        temperature_K=temperature_K, case_dir=case_dir, moose_input=None,
        timeout_s=timeout_s, max_coupling_iterations=max_coupling_iterations,
        execution=execution, require_calibration=False)
    convergence, outputs = run['convergence'], run['outputs']
    numeric_failures = [
        name for name, passed in (
            ('converged', convergence.get('converged') is True),
            ('mesh_independent', convergence.get('mesh_independent') is True),
            ('conservation_residuals',
             convergence.get('conservation_satisfied') is True),
            ('required_outputs', REQUIRED_OUTPUTS[reactor_type].issubset(outputs)),
        ) if not passed]
    if not numeric_failures:
        for name in ('CH4_conversion', 'H2_selectivity', 'solid_C_selectivity',
                     'faradaic_efficiency_H2'):
            if name not in outputs:
                continue
            value = outputs[name]
            if (isinstance(value, bool) or not isinstance(value, (int, float)) or
                    not math.isfinite(float(value)) or not 0 <= float(value) <= 1):
                numeric_failures.append(f'bounded_{name}')
    if numeric_failures:
        raise RuntimeError('screening solver checks failed: ' +
                           ', '.join(numeric_failures))
    _validate_cantera_coupling(
        run, candidate_id=candidate_id, temperature_K=temperature_K,
        reactor_type=reactor_type, execution=execution)
    case, metadata, input_file = run['case'], run['metadata'], run['input_file']
    physical_case = run['physical_case']
    result = {
        'schema_version': 1,
        'valid': True,
        'tier': SCREENING_TIER,
        'calibrated': False,
        'candidate_exclusion_authorized': False,
        'candidate_id': candidate_id,
        'pathway_mode': mode,
        'reactor_type': reactor_type,
        'temperature_K': float(temperature_K),
        'outputs': outputs,
        'convergence': convergence,
        'backend_solvers': run['versions'],
        'solver_coupling': metadata['solver_coupling'],
        'mechanism': metadata.get('mechanism'),
        'case_dir': str(case),
        'input_sha256': run['input_digest'],
        'moose_input_sha256': sha256(input_file) if input_file else None,
        'hydrodynamic_handoff': run['handoff'],
        'observed_coupling_iterations': run['observed_coupling'],
        'kinetics_source': physical_case['kinetics']['source'],
        'feed_source': physical_case['feed']['source'],
        'limitations': ['uncalibrated_representative_operating_point',
                        'no_holdout_validation'],
    }
    if reactor_type == 'Electrochemical':
        result['electrolyte_phase'] = metadata.get(
            'electrolyte_phase', physical_case.get('electrolyte_phase'))
    (case / SCREENING_RESULT).write_text(
        json.dumps(result, indent=1, sort_keys=True) + '\n')
    return result


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
