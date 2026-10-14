"""MOOSE backend discovery and command construction for NTEC/electrochemistry.

MOOSE solves the spatial fields (potential, species, and for NTEC the liquid
flow) from repository-owned input templates in ``moose/inputs``. Case values
reach MOOSE as command-line overrides of top-level HIT variables, so the
input file itself stays byte-identical and its SHA-256 identifies the model.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path
from typing import Mapping

from pipeline.simulation.executables import resolve_executable

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
MOOSE_INPUT_DIR = REPOSITORY_ROOT / 'moose' / 'inputs'
DEFAULT_MOOSE_INPUTS = {
    'NTEC': MOOSE_INPUT_DIR / 'ntec_shear_layer.i',
    'Electrochemical': MOOSE_INPUT_DIR / 'electrochemical_cell.i',
}
# Executable names shipped by MOOSE distributions, in preference order.
MOOSE_EXECUTABLE_NAMES = ('moose-opt', 'combined-opt', 'moose')
MOOSE_CONDA_ENV = 'moose-env'
PARAMETER_RECORD = 'hydrogen_moose_parameters.json'


def resolve_moose_executable(required: bool = True) -> str | None:
    """Locate a MOOSE modules executable.

    Resolution order: ``MOOSE_EXECUTABLE``, then each known executable name on
    ``PATH``, then the documented ``moose-env`` Conda environment.

    Args:
        required: Raise instead of returning ``None`` when nothing is found.

    Returns:
        Absolute executable path, or ``None`` when optional and unavailable.
    """
    found = resolve_executable(
        MOOSE_EXECUTABLE_NAMES[0], env_var='MOOSE_EXECUTABLE',
        conda_env=MOOSE_CONDA_ENV, required=False)
    if found:
        return found
    for name in MOOSE_EXECUTABLE_NAMES[1:]:
        found = resolve_executable(name, conda_env=MOOSE_CONDA_ENV, required=False)
        if found:
            return found
    if required:
        raise RuntimeError(
            'MOOSE is required for NTEC/electrochemical modes; set '
            'MOOSE_EXECUTABLE or install environments/environment-moose.yml')
    return None


def moose_environment(executable: str | None) -> dict[str, str]:
    """Return a subprocess environment that runs MOOSE with its own toolchain.

    MOOSE JIT-compiles parsed expressions (``ADParsedMaterial``) with the first
    ``mpicxx`` on ``PATH``. When another MPI stack is active, e.g. the NVIDIA
    HPC SDK sourced for GPU Quantum ESPRESSO, that compiler produces objects
    MOOSE cannot load, so the executable's own ``bin`` directory goes first.

    Args:
        executable: Resolved MOOSE executable, or None to leave PATH unchanged.

    Returns:
        A copy of ``os.environ`` with the MOOSE toolchain first on ``PATH``.
    """
    environment = os.environ.copy()
    if executable:
        toolchain = str(Path(executable).resolve().parent)
        environment['PATH'] = os.pathsep.join(
            [toolchain, *[p for p in environment.get('PATH', '').split(os.pathsep)
                          if p and p != toolchain]])
    return environment


def moose_version(executable: str) -> str:
    """Return the MOOSE build identity reported by ``--version``.

    Args:
        executable: MOOSE executable path.

    Returns:
        A non-empty version string for provenance.
    """
    probe = subprocess.run(
        [executable, '--version'], capture_output=True, text=True, timeout=60)
    lines = [line.strip() for line in (probe.stdout + probe.stderr).splitlines()
             if line.strip()]
    if probe.returncode or not lines:
        raise RuntimeError(f'MOOSE executable did not report a version: {executable}')
    return ' '.join(lines[:2])


def default_moose_input(reactor_type: str) -> Path:
    """Return the repository-owned input template for a reactor type.

    Args:
        reactor_type: ``NTEC`` or ``Electrochemical``.

    Returns:
        Path of the MOOSE input template.
    """
    try:
        return DEFAULT_MOOSE_INPUTS[reactor_type]
    except KeyError as exc:
        raise ValueError(f'{reactor_type} has no MOOSE model') from exc


# Case sections whose numeric fields become MOOSE variables of the same name.
CASE_PARAMETER_SECTIONS = ('geometry', 'operating', 'properties', 'model_parameters')
COUPLING_FEEDBACK = 'hydrogen_interface_rates.json'


def moose_case_parameters(physical_case: Mapping, case: Path,
                          iteration: int) -> dict[str, object]:
    """Derive MOOSE variables from a validated physical case.

    Every finite numeric field in ``geometry``, ``operating``, ``properties``
    and the optional ``model_parameters`` section is passed under its own
    unit-bearing name (for example ``electrode_area_m2``). Duplicate names
    across sections are rejected rather than silently overwritten. The
    Cantera interface rates from the previous iteration are passed as scalar
    variables when present, so iteration 1 runs with the declared initial
    guess and later iterations consume the exchanged kinetics.

    Args:
        physical_case: Case returned by ``load_physical_case``.
        case: Case directory holding the coupling feedback file.
        iteration: One-based coupling iteration.

    Returns:
        Variable name to value mapping for ``build_moose_command``.
    """
    parameters: dict[str, object] = {}
    for section in CASE_PARAMETER_SECTIONS:
        for name, value in dict(physical_case.get(section, {})).items():
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                continue
            if name in parameters:
                raise ValueError(f'case field {name!r} is declared twice')
            parameters[name] = float(value)
    parameters['coupling_iteration'] = int(iteration)
    feedback = case / COUPLING_FEEDBACK
    if iteration > 1:
        if not feedback.is_file():
            raise RuntimeError(f'coupling feedback missing before iteration {iteration}')
        exchanged = json.loads(feedback.read_text()).get('moose_parameters', {})
        for name, value in exchanged.items():
            if name in parameters:
                raise ValueError(f'coupling feedback overrides case field {name!r}')
            parameters[name] = float(value)
    return parameters


def _format_value(value) -> str:
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, (int, float)):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError('MOOSE parameters must be finite')
        return repr(number)
    text = str(value)
    if not text or any(char.isspace() for char in text) or '=' in text:
        raise ValueError(f'MOOSE string parameter is not a single token: {text!r}')
    return text


def moose_overrides(parameters: Mapping[str, object]) -> list[str]:
    """Convert named parameters into deterministic ``name=value`` overrides.

    Args:
        parameters: Top-level HIT variable names mapped to values.

    Returns:
        Sorted command-line override tokens.
    """
    tokens = []
    for name in sorted(parameters):
        if not name.replace('_', '').isalnum():
            raise ValueError(f'invalid MOOSE parameter name: {name!r}')
        tokens.append(f'{name}={_format_value(parameters[name])}')
    return tokens


def build_moose_command(executable: str, input_file: Path, case: Path,
                        parameters: Mapping[str, object],
                        processes: int = 1) -> list[str]:
    """Build a shell-free MOOSE command and record its parameters in the case.

    Args:
        executable: MOOSE executable path.
        input_file: Byte-stable MOOSE input template.
        case: Case directory used as the working directory.
        parameters: Top-level HIT variable overrides derived from the case.
        processes: MPI ranks; ``1`` runs serially without ``mpiexec``.

    Returns:
        Command tokens.
    """
    overrides = moose_overrides(parameters)
    # Cases legitimately carry fields a given model does not consume. MOOSE
    # resolves a relative file_base against the input file, so use the case.
    command = [executable, '-i', str(input_file), '--color', 'off', '--allow-unused',
               f'Outputs/file_base={case / "hydrogen_moose"}', *overrides]
    if processes > 1:
        mpiexec = resolve_executable('mpiexec', conda_env=MOOSE_CONDA_ENV)
        command = [mpiexec, '-n', str(int(processes)), *command]
    (case / PARAMETER_RECORD).write_text(json.dumps(
        {'schema_version': 1, 'input_file': input_file.name,
         'overrides': overrides, 'command': command}, sort_keys=True) + '\n')
    return command
