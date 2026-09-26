"""Portable discovery of external executables used by pipeline stages."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


def _usable(candidate: str | None) -> str | None:
    if not candidate:
        return None
    found = shutil.which(candidate)
    if found:
        # Preserve the invoked symlink name. Multicall launchers such as the
        # NVIDIA HPC-X ``mpirun -> env.sh`` wrapper dispatch through argv[0]
        # and stop working if canonicalized to their implementation target.
        return str(Path(found).absolute())
    path = Path(candidate).expanduser()
    if path.is_file() and os.access(path, os.X_OK):
        return str(path.absolute())
    return None


def resolve_executable(name: str, *, env_var: str | None = None,
                       conda_env: str | None = None,
                       required: bool = True) -> str | None:
    """Resolve without assuming where Python, Conda, or environments live.

        Resolution order is an explicit environment override, the caller's PATH,
        then a PATH-resolved Conda executable querying a documented environment.

    Args:
        name: Name used by this operation.
        env_var: Env var used by this operation.
        conda_env: Conda env used by this operation.
        required: Whether to enable required.

    Returns:
        Computed `str | None` result.
    """
    override = os.environ.get(env_var, '') if env_var else ''
    if override:
        resolved = _usable(override)
        if resolved:
            return resolved
        raise RuntimeError(
            f'{env_var} is set but is not an executable file or PATH command: '
            f'{override!r}')

    resolved = _usable(name)
    if resolved:
        return resolved

    conda = shutil.which('conda')
    if conda and conda_env:
        probe = subprocess.run(
            [conda, 'run', '-n', conda_env, 'which', name],
            capture_output=True, text=True, timeout=30)
        if probe.returncode == 0:
            lines = [line.strip() for line in probe.stdout.splitlines()
                     if line.strip()]
            if lines:
                resolved = _usable(lines[-1])
                if resolved:
                    return resolved

    if not required:
        return None
    override_help = f' set {env_var},' if env_var else ''
    conda_help = (f' or install it in the documented {conda_env!r} Conda environment'
                  if conda_env else '')
    raise RuntimeError(
        f'Unable to locate {name!r};{override_help} add it to PATH{conda_help}. '
        'See README.md Environment Setup.')


def resolve_qe_executable(name: str) -> str:
    """Resolve only an explicitly activated, pinned GPU QE executable.

    Production QE intentionally does not fall back to ``PATH`` or Conda. Those
    locations can contain a CPU-only build with the same executable name. The
    Blackwell installer writes ``PW_X`` and ``NEB_X`` into its activation file,
    making the selected native toolchain explicit and auditable.

    Args:
        name: Human-readable identifier used in diagnostics and output.

    Returns:
        A `str` containing the resolve qe executable result.
    """
    variables = {'pw.x': 'PW_X', 'neb.x': 'NEB_X'}
    variable = variables.get(name)
    if variable is None:
        raise ValueError(f'unsupported Quantum ESPRESSO executable: {name!r}')
    override = os.environ.get(variable, '')
    if not override:
        raise RuntimeError(
            f'{variable} is not set. Run scripts/install_qe_gpu.sh and source '
            'the generated activate.sh; generic PATH and Conda QE builds are '
            'not accepted for production.')
    resolved = _usable(override)
    if not resolved:
        raise RuntimeError(
            f'{variable} does not identify an executable: {override!r}')
    return resolved
