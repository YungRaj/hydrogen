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
        return str(Path(found).resolve())
    path = Path(candidate).expanduser()
    if path.is_file() and os.access(path, os.X_OK):
        return str(path.resolve())
    return None


def resolve_executable(name: str, *, env_var: str | None = None,
                       conda_env: str | None = None,
                       required: bool = True) -> str | None:
    """Resolve without assuming where Python, Conda, or environments live.

    Resolution order is an explicit environment override, the caller's PATH,
    then a PATH-resolved Conda executable querying a documented environment.
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
    variables = {'pw.x': 'PW_X', 'neb.x': 'NEB_X'}
    return resolve_executable(
        name, env_var=variables.get(name), conda_env='qe-env', required=True)
