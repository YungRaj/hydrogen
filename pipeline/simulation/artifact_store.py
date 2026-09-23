"""Persistence boundary for fail-closed multiphysics artifacts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Mapping


def persist_validated_artifact(
        artifact: Mapping, *, results_dir: str | Path, candidate_id: str,
        mode: str, reactor_type: str, temperature_K: float,
        path_builder: Callable, validator: Callable) -> Path:
    """Write by atomic replacement and retain only validator-accepted artifacts.

        The target is removed if validation fails, preventing an invalid partial
        result from being mistaken for reusable solver evidence.

    Args:
        artifact: Mapping supplying artifact.
        results_dir: Directory containing or receiving calculation results.
        candidate_id: Stable candidate identifier.
        mode: Configured methane-conversion pathway.
        reactor_type: Physical reactor implementation identifier.
        temperature_K: Absolute temperature in kelvin.
        path_builder: Injected callable used to perform path builder.
        validator: Injected callable used to perform validator.

    Returns:
        Filesystem path produced or resolved by the operation.
    """
    target = path_builder(
        results_dir, candidate_id, mode, reactor_type, temperature_K)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + '.tmp')
    temporary.write_text(json.dumps(dict(artifact), indent=2, sort_keys=True) + '\n')
    temporary.replace(target)
    validated = validator(
        results_dir, candidate_id, mode, reactor_type, temperature_K)
    if not validated['valid']:
        target.unlink(missing_ok=True)
        failed = ', '.join(validated.get('failed_checks', ()))
        raise RuntimeError(
            'solver output did not satisfy the artifact contract: ' +
            (failed or validated['reason']))
    return target
