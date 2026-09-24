"""Persisted top-level pipeline state and phase-record shapes."""

from __future__ import annotations

import math
from typing import Any, Mapping, TypedDict, cast


class _PipelineStateOptional(TypedDict, total=False):
    phase1: dict[str, Any]
    phase2: dict[str, Any]
    phase3: dict[str, Any]
    phase4: dict[str, Any]
    phase5: dict[str, Any]
    phase6: dict[str, Any]
    total_elapsed_s: float


class PipelineStateDocument(_PipelineStateOptional):
    """Versionable persisted summary assembled by the master coordinator."""


def validate_pipeline_state(value: object) -> PipelineStateDocument:
    """Validate a loaded state document before orchestration mutates it.

    Args:
        value: Deserialized state read from the configured persistence adapter.

    Returns:
        The original mapping copied into a typed mutable state document.
    """
    if not isinstance(value, Mapping):
        raise ValueError('pipeline state must be a mapping')
    state = dict(cast(Mapping[str, Any], value))
    for phase in ('phase1', 'phase2', 'phase3', 'phase4', 'phase5', 'phase6'):
        if phase in state and not isinstance(state[phase], Mapping):
            raise ValueError(f'{phase} pipeline state must be a mapping')
    elapsed = state.get('total_elapsed_s')
    if elapsed is not None:
        if (not isinstance(elapsed, (int, float)) or isinstance(elapsed, bool)
                or not math.isfinite(float(elapsed)) or float(elapsed) < 0):
            raise ValueError('pipeline total_elapsed_s must be finite and nonnegative')
        state['total_elapsed_s'] = float(elapsed)
    return cast(PipelineStateDocument, state)
