"""Persisted top-level pipeline state and phase-record shapes."""

from __future__ import annotations

import math
from typing import Any, Mapping, TypedDict, cast

from pipeline.data_models.stages import (
    DFTState, DiscoveryState, FuelCellState, ReactorBatchState, ReportState,
    VQEState)


PIPELINE_STATE_SCHEMA_VERSION = 1


class DiscoveryPhaseRecord(DiscoveryState):
    """Persisted discovery summary plus coordinator timing."""
    elapsed_s: float


class ReactorPhaseRecord(ReactorBatchState):
    """Persisted reactor summary plus coordinator timing."""
    elapsed_s: float


class DFTPhaseRecord(DFTState):
    """Persisted DFT summary plus coordinator timing."""
    elapsed_s: float


class VQEPhaseRecord(VQEState):
    """Persisted VQE summary plus coordinator timing."""
    elapsed_s: float


class FuelCellPhaseRecord(FuelCellState):
    """Persisted fuel-cell summary plus coordinator timing."""
    elapsed_s: float


class ReportPhaseRecord(ReportState):
    """Persisted report summary plus coordinator timing."""
    elapsed_s: float


class _PipelineStateOptional(TypedDict, total=False):
    phase1: DiscoveryPhaseRecord
    phase2: ReactorPhaseRecord
    phase3: DFTPhaseRecord
    phase4: VQEPhaseRecord
    phase5: FuelCellPhaseRecord
    phase6: ReportPhaseRecord
    total_elapsed_s: float


class PipelineStateDocument(_PipelineStateOptional):
    """Versionable persisted summary assembled by the master coordinator."""
    schema_version: int


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
    version = state.get('schema_version', 0)
    if not isinstance(version, int) or isinstance(version, bool) or version < 0:
        raise ValueError('pipeline schema_version must be a nonnegative integer')
    if version > PIPELINE_STATE_SCHEMA_VERSION:
        raise ValueError(
            f'pipeline state schema {version} is newer than supported '
            f'version {PIPELINE_STATE_SCHEMA_VERSION}')
    # Version 0 is the historical unversioned shape. Its phase dictionaries
    # already match version 1, so migration is explicit and lossless.
    state['schema_version'] = PIPELINE_STATE_SCHEMA_VERSION
    for phase in ('phase1', 'phase2', 'phase3', 'phase4', 'phase5', 'phase6'):
        if phase in state and not isinstance(state[phase], Mapping):
            raise ValueError(f'{phase} pipeline state must be a mapping')
        record = state.get(phase)
        if isinstance(record, Mapping) and 'elapsed_s' in record:
            typed_record = cast(Mapping[str, object], record)
            phase_elapsed = typed_record['elapsed_s']
            if (not isinstance(phase_elapsed, (int, float)) or
                    isinstance(phase_elapsed, bool) or
                    not math.isfinite(float(phase_elapsed)) or
                    float(phase_elapsed) < 0):
                raise ValueError(
                    f'{phase} elapsed_s must be finite and nonnegative')
    elapsed = state.get('total_elapsed_s')
    if elapsed is not None:
        if (not isinstance(elapsed, (int, float)) or isinstance(elapsed, bool)
                or not math.isfinite(float(elapsed)) or float(elapsed) < 0):
            raise ValueError('pipeline total_elapsed_s must be finite and nonnegative')
        state['total_elapsed_s'] = float(elapsed)
    return cast(PipelineStateDocument, state)
