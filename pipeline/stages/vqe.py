"""Independent CUDA-Q VQE validation stage."""

from __future__ import annotations

from typing import Callable

from pipeline.stages.contracts import StageOutcome


def run_vqe_stage(*, top_k: int, execute_quantum: bool,
                  validator: Callable | None = None) -> StageOutcome:
    """Run candidate VQE validations with an injectable solver boundary.

    Args:
        top_k: Bound controlling top k.
        execute_quantum: Whether to enable execute quantum.
        validator: Injected callable used to perform validator.

    Returns:
        Computed `StageOutcome` result.
    """
    if not isinstance(top_k, int) or top_k < 0:
        raise ValueError('top_k must be a nonnegative integer')
    if validator is None:
        from pipeline.validation.vqe_transition_state import (
            validate_transition_state)
        validator = validate_transition_state
    target = 'nvidia' if execute_quantum else 'default'
    results = [
        validator(f'champion_{index}', 'CH_split', target=target)
        for index in range(min(top_k, 3))
    ]
    return StageOutcome(
        state={'n_vqe_runs': len(results)},
        products={'vqe_results': results})
