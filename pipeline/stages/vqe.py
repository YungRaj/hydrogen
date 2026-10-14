"""Independent CUDA-Q VQE validation stage."""

from __future__ import annotations

from pathlib import Path
from typing import Callable

from pipeline.data_models.stages import VQEProducts, VQEState
from pipeline.stages.contracts import StageOutcome


def _candidate_rows(candidates, limit: int):
    if candidates is None or limit <= 0:
        return []
    if hasattr(candidates, 'head'):
        return [row for _, row in candidates.head(limit).iterrows()]
    return list(candidates)[:limit]


def _candidate_key(row) -> str:
    """Return the path-safe candidate identity used to locate its FCIDUMP."""
    value = row.get('candidate_id') if hasattr(row, 'get') else None
    if value is None or (isinstance(value, float) and value != value):
        raise ValueError('VQE candidate row has no candidate_id')
    key = str(value)
    if not key or Path(key).name != key or key in {'.', '..'}:
        raise ValueError(f'candidate_id is not a safe path component: {key!r}')
    return key


def _default_hamiltonian_root() -> Path:
    from pipeline.utils import VQE_DIR

    return VQE_DIR / 'hamiltonians'


def run_vqe_stage(
    candidates=None,
    *,
    top_k: int,
    execute_quantum: bool,
    validator: Callable | None = None,
    hamiltonian_root: str | Path | None = None,
) -> StageOutcome[VQEState, VQEProducts]:
    """Run candidate VQE validations with an injectable solver boundary.

    Only candidates with a provenance-bound, candidate-specific Hamiltonian at
    ``<hamiltonian_root>/<candidate_id>/FCIDUMP`` are run. Without one the
    validator would fall back to a toy benchmark Hamiltonian that carries no
    candidate identity, so those candidates are reported as missing inputs
    instead of producing placeholder ``champion_<i>`` results.

    Args:
        candidates: Ranked candidate table or sequence of mappings carrying
            ``candidate_id``; ``None`` means no candidates are available.
        top_k: Bound controlling how many leading candidates are considered.
        execute_quantum: Whether to use the GPU CUDA-Q target.
        validator: Injected callable used to perform validator.
        hamiltonian_root: Directory holding ``<candidate_id>/FCIDUMP`` inputs;
            defaults to ``<results>/vqe/hamiltonians``.

    Returns:
        Computed `StageOutcome` result.
    """
    if not isinstance(top_k, int) or top_k < 0:
        raise ValueError('top_k must be a nonnegative integer')
    if validator is None:
        from pipeline.validation.vqe_transition_state import validate_transition_state

        validator = validate_transition_state
    root = (Path(hamiltonian_root) if hamiltonian_root is not None
            else _default_hamiltonian_root())
    target = 'nvidia' if execute_quantum else 'default'
    ready, missing = [], []
    for row in _candidate_rows(candidates, top_k):
        key = _candidate_key(row)
        hamiltonian = root / key / 'FCIDUMP'
        (ready if hamiltonian.is_file() else missing).append((key, hamiltonian))
    results = [
        validator(key, 'CH_split', target=target,
                  candidate_hamiltonian=str(hamiltonian))
        for key, hamiltonian in ready
    ]
    state: VQEState = {
        'n_vqe_runs': len(results),
        'candidates': [key for key, _ in ready],
        'n_missing_hamiltonians': len(missing),
    }
    products: VQEProducts = {
        'vqe_results': results,
        'missing_hamiltonians': [str(path) for _, path in missing],
    }
    return StageOutcome(state=state, products=products)
