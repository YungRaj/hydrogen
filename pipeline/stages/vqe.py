"""Independent CUDA-Q VQE validation stage."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

from pipeline.data_models.quantum import VQEResult
from pipeline.data_models.stages import VQEProducts, VQEState
from pipeline.stages.contracts import StageOutcome


def _candidate_rows(candidates: Any, limit: int) -> list[Any]:
    if candidates is None or limit <= 0:
        return []
    if hasattr(candidates, 'head'):
        return [row for _, row in candidates.head(limit).iterrows()]
    return list(candidates)[:limit]


def _candidate_key(row: Any) -> str:
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


def _extract_row_genome(row: Any) -> tuple[object, ...] | None:
    """Return a parsed genome tuple from a candidate row if present."""
    import ast

    raw_g = row.get('genome') if hasattr(row, 'get') else None
    if raw_g is None:
        return None
    if isinstance(raw_g, (tuple, list)):
        return tuple(tuple(x) if isinstance(x, list) else x for x in raw_g)
    if isinstance(raw_g, str) and raw_g.strip():
        try:
            parsed = ast.literal_eval(raw_g.strip())
            if isinstance(parsed, (tuple, list)) and parsed:
                return tuple(
                    tuple(x) if isinstance(x, list) else x for x in parsed
                )
        except (ValueError, SyntaxError):
            return None
    return None


def run_vqe_stage(
    candidates: Any = None,
    *,
    top_k: int,
    execute_quantum: bool,
    validator: Callable[..., VQEResult] | None = None,
    hamiltonian_root: str | Path | None = None,
    auto_generate_hamiltonian: bool | None = None,
) -> StageOutcome[VQEState, VQEProducts]:
    """Run candidate VQE validations with an injectable solver boundary.

    Only candidates with a provenance-bound, candidate-specific Hamiltonian at
    ``<hamiltonian_root>/<candidate_id>/FCIDUMP`` are run. When a candidate row
    carries a valid ``genome`` and ``auto_generate_hamiltonian`` is enabled
    (defaults to ``execute_quantum``), ``generate_catalyst_ts_fcidump`` is
    invoked automatically to construct the 8-qubit ``CASSCF(4,4)`` transition-state
    ``FCIDUMP`` and ``FCIDUMP.json`` provenance sidecar on disk before VQE runs.

    Args:
        candidates: Ranked candidate table or sequence of mappings carrying
            ``candidate_id``; ``None`` means no candidates are available.
        top_k: Bound controlling how many leading candidates are considered.
        execute_quantum: Whether to use the GPU CUDA-Q target.
        validator: Injected callable used to perform validator.
        hamiltonian_root: Directory holding ``<candidate_id>/FCIDUMP`` inputs;
            defaults to ``<results>/vqe/hamiltonians``.
        auto_generate_hamiltonian: Optional override controlling automatic
            PySCF ``CASSCF(4,4)`` ``FCIDUMP`` generation when a candidate row
            includes a valid ``genome``; defaults to ``execute_quantum``.

    Returns:
        Computed `StageOutcome` result.
    """
    if not isinstance(top_k, int) or top_k < 0:  # pyright: ignore[reportUnnecessaryIsInstance]
        raise ValueError('top_k must be a nonnegative integer')
    active_validator: Callable[..., VQEResult]
    if validator is None:
        from pipeline.validation.vqe_transition_state import (
            validate_transition_state,  # pyright: ignore[reportUnknownVariableType]
        )

        active_validator = cast(Callable[..., VQEResult], validate_transition_state)
    else:
        active_validator = validator
    root = (
        Path(hamiltonian_root)
        if hamiltonian_root is not None
        else _default_hamiltonian_root()
    )
    should_auto_generate = (
        bool(execute_quantum)
        if auto_generate_hamiltonian is None
        else bool(auto_generate_hamiltonian)
    )
    # Use GPU statevector ('nvidia') when quantum execution is enabled, else CPU ('default')
    target = 'nvidia' if execute_quantum else 'default'
    ready: list[tuple[str, Path]] = []
    missing: list[tuple[str, Path]] = []
    for row in _candidate_rows(candidates, top_k):
        key = _candidate_key(row)
        hamiltonian = root / key / 'FCIDUMP'
        if not hamiltonian.is_file() and should_auto_generate:
            parsed_genome = _extract_row_genome(row)
            if parsed_genome is not None:
                try:
                    from pipeline.validation.candidate_hamiltonian import (
                        generate_catalyst_ts_fcidump,
                    )

                    generate_catalyst_ts_fcidump(
                        key,
                        parsed_genome,
                        output_dir=root / key,
                    )
                except Exception:
                    pass
        (ready if hamiltonian.is_file() else missing).append((key, hamiltonian))
    results: list[VQEResult] = [
        active_validator(
            key,
            'CH_split',  # Rate-determining first C-H cleavage transition state active space
            target=target,
            candidate_hamiltonian=str(hamiltonian),
        )
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
