"""Independent Quantum ESPRESSO candidate-validation stage."""

from __future__ import annotations

import ast
from collections.abc import Callable, Iterator, Mapping, MutableMapping
from typing import Any, cast

from pipeline.data_models.quantum import DFTResult
from pipeline.data_models.stages import DFTProducts, DFTState
from pipeline.stages.contracts import StageOutcome


def _indexed_rows(candidates: Any) -> Iterator[tuple[Any, Any]]:
    if hasattr(candidates, 'iterrows'):
        yield from candidates.iterrows()
    else:
        yield from enumerate(candidates)


def _row_as_mapping(row: Any) -> Mapping[str, Any] | None:
    """Accept a dict, Mapping, or pandas Series. A Series is not a Mapping."""
    if isinstance(row, (MutableMapping, Mapping)):
        return cast(Mapping[str, Any], row)
    if hasattr(row, 'to_dict'):
        return cast(Mapping[str, Any], row.to_dict())
    return None


def _genome_from_row(row: Any) -> Any:
    mapping = _row_as_mapping(row)
    if mapping is None:
        return row
    if 'genome' in mapping:
        return mapping['genome']
    raise KeyError('candidate row has no genome')


def _job_identity(genome: Any) -> str:
    """Return a stable, path-safe job identity for one candidate.

    Positional names (``dft_cat_0``) collide across reruns with a different
    ranking; the content-addressed candidate ID (the same hash recorded by
    ``validate_catalyst`` as ``candidate_id``) does not.
    """
    from pipeline.search.discovery import (
        candidate_id,  # pyright: ignore[reportUnknownVariableType]
    )

    return str(candidate_id(genome))


def run_dft_stage(
    candidates: Any,
    *,
    top_k: int,
    execute_dft: bool,
    validator: Callable[..., DFTResult] | None = None,
    name_prefix: str = 'dft_cat',
    error_sink: Callable[[str], None] | None = None,
) -> StageOutcome[DFTState, DFTProducts]:
    """Validate candidate genomes through an injectable QE workflow boundary.

    Args:
        candidates: Candidate records to process.
        top_k: Bound controlling top k.
        execute_dft: Whether to enable execute dft.
        validator: Injected callable used to perform validator.
        name_prefix: Name prefix used by this operation.
        error_sink: Injected callable used to perform error sink.

    Returns:
        Computed `StageOutcome` result.
    """
    if not isinstance(top_k, int) or top_k < 0:  # pyright: ignore[reportUnnecessaryIsInstance]
        raise ValueError('top_k must be a nonnegative integer')
    active_validator: Callable[..., DFTResult]
    if validator is None:
        from pipeline.validation.dft_validator import (
            validate_catalyst,  # pyright: ignore[reportUnknownVariableType]
        )

        active_validator = cast(Callable[..., DFTResult], validate_catalyst)
    else:
        active_validator = validator
    selected = (
        candidates.head(top_k)
        if hasattr(candidates, 'head')
        else list(candidates)[:top_k]
    )
    results: list[DFTResult] = []
    failures: list[dict[str, str]] = []
    for index, row in _indexed_rows(selected):
        name = f'{name_prefix}_{index}'
        try:
            genome = _genome_from_row(row)
            if isinstance(genome, str):
                genome = ast.literal_eval(genome)
            name = f'{name_prefix}_{_job_identity(genome)}'
            result = active_validator(name, genome, run_dft=execute_dft)
            results.append(result)
        except Exception as exc:
            message = f'DFT failed for {name}: {exc}'
            failures.append(
                {
                    'candidate': name,
                    'error_type': type(exc).__name__,
                    'error': str(exc),
                }
            )
            if error_sink:
                error_sink(message)
    state: DFTState = {
        'n_validated': len(results),
        'n_converged': sum(1 for result in results if result.get('converged', False)),
        'n_failed': len(failures),
    }
    products: DFTProducts = {'dft_results': results, 'failures': failures}
    return StageOutcome(state=state, products=products)
