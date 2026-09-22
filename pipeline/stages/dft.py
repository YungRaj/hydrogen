"""Independent Quantum ESPRESSO candidate-validation stage."""

from __future__ import annotations

import ast
from typing import Callable, Iterable, Mapping, MutableMapping

from pipeline.stages.contracts import StageOutcome


def _indexed_rows(candidates):
    if hasattr(candidates, 'iterrows'):
        yield from candidates.iterrows()
    else:
        yield from enumerate(candidates)


def _row_as_mapping(row):
    """Accept a dict, Mapping, or pandas Series. A Series is not a Mapping."""
    if isinstance(row, MutableMapping) or isinstance(row, Mapping):
        return row
    if hasattr(row, 'to_dict'):
        return row.to_dict()
    return None


def _genome_from_row(row):
    mapping = _row_as_mapping(row)
    if mapping is None:
        return row
    if 'genome' in mapping:
        return mapping['genome']
    raise KeyError('candidate row has no genome')


def run_dft_stage(candidates, *, top_k: int, execute_dft: bool,
                  validator: Callable | None = None,
                  name_prefix: str = 'dft_cat',
                  error_sink: Callable[[str], None] | None = None) -> StageOutcome:
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
    if not isinstance(top_k, int) or top_k < 0:
        raise ValueError('top_k must be a nonnegative integer')
    if validator is None:
        from pipeline.validation.dft_validator import validate_catalyst
        validator = validate_catalyst
    selected = candidates.head(top_k) if hasattr(candidates, 'head') \
        else list(candidates)[:top_k]
    results, failures = [], []
    for index, row in _indexed_rows(selected):
        try:
            genome = _genome_from_row(row)
            if isinstance(genome, str):
                genome = ast.literal_eval(genome)
            result = validator(
                f'{name_prefix}_{index}', genome, run_dft=execute_dft)
            results.append(result)
        except Exception as exc:
            message = f'DFT failed for {name_prefix}_{index}: {exc}'
            failures.append({
                'candidate': f'{name_prefix}_{index}',
                'error_type': type(exc).__name__, 'error': str(exc)})
            if error_sink:
                error_sink(message)
    return StageOutcome(
        state={
            'n_validated': len(results),
            'n_converged': sum(
                1 for result in results if result.get('converged', False)),
            'n_failed': len(failures),
        },
        products={'dft_results': results, 'failures': failures})
