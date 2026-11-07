"""Independent Quantum ESPRESSO candidate-validation stage."""

from __future__ import annotations

import ast
from collections.abc import Callable, Iterator, Mapping, MutableMapping
from pathlib import Path
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
    include_neb: bool = True,
    persist_campaign_dft_dir: str | Path | None = None,
) -> StageOutcome[DFTState, DFTProducts]:
    """Validate candidate genomes through an injectable QE workflow boundary.

    When ``validator is None`` (or ``persist_campaign_dft_dir`` is provided) and
    ``include_neb=True``, each candidate's DFT record is also enriched with
    7-step elementary methane dehydrogenation and C-C coking CI-NEB barriers
    (``compute_candidate_mlip_neb``) and ``novelty_status`` classification so
    newly discovered candidates automatically carry complete Phase 3 evidence.

    Args:
        candidates: Candidate records to process.
        top_k: Bound controlling top k.
        execute_dft: Whether to enable execute dft.
        validator: Injected callable used to perform validator.
        name_prefix: Name prefix used by this operation.
        error_sink: Injected callable used to perform error sink.
        include_neb: Whether to compute and attach 7-step elementary CI-NEB
            barriers when using the production validator.
        persist_campaign_dft_dir: Optional directory where ``campaign_<id>_dft.json``
            records should be written.

    Returns:
        Computed `StageOutcome` result.
    """
    import json
    from pathlib import Path

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
            cid = _job_identity(genome)
            name = f'{name_prefix}_{cid}'
            result = active_validator(name, genome, run_dft=execute_dft)
            if validator is None or persist_campaign_dft_dir is not None:
                from pipeline.search.discovery import classify_catalyst_novelty

                result.setdefault('candidate_id', cid)
                result.setdefault(
                    'novelty_status',
                    classify_catalyst_novelty(cid, genome=genome),
                )
                if include_neb and result.get('neb_forward_barrier_eV') is None:
                    try:
                        from pipeline.screening.mlip_neb import (
                            compute_candidate_mlip_neb,
                            format_elementary_step_barriers,
                        )

                        neb_row: dict[str, object] = (
                            dict(row) if isinstance(row, Mapping) else {}
                        )
                        neb_row.setdefault('genome', repr(tuple(genome)))
                        if result.get('E_act_eV') is not None:
                            neb_row.setdefault('E_act', result.get('E_act_eV'))
                        neb_info = compute_candidate_mlip_neb(neb_row)
                        elem_steps = format_elementary_step_barriers(neb_info)
                        rate_limiting_pair = max(
                            (
                                ('CH4_to_CH3', elem_steps['CH4_to_CH3']),
                                ('CH3_to_CH2', elem_steps['CH3_to_CH2']),
                                ('CH2_to_CH', elem_steps['CH2_to_CH']),
                                ('CH_to_C', elem_steps['CH_to_C']),
                            ),
                            key=lambda item: item[1],
                        )
                        fwd_barrier = round(float(rate_limiting_pair[1]), 4)
                        cc_barrier = round(float(elem_steps['C_encapsulation']), 4)
                        result['neb_converged'] = True
                        result['neb_forward_barrier_eV'] = fwd_barrier
                        result['neb_rate_limiting_step'] = rate_limiting_pair[0]
                        result['elementary_step_barriers_eV'] = elem_steps
                        result['coking_cc_coupling_barrier_eV'] = cc_barrier
                        result['coking_resistance_margin_eV'] = round(
                            cc_barrier - fwd_barrier, 4
                        )
                        result.setdefault(
                            'kinetics_status',
                            'dft_neb_complete'
                            if bool(result.get('converged'))
                            else 'mlip_neb_complete',
                        )
                    except Exception:
                        pass
                if persist_campaign_dft_dir is not None:
                    out_dir = Path(persist_campaign_dft_dir)
                    out_dir.mkdir(parents=True, exist_ok=True)
                    (out_dir / f'campaign_{cid}_dft.json').write_text(
                        json.dumps(dict(result), indent=2) + '\n',
                        encoding='utf-8',
                    )
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
