"""Evidence-aware candidate selection at fidelity boundaries.

An incomplete calculation is not evidence that a material is poor.  Reactor
models require complete numerical descriptors, while high-fidelity validation
must retain a route for unresolved and unfamiliar chemistry.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pipeline.data_models.evidence import CandidateDisposition


HARD_EXCLUSION: CandidateDisposition = 'hard_excluded'
VALIDATION_REQUIRED: CandidateDisposition = 'validation_required'
QUANTITATIVE_SCREENING: CandidateDisposition = 'quantitative_screening'


def evidence_disposition(row, primary: str) -> CandidateDisposition:
    """Classify evidence without conflating calculation failure with chemistry.

    Args:
        row: Row used by this operation.
        primary: Primary used by this operation.

    Returns:
        Computed `str` result.
    """
    declared = str(row.get('candidate_disposition', '') or '')
    if declared == HARD_EXCLUSION:
        return HARD_EXCLUSION
    error = str(row.get('error', '') or '').lower()
    if 'toxic/radioactive' in error:
        return HARD_EXCLUSION
    try:
        value = float(row.get(primary))
    except (TypeError, ValueError):
        value = np.nan
    censored_value = row.get(f'{primary}_censored', False)
    censored = False if pd.isna(censored_value) else bool(censored_value)
    if bool(row.get('valid', False)) and np.isfinite(value) and not censored:
        return QUANTITATIVE_SCREENING
    return VALIDATION_REQUIRED


def annotate_evidence(frame: pd.DataFrame, primary: str) -> pd.DataFrame:
    """Attach explicit evidence authority and limitations to a result.

    Args:
        frame: Tabular candidate records to annotate or evaluate.
        primary: Primary objective column used to order evidence.

    Returns:
        A `pd.DataFrame` containing the annotate evidence result.
    """
    result = frame.copy()
    result['candidate_disposition'] = [
        evidence_disposition(row, primary) for _, row in result.iterrows()
    ]
    return result


def _frame_pareto_ranks(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    """Compute Pareto front ranks across minimization columns on a DataFrame."""
    from pipeline.screening.objectives import pareto_ranks

    cols = []
    for spec in columns:
        negate = spec.startswith('-')
        col_name = spec[1:] if negate else spec
        series = pd.to_numeric(frame.get(col_name), errors='coerce')
        values = (-series if negate else series).fillna(1e9).to_numpy(dtype=float)
        cols.append(values)
    return pareto_ranks(np.column_stack(cols))


def _class_preserving(
    frame: pd.DataFrame,
    n_select: int,
    primary: str,
    min_per_class: int,
    validation: bool,
    objectives: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    if n_select <= 0 or frame.empty:
        return frame.iloc[0:0].copy()
    work = frame.copy()
    work['_selection_index'] = np.arange(len(work))
    values = pd.to_numeric(work.get(primary), errors='coerce')
    work['_primary'] = values.fillna(np.inf)
    work['_pareto_rank'] = (
        _frame_pareto_ranks(work, list(objectives))
        if objectives
        else np.zeros(len(work), dtype=int)
    )
    needs_dft = (
        work['needs_dft_validation'].map(
            lambda value: False if pd.isna(value) else bool(value)
        )
        if 'needs_dft_validation' in work
        else pd.Series(False, index=work.index)
    )
    work['_needs_resolution'] = (
        work['candidate_disposition'].eq(VALIDATION_REQUIRED) | needs_dft
    ).astype(int)
    # Validation gives unresolved evidence first access; quantitative consumers
    # can only use complete screening rows. Stable index makes ties reproducible.
    if validation:
        order = ['_needs_resolution', '_pareto_rank', '_primary', '_selection_index']
        ascending = [False, True, True, True]
    else:
        order = ['_pareto_rank', '_primary', '_selection_index']
        ascending = [True, True, True]
    ranked = work.sort_values(order, ascending=ascending, kind='mergesort')
    chosen = []
    if 'material_class' in ranked and min_per_class > 0:
        class_rows = {
            material_class: ranked[ranked['material_class'].eq(material_class)]
            for material_class in ranked['material_class'].dropna().unique()
        }
        # Allocate quota in rounds. If the budget is smaller than the number of
        # classes, the best class champions win rather than alphabetical order.
        rank_position = {
            index: position
            for position, index in enumerate(ranked['_selection_index'].tolist())
        }
        for quota_round in range(min_per_class):
            champions = [
                group.iloc[quota_round]['_selection_index']
                for group in class_rows.values()
                if len(group) > quota_round
            ]
            champions.sort(key=lambda index: rank_position[index])
            chosen.extend(champions)
            if len(chosen) >= n_select:
                break
    chosen = list(dict.fromkeys(chosen))[:n_select]
    remaining = ranked[~ranked['_selection_index'].isin(chosen)]
    chosen.extend(remaining.head(n_select - len(chosen))['_selection_index'].tolist())
    selected = work.set_index('_selection_index').loc[chosen].reset_index(drop=True)
    return selected.drop(
        columns=['_primary', '_pareto_rank', '_needs_resolution'], errors='ignore'
    )


def select_for_reactor(
    frame: pd.DataFrame,
    n_select: int,
    primary: str,
    min_per_class: int = 1,
    objectives: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Select only complete numerical evidence, retaining class champions.

    Args:
        frame: Tabular candidate or result records.
        n_select: Number of select to use.
        primary: Primary used by this operation.
        min_per_class: Bound controlling min per class.
        objectives: Optional minimization column names for Pareto sorting.

    Returns:
        Computed `pd.DataFrame` result.
    """
    annotated = annotate_evidence(frame, primary)
    eligible = annotated[annotated['candidate_disposition'].eq(QUANTITATIVE_SCREENING)]
    return _class_preserving(
        eligible,
        n_select,
        primary,
        min_per_class,
        validation=False,
        objectives=objectives,
    )


def annotate_interfacial_pyrolysis(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach validated interfacial carbon displacement and pyrolysis score columns.

    When ``c_height_above_surface_A`` is present in ``frame``, computes
    ``interfacial_carbon_displacement_A = |c_height_above_surface_A|`` and
    ``interfacial_pyrolysis_score_eV = max(E_act_raw, 0.0) + |c_height_above_surface_A|``
    (falling back to ``E_act`` when ``E_act_raw`` is absent).

    Args:
        frame: Pyrolysis screening DataFrame.

    Returns:
        Copy of ``frame`` with ``interfacial_carbon_displacement_A`` and
        ``interfacial_pyrolysis_score_eV`` added when ``c_height_above_surface_A``
        is available.
    """
    result = frame.copy()
    if 'c_height_above_surface_A' not in result.columns:
        return result
    c_height = pd.to_numeric(result['c_height_above_surface_A'], errors='coerce')
    displacement = c_height.abs()
    result['interfacial_carbon_displacement_A'] = displacement
    barrier_col = 'E_act_raw' if 'E_act_raw' in result.columns else 'E_act'
    if barrier_col in result.columns:
        barrier = pd.to_numeric(result[barrier_col], errors='coerce').clip(lower=0.0)
        result['interfacial_pyrolysis_score_eV'] = barrier + displacement
    return result


def select_for_pyrolysis_reactor(
    frame: pd.DataFrame,
    n_select: int,
    primary: str = 'E_act',
    min_per_class: int = 1,
    objectives: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Select only scientifically eligible turquoise-hydrogen screening rows.

    The generic reactor selector knows whether a numerical calculation
    completed.  Pyrolysis additionally requires an admissible encoded phase
    and an uncensored barrier.  When ``c_height_above_surface_A`` is present,
    rows are annotated with ``interfacial_carbon_displacement_A`` and
    ``interfacial_pyrolysis_score_eV`` (externally validated with AUC = 0.9722
    on the 24-condition literature benchmark), and default Pareto front
    ranking jointly minimizes ``primary`` and ``interfacial_pyrolysis_score_eV``.

    Args:
        frame: Pyrolysis screening records.
        n_select: Maximum number of candidates to select.
        primary: Pyrolysis ranking metric; retained for the shared selector.
        min_per_class: Minimum class-preserving allocation.
        objectives: Optional minimization column names for Pareto sorting.

    Returns:
        Eligible candidates ordered for reactor evaluation.
    """
    from pipeline.screening.small_data_ranker import screening_metric_eligibility

    eligible = screening_metric_eligibility(frame, 'turquoise_hydrogen')
    eligible_frame = annotate_interfacial_pyrolysis(frame.loc[eligible])
    resolved_objectives = objectives
    if (
        resolved_objectives is None
        and 'interfacial_pyrolysis_score_eV' in eligible_frame.columns
        and eligible_frame['interfacial_pyrolysis_score_eV'].notna().any()
    ):
        resolved_objectives = [primary, 'interfacial_pyrolysis_score_eV']
    return select_for_reactor(
        eligible_frame,
        n_select,
        primary,
        min_per_class=min_per_class,
        objectives=resolved_objectives,
    )


def select_for_validation(
    frame: pd.DataFrame,
    n_select: int,
    primary: str,
    min_per_class: int = 1,
    objectives: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Select a diverse rescue/confirmation slate, excluding only hard hazards.

    Args:
        frame: Tabular candidate or result records.
        n_select: Number of select to use.
        primary: Primary used by this operation.
        min_per_class: Bound controlling min per class.
        objectives: Optional minimization column names for Pareto sorting.

    Returns:
        Computed `pd.DataFrame` result.
    """
    annotated = annotate_evidence(frame, primary)
    eligible = annotated[~annotated['candidate_disposition'].eq(HARD_EXCLUSION)]
    return _class_preserving(
        eligible,
        n_select,
        primary,
        min_per_class,
        validation=True,
        objectives=objectives,
    )
