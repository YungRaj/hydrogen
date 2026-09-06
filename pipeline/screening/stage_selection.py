"""Evidence-aware candidate selection at fidelity boundaries.

An incomplete calculation is not evidence that a material is poor.  Reactor
models require complete numerical descriptors, while high-fidelity validation
must retain a route for unresolved and unfamiliar chemistry.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


HARD_EXCLUSION = 'hard_excluded'
VALIDATION_REQUIRED = 'validation_required'
QUANTITATIVE_SCREENING = 'quantitative_screening'


def evidence_disposition(row, primary: str) -> str:
    """Classify evidence without conflating calculation failure with chemistry."""
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
    result = frame.copy()
    result['candidate_disposition'] = [
        evidence_disposition(row, primary) for _, row in result.iterrows()]
    return result


def _class_preserving(frame: pd.DataFrame, n_select: int, primary: str,
                      min_per_class: int, validation: bool) -> pd.DataFrame:
    if n_select <= 0 or frame.empty:
        return frame.iloc[0:0].copy()
    work = frame.copy()
    work['_selection_index'] = np.arange(len(work))
    values = pd.to_numeric(work.get(primary), errors='coerce')
    work['_primary'] = values.fillna(np.inf)
    needs_dft = (work['needs_dft_validation'].map(
                    lambda value: False if pd.isna(value) else bool(value))
                 if 'needs_dft_validation' in work else
                 pd.Series(False, index=work.index))
    work['_needs_resolution'] = (
        work['candidate_disposition'].eq(VALIDATION_REQUIRED) |
        needs_dft).astype(int)
    # Validation gives unresolved evidence first access; quantitative consumers
    # can only use complete screening rows. Stable index makes ties reproducible.
    if validation:
        order = ['_needs_resolution', '_primary', '_selection_index']
        ascending = [False, True, True]
    else:
        order = ['_primary', '_selection_index']
        ascending = [True, True]
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
            index: position for position, index in
            enumerate(ranked['_selection_index'].tolist())}
        for quota_round in range(min_per_class):
            champions = [
                group.iloc[quota_round]['_selection_index']
                for group in class_rows.values() if len(group) > quota_round]
            champions.sort(key=lambda index: rank_position[index])
            chosen.extend(champions)
            if len(chosen) >= n_select:
                break
    chosen = list(dict.fromkeys(chosen))[:n_select]
    remaining = ranked[~ranked['_selection_index'].isin(chosen)]
    chosen.extend(remaining.head(n_select - len(chosen))['_selection_index'].tolist())
    selected = work.set_index('_selection_index').loc[chosen].reset_index(drop=True)
    return selected.drop(columns=['_primary', '_needs_resolution'], errors='ignore')


def select_for_reactor(frame: pd.DataFrame, n_select: int, primary: str,
                       min_per_class: int = 1) -> pd.DataFrame:
    """Select only complete numerical evidence, retaining class champions."""
    annotated = annotate_evidence(frame, primary)
    eligible = annotated[
        annotated['candidate_disposition'].eq(QUANTITATIVE_SCREENING)]
    return _class_preserving(
        eligible, n_select, primary, min_per_class, validation=False)


def select_for_validation(frame: pd.DataFrame, n_select: int, primary: str,
                          min_per_class: int = 1) -> pd.DataFrame:
    """Select a diverse rescue/confirmation slate, excluding only hard hazards."""
    annotated = annotate_evidence(frame, primary)
    eligible = annotated[~annotated['candidate_disposition'].eq(HARD_EXCLUSION)]
    return _class_preserving(
        eligible, n_select, primary, min_per_class, validation=True)
