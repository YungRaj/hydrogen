"""Persistence adapter for candidate handoff between discovery and later stages."""

from __future__ import annotations

from pathlib import Path


def load_selected_candidates(path: str | Path, *, top_k_reactor: int,
                             top_k_dft: int):
    """Load a screening table and reproduce the two distinct admission routes.

    Args:
        path: Filesystem path to the input or output artifact.
        top_k_reactor: Bound controlling top k reactor.
        top_k_dft: Bound controlling top k dft.

    Returns:
        Computed result described above.
    """
    import pandas as pd
    from pipeline.common.application_scope import scope_pyrolysis_pool
    from pipeline.screening.stage_selection import (
        select_for_reactor, select_for_validation)

    source = Path(path)
    if not source.is_file():
        return None
    database = pd.read_csv(source)
    # Same admissibility pool as the live discovery stage (ADR 0001), so a
    # phase-2-only restart draws its slates from the same rows. Both routes
    # stay distinct: the validation route may still rescue invalid rows.
    pool, admissibility = scope_pyrolysis_pool(database)
    return {
        'screening_database': database,
        'admissibility': admissibility,
        'top_catalysts': select_for_reactor(
            pool, top_k_reactor, 'E_act', min_per_class=1),
        'dft_candidates': select_for_validation(
            pool, top_k_dft, 'E_act', min_per_class=1),
    }
