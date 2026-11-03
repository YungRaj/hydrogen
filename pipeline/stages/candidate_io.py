"""Persistence adapter for candidate handoff between discovery and later stages."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

from pipeline.data_models.stages import SelectedCandidates, TableLike

# Screening evidence persisted by the branch discovery workflow
# (BranchDiscoveryWorkflow.evidence_db for turquoise hydrogen). The retired
# genetic search wrote ``ga_full_database.csv``; nothing writes that any more.
SCREENING_EVIDENCE_DB = 'branch_ranker_evidence.csv'


def load_reactor_reference(candidate_id: str) -> Any | None:
    """Load a tracked literature reference for the reactor batch only.

    Reference rows are separate from discovery and DFT candidate selection.
    Unknown names may still identify a judge already supplied by the caller.

    Args:
        candidate_id: Input controlling candidate id.

    Returns:
        The matching reference screening row, or ``None`` if not tracked.
    """
    if candidate_id != 'ni_np_lit':
        return None
    import pandas as pd  # pyright: ignore[reportUnknownVariableType]

    from pipeline.utils import BASE_DIR

    path = BASE_DIR / 'sweeps' / 'ni_np_lit_screening_row.csv'
    table: Any = pd.read_csv(path)  # pyright: ignore[reportUnknownMemberType]
    matches: Any = table.loc[table['candidate_id'] == candidate_id]
    if len(matches) != 1:
        raise ValueError(f'expected one {candidate_id} reference row in {path}')
    return matches.iloc[0]


def load_selected_candidates(
    path: str | Path, *, top_k_reactor: int, top_k_dft: int
) -> SelectedCandidates | None:
    """Load a screening table and reproduce the two distinct admission routes.

    Args:
        path: Filesystem path to the input or output artifact.
        top_k_reactor: Bound controlling top k reactor.
        top_k_dft: Bound controlling top k dft.

    Returns:
        Computed result described above.
    """
    import pandas as pd  # pyright: ignore[reportUnknownVariableType]

    from pipeline.screening.stage_selection import (
        select_for_pyrolysis_reactor,  # pyright: ignore[reportUnknownVariableType]
        select_for_validation,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.search.discovery import (
        add_discovery_metadata,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.search.scope import scope_pyrolysis_pool

    source = Path(path)
    if not source.is_file():
        return None
    # The persisted evidence table is written before discovery adds the
    # content-addressed candidate_id; restore it so restarted reactor/DFT/VQE
    # stages key artifacts exactly as the live run does.
    database: Any = pd.read_csv(source)  # pyright: ignore[reportUnknownMemberType]
    if 'candidate_id' not in database.columns:
        database = cast(Any, add_discovery_metadata(database))
    # Same admissibility pool as the live discovery stage (ADR 0001), so a
    # phase-2-only restart draws its slates from the same rows. Both routes
    # stay distinct: the validation route may still rescue invalid rows.
    pool, raw_admissibility = scope_pyrolysis_pool(database)
    reactor_slate = cast(
        TableLike,
        select_for_pyrolysis_reactor(pool, top_k_reactor, 'E_act', min_per_class=1),
    )
    dft_slate = cast(
        TableLike,
        select_for_validation(pool, top_k_dft, 'E_act', min_per_class=1),
    )
    selected: SelectedCandidates = {
        'screening_database': cast(TableLike, database),
        'admissibility': raw_admissibility,
        'top_catalysts': reactor_slate,
        'dft_candidates': dft_slate,
    }
    return selected
