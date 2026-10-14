"""Load deduplicated legacy outcomes used by historical benchmark evidence."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.search.discovery import candidate_id


@dataclass(frozen=True)
class PilotSpec:
    """Identify legacy outcome tables and their application metric."""

    application: str
    paths: tuple[str, ...]
    outcome: str


def _genome(value: object) -> tuple:
    genome = value if isinstance(value, tuple) else ast.literal_eval(str(value))
    if not isinstance(genome, tuple) or not genome:
        raise ValueError("invalid genome")
    return genome


def _identity(genome: tuple) -> str:
    # Canonical repository identity: representational equivalences (e.g.
    # SolidCatalyst dopant order) collapse to one candidate.
    return candidate_id(genome)


def load_legacy_outcomes(spec: PilotSpec) -> pd.DataFrame:
    """Load and candidate-deduplicate valid legacy outcomes.

    Args:
        spec: Input tables, application, and outcome column.

    Returns:
        Valid rows collapsed to one median outcome per catalyst genome.

    Raises:
        FileNotFoundError: If none of the input tables exist.
        ValueError: If required columns are missing or no valid, finite,
            parseable outcome rows remain.
    """
    frames = [pd.read_csv(path) for path in spec.paths if Path(path).is_file()]
    if not frames:
        raise FileNotFoundError(f"no pilot inputs found for {spec.application}")
    frame = pd.concat(frames, ignore_index=True)
    missing = {"genome", "valid", spec.outcome} - set(frame.columns)
    if missing:
        raise ValueError(f"missing pilot columns: {sorted(missing)}")
    frame = frame[frame.valid.eq(True)].copy()
    frame[spec.outcome] = pd.to_numeric(frame[spec.outcome], errors="coerce")
    frame = frame[np.isfinite(frame[spec.outcome])]
    parsed_rows: list[dict] = []
    for _, row in frame.iterrows():
        try:
            genome = _genome(row.genome)
            identity = _identity(genome)
        except (ValueError, SyntaxError, TypeError, IndexError):
            continue
        record = row.to_dict()
        record.update(parsed_genome=genome, candidate_id=identity)
        parsed_rows.append(record)
    if not parsed_rows:
        raise ValueError(
            f"no valid, finite, parseable {spec.outcome} rows for {spec.application}"
        )
    parsed = pd.DataFrame(parsed_rows)
    grouped = []
    for _, rows in parsed.groupby("candidate_id", sort=True):
        first = rows.iloc[0].copy()
        first[spec.outcome] = float(rows[spec.outcome].median())
        first["replicates"] = int(len(rows))
        grouped.append(first)
    return pd.DataFrame(grouped).reset_index(drop=True)


