"""Load deduplicated legacy outcomes used by historical benchmark evidence."""

from __future__ import annotations

import ast
import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PilotSpec:
    """Identify legacy outcome tables and their application metric."""

    application: str
    paths: tuple[str, ...]
    outcome: str
    folds: int = 5
    selection_fraction: float = 0.20
    random_trials: int = 20_000


def _genome(value: object) -> tuple:
    genome = value if isinstance(value, tuple) else ast.literal_eval(str(value))
    if not isinstance(genome, tuple) or not genome:
        raise ValueError("invalid genome")
    return genome


def _identity(genome: tuple) -> str:
    return hashlib.sha256(repr(genome).encode()).hexdigest()


def load_legacy_outcomes(spec: PilotSpec) -> pd.DataFrame:
    """Load and candidate-deduplicate valid legacy outcomes.

    Args:
        spec: Input tables, application, and outcome column.

    Returns:
        Valid rows collapsed to one median outcome per catalyst genome.
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
        except (ValueError, SyntaxError):
            continue
        record = row.to_dict()
        record.update(parsed_genome=genome, candidate_id=_identity(genome))
        parsed_rows.append(record)
    parsed = pd.DataFrame(parsed_rows)
    grouped = []
    for _, rows in parsed.groupby("candidate_id", sort=True):
        first = rows.iloc[0].copy()
        first[spec.outcome] = float(rows[spec.outcome].median())
        first["replicates"] = int(len(rows))
        grouped.append(first)
    return pd.DataFrame(grouped).reset_index(drop=True)


def default_specs() -> tuple[PilotSpec, PilotSpec]:
    """Return the legacy input tables required by prospective campaign training.

    Returns:
        Turquoise-hydrogen and fuel-cell legacy outcome specifications.
    """
    return (
        PilotSpec(
            "turquoise_hydrogen",
            (
                "results/screening/ga_initial_screening.csv",
                "results/screening/ga_fairchem_gen1.csv",
                "results/screening/ga_fairchem_gen2.csv",
                "results/screening/ga_mace_gen1.csv",
                "results/screening/ga_mace_gen2.csv",
            ),
            "E_act",
        ),
        PilotSpec(
            "fuel_cell_orr",
            (
                "results/fuel_cell/cathode_screening.csv",
                "results/fuel_cell/fc_initial_screening.csv",
                "results/fuel_cell/fc_fairchem_gen1.csv",
                "results/fuel_cell/fc_mace_gen1.csv",
            ),
            "orr_overpotential_V",
        ),
    )
