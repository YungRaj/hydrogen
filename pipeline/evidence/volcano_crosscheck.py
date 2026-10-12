"""Descriptor-only ORR volcano cross-check.

This is intentionally independent of the screening model.  It provides the
standard one-dimensional OH binding-energy volcano and a CSV export that can
be consumed by CatMAP when that optional dependency is available.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, Mapping


def volcano_rate_log10(dG_OH_eV: float, *, temperature_K: float = 298.15) -> float:
    """Return a dimensionless log10 relative rate from a simple OH volcano.

    The optimum is 0.86 eV (the conventional ORR OH descriptor); deviations
    are penalized symmetrically as a descriptor sanity check, not as a fitted
    screening score.
    """
    if temperature_K <= 0 or not math.isfinite(dG_OH_eV):
        raise ValueError('finite descriptor and positive temperature required')
    activation_eV = abs(float(dG_OH_eV) - 0.86)
    return -activation_eV / (8.617333262e-5 * temperature_K * math.log(10.0))


def rank_descriptors(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    """Annotate and rank rows without changing their source values."""
    ranked = []
    for row in rows:
        value = float(row['dG_OH_eV'])
        ranked.append({**row, 'volcano_log10_rate': volcano_rate_log10(value)})
    return sorted(ranked, key=lambda row: float(row['volcano_log10_rate']), reverse=True)


def write_catmap_input(rows: Iterable[Mapping[str, object]], path: str | Path) -> None:
    """Write the minimal descriptor table accepted by a CatMAP preprocessor."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=('candidate_id', 'dG_OH_eV'))
        writer.writeheader()
        for row in rows:
            writer.writerow({
                'candidate_id': str(row.get('candidate_id', 'unknown')),
                'dG_OH_eV': float(row['dG_OH_eV']),
            })
