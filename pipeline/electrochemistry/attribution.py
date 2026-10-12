"""Mechanism-on/off attribution for native screening results."""

from __future__ import annotations

from typing import Any, Mapping


def compare_regeneration(
    baseline: Mapping[str, float], regeneration: Mapping[str, float]
) -> dict[str, Any]:
    """Return carbon-removal deltas without changing either result.

    Args:
        baseline: Reactor outputs without carbon regeneration.
        regeneration: Reactor outputs with carbon regeneration.
    """
    fields = (
        'CH4_conversion', 'H2_selectivity', 'solid_C_selectivity',
        'specific_energy_kWh_kg_H2',
    )
    missing = [field for field in fields if field not in baseline or field not in regeneration]
    if missing:
        raise ValueError(f'attribution results missing fields: {missing}')
    return {
        'baseline': dict(baseline),
        'regeneration': dict(regeneration),
        'delta': {field: float(regeneration[field]) - float(baseline[field]) for field in fields},
        'mechanism': 'carbon_regeneration',
    }
