"""Shared data contracts for independently runnable pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class StageOutcome:
    """Explicit persisted state and in-process products from one stage."""

    state: Mapping[str, Any]
    products: Mapping[str, Any] = field(default_factory=dict)


def require_stage_outcome(value, *, stage: str,
                          required_products: tuple[str, ...] = ()) -> StageOutcome:
    """Validate a replaceable stage at its orchestration boundary."""
    if not isinstance(value, StageOutcome):
        raise TypeError(f'{stage} stage must return StageOutcome')
    if not isinstance(value.state, Mapping) or not isinstance(value.products, Mapping):
        raise TypeError(f'{stage} StageOutcome mappings are invalid')
    missing = set(required_products).difference(value.products)
    if missing:
        raise ValueError(
            f'{stage} stage omitted required products: {sorted(missing)}')
    return value
