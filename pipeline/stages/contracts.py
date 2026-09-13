"""Shared data contracts for independently runnable pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class StageOutcome:
    """Explicit persisted state and in-process products from one stage."""

    state: Mapping[str, Any]
    products: Mapping[str, Any] = field(default_factory=dict)
