"""Shared data contracts for independently runnable pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Generic, Mapping, TypeVar, cast


StateT = TypeVar("StateT", bound=Mapping[str, Any])
ProductsT = TypeVar("ProductsT", bound=Mapping[str, Any])


@dataclass(frozen=True, slots=True)
class StageOutcome(Generic[StateT, ProductsT]):
    """Typed persisted state and in-process products from one stage.

    The generic parameters retain each stage's concrete mapping shape for
    static analysis.  Runtime values remain mappings to preserve existing
    artifacts, tests, and third-party stage implementations.
    """

    state: StateT
    products: ProductsT


def require_stage_outcome(
        value: StageOutcome[StateT, ProductsT], *, stage: str,
        required_products: tuple[str, ...] = (),
        ) -> StageOutcome[StateT, ProductsT]:
    """Validate a replaceable stage at its orchestration boundary.

    Args:
        value: Value used by this operation.
        stage: Stage used by this operation.
        required_products: Ordered values supplying required products.

    Returns:
        Computed `StageOutcome` result.
    """
    # Replacements can violate the annotation at runtime, so retain this
    # guard even though static callers already see the generic contract.
    if not isinstance(value, StageOutcome):  # pyright: ignore[reportUnnecessaryIsInstance]
        raise TypeError(f'{stage} stage must return StageOutcome')
    unchecked = cast(StageOutcome[Any, Any], value)
    if (not isinstance(unchecked.state, Mapping) or
            not isinstance(unchecked.products, Mapping)):
        raise TypeError(f'{stage} StageOutcome mappings are invalid')
    outcome = cast(StageOutcome[StateT, ProductsT], unchecked)
    missing = set(required_products).difference(outcome.products)
    if missing:
        raise ValueError(
            f'{stage} stage omitted required products: {sorted(missing)}')
    return outcome
