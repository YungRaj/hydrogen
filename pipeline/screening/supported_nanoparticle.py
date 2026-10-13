"""Versioned supported-nanoparticle screening adapter."""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterator

from pipeline.screening import surface_screener

SUPPORTED_NANOPARTICLE_PROTOCOL = 'supported-nanoparticle-v1'
_SUPPORTED_CLASSES = frozenset({'SupportedNanoparticle'})


@contextmanager
def _lift_supported_adsorption_site() -> Iterator[None]:
    """Apply the supported-particle placement correction during evaluation."""
    original = surface_screener.screening_adsorption_site

    def corrected(structure: Any, active: list[int], material: str, distance: float):
        position, direction = original(structure, active, material, distance)
        if material in _SUPPORTED_CLASSES:
            position = position + 0.35 * direction
        return position, direction

    surface_screener.screening_adsorption_site = corrected
    try:
        yield
    finally:
        surface_screener.screening_adsorption_site = original


def evaluate_supported_nanoparticle(
    genome: tuple[Any, ...], calculator: Any, references: dict[str, float],
    *, structure_seed: int = 0,
) -> dict[str, Any]:
    """Evaluate one supported nanoparticle with the v1 placement protocol.

    Args:
        genome: SupportedNanoparticle genome tuple.
        calculator: Atomic calculator used by the shared screener.
        references: Gas-phase reference energies.
        structure_seed: Deterministic structure replica index.

    Returns:
        Screening result tagged with the supported-nanoparticle protocol.
    """
    if not genome or genome[0] not in _SUPPORTED_CLASSES:
        raise ValueError('supported-nanoparticle adapter received another class')
    with _lift_supported_adsorption_site():
        result = surface_screener.evaluate_candidate(
            genome, calculator, references, structure_seed=structure_seed
        )
    result['screening_protocol'] = SUPPORTED_NANOPARTICLE_PROTOCOL
    return result
