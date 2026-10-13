"""Versioned supported-nanoparticle screening adapter."""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterator

from ase.constraints import FixBondLength
from pipeline.screening import surface_screener

SUPPORTED_NANOPARTICLE_PROTOCOL = 'supported-nanoparticle-v1'
SUPPORTED_ADSORBATE_LIFTS_A = (0.8, 1.2, 1.6)
_SUPPORTED_CLASSES = frozenset({'SupportedNanoparticle'})


@contextmanager
def _lift_supported_adsorption_site(lift_A: float) -> Iterator[None]:
    """Apply the supported-particle placement correction during evaluation."""
    original = surface_screener.screening_adsorption_site
    original_relax = surface_screener.require_relaxation

    def corrected(structure: Any, active: list[int], material: str, distance: float):
        position, direction = original(structure, active, material, distance)
        if material in _SUPPORTED_CLASSES:
            position = position + lift_A * direction
        return position, direction

    surface_screener.screening_adsorption_site = corrected
    def constrained_relax(result: dict[str, Any], atoms: Any, label: str,
                          fmax: float, steps: int) -> bool:
        if label == 'ch3' and len(atoms) >= 4:
            start = len(atoms) - 3
            atoms.set_constraint(list(atoms.constraints) + [
                FixBondLength(len(atoms) - 4, index)
                for index in range(start, len(atoms))
            ])
        return original_relax(result, atoms, label, fmax, steps)

    surface_screener.require_relaxation = constrained_relax
    try:
        yield
    finally:
        surface_screener.screening_adsorption_site = original
        surface_screener.require_relaxation = original_relax


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
    result: dict[str, Any] = {}
    for lift_A in SUPPORTED_ADSORBATE_LIFTS_A:
        with _lift_supported_adsorption_site(lift_A):
            result = surface_screener.evaluate_candidate(
                genome, calculator, references, structure_seed=structure_seed
            )
        if result.get('valid') or result.get('error') != 'CH3 adsorbate lost molecular integrity':
            break
    result['screening_protocol'] = SUPPORTED_NANOPARTICLE_PROTOCOL
    return result
