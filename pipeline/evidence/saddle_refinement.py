"""Optional saddle-point refinement for already-built NEB transition states."""

from __future__ import annotations

from typing import Any

from ase import Atoms


def refine_saddle(
    atoms: Atoms,
    *,
    fmax_eV_A: float = 0.05,
    steps: int = 500,
) -> dict[str, Any]:
    """Refine one NEB maximum with Sella and return an auditable summary.

    Args:
        atoms: Transition-state candidate with an attached calculator.
        fmax_eV_A: Maximum force convergence threshold in eV/Å.
        steps: Maximum Sella optimization steps.

    Returns:
        Convergence, force, energy, and step metadata.

    Sella is deliberately optional: importing this function does not require
    the package, while invoking it fails with a clear installation error.
    """
    if atoms.calc is None:
        raise ValueError('saddle refinement requires an attached calculator')
    if fmax_eV_A <= 0.0 or steps <= 0:
        raise ValueError('fmax_eV_A and steps must be positive')
    try:
        from sella import Sella
    except ImportError as exc:
        raise RuntimeError(
            'Sella is required for saddle refinement; install the atomistic '
            'validation dependencies'
        ) from exc

    optimizer = Sella(atoms, order=1)
    optimizer.run(fmax=fmax_eV_A, steps=steps)
    return {
        'converged': bool(optimizer.converged()),
        'steps': int(getattr(optimizer, 'nsteps', steps)),
        'energy_eV': float(atoms.get_potential_energy()),
        'max_force_eV_A': float(max(abs(value) for value in atoms.get_forces().ravel())),
    }
