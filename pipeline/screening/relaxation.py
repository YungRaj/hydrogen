"""Fail-closed, auditable geometry relaxation helpers."""

from __future__ import annotations

import hashlib
import json

import numpy as np
from ase.optimize import BFGS, FIRE


def geometry_digest(atoms) -> str:
    """Stable digest of the final chemical identity, cell, PBC, and positions."""
    digest = hashlib.sha256()
    digest.update(np.asarray(atoms.numbers, dtype=np.int16).tobytes())
    digest.update(np.asarray(atoms.positions, dtype=np.float64).tobytes())
    digest.update(np.asarray(atoms.cell, dtype=np.float64).tobytes())
    digest.update(np.asarray(atoms.pbc, dtype=np.bool_).tobytes())
    return digest.hexdigest()


def _max_force(atoms) -> float:
    forces = np.asarray(atoms.get_forces(), dtype=float)
    return float(np.linalg.norm(forces, axis=1).max()) if len(forces) else 0.0


def _initial_geometry_error(atoms) -> str | None:
    if len(atoms) == 0:
        return 'empty_structure'
    if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell).all():
        return 'non_finite_geometry'
    if np.any(atoms.pbc) and (not np.all(atoms.pbc) or abs(atoms.cell.volume) < 1e-8):
        return 'invalid_periodic_cell'
    if len(atoms) > 1:
        distances = np.asarray(atoms.get_all_distances(mic=bool(np.all(atoms.pbc))))
        distances[np.eye(len(atoms), dtype=bool)] = np.inf
        if float(distances.min()) < 0.35:
            return 'atomic_overlap'
    return None


def _attempt(optimizer, atoms, fmax: float, steps: int, name: str) -> dict:
    started = geometry_digest(atoms)
    try:
        ase_converged = bool(optimizer.run(fmax=fmax, steps=steps))
        max_force = _max_force(atoms)
        finite = np.isfinite(max_force) and np.isfinite(atoms.positions).all()
        converged = ase_converged and finite and max_force <= fmax
        return {'name': name, 'converged': bool(converged),
                'ase_converged': ase_converged, 'steps': int(optimizer.nsteps),
                'max_force_eV_A': max_force, 'start_geometry_sha256': started,
                'end_geometry_sha256': geometry_digest(atoms),
                'failure_class': None if converged else (
                    'numerical_failure' if not finite else 'force_threshold')}
    except Exception as exc:
        return {'name': name, 'converged': False,
                'ase_converged': False, 'steps': int(getattr(optimizer, 'nsteps', 0)),
                'max_force_eV_A': None, 'start_geometry_sha256': started,
                'end_geometry_sha256': geometry_digest(atoms),
                'failure_class': 'optimizer_exception',
                'error': f'{type(exc).__name__}: {str(exc)[:300]}'}


def relax_with_record(atoms, label: str, fmax: float, steps: int,
                      recovery: bool = True) -> dict:
    """Relax atoms and return explicit convergence evidence.

    ASE's return value is intentionally checked; exhausting the step budget is
    not accepted as convergence.  Forces are evaluated once at the final
    geometry so the recorded maximum force is the quantity behind the gate.
    """
    prefix = f'relax_{label}'
    geometry_error = _initial_geometry_error(atoms)
    if geometry_error:
        return {
            f'{prefix}_converged': False,
            f'{prefix}_max_force_eV_A': np.nan,
            f'{prefix}_steps': 0,
            f'{prefix}_step_limit': int(steps),
            f'{prefix}_fmax_eV_A': float(fmax),
            f'{prefix}_termination': 'invalid_initial_geometry',
            f'{prefix}_failure_class': geometry_error,
            f'{prefix}_selected_attempt': 'none',
            f'{prefix}_attempts_json': '[]',
            f'{prefix}_geometry_sha256': geometry_digest(atoms),
        }

    initial_positions = atoms.positions.copy()
    attempts = []
    best_positions = initial_positions.copy()
    best_force = np.inf

    def retain(attempt):
        nonlocal best_force, best_positions
        attempts.append(attempt)
        force = attempt.get('max_force_eV_A')
        if force is not None and np.isfinite(force) and force < best_force:
            best_force = force
            best_positions = atoms.positions.copy()

    retain(_attempt(BFGS(atoms, logfile=None, maxstep=0.20), atoms,
                    fmax, steps, 'bfgs_standard'))

    if recovery and not attempts[-1]['converged']:
        atoms.set_positions(initial_positions)
        precondition = _attempt(
            FIRE(atoms, logfile=None, dt=0.05, dtmax=0.5, maxstep=0.05),
            atoms, max(0.25, 3.0 * fmax), max(20, steps // 2),
            'fire_precondition')
        # This stage uses a deliberately loose force target and only prepares
        # the geometry; it can never satisfy the final evidence gate itself.
        precondition['qualifying'] = False
        retain(precondition)
        retain(_attempt(BFGS(atoms, logfile=None, maxstep=0.08), atoms,
                        fmax, 2 * steps, 'bfgs_after_fire'))

    if recovery and not attempts[-1]['converged']:
        atoms.set_positions(best_positions)
        retain(_attempt(
            FIRE(atoms, logfile=None, dt=0.025, dtmax=0.25, maxstep=0.03,
                 downhill_check=True), atoms, fmax, 3 * steps,
            'fire_small_step'))

    successful = next((attempt for attempt in attempts
                       if attempt['converged'] and attempt.get('qualifying', True)), None)
    if successful is None:
        atoms.set_positions(best_positions)
        selected = min(
            (attempt for attempt in attempts
             if attempt.get('max_force_eV_A') is not None and
             np.isfinite(attempt['max_force_eV_A'])),
            key=lambda attempt: attempt['max_force_eV_A'], default=attempts[-1])
    else:
        selected = successful
    max_force = selected.get('max_force_eV_A')
    converged = successful is not None
    total_steps = sum(attempt['steps'] for attempt in attempts)
    failure_class = None if converged else selected.get('failure_class', 'force_threshold')
    return {
        f'{prefix}_converged': converged,
        f'{prefix}_max_force_eV_A': max_force,
        f'{prefix}_steps': total_steps,
        f'{prefix}_step_limit': int(steps),
        f'{prefix}_fmax_eV_A': float(fmax),
        f'{prefix}_termination': 'converged' if converged else 'recovery_exhausted',
        f'{prefix}_failure_class': failure_class,
        f'{prefix}_selected_attempt': selected['name'],
        f'{prefix}_attempts_json': json.dumps(attempts, sort_keys=True),
        f'{prefix}_geometry_sha256': geometry_digest(atoms),
    }


def require_relaxation(result: dict, atoms, label: str,
                       fmax: float, steps: int, recovery: bool = True) -> bool:
    """Add a relaxation record and fail the candidate closed when incomplete."""
    record = relax_with_record(atoms, label, fmax, steps, recovery=recovery)
    result.update(record)
    if record[f'relax_{label}_converged']:
        return True
    result['valid'] = False
    if record[f'relax_{label}_termination'] == 'invalid_initial_geometry':
        result['error'] = (
            f"Invalid {label} initial geometry: "
            f"{record[f'relax_{label}_failure_class']}")
    else:
        result['error'] = (
            f"Unconverged {label} relaxation: max_force="
            f"{record[f'relax_{label}_max_force_eV_A']:.6g} eV/A after "
            f"{record[f'relax_{label}_steps']} cumulative recovery steps")
    result['needs_dft_validation'] = True
    result['candidate_disposition'] = 'validation_required'
    return False
