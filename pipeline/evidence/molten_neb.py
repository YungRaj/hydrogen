"""Fail-closed eSEN CI-NEB validation for molten-metal methane activation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from ase import Atom, Atoms
from ase.build import molecule
from ase.mep import NEB
from ase.optimize import FIRE

from pipeline.screening.relaxation import geometry_digest, relax_with_record
from pipeline.screening.surface_calculator import get_ocp_calculator
from pipeline.screening.surface_screener import generate_structure


MODEL_NAME = 'esen-sm-conserving-all-oc25'


@dataclass(frozen=True)
class NEBConfig:
    """Numerical controls frozen before external barrier evaluation."""

    images: int = 7
    endpoint_fmax_eV_A: float = 0.08
    endpoint_steps: int = 200
    path_fmax_eV_A: float = 0.08
    path_steps: int = 400
    spring_eV_A2: float = 0.10
    max_image_displacement_A: float = 1.5


def _active_site(structure: Atoms, indices: list[int], promoter: str) -> int:
    """Choose a promoter surface atom, or the most central host site."""
    candidates = [index for index in indices if 0 <= index < len(structure)]
    promoted = [index for index in candidates if structure[index].symbol == promoter]
    eligible = promoted or candidates
    if not eligible:
        raise ValueError('molten NEB requires a surface active site')
    center = 0.5 * (structure.cell[0] + structure.cell[1])
    return int(min(
        eligible,
        key=lambda index: float(
            np.linalg.norm((structure[index].position - center)[:2])
        ),
    ))


def methane_endpoints(
    structure: Atoms, active_index: int, calculator: Any, config: NEBConfig
) -> tuple[Atoms, Atoms, dict[str, Any]]:
    """Build cleanly ordered CH4 and CH3*+H* endpoints and relax the product.

    Args:
        structure: Relaxed clean catalyst surface.
        active_index: Surface atom on which methane is activated.
        calculator: Shared energy-conserving MLIP calculator.
        config: Frozen endpoint relaxation controls.

    Returns:
        Initial state, relaxed product state, and product relaxation evidence.
    """
    site = structure[active_index].position.copy()
    initial = structure.copy()
    methane = molecule('CH4')
    methane.translate(site + np.array([0.0, 0.0, 3.5]) - methane[0].position)
    initial += methane

    final = structure.copy()
    carbon = site + np.array([0.0, 0.0, 2.0])
    final.append(Atom('C', carbon))
    final.extend(Atoms('H3', positions=(
        carbon + np.array([0.0, 1.02, 0.35]),
        carbon + np.array([-0.88, -0.51, 0.35]),
        carbon + np.array([0.88, -0.51, 0.35]),
    )))
    final.append(Atom('H', site + np.array([1.5, 0.0, 1.2])))
    for endpoint in (initial, final):
        endpoint.set_cell(structure.cell)
        endpoint.set_pbc(structure.pbc)
        endpoint.set_constraint(structure.constraints)
        endpoint.calc = calculator
    relaxation = relax_with_record(
        final,
        'neb_product',
        config.endpoint_fmax_eV_A,
        config.endpoint_steps,
    )
    if not relaxation['relax_neb_product_converged']:
        raise RuntimeError('MLIP-NEB product endpoint did not converge')
    return initial, final, relaxation


def _path_continuous(images: list[Atoms], maximum_A: float) -> bool:
    return all(
        float(np.linalg.norm(right.positions - left.positions, axis=1).max())
        <= maximum_A
        for left, right in zip(images, images[1:])
    )


def run_replica(
    genome: tuple[Any, ...], structure_seed: int, device: str, config: NEBConfig
) -> dict[str, Any]:
    """Run one climbing-image path and return explicit convergence evidence.

    Args:
        genome: Molten-metal catalyst encoding.
        structure_seed: Deterministic surface configuration index.
        device: Calculator device selected by the caller.
        config: Frozen endpoint and path controls.

    Returns:
        Barrier, endpoint, convergence, continuity, and provenance evidence.
    """
    calculator = get_ocp_calculator(MODEL_NAME, device=device)
    if calculator is None:
        raise RuntimeError(f'could not load {MODEL_NAME}')
    structure, active_indices, material_class = generate_structure(genome, structure_seed)
    if material_class != 'MoltenMetal':
        raise ValueError('MLIP-NEB validation accepts only MoltenMetal genomes')
    structure.pbc = True
    structure.calc = calculator
    clean = relax_with_record(
        structure,
        'neb_clean',
        config.endpoint_fmax_eV_A,
        config.endpoint_steps,
    )
    if not clean['relax_neb_clean_converged']:
        raise RuntimeError('MLIP-NEB clean endpoint did not converge')
    active_index = _active_site(structure, active_indices, str(genome[2]))
    initial, final, product = methane_endpoints(
        structure, active_index, calculator, config
    )
    images = [initial] + [initial.copy() for _ in range(config.images - 2)] + [final]
    neb = NEB(
        images,
        k=config.spring_eV_A2,
        climb=True,
        method='improvedtangent',
        allow_shared_calculator=True,
    )
    neb.interpolate(method='idpp', apply_constraint=True)
    for image in images:
        image.calc = calculator
    optimizer = FIRE(neb, logfile=None, dt=0.03, dtmax=0.2, maxstep=0.05)
    converged = bool(
        optimizer.run(fmax=config.path_fmax_eV_A, steps=config.path_steps)
    )
    energies = [float(image.get_potential_energy()) for image in images]
    forces = np.asarray(neb.get_forces(), dtype=float)
    max_force = float(np.linalg.norm(forces, axis=1).max())
    continuous = _path_continuous(images, config.max_image_displacement_A)
    valid = bool(
        converged
        and math.isfinite(max_force)
        and max_force <= config.path_fmax_eV_A
        and continuous
        and max(energies) not in (energies[0], energies[-1])
    )
    return {
        'genome': repr(genome),
        'structure_seed': structure_seed,
        'model': MODEL_NAME,
        'config': asdict(config),
        'active_index': active_index,
        'active_element': structure[active_index].symbol,
        'clean_geometry_sha256': geometry_digest(structure),
        'initial_geometry_sha256': geometry_digest(initial),
        'final_geometry_sha256': geometry_digest(final),
        'endpoint_relaxation': {'clean': clean, 'product': product},
        'energies_eV': energies,
        'relative_energies_eV': [value - energies[0] for value in energies],
        'barrier_eV': max(energies) - energies[0],
        'reaction_energy_eV': energies[-1] - energies[0],
        'path_max_force_eV_A': max_force,
        'path_steps': int(optimizer.nsteps),
        'path_continuous': continuous,
        'valid': valid,
    }


def lock_digest(path: str | Path) -> str:
    """Return the SHA-256 identity of a preregistration lock.

    Args:
        path: Preregistration lock file.

    Returns:
        Lowercase SHA-256 digest of the exact lock bytes.
    """
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main() -> None:
    """Run locked replicas without reading held-out experimental labels."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lock', default='docs/evidence/molten_neb_validation_lock.json')
    parser.add_argument('--output', default='results/evidence/molten_neb_validation.json')
    parser.add_argument('--device', default='cuda:0')
    args = parser.parse_args()
    lock_path = Path(args.lock)
    lock = json.loads(lock_path.read_text())
    implementation = Path(__file__)
    actual_implementation = hashlib.sha256(implementation.read_bytes()).hexdigest()
    if actual_implementation != lock.get('implementation_sha256'):
        parser.error('molten NEB implementation does not match the preregistration lock')
    config = NEBConfig(**lock['neb_config'])
    records = [
        run_replica(tuple(case['genome']), seed, args.device, config)
        for case in lock['external_cases']
        for seed in case['structure_seeds']
    ]
    valid = [record for record in records if record['valid']]
    references = [
        float(value)
        for case in lock['external_cases']
        for value in case['reference_barriers_eV']
    ]
    predicted = [record['barrier_eV'] for record in valid]
    minimum_valid = int(lock['acceptance']['minimum_valid_paths'])
    median_error = (
        abs(float(np.median(predicted)) - float(np.median(references)))
        if predicted else None
    )
    passed = bool(
        len(valid) >= minimum_valid
        and median_error is not None
        and median_error <= float(lock['acceptance']['maximum_median_absolute_error_eV'])
    )
    report = {
        'schema_version': 1,
        'lock_sha256': lock_digest(lock_path),
        'valid_paths': len(valid),
        'attempted_paths': len(records),
        'median_predicted_barrier_eV': float(np.median(predicted)) if predicted else None,
        'median_reference_barrier_eV': float(np.median(references)),
        'median_absolute_error_eV': median_error,
        'acceptance_passed': passed,
        'held_out_conditions_consumed': False,
        'records': records,
    }
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(target)
    print(json.dumps({key: value for key, value in report.items() if key != 'records'}, indent=2))


if __name__ == '__main__':
    main()
