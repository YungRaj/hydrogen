"""Fail-closed staged eSEN NEB validation for methane activation."""

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
from ase.constraints import FixAtoms, Hookean
from ase.geometry import find_mic
from ase.mep import NEB
from ase.optimize import FIRE

from pipeline.screening.relaxation import geometry_digest, relax_with_record
from pipeline.screening.surface_calculator import get_ocp_calculator
from pipeline.screening.surface_screener import generate_structure

MODEL_NAME = 'esen-sm-conserving-all-oc25'


@dataclass(frozen=True)
class NEBConfig:
    """Frozen endpoint and staged-path numerical controls."""

    images: int = 11
    endpoint_fmax_eV_A: float = 0.08
    endpoint_steps: int = 400
    path_fmax_eV_A: float = 0.08
    plain_path_steps: int = 1000
    climbing_path_steps: int = 1000
    climbing_maxstep: float = 0.02
    spring_eV_A2: float = 0.70
    max_image_displacement_A: float = 1.5
    local_radius_A: float = 4.5
    reactant_height_A: float = 2.30
    product_CH_restraint_A: float = 2.25
    maximum_product_CH_A: float = 2.5


def _active_site(structure: Atoms, indices: list[int], promoter: str) -> int:
    candidates = [index for index in indices if 0 <= index < len(structure)]
    promoted = [index for index in candidates if structure[index].symbol == promoter]
    eligible = promoted or candidates
    if not eligible:
        raise ValueError('molten NEB requires a surface active site')
    center = 0.5 * (structure.cell[0] + structure.cell[1])
    return int(min(eligible, key=lambda index: float(
        np.linalg.norm((structure[index].position - center)[:2])
    )))


def _local_constraint(
    structure: Atoms, active_index: int, local_radius_A: float
) -> FixAtoms:
    fixed = {
        int(index)
        for constraint in structure.constraints
        if hasattr(constraint, 'get_indices')
        for index in constraint.get_indices()
    }
    fixed.update(
        index for index in range(len(structure))
        if structure.get_distance(active_index, index, mic=True) > local_radius_A
    )
    return FixAtoms(indices=sorted(fixed))


def _tetrahedral_vectors(direction: np.ndarray, bond_A: float) -> list[np.ndarray]:
    axis = np.asarray(direction, dtype=float)
    norm = float(np.linalg.norm(axis))
    if norm == 0.0:
        raise ValueError('breaking-bond direction must be nonzero')
    axis /= norm
    trial = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(axis, trial))) > 0.9:
        trial = np.array([1.0, 0.0, 0.0])
    first = np.cross(axis, trial)
    first /= np.linalg.norm(first)
    second = np.cross(axis, first)
    radial = math.sqrt(8.0 / 9.0)
    others = [
        -axis / 3.0
        + radial * (math.cos(angle) * first + math.sin(angle) * second)
        for angle in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)
    ]
    return [bond_A * vector for vector in others] + [bond_A * axis]


def _path_continuous(images: list[Atoms], maximum_A: float) -> bool:
    for left, right in zip(images, images[1:]):
        displacement, _ = find_mic(
            right.positions - left.positions,
            cell=left.cell,
            pbc=left.pbc,
        )
        if float(np.linalg.norm(displacement, axis=1).max()) > maximum_A:
            return False
    return True


def _path_jumps(images: list[Atoms]) -> list[dict[str, float | int]]:
    """Report the largest atom displacement for every adjacent image pair."""
    jumps: list[dict[str, float | int]] = []
    for image_index, (left, right) in enumerate(zip(images, images[1:])):
        displacement, _ = find_mic(
            right.positions - left.positions,
            cell=left.cell,
            pbc=left.pbc,
        )
        distances = np.linalg.norm(displacement, axis=1)
        atom_index = int(np.argmax(distances))
        jumps.append({
            'left_image': image_index,
            'right_image': image_index + 1,
            'atom_index': atom_index,
            'displacement_A': float(distances[atom_index]),
        })
    return jumps


def _relaxed_methane_bond(
    calculator: Any, config: NEBConfig
) -> tuple[float, float, dict[str, Any]]:
    methane = molecule('CH4')
    methane.set_cell((10.0, 10.0, 10.0))
    methane.center()
    methane.pbc = True
    methane.calc = calculator
    relaxation = relax_with_record(
        methane, 'neb_gas_ch4', config.endpoint_fmax_eV_A, config.endpoint_steps
    )
    if not relaxation['relax_neb_gas_ch4_converged']:
        raise RuntimeError('MLIP-NEB gas-phase methane did not converge')
    bonds = [methane.get_distance(0, index) for index in range(1, 5)]
    return (
        float(np.mean(bonds)),
        float(methane.get_potential_energy()),
        relaxation,
    )


def _path_result(
    images: list[Atoms], neb: NEB, plain: FIRE, climbing: FIRE | None,
    config: NEBConfig, converged: bool, reference_energy_eV: float,
) -> dict[str, Any]:
    energies = [float(image.get_potential_energy()) for image in images]
    relative = [value - reference_energy_eV for value in energies]
    forces = np.asarray(neb.get_forces(), dtype=float)
    max_force = float(np.linalg.norm(forces, axis=1).max())
    jumps = _path_jumps(images)
    continuous = all(
        float(jump['displacement_A']) <= config.max_image_displacement_A
        for jump in jumps
    )
    interior_maximum = int(np.argmax(energies)) not in (0, len(energies) - 1)
    valid = bool(
        converged and math.isfinite(max_force)
        and max_force <= config.path_fmax_eV_A
        and continuous and interior_maximum
    )
    apparent = max(relative)
    return {
        'energies_eV': energies,
        'reference_energy_eV': reference_energy_eV,
        'relative_energies_eV': relative,
        'initial_interaction_eV': energies[0] - reference_energy_eV,
        'barrier_eV': apparent if valid else None,
        'apparent_barrier_eV': apparent,
        'reaction_energy_eV': energies[-1] - reference_energy_eV,
        'path_max_force_eV_A': max_force,
        'plain_path_steps': int(plain.nsteps),
        'climbing_path_steps': int(climbing.nsteps) if climbing else 0,
        'plain_converged': climbing is not None,
        'climbing_converged': bool(converged),
        'path_continuous': continuous,
        'path_jumps': jumps,
        'interior_maximum': interior_maximum,
        'valid': valid,
    }


def _run_path(
    initial: Atoms, final: Atoms, calculator: Any, config: NEBConfig,
    reference_energy_eV: float,
) -> dict[str, Any]:
    images = [initial] + [initial.copy() for _ in range(config.images - 2)] + [final]
    neb = NEB(
        images, k=config.spring_eV_A2, climb=False,
        method='improvedtangent', allow_shared_calculator=True,
    )
    # Direct interpolation preserves the deliberately matched C/H reaction
    # coordinate. IDPP can move the carbon through a lower-energy surface site
    # before NEB starts, creating a nonphysical jump unrelated to C-H breaking.
    neb.interpolate(method='linear', apply_constraint=True)
    for image in images:
        image.calc = calculator
    plain = FIRE(neb, logfile=None, dt=0.03, dtmax=0.2, maxstep=0.05)
    if not plain.run(fmax=config.path_fmax_eV_A, steps=config.plain_path_steps):
        return _path_result(
            images, neb, plain, None, config, False, reference_energy_eV
        )
    plain_result = _path_result(
        images, neb, plain, None, config, True, reference_energy_eV
    )
    plain_result['plain_converged'] = True
    neb.climb = True
    climbing = FIRE(
        neb, logfile=None, dt=0.01, dtmax=0.05,
        maxstep=config.climbing_maxstep,
    )
    converged = bool(climbing.run(
        fmax=config.path_fmax_eV_A, steps=config.climbing_path_steps
    ))
    result = _path_result(
        images, neb, plain, climbing, config, converged, reference_energy_eV
    )
    result['plain_stage'] = {
        key: value for key, value in plain_result.items()
        if key not in {'energies_eV'}
    }
    return result


def methane_endpoints(
    structure: Atoms, active_index: int, calculator: Any, config: NEBConfig
) -> tuple[Atoms, Atoms, dict[str, Any]]:
    """Build atom-matched molten endpoints and relax both consistently.

    Args:
        structure: Relaxed clean catalyst surface.
        active_index: Surface atom on which methane is activated.
        calculator: Shared energy-conserving MLIP calculator.
        config: Frozen endpoint relaxation controls.

    Returns:
        Relaxed reactant and product states plus endpoint evidence.
    """
    site = structure[active_index].position.copy()
    final_h = site + np.array([1.5, 0.0, 1.2])
    initial_c = site + np.array([0.0, 0.0, config.reactant_height_A])
    methane_bond_A, gas_energy_eV, gas_relaxation = _relaxed_methane_bond(
        calculator, config
    )
    initial = structure.copy()
    initial.extend(Atoms(
        'CH4', positions=[initial_c] + [initial_c + vector for vector in
        _tetrahedral_vectors(final_h - initial_c, methane_bond_A)]
    ))
    final = structure.copy()
    carbon = site + np.array([0.0, 0.0, 2.0])
    final.append(Atom('C', carbon))
    final.extend(Atoms('H3', positions=(
        carbon + np.array([0.0, 1.02, 0.35]),
        carbon + np.array([-0.88, -0.51, 0.35]),
        carbon + np.array([0.88, -0.51, 0.35]),
    )))
    final.append(Atom('H', final_h))
    slab_constraint = _local_constraint(structure, active_index, config.local_radius_A)
    for endpoint in (initial, final):
        endpoint.set_cell(structure.cell)
        endpoint.set_pbc(structure.pbc)
        endpoint.calc = calculator
    initial.set_constraint(slab_constraint)
    carbon_index = len(structure)
    initial.set_constraint([
        slab_constraint,
        Hookean(
            carbon_index,
            active_index,
            k=5.0,
            rt=config.reactant_height_A,
        ),
    ])
    reactant = relax_with_record(
        initial, 'neb_reactant', config.endpoint_fmax_eV_A, config.endpoint_steps
    )
    initial.set_constraint(slab_constraint)
    if not reactant['relax_neb_reactant_converged']:
        raise RuntimeError('MLIP-NEB reactant endpoint did not converge')
    final.set_constraint([
        slab_constraint,
        Hookean(
            carbon_index,
            len(final) - 1,
            k=5.0,
            rt=config.product_CH_restraint_A,
        ),
    ])
    product = relax_with_record(
        final, 'neb_product', config.endpoint_fmax_eV_A, config.endpoint_steps
    )
    final.set_constraint(slab_constraint)
    if not product['relax_neb_product_converged']:
        raise RuntimeError('MLIP-NEB product endpoint did not converge')
    detached_ch_A = float(final.get_distance(carbon_index, len(final) - 1, mic=True))
    if detached_ch_A > config.maximum_product_CH_A:
        raise RuntimeError('MLIP-NEB product H diffused outside the reaction site')
    return initial, final, {
        'reactant_components': {
            'clean_surface_pre_relaxed': True,
            'gas_ch4': gas_relaxation,
            'gas_ch4_bond_A': methane_bond_A,
            'gas_ch4_energy_eV': gas_energy_eV,
        },
        'reactant': reactant,
        'product': product,
        'product_detached_CH_A': detached_ch_A,
    }


def run_replica(
    genome: tuple[Any, ...], structure_seed: int, device: str, config: NEBConfig
) -> dict[str, Any]:
    """Run one staged molten path and return convergence evidence.

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
        structure, 'neb_clean', config.endpoint_fmax_eV_A, config.endpoint_steps
    )
    if not clean['relax_neb_clean_converged']:
        raise RuntimeError('MLIP-NEB clean endpoint did not converge')
    active_index = _active_site(structure, active_indices, str(genome[2]))
    initial, final, endpoints = methane_endpoints(
        structure, active_index, calculator, config
    )
    reference_energy_eV = float(structure.get_potential_energy()) + float(
        endpoints['reactant_components']['gas_ch4_energy_eV']
    )
    return _run_path(
        initial, final, calculator, config, reference_energy_eV
    ) | {
        'genome': repr(genome), 'structure_seed': structure_seed,
        'model': MODEL_NAME, 'config': asdict(config),
        'active_index': active_index,
        'active_element': structure[active_index].symbol,
        'clean_geometry_sha256': geometry_digest(structure),
        'initial_geometry_sha256': geometry_digest(initial),
        'final_geometry_sha256': geometry_digest(final),
        'endpoint_relaxation': {'clean': clean} | endpoints,
    }


def run_ni111_reference(device: str, config: NEBConfig) -> dict[str, Any]:
    """Evaluate staged eSEN NEB on the independent Bengaard Ni(111) case.

    Args:
        device: Calculator device selected by the caller.
        config: Frozen endpoint and path controls.

    Returns:
        Complete convergence and barrier evidence for the solid reference.
    """
    from pipeline.validation.ni111_methane_reference import (
        REFERENCE_BARRIER_EV, build_ni111_methane_endpoints,
    )
    calculator = get_ocp_calculator(MODEL_NAME, device=device)
    if calculator is None:
        raise RuntimeError(f'could not load {MODEL_NAME}')
    initial, final = build_ni111_methane_endpoints()
    carbon = len(initial) - 5
    fixed = initial.constraints[0]
    clean = initial[:carbon]
    clean.set_constraint(fixed)
    clean.calc = calculator
    clean_relaxation = relax_with_record(
        clean, 'ni111_clean', config.endpoint_fmax_eV_A, config.endpoint_steps
    )
    if not clean_relaxation['relax_ni111_clean_converged']:
        raise RuntimeError('Ni(111) clean surface did not converge')
    initial.positions[:carbon] = clean.positions
    final.positions[:carbon] = clean.positions
    methane_bond_A, gas_energy_eV, gas_relaxation = _relaxed_methane_bond(
        calculator, config
    )
    direction = initial.positions[-1] - initial.positions[carbon]
    initial.positions[carbon:, 2] += (
        float(clean.positions[:, 2].max()) + config.reactant_height_A
        - float(initial.positions[carbon, 2])
    )
    initial.positions[carbon + 1:] = [
        initial.positions[carbon] + vector
        for vector in _tetrahedral_vectors(direction, methane_bond_A)
    ]
    initial.set_constraint(fixed)
    for endpoint in (initial, final):
        endpoint.calc = calculator
    final.set_constraint([
        final.constraints[0],
        Hookean(carbon, len(final) - 1, k=5.0, rt=config.product_CH_restraint_A),
    ])
    product = relax_with_record(
        final, 'ni111_product', config.endpoint_fmax_eV_A, config.endpoint_steps
    )
    final.set_constraint(fixed)
    if not product['relax_ni111_product_converged']:
        raise RuntimeError('Ni(111) product endpoint did not converge')
    detached_ch_A = float(final.get_distance(carbon, len(final) - 1, mic=True))
    if detached_ch_A > config.maximum_product_CH_A:
        raise RuntimeError('Ni(111) product H diffused outside the reaction site')
    reference_energy_eV = float(clean.get_potential_energy()) + gas_energy_eV
    result = _run_path(
        initial, final, calculator, config, reference_energy_eV
    )
    error = (
        abs(float(result['barrier_eV']) - REFERENCE_BARRIER_EV)
        if result['valid'] else None
    )
    return result | {
        'reference_id': 'bengaard-ni111-ch4-barrier',
        'reference_barrier_eV': REFERENCE_BARRIER_EV,
        'absolute_error_eV': error,
        'endpoint_relaxation': {
            'reactant_components': {
                'clean_surface': clean_relaxation,
                'gas_ch4': gas_relaxation,
                'gas_ch4_bond_A': methane_bond_A,
                'gas_ch4_energy_eV': gas_energy_eV,
            },
            'product': product,
            'product_detached_CH_A': detached_ch_A,
        },
    }


def lock_digest(path: str | Path) -> str:
    """Return the SHA-256 identity of a preregistration lock.

    Args:
        path: Preregistration lock file.

    Returns:
        Lowercase SHA-256 digest of the exact lock bytes.
    """
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _write_report(path: str | Path, report: dict[str, Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(target)


def main() -> None:
    """Run a locked reference or molten validation without reading holdouts."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lock', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--device', default='cuda')
    args = parser.parse_args()
    lock_path = Path(args.lock)
    lock = json.loads(lock_path.read_text())
    actual = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if actual != lock.get('implementation_sha256'):
        parser.error('NEB implementation does not match the preregistration lock')
    config = NEBConfig(**lock['neb_config'])
    if lock['kind'] == 'ni111_reference':
        record = run_ni111_reference(args.device, config)
        tolerance = float(lock['acceptance']['maximum_absolute_error_eV'])
        passed = bool(
            record['valid'] and record['absolute_error_eV'] is not None
            and record['absolute_error_eV'] <= tolerance
        )
        report = {
            'schema_version': 1, 'kind': lock['kind'],
            'lock_sha256': lock_digest(lock_path),
            'acceptance_passed': passed,
            'held_out_conditions_consumed': False, 'record': record,
        }
    elif lock['kind'] == 'molten_external':
        records = [
            run_replica(tuple(case['genome']), seed, args.device, config)
            for case in lock['external_cases'] for seed in case['structure_seeds']
        ]
        valid = [record for record in records if record['valid']]
        references = [
            float(value) for case in lock['external_cases']
            for value in case['reference_barriers_eV']
        ]
        predicted = [float(record['barrier_eV']) for record in valid]
        error = (
            abs(float(np.median(predicted)) - float(np.median(references)))
            if predicted else None
        )
        passed = bool(
            len(valid) >= int(lock['acceptance']['minimum_valid_paths'])
            and error is not None
            and error <= float(lock['acceptance']['maximum_median_absolute_error_eV'])
        )
        report = {
            'schema_version': 1, 'kind': lock['kind'],
            'lock_sha256': lock_digest(lock_path),
            'valid_paths': len(valid), 'attempted_paths': len(records),
            'median_predicted_barrier_eV': (
                float(np.median(predicted)) if predicted else None
            ),
            'median_reference_barrier_eV': float(np.median(references)),
            'median_absolute_error_eV': error, 'acceptance_passed': passed,
            'held_out_conditions_consumed': False, 'records': records,
        }
    else:
        parser.error(f"unsupported validation kind: {lock['kind']}")
    _write_report(args.output, report)
    print(json.dumps({
        key: value for key, value in report.items()
        if key not in {'record', 'records'}
    }, indent=2))


if __name__ == '__main__':
    main()
