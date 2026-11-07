"""Deterministic 3D ASE structure builders for all 14 catalyst material classes."""

from __future__ import annotations

import hashlib
import random
from typing import Tuple

import numpy as np
from ase import Atom, Atoms
from ase.build import bcc110, fcc100, fcc111, hcp0001
from ase.constraints import FixAtoms

# Shared lattice constants (Å) across FCC, BCC, HCP, and semimetal hosts.
LATTICE_CONSTANTS: dict[str, float] = {
    # FCC metals
    'Ni': 3.52,
    'Cu': 3.61,
    'Ag': 4.09,
    'Au': 4.08,
    'Al': 4.05,
    'Pd': 3.89,
    'Pt': 3.92,
    'Rh': 3.80,
    'Ir': 3.84,
    'Pb': 4.95,
    'Ca': 5.58,
    'Sr': 6.08,
    'Ce': 5.16,
    'Yb': 5.49,
    'Th': 5.08,
    # BCC metals
    'Fe': 2.87,
    'Cr': 2.88,
    'V': 3.03,
    'Nb': 3.30,
    'Mo': 3.15,
    'W': 3.16,
    'Ta': 3.30,
    'Na': 4.29,
    'K': 5.33,
    'Ba': 5.02,
    'Li': 3.51,
    'Cs': 6.14,
    'Rb': 5.59,
    # HCP metals (using equivalent FCC a)
    'Ti': 2.95,
    'Zr': 3.23,
    'Hf': 3.19,
    'Co': 2.51,
    'Ru': 2.71,
    'Os': 2.74,
    'Re': 2.76,
    'Sc': 3.31,
    'Y': 3.65,
    'La': 3.75,
    'Mg': 3.21,
    'Zn': 2.66,
    'Cd': 2.98,
    # Non-standard / semimetals (use nearest FCC equivalent)
    'Sn': 5.83,
    'In': 4.60,
    'Ga': 4.52,
    'Sb': 4.31,
    'Bi': 4.75,
    'Ge': 5.66,
    'Si': 5.43,
    'Te': 4.45,
    'Se': 4.36,
    'As': 4.13,
    'Mn': 3.50,
    'Pr': 5.16,
    'Nd': 5.08,
    'Gd': 3.63,
    'Tb': 3.60,
    'Dy': 3.59,
    'Ho': 3.58,
    'Er': 3.56,
    'Tm': 3.54,
    'Lu': 3.50,
    'Sm': 3.62,
    'Eu': 4.58,
    'Tl': 3.46,
    'Po': 3.35,
}

_SLAB_FACETS = ('fcc111', 'fcc100', 'bcc110', 'hcp0001')
_MILLER_FACETS = {'111': 'fcc111', '100': 'fcc100', '110': 'bcc110', '211': 'fcc111'}
_MAX_FACETS = {'basal_0001': 'hcp0001', 'edge_1010': 'bcc110', 'edge_1120': 'fcc100'}
_SPINEL_MORPHOLOGY_DEFECTS = {
    'mesoporous': 'O_vacancy',
    'hollow_sphere': 'A_vacancy',
    'nanorod': 'B_vacancy',
}
_MXENE_TERMINATIONS = {
    'OH': ['O', 'H'],
    'O': ['O'],
    'F': ['F'],
    'Cl': ['Cl'],
    'S': ['S'],
    'mixed_OH_O': ['O', 'H', 'O'],
    'mixed_OH_F': ['O', 'H', 'F'],
    'bare': [],
}
_CARBON_N_CAVITIES = {
    'pyridinic': 'N3C',
    'pyrrolic': 'N4_pyrrole',
    'graphitic': 'N4',
    'oxidized': 'N3O',
    'mixed_pyridinic_graphitic': 'N4C2',
}


def _stable_seed(*parts: object) -> int:
    """Return a process-independent seed for one physical realization."""
    payload = repr(parts).encode('utf-8')
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], 'big')


def _label_fraction(label: object) -> float:
    """Stable [0, 1) coordinate used to distinguish categorical realizations."""
    return (_stable_seed('label-coordinate-v1', label) % 10_000) / 10_000.0


def _representative_element(label: str, default: str = 'C') -> str:
    """Map support/linker labels to a chemically relevant explicit atom."""
    text = str(label)
    for token in (
        'Mo',
        'Ti',
        'Al',
        'Si',
        'Mg',
        'Zr',
        'Ce',
        'La',
        'Fe',
        'Cl',
        'Se',
        'Zn',
        'Ni',
        'Co',
        'Cu',
        'Ca',
        'Na',
        'B',
        'N',
        'O',
        'S',
        'P',
        'F',
        'C',
    ):
        if token in text:
            return token
    return default


def _append_marker(
    atoms: Atoms, symbol: str, offset: tuple[float, float, float]
) -> int:
    """Add one explicit environment atom without overlapping the active site."""
    center = atoms.get_center_of_mass()
    atoms.append(Atom(symbol, position=center + np.asarray(offset, dtype=float)))
    return len(atoms) - 1


def _add_surface_markers(atoms: Atoms, symbols: list[str]) -> list[int]:
    """Add sparse subsurface/backside environment atoms and restore bottom constraint."""
    if not symbols:
        return []
    atoms.set_constraint()
    z_min = float(atoms.positions[:, 2].min())
    cell = atoms.cell.lengths()
    added = []
    for index, symbol in enumerate(symbols):
        x = (index + 1) * float(cell[0]) / (len(symbols) + 1)
        y = (index % 2 + 1) * float(cell[1]) / 3.0
        atoms.append(
            Atom(symbol, position=(x, y, max(0.6, z_min - 1.6 - 0.55 * index)))
        )
        added.append(len(atoms) - 1)
    z = atoms.positions[:, 2]
    atoms.set_constraint(FixAtoms(mask=z < z_min + 4.0))
    return added


def generate_alloy_slab(
    host: str,
    facet: str,
    strain: float,
    dopants: tuple,
    n_sub: int,
    n_vac: int,
    size: Tuple[int, int, int] = (3, 3, 4),
    vacuum: float = 12.0,
    seed_key=None,
    environment_label: str | None = None,
) -> Tuple[Atoms, list]:
    """Generate a strained alloy slab with dopant substitutions and vacancies.

    Args:
        host: Host metal element symbol.
        facet: Surface facet identifier ('fcc111', 'fcc100', 'bcc110', 'hcp0001').
        strain: Biaxial lattice strain fraction.
        dopants: Ordered dopant element symbols.
        n_sub: Number of top-layer dopant substitutions.
        n_vac: Number of top-layer vacancies.
        size: Supercell dimensions (nx, ny, nz).
        vacuum: Vacuum padding along z in angstrom.
        seed_key: Deterministic seed key for dopant/vacancy placement.
        environment_label: Optional support/environment label.

    Returns:
        Tuple of `(slab, top_indices)`.
    """
    builders = {
        'fcc111': fcc111,
        'fcc100': fcc100,
        'bcc110': bcc110,
        'hcp0001': hcp0001,
    }
    builder = builders.get(facet, fcc111)
    a = LATTICE_CONSTANTS.get(host, 3.60)
    try:
        slab = builder(host, size=size, vacuum=vacuum, a=a)
    except Exception:
        slab = fcc111('Ni', size=size, vacuum=vacuum, a=a)
        for atom in slab:
            atom.symbol = host

    z_coords = slab.positions[:, 2]
    z_max = z_coords.max()
    top_mask = z_coords > (z_max - 2.5)
    top_indices = np.where(top_mask)[0]

    rng = random.Random(
        _stable_seed(
            'alloy-slab-v2',
            host,
            facet,
            tuple(dopants),
            int(n_sub),
            int(n_vac),
            tuple(size),
            float(vacuum),
            seed_key,
        )
    )
    shuffled = list(top_indices)
    rng.shuffle(shuffled)

    sub_count = (
        min(max(n_sub, len(dopants)), len(shuffled))
        if (dopants and n_sub > 0)
        else 0
    )
    for i in range(sub_count):
        slab[shuffled[i]].symbol = dopants[i % len(dopants)]

    vac_count = min(n_vac, len(shuffled) - sub_count)
    if vac_count > 0:
        del_indices = sorted(shuffled[sub_count : sub_count + vac_count], reverse=True)
        for idx in del_indices:
            del slab[idx]
        z_coords = slab.positions[:, 2]
        top_mask = z_coords > (z_coords.max() - 2.5)
        top_indices = np.where(top_mask)[0]
    elif sub_count > 0:
        sub_set = set(int(idx) for idx in shuffled[:sub_count])
        top_indices = [int(idx) for idx in shuffled[:sub_count]] + [
            int(idx) for idx in top_indices if int(idx) not in sub_set
        ]

    if environment_label and str(environment_label).lower() != 'none':
        bottom = int(np.argmin(slab.positions[:, 2]))
        slab[bottom].symbol = _representative_element(environment_label, 'C')
        slab.positions[bottom, 0] += 0.05 * _label_fraction(environment_label)

    cell = slab.get_cell()
    cell[0] *= 1.0 + strain
    cell[1] *= 1.0 + strain
    slab.set_cell(cell, scale_atoms=True)

    z_coords = slab.positions[:, 2]
    z_min = z_coords.min()
    fix_mask = z_coords < (z_min + 4.0)
    slab.set_constraint(FixAtoms(mask=fix_mask))

    return slab, [int(idx) for idx in top_indices]


def generate_porphyrin_cluster(
    metal: str,
    cavity: str,
    d_metal_n: float = 2.0,
    substrate: str = 'N-graphene',
    axial_ligand: str = 'none',
    linker: str | None = None,
    pore_size: float | None = None,
) -> Atoms:
    """Generate a metal-porphyrin active site cluster for SAC/MOF evaluation.

    Args:
        metal: Central metal element symbol (or 'None' for metal-free).
        cavity: Coordination motif identifier (e.g. 'N4', 'N3C').
        d_metal_n: Metal-ligand bond distance in angstrom.
        substrate: Support/substrate identifier.
        axial_ligand: Axial ligand identifier.
        linker: Optional framework linker label.
        pore_size: Optional framework pore size in angstrom.

    Returns:
        Periodic ASE `Atoms` cluster representation.
    """
    atoms = Atoms()
    if metal != 'None':
        atoms.append(Atom(metal, position=(0.0, 0.0, 0.0)))
    else:
        atoms.append(Atom('C', position=(0.0, 0.0, 0.0)))

    coord_map = {
        'N4': ['N', 'N', 'N', 'N'],
        'N3C': ['N', 'N', 'N', 'C'],
        'N2C2': ['N', 'N', 'C', 'C'],
        'N3B': ['N', 'N', 'N', 'B'],
        'N3S': ['N', 'N', 'N', 'S'],
        'N3P': ['N', 'N', 'N', 'P'],
        'N2S2': ['N', 'N', 'S', 'S'],
        'N2O2': ['N', 'N', 'O', 'O'],
        'N3O': ['N', 'N', 'N', 'O'],
        'O4': ['O', 'O', 'O', 'O'],
        'N4_pyridine': ['N', 'N', 'N', 'N'],
        'N4_pyrrole': ['N', 'N', 'N', 'N'],
        'N6': ['N', 'N', 'N', 'N', 'N', 'N'],
        'N8': ['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'],
        'N4C2': ['N', 'N', 'N', 'N', 'C', 'C'],
        'N3SN3': ['N', 'N', 'N', 'S', 'N', 'N'],
        'N4N4': ['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'],
        'N2P2': ['N', 'N', 'P', 'P'],
    }
    coord_atoms = coord_map.get(cavity, ['N', 'N', 'N', 'N'])

    if len(coord_atoms) <= 4:
        positions = [
            (-d_metal_n, 0.0, 0.0),
            (d_metal_n, 0.0, 0.0),
            (0.0, -d_metal_n, 0.0),
            (0.0, d_metal_n, 0.0),
        ]
    else:
        n_coord = len(coord_atoms)
        positions = [
            (
                d_metal_n * np.cos(2 * np.pi * k / n_coord),
                d_metal_n * np.sin(2 * np.pi * k / n_coord),
                0.0,
            )
            for k in range(n_coord)
        ]

    for sym, pos in zip(coord_atoms, positions[: len(coord_atoms)]):
        atoms.append(Atom(sym, position=pos))
    if len(atoms) > 1:
        atoms.positions[1, 2] += 0.08 * _label_fraction(cavity)

    n_c = 8
    pore_scale = (
        1.0 if pore_size is None else np.clip(float(pore_size) / 12.0, 0.75, 1.75)
    )
    for k in range(n_c):
        angle = 2 * np.pi * k / n_c + np.pi / n_c
        r = (d_metal_n + 1.2) * pore_scale
        atoms.append(Atom('C', position=(r * np.cos(angle), r * np.sin(angle), 0.0)))

    for k in range(n_c):
        angle = 2 * np.pi * k / n_c + np.pi / n_c
        r = (d_metal_n + 2.2) * pore_scale
        atoms.append(Atom('H', position=(r * np.cos(angle), r * np.sin(angle), 0.0)))

    substrate_element = _representative_element(substrate, 'C')
    atoms[5].symbol = substrate_element
    atoms.positions[5, 2] += 0.08 * _label_fraction(substrate)
    if linker is not None:
        atoms[6].symbol = _representative_element(linker, 'C')
        atoms.positions[6, 2] += 0.08 * _label_fraction(linker)

    ligand = str(axial_ligand).lower()
    if ligand not in ('none', '') and metal != 'None':
        ligand_symbol = _representative_element(axial_ligand, 'O')
        _append_marker(
            atoms,
            ligand_symbol,
            (0.0, 0.0, -(1.9 + 0.3 * _label_fraction(axial_ligand))),
        )

    atoms.set_cell((15.0, 15.0, 15.0))
    atoms.center()
    atoms.pbc = True
    return atoms


def adsorption_site(
    atoms: Atoms, active_indices: list[int], distance: float
) -> tuple[np.ndarray, np.ndarray]:
    """Choose the less crowded side of a surface or molecular active site.

    Args:
        atoms: Relaxed catalyst structure.
        active_indices: Atom indices defining the active site.
        distance: Adsorbate distance from the active-site centroid in angstrom.

    Returns:
        Adsorbate position and outward unit direction.
    """
    indices = [index for index in active_indices if 0 <= index < len(atoms)]
    base = (
        atoms.positions[indices].mean(axis=0)
        if indices
        else atoms.positions.mean(axis=0)
    )
    directions = (np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0]))

    def clearance(direction: np.ndarray) -> float:
        position = base + distance * direction
        return float(np.linalg.norm(atoms.positions - position, axis=1).min())

    direction = max(directions, key=clearance)
    return base + distance * direction, direction


def screening_adsorption_site(
    atoms: Atoms,
    active_indices: list[int],
    material_class: str,
    distance: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Preserve legacy slab atop sites while avoiding SAC/DAC collisions.

    Args:
        atoms: Relaxed catalyst structure.
        active_indices: Atom indices defining the active site.
        material_class: Encoded catalyst material class.
        distance: Adsorbate distance from the active site in angstrom.

    Returns:
        Adsorbate position and outward unit direction.
    """
    if material_class in ('SAC', 'DAC'):
        return adsorption_site(atoms, active_indices, distance)
    if active_indices and 0 <= active_indices[0] < len(atoms):
        base = atoms[active_indices[0]].position.copy()
    else:
        base = atoms.positions.mean(axis=0)
    direction = np.array([0.0, 0.0, 1.0])
    return base + distance * direction, direction


def _slab_facet_key(facet: str) -> str:
    """Map an extended facet label to an available slab builder facet."""
    facet_key = facet.split('_')[0] if '_' in facet else facet
    return facet_key if facet_key in _SLAB_FACETS else 'fcc111'


def _top_layer_indices(
    slab: Atoms, preferred_symbols: tuple[str, ...] = ()
) -> list[int]:
    """Return indices of atoms within 3 angstrom of the slab top, prioritizing exposed active cations."""
    z = slab.positions[:, 2]
    z_max = float(z.max())
    raw = [int(i) for i in np.where(z > z_max - 3.0)[0]]
    pref = {s for s in preferred_symbols if s and s != 'None'}

    def _sort_key(idx: int) -> tuple[int, float, int]:
        sym = slab[idx].symbol
        in_outer_plane = float(z[idx]) >= z_max - 0.85
        is_pref = sym in pref
        tier = (
            0
            if (in_outer_plane and is_pref)
            else (
                1
                if (in_outer_plane and sym != 'O')
                else (2 if is_pref else (3 if sym != 'O' else 4))
            )
        )
        return (tier, -float(z[idx]), idx)

    return sorted(raw, key=_sort_key)


# Effective Shannon 6/12-coordinate ionic radii (Å) for ABO3 / AB2O4 lattice scaling
_CATION_RADII_A: dict[str, float] = {
    'La': 1.36,
    'Sr': 1.44,
    'Ba': 1.61,
    'Ca': 1.34,
    'Ce': 1.14,
    'Pr': 1.18,
    'Nd': 1.16,
    'Sm': 1.13,
    'Gd': 1.10,
    'Y': 1.08,
    'Mg': 0.72,
    'Zn': 0.74,
    'Ti': 0.605,
    'V': 0.64,
    'Cr': 0.615,
    'Mn': 0.645,
    'Fe': 0.645,
    'Co': 0.61,
    'Ni': 0.69,
    'Cu': 0.73,
    'Zr': 0.72,
    'Nb': 0.64,
    'Mo': 0.65,
    'Ru': 0.62,
    'Rh': 0.665,
    'Al': 0.535,
    'Ga': 0.62,
    'In': 0.80,
    'Sn': 0.69,
}


def _generate_perovskite_slab(
    A: str, B: str, dopant: str, frac: float, defect: str, seed_key=None
) -> Atoms:
    """Generate a BO2-terminated Pm-3m ABO3(001) perovskite slab with composition-scaled lattice constant."""
    r_a = _CATION_RADII_A.get(A, 1.34)
    r_b = _CATION_RADII_A.get(B, 0.64)
    r_dop = _CATION_RADII_A.get(dopant, r_b) if dopant != 'None' else r_b
    r_b_eff = (1.0 - float(frac)) * r_b + float(frac) * r_dop
    # Cubic Pm-3m lattice constant from Goldschmidt A-O and B-O bond lengths (r_O = 1.40 Å)
    a_ao = (r_a + 1.40) * np.sqrt(2.0)
    a_bo = 2.0 * (r_b_eff + 1.40)
    a = float(np.clip(0.5 * (a_ao + a_bo), 3.76, 4.12))
    atoms = Atoms()
    a_sites, b_sites, o_sites = [], [], []
    for ix in range(2):
        for iy in range(2):
            for iz in range(3):
                base = np.array([ix * a, iy * a, iz * a])
                atoms.append(Atom(A, position=base))
                a_sites.append(len(atoms) - 1)
                atoms.append(Atom(B, position=base + np.array([a / 2, a / 2, a / 2])))
                b_sites.append(len(atoms) - 1)
                atoms.append(Atom('O', position=base + np.array([a / 2, a / 2, 0])))
                o_sites.append(len(atoms) - 1)
                atoms.append(Atom('O', position=base + np.array([a / 2, 0, a / 2])))
                o_sites.append(len(atoms) - 1)
                atoms.append(Atom('O', position=base + np.array([0, a / 2, a / 2])))
                o_sites.append(len(atoms) - 1)

    rng = random.Random(
        _stable_seed('perovskite-v2', A, B, dopant, float(frac), defect, seed_key)
    )
    if dopant != 'None' and frac > 0:
        count = max(1, min(len(b_sites), round(float(frac) * len(b_sites))))
        # Ensure at least one top-layer B-site is substituted when doped so the exposed
        # BO2(001) termination probes the dopant-modified active site.
        top_b = [idx for idx in b_sites if atoms.positions[idx, 2] > 2.0 * a]
        chosen = []
        if top_b:
            chosen.append(rng.choice(top_b))
        rem_pool = [idx for idx in b_sites if idx not in chosen]
        if count > len(chosen) and rem_pool:
            chosen.extend(rng.sample(rem_pool, min(count - len(chosen), len(rem_pool))))
        for index in chosen:
            atoms[index].symbol = dopant

    vacancy_pool = {
        'A_vacancy': a_sites,
        'B_vacancy': b_sites,
        'O_vacancy': o_sites,
    }.get(defect)
    if vacancy_pool:
        del atoms[rng.choice(vacancy_pool)]
    elif defect == 'A_excess':
        atoms.append(
            Atom(A, position=(a / 2, a / 2, float(atoms.positions[:, 2].max()) + 1.8))
        )

    cell = [2 * a, 2 * a, 3 * a + 12.0]
    atoms.set_cell(cell)
    atoms.pbc = True
    z = atoms.positions[:, 2]
    atoms.set_constraint(FixAtoms(mask=z < z.min() + 2.0))
    return atoms


def _generate_spinel_slab(
    A: str, B: str, dopant: str, defect: str, support: str, seed_key=None
) -> Atoms:
    """Build an octahedral-B-terminated AB2O4 spinel slab with composition-scaled lattice constant."""
    atoms = Atoms()
    r_a = _CATION_RADII_A.get(A, 0.72)
    r_b = _CATION_RADII_A.get(B, 0.62)
    a = float(np.clip(4.02 + 0.22 * (r_a - 0.70) + 0.30 * (r_b - 0.60), 3.95, 4.30))
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                base = np.array([ix * a, iy * a, iz * a])
                atoms.append(Atom(A, position=base))
                atoms.append(Atom(B, position=base + (a / 4, a / 4, a / 4)))
                atoms.append(Atom(B, position=base + (3 * a / 4, 3 * a / 4, 3 * a / 4)))
                for offset in (
                    (a / 2, 0, 0),
                    (0, a / 2, 0),
                    (0, 0, a / 2),
                    (a / 2, a / 2, a / 2),
                ):
                    atoms.append(Atom('O', position=base + offset))
    rng = random.Random(
        _stable_seed('spinel-v2', A, B, dopant, defect, support, seed_key)
    )
    if dopant != 'None':
        b_sites = [i for i, atom in enumerate(atoms) if atom.symbol == B]
        top_b = [i for i in b_sites if atoms.positions[i, 2] > 1.5 * a]
        atoms[rng.choice(top_b if top_b else b_sites)].symbol = dopant
    vacancy_symbol = {'A_vacancy': A, 'B_vacancy': B, 'O_vacancy': 'O'}.get(defect)
    if vacancy_symbol:
        choices = [i for i, atom in enumerate(atoms) if atom.symbol == vacancy_symbol]
        if choices:
            del atoms[rng.choice(choices)]
    if support != 'none':
        atoms[0].symbol = _representative_element(support, 'C')
    atoms.set_cell((2 * a, 2 * a, 2 * a + 12.0))
    atoms.pbc = True
    z = atoms.positions[:, 2]
    atoms.set_constraint(FixAtoms(mask=z < z.min() + 2.0))
    return atoms


def _generate_hydride_slab(
    metal: str,
    second: str,
    h_type: str = 'simple',
    additive: str = 'None',
    temperature_K: float = 300,
    seed_key=None,
) -> Atoms:
    """Generate a metal-hydride slab: metal FCC + interstitial H."""
    a = LATTICE_CONSTANTS.get(metal, 3.60)
    try:
        slab = fcc111(metal, size=(3, 3, 3), vacuum=12.0, a=a)
    except Exception:
        slab = fcc111('Ni', size=(3, 3, 3), vacuum=12.0, a=a)
        for atom in slab:
            atom.symbol = metal

    z = slab.positions[:, 2]
    top_z = z.max()
    h_positions = []
    for atom in slab:
        if atom.position[2] > top_z - 2.5:
            h_pos = atom.position.copy()
            h_pos[2] += 1.0
            h_positions.append(h_pos)
    for h_pos in h_positions:
        slab.append(Atom('H', position=h_pos))

    extra_h = {
        'simple': 0,
        'intermetallic_AB': 1,
        'intermetallic_AB2': 2,
        'intermetallic_AB5': 3,
        'intermetallic_A2B': 2,
        'complex_alanate': 4,
        'complex_borohydride': 4,
        'complex_amide': 3,
        'perovskite_hydride': 2,
    }.get(h_type, 0)
    for index in range(extra_h):
        slab.append(
            Atom(
                'H',
                position=(
                    (index + 1) * slab.cell.lengths()[0] / (extra_h + 1),
                    (index % 2 + 1) * slab.cell.lengths()[1] / 3.0,
                    top_z + 2.5 + 0.35 * index,
                ),
            )
        )

    if second != 'None' and second != metal:
        indices = [
            i
            for i in range(len(slab))
            if slab[i].symbol == metal and slab[i].position[2] > top_z - 2.5
        ]
        for i in indices[:2]:
            slab[i].symbol = second

    if additive != 'None':
        added = _add_surface_markers(slab, [_representative_element(additive, 'C')])
        slab.positions[added[0], 0] += 0.1 * _label_fraction(additive)

    thermal_scale = 1.0 + (float(temperature_K) - 300.0) * 1.0e-5
    slab.set_cell(slab.cell * [thermal_scale, thermal_scale, 1.0], scale_atoms=True)

    z_all = slab.positions[:, 2]
    slab.set_constraint(FixAtoms(mask=z_all < z_all.min() + 3.0))
    return slab


def _build_supported_nanoparticle(genome, structure_seed):
    _, metal, support, facet, diameter_nm, loading_wt_pct, dopants, n_sub, n_vac = genome
    d_nm = max(float(diameter_nm), 0.5)
    load_pct = min(max(float(loading_wt_pct), 0.0), 100.0)
    # Spherical nanoparticle dispersion D ~ 1.1 / d_nm (Bergeret & Gallezot) and
    # cuboctahedral corner/edge step-site fraction f_step ~ (2.5 / d_nm)^2 (Van Hardeveld & Hartog).
    dispersion = min(1.0, 1.1 / d_nm)
    step_fraction = min(1.0, (2.5 / d_nm) ** 2)
    # Young-Laplace capillary compression (-1.5% at 1 nm) plus support-epitaxy tensile offset.
    capillary_strain = -0.015 / max(d_nm, 1.0)
    loading_strain = 0.002 * (1.0 - min(load_pct, 30.0) / 30.0)
    effective_strain = float(np.clip(capillary_strain + loading_strain, -0.03, 0.03))
    effective_n_vac = max(int(n_vac), 1 if d_nm < 5.0 else 0)
    slab, top_idx = generate_alloy_slab(
        metal,
        _slab_facet_key(facet),
        effective_strain,
        tuple(dopants),
        int(n_sub),
        effective_n_vac,
        seed_key=(genome, int(structure_seed)),
        environment_label=support,
    )
    slab.info.update(
        {
            'particle_diameter_nm': float(diameter_nm),
            'loading_wt_pct': float(loading_wt_pct),
            'support': str(support),
            'dispersion': float(dispersion),
            'step_site_fraction': float(step_fraction),
            'effective_strain': effective_strain,
            'effective_surface_vacancies': int(effective_n_vac),
        }
    )
    return slab, top_idx


def _build_solid_catalyst(genome, structure_seed):
    _, metal, support, facet, strain, dopants, n_sub, n_vac = genome
    return generate_alloy_slab(
        metal,
        _slab_facet_key(facet),
        strain,
        dopants,
        n_sub,
        n_vac,
        seed_key=genome,
        environment_label=support,
    )


def _build_molten_metal(genome, structure_seed):
    _, host, promoter, at_pct, temp = genome
    # Map atomic percent onto the 9 sites of the 3x3 top surface layer; apply linear
    # thermal expansion alpha = 1.0e-5 K^-1 relative to the 1000 K reference state.
    dopants = (promoter,) if (promoter != 'None' and float(at_pct) > 0.0) else ()
    n_sub = max(1, round(float(at_pct) / 100.0 * 9)) if dopants else 0
    thermal_strain = (float(temp) - 1000.0) * 1.0e-5
    facet = 'hcp0001' if host in {'Mg', 'Zn', 'Cd'} else 'fcc111'
    slab, top_idx = generate_alloy_slab(
        host,
        facet,
        thermal_strain,
        dopants,
        n_sub,
        0,
        seed_key=('MoltenMetal', host, promoter, float(at_pct), structure_seed),
    )
    # Apply deterministic liquid-metal positional disorder (capillary/thermal
    # corrugation amplitude 0.08 A * sqrt(T / 1000 K)) and Gibbs surface-tension
    # relaxation where higher-surface-tension transition-metal promoters relax 0.14 A
    # inward (-z) relative to low-surface-tension liquid p-block hosts (Upham et al. 2017).
    rng_liq = random.Random(
        _stable_seed(
            'molten-liquid-corrugation-v1',
            host,
            promoter,
            float(at_pct),
            float(temp),
            structure_seed,
        )
    )
    amp = 0.08 * float(np.sqrt(max(float(temp), 300.0) / 1000.0))
    high_gamma_promoters = {
        'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
        'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ta', 'W', 'Re', 'Pt',
    }
    for idx in top_idx:
        dx = rng_liq.uniform(-amp, amp)
        dy = rng_liq.uniform(-amp, amp)
        dz = rng_liq.uniform(-0.5 * amp, 0.5 * amp)
        if slab[int(idx)].symbol in high_gamma_promoters and host not in high_gamma_promoters:
            dz -= 0.14
        slab.positions[int(idx)] += np.array([dx, dy, dz], dtype=float)
    return slab, top_idx


def _build_sac(genome, structure_seed):
    cluster = generate_porphyrin_cluster(
        genome[1], genome[2], substrate=genome[3], axial_ligand=genome[4]
    )
    return cluster, [0]


def _build_dac(genome, structure_seed):
    _, m1, m2, coord, substrate = genome
    cluster = generate_porphyrin_cluster(m1, coord, substrate=substrate)
    separation = 2.5
    cluster[0].position[0] -= separation / 2.0
    second = cluster[0].position.copy()
    second[0] += separation
    midpoint_x = 0.5 * (cluster[0].position[0] + second[0])
    coordination_count = 8 if coord in ('N8', 'N4N4') else 6
    for index in range(1, 1 + coordination_count):
        atom = cluster[index]
        if (
            min(
                np.linalg.norm(atom.position - cluster[0].position),
                np.linalg.norm(atom.position - second),
            )
            < 1.2
        ):
            atom.position[0] = (
                second[0] + 2.0
                if atom.position[0] >= midpoint_x
                else cluster[0].position[0] - 2.0
            )
    cluster.append(Atom(m2, position=second))
    return cluster, [0, len(cluster) - 1]


def _build_framework(genome, structure_seed):
    _, metal, linker, cavity, pore = genome
    cluster = generate_porphyrin_cluster(
        metal, cavity, substrate='framework', linker=linker, pore_size=pore
    )
    return cluster, [0]


def _build_perovskite(genome, structure_seed):
    _, A, B, dopant, frac, defect = genome
    slab = _generate_perovskite_slab(A, B, dopant, frac, defect, seed_key=genome)
    return slab, _top_layer_indices(slab, preferred_symbols=(str(dopant), str(B)))


def _build_metal_hydride(genome, structure_seed):
    _, metal, h_type, second, additive, temp = genome
    slab = _generate_hydride_slab(
        metal,
        second,
        h_type=h_type,
        additive=additive,
        temperature_K=temp,
        seed_key=genome,
    )
    return slab, _top_layer_indices(slab, preferred_symbols=(str(second), str(metal)))


def _insert_octahedral_interstitial_sublayers(
    slab: Atoms, x_symbols: list[str]
) -> list[int]:
    """Insert octahedral interstitial C/N atoms at interlayer midplanes of an M_{n+1}X_n slab."""
    if not x_symbols:
        return []
    slab.set_constraint()
    z_vals = np.sort(np.unique(np.round(slab.positions[:, 2], 2)))
    cell_len = slab.cell.lengths()
    lx, ly = float(cell_len[0]), float(cell_len[1])
    mid_zs = (
        [0.5 * float(z_vals[k] + z_vals[k + 1]) for k in range(len(z_vals) - 1)]
        if len(z_vals) >= 2
        else [float(slab.positions[:, 2].mean())]
    )
    # Candidate hollow-site fractional (fx, fy) grid in the lateral unit cell
    hollow_grid = [
        (1.0 / 3.0, 2.0 / 3.0),
        (2.0 / 3.0, 1.0 / 3.0),
        (0.50, 0.50),
        (0.25, 0.75),
        (0.75, 0.25),
        (0.18, 0.42),
        (0.82, 0.58),
    ]
    added: list[int] = []
    for idx_s, sym in enumerate(x_symbols):
        z_mid = mid_zs[idx_s % len(mid_zs)]
        best_pos = None
        best_min_d = -1.0
        for shift in range(len(hollow_grid)):
            fx, fy = hollow_grid[(idx_s + shift) % len(hollow_grid)]
            cand = np.array([fx * lx, fy * ly, z_mid], dtype=float)
            min_d = float(np.linalg.norm(slab.positions - cand, axis=1).min())
            if min_d > best_min_d:
                best_min_d = min_d
                best_pos = cand
        if best_pos is not None and best_min_d >= 1.15:
            slab.append(Atom(sym, position=best_pos))
            added.append(len(slab) - 1)
        else:
            added.extend(_add_surface_markers(slab, [sym]))
    z_all = slab.positions[:, 2]
    z_min = float(z_all.min())
    slab.set_constraint(FixAtoms(mask=z_all < z_min + 3.5))
    return added


def _build_max_phase(genome, structure_seed):
    _, M, A, X, n_val, dopant, facet = genome
    n_int = max(1, int(n_val))
    # On basal_0001, the top catalytic layer is M (substituted by dopant if present) while
    # the A-group nanolaminate layer resides in the subsurface plane; on prism edge facets
    # (edge_1010, edge_1120), both M/dopant and A are exposed at the surface step.
    if str(facet) == 'basal_0001':
        top_dopants = (dopant,) if str(dopant) != 'None' else ()
    else:
        top_dopants = tuple(x for x in (dopant, A) if str(x) != 'None')
    slab, top_idx = generate_alloy_slab(
        M,
        _MAX_FACETS.get(facet, 'hcp0001'),
        0.0,
        top_dopants,
        min(len(top_dopants), 2),
        0,
        size=(3, 3, n_int + 1),
        seed_key=genome,
    )
    if str(facet) == 'basal_0001' and str(A) != 'None':
        # Substitute A-group atoms into the subsurface nanolaminate plane
        z_coords = slab.positions[:, 2]
        z_max = float(z_coords.max())
        sub_indices = [
            int(i)
            for i in range(len(slab))
            if int(i) not in set(top_idx) and float(z_coords[i]) < z_max - 1.5
        ]
        for sub_i in sub_indices[:3]:
            slab[sub_i].symbol = str(A)
    _insert_octahedral_interstitial_sublayers(slab, [str(X)] * n_int)
    return slab, top_idx


def _build_hea(genome, structure_seed):
    _, components, structure, facet_str, temp = genome
    host = components[0]
    dopants = components[1:]
    facet_key = {'bcc': 'bcc110', 'hcp': 'hcp0001'}.get(
        structure, _MILLER_FACETS.get(facet_str, 'fcc111')
    )
    thermal_strain = (float(temp) - 1000.0) * 5.0e-6
    if structure == 'amorphous':
        thermal_strain += 0.003
    elif structure == 'fcc_bcc_dual':
        thermal_strain -= 0.003
    return generate_alloy_slab(
        host,
        facet_key,
        thermal_strain,
        dopants,
        min(len(dopants), 6),
        0,
        seed_key=genome,
    )


def _build_spinel(genome, structure_seed):
    _, A, B, dopant, morph, support = genome
    slab = _generate_spinel_slab(
        A,
        B,
        dopant,
        _SPINEL_MORPHOLOGY_DEFECTS.get(morph, 'none'),
        support,
        seed_key=genome,
    )
    morphology_scale = 0.98 + 0.04 * _label_fraction(morph)
    slab.set_cell(slab.cell * [morphology_scale, 1.0, 1.0], scale_atoms=True)
    return slab, _top_layer_indices(slab, preferred_symbols=(str(dopant), str(B)))


def _build_mxene(genome, structure_seed):
    _, M, X_elem, n_val, term, sac_metal = genome
    n_int = max(1, int(n_val))
    slab, top_idx = generate_alloy_slab(
        M,
        'hcp0001',
        0.0,
        (sac_metal,) if sac_metal != 'None' else (),
        1 if sac_metal != 'None' else 0,
        0,
        size=(3, 3, n_int + 1),
        seed_key=genome,
    )
    x_symbols = list(('C', 'N') if X_elem == 'CN' else (X_elem,))
    _insert_octahedral_interstitial_sublayers(slab, x_symbols * n_int)
    term_syms = _MXENE_TERMINATIONS.get(term, [])
    if term_syms:
        _add_surface_markers(slab, list(term_syms))
    return slab, top_idx


def _build_saa(genome, structure_seed):
    _, trace, host, facet_str, loading = genome
    lateral = 4 if loading <= 100 else (3 if loading <= 1000 else 2)
    slab, top_idx = generate_alloy_slab(
        host,
        _MILLER_FACETS.get(facet_str, 'fcc111'),
        0.0,
        (trace,),
        1,
        0,
        size=(lateral, lateral, 4),
        seed_key=genome,
    )
    loading_scale = 1.0 + min(float(loading), 10000.0) * 1.0e-7
    slab.set_cell(slab.cell * [loading_scale, loading_scale, 1.0], scale_atoms=True)
    return slab, top_idx


def _build_metal_free_carbon(genome, structure_seed):
    _, n_type, n_frac, defect, substrate, co_dop = genome
    cluster = generate_porphyrin_cluster(
        'N',
        _CARBON_N_CAVITIES.get(n_type, 'N4'),
        d_metal_n=1.8 + float(n_frac),
        substrate=substrate,
        linker=defect,
        pore_size=8.0 + 20.0 * float(n_frac),
    )
    if co_dop.lower() != 'none':
        _append_marker(cluster, co_dop, (2.5, 0.0, 0.5))
    return cluster, [0]


_STRUCTURE_BUILDERS = {
    'SupportedNanoparticle': _build_supported_nanoparticle,
    'SolidCatalyst': _build_solid_catalyst,
    'MoltenMetal': _build_molten_metal,
    'SAC': _build_sac,
    'DAC': _build_dac,
    'MOF': _build_framework,
    'COF': _build_framework,
    'Perovskite': _build_perovskite,
    'MetalHydride': _build_metal_hydride,
    'MAXPhase': _build_max_phase,
    'HEA': _build_hea,
    'Spinel': _build_spinel,
    'MXene': _build_mxene,
    'SAA': _build_saa,
    'MetalFreeCarbon': _build_metal_free_carbon,
}


def generate_structure(
    genome: tuple, structure_seed: int = 0
) -> Tuple[Atoms, list, str]:
    """Generate an atomic structure from a catalyst genome.

    Args:
        genome: Encoded catalyst composition and structural configuration.
        structure_seed: Deterministic structural replica index.

    Returns:
        Tuple of `(atoms, active_site_indices, material_class)`.
    """
    mat_class = genome[0]
    builder = _STRUCTURE_BUILDERS.get(mat_class) if isinstance(mat_class, str) else None
    if builder is None:
        raise ValueError(f"Unknown material class: {mat_class}")
    atoms, active_indices = builder(genome, structure_seed)
    return atoms, active_indices, mat_class
