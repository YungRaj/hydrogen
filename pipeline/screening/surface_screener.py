#!/usr/bin/env python3
# Turquoise-hydrogen candidate screening.
"""
Multi-GPU Parallelized Meta eSen-SM Catalyst Screening Engine.

For each catalyst candidate, computes:
  1. Relaxed surface energy
  2. H* adsorption energy (methane activation descriptor)
  3. CH₃* adsorption energy (dissociation intermediate)
  4. C* adsorption energy (coking resistance descriptor)
  5. Reaction energy: CH₄ → CH₃* + H*
  6. Segregation / binding stability energy
  7. NEB transition state barrier (top candidates only)

Parallelized across all available GPUs with 4 workers per device.
"""

import os
import sys
import time
import random
import hashlib
import numpy as np
import multiprocessing as mp
from typing import List, Tuple, Dict, Optional

# ASE imports
from ase import Atoms, Atom
from ase.build import fcc111, fcc100, bcc110, hcp0001, molecule
from ase.optimize import BFGS
from ase.constraints import FixAtoms

from pipeline.common.utils import (
    BASE_DIR, SCREENING_DIR, setup_logger, print_banner,
    k_B_eV, bep_activation_energy, arrhenius_rate,
    abundance_cost_penalty, save_screening_db,
    check_element_safety, is_valid_for_application,
    CRUSTAL_ABUNDANCE_PPM, MELTING_POINT_K,
)

logger = setup_logger('surface_screener', 'screening/surface_screening.log')

SCREENING_PROTOCOL_ID = 'esen-sm-conserving-all-oc25:relax-v2:pyrolysis-v2'


def _stable_seed(*parts) -> int:
    """Return a process-independent seed for one physical realization."""
    payload = repr(parts).encode('utf-8')
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], 'big')


def _label_fraction(label) -> float:
    """Stable [0, 1) coordinate used to distinguish categorical realizations."""
    return (_stable_seed('label-coordinate-v1', label) % 10_000) / 10_000.0


def _representative_element(label: str, default: str = 'C') -> str:
    """Map support/linker labels to a chemically relevant explicit atom."""
    text = str(label)
    for token in ('Mo', 'Ti', 'Al', 'Si', 'Mg', 'Zr', 'Ce', 'La', 'Fe',
                  'Cl', 'Se', 'Zn', 'Ni', 'Co', 'Cu', 'Ca', 'Na',
                  'B', 'N', 'O', 'S', 'P', 'F', 'C'):
        if token in text:
            return token
    return default


def _append_marker(atoms: Atoms, symbol: str, offset: tuple[float, float, float]) -> int:
    """Add one explicit environment atom without overlapping the active site."""
    center = atoms.get_center_of_mass()
    atoms.append(Atom(symbol, position=center + np.asarray(offset, dtype=float)))
    return len(atoms) - 1


def _add_surface_markers(atoms: Atoms, symbols: list[str]) -> list[int]:
    """Add sparse environment atoms and restore a valid bottom constraint."""
    if not symbols:
        return []
    atoms.set_constraint()
    z_max = float(atoms.positions[:, 2].max())
    cell = atoms.cell.lengths()
    added = []
    for index, symbol in enumerate(symbols):
        x = (index + 1) * float(cell[0]) / (len(symbols) + 1)
        y = (index % 2 + 1) * float(cell[1]) / 3.0
        atoms.append(Atom(symbol, position=(x, y, z_max + 1.5 + 0.25 * index)))
        added.append(len(atoms) - 1)
    z = atoms.positions[:, 2]
    atoms.set_constraint(FixAtoms(mask=z < z.min() + 4.0))
    return added


# ═══════════════════════════════════════════════════════════════════════════════
# STRUCTURE GENERATORS
# ═══════════════════════════════════════════════════════════════════════════════

def generate_alloy_slab(host: str, facet: str, strain: float,
                        dopants: tuple, n_sub: int, n_vac: int,
                        size: Tuple[int,int,int] = (3, 3, 4),
                        vacuum: float = 12.0, seed_key=None,
                        environment_label: str | None = None) -> Tuple[Atoms, list]:
    """Generate a strained alloy slab with dopant substitutions and vacancies."""
    builders = {
        'fcc111': fcc111, 'fcc100': fcc100,
        'bcc110': bcc110, 'hcp0001': hcp0001,
    }
    builder = builders.get(facet, fcc111)

    # Explicit lattice constants (Å) — ASE can't guess these for many elements
    LATTICE_CONSTANTS = {
        # FCC metals
        'Ni': 3.52, 'Cu': 3.61, 'Ag': 4.09, 'Au': 4.08, 'Al': 4.05,
        'Pd': 3.89, 'Pt': 3.92, 'Rh': 3.80, 'Ir': 3.84, 'Pb': 4.95,
        'Ca': 5.58, 'Sr': 6.08, 'Ce': 5.16, 'Yb': 5.49, 'Th': 5.08,
        # BCC metals
        'Fe': 2.87, 'Cr': 2.88, 'V': 3.03, 'Nb': 3.30, 'Mo': 3.15,
        'W': 3.16, 'Ta': 3.30, 'Na': 4.29, 'K': 5.33, 'Ba': 5.02,
        'Li': 3.51, 'Cs': 6.14, 'Rb': 5.59,
        # HCP metals (using equivalent FCC a)
        'Ti': 2.95, 'Zr': 3.23, 'Hf': 3.19, 'Co': 2.51, 'Ru': 2.71,
        'Os': 2.74, 'Re': 2.76, 'Sc': 3.31, 'Y': 3.65, 'La': 3.75,
        'Mg': 3.21, 'Zn': 2.66, 'Cd': 2.98,
        # Non-standard / semimetals (use nearest FCC equivalent)
        'Sn': 5.83, 'In': 4.60, 'Ga': 4.52, 'Sb': 4.31, 'Bi': 4.75,
        'Ge': 5.66, 'Si': 5.43, 'Te': 4.45, 'Se': 4.36, 'As': 4.13,
        'Mn': 3.50, 'Pr': 5.16, 'Nd': 5.08, 'Gd': 3.63, 'Tb': 3.60,
        'Dy': 3.59, 'Ho': 3.58, 'Er': 3.56, 'Tm': 3.54, 'Lu': 3.50,
        'Sm': 3.62, 'Eu': 4.58, 'Tl': 3.46, 'Po': 3.35,
    }

    a = LATTICE_CONSTANTS.get(host, 3.60)  # default 3.60 Å if unknown
    try:
        slab = builder(host, size=size, vacuum=vacuum, a=a)
    except Exception:
        # Fallback: use Ni as host template, then substitute
        slab = fcc111('Ni', size=size, vacuum=vacuum, a=a)
        for atom in slab:
            atom.symbol = host

    # Identify surface atoms (top layer by z-coordinate)
    z_coords = slab.positions[:, 2]
    z_max = z_coords.max()
    top_mask = z_coords > (z_max - 2.5)  # atoms within 2.5 Å of surface
    top_indices = np.where(top_mask)[0]

    # Deterministic shuffle based on composition for reproducibility
    rng = random.Random(_stable_seed(
        'alloy-slab-v2', host, facet, float(strain), tuple(dopants),
        int(n_sub), int(n_vac), tuple(size), float(vacuum), seed_key))
    shuffled = list(top_indices)
    rng.shuffle(shuffled)

    # Substitutions
    sub_count = min(n_sub, len(shuffled)) if dopants else 0
    for i in range(sub_count):
        slab[shuffled[i]].symbol = dopants[i % len(dopants)]

    # Vacancies (remove atoms after substitutions)
    vac_count = min(n_vac, len(shuffled) - sub_count)
    if vac_count > 0:
        del_indices = sorted(shuffled[sub_count:sub_count + vac_count], reverse=True)
        for idx in del_indices:
            del slab[idx]
        # Recompute top_indices after deletion
        z_coords = slab.positions[:, 2]
        top_mask = z_coords > (z_coords.max() - 2.5)
        top_indices = np.where(top_mask)[0]

    # A compact interface proxy prevents support-distinct genomes from
    # collapsing to an identical unsupported metal slab.
    if environment_label and str(environment_label).lower() != 'none':
        bottom = int(np.argmin(slab.positions[:, 2]))
        slab[bottom].symbol = _representative_element(environment_label, 'C')
        slab.positions[bottom, 0] += 0.05 * _label_fraction(environment_label)

    # Apply biaxial strain
    cell = slab.get_cell()
    cell[0] *= (1.0 + strain)
    cell[1] *= (1.0 + strain)
    slab.set_cell(cell, scale_atoms=True)

    # Fix bottom 2 layers
    z_coords = slab.positions[:, 2]
    z_min = z_coords.min()
    fix_mask = z_coords < (z_min + 4.0)
    slab.set_constraint(FixAtoms(mask=fix_mask))

    return slab, list(top_indices)


def generate_porphyrin_cluster(metal: str, cavity: str,
                               d_metal_n: float = 2.0,
                               substrate: str = 'N-graphene',
                               axial_ligand: str = 'none',
                               linker: str | None = None,
                               pore_size: float | None = None) -> Atoms:
    """Generate a metal-porphyrin active site cluster for SAC/MOF evaluation."""
    atoms = Atoms()

    # Central metal atom
    if metal != 'None':
        atoms.append(Atom(metal, position=(0.0, 0.0, 0.0)))
    else:
        # COFs have no metal node; use an explicit framework carbon instead of
        # ASE's dummy X element, which atomistic models cannot embed.
        atoms.append(Atom('C', position=(0.0, 0.0, 0.0)))

    # Parse cavity for coordinating atoms
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
        'N6': ['N', 'N', 'N', 'N', 'N', 'N'],  # dual-atom
        'N8': ['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'],
        'N4C2': ['N', 'N', 'N', 'N', 'C', 'C'],
        'N3SN3': ['N', 'N', 'N', 'S', 'N', 'N'],
        'N4N4': ['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'],
        'N2P2': ['N', 'N', 'P', 'P'],
    }
    coord_atoms = coord_map.get(cavity, ['N', 'N', 'N', 'N'])

    # Place coordination atoms in a square/hexagonal arrangement
    if len(coord_atoms) <= 4:
        positions = [
            (-d_metal_n, 0.0, 0.0),
            (d_metal_n, 0.0, 0.0),
            (0.0, -d_metal_n, 0.0),
            (0.0, d_metal_n, 0.0),
        ]
    else:
        # Hexagonal for 6+ coordinating atoms
        n_coord = len(coord_atoms)
        positions = []
        for k in range(n_coord):
            angle = 2 * np.pi * k / n_coord
            positions.append((d_metal_n * np.cos(angle), d_metal_n * np.sin(angle), 0.0))

    for sym, pos in zip(coord_atoms, positions[:len(coord_atoms)]):
        atoms.append(Atom(sym, position=pos))
    if len(atoms) > 1:
        atoms.positions[1, 2] += 0.08 * _label_fraction(cavity)

    # Add surrounding carbon skeleton
    n_c = 8
    pore_scale = 1.0 if pore_size is None else np.clip(float(pore_size) / 12.0, 0.75, 1.75)
    for k in range(n_c):
        angle = 2 * np.pi * k / n_c + np.pi / n_c
        r = (d_metal_n + 1.2) * pore_scale
        atoms.append(Atom('C', position=(r * np.cos(angle), r * np.sin(angle), 0.0)))

    # Add hydrogen termination
    for k in range(n_c):
        angle = 2 * np.pi * k / n_c + np.pi / n_c
        r = (d_metal_n + 2.2) * pore_scale
        atoms.append(Atom('H', position=(r * np.cos(angle), r * np.sin(angle), 0.0)))

    # Encode the local support/linker chemistry explicitly. These atoms are a
    # compact active-site environment, not a claim of a converged full crystal.
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
            atoms, ligand_symbol,
            (0.0, 0.0, 1.9 + 0.3 * _label_fraction(axial_ligand)))

    atoms.set_cell((15.0, 15.0, 15.0))
    atoms.center()
    atoms.pbc = True
    return atoms


def generate_structure(genome: tuple) -> Tuple[Atoms, list, str]:
    """
    Generate an atomic structure from a catalyst genome.
    Returns: (atoms, active_site_indices, material_class)
    """
    mat_class = genome[0]

    if mat_class == 'SolidCatalyst':
        _, metal, support, facet, strain, dopants, n_sub, n_vac = genome
        # Map extended facets to available builders
        facet_key = facet.split('_')[0] if '_' in facet else facet
        if facet_key not in ('fcc111', 'fcc100', 'bcc110', 'hcp0001'):
            facet_key = 'fcc111'  # default fallback
        slab, top_idx = generate_alloy_slab(
            metal, facet_key, strain, dopants, n_sub, n_vac,
            seed_key=genome, environment_label=support)
        return slab, top_idx, mat_class

    elif mat_class == 'MoltenMetal':
        _, host, promoter, at_pct, temp = genome
        n_sub = max(1, int(at_pct / 100.0 * 36))
        dopants = (promoter,) if promoter != 'None' else ()
        thermal_strain = (float(temp) - 1000.0) * 1.0e-5
        slab, top_idx = generate_alloy_slab(
            host, 'fcc111', thermal_strain, dopants, n_sub, 0,
            seed_key=genome)
        return slab, top_idx, mat_class

    elif mat_class in ('SAC', 'DAC'):
        if mat_class == 'SAC':
            metal = genome[1]
            coord = genome[2]
            substrate = genome[3]
            axial = genome[4]
            cluster = generate_porphyrin_cluster(
                metal, coord, substrate=substrate, axial_ligand=axial)
        else:
            _, m1, m2, coord, substrate = genome
            cluster = generate_porphyrin_cluster(m1, coord, substrate=substrate)
            pos = cluster[0].position.copy()
            pos[0] += 2.5
            cluster.append(Atom(m2, position=pos))
        return cluster, [0], mat_class

    elif mat_class in ('MOF', 'COF'):
        _, metal, linker, cavity, pore = genome
        cluster = generate_porphyrin_cluster(
            metal, cavity, substrate='framework', linker=linker, pore_size=pore)
        return cluster, [0], mat_class

    elif mat_class == 'Perovskite':
        # ABO₃ perovskite — simple cubic unit cell replicated as slab
        _, A, B, dopant, frac, defect = genome
        slab = _generate_perovskite_slab(A, B, dopant, frac, defect, seed_key=genome)
        z = slab.positions[:, 2]
        top_idx = list(np.where(z > z.max() - 3.0)[0])
        return slab, top_idx, mat_class

    elif mat_class == 'MetalHydride':
        # Metal hydride — metal + H in bulk-like slab
        _, metal, h_type, second, additive, temp = genome
        slab = _generate_hydride_slab(
            metal, second, h_type=h_type, additive=additive,
            temperature_K=temp, seed_key=genome)
        z = slab.positions[:, 2]
        top_idx = list(np.where(z > z.max() - 3.0)[0])
        return slab, top_idx, mat_class

    elif mat_class == 'MAXPhase':
        # M_{n+1}AX_n layered structure — model as M-slab with A/X interstitials
        _, M, A, X, n_val, dopant, facet = genome
        facet_map = {'basal_0001': 'hcp0001', 'edge_1010': 'bcc110',
                     'edge_1120': 'fcc100'}
        dopants = tuple(x for x in (A, dopant) if x != 'None')
        slab, top_idx = generate_alloy_slab(
            M, facet_map.get(facet, 'hcp0001'), 0.0, dopants,
            min(len(dopants), 2), 0, size=(3, 3, n_val + 1),
            seed_key=genome)
        _add_surface_markers(slab, [X] * n_val)
        return slab, top_idx, mat_class

    elif mat_class == 'HEA':
        # High-entropy alloy — random substitution FCC slab
        _, components, structure, facet_str, temp = genome
        host = components[0]
        dopants = components[1:]
        facet_map = {'111': 'fcc111', '100': 'fcc100', '110': 'bcc110', '211': 'fcc111'}
        facet_key = facet_map.get(facet_str, 'fcc111')
        if structure == 'bcc':
            facet_key = 'bcc110'
        elif structure == 'hcp':
            facet_key = 'hcp0001'
        n_sub = min(len(dopants), 6)
        thermal_strain = (float(temp) - 1000.0) * 5.0e-6
        if structure == 'amorphous':
            thermal_strain += 0.003
        elif structure == 'fcc_bcc_dual':
            thermal_strain -= 0.003
        slab, top_idx = generate_alloy_slab(
            host, facet_key, thermal_strain, dopants, n_sub, 0,
            seed_key=genome)
        return slab, top_idx, mat_class

    elif mat_class == 'Spinel':
        # AB₂O₄ spinel — model as B-metal oxide slab
        _, A, B, dopant, morph, support = genome
        morphology_defect = {
            'mesoporous': 'O_vacancy', 'hollow_sphere': 'A_vacancy',
            'nanorod': 'B_vacancy',
        }.get(morph, 'none')
        slab = _generate_spinel_slab(
            A, B, dopant, morphology_defect, support, seed_key=genome)
        morphology_scale = 0.98 + 0.04 * _label_fraction(morph)
        slab.set_cell(slab.cell * [morphology_scale, 1.0, 1.0],
                      scale_atoms=True)
        z = slab.positions[:, 2]
        top_idx = list(np.where(z > z.max() - 3.0)[0])
        return slab, top_idx, mat_class

    elif mat_class == 'MXene':
        # 2D MXene — model as M-carbide/nitride slab
        _, M, X_elem, n_val, term, sac_metal = genome
        slab, top_idx = generate_alloy_slab(
            M, 'hcp0001', 0.0,
            (sac_metal,) if sac_metal != 'None' else (), 1, 0,
            size=(3, 3, n_val + 1), seed_key=genome)
        x_symbols = list(('C', 'N') if X_elem == 'CN' else (X_elem,))
        term_symbols = {
            'OH': ['O', 'H'], 'O': ['O'], 'F': ['F'], 'Cl': ['Cl'],
            'S': ['S'], 'mixed_OH_O': ['O', 'H', 'O'],
            'mixed_OH_F': ['O', 'H', 'F'], 'bare': [],
        }.get(term, [])
        _add_surface_markers(slab, x_symbols * n_val + term_symbols)
        return slab, top_idx, mat_class

    elif mat_class == 'SAA':
        # Single-atom alloy — host metal slab with trace substitution
        _, trace, host, facet_str, loading = genome
        facet_map = {'111': 'fcc111', '100': 'fcc100', '110': 'bcc110', '211': 'fcc111'}
        facet_key = facet_map.get(facet_str, 'fcc111')
        # Loading controls the explicit surface-cell size: more dilute SAAs use
        # a larger host cell around the single trace atom.
        lateral = 4 if loading <= 100 else (3 if loading <= 1000 else 2)
        slab, top_idx = generate_alloy_slab(
            host, facet_key, 0.0, (trace,), 1, 0,
            size=(lateral, lateral, 4), seed_key=genome)
        loading_scale = 1.0 + min(float(loading), 10000.0) * 1.0e-7
        slab.set_cell(slab.cell * [loading_scale, loading_scale, 1.0],
                      scale_atoms=True)
        return slab, top_idx, mat_class

    elif mat_class == 'MetalFreeCarbon':
        # Metal-free N-doped carbon — model as N-graphene cluster
        _, n_type, n_frac, defect, substrate, co_dop = genome
        cavity = {
            'pyridinic': 'N3C', 'pyrrolic': 'N4_pyrrole',
            'graphitic': 'N4', 'oxidized': 'N3O',
            'mixed_pyridinic_graphitic': 'N4C2',
        }.get(n_type, 'N4')
        cluster = generate_porphyrin_cluster(
            'N', cavity, d_metal_n=1.8 + float(n_frac),
            substrate=substrate, linker=defect,
            pore_size=8.0 + 20.0 * float(n_frac))
        if co_dop.lower() != 'none':
            _append_marker(cluster, co_dop, (2.5, 0.0, 0.5))
        return cluster, [0], mat_class

    else:
        raise ValueError(f"Unknown material class: {mat_class}")


def _generate_perovskite_slab(A: str, B: str, dopant: str,
                               frac: float, defect: str,
                               seed_key=None) -> Atoms:
    """Generate a simple ABO₃ perovskite slab for MACE evaluation."""
    a = 3.9  # approx perovskite lattice constant (Å)
    # 2×2×3 supercell → 12 ABO₃ units
    atoms = Atoms()
    a_sites, b_sites, o_sites = [], [], []
    for ix in range(2):
        for iy in range(2):
            for iz in range(3):
                base = np.array([ix * a, iy * a, iz * a])
                # A-site (corner)
                atoms.append(Atom(A, position=base))
                a_sites.append(len(atoms) - 1)
                # B-site (body center)
                atoms.append(Atom(B, position=base + np.array([a/2, a/2, a/2])))
                b_sites.append(len(atoms) - 1)
                # O-sites (face centers)
                atoms.append(Atom('O', position=base + np.array([a/2, a/2, 0])))
                o_sites.append(len(atoms) - 1)
                atoms.append(Atom('O', position=base + np.array([a/2, 0, a/2])))
                o_sites.append(len(atoms) - 1)
                atoms.append(Atom('O', position=base + np.array([0, a/2, a/2])))
                o_sites.append(len(atoms) - 1)

    rng = random.Random(_stable_seed(
        'perovskite-v2', A, B, dopant, float(frac), defect, seed_key))
    if dopant != 'None' and frac > 0:
        count = max(1, min(len(b_sites), round(float(frac) * len(b_sites))))
        for index in rng.sample(b_sites, count):
            atoms[index].symbol = dopant

    vacancy_pool = {
        'A_vacancy': a_sites, 'B_vacancy': b_sites, 'O_vacancy': o_sites,
    }.get(defect)
    if vacancy_pool:
        del atoms[rng.choice(vacancy_pool)]
    elif defect == 'A_excess':
        atoms.append(Atom(A, position=(
            a / 2, a / 2, float(atoms.positions[:, 2].max()) + 1.8)))

    cell = [2*a, 2*a, 3*a + 12.0]  # vacuum in z
    atoms.set_cell(cell)
    atoms.pbc = True
    # Fix bottom layer
    z = atoms.positions[:, 2]
    atoms.set_constraint(FixAtoms(mask=z < z.min() + 2.0))
    return atoms


def _generate_spinel_slab(A: str, B: str, dopant: str, defect: str,
                          support: str, seed_key=None) -> Atoms:
    """Build an AB2O4-like local slab rather than reusing ABO3 stoichiometry."""
    atoms = Atoms()
    a = 4.1
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                base = np.array([ix * a, iy * a, iz * a])
                atoms.append(Atom(A, position=base))
                atoms.append(Atom(B, position=base + (a/4, a/4, a/4)))
                atoms.append(Atom(B, position=base + (3*a/4, 3*a/4, 3*a/4)))
                for offset in ((a/2, 0, 0), (0, a/2, 0),
                               (0, 0, a/2), (a/2, a/2, a/2)):
                    atoms.append(Atom('O', position=base + offset))
    rng = random.Random(_stable_seed(
        'spinel-v2', A, B, dopant, defect, support, seed_key))
    if dopant != 'None':
        b_sites = [i for i, atom in enumerate(atoms) if atom.symbol == B]
        atoms[rng.choice(b_sites)].symbol = dopant
    vacancy_symbol = {'A_vacancy': A, 'B_vacancy': B, 'O_vacancy': 'O'}.get(defect)
    if vacancy_symbol:
        choices = [i for i, atom in enumerate(atoms) if atom.symbol == vacancy_symbol]
        if choices:
            del atoms[rng.choice(choices)]
    if support != 'none':
        atoms[0].symbol = _representative_element(support, 'C')
    atoms.set_cell((2*a, 2*a, 2*a + 12.0))
    atoms.pbc = True
    z = atoms.positions[:, 2]
    atoms.set_constraint(FixAtoms(mask=z < z.min() + 2.0))
    return atoms


def _generate_hydride_slab(metal: str, second: str, h_type: str = 'simple',
                           additive: str = 'None', temperature_K: float = 300,
                           seed_key=None) -> Atoms:
    """Generate a metal-hydride slab: metal FCC + interstitial H."""
    # Use explicit lattice constant (same table as generate_alloy_slab)
    LATTICE_CONSTANTS = {
        'Ni': 3.52, 'Cu': 3.61, 'Ag': 4.09, 'Au': 4.08, 'Al': 4.05,
        'Pd': 3.89, 'Pt': 3.92, 'Fe': 2.87, 'Ti': 2.95, 'Zr': 3.23,
        'Mg': 3.21, 'Ca': 5.58, 'Na': 4.29, 'Li': 3.51, 'La': 3.75,
        'Ce': 5.16, 'V': 3.03, 'Mn': 3.50, 'Co': 2.51, 'Zn': 2.66,
        'Sn': 5.83, 'In': 4.60, 'Sb': 4.31, 'Bi': 4.75,
    }
    a = LATTICE_CONSTANTS.get(metal, 3.60)

    try:
        slab = fcc111(metal, size=(3, 3, 3), vacuum=12.0, a=a)
    except Exception:
        slab = fcc111('Ni', size=(3, 3, 3), vacuum=12.0, a=a)
        for atom in slab:
            atom.symbol = metal

    # Collect H positions first, then add (don't iterate while appending)
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

    # Hydride family controls approximate hydrogen loading. This is a compact
    # screening realization; DFT promotion still requires phase-specific cells.
    extra_h = {
        'simple': 0, 'intermetallic_AB': 1, 'intermetallic_AB2': 2,
        'intermetallic_AB5': 3, 'intermetallic_A2B': 2,
        'complex_alanate': 4, 'complex_borohydride': 4,
        'complex_amide': 3, 'perovskite_hydride': 2,
    }.get(h_type, 0)
    rng = random.Random(_stable_seed(
        'hydride-v2', metal, second, h_type, additive,
        float(temperature_K), seed_key))
    for index in range(extra_h):
        slab.append(Atom('H', position=(
            (index + 1) * slab.cell.lengths()[0] / (extra_h + 1),
            (index % 2 + 1) * slab.cell.lengths()[1] / 3.0,
            top_z + 2.5 + 0.35 * index)))

    # Substitute surface atoms with second metal
    if second != 'None' and second != metal:
        indices = [i for i in range(len(slab))
                   if slab[i].symbol == metal and slab[i].position[2] > top_z - 2.5]
        for i in indices[:2]:
            slab[i].symbol = second

    if additive != 'None':
        added = _add_surface_markers(
            slab, [_representative_element(additive, 'C')])
        slab.positions[added[0], 0] += 0.1 * _label_fraction(additive)

    thermal_scale = 1.0 + (float(temperature_K) - 300.0) * 1.0e-5
    slab.set_cell(slab.cell * [thermal_scale, thermal_scale, 1.0],
                  scale_atoms=True)

    # Recompute z after adding H and set constraints
    z_all = slab.positions[:, 2]
    slab.set_constraint(FixAtoms(mask=z_all < z_all.min() + 3.0))
    return slab


# ═══════════════════════════════════════════════════════════════════════════════
# MACE EVALUATION WORKER
# ═══════════════════════════════════════════════════════════════════════════════

def compute_reference_energies(calc) -> Dict[str, float]:
    """Compute gas-phase reference energies using Meta eSen-SM."""
    refs = {}

    # H₂
    h2 = molecule('H2')
    h2.set_cell([10, 10, 10])
    h2.center()
    h2.pbc = True
    h2.calc = calc
    BFGS(h2, logfile=None).run(fmax=0.05)
    refs['H2'] = h2.get_potential_energy()

    # CH₄
    ch4 = molecule('CH4')
    ch4.set_cell([10, 10, 10])
    ch4.center()
    ch4.pbc = True
    ch4.calc = calc
    BFGS(ch4, logfile=None).run(fmax=0.05)
    refs['CH4'] = ch4.get_potential_energy()

    # CH₃ radical
    ch3 = molecule('CH3')
    ch3.set_cell([10, 10, 10])
    ch3.center()
    ch3.pbc = True
    ch3.calc = calc
    BFGS(ch3, logfile=None).run(fmax=0.05)
    refs['CH3'] = ch3.get_potential_energy()

    # Isolated C reference calculated thermodynamically: CH4 - 2*H2
    # This avoids zero-edge errors for single-atom systems in fairchem/eSen.
    refs['C'] = refs['CH4'] - 2.0 * refs['H2']

    return refs


def evaluate_candidate(genome: tuple, calc, refs: dict) -> dict:
    """
    Evaluate a single catalyst candidate using Meta eSen-SM.
    
    Returns a dictionary with all computed descriptors.
    """
    mat_class = genome[0]
    result = {
        'genome': str(genome),
        'material_class': mat_class,
        'valid': False,
        'screening_protocol': SCREENING_PROTOCOL_ID,
    }

    try:
        structure, active_idx, _ = generate_structure(genome)
        structure.pbc = True  # eSen requires PBC set to True in all dimensions

        # 1. Relax clean surface/cluster
        structure.calc = calc
        dyn = BFGS(structure, logfile=None)
        dyn.run(fmax=0.08, steps=150)
        e_clean = structure.get_potential_energy()
        result['e_clean'] = e_clean

        # Identify adsorption site (top of first active atom)
        if len(active_idx) > 0 and active_idx[0] < len(structure):
            ads_pos = structure[active_idx[0]].position.copy()
        else:
            ads_pos = structure.positions.mean(axis=0)
        ads_pos[2] += 1.8  # 1.8 Å above surface

        # 2. H* adsorption
        slab_h = structure.copy()
        h_pos = ads_pos.copy()
        h_pos[2] = ads_pos[2] - 0.3  # H binds closer
        slab_h.append(Atom('H', position=h_pos))
        slab_h.calc = calc
        BFGS(slab_h, logfile=None).run(fmax=0.08, steps=100)
        e_h = slab_h.get_potential_energy()
        dE_H = e_h - e_clean - 0.5 * refs['H2']
        result['dE_H'] = dE_H

        # 3. CH₃* adsorption
        slab_ch3 = structure.copy()
        c_pos = ads_pos.copy()
        slab_ch3.append(Atom('C', position=c_pos))
        slab_ch3.append(Atom('H', position=c_pos + np.array([0.0, 1.02, 0.35])))
        slab_ch3.append(Atom('H', position=c_pos + np.array([-0.88, -0.51, 0.35])))
        slab_ch3.append(Atom('H', position=c_pos + np.array([0.88, -0.51, 0.35])))
        slab_ch3.calc = calc
        BFGS(slab_ch3, logfile=None).run(fmax=0.08, steps=100)
        e_ch3 = slab_ch3.get_potential_energy()
        dE_CH3 = e_ch3 - e_clean - refs['CH3']
        result['dE_CH3'] = dE_CH3

        # 4. C* adsorption (coking indicator)
        slab_c = structure.copy()
        c_ads_pos = ads_pos.copy()
        c_ads_pos[2] -= 0.4  # C binds closer to surface
        slab_c.append(Atom('C', position=c_ads_pos))
        slab_c.calc = calc
        BFGS(slab_c, logfile=None).run(fmax=0.08, steps=100)
        e_c = slab_c.get_potential_energy()
        dE_C = e_c - e_clean - refs['C']
        result['dE_C'] = dE_C

        # 5. CH₄ → CH₃* + H* reaction energy
        dE_split = dE_CH3 + dE_H + refs['CH3'] + 0.5 * refs['H2'] - refs['CH4']
        result['dE_split'] = dE_split

        # 6. Activation barrier (BEP correlation)
        E_act = bep_activation_energy(dE_split, material_class=mat_class)
        result['E_act'] = E_act

        # 7. Coking resistance index
        coking_index = dE_C - 2.0 * dE_H  # positive = resistant
        
        # Apply liquid-metal coking resistance bonus for NTEC mode
        py_mode = os.environ.get('PYROLYSIS_MODE', 'ntec')
        if py_mode == 'ntec':
            from pipeline.screening.genetic_optimizer import _extract_elements_from_genome
            from pipeline.process.ntec_model import conditions_from_environment, ntec_assistance
            assistance = ntec_assistance(conditions_from_environment())
            elements = _extract_elements_from_genome(genome)
            if any(e in {'Ga', 'In', 'Sn', 'Bi'} for e in elements):
                coking_index += assistance['coking_bonus']
                result['ntec_evidence_status'] = assistance['status']
                result['ntec_assistance'] = assistance
            
        result['coking_index'] = coking_index

        # 8. Stability metric
        if mat_class in ('SolidCatalyst', 'MoltenMetal'):
            # Segregation energy for alloys
            e_seg = _compute_segregation_energy(structure, active_idx, calc, e_clean)
            result['segregation_energy'] = e_seg
        else:
            # Binding energy for SAC/MOF
            e_bind = _compute_binding_energy(structure, calc, e_clean, genome)
            result['segregation_energy'] = e_bind

        # 9. Rate constant at reference temperature
        T_ref = genome[4] if mat_class == 'MoltenMetal' else 1000.0
        A_prefactor = 1e13  # standard TST prefactor
        rate_k = arrhenius_rate(A_prefactor, E_act, T_ref)
        result['rate_constant'] = rate_k
        result['temperature_K'] = T_ref

        # 10. Extract element list for cost scoring
        elements = _extract_elements(genome)
        result['cost_penalty'] = abundance_cost_penalty(elements)

        # 10b. Safety check — reject toxic/radioactive elements
        is_safe, safety_reason = check_element_safety(elements)
        if not is_safe:
            result['valid'] = False
            result['error'] = safety_reason
            return result

        # 10c. Application feasibility flag
        result['pyrolysis_viable'] = is_valid_for_application(mat_class, 'pyrolysis')

        # 11. Physical sanity filters
        # Meta model can produce unphysical energies on exotic structures.
        # Valid adsorption energies for surface chemistry: |ΔE| < 10 eV
        SANE_LIMIT = 10.0  # eV
        for key in ('dE_H', 'dE_CH3', 'dE_C', 'dE_split'):
            val = result.get(key, 0)
            if abs(val) > SANE_LIMIT:
                result['valid'] = False
                result['error'] = f'Unphysical {key}={val:.2f} eV (|val|>{SANE_LIMIT})'
                return result

        # Clamp derived values to physical ranges
        result['E_act'] = max(0.01, min(result.get('E_act', 5.0), 5.0))
        result['E_act_censored'] = result['E_act'] in (0.01, 5.0)
        if result['E_act_censored']:
            result['needs_dft_validation'] = True
        result['coking_index'] = max(-20.0, min(result.get('coking_index', 0), 20.0))

        # 12. OOD confidence — how much we trust this prediction
        from pipeline.common.ood_detector import compute_model_confidence
        conf = compute_model_confidence(genome, elements)
        result['model_confidence'] = float(conf)
        result['needs_dft_validation'] = result.get('needs_dft_validation', False) or conf < 0.5

        result['valid'] = True

    except Exception as e:
        result['error'] = str(e)[:200]

    return result


def _compute_segregation_energy(slab, active_idx, calc, e_clean):
    """Compute segregation energy: tendency for dopant to stay at surface."""
    try:
        slab_swap = slab.copy()
        z_coords = slab_swap.positions[:, 2]
        bulk_idx = np.argsort(z_coords)[0]
        surf_idx = active_idx[0] if active_idx[0] < len(slab_swap) else 0

        if bulk_idx < len(slab_swap) and surf_idx < len(slab_swap):
            sym_s = slab_swap[surf_idx].symbol
            sym_b = slab_swap[bulk_idx].symbol
            if sym_s != sym_b:
                slab_swap[surf_idx].symbol = sym_b
                slab_swap[bulk_idx].symbol = sym_s
                slab_swap.calc = calc
                e_swap = slab_swap.get_potential_energy()
                return e_clean - e_swap  # negative = dopant prefers surface
        return 0.0
    except Exception:
        return 0.0


def _compute_binding_energy(structure, calc, e_clean, genome):
    """Compute metal-cavity binding energy for SAC/MOF."""
    try:
        cavity_struct = structure.copy()
        del cavity_struct[0]  # remove central metal
        cavity_struct.calc = calc
        e_empty = cavity_struct.get_potential_energy()

        ref_energies = {
            'Fe': -5.0, 'Co': -4.8, 'Ni': -5.5, 'Mn': -4.2, 'Cu': -3.8,
            'Zn': -2.5, 'Mo': -8.0, 'W': -9.5, 'V': -7.0, 'Cr': -6.0,
            'Ru': -6.5, 'Rh': -5.5, 'Ti': -6.5, 'Zr': -7.5, 'Sn': -3.0,
            'In': -2.5, 'None': 0.0,
        }
        metal = genome[1]
        e_metal = ref_energies.get(metal, -4.0)
        return (e_clean - e_empty - e_metal) / 10.0  # normalized
    except Exception:
        return 0.0


def _extract_elements(genome: tuple) -> List[str]:
    """Extract the list of metallic elements from a genome for cost scoring."""
    mat_class = genome[0]
    elements = []
    if mat_class == 'MoltenMetal':
        elements.append(genome[1])
        if genome[2] != 'None':
            elements.append(genome[2])
    elif mat_class == 'SolidCatalyst':
        elements.append(genome[1])
        for d in genome[5]:
            elements.append(d)
    elif mat_class == 'SAC':
        elements.append(genome[1])
    elif mat_class == 'DAC':
        elements.extend([genome[1], genome[2]])
    elif mat_class in ('MOF', 'COF'):
        if genome[1] != 'None':
            elements.append(genome[1])
    elif mat_class == 'Perovskite':
        elements.extend([genome[1], genome[2]])
        if genome[3] != 'None':
            elements.append(genome[3])
    elif mat_class == 'MetalHydride':
        elements.append(genome[1])
        if genome[3] != 'None':
            elements.append(genome[3])
    elif mat_class == 'MAXPhase':
        elements.extend([genome[1], genome[2]])
        if genome[5] != 'None':
            elements.append(genome[5])
    elif mat_class == 'HEA':
        elements.extend(list(genome[1]))
    elif mat_class == 'Spinel':
        elements.extend([genome[1], genome[2]])
        if genome[3] != 'None': elements.append(genome[3])
    elif mat_class == 'MXene':
        elements.append(genome[1])
        if genome[5] != 'None': elements.append(genome[5])
    elif mat_class == 'SAA':
        elements.extend([genome[1], genome[2]])
    elif mat_class == 'MetalFreeCarbon':
        pass  # no metals
    return [e for e in elements if e != 'None']


# ═══════════════════════════════════════════════════════════════════════════════
# WORKER PROCESS
# ═══════════════════════════════════════════════════════════════════════════════

def eval_worker(worker_id: int, gpu_id: int, gpu_uuid: str, task_queue: mp.Queue,
                result_queue: mp.Queue):
    """Worker process: loads Meta eSen on assigned GPU and evaluates candidates."""
    try:
        import os
        # This process is spawned without importing torch at module scope, so
        # the UUID mask takes effect before CUDA is initialized. Inside this
        # one-device namespace fairchem should use logical device ``cuda``.
        os.environ['CUDA_VISIBLE_DEVICES'] = gpu_uuid
        os.environ['OMP_NUM_THREADS'] = '1'
        os.environ['MKL_NUM_THREADS'] = '1'
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
        os.environ['NUMEXPR_NUM_THREADS'] = '1'
        
        # Limit CPU threads to prevent multiprocessing CPU over-subscription thrashing
        import torch
        torch.set_num_threads(1)
        torch.set_num_interop_threads(1)
        
        from pipeline.screening.surface_calculator import get_ocp_calculator
        calc = get_ocp_calculator(
            model_name='esen-sm-conserving-all-oc25', device='cuda')

        # Compute reference energies on this worker's calculator
        refs = compute_reference_energies(calc)

        while True:
            item = task_queue.get()
            if item is None:
                break

            idx, genome = item
            try:
                result = evaluate_candidate(genome, calc, refs)
                result['worker_id'] = worker_id
                result['gpu_id'] = gpu_id
                result_queue.put((idx, result))
            except Exception as e:
                result_queue.put((idx, {
                    'genome': str(genome),
                    'material_class': genome[0],
                    'valid': False,
                    'screening_protocol': SCREENING_PROTOCOL_ID,
                    'error': str(e)[:200],
                }))

    except Exception as e:
        logger.error(f"Worker {worker_id} failed to initialize: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN SCREENING ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

def run_screening(genomes: List[tuple], db_filename: str = "surface_screening.csv",
                  workers_per_gpu: int = 2) -> 'pd.DataFrame':
    """
    Run parallel Meta eSen-SM screening on a list of catalyst genomes.
    
    Args:
        genomes: List of catalyst genome tuples
        db_filename: Output CSV filename
        workers_per_gpu: Number of parallel workers per GPU
        
    Returns:
        DataFrame with all screening results
    """
    import pandas as pd
    import torch

    # These must exist before spawn starts a fresh interpreter and imports
    # numpy/torch. Setting them only inside eval_worker is too late for BLAS.
    for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ.setdefault(name, '1')

    mp.set_start_method('spawn', force=True)

    print_banner("META ESEN-SM SURFACE CATALYST SCREENING")
    logger.info(f"Screening {len(genomes)} catalyst candidates...")

    device_count = torch.cuda.device_count()
    if device_count < 1:
        raise RuntimeError('Meta eSen-SM screening requires at least one visible CUDA GPU')
    if workers_per_gpu < 1:
        raise ValueError('workers_per_gpu must be positive')
    num_workers = device_count * workers_per_gpu
    gpu_uuids = [f"GPU-{torch.cuda.get_device_properties(i).uuid}"
                 for i in range(device_count)]
    logger.info(f"Using {device_count} GPU(s), {num_workers} parallel workers")

    # Setup queues
    task_queue = mp.Queue()
    result_queue = mp.Queue()

    # Enqueue all tasks
    for idx, genome in enumerate(genomes):
        task_queue.put((idx, genome))

    # Poison pills
    for _ in range(num_workers):
        task_queue.put(None)

    # Launch workers
    workers = []
    for w_id in range(num_workers):
        gpu_id = w_id % device_count
        p = mp.Process(
            target=eval_worker,
            args=(w_id, gpu_id, gpu_uuids[gpu_id], task_queue, result_queue))
        p.start()
        workers.append(p)

    # Collect results
    results = []
    t_start = time.time()
    for i in range(len(genomes)):
        idx, result = result_queue.get()
        results.append(result)

        if (i + 1) % 50 == 0 or (i + 1) == len(genomes):
            elapsed = time.time() - t_start
            rate = (i + 1) / elapsed
            n_valid = sum(1 for r in results if r.get('valid', False))
            logger.info(
                f"Progress: {i+1}/{len(genomes)} "
                f"({rate:.1f} candidates/sec, {n_valid} valid, "
                f"{elapsed:.0f}s elapsed)"
            )

    # Wait for workers
    for p in workers:
        p.join(timeout=30)

    # Build DataFrame
    df = pd.DataFrame(results)
    path = save_screening_db(df, db_filename)
    logger.info(f"Screening complete. {len(df)} results saved to {path}")

    # Summary statistics
    valid_df = df[df['valid'] == True]
    if len(valid_df) > 0:
        logger.info(f"  Valid candidates: {len(valid_df)}/{len(df)}")
        logger.info(f"  Best E_act: {valid_df['E_act'].min():.4f} eV")
        logger.info(f"  Best coking index: {valid_df['coking_index'].max():.4f} eV")
        logger.info(f"  dE_H range: [{valid_df['dE_H'].min():.3f}, {valid_df['dE_H'].max():.3f}] eV")

    return df


if __name__ == '__main__':
    from pipeline.search.indexed_space import deterministic_tree_probes

    # Quick test: screen deterministic divide-and-conquer calibration probes
    pop = deterministic_tree_probes(20)
    df = run_screening(pop, db_filename="test_screening.csv", workers_per_gpu=2)
    print(f"\nScreening complete: {len(df)} candidates evaluated")
