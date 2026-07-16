#!/usr/bin/env python3
"""
Catalyst Design Space Definitions for Turquoise Hydrogen Production.

Defines the complete, physically grounded chemical space of catalysts
organized by material class. Each class has:
  - A genome specification (tuple structure for the genetic algorithm)
  - A random genome generator
  - A feature encoder for surrogate model training
  - Constraints that ensure only lab-synthesizable candidates are generated

Material Classes:
  A. MoltenMetal    — Low-melting-point alloys for bubble column reactors
  B. SolidCatalyst  — Supported metal catalysts for packed/fluidized bed reactors
  C. SAC / DAC      — Single/Dual-atom catalysts on N-doped carbon
  D. MOF / COF      — Porous framework catalysts
"""

import random
import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass, field

from pipeline.utils import (
    CRUSTAL_ABUNDANCE_PPM, MELTING_POINT_K, METAL_PRICE_USD_KG,
    is_molten_at_temperature
)


# ═══════════════════════════════════════════════════════════════════════════════
# A. MOLTEN METAL ALLOYS
# ═══════════════════════════════════════════════════════════════════════════════
# For bubble column reactors. The host metal must be liquid at operating temperature.
# A catalytic promoter (transition metal) dissolves into the melt at low concentration.

MOLTEN_HOSTS = ['Sn', 'Bi', 'In', 'Ga', 'Pb', 'Sb', 'Te']

MOLTEN_PROMOTERS = [
    'Ni', 'Cu', 'Fe', 'Co', 'Mn', 'Pd', 'Pt', 'Ru', 'Rh',
    'Mo', 'W', 'V', 'Cr', 'Ti', 'Zr', 'Nb', 'La', 'Ce',
    'None',  # Pure molten metal (no promoter)
]

MOLTEN_PROMOTER_AT_PCT = [0.0, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]

MOLTEN_TEMPERATURES_K = [
    700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200
]


def generate_molten_metal_genome() -> tuple:
    """Generate a random molten metal catalyst genome."""
    host = random.choice(MOLTEN_HOSTS)
    promoter = random.choice(MOLTEN_PROMOTERS)
    at_pct = random.choice(MOLTEN_PROMOTER_AT_PCT) if promoter != 'None' else 0.0
    temp = random.choice(MOLTEN_TEMPERATURES_K)
    # Ensure host is molten at operating temperature
    if not is_molten_at_temperature(host, temp):
        temp = max(MOLTEN_TEMPERATURES_K)  # push to highest T
    return ('MoltenMetal', host, promoter, at_pct, temp)


def validate_molten_metal(genome: tuple) -> bool:
    """Check if a molten metal genome is physically feasible."""
    _, host, promoter, at_pct, temp = genome
    if not is_molten_at_temperature(host, temp):
        return False
    if promoter != 'None' and at_pct <= 0:
        return False
    if promoter == 'None' and at_pct > 0:
        return False
    if at_pct > 30.0:
        return False
    # Promoter must not form immiscible phase at concentration
    # (simplified check — real phase diagrams are more complex)
    if promoter in ['W', 'Mo', 'Nb', 'Zr', 'Ti'] and at_pct > 10.0:
        return False  # refractory metals have limited solubility
    return True


# ═══════════════════════════════════════════════════════════════════════════════
# B. SOLID CATALYSTS ON SUPPORTS
# ═══════════════════════════════════════════════════════════════════════════════
# For packed-bed and fluidized-bed reactors.

SOLID_ACTIVE_METALS = ['Ni', 'Fe', 'Co', 'Cu', 'Mo', 'W', 'Mn', 'V', 'Cr']

SOLID_SUPPORTS = [
    'Al2O3', 'SiO2', 'MgO', 'CeO2', 'TiO2', 'ZrO2',
    'Carbon',  # Activated carbon / carbon black
    'Graphene',
]

SOLID_FACETS = ['fcc111', 'fcc100', 'bcc110', 'hcp0001']

SOLID_DOPANTS = [
    'Li', 'Na', 'Mg', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
    'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'In', 'Sn',
    'Sb', 'Hf', 'Ta', 'W', 'Re', 'Ir', 'Pt', 'Au', 'Pb', 'Bi', 'B', 'C',
    'Si', 'Ge', 'P', 'S', 'Ga', 'Al', 'La', 'Ce',
]

SOLID_STRAIN_RANGE = (-0.08, 0.08)  # ±8%


def generate_solid_catalyst_genome() -> tuple:
    """Generate a random solid catalyst genome."""
    metal = random.choice(SOLID_ACTIVE_METALS)
    support = random.choice(SOLID_SUPPORTS)
    facet = random.choice(SOLID_FACETS)
    strain = random.uniform(*SOLID_STRAIN_RANGE)
    num_dop = random.randint(0, 3)
    dopants = tuple(random.choice(SOLID_DOPANTS) for _ in range(num_dop))
    num_sub = random.randint(1, min(4, max(1, num_dop)))
    num_vac = random.randint(0, 2)
    return ('SolidCatalyst', metal, support, facet, strain, dopants, num_sub, num_vac)


# ═══════════════════════════════════════════════════════════════════════════════
# C. SINGLE-ATOM & DUAL-ATOM CATALYSTS (SAC / DAC)
# ═══════════════════════════════════════════════════════════════════════════════
# Metal atoms atomically dispersed in nitrogen-doped carbon matrices.

SAC_METALS = ['Fe', 'Co', 'Ni', 'Cu', 'Mn', 'Mo', 'W', 'Ru', 'Rh', 'V', 'Cr', 'Zn', 'Sn', 'In']

SAC_COORDINATIONS = [
    'N4', 'N3C', 'N2C2', 'N3B', 'N3S', 'N3P', 'N2S2', 'N2O2',
    'N3O', 'N4_pyridine', 'N4_pyrrole',
]

SAC_SUBSTRATES = ['N-graphene', 'N-CNT', 'N-carbon_black']


def generate_sac_genome() -> tuple:
    """Generate a random SAC genome."""
    metal = random.choice(SAC_METALS)
    coord = random.choice(SAC_COORDINATIONS)
    substrate = random.choice(SAC_SUBSTRATES)
    return ('SAC', metal, coord, substrate)


DAC_METALS_1 = ['Fe', 'Co', 'Ni', 'Cu', 'Mn', 'Mo', 'W']
DAC_METALS_2 = ['Fe', 'Co', 'Ni', 'Cu', 'Mn', 'Zn', 'Sn', 'In', 'V', 'Cr']

DAC_COORDINATIONS = ['N6', 'N8', 'N4C2', 'N3SN3', 'N4N4']


def generate_dac_genome() -> tuple:
    """Generate a random DAC genome."""
    m1 = random.choice(DAC_METALS_1)
    m2 = random.choice(DAC_METALS_2)
    coord = random.choice(DAC_COORDINATIONS)
    substrate = random.choice(SAC_SUBSTRATES)
    return ('DAC', m1, m2, coord, substrate)


# ═══════════════════════════════════════════════════════════════════════════════
# D. METAL-ORGANIC FRAMEWORKS (MOFs) & COVALENT ORGANIC FRAMEWORKS (COFs)
# ═══════════════════════════════════════════════════════════════════════════════

MOF_METAL_NODES = ['Fe', 'Co', 'Ni', 'Mn', 'Cu', 'Zn', 'Zr', 'Ti', 'V', 'Cr', 'Mo']

MOF_LINKERS = ['BDC', 'BTC', 'Porphyrin', 'Phthalocyanine', 'NDC', 'BPDC', 'Pyrazole']

MOF_CAVITIES = ['N4', 'N2S2', 'N2O2', 'O4', 'N3S', 'N3P', 'N2P2']

MOF_PORE_SIZES = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 25.0]

COF_LINKAGES = ['Imine', 'Triazine', 'Boroxine', 'Phenazine', 'Olefin', 'Azine']


def generate_mof_genome() -> tuple:
    """Generate a random MOF genome."""
    metal = random.choice(MOF_METAL_NODES)
    linker = random.choice(MOF_LINKERS)
    cavity = random.choice(MOF_CAVITIES)
    pore = random.choice(MOF_PORE_SIZES)
    return ('MOF', metal, linker, cavity, pore)


def generate_cof_genome() -> tuple:
    """Generate a random COF genome."""
    metal = random.choice(MOF_METAL_NODES + ['None'])
    linkage = random.choice(COF_LINKAGES)
    cavity = random.choice(MOF_CAVITIES)
    pore = random.choice(MOF_PORE_SIZES)
    return ('COF', metal, linkage, cavity, pore)


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

ALL_MATERIAL_CLASSES = ['MoltenMetal', 'SolidCatalyst', 'SAC', 'DAC', 'MOF', 'COF']

# Weights for random generation — bias toward industrially relevant classes
CLASS_WEIGHTS = {
    'MoltenMetal': 0.30,     # Primary focus for turquoise H₂
    'SolidCatalyst': 0.30,   # Well-established technology
    'SAC': 0.15,             # Emerging, high potential
    'DAC': 0.10,             # Cutting-edge
    'MOF': 0.10,             # Exploratory
    'COF': 0.05,             # Exploratory
}

GENERATORS = {
    'MoltenMetal': generate_molten_metal_genome,
    'SolidCatalyst': generate_solid_catalyst_genome,
    'SAC': generate_sac_genome,
    'DAC': generate_dac_genome,
    'MOF': generate_mof_genome,
    'COF': generate_cof_genome,
}


def generate_random_genome(material_class: Optional[str] = None) -> tuple:
    """
    Generate a random catalyst genome.
    If material_class is None, randomly select one based on CLASS_WEIGHTS.
    """
    if material_class is None:
        classes = list(CLASS_WEIGHTS.keys())
        weights = list(CLASS_WEIGHTS.values())
        material_class = random.choices(classes, weights=weights, k=1)[0]
    return GENERATORS[material_class]()


def generate_population(pop_size: int, material_class: Optional[str] = None) -> List[tuple]:
    """Generate a diverse initial population of catalyst genomes."""
    population = []
    for _ in range(pop_size):
        genome = generate_random_genome(material_class)
        population.append(genome)
    return population


# ═══════════════════════════════════════════════════════════════════════════════
# FEATURE ENCODING FOR SURROGATE MODEL
# ═══════════════════════════════════════════════════════════════════════════════

# Build lookup dictionaries for O(1) encoding
_ALL_METALS = sorted(set(
    MOLTEN_HOSTS + MOLTEN_PROMOTERS + SOLID_ACTIVE_METALS +
    SAC_METALS + DAC_METALS_1 + DAC_METALS_2 + MOF_METAL_NODES +
    ['None']
))
_METAL_IDX = {m: i for i, m in enumerate(_ALL_METALS)}

_ALL_SUPPORTS = sorted(set(SOLID_SUPPORTS + SAC_SUBSTRATES))
_SUPPORT_IDX = {s: i for i, s in enumerate(_ALL_SUPPORTS)}

_ALL_FACETS_LINKERS = sorted(set(
    SOLID_FACETS + MOF_LINKERS + COF_LINKAGES
))
_FACET_LINKER_IDX = {f: i for i, f in enumerate(_ALL_FACETS_LINKERS)}

_ALL_COORDS = sorted(set(SAC_COORDINATIONS + DAC_COORDINATIONS + MOF_CAVITIES))
_COORD_IDX = {c: i for i, c in enumerate(_ALL_COORDS)}

_ALL_DOPANTS = sorted(set(SOLID_DOPANTS))
_DOPANT_IDX = {d: i for i, d in enumerate(_ALL_DOPANTS)}

# Total feature vector dimension
N_CLASSES = len(ALL_MATERIAL_CLASSES)
N_METALS = len(_ALL_METALS)
N_SUPPORTS = len(_ALL_SUPPORTS)
N_FACETS = len(_ALL_FACETS_LINKERS)
N_COORDS = len(_ALL_COORDS)
N_DOPANTS = len(_ALL_DOPANTS)
N_CONTINUOUS = 4  # strain, at_pct/pore_size, num_sub, num_vac/temperature

FEATURE_DIM = N_CLASSES + 2 * N_METALS + N_SUPPORTS + N_FACETS + N_COORDS + N_DOPANTS + N_CONTINUOUS


def encode_genome(genome: tuple) -> np.ndarray:
    """
    Encode a catalyst genome into a fixed-length feature vector
    suitable for surrogate model training.
    
    Returns: np.ndarray of shape (FEATURE_DIM,)
    """
    mat_class = genome[0]

    # One-hot: material class
    cls_vec = np.zeros(N_CLASSES)
    cls_idx = ALL_MATERIAL_CLASSES.index(mat_class) if mat_class in ALL_MATERIAL_CLASSES else 0
    cls_vec[cls_idx] = 1.0

    # One-hot: primary metal
    metal1_vec = np.zeros(N_METALS)
    # One-hot: secondary metal (promoter / second DAC metal)
    metal2_vec = np.zeros(N_METALS)
    # One-hot: support / substrate
    support_vec = np.zeros(N_SUPPORTS)
    # One-hot: facet / linker
    facet_vec = np.zeros(N_FACETS)
    # One-hot: coordination / cavity
    coord_vec = np.zeros(N_COORDS)
    # Multi-hot: dopants
    dop_vec = np.zeros(N_DOPANTS)
    # Continuous features
    cont = np.zeros(N_CONTINUOUS)

    if mat_class == 'MoltenMetal':
        _, host, promoter, at_pct, temp = genome
        if host in _METAL_IDX:
            metal1_vec[_METAL_IDX[host]] = 1.0
        if promoter in _METAL_IDX:
            metal2_vec[_METAL_IDX[promoter]] = 1.0
        cont[0] = 0.0  # no strain concept
        cont[1] = at_pct / 30.0  # normalized
        cont[2] = 0.0
        cont[3] = (temp - 700.0) / 500.0  # normalized temperature

    elif mat_class == 'SolidCatalyst':
        _, metal, support, facet, strain, dopants, n_sub, n_vac = genome
        if metal in _METAL_IDX:
            metal1_vec[_METAL_IDX[metal]] = 1.0
        if support in _SUPPORT_IDX:
            support_vec[_SUPPORT_IDX[support]] = 1.0
        if facet in _FACET_LINKER_IDX:
            facet_vec[_FACET_LINKER_IDX[facet]] = 1.0
        for d in dopants:
            if d in _DOPANT_IDX:
                dop_vec[_DOPANT_IDX[d]] = 1.0
        cont[0] = strain / 0.08  # normalized
        cont[1] = 0.0
        cont[2] = n_sub / 4.0
        cont[3] = n_vac / 2.0

    elif mat_class == 'SAC':
        _, metal, coord, substrate = genome
        if metal in _METAL_IDX:
            metal1_vec[_METAL_IDX[metal]] = 1.0
        if coord in _COORD_IDX:
            coord_vec[_COORD_IDX[coord]] = 1.0
        if substrate in _SUPPORT_IDX:
            support_vec[_SUPPORT_IDX[substrate]] = 1.0

    elif mat_class == 'DAC':
        _, m1, m2, coord, substrate = genome
        if m1 in _METAL_IDX:
            metal1_vec[_METAL_IDX[m1]] = 1.0
        if m2 in _METAL_IDX:
            metal2_vec[_METAL_IDX[m2]] = 1.0
        if coord in _COORD_IDX:
            coord_vec[_COORD_IDX[coord]] = 1.0
        if substrate in _SUPPORT_IDX:
            support_vec[_SUPPORT_IDX[substrate]] = 1.0

    elif mat_class in ('MOF', 'COF'):
        _, metal, linker, cavity, pore = genome
        if metal in _METAL_IDX:
            metal1_vec[_METAL_IDX[metal]] = 1.0
        if linker in _FACET_LINKER_IDX:
            facet_vec[_FACET_LINKER_IDX[linker]] = 1.0
        if cavity in _COORD_IDX:
            coord_vec[_COORD_IDX[cavity]] = 1.0
        cont[1] = pore / 25.0  # normalized pore size

    return np.concatenate([
        cls_vec, metal1_vec, metal2_vec, support_vec,
        facet_vec, coord_vec, dop_vec, cont
    ])


def encode_population(population: List[tuple]) -> np.ndarray:
    """Encode a list of genomes into a feature matrix."""
    return np.stack([encode_genome(g) for g in population])


# ═══════════════════════════════════════════════════════════════════════════════
# CROSSOVER & MUTATION OPERATORS
# ═══════════════════════════════════════════════════════════════════════════════

def crossover(parent1: tuple, parent2: tuple) -> tuple:
    """
    Uniform crossover between two genomes of the same material class.
    If classes differ, randomly pick one parent's class and build child
    by mixing compatible genes.
    """
    if parent1[0] != parent2[0]:
        # Cross-class: just pick one parent (class-preserving)
        return parent1 if random.random() < 0.5 else parent2

    mat_class = parent1[0]
    child_genes = [mat_class]

    # Gene-wise crossover (skip index 0 which is the class tag)
    for i in range(1, min(len(parent1), len(parent2))):
        gene = parent1[i] if random.random() < 0.5 else parent2[i]
        child_genes.append(gene)

    return tuple(child_genes)


def mutate(genome: tuple, rate: float = 0.2) -> tuple:
    """
    Mutate a genome with probability `rate` per gene.
    Respects material-class-specific constraints.
    """
    if random.random() > rate:
        return genome

    mat_class = genome[0]
    genes = list(genome)

    if mat_class == 'MoltenMetal':
        gene_idx = random.randint(1, 4)
        if gene_idx == 1:
            genes[1] = random.choice(MOLTEN_HOSTS)
        elif gene_idx == 2:
            genes[2] = random.choice(MOLTEN_PROMOTERS)
        elif gene_idx == 3:
            genes[3] = random.choice(MOLTEN_PROMOTER_AT_PCT)
        elif gene_idx == 4:
            genes[4] = random.choice(MOLTEN_TEMPERATURES_K)

    elif mat_class == 'SolidCatalyst':
        gene_idx = random.randint(1, 7)
        if gene_idx == 1:
            genes[1] = random.choice(SOLID_ACTIVE_METALS)
        elif gene_idx == 2:
            genes[2] = random.choice(SOLID_SUPPORTS)
        elif gene_idx == 3:
            genes[3] = random.choice(SOLID_FACETS)
        elif gene_idx == 4:
            genes[4] = random.uniform(*SOLID_STRAIN_RANGE)
        elif gene_idx == 5:
            n = random.randint(0, 3)
            genes[5] = tuple(random.choice(SOLID_DOPANTS) for _ in range(n))
        elif gene_idx == 6:
            genes[6] = random.randint(1, 4)
        elif gene_idx == 7:
            genes[7] = random.randint(0, 2)

    elif mat_class == 'SAC':
        gene_idx = random.randint(1, 3)
        if gene_idx == 1:
            genes[1] = random.choice(SAC_METALS)
        elif gene_idx == 2:
            genes[2] = random.choice(SAC_COORDINATIONS)
        elif gene_idx == 3:
            genes[3] = random.choice(SAC_SUBSTRATES)

    elif mat_class == 'DAC':
        gene_idx = random.randint(1, 4)
        if gene_idx == 1:
            genes[1] = random.choice(DAC_METALS_1)
        elif gene_idx == 2:
            genes[2] = random.choice(DAC_METALS_2)
        elif gene_idx == 3:
            genes[3] = random.choice(DAC_COORDINATIONS)
        elif gene_idx == 4:
            genes[4] = random.choice(SAC_SUBSTRATES)

    elif mat_class in ('MOF', 'COF'):
        gene_idx = random.randint(1, 4)
        if gene_idx == 1:
            genes[1] = random.choice(MOF_METAL_NODES + (['None'] if mat_class == 'COF' else []))
        elif gene_idx == 2:
            genes[2] = random.choice(MOF_LINKERS if mat_class == 'MOF' else COF_LINKAGES)
        elif gene_idx == 3:
            genes[3] = random.choice(MOF_CAVITIES)
        elif gene_idx == 4:
            genes[4] = random.choice(MOF_PORE_SIZES)

    return tuple(genes)


# ═══════════════════════════════════════════════════════════════════════════════
# DESIGN SPACE STATISTICS
# ═══════════════════════════════════════════════════════════════════════════════

def estimate_design_space_size() -> Dict[str, int]:
    """Estimate the combinatorial size of each material class."""
    sizes = {
        'MoltenMetal': (
            len(MOLTEN_HOSTS) * len(MOLTEN_PROMOTERS) *
            len(MOLTEN_PROMOTER_AT_PCT) * len(MOLTEN_TEMPERATURES_K)
        ),
        'SolidCatalyst': (
            len(SOLID_ACTIVE_METALS) * len(SOLID_SUPPORTS) *
            len(SOLID_FACETS) * 16 *  # strain discretized to 16 bins
            (len(SOLID_DOPANTS) ** 2) * 4 * 3  # dopant combos × subs × vacs
        ),
        'SAC': len(SAC_METALS) * len(SAC_COORDINATIONS) * len(SAC_SUBSTRATES),
        'DAC': (
            len(DAC_METALS_1) * len(DAC_METALS_2) *
            len(DAC_COORDINATIONS) * len(SAC_SUBSTRATES)
        ),
        'MOF': (
            len(MOF_METAL_NODES) * len(MOF_LINKERS) *
            len(MOF_CAVITIES) * len(MOF_PORE_SIZES)
        ),
        'COF': (
            (len(MOF_METAL_NODES) + 1) * len(COF_LINKAGES) *
            len(MOF_CAVITIES) * len(MOF_PORE_SIZES)
        ),
    }
    sizes['TOTAL'] = sum(sizes.values())
    return sizes


if __name__ == '__main__':
    print("=" * 70)
    print("  CATALYST DESIGN SPACE SUMMARY")
    print("=" * 70)

    sizes = estimate_design_space_size()
    for cls, size in sizes.items():
        print(f"  {cls:20s}: {size:>15,} configurations")

    print(f"\n  Feature vector dimension: {FEATURE_DIM}")
    print(f"\n  Generating 5 random genomes:")
    for _ in range(5):
        g = generate_random_genome()
        print(f"    {g}")
