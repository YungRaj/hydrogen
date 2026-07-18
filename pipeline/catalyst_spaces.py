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
    is_molten_at_temperature, TOXIC_ELEMENTS,
)


def _safe(elements: list) -> list:
    """Remove toxic/radioactive elements from a design space list."""
    return [e for e in elements if e not in TOXIC_ELEMENTS]


# ═══════════════════════════════════════════════════════════════════════════════
# A. MOLTEN METAL ALLOYS
# ═══════════════════════════════════════════════════════════════════════════════
# For bubble column reactors. The host metal must be liquid at operating temperature.
# A catalytic promoter (transition metal) dissolves into the melt at low concentration.

MOLTEN_HOSTS = _safe([
    # Classic low-melting hosts
    'Sn', 'Bi', 'In', 'Ga', 'Pb', 'Sb', 'Te',
    # Extended hosts (higher melting but viable at elevated T)
    'Zn', 'Al', 'Ag', 'Au', 'Cd', 'Tl',
])

MOLTEN_PROMOTERS = _safe([
    # 3d transition metals
    'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    # 4d transition metals
    'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag',
    # 5d transition metals
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    # Rare earth metals
    'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Gd', 'Dy', 'Er', 'Yb',
    # Post-transition metals & metalloids
    'Al', 'Ga', 'Ge', 'Se', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Tl', 'Pb', 'Bi',
    # Light elements
    'Li', 'Na', 'K', 'Mg', 'Ca', 'Sr', 'Ba',
    'None',  # Pure molten metal (no promoter)
])

MOLTEN_PROMOTER_AT_PCT = [
    0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5,
    10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0,
]

MOLTEN_TEMPERATURES_K = [
    600, 650, 700, 750, 800, 850, 900, 950,
    1000, 1050, 1100, 1150, 1200, 1250, 1300,
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

SOLID_ACTIVE_METALS = _safe([
    # 3d TMs
    'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    # 4d TMs
    'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag',
    # 5d TMs
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    # Rare earths
    'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Gd',
    # Post-transition metals
    'Al', 'Ga', 'In', 'Sn', 'Sb', 'Bi', 'Pb',
])

SOLID_SUPPORTS = [
    # Oxides
    'Al2O3', 'SiO2', 'MgO', 'CeO2', 'TiO2', 'ZrO2',
    'Fe2O3', 'Cr2O3', 'V2O5', 'Nb2O5', 'La2O3', 'Y2O3',
    'WO3', 'MoO3', 'SnO2', 'Bi2O3', 'BaO', 'SrTiO3',
    # Mixed oxides
    'CeZrO2', 'MgAl2O4', 'BaTiO3', 'LaAlO3',
    # Carbon-based
    'Carbon', 'Graphene', 'CNT', 'Graphite', 'Diamond',
    'Fullerene', 'CarbonNitride',
    # Nitrides / Carbides / Borides
    'BN', 'Si3N4', 'AlN', 'TiN', 'WC', 'SiC', 'MoC', 'TiC', 'TaC',
    'TiB2', 'ZrB2', 'HfB2',
    # Sulfides / Phosphides
    'MoS2', 'WS2', 'NiP', 'FeP', 'CoP',
    # Zeolites
    'ZSM5', 'Beta', 'Mordenite', 'Y_zeolite', 'SAPO34',
]

SOLID_FACETS = [
    'fcc111', 'fcc100', 'fcc110', 'fcc211', 'fcc311',
    'bcc110', 'bcc100', 'bcc111', 'bcc211',
    'hcp0001', 'hcp1010', 'hcp1120',
]

SOLID_DOPANTS = _safe([
    # All accessible dopants in the periodic table
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
    'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag',
    'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    'Tl', 'Pb', 'Bi',
])

SOLID_STRAIN_RANGE = (-0.10, 0.10)  # ±10%


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

SAC_METALS = _safe([
    # All viable single-atom metals
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Gd',
    'Al', 'Ga', 'In', 'Sn', 'Sb', 'Bi', 'Pb',
])

SAC_COORDINATIONS = [
    'N4', 'N3C', 'N2C2', 'N3B', 'N3S', 'N3P', 'N2S2', 'N2O2',
    'N3O', 'N4_pyridine', 'N4_pyrrole', 'N2P2', 'N2Se2',
    'O4', 'S4', 'N2B2', 'C4', 'N3Se',
]

SAC_SUBSTRATES = [
    'N-graphene', 'N-CNT', 'N-carbon_black',
    'N-graphdiyne', 'N-porous_carbon', 'N-mesoporous_carbon',
    'BN_sheet', 'MoS2_sheet', 'Phosphorene',
]


def generate_sac_genome() -> tuple:
    """Generate a random SAC genome."""
    metal = random.choice(SAC_METALS)
    coord = random.choice(SAC_COORDINATIONS)
    substrate = random.choice(SAC_SUBSTRATES)
    return ('SAC', metal, coord, substrate)


DAC_METALS_1 = SAC_METALS  # Full metal set for dual-atom site 1
DAC_METALS_2 = SAC_METALS  # Full metal set for dual-atom site 2

DAC_COORDINATIONS = [
    'N6', 'N8', 'N4C2', 'N3SN3', 'N4N4', 'N6C2', 'N4S2',
    'N4O2', 'N4P2', 'N2S2N2', 'N3BN3', 'C4N4',
]


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

MOF_METAL_NODES = _safe([
    'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    'La', 'Ce', 'Nd', 'Gd', 'Al', 'Ga', 'In', 'Sn', 'Sb', 'Bi',
    'Mg', 'Ca', 'Sr', 'Ba', 'Cd', 'Pb',
])

MOF_LINKERS = [
    'BDC', 'BTC', 'Porphyrin', 'Phthalocyanine', 'NDC', 'BPDC', 'Pyrazole',
    'Imidazolate', 'Triazolate', 'Tetrazolate', 'Oxalate', 'Fumarate',
    'DOBDC', 'TCPP', 'HHTP', 'HITP', 'BTB',
]

MOF_CAVITIES = [
    'N4', 'N2S2', 'N2O2', 'O4', 'N3S', 'N3P', 'N2P2',
    'O6', 'S4', 'N6', 'N4O2', 'N2Se2', 'P4',
]

MOF_PORE_SIZES = [4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 25.0, 30.0, 40.0]

COF_LINKAGES = [
    'Imine', 'Triazine', 'Boroxine', 'Phenazine', 'Olefin', 'Azine',
    'Hydrazone', 'Squaraine', 'Beta_ketoenamine', 'Dioxin', 'Amide', 'Imide',
]


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
# E. PEROVSKITE CATALYSTS (ABO₃)
# ═══════════════════════════════════════════════════════════════════════════════

PEROVSKITE_A_SITE = _safe([
    'La', 'Sr', 'Ba', 'Ca', 'Pr', 'Nd', 'Sm', 'Gd', 'Y', 'Ce',
    'K', 'Na', 'Li', 'Bi', 'Pb', 'Ag',
])
PEROVSKITE_B_SITE = _safe([
    'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd',
    'Hf', 'Ta', 'W', 'Re', 'Ir', 'Pt',
    'Al', 'Ga', 'In', 'Sn', 'Sb',
])
PEROVSKITE_DOPANT_FRAC = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50]
PEROVSKITE_DEFECTS = ['none', 'A_vacancy', 'B_vacancy', 'O_vacancy', 'A_excess']

def generate_perovskite_genome() -> tuple:
    A = random.choice(PEROVSKITE_A_SITE)
    B = random.choice(PEROVSKITE_B_SITE)
    dopant = random.choice(PEROVSKITE_B_SITE + ['None'])
    frac = random.choice(PEROVSKITE_DOPANT_FRAC) if dopant != 'None' else 0.0
    defect = random.choice(PEROVSKITE_DEFECTS)
    return ('Perovskite', A, B, dopant, frac, defect)


# ═══════════════════════════════════════════════════════════════════════════════
# F. METAL HYDRIDES (for H₂ storage & catalytic decomposition)
# ═══════════════════════════════════════════════════════════════════════════════

HYDRIDE_METALS = _safe([
    # Simple hydrides
    'Li', 'Na', 'K', 'Mg', 'Ca', 'Sr', 'Ba', 'Al', 'Ti', 'Zr', 'Hf',
    'V', 'Nb', 'Ta', 'La', 'Ce', 'Pr', 'Nd', 'Y', 'Sc',
    # Intermetallic hydrides
    'Fe', 'Co', 'Ni', 'Cu', 'Mn', 'Pd', 'Pt',
])
HYDRIDE_TYPES = [
    'simple', 'complex_alanate', 'complex_borohydride', 'complex_amide',
    'intermetallic_AB5', 'intermetallic_AB2', 'intermetallic_AB',
    'intermetallic_A2B', 'perovskite_hydride',
]
HYDRIDE_SECOND_METAL = HYDRIDE_METALS + ['None']
HYDRIDE_ADDITIVES = _safe([
    'None', 'TiCl3', 'TiF3', 'VCl3', 'NbF5', 'ZrCl4', 'CeO2',
    'MgO', 'Carbon', 'Graphene', 'CNT', 'TiO2', 'Fe2O3',
])

def generate_hydride_genome() -> tuple:
    metal = random.choice(HYDRIDE_METALS)
    h_type = random.choice(HYDRIDE_TYPES)
    second = random.choice(HYDRIDE_SECOND_METAL)
    additive = random.choice(HYDRIDE_ADDITIVES)
    temp = random.choice([300, 350, 400, 450, 500, 550, 600, 700, 800])
    return ('MetalHydride', metal, h_type, second, additive, temp)


# ═══════════════════════════════════════════════════════════════════════════════
# G. MAX PHASES (M_{n+1}AX_n) — layered ternary carbides/nitrides
# ═══════════════════════════════════════════════════════════════════════════════

MAX_M_ELEMENTS = _safe([
    'Ti', 'V', 'Cr', 'Zr', 'Nb', 'Mo', 'Hf', 'Ta', 'W',
    'Sc', 'Mn', 'Fe', 'Co', 'Ni',
])
MAX_A_ELEMENTS = _safe([
    'Al', 'Si', 'P', 'S', 'Ga', 'Ge', 'As', 'Cd', 'In', 'Sn',
    'Tl', 'Pb', 'Bi',
])
MAX_X_ELEMENTS = ['C', 'N']  # carbide or nitride
MAX_N_VALUES = [1, 2, 3]  # n in M_{n+1}AX_n → 211, 312, 413 phases

def generate_max_genome() -> tuple:
    M = random.choice(MAX_M_ELEMENTS)
    A = random.choice(MAX_A_ELEMENTS)
    X = random.choice(MAX_X_ELEMENTS)
    n = random.choice(MAX_N_VALUES)
    dopant = random.choice(MAX_M_ELEMENTS + ['None'])
    facet = random.choice(['basal_0001', 'edge_1010', 'edge_1120'])
    return ('MAXPhase', M, A, X, n, dopant, facet)


# ═══════════════════════════════════════════════════════════════════════════════
# H. HIGH-ENTROPY ALLOYS (HEAs) — 4-6 component equimolar alloys
# ═══════════════════════════════════════════════════════════════════════════════

HEA_ELEMENTS = _safe([
    'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    'Al', 'Si', 'Ga', 'Ge', 'Sn', 'Sb', 'Bi',
    'La', 'Ce', 'Y', 'Sc',
])
HEA_STRUCTURES = ['fcc', 'bcc', 'hcp', 'fcc_bcc_dual', 'amorphous']

def generate_hea_genome() -> tuple:
    n_components = random.randint(4, 6)
    components = tuple(sorted(random.sample(HEA_ELEMENTS, n_components)))
    structure = random.choice(HEA_STRUCTURES)
    facet = random.choice(['111', '100', '110', '211'])
    temp = random.choice([800, 900, 1000, 1100, 1200, 1300])
    return ('HEA', components, structure, facet, temp)


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

ALL_MATERIAL_CLASSES = [
    'MoltenMetal', 'SolidCatalyst', 'SAC', 'DAC',
    'MOF', 'COF', 'Perovskite', 'MetalHydride', 'MAXPhase', 'HEA',
]

# Weights for random generation — balanced across all classes
CLASS_WEIGHTS = {
    'MoltenMetal': 0.15,     # Primary focus for turquoise H₂
    'SolidCatalyst': 0.15,   # Well-established technology
    'SAC': 0.12,             # Emerging, high potential
    'DAC': 0.10,             # Cutting-edge
    'MOF': 0.10,             # Porous frameworks
    'COF': 0.08,             # Covalent frameworks
    'Perovskite': 0.10,      # ABO₃ oxides
    'MetalHydride': 0.08,    # H₂ storage & decomposition
    'MAXPhase': 0.06,        # Layered ternary ceramics
    'HEA': 0.06,             # High-entropy alloys
}

GENERATORS = {
    'MoltenMetal': generate_molten_metal_genome,
    'SolidCatalyst': generate_solid_catalyst_genome,
    'SAC': generate_sac_genome,
    'DAC': generate_dac_genome,
    'MOF': generate_mof_genome,
    'COF': generate_cof_genome,
    'Perovskite': generate_perovskite_genome,
    'MetalHydride': generate_hydride_genome,
    'MAXPhase': generate_max_genome,
    'HEA': generate_hea_genome,
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
    PEROVSKITE_A_SITE + PEROVSKITE_B_SITE +
    HYDRIDE_METALS + MAX_M_ELEMENTS + MAX_A_ELEMENTS +
    HEA_ELEMENTS + ['None']
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

    elif mat_class == 'Perovskite':
        _, A, B, dopant, frac, defect = genome
        if A in _METAL_IDX:
            metal1_vec[_METAL_IDX[A]] = 1.0
        if B in _METAL_IDX:
            metal2_vec[_METAL_IDX[B]] = 1.0
        if dopant in _DOPANT_IDX:
            dop_vec[_DOPANT_IDX[dopant]] = 1.0
        cont[0] = frac / 0.50  # normalized dopant fraction
        # Encode defect type as continuous
        _defect_map = {'none': 0.0, 'A_vacancy': 0.2, 'B_vacancy': 0.4,
                       'O_vacancy': 0.6, 'A_excess': 0.8}
        cont[1] = _defect_map.get(defect, 0.0)

    elif mat_class == 'MetalHydride':
        _, metal, h_type, second, additive, temp = genome
        if metal in _METAL_IDX:
            metal1_vec[_METAL_IDX[metal]] = 1.0
        if second in _METAL_IDX:
            metal2_vec[_METAL_IDX[second]] = 1.0
        # Encode hydride type as normalized index
        _hydride_types = [
            'simple', 'complex_alanate', 'complex_borohydride', 'complex_amide',
            'intermetallic_AB5', 'intermetallic_AB2', 'intermetallic_AB',
            'intermetallic_A2B', 'perovskite_hydride',
        ]
        cont[0] = _hydride_types.index(h_type) / len(_hydride_types) if h_type in _hydride_types else 0.0
        cont[3] = (temp - 400.0) / 400.0  # normalized temperature

    elif mat_class == 'MAXPhase':
        _, M, A, X, n, dopant, facet = genome
        if M in _METAL_IDX:
            metal1_vec[_METAL_IDX[M]] = 1.0
        if A in _METAL_IDX:
            metal2_vec[_METAL_IDX[A]] = 1.0
        if dopant in _DOPANT_IDX:
            dop_vec[_DOPANT_IDX[dopant]] = 1.0
        cont[0] = 1.0 if X == 'N' else 0.0  # carbide vs nitride
        cont[1] = n / 3.0  # normalized phase order
        # Encode facet
        _max_facet_map = {'basal_0001': 0.0, 'edge_1010': 0.5, 'edge_1120': 1.0}
        cont[2] = _max_facet_map.get(facet, 0.0)

    elif mat_class == 'HEA':
        _, components, structure, facet, temp = genome
        # Multi-hot encode all component elements
        for comp in components:
            if comp in _METAL_IDX:
                metal1_vec[_METAL_IDX[comp]] = 1.0
        # Encode structure as normalized index
        _hea_structs = ['fcc', 'bcc', 'hcp', 'fcc_bcc_dual', 'amorphous']
        cont[0] = _hea_structs.index(structure) / len(_hea_structs) if structure in _hea_structs else 0.0
        # Encode facet
        _hea_facet_map = {'111': 0.0, '100': 0.33, '110': 0.67, '211': 1.0}
        cont[1] = _hea_facet_map.get(facet, 0.0)
        cont[3] = (temp - 900.0) / 400.0  # normalized temperature

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

    elif mat_class == 'Perovskite':
        gene_idx = random.randint(1, 5)
        if gene_idx == 1:
            genes[1] = random.choice(PEROVSKITE_A_SITE)
        elif gene_idx == 2:
            genes[2] = random.choice(PEROVSKITE_B_SITE)
        elif gene_idx == 3:
            genes[3] = random.choice(PEROVSKITE_B_SITE + ['None'])
        elif gene_idx == 4:
            genes[4] = random.choice(PEROVSKITE_DOPANT_FRAC)
        elif gene_idx == 5:
            genes[5] = random.choice(PEROVSKITE_DEFECTS)

    elif mat_class == 'MetalHydride':
        gene_idx = random.randint(1, 5)
        if gene_idx == 1:
            genes[1] = random.choice(HYDRIDE_METALS)
        elif gene_idx == 2:
            genes[2] = random.choice(HYDRIDE_TYPES)
        elif gene_idx == 3:
            genes[3] = random.choice(HYDRIDE_SECOND_METAL)
        elif gene_idx == 4:
            genes[4] = random.choice(HYDRIDE_ADDITIVES)
        elif gene_idx == 5:
            genes[5] = random.choice([300, 350, 400, 450, 500, 550, 600, 700, 800])

    elif mat_class == 'MAXPhase':
        gene_idx = random.randint(1, 6)
        if gene_idx == 1:
            genes[1] = random.choice(MAX_M_ELEMENTS)
        elif gene_idx == 2:
            genes[2] = random.choice(MAX_A_ELEMENTS)
        elif gene_idx == 3:
            genes[3] = random.choice(MAX_X_ELEMENTS)
        elif gene_idx == 4:
            genes[4] = random.choice(MAX_N_VALUES)
        elif gene_idx == 5:
            genes[5] = random.choice(MAX_M_ELEMENTS + ['None'])
        elif gene_idx == 6:
            genes[6] = random.choice(['basal_0001', 'edge_1010', 'edge_1120'])

    elif mat_class == 'HEA':
        gene_idx = random.randint(1, 4)
        if gene_idx == 1:
            n = random.randint(4, 6)
            genes[1] = tuple(sorted(random.sample(HEA_ELEMENTS, n)))
        elif gene_idx == 2:
            genes[2] = random.choice(HEA_STRUCTURES)
        elif gene_idx == 3:
            genes[3] = random.choice(['111', '100', '110', '211'])
        elif gene_idx == 4:
            genes[4] = random.choice([800, 900, 1000, 1100, 1200, 1300])

    return tuple(genes)


# ═══════════════════════════════════════════════════════════════════════════════
# DESIGN SPACE STATISTICS
# ═══════════════════════════════════════════════════════════════════════════════

def estimate_design_space_size() -> Dict[str, int]:
    """Estimate the combinatorial size of each material class."""
    from math import comb
    sizes = {
        'MoltenMetal': (
            len(MOLTEN_HOSTS) * len(MOLTEN_PROMOTERS) *
            len(MOLTEN_PROMOTER_AT_PCT) * len(MOLTEN_TEMPERATURES_K)
        ),
        'SolidCatalyst': (
            len(SOLID_ACTIVE_METALS) * len(SOLID_SUPPORTS) *
            len(SOLID_FACETS) * 20 *  # strain discretized to 20 bins
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
        'Perovskite': (
            len(PEROVSKITE_A_SITE) * len(PEROVSKITE_B_SITE) *
            (len(PEROVSKITE_B_SITE) + 1) * len(PEROVSKITE_DOPANT_FRAC) *
            len(PEROVSKITE_DEFECTS)
        ),
        'MetalHydride': (
            len(HYDRIDE_METALS) * len(HYDRIDE_TYPES) *
            len(HYDRIDE_SECOND_METAL) * len(HYDRIDE_ADDITIVES) * 9  # 9 temps
        ),
        'MAXPhase': (
            len(MAX_M_ELEMENTS) * len(MAX_A_ELEMENTS) *
            len(MAX_X_ELEMENTS) * len(MAX_N_VALUES) *
            (len(MAX_M_ELEMENTS) + 1) * 3  # dopant + facets
        ),
        'HEA': (
            # C(35,4) + C(35,5) + C(35,6) component combos × structure × facet × temp
            (comb(len(HEA_ELEMENTS), 4) + comb(len(HEA_ELEMENTS), 5) +
             comb(len(HEA_ELEMENTS), 6)) *
            len(HEA_STRUCTURES) * 4 * 6
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
