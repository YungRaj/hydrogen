#!/usr/bin/env python3
# Shared catalyst design-space definitions.
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

from pipeline.common.utils import (
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
    # In-plane coordination environments
    'N4', 'N3C', 'N2C2', 'N3B', 'N3S', 'N3P', 'N2S2', 'N2O2',
    'N3O', 'N4_pyridine', 'N4_pyrrole', 'N2P2', 'N2Se2',
    'O4', 'S4', 'N2B2', 'C4', 'N3Se',
]

# Axial ligand variants — 2026 "dual-modulation" breakthrough
# (Jaouen et al., July 2026: in-plane + axial coordination decoration)
# Axial ligands modify the d-band center and spin state of the metal center
SAC_AXIAL_LIGANDS = [
    'none',       # bare M-Nx (standard)
    'OH',         # hydroxyl — E₁/₂ = 0.91V for Fe-N4-OH (2026 record)
    'O2',         # superoxo — stabilizes high-spin Fe³⁺
    'Cl',         # chloride — from pyrolysis precursor
    'NH3',        # amino — from ammonia treatment
    'CO',         # carbonyl — from CO₂ activation
    'H2O',        # aqua — operando coordination
    'pyridine',   # N-heterocycle — second coordination shell effect
    'imidazole',  # histidine-like — biomimetic coordination
]

SAC_SUBSTRATES = [
    'N-graphene', 'N-CNT', 'N-carbon_black',
    'N-graphdiyne', 'N-porous_carbon', 'N-mesoporous_carbon',
    'BN_sheet', 'MoS2_sheet', 'Phosphorene',
]


def generate_sac_genome() -> tuple:
    """Generate a random SAC genome.
    
    Genome: ('SAC', metal, coordination, substrate, axial_ligand)
    """
    metal = random.choice(SAC_METALS)
    coord = random.choice(SAC_COORDINATIONS)
    substrate = random.choice(SAC_SUBSTRATES)
    axial = random.choice(SAC_AXIAL_LIGANDS)
    return ('SAC', metal, coord, substrate, axial)


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
# I. SPINEL OXIDES (AB₂O₄) — bifunctional ORR/OER, earth-abundant
# ═══════════════════════════════════════════════════════════════════════════════

SPINEL_A_METALS = _safe(['Co', 'Mn', 'Ni', 'Fe', 'Cu', 'Zn', 'Mg', 'Li'])
SPINEL_B_METALS = _safe(['Co', 'Mn', 'Fe', 'Cr', 'Al', 'Ti', 'V', 'Ni', 'Cu'])
SPINEL_DOPANTS = _safe(['Ce', 'La', 'Y', 'Nb', 'Mo', 'W', 'Zr', 'None'])
SPINEL_MORPHOLOGIES = ['nanoparticle', 'nanorod', 'nanosheet', 'mesoporous', 'hollow_sphere']
SPINEL_SUPPORT_CARBONS = ['N-graphene', 'N-CNT', 'carbon_black', 'graphene_oxide', 'none']

def generate_spinel_genome() -> tuple:
    """Generate a random Spinel (AB₂O₄) genome.
    
    Genome: ('Spinel', A_metal, B_metal, dopant, morphology, support)
    """
    A = random.choice(SPINEL_A_METALS)
    B = random.choice(SPINEL_B_METALS)
    dopant = random.choice(SPINEL_DOPANTS)
    morph = random.choice(SPINEL_MORPHOLOGIES)
    support = random.choice(SPINEL_SUPPORT_CARBONS)
    return ('Spinel', A, B, dopant, morph, support)


# ═══════════════════════════════════════════════════════════════════════════════
# J. MXENES (Mn+1XnTx) — 2D carbides/nitrides with surface terminations
# ═══════════════════════════════════════════════════════════════════════════════

MXENE_M_ELEMENTS = _safe(['Ti', 'V', 'Cr', 'Mo', 'Nb', 'Ta', 'Zr', 'Hf', 'W', 'Mn'])
MXENE_X_ELEMENTS = ['C', 'N', 'CN']  # carbide, nitride, carbonitride
MXENE_TERMINATIONS = ['OH', 'O', 'F', 'Cl', 'S', 'mixed_OH_O', 'mixed_OH_F', 'bare']
MXENE_N_VALUES = [1, 2, 3]  # M2X, M3X2, M4X3
MXENE_SAC_METALS = _safe(['Fe', 'Co', 'Ni', 'Cu', 'Mn', 'Pt', 'Pd', 'Ru', 'None'])

def generate_mxene_genome() -> tuple:
    """Generate a random MXene genome.
    
    Genome: ('MXene', M_element, X_element, n, termination, sac_metal)
    MXenes are 2D sheets derived from MAX phases with distinct surface chemistry.
    Optionally decorated with single-atom metal sites.
    """
    M = random.choice(MXENE_M_ELEMENTS)
    X = random.choice(MXENE_X_ELEMENTS)
    n = random.choice(MXENE_N_VALUES)
    term = random.choice(MXENE_TERMINATIONS)
    sac = random.choice(MXENE_SAC_METALS)
    return ('MXene', M, X, n, term, sac)


# ═══════════════════════════════════════════════════════════════════════════════
# K. SINGLE-ATOM ALLOYS (SAAs) — trace PGM in abundant host
# ═══════════════════════════════════════════════════════════════════════════════

SAA_TRACE_METALS = _safe(['Pt', 'Pd', 'Rh', 'Ir', 'Ru', 'Au', 'Ag'])  # dispersed single atoms
SAA_HOST_METALS = _safe(['Cu', 'Ni', 'Co', 'Fe', 'Ag', 'Au', 'Sn', 'In', 'Ga', 'Al'])
SAA_FACETS = ['111', '100', '110', '211']
SAA_LOADINGS_PPM = [100, 500, 1000, 2000, 5000, 10000]  # trace metal loading

def generate_saa_genome() -> tuple:
    """Generate a random Single-Atom Alloy genome.
    
    Genome: ('SAA', trace_metal, host_metal, facet, loading_ppm)
    SAAs disperse isolated atoms of one metal in a host of another.
    Distinct from SACs (which use N-carbon supports).
    """
    trace = random.choice(SAA_TRACE_METALS)
    host = random.choice(SAA_HOST_METALS)
    facet = random.choice(SAA_FACETS)
    loading = random.choice(SAA_LOADINGS_PPM)
    return ('SAA', trace, host, facet, loading)


# ═══════════════════════════════════════════════════════════════════════════════
# L. METAL-FREE N-CARBON — zero metal cost, intrinsic ORR activity
# ═══════════════════════════════════════════════════════════════════════════════

MFC_N_TYPES = ['pyridinic', 'pyrrolic', 'graphitic', 'oxidized', 'mixed_pyridinic_graphitic']
MFC_N_FRACTIONS = [0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20]
MFC_DEFECT_TYPES = ['none', 'Stone_Wales', 'divacancy', 'edge_zigzag', 'edge_armchair', 'pentagon_heptagon']
MFC_SUBSTRATES = ['graphene', 'CNT', 'graphdiyne', 'porous_carbon', 'carbon_black', 'graphene_nanoribbon']
MFC_DOPANTS = ['B', 'S', 'P', 'F', 'none']  # co-dopants with nitrogen

def generate_mfc_genome() -> tuple:
    """Generate a random Metal-Free N-Carbon genome.
    
    Genome: ('MetalFreeCarbon', n_type, n_fraction, defect, substrate, co_dopant)
    Pure N-doped carbon without any metal center — zero catalyst cost.
    """
    n_type = random.choice(MFC_N_TYPES)
    n_frac = random.choice(MFC_N_FRACTIONS)
    defect = random.choice(MFC_DEFECT_TYPES)
    substrate = random.choice(MFC_SUBSTRATES)
    co_dop = random.choice(MFC_DOPANTS)
    return ('MetalFreeCarbon', n_type, n_frac, defect, substrate, co_dop)



# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

ALL_MATERIAL_CLASSES = [
    'MoltenMetal', 'SolidCatalyst', 'SAC', 'DAC',
    'MOF', 'COF', 'Perovskite', 'MetalHydride', 'MAXPhase', 'HEA',
    'Spinel', 'MXene', 'SAA', 'MetalFreeCarbon',
]

# Weights for random generation — balanced across all classes
CLASS_WEIGHTS = {
    'MoltenMetal': 0.12,       # Primary focus for turquoise H₂
    'SolidCatalyst': 0.10,     # Well-established technology
    'SAC': 0.12,               # 2026 breakthrough class (with axial ligands)
    'DAC': 0.08,               # Cutting-edge dual-atom
    'MOF': 0.06,               # Porous frameworks
    'COF': 0.05,               # Covalent frameworks
    'Perovskite': 0.06,        # ABO₃ oxides
    'MetalHydride': 0.05,      # H₂ storage & decomposition
    'MAXPhase': 0.04,          # Layered ternary ceramics
    'HEA': 0.06,               # High-entropy alloys
    'Spinel': 0.08,            # AB₂O₄ bifunctional — proven ORR/OER
    'MXene': 0.06,             # 2D carbides — high surface area
    'SAA': 0.06,               # Single-atom alloys — ultra-low PGM
    'MetalFreeCarbon': 0.06,   # Zero-cost metal-free ORR
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
    'Spinel': generate_spinel_genome,
    'MXene': generate_mxene_genome,
    'SAA': generate_saa_genome,
    'MetalFreeCarbon': generate_mfc_genome,
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
    HEA_ELEMENTS + SPINEL_A_METALS + SPINEL_B_METALS +
    MXENE_M_ELEMENTS + MXENE_SAC_METALS +
    SAA_TRACE_METALS + SAA_HOST_METALS + ['None']
))
_METAL_IDX = {m: i for i, m in enumerate(_ALL_METALS)}

_ALL_SUPPORTS = sorted(set(SOLID_SUPPORTS + SAC_SUBSTRATES))
_SUPPORT_IDX = {s: i for i, s in enumerate(_ALL_SUPPORTS)}

_ALL_FACETS_LINKERS = sorted(set(
    SOLID_FACETS + MOF_LINKERS + COF_LINKAGES
))
_FACET_LINKER_IDX = {f: i for i, f in enumerate(_ALL_FACETS_LINKERS)}

_ALL_COORDS = sorted(set(
    SAC_COORDINATIONS + DAC_COORDINATIONS + MOF_CAVITIES +
    SAC_AXIAL_LIGANDS + MXENE_TERMINATIONS +
    MFC_N_TYPES + MFC_DEFECT_TYPES + SPINEL_MORPHOLOGIES
))
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
        metal = genome[1]
        coord = genome[2]
        substrate = genome[3]
        axial = genome[4] if len(genome) > 4 else 'none'
        if metal in _METAL_IDX:
            metal1_vec[_METAL_IDX[metal]] = 1.0
        if coord in _COORD_IDX:
            coord_vec[_COORD_IDX[coord]] = 1.0
        if substrate in _SUPPORT_IDX:
            support_vec[_SUPPORT_IDX[substrate]] = 1.0
        # Encode axial ligand as a second coordination feature
        if axial in _COORD_IDX:
            coord_vec[_COORD_IDX[axial]] = 0.5  # half weight to distinguish from in-plane

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

    elif mat_class == 'Spinel':
        _, A, B, dopant, morph, support = genome
        if A in _METAL_IDX:
            metal1_vec[_METAL_IDX[A]] = 1.0
        if B in _METAL_IDX:
            metal2_vec[_METAL_IDX[B]] = 1.0
        if dopant in _DOPANT_IDX:
            dop_vec[_DOPANT_IDX[dopant]] = 1.0
        if morph in _COORD_IDX:
            coord_vec[_COORD_IDX[morph]] = 1.0
        if support in _SUPPORT_IDX:
            support_vec[_SUPPORT_IDX[support]] = 1.0

    elif mat_class == 'MXene':
        _, M, X, n, term, sac_metal = genome
        if M in _METAL_IDX:
            metal1_vec[_METAL_IDX[M]] = 1.0
        if sac_metal in _METAL_IDX:
            metal2_vec[_METAL_IDX[sac_metal]] = 1.0
        if term in _COORD_IDX:
            coord_vec[_COORD_IDX[term]] = 1.0
        cont[0] = 1.0 if X == 'N' else (0.5 if X == 'CN' else 0.0)
        cont[1] = n / 3.0

    elif mat_class == 'SAA':
        _, trace, host, facet, loading = genome
        if trace in _METAL_IDX:
            metal1_vec[_METAL_IDX[trace]] = 1.0
        if host in _METAL_IDX:
            metal2_vec[_METAL_IDX[host]] = 1.0
        _saa_facet_map = {'111': 0.0, '100': 0.33, '110': 0.67, '211': 1.0}
        cont[0] = _saa_facet_map.get(facet, 0.0)
        cont[1] = loading / 10000.0  # normalized loading

    elif mat_class == 'MetalFreeCarbon':
        _, n_type, n_frac, defect, substrate, co_dop = genome
        # No metal encoding — this is the point
        if n_type in _COORD_IDX:
            coord_vec[_COORD_IDX[n_type]] = 1.0
        if defect in _COORD_IDX:
            coord_vec[_COORD_IDX[defect]] = 0.5
        if substrate in _SUPPORT_IDX:
            support_vec[_SUPPORT_IDX[substrate]] = 1.0
        cont[0] = n_frac / 0.20  # normalized N fraction
        _codop_map = {'B': 0.2, 'S': 0.4, 'P': 0.6, 'F': 0.8, 'none': 0.0}
        cont[1] = _codop_map.get(co_dop, 0.0)

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
        gene_idx = random.randint(1, 4)
        if gene_idx == 1:
            genes[1] = random.choice(SAC_METALS)
        elif gene_idx == 2:
            genes[2] = random.choice(SAC_COORDINATIONS)
        elif gene_idx == 3:
            genes[3] = random.choice(SAC_SUBSTRATES)
        elif gene_idx == 4:
            if len(genes) > 4:
                genes[4] = random.choice(SAC_AXIAL_LIGANDS)
            else:
                genes.append(random.choice(SAC_AXIAL_LIGANDS))

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

    elif mat_class == 'Spinel':
        gene_idx = random.randint(1, 5)
        if gene_idx == 1:
            genes[1] = random.choice(SPINEL_A_METALS)
        elif gene_idx == 2:
            genes[2] = random.choice(SPINEL_B_METALS)
        elif gene_idx == 3:
            genes[3] = random.choice(SPINEL_DOPANTS)
        elif gene_idx == 4:
            genes[4] = random.choice(SPINEL_MORPHOLOGIES)
        elif gene_idx == 5:
            genes[5] = random.choice(SPINEL_SUPPORT_CARBONS)

    elif mat_class == 'MXene':
        gene_idx = random.randint(1, 5)
        if gene_idx == 1:
            genes[1] = random.choice(MXENE_M_ELEMENTS)
        elif gene_idx == 2:
            genes[2] = random.choice(MXENE_X_ELEMENTS)
        elif gene_idx == 3:
            genes[3] = random.choice(MXENE_N_VALUES)
        elif gene_idx == 4:
            genes[4] = random.choice(MXENE_TERMINATIONS)
        elif gene_idx == 5:
            genes[5] = random.choice(MXENE_SAC_METALS)

    elif mat_class == 'SAA':
        gene_idx = random.randint(1, 4)
        if gene_idx == 1:
            genes[1] = random.choice(SAA_TRACE_METALS)
        elif gene_idx == 2:
            genes[2] = random.choice(SAA_HOST_METALS)
        elif gene_idx == 3:
            genes[3] = random.choice(SAA_FACETS)
        elif gene_idx == 4:
            genes[4] = random.choice(SAA_LOADINGS_PPM)

    elif mat_class == 'MetalFreeCarbon':
        gene_idx = random.randint(1, 5)
        if gene_idx == 1:
            genes[1] = random.choice(MFC_N_TYPES)
        elif gene_idx == 2:
            genes[2] = random.choice(MFC_N_FRACTIONS)
        elif gene_idx == 3:
            genes[3] = random.choice(MFC_DEFECT_TYPES)
        elif gene_idx == 4:
            genes[4] = random.choice(MFC_SUBSTRATES)
        elif gene_idx == 5:
            genes[5] = random.choice(MFC_DOPANTS)

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
        'SAC': (
            len(SAC_METALS) * len(SAC_COORDINATIONS) *
            len(SAC_SUBSTRATES) * len(SAC_AXIAL_LIGANDS)
        ),
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
        'Spinel': (
            len(SPINEL_A_METALS) * len(SPINEL_B_METALS) *
            len(SPINEL_DOPANTS) * len(SPINEL_MORPHOLOGIES) *
            len(SPINEL_SUPPORT_CARBONS)
        ),
        'MXene': (
            len(MXENE_M_ELEMENTS) * len(MXENE_X_ELEMENTS) *
            len(MXENE_N_VALUES) * len(MXENE_TERMINATIONS) *
            len(MXENE_SAC_METALS)
        ),
        'SAA': (
            len(SAA_TRACE_METALS) * len(SAA_HOST_METALS) *
            len(SAA_FACETS) * len(SAA_LOADINGS_PPM)
        ),
        'MetalFreeCarbon': (
            len(MFC_N_TYPES) * len(MFC_N_FRACTIONS) *
            len(MFC_DEFECT_TYPES) * len(MFC_SUBSTRATES) *
            len(MFC_DOPANTS)
        ),
    }
    sizes['TOTAL'] = sum(sizes.values())
    return sizes


def generate_hierarchical_htvs_pool(pool_size: int, scorer=None,
                                    campaign_round: int = 0) -> List[tuple]:
    """
    Generate a deterministic, non-random pool of catalyst candidates using hierarchical screening.
    
    1. Exhaustively builds grids for small classes and core configurations for large classes.
    2. Screens the cores with the scorer callback (if provided).
    3. Expands the best-performing cores with local modifications (dopants, strain, defects).
    4. Selects and returns the best candidates.
    """
    # Rotate every strided dimension between rounds.  Unlike a fixed [::2]
    # slice, successive reinjection rounds cover the complementary half.
    # Each categorical dimension uses a different bit of campaign_round.  Over
    # 2**N rounds this visits every even/odd Cartesian shard, including mixed
    # combinations (not merely the all-even and all-odd diagonals).
    stride2 = lambda values, dimension: values[(campaign_round >> dimension) & 1::2]

    # A. Generate core configurations deterministically
    cores = []

    # 1. MoltenMetal (stride to ~5000)
    for host in MOLTEN_HOSTS:
        for promoter in stride2(MOLTEN_PROMOTERS, 0):
            for at_pct in stride2(MOLTEN_PROMOTER_AT_PCT, 1):
                for temp in stride2(MOLTEN_TEMPERATURES_K, 2):
                    cores.append(('MoltenMetal', host, promoter, at_pct, temp))

    # 2. SAC (100% exhaustive: 6,318)
    for metal in SAC_METALS:
        for coord in SAC_COORDINATIONS:
            for substrate in SAC_SUBSTRATES:
                for axial in SAC_AXIAL_LIGANDS:
                    cores.append(('SAC', metal, coord, substrate, axial))

    # 3. DAC (strided to ~10,000)
    for m1 in stride2(DAC_METALS_1, 0):
        for m2 in stride2(DAC_METALS_2, 1):
            for coord in stride2(DAC_COORDINATIONS, 2):
                for substrate in SAC_SUBSTRATES:
                    cores.append(('DAC', m1, m2, coord, substrate))

    # 4. MOF (strided to ~5,000)
    for metal in stride2(MOF_METAL_NODES, 0):
        for linker in stride2(MOF_LINKERS, 1):
            for cavity in stride2(MOF_CAVITIES, 2):
                for pore in stride2(MOF_PORE_SIZES, 3):
                    cores.append(('MOF', metal, linker, cavity, pore))

    # 5. COF (strided to ~5,000)
    for metal in stride2(MOF_METAL_NODES + ['None'], 0):
        for linkage in stride2(COF_LINKAGES, 1):
            for cavity in stride2(MOF_CAVITIES, 2):
                for pore in stride2(MOF_PORE_SIZES, 3):
                    cores.append(('COF', metal, linkage, cavity, pore))

    # 6. Perovskite (strided to ~3,000)
    for A in stride2(PEROVSKITE_A_SITE, 0):
        for B in stride2(PEROVSKITE_B_SITE, 1):
            for dopant in stride2(PEROVSKITE_B_SITE + ['None'], 2):
                frac = 0.0 if dopant == 'None' else PEROVSKITE_DOPANT_FRAC[(campaign_round % (len(PEROVSKITE_DOPANT_FRAC)-1))+1]
                defect = PEROVSKITE_DEFECTS[campaign_round % len(PEROVSKITE_DEFECTS)]
                cores.append(('Perovskite', A, B, dopant, frac, defect))

    # 7. MetalHydride (strided to ~5,000)
    for metal in stride2(HYDRIDE_METALS, 0):
        for h_type in stride2(HYDRIDE_TYPES, 1):
            for second in stride2(HYDRIDE_SECOND_METAL, 2):
                additive = HYDRIDE_ADDITIVES[campaign_round % len(HYDRIDE_ADDITIVES)]
                temp = [300, 350, 400, 450, 500, 550, 600, 700, 800][campaign_round % 9]
                cores.append(('MetalHydride', metal, h_type, second, additive, temp))

    # 8. MAXPhase (100% exhaustive: 49,140)
    for M in MAX_M_ELEMENTS:
        for A in MAX_A_ELEMENTS:
            for X in MAX_X_ELEMENTS:
                for N in MAX_N_VALUES:
                    cores.append(('MAXPhase', M, A, X, N, 'None', 'basal_0001'))

    # 9. HEA (strided to ~10,000)
    from itertools import combinations
    n_elements = len(HEA_ELEMENTS)
    step = max(1, n_elements // 10)
    subset_elements = [HEA_ELEMENTS[i] for i in range(0, n_elements, step)]
    for combi in combinations(subset_elements, 4):
        for struct in HEA_STRUCTURES[:2]:
            for facet in ['111']:
                cores.append(('HEA', tuple(sorted(combi)), struct, facet, 1000))

    # 10. Spinel (100% exhaustive: 14,400)
    for A in SPINEL_A_METALS:
        for B in SPINEL_B_METALS:
            for dopant in SPINEL_DOPANTS:
                for morph in SPINEL_MORPHOLOGIES:
                    for support in SPINEL_SUPPORT_CARBONS:
                        cores.append(('Spinel', A, B, dopant, morph, support))

    # 11. MXene (100% exhaustive: 2,400)
    for M in MXENE_M_ELEMENTS:
        for X in MXENE_X_ELEMENTS:
            for N in MXENE_N_VALUES:
                for term in MXENE_TERMINATIONS:
                    for sac in MXENE_SAC_METALS:
                        cores.append(('MXene', M, X, N, term, sac))

    # 12. SAA (100% exhaustive: 1,500)
    for trace in SAA_TRACE_METALS:
        for host in SAA_HOST_METALS:
            for facet in SAA_FACETS:
                for loading in SAA_LOADINGS_PPM:
                    cores.append(('SAA', trace, host, facet, loading))

    # 13. MetalFreeCarbon (100% exhaustive: 1,575)
    for n_type in MFC_N_TYPES:
        for n_frac in MFC_N_FRACTIONS:
            for defect in MFC_DEFECT_TYPES:
                for substrate in MFC_SUBSTRATES:
                    for dopant in MFC_DOPANTS:
                        cores.append(('MetalFreeCarbon', n_type, n_frac, defect, substrate, dopant))

    # 14. SolidCatalyst Cores (all 25,800 combinations)
    for metal in SOLID_ACTIVE_METALS:
        for support in SOLID_SUPPORTS:
            for facet in SOLID_FACETS:
                cores.append(('SolidCatalyst', metal, support, facet, 0.0, (), 1, 0))

    if scorer is None:
        step = max(1, len(cores) // pool_size)
        return cores[::step][:pool_size]

    # B. Screen cores to select the most promising ones
    acquisition_scores = scorer(cores)
    best_indices = np.argsort(acquisition_scores)

    # C. Perform neighborhood expansion for the top cores
    expanded_pool = []
    solid_expanded_count = 0
    common_dopants = ['B', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Si', 'Al']

    for idx in best_indices:
        g = cores[idx]
        mat_class = g[0]
        if mat_class == 'SolidCatalyst':
            if solid_expanded_count >= 500:
                expanded_pool.append(g)
                continue
            _, metal, support, facet, _, _, _, _ = g
            for strain in [-0.05, 0.05]:
                for vac in [0, 1]:
                    expanded_pool.append(('SolidCatalyst', metal, support, facet, strain, (), 1, vac))
                    for d in common_dopants:
                        expanded_pool.append(('SolidCatalyst', metal, support, facet, strain, (d,), 1, vac))
            solid_expanded_count += 1
        elif mat_class == 'Perovskite':
            _, A, B, _, _, _ = g
            for dopant in PEROVSKITE_B_SITE[:3]:
                for frac in [0.05, 0.1]:
                    for defect in PEROVSKITE_DEFECTS:
                        expanded_pool.append(('Perovskite', A, B, dopant, frac, defect))
        elif mat_class == 'MetalHydride':
            _, metal, h_type, second, _, _ = g
            for additive in HYDRIDE_ADDITIVES[1:4]:
                for temp in [400, 600]:
                    expanded_pool.append(('MetalHydride', metal, h_type, second, additive, temp))
        else:
            expanded_pool.append(g)

    full_pool = list(set(cores + expanded_pool))

    # D. Final screen of the combined pool
    scores = scorer(full_pool)
    best_full_indices = np.argsort(scores)[:pool_size]
    return [full_pool[i] for i in best_full_indices]


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
