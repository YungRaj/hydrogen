#!/usr/bin/env python3
"""
Shared utilities, physical constants, logging, and helper functions
for the turquoise hydrogen → fuel cell simulation pipeline.
"""

import os
import sys
import time
import json
import hashlib
import logging
import numpy as np
try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

# ─── Project Paths ──────────────────────────────────────────────────────────────

BASE_DIR = Path(__file__).resolve().parent.parent
PIPELINE_DIR = BASE_DIR / "pipeline"
MECHANISMS_DIR = BASE_DIR / "mechanisms"
RESULTS_DIR = BASE_DIR / "results"
SCREENING_DIR = RESULTS_DIR / "screening"
REACTOR_DIR = RESULTS_DIR / "reactor"
DFT_DIR = RESULTS_DIR / "dft"
VQE_DIR = RESULTS_DIR / "vqe"
FUEL_CELL_DIR = RESULTS_DIR / "fuel_cell"
REPORTS_DIR = RESULTS_DIR / "reports"
QE_PSEUDO_DIR = BASE_DIR / "quantum_espresso" / "pseudo"

# ─── Conda Environments ─────────────────────────────────────────────────────────

ENV_MACE = "deepmd-env"         # MACE-MP-0, PyTorch, ASE
ENV_BATTERY = "battery-env"     # ASE, Pymatgen, Pandas
ENV_CANTERA = "cp2k-env"        # Cantera 3.2
ENV_QE = "qe-env"               # Quantum ESPRESSO
ENV_QUANTUM = "quantum-env"     # CUDA-Q

# ─── Physical Constants ─────────────────────────────────────────────────────────

# Fundamental
k_B_eV = 8.617333262145e-5      # Boltzmann constant (eV/K)
k_B_J = 1.380649e-23            # Boltzmann constant (J/K)
h_eV = 4.135667696e-15          # Planck constant (eV·s)
h_J = 6.62607015e-34            # Planck constant (J·s)
R_gas = 8.314462618             # Universal gas constant (J/(mol·K))
F_const = 96485.33212           # Faraday constant (C/mol)
eV_to_J = 1.602176634e-19       # eV to Joules
Ry_to_eV = 13.605693122994      # Rydberg to eV
Ha_to_eV = 27.211386245988      # Hartree to eV
amu_to_kg = 1.66053906660e-27   # Atomic mass unit to kg

# Thermodynamic reference values
E_H2O_G = -2.4583               # Free energy correction for H₂O(l) at 298 K (eV)
ZPE_H2 = 0.27                   # Zero-point energy of H₂ (eV)
TS_H2 = 0.40                    # T·S for H₂ at 298 K, 1 bar (eV)
G_H2_correction = ZPE_H2 - TS_H2  # Free energy correction for H₂

# Standard electrode potentials
E_SHE = 0.0                     # Standard hydrogen electrode (V)
E_ORR_eq = 1.229                # Equilibrium potential for ORR at 298 K (V)

# ─── Atomic Data ─────────────────────────────────────────────────────────────────

ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Na': 22.990,
    'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.06,
    'Cl': 35.45, 'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867,
    'V': 50.942, 'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933,
    'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.63,
    'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Y': 88.906, 'Zr': 91.224,
    'Nb': 92.906, 'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.906, 'Pd': 106.42,
    'Ag': 107.868, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.710,
    'Sb': 121.760, 'Te': 127.60, 'La': 138.905, 'Ce': 140.116,
    'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
    'Pt': 195.084, 'Au': 196.967, 'Pb': 207.2, 'Bi': 208.980,
}

# USGS crustal abundance (ppm by mass) — used for cost/availability scoring
CRUSTAL_ABUNDANCE_PPM = {
    # Transition metals (core design space)
    'Ni': 84.0, 'Pt': 0.005, 'Cu': 60.0, 'Pd': 0.015, 'Fe': 56300.0,
    'Co': 25.0, 'Ru': 0.001, 'Rh': 0.001, 'Ir': 0.001, 'Au': 0.004,
    'Ag': 0.075, 'W': 1.25, 'Mo': 1.2, 'Ti': 5650.0, 'V': 120.0,
    'Cr': 102.0, 'Mn': 950.0, 'Zn': 70.0, 'Zr': 165.0, 'Nb': 20.0,
    'Ta': 2.0, 'Re': 0.0007, 'Hf': 3.0, 'Sc': 22.0, 'Y': 33.0,
    # Post-transition metals & metalloids
    'Al': 82300.0, 'Ga': 19.0, 'In': 0.25, 'Sn': 2.3, 'Bi': 0.008,
    'Sb': 0.2, 'Pb': 14.0, 'Ge': 1.5, 'Te': 0.001, 'Se': 0.05,
    'Cd': 0.15, 'Os': 0.0015,
    # Alkali & alkaline earth (perovskite A-sites, spinel)
    'Li': 20.0, 'Na': 23600.0, 'K': 20900.0, 'Rb': 90.0, 'Cs': 3.0,
    'Mg': 23300.0, 'Ca': 41500.0, 'Sr': 370.0, 'Ba': 425.0,
    # Lanthanides (perovskite A-sites)
    'La': 39.0, 'Ce': 66.5, 'Pr': 9.2, 'Nd': 41.5, 'Sm': 7.05,
    'Eu': 2.0, 'Gd': 6.2, 'Dy': 5.2, 'Er': 3.5, 'Yb': 3.2,
    # Non-metals (MXene terminations, N-doped carbon, linkers)
    'B': 10.0, 'C': 200.0, 'N': 19.0, 'O': 461000.0, 'F': 585.0,
    'P': 1050.0, 'S': 350.0, 'Si': 282000.0,
    'Cl': 145.0, 'Br': 2.4, 'I': 0.45,
}

# Approximate bulk commodity price ($/kg, 2025 est.)
METAL_PRICE_USD_KG = {
    # PGMs and precious metals
    'Pt': 31000.0, 'Pd': 40000.0, 'Rh': 145000.0, 'Ir': 150000.0,
    'Ru': 14000.0, 'Au': 75000.0, 'Ag': 850.0, 'Re': 4500.0,
    # Base transition metals
    'Ni': 16.0, 'Cu': 9.0, 'Fe': 0.10, 'Co': 33.0, 'Ti': 11.0,
    'V': 30.0, 'Cr': 10.0, 'Mn': 2.0, 'Zn': 2.7, 'Zr': 35.0,
    'Nb': 73.0, 'Mo': 45.0, 'W': 35.0, 'Hf': 900.0, 'Ta': 300.0,
    'Sc': 3500.0,
    # Post-transition metals
    'Al': 2.5, 'Ga': 300.0, 'In': 300.0, 'Sn': 25.0, 'Bi': 8.0,
    'Sb': 12.0, 'Te': 80.0, 'Pb': 2.1, 'Ge': 1200.0, 'Se': 45.0,
    # Alkali & alkaline earth
    'Li': 70.0, 'Na': 0.20, 'K': 0.80, 'Rb': 12000.0, 'Cs': 60000.0,
    'Mg': 2.5, 'Ca': 2.0, 'Sr': 6.5, 'Ba': 0.30,
    # Lanthanides
    'La': 5.0, 'Ce': 5.0, 'Y': 35.0, 'Pr': 100.0, 'Nd': 75.0,
    'Sm': 15.0, 'Eu': 500.0, 'Gd': 55.0, 'Dy': 350.0, 'Er': 70.0,
    'Yb': 60.0,
    # Non-metals (commodity chemicals, very cheap per kg)
    'B': 5.0, 'C': 0.10, 'N': 0.50, 'O': 0.20, 'F': 2.0,
    'P': 1.0, 'S': 0.10, 'Si': 2.0, 'Cl': 0.30, 'Br': 4.0, 'I': 35.0,
}

# Melting points (K) — critical for molten metal catalyst feasibility
MELTING_POINT_K = {
    'Sn': 505.1, 'Bi': 544.6, 'In': 429.7, 'Ga': 302.9, 'Pb': 600.6,
    'Sb': 903.8, 'Te': 722.7, 'Ni': 1728.0, 'Cu': 1358.0, 'Fe': 1811.0,
    'Co': 1768.0, 'Mn': 1519.0, 'Pd': 1828.0, 'Pt': 2041.0, 'Al': 933.5,
    'Zn': 692.7, 'Ag': 1234.9, 'Au': 1337.3, 'Mo': 2896.0, 'W': 3695.0,
}

# ─── Logging ─────────────────────────────────────────────────────────────────────

def setup_logger(name: str, log_file: Optional[str] = None, level=logging.INFO) -> logging.Logger:
    """Create a configured logger with both console and file output."""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers = []  # Clear any existing handlers

    formatter = logging.Formatter(
        '[%(asctime)s] %(name)s | %(levelname)s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File handler
    if log_file:
        log_path = RESULTS_DIR / log_file
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(str(log_path), mode='a')
        fh.setLevel(level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def print_banner(title: str, width: int = 80):
    """Print a formatted section banner."""
    print("=" * width)
    print(f" {title.center(width - 2)} ")
    print("=" * width)


# ─── Genome Hashing & Serialization ─────────────────────────────────────────────

def genome_hash(genome: tuple) -> str:
    """Create a deterministic hash for a catalyst genome for deduplication."""
    serialized = json.dumps(genome, sort_keys=True, default=str)
    return hashlib.md5(serialized.encode()).hexdigest()[:12]


def genome_to_dict(genome: tuple, material_class: str) -> Dict[str, Any]:
    """Convert a genome tuple to a labeled dictionary based on material class."""
    if material_class == 'MoltenMetal':
        return {
            'material_class': material_class,
            'host_metal': genome[0],
            'promoter': genome[1],
            'promoter_at_pct': genome[2],
            'temperature_K': genome[3],
        }
    elif material_class == 'SolidCatalyst':
        return {
            'material_class': material_class,
            'active_metal': genome[0],
            'support': genome[1],
            'facet': genome[2],
            'strain': genome[3],
            'dopants': genome[4],
            'num_substitutions': genome[5],
            'num_vacancies': genome[6],
        }
    elif material_class in ('SAC', 'DAC'):
        return {
            'material_class': material_class,
            'metal_1': genome[0],
            'metal_2': genome[1] if material_class == 'DAC' else None,
            'coordination': genome[2] if material_class == 'DAC' else genome[1],
            'substrate': genome[3] if material_class == 'DAC' else genome[2],
        }
    elif material_class in ('MOF', 'COF'):
        return {
            'material_class': material_class,
            'metal_node': genome[0],
            'linker': genome[1],
            'cavity': genome[2],
            'pore_size_A': genome[3],
        }
    else:
        return {'material_class': material_class, 'genome': genome}


# ─── Cost & Feasibility Scoring ──────────────────────────────────────────────────

def abundance_cost_penalty(elements: List[str]) -> float:
    """
    Compute a log-scaled cost penalty based on crustal abundance.
    More abundant elements get lower (better) penalties.
    Returns a value in [−2, 0] where 0 = abundant, −2 = very rare.

    Uses geometric mean to prevent a single abundant element from
    masking rare/expensive components (e.g. Fe+Ir should still penalize Ir).
    """
    if not elements:
        return 0.0
    abundances = [max(CRUSTAL_ABUNDANCE_PPM.get(e, 1.0), 1e-4) for e in elements]
    # Geometric mean: exp(mean(log(x))) — sensitive to any rare component
    geo_mean = np.exp(np.mean(np.log(abundances)))
    if geo_mean < 1.0:
        return -2.0
    elif geo_mean < 10.0:
        return 2.0 * (np.log10(geo_mean) - 1.0)  # range: [-2, 0]
    else:
        return 0.0



# ─── Safety & Feasibility Filters ─────────────────────────────────────────────

# Elements that should never appear in a recommended catalyst
TOXIC_ELEMENTS = {
    'Os',  # OsO₄ is extremely volatile and toxic
    'Tl',  # Lethal poison (thallium)
    'Po',  # Radioactive
    'Tc',  # Radioactive (technetium)
    'Pm',  # Radioactive (promethium)
    'Cd',  # Carcinogenic heavy metal
    'Hg',  # Neurotoxic
    'Be',  # Berylliosis risk
    'As',  # Carcinogenic
}

# Material classes valid for each application
# All classes are included to ensure exhaustive exploration.
# Physical viability is already handled by OOD confidence and
# the multi-fidelity validation pipeline (Tier 3 DFT for champions).
VALID_CLASSES_PYROLYSIS = {
    'SolidCatalyst', 'MoltenMetal', 'HEA', 'MAXPhase',
    'Perovskite', 'MetalHydride',
    'Spinel',          # Earth-abundant oxides, stable at pyrolysis temps
    'MXene',           # Ti₃C₂ stable to ~800°C, catalytically active
    'SAA',             # Dilute alloys, same stability as SolidCatalyst
    'SAC',             # Carbon support marginal at >700°C — flag in report
    'DAC',             # Same caveat as SAC
    'MetalFreeCarbon', # N-carbon stable to ~600°C — flag in report
    # MOF/COF: decompose at >300°C but included for completeness — OOD handles
    'MOF', 'COF',
}

VALID_CLASSES_FUEL_CELL = {
    'SolidCatalyst', 'SAC', 'DAC', 'HEA', 'MAXPhase',
    'Perovskite',
    'Spinel',          # NiFe₂O₄ is a proven bifunctional ORR/OER catalyst
    'MXene',           # Promising ORR activity in alkaline media
    'SAA',             # Pt₁Cu etc. are frontier FC catalysts
    'MetalFreeCarbon', # N-doped graphene is a benchmark PGM-free ORR catalyst
    # MOF/COF: unstable in acid PEMFC but viable in AEM fuel cells
    'MOF', 'COF',
    # MetalHydride/MoltenMetal: non-standard for FC, but included to avoid pruning
    'MetalHydride', 'MoltenMetal',
}


def check_element_safety(elements: List[str]) -> Tuple[bool, str]:
    """
    Check if a catalyst's elements are safe for recommendation.
    Returns (is_safe, reason).
    """
    toxic_found = [e for e in elements if e in TOXIC_ELEMENTS]
    if toxic_found:
        return False, f"Contains toxic/radioactive elements: {', '.join(toxic_found)}"
    return True, "OK"


def is_valid_for_application(material_class: str, application: str = 'pyrolysis') -> bool:
    """Check if a material class is physically viable for the target application."""
    if application == 'pyrolysis':
        return material_class in VALID_CLASSES_PYROLYSIS
    elif application in ('fuel_cell', 'orr'):
        return material_class in VALID_CLASSES_FUEL_CELL
    return True  # unknown application — allow all


def material_cost_usd_per_kg(elements: List[str], fractions: Optional[List[float]] = None) -> float:
    """Estimate raw material cost in $/kg for a multi-component catalyst."""
    if fractions is None:
        fractions = [1.0 / len(elements)] * len(elements)
    cost = sum(
        METAL_PRICE_USD_KG.get(e, 50.0) * f
        for e, f in zip(elements, fractions)
    )
    return cost


def is_molten_at_temperature(metal: str, temperature_K: float) -> bool:
    """Check if a metal is molten at the given temperature."""
    mp = MELTING_POINT_K.get(metal, 3000.0)
    return temperature_K > mp


# ─── Arrhenius & Kinetics ────────────────────────────────────────────────────────

def arrhenius_rate(A: float, E_act_eV: float, T_K: float) -> float:
    """Compute Arrhenius rate constant k = A * exp(-Eₐ / kT)."""
    if T_K <= 0:
        return 0.0
    return A * np.exp(-E_act_eV / (k_B_eV * T_K))


def tst_prefactor(T_K: float, delta_S_eV_K: float = 0.0) -> float:
    """
    Transition State Theory pre-exponential factor:
    A = (kT/h) * exp(ΔS‡/k)
    """
    A_tst = (k_B_eV * T_K) / h_eV
    if delta_S_eV_K != 0.0:
        A_tst *= np.exp(delta_S_eV_K / k_B_eV)
    return A_tst


def bep_activation_energy(delta_E_rxn: float, alpha: float = 0.87,
                          beta: float = 0.75, material_class: str = None) -> float:
    """
    Brønsted-Evans-Polanyi (BEP) correlation for activation energy:
    Eₐ = alpha + beta * ΔE_rxn  (for exothermic reactions, ΔE < 0)

    Class-specific BEP parameters from literature:
      - Transition metals: Nørskov et al., J. Catal. 2002
      - Metal oxides: Vojvodic et al., Chem. Rev. 2014
      - Zeolites: Bligaard et al., J. Catal. 2004
      - Molten metals: Upham et al., Science 2017
      - SAC/DAC: Li et al., Nat. Catal. 2019 (approximate)
      - Others: Use metal defaults with uncertainty flag

    Returns Eₐ in eV, clamped to [0.01, 5.0].
    """
    # Class-specific BEP parameters (alpha = intercept, beta = slope)
    BEP_PARAMS = {
        'SolidCatalyst': (0.87, 0.75),   # Transition metal surfaces
        'HEA':           (0.90, 0.78),   # Multi-component alloys (slightly higher barrier)
        'MoltenMetal':   (0.50, 0.60),   # Liquid metal catalysis (Upham et al.)
        'SAC':           (0.95, 0.70),   # Single-atom (stronger binding, modified scaling)
        'DAC':           (0.92, 0.72),   # Dual-atom (intermediate)
        'Perovskite':    (1.50, 0.55),   # Oxide surfaces (higher intercept)
        'MAXPhase':      (1.00, 0.70),   # Carbide/nitride surfaces
        'MetalHydride':  (0.80, 0.65),   # Hydride surfaces
        'MOF':           (2.00, 0.40),   # Framework catalysts (very different scaling)
        'COF':           (2.00, 0.40),   # Same as MOF (limited data)
    }

    if material_class and material_class in BEP_PARAMS:
        alpha, beta = BEP_PARAMS[material_class]

    E_act = alpha + beta * delta_E_rxn
    return float(np.clip(E_act, 0.01, 5.0))


# ─── Electrochemistry ───────────────────────────────────────────────────────────

def orr_overpotential(dG_OH: float, dG_O: float, dG_OOH: float) -> Tuple[float, str]:
    """
    Compute ORR overpotential using the computational hydrogen electrode (CHE).
    
    The 4-electron ORR pathway:
      O₂ + * + H⁺ + e⁻ → OOH*      ΔG₁ = dG_OOH
      OOH* + H⁺ + e⁻ → O* + H₂O    ΔG₂ = dG_O - dG_OOH + 3.33 (water correction)
      O* + H⁺ + e⁻ → OH*            ΔG₃ = dG_OH - dG_O
      OH* + H⁺ + e⁻ → H₂O + *      ΔG₄ = -dG_OH

    η = max(ΔGᵢ)/e − 1.23 V
    
    Returns (overpotential_V, rate_determining_step).
    """
    dG1 = dG_OOH - 4.92  # relative to O₂ + 2H₂O reference
    dG2 = dG_O - dG_OOH
    dG3 = dG_OH - dG_O
    dG4 = -dG_OH

    # At U = 0, all steps must be downhill. At U = U_eq, the potential-determining
    # step is the one with the largest (most positive) ΔG.
    # ΔG at applied potential U: ΔGᵢ(U) = ΔGᵢ + eU
    # Limiting potential: U_L = −max(ΔGᵢ) / e
    # Overpotential: η = 1.23 − U_L

    steps = {'step_1_OOH': dG1, 'step_2_O': dG2, 'step_3_OH': dG3, 'step_4_H2O': dG4}
    rds_name = max(steps, key=steps.get)
    rds_dG = steps[rds_name]

    U_limiting = -rds_dG
    eta = E_ORR_eq - U_limiting
    return max(eta, 0.0), rds_name


def butler_volmer_current(j0: float, eta: float, alpha_a: float = 0.5,
                          alpha_c: float = 0.5, T_K: float = 353.0) -> float:
    """
    Butler-Volmer equation for electrode kinetics.
    j = j₀ * [exp(αₐFη/RT) − exp(−αcFη/RT)]
    
    Args:
        j0: Exchange current density (A/cm²)
        eta: Overpotential (V), positive for anodic
        alpha_a: Anodic transfer coefficient
        alpha_c: Cathodic transfer coefficient  
        T_K: Temperature (K)
    Returns:
        Current density (A/cm²)
    """
    f = F_const / (R_gas * T_K)
    return j0 * (np.exp(alpha_a * f * eta) - np.exp(-alpha_c * f * eta))


# ─── Data I/O ────────────────────────────────────────────────────────────────────

def save_screening_db(df, filename: str, subdir: str = "screening"):
    """Save a screening database as CSV."""
    path = RESULTS_DIR / subdir / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def load_screening_db(filename: str, subdir: str = "screening"):
    """Load a screening database if it exists."""
    path = RESULTS_DIR / subdir / filename
    if path.exists():
        return pd.read_csv(path)
    return None


def save_json(data: dict, filename: str, subdir: str = "reports"):
    """Save a dictionary as JSON."""
    path = RESULTS_DIR / subdir / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, default=str)
    return path


def load_json(filename: str, subdir: str = "reports") -> Optional[dict]:
    """Load a JSON file if it exists."""
    path = RESULTS_DIR / subdir / filename
    if path.exists():
        with open(path, 'r') as f:
            return json.load(f)
    return None


# ─── Subprocess Helpers ──────────────────────────────────────────────────────────

def conda_run_cmd(env_name: str, python_cmd: str, cwd: Optional[str] = None) -> str:
    """Build a conda run command string."""
    cmd = f"conda run -n {env_name} python {python_cmd}"
    return cmd


def run_in_env(env_name: str, script_path: str, args: str = "",
               cwd: Optional[str] = None, check: bool = True):
    """Execute a Python script in a specific conda environment."""
    import subprocess
    cmd = f"conda run -n {env_name} python {script_path} {args}"
    work_dir = cwd or str(BASE_DIR)
    result = subprocess.run(cmd, shell=True, cwd=work_dir, 
                           capture_output=True, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed in {env_name}:\n"
            f"CMD: {cmd}\n"
            f"STDERR: {result.stderr[-2000:]}"
        )
    return result
