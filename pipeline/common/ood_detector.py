#!/usr/bin/env python3
# Shared out-of-distribution confidence policies.
"""
Out-of-Distribution (OOD) Detection for GNN Screening.

The eSen-SM model was trained on ~260M DFT calculations from OC20/OC22/OC25,
covering solid-gas and solid-liquid interfaces for ~88 elements. However,
our 14-class design space includes material types NOT in that training data
(MOFs, metal-free carbon, metal hydrides, etc.).

This module provides three layers of OOD detection:

  Layer 1: CLASS_CONFIDENCE — static prior based on training set coverage.
           Free, applied to every genome.

  Layer 2: element_coverage_score() — checks if the elements in a genome
           were well-represented in OC20/OC25 surface calculations.
           Cheap, applied to every genome.

  Layer 3: dual_model_disagreement() — runs both eSen-SM and MACE-MP-0
           and measures energy disagreement. Expensive, applied only
           to top-k candidates during periodic GA validation.

The combined confidence score is used to discount screening results:
  - confidence > 0.7  → trust the prediction
  - confidence 0.4-0.7 → flag for Tier 3 DFT validation
  - confidence < 0.4  → heavily penalize in NSGA-II ranking
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger('ood_detector')


# ═══════════════════════════════════════════════════════════════════════════════
# LAYER 1: CLASS-LEVEL CONFIDENCE (cost = 0)
# ═══════════════════════════════════════════════════════════════════════════════
# Based on what fraction of each class's bonding environments appear
# in OC20 (gas-surface) + OC22 (oxide surfaces) + OC25 (electrochemical).
#
# Scored 0.0 to 1.0:
#   1.0 = core training data (metal FCC/BCC/HCP slabs)
#   0.0 = completely absent from training

CLASS_CONFIDENCE = {
    'SolidCatalyst':   0.95,  # Metal slabs are the core of OC20
    'HEA':             0.85,  # Multi-component alloys well-covered
    'SAA':             0.85,  # Dilute substitutional alloys in OC20
    'Spinel':          0.60,  # OC22 has metal oxides but spinels are sparse
    'Perovskite':      0.55,  # Some ABO₃ in OC22, not comprehensive
    'MXene':           0.50,  # Ti/Mo carbides partially in OC20
    'MAXPhase':        0.55,  # Carbide surfaces partially in OC20, similar to MXene
    'MoltenMetal':     0.45,  # Modeled as solid slab — liquid structure missing
    'SAC':             0.40,  # Fe-N₄/graphene partially in OC20, but
                              # axial ligand effects are NOT in training
    'DAC':             0.35,  # Dual-atom sites barely in training data
    'MOF':             0.20,  # Metal-organic frameworks NOT in OC20
    'COF':             0.15,  # Covalent organic frameworks NOT in OC20
    'MetalHydride':    0.20,  # Interstitial H in bulk NOT in OC20
    'MetalFreeCarbon': 0.15,  # No pure N-carbon surfaces in OC20
}


# ═══════════════════════════════════════════════════════════════════════════════
# LAYER 2: ELEMENT COVERAGE (cost = O(1) per genome)
# ═══════════════════════════════════════════════════════════════════════════════
# OC20/OC22 training data element frequency. Elements with >10k slab
# calculations in OC20 get score 1.0, elements absent get 0.0.
# From OC20 dataset statistics (Chanussot et al., ACS Catal., 2021).

# Tier 1: Extensively sampled in OC20 (>50k slabs)
_TIER1_ELEMENTS = {
    'Pt', 'Pd', 'Ni', 'Cu', 'Au', 'Ag', 'Rh', 'Ir', 'Ru', 'Co',
    'Fe', 'Mn', 'Mo', 'W', 'Re', 'Ti', 'Zr', 'V', 'Cr', 'Nb',
    'Ta', 'Hf', 'Os', 'Zn', 'Al', 'Si', 'Ga', 'Ge', 'Sn', 'Sb',
    'In', 'Bi', 'Sc', 'Y',
}

# Tier 2: Moderately sampled (1k-50k slabs or via OC22 oxides)
_TIER2_ELEMENTS = {
    'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Ca', 'Sr', 'Ba', 'Mg', 'Li',
    'Na', 'K', 'Pb', 'Te', 'Se', 'Cd',
    'Rb', 'Cs',  # Alkali metals in some OC22 oxides
    'Eu', 'Gd', 'Dy', 'Er', 'Yb',  # Rare earths in OC22 oxide surfaces
}

# Element coverage score
_ELEMENT_COVERAGE = {}
for e in _TIER1_ELEMENTS:
    _ELEMENT_COVERAGE[e] = 1.0
for e in _TIER2_ELEMENTS:
    _ELEMENT_COVERAGE[e] = 0.6
# Everything else defaults to 0.1 (barely or never seen)


def element_coverage_score(elements: List[str]) -> float:
    """
    Score how well the elements in a genome are covered by OC20/OC25.
    
    Returns: float in [0, 1], where 1.0 = all elements extensively trained,
             0.1 = elements never in training data.
    """
    if not elements:
        return 0.5  # no metals (e.g., MetalFreeCarbon) — moderate uncertainty
    
    scores = [_ELEMENT_COVERAGE.get(e, 0.1) for e in elements if e != 'None']
    if not scores:
        return 0.5
    
    # Geometric mean: one bad element tanks the whole score
    return float(np.exp(np.mean(np.log(np.array(scores) + 1e-10))))


# ═══════════════════════════════════════════════════════════════════════════════
# LAYER 3: DUAL-MODEL DISAGREEMENT (cost = 2× per candidate)
# ═══════════════════════════════════════════════════════════════════════════════

def dual_model_disagreement(atoms, primary_calc, fallback_calc,
                            property_name: str = 'energy') -> Dict:
    """
    Compare predictions from two calculators on the same structure.
    
    Returns dict with:
      - 'primary_energy': energy from primary model (eSen-SM)
      - 'fallback_energy': energy from fallback model (MACE/EquiformerV2)
      - 'disagreement_eV': absolute difference
      - 'relative_disagreement': |diff| / mean(|energies|)
      - 'is_ood': True if disagreement exceeds threshold
    
    Should only be called on top-k candidates during validation rounds.
    """
    from ase.optimize import BFGS
    
    result = {}
    
    try:
        # Primary model (eSen-SM)
        atoms_p = atoms.copy()
        atoms_p.calc = primary_calc
        BFGS(atoms_p, logfile=None).run(fmax=0.1, steps=50)
        e_primary = atoms_p.get_potential_energy()
        result['primary_energy'] = float(e_primary)
    except Exception as e:
        result['primary_error'] = str(e)[:100]
        e_primary = None

    try:
        # Fallback model (MACE-MP-0 or second eSen variant)
        atoms_f = atoms.copy()
        atoms_f.calc = fallback_calc
        BFGS(atoms_f, logfile=None).run(fmax=0.1, steps=50)
        e_fallback = atoms_f.get_potential_energy()
        result['fallback_energy'] = float(e_fallback)
    except Exception as e:
        result['fallback_error'] = str(e)[:100]
        e_fallback = None

    if e_primary is not None and e_fallback is not None:
        diff = abs(e_primary - e_fallback)
        mean_abs = (abs(e_primary) + abs(e_fallback)) / 2.0
        rel_diff = diff / max(mean_abs, 0.01)
        
        result['disagreement_eV'] = float(diff)
        result['relative_disagreement'] = float(rel_diff)
        # OOD threshold: >2 eV absolute OR >50% relative
        result['is_ood'] = bool(diff > 2.0 or rel_diff > 0.5)
    else:
        result['disagreement_eV'] = float('inf')
        result['relative_disagreement'] = float('inf')
        result['is_ood'] = True  # can't verify → assume OOD

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# COMBINED CONFIDENCE SCORE
# ═══════════════════════════════════════════════════════════════════════════════

def compute_model_confidence(genome: tuple, elements: List[str],
                             dual_result: Optional[Dict] = None) -> float:
    """
    Compute combined confidence score for a screening prediction.
    
    Combines:
      - Layer 1: class-level prior (40% weight)
      - Layer 2: element coverage (30% weight)
      - Layer 3: dual-model disagreement if available (30% weight)
    
    Returns: float in [0, 1], where 1.0 = high confidence, 0.0 = don't trust.
    """
    mat_class = genome[0]
    
    # Layer 1: class prior
    cls_conf = CLASS_CONFIDENCE.get(mat_class, 0.3)
    
    # Layer 2: element coverage
    elem_conf = element_coverage_score(elements)
    
    if dual_result is not None and 'relative_disagreement' in dual_result:
        # Layer 3 available: use all three
        rel_dis = dual_result['relative_disagreement']
        # Convert disagreement to confidence: 0% disagreement → 1.0, 100% → 0.0
        dual_conf = max(0.0, 1.0 - rel_dis)
        
        confidence = 0.50 * cls_conf + 0.25 * elem_conf + 0.25 * dual_conf
    else:
        # No dual model: class prior dominates (it captures structural coverage)
        confidence = 0.65 * cls_conf + 0.35 * elem_conf
    
    return float(np.clip(confidence, 0.0, 1.0))


def confidence_penalty(confidence: float) -> float:
    """
    Convert model confidence to NSGA-II objective penalty.
    
    Applied ADDITIVELY to predicted properties:
      confidence 1.0 → penalty 0.0 (no change)
      confidence 0.5 → penalty 0.5 (shifted by 0.5 V)
      confidence 0.2 → penalty 0.8 (shifted by 0.8 V)
      confidence 0.0 → penalty 1.0 (maximum shift)
    
    Additive penalty ensures OOD candidates get worse rankings regardless
    of whether the surrogate prediction is positive or negative.
    This makes the GA prefer candidates where the model is trustworthy,
    without completely eliminating OOD candidates from exploration.
    """
    return 1.0 - confidence
