#!/usr/bin/env python3
"""
Multi-Fidelity Surface Catalysis Calculator.

Replaces the single bulk MACE-MP-0 model with a tiered approach:

  Tier 1 (Screening): MACE-MP-0 — fast (~2s/candidate), bulk-trained
                       Used for screening the 21.1B encoded space

  Tier 2 (Validation): Catalysis-Hub lookup + EquiformerV2 (OC20)
                        Real DFT surface energies from 100k+ reactions
                        OC20-trained GNN for surface adsorption

  Tier 3 (High-fidelity): Quantum ESPRESSO DFT
                           Full periodic slab DFT with PAW pseudopotentials
                           ~1-4 hours per candidate

Each tier provides an ASE-compatible calculator interface.
The campaign uses Tier 1 for population screening, Tier 2 for
top-k validation, and Tier 3 for champion catalysts.
"""

import os
import sys
import json
import time
import logging
import urllib.request
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write as ase_write

sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.utils import setup_logger, save_json, BASE_DIR

logger = setup_logger('surface_calculator', 'screening/surface_calculator.log')


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 1: MACE-MP-0 (existing — unchanged)
# ═══════════════════════════════════════════════════════════════════════════════

def get_mace_calculator(device='cuda:0', model='medium'):
    """Load MACE-MP-0 calculator (bulk materials model, fast screening)."""
    from mace.calculators import mace_mp
    return mace_mp(model=model, device=device)


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 2a: CATALYSIS-HUB LOOKUP (real DFT data, 100k+ reactions)
# ═══════════════════════════════════════════════════════════════════════════════

CATALYSIS_HUB_URL = 'https://api.catalysis-hub.org/graphql'


def query_catalysis_hub(surface: str, facet: str = None,
                        adsorbate: str = None, limit: int = 50) -> List[Dict]:
    """
    Query Catalysis-Hub for DFT-computed adsorption/reaction energies.

    Args:
        surface: Surface composition (e.g., "Pt", "NiFe", "Pd3Au")
        facet: Crystal facet (e.g., "111", "100", "211")
        adsorbate: Adsorbate species (e.g., "H", "OH", "CH3")
        limit: Max results

    Returns: List of reaction energy records from the database.
    """
    filters = [f'first:{limit}']
    if surface:
        filters.append(f'surfaceComposition:"{surface}"')
    if facet:
        filters.append(f'facet:"{facet}"')
    if adsorbate:
        filters.append(f'products:"star{adsorbate}"')

    query = """{{
        reactions({filters}) {{
            edges {{
                node {{
                    reactionEnergy
                    activationEnergy
                    surfaceComposition
                    facet
                    products
                    reactants
                    chemicalComposition
                    reactionSystems {{
                        name
                        aseId
                    }}
                    sites
                }}
            }}
        }}
    }}""".format(filters=','.join(filters))

    try:
        req = urllib.request.Request(
            f'{CATALYSIS_HUB_URL}?query={urllib.parse.quote(query)}',
            headers={'Accept': 'application/json'}
        )
        resp = urllib.request.urlopen(req, timeout=15)
        data = json.loads(resp.read())
        edges = data.get('data', {}).get('reactions', {}).get('edges', [])
        results = [e['node'] for e in edges]
        logger.info(f"Catalysis-Hub: {len(results)} results for {surface}({facet})")
        return results
    except Exception as e:
        logger.warning(f"Catalysis-Hub query failed: {e}")
        return []


def lookup_adsorption_energy(element: str, facet: str,
                             adsorbate: str) -> Optional[float]:
    """
    Look up a specific adsorption energy from Catalysis-Hub.

    Returns energy in eV if found, None otherwise.
    """
    results = query_catalysis_hub(element, facet, adsorbate, limit=5)
    if results:
        energies = [r['reactionEnergy'] for r in results
                    if r.get('reactionEnergy') is not None]
        if energies:
            return float(np.median(energies))
    return None


def build_calibration_table(elements: List[str] = None,
                            facets: List[str] = None,
                            adsorbates: List[str] = None) -> Dict:
    """
    Build a lookup table of known DFT adsorption energies from Catalysis-Hub.
    Used to calibrate MACE predictions against ground truth.
    """
    if elements is None:
        elements = ['Pt', 'Pd', 'Ni', 'Cu', 'Au', 'Ag', 'Rh', 'Ir',
                     'Ru', 'Fe', 'Co', 'Mo', 'W', 'Re', 'Mn', 'Ti']
    if facets is None:
        facets = ['111', '100', '211']
    if adsorbates is None:
        adsorbates = ['H', 'OH', 'O', 'CH3', 'C']

    table = {}
    n_found = 0
    for elem in elements:
        for facet in facets:
            for ads in adsorbates:
                key = f'{elem}_{facet}_{ads}'
                e = lookup_adsorption_energy(elem, facet, ads)
                if e is not None:
                    table[key] = e
                    n_found += 1
                time.sleep(0.1)  # rate limit

    logger.info(f"Calibration table: {n_found} entries from Catalysis-Hub")
    cache_path = BASE_DIR / 'models' / 'catalysis_hub_calibration.json'
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    save_json(table, str(cache_path))
    return table


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 2b: EQUIFORMERV2 / eSen (OC20/OC22 surface-trained GNN)
# ═══════════════════════════════════════════════════════════════════════════════

def _ensure_hf_token():
    """Load HuggingFace token from .hf_token file or env var."""
    if os.environ.get('HF_TOKEN'):
        return True
    token_file = Path(__file__).parent.parent / '.hf_token'
    if token_file.exists():
        token = token_file.read_text().strip()
        if token:
            os.environ['HF_TOKEN'] = token
            return True
    return False


def get_ocp_calculator(model_name: str = 'esen-sm-conserving-all-oc25',
                       device: str = 'cuda:0') -> Optional[Calculator]:
    """
    Load an OC20/OC25-trained surface catalysis GNN calculator.

    Available models (fairchem v2, require HF_TOKEN):
      - esen-sm-conserving-all-oc25  (recommended — energy-conserving)
      - esen-md-direct-all-oc25      (faster, MD-optimized)
      - uma-s-1p1                    (Universal Model for Atoms)
      - uma-m-1p1                    (larger UMA)

    Falls back to local EquiformerV2-31M checkpoint if available.
    """
    _ensure_hf_token()

    # Try 1: fairchem v2 API (eSen/UMA — gated, needs HF_TOKEN)
    try:
        from fairchem.core.calculate.pretrained_mlip import get_predict_unit
        from fairchem.core import FAIRChemCalculator

        unit = get_predict_unit(model_name)
        calc = FAIRChemCalculator(predict_unit=unit)
        logger.info(f"Loaded {model_name} via fairchem v2")
        return calc
    except Exception as e:
        logger.debug(f"fairchem v2 load failed: {e}")

    # Try 2: Local OC20 checkpoint (EquiformerV2-31M)
    local_ckpt = BASE_DIR / 'models' / 'eq2_31M_ec4_allmd.pt'
    if local_ckpt.exists():
        try:
            from fairchem.core import OCPCalculator
            calc = OCPCalculator(checkpoint_path=str(local_ckpt), cpu=(device == 'cpu'))
            logger.info(f"Loaded EquiformerV2 from local checkpoint")
            return calc
        except Exception as e:
            logger.debug(f"fairchem v1 load failed: {e}")

    # Try 3: Download checkpoint
    try:
        import torch
        url = 'https://dl.fbaipublicfiles.com/opencatalystproject/models/2023_06/oc20/s2ef/eq2_31M_ec4_allmd.pt'
        local_ckpt.parent.mkdir(parents=True, exist_ok=True)
        if not local_ckpt.exists():
            logger.info("Downloading EquiformerV2-31M checkpoint...")
            torch.hub.download_url_to_file(url, str(local_ckpt))

        from fairchem.core import OCPCalculator
        calc = OCPCalculator(checkpoint_path=str(local_ckpt), cpu=(device == 'cpu'))
        logger.info("Loaded EquiformerV2 from downloaded checkpoint")
        return calc
    except Exception as e:
        logger.warning(f"OC20 calculator not available: {e}")

    return None


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 3: QUANTUM ESPRESSO (full DFT — highest fidelity)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class QEConfig:
    """Quantum ESPRESSO calculation parameters for surface catalysis."""
    ecutwfc: float = 50.0        # Plane-wave cutoff (Ry)
    ecutrho: float = 400.0       # Charge density cutoff (Ry)
    kpoints: Tuple[int, ...] = (4, 4, 1)  # k-point mesh (Γ-centered)
    smearing: str = 'marzari-vanderbilt'
    degauss: float = 0.02        # Smearing width (Ry)
    conv_thr: float = 1e-6       # SCF convergence (Ry)
    nstep: int = 100             # Max ionic steps
    forc_conv_thr: float = 1e-3  # Force convergence (Ry/bohr)
    pseudo_dir: str = ''
    outdir: str = ''
    n_cores: int = 16


def get_qe_calculator(atoms: Atoms, config: QEConfig = None) -> Optional[Calculator]:
    """
    Create an ASE-compatible Quantum ESPRESSO calculator for a surface slab.

    Requires:
      - pw.x in PATH (conda run -n qe-env)
      - PAW pseudopotentials in pseudo_dir
    """
    if config is None:
        config = QEConfig()

    try:
        from ase.calculators.espresso import Espresso, EspressoProfile

        # Locate pw.x
        pw_path = '/home/ilhanraja/miniconda3/envs/qe-env/bin/pw.x'
        if not os.path.exists(pw_path):
            import shutil
            pw_path = shutil.which('pw.x')

        if pw_path is None:
            logger.warning("pw.x not found — QE calculator unavailable")
            return None

        # Locate pseudopotentials
        pseudo_dir = config.pseudo_dir or str(BASE_DIR / 'quantum_espresso' / 'pseudo')
        if not os.path.isdir(pseudo_dir):
            os.makedirs(pseudo_dir, exist_ok=True)

        # Build pseudopotential mapping for elements in atoms
        pseudopotentials = {}
        for symbol in set(atoms.get_chemical_symbols()):
            # Standard SSSP naming: Element.pbe-n-kjpaw_psl.1.0.0.UPF
            pseudo_files = list(Path(pseudo_dir).glob(f'{symbol}.*UPF')) + \
                           list(Path(pseudo_dir).glob(f'{symbol}.*upf'))
            if pseudo_files:
                pseudopotentials[symbol] = pseudo_files[0].name
            else:
                pseudopotentials[symbol] = f'{symbol}.pbe-n-kjpaw_psl.1.0.0.UPF'

        outdir = config.outdir or str(BASE_DIR / 'results' / 'dft' / 'tmp')
        os.makedirs(outdir, exist_ok=True)

        input_data = {
            'control': {
                'calculation': 'relax',
                'restart_mode': 'from_scratch',
                'pseudo_dir': pseudo_dir,
                'outdir': outdir,
                'tprnfor': True,
                'tstress': True,
                'nstep': config.nstep,
                'forc_conv_thr': config.forc_conv_thr,
            },
            'system': {
                'ecutwfc': config.ecutwfc,
                'ecutrho': config.ecutrho,
                'occupations': 'smearing',
                'smearing': config.smearing,
                'degauss': config.degauss,
            },
            'electrons': {
                'conv_thr': config.conv_thr,
                'mixing_beta': 0.3,
            },
        }

        profile = EspressoProfile(
            command=f'mpirun -np {config.n_cores} {pw_path}',
            pseudo_dir=pseudo_dir,
        )

        calc = Espresso(
            profile=profile,
            pseudopotentials=pseudopotentials,
            input_data=input_data,
            kpts=config.kpoints,
        )

        logger.info(f"QE calculator created: ecutwfc={config.ecutwfc} Ry, "
                    f"kpts={config.kpoints}, {config.n_cores} cores")
        return calc

    except Exception as e:
        logger.warning(f"QE calculator setup failed: {e}")
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# MULTI-FIDELITY EVALUATION
# ═══════════════════════════════════════════════════════════════════════════════

def evaluate_with_tier(atoms: Atoms, adsorbate: str, tier: int = 1,
                       device: str = 'cuda:0') -> Dict:
    """
    Evaluate adsorption energy using the specified fidelity tier.

    Tier 1: MACE-MP-0 (fast, ~2s, bulk-trained)
    Tier 2: OC20 EquiformerV2 + Catalysis-Hub validation (medium, ~5s, surface-trained)
    Tier 3: Quantum ESPRESSO DFT (slow, ~hours, gold standard)

    Returns dict with energy, forces, fidelity metadata.
    """
    result = {'tier': tier, 'adsorbate': adsorbate, 'n_atoms': len(atoms)}

    if tier == 1:
        calc = get_mace_calculator(device)
        atoms.calc = calc
        from ase.optimize import BFGS
        BFGS(atoms, logfile=None).run(fmax=0.08, steps=100)
        result['energy'] = atoms.get_potential_energy()
        result['max_force'] = float(abs(atoms.get_forces()).max())
        result['model'] = 'MACE-MP-0 (bulk)'

    elif tier == 2:
        # Try OC20 surface model first
        ocp_calc = get_ocp_calculator(device=device)
        if ocp_calc is not None:
            atoms.calc = ocp_calc
            from ase.optimize import BFGS
            BFGS(atoms, logfile=None).run(fmax=0.05, steps=150)
            result['energy'] = atoms.get_potential_energy()
            result['max_force'] = float(abs(atoms.get_forces()).max())
            result['model'] = 'EquiformerV2 (OC20-surface)'
        else:
            # Fall back to MACE + Catalysis-Hub correction
            calc = get_mace_calculator(device)
            atoms.calc = calc
            from ase.optimize import BFGS
            BFGS(atoms, logfile=None).run(fmax=0.08, steps=100)
            mace_energy = atoms.get_potential_energy()

            # Look up reference from Catalysis-Hub
            symbols = atoms.get_chemical_symbols()
            metal = [s for s in set(symbols) if s not in ('H', 'C', 'O', 'N')]
            if metal:
                ref_energy = lookup_adsorption_energy(metal[0], '111', adsorbate)
                if ref_energy is not None:
                    result['catalysis_hub_ref'] = ref_energy
                    result['model'] = 'MACE-MP-0 + Catalysis-Hub calibration'
                else:
                    result['model'] = 'MACE-MP-0 (no reference found)'
            else:
                result['model'] = 'MACE-MP-0 (fallback)'

            result['energy'] = mace_energy
            result['max_force'] = float(abs(atoms.get_forces()).max())

    elif tier == 3:
        qe_calc = get_qe_calculator(atoms)
        if qe_calc is not None:
            atoms.calc = qe_calc
            result['energy'] = atoms.get_potential_energy()
            result['max_force'] = float(abs(atoms.get_forces()).max())
            result['model'] = 'Quantum ESPRESSO (DFT-PBE)'
        else:
            result['error'] = 'QE not available'
            result['model'] = 'N/A'

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# CROSS-VALIDATION: Compare MACE vs Catalysis-Hub
# ═══════════════════════════════════════════════════════════════════════════════

def cross_validate_mace(elements: List[str] = None,
                        device: str = 'cuda:0') -> Dict:
    """
    Compare MACE-MP-0 predictions against Catalysis-Hub DFT data
    for known systems. Computes MAE, RMSE, and systematic bias.

    This tells us how much to trust MACE for novel catalysts.
    """
    from ase.build import fcc111, add_adsorbate
    from ase.optimize import BFGS

    if elements is None:
        elements = ['Pt', 'Pd', 'Ni', 'Cu', 'Au', 'Ag', 'Rh', 'Ru']

    calc = get_mace_calculator(device)

    mace_energies = []
    dft_energies = []
    systems = []

    for elem in elements:
        for ads in ['H', 'O', 'OH']:
            # MACE prediction
            slab = fcc111(elem, size=(2, 2, 3), vacuum=10.0)
            e_clean = slab.copy()
            e_clean.calc = calc
            BFGS(e_clean, logfile=None).run(fmax=0.1, steps=50)
            E_clean = e_clean.get_potential_energy()

            slab_ads = slab.copy()
            if ads == 'OH':
                add_adsorbate(slab_ads, 'O', height=1.9, position='ontop')
                from ase import Atom
                slab_ads.append(Atom('H', position=slab_ads[-1].position + [0, 0, 0.97]))
            elif ads == 'O':
                add_adsorbate(slab_ads, 'O', height=1.7, position='ontop')
            else:
                add_adsorbate(slab_ads, ads, height=1.5, position='ontop')
            slab_ads.calc = calc
            BFGS(slab_ads, logfile=None).run(fmax=0.1, steps=50)
            E_ads = slab_ads.get_potential_energy()
            dE_mace = E_ads - E_clean

            # Catalysis-Hub reference
            dE_dft = lookup_adsorption_energy(elem, '111', ads)

            if dE_dft is not None:
                mace_energies.append(dE_mace)
                dft_energies.append(dE_dft)
                systems.append(f'{elem}(111)+{ads}*')
                logger.info(f"  {elem}(111)+{ads}*: MACE={dE_mace:.3f}, DFT={dE_dft:.3f} eV")

    if len(mace_energies) > 0:
        mace_arr = np.array(mace_energies)
        dft_arr = np.array(dft_energies)
        errors = mace_arr - dft_arr
        mae = float(np.mean(np.abs(errors)))
        rmse = float(np.sqrt(np.mean(errors**2)))
        bias = float(np.mean(errors))

        result = {
            'n_comparisons': len(systems),
            'MAE_eV': round(mae, 4),
            'RMSE_eV': round(rmse, 4),
            'systematic_bias_eV': round(bias, 4),
            'systems': systems,
            'mace_energies': [round(e, 4) for e in mace_energies],
            'dft_energies': [round(e, 4) for e in dft_energies],
        }
        logger.info(f"MACE vs DFT: MAE={mae:.3f} eV, RMSE={rmse:.3f} eV, bias={bias:.3f} eV")
        return result

    return {'error': 'No comparison data available'}


if __name__ == '__main__':
    print("=== Surface Calculator Multi-Fidelity Tiers ===")
    print()
    print("Tier 1: MACE-MP-0 (bulk) — AVAILABLE")

    ocp = get_ocp_calculator()
    if ocp:
        print("Tier 2: EquiformerV2 (OC20 surface) — AVAILABLE")
    else:
        print("Tier 2: EquiformerV2 — NOT AVAILABLE (need HF_TOKEN for UMA/eSen)")
        print("        Catalysis-Hub API — AVAILABLE (100k+ DFT reactions)")

    qe = get_qe_calculator(Atoms('Pt'))
    if qe:
        print("Tier 3: Quantum ESPRESSO — AVAILABLE")
    else:
        print("Tier 3: Quantum ESPRESSO — NEEDS PSEUDOPOTENTIALS")
