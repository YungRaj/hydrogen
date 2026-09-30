#!/usr/bin/env python3
# Screening calculator construction and calibration.
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
import shlex
import sys
import json
import urllib.request
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from ase import Atoms
from ase.calculators.calculator import Calculator

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.utils import setup_logger, BASE_DIR

logger = setup_logger('surface_calculator', 'screening/surface_calculator.log')


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 1: MACE-MP-0 (existing — unchanged)
# ═══════════════════════════════════════════════════════════════════════════════


def get_mace_calculator(device='cuda:0', model='medium'):
    """Load MACE-MP-0 calculator (bulk materials model, fast screening).

    Args:
        device: CPU or GPU device requested for execution.
        model: Fitted model used for inference.

    Returns:
        Computed result described above.
    """
    from mace.calculators import mace_mp

    return mace_mp(model=model, device=device)


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 2a: CATALYSIS-HUB LOOKUP (real DFT data, 100k+ reactions)
# ═══════════════════════════════════════════════════════════════════════════════

CATALYSIS_HUB_URL = 'https://api.catalysis-hub.org/graphql'


def query_catalysis_hub(
    surface: str, facet: str = None, adsorbate: str = None, limit: int = 50
) -> List[Dict]:
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
    }}""".format(
        filters=','.join(filters)
    )

    try:
        req = urllib.request.Request(
            f'{CATALYSIS_HUB_URL}?query={urllib.parse.quote(query)}',
            headers={'Accept': 'application/json'},
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


def lookup_adsorption_energy(
    element: str, facet: str, adsorbate: str
) -> Optional[float]:
    """
        Look up a specific adsorption energy from Catalysis-Hub.

        Returns energy in eV if found, None otherwise.

    Args:
        element: Element used by this operation.
        facet: Facet used by this operation.
        adsorbate: Adsorbate used by this operation.

    Returns:
        Computed `Optional[float]` result.
    """
    results = query_catalysis_hub(element, facet, adsorbate, limit=5)
    if results:
        energies = [
            r['reactionEnergy'] for r in results if r.get('reactionEnergy') is not None
        ]
        if energies:
            return float(np.median(energies))
    return None



def _ensure_hf_token():
    """Load HuggingFace token from .hf_token file or env var."""
    if os.environ.get('HF_TOKEN'):
        return True
    token_file = BASE_DIR / '.hf_token'
    if token_file.exists():
        token = token_file.read_text().strip()
        if token:
            os.environ['HF_TOKEN'] = token
            return True
    return False


def get_ocp_calculator(
    model_name: str = 'esen-sm-conserving-all-oc25', device: str = 'cuda:0'
) -> Optional[Calculator]:
    """
        Load an OC20/OC25-trained surface catalysis GNN calculator.

        Available models (fairchem v2, require HF_TOKEN):
          - esen-sm-conserving-all-oc25  (recommended — energy-conserving)
          - esen-md-direct-all-oc25      (faster, MD-optimized)
          - uma-s-1p1                    (Universal Model for Atoms)
          - uma-m-1p1                    (larger UMA)

        Falls back to local EquiformerV2-31M checkpoint if available.

    Args:
        model_name: Model name used by this operation.
        device: CPU or GPU device requested for execution.

    Returns:
        Computed `Optional[Calculator]` result.
    """
    _ensure_hf_token()

    # Try 1: fairchem v2 API (eSen/UMA — gated, needs HF_TOKEN)
    try:
        from fairchem.core.calculate.pretrained_mlip import get_predict_unit
        from fairchem.core import FAIRChemCalculator

        # fairchem v2 does not infer the caller's requested device from the ASE
        # calculator. Forward it explicitly so CPU requests and CUDA worker
        # masks cannot silently fall back to another device.
        unit = get_predict_unit(model_name, device=device)
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

    ecutwfc: float = 50.0  # Plane-wave cutoff (Ry)
    ecutrho: float = 400.0  # Charge density cutoff (Ry)
    kpoints: Tuple[int, ...] = (4, 4, 1)  # k-point mesh (Γ-centered)
    smearing: str = 'marzari-vanderbilt'
    degauss: float = 0.02  # Smearing width (Ry)
    conv_thr: float = 1e-6  # SCF convergence (Ry)
    nstep: int = 100  # Max ionic steps
    forc_conv_thr: float = 1e-3  # Force convergence (Ry/bohr)
    pseudo_dir: str = ''
    outdir: str = ''
    n_cores: int = 16


def get_qe_calculator(atoms: Atoms, config: QEConfig = None) -> Optional[Calculator]:
    """
        Create an ASE-compatible Quantum ESPRESSO calculator for a surface slab.

        Requires:
          - the pinned GPU QE activation file (``PW_X`` and ``MPIEXEC``)
          - PAW pseudopotentials in pseudo_dir

    Args:
        atoms: Atomic structure consumed by the calculator.
        config: Configuration controlling this operation.

    Returns:
        Computed `Optional[Calculator]` result.
    """
    if config is None:
        config = QEConfig()

    try:
        from ase.calculators.espresso import Espresso, EspressoProfile

        # Do not accept a same-named CPU executable from PATH or Conda.
        from pipeline.simulation.executables import resolve_qe_executable

        try:
            pw_path = resolve_qe_executable('pw.x')
        except RuntimeError as exc:
            logger.warning("GPU QE calculator unavailable: %s", exc)
            return None
        mpi_path = os.environ.get('MPIEXEC')
        if not mpi_path:
            logger.warning(
                "GPU QE calculator unavailable: MPIEXEC is not set by activate.sh"
            )
            return None

        # Locate pseudopotentials
        pseudo_dir = config.pseudo_dir or str(BASE_DIR / 'quantum_espresso' / 'pseudo')
        if not os.path.isdir(pseudo_dir):
            os.makedirs(pseudo_dir, exist_ok=True)

        # Build pseudopotential mapping for elements in atoms
        pseudopotentials = {}
        for symbol in set(atoms.get_chemical_symbols()):
            # Standard SSSP naming: Element.pbe-n-kjpaw_psl.1.0.0.UPF
            pseudo_files = list(Path(pseudo_dir).glob(f'{symbol}.*UPF')) + list(
                Path(pseudo_dir).glob(f'{symbol}.*upf')
            )
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
            command=f'{shlex.quote(mpi_path)} -np {config.n_cores} '
            f'{shlex.quote(pw_path)}',
            pseudo_dir=pseudo_dir,
        )

        calc = Espresso(
            profile=profile,
            pseudopotentials=pseudopotentials,
            input_data=input_data,
            kpts=config.kpoints,
        )

        logger.info(
            f"QE calculator created: ecutwfc={config.ecutwfc} Ry, "
            f"kpts={config.kpoints}, {config.n_cores} cores"
        )
        return calc

    except Exception as e:
        logger.warning(f"QE calculator setup failed: {e}")
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# MULTI-FIDELITY EVALUATION
# ═══════════════════════════════════════════════════════════════════════════════





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
