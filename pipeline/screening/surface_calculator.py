#!/usr/bin/env python3
# Screening calculator construction and calibration.
"""
Multi-Fidelity Surface Catalysis Calculator.

  Tier 1 (Screening): eSEN/UMA surface GNN
  Tier 2 (High-fidelity): Quantum ESPRESSO DFT
                           Full periodic slab DFT with PAW pseudopotentials
                           ~1-4 hours per candidate
"""

import os
import hashlib
import shlex
import sys
from pathlib import Path
from typing import Optional, Tuple
from dataclasses import dataclass

from ase import Atoms
from ase.calculators.calculator import Calculator

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.utils import setup_logger, BASE_DIR

logger = setup_logger('surface_calculator', 'screening/surface_calculator.log')


def _verify_checkpoint(path: Path) -> None:
    """Verify an optional checkpoint digest before loading model weights.

    The expected digest is intentionally supplied by the deployment environment
    because gated Fairchem weights are not stored in this repository.
    """
    expected = os.environ.get('FAIRCHEM_CHECKPOINT_SHA256', '').strip().lower()
    if not expected:
        return
    if len(expected) != 64 or any(char not in '0123456789abcdef' for char in expected):
        raise RuntimeError('FAIRCHEM_CHECKPOINT_SHA256 must be a 64-character hex digest')
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    actual = digest.hexdigest()
    if actual != expected:
        raise RuntimeError(
            f'Fairchem checkpoint digest mismatch for {path}: '
            f'expected {expected}, got {actual}'
        )


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 1: SURFACE GNN
# ═══════════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════════════

def _token_files() -> list[Path]:
    """Token locations in priority order: HF cache standard, then repo file."""
    hf_home = Path(os.environ.get('HF_HOME', Path.home() / '.cache' / 'huggingface'))
    return [hf_home / 'token', BASE_DIR / '.hf_token']


def _ensure_hf_token():
    """Load the HuggingFace token from the environment or a token file.

    Checks ``HF_TOKEN``/``HUGGING_FACE_HUB_TOKEN``, then ``$HF_HOME/token``
    (where ``huggingface-cli login`` writes it), then the repo ``.hf_token``.
    """
    for name in ('HF_TOKEN', 'HUGGING_FACE_HUB_TOKEN'):
        if os.environ.get(name):
            os.environ['HF_TOKEN'] = os.environ[name]
            return True
    for token_file in _token_files():
        if token_file.is_file():
            token = token_file.read_text().strip()
            if token:
                os.environ['HF_TOKEN'] = token
                return True
    return False


_DEFAULT_MODEL_NAME = 'esen-sm-conserving-all-oc25'


def screening_model_name() -> str:
    """Return the fairchem model used for production screening.

    ``HYDROGEN_SCREENING_MODEL`` selects an alternative checkpoint (e.g.
    ``uma-m-1p1``) for model-comparison studies; results record the model and
    the result cache keys on it, so outputs from different models never mix.

    Returns:
        Fairchem model identifier.
    """
    return os.environ.get('HYDROGEN_SCREENING_MODEL', _DEFAULT_MODEL_NAME)


def get_ocp_calculator(
    model_name: Optional[str] = None, device: str = 'cuda:0'
) -> Optional[Calculator]:
    """
        Load an OC20/OC25-trained surface catalysis GNN calculator.

        Available models (fairchem v2, require HF_TOKEN):
          - esen-sm-conserving-all-oc25  (recommended — energy-conserving)
          - esen-md-direct-all-oc25      (faster, MD-optimized)
          - uma-s-1p1                    (Universal Model for Atoms)
          - uma-m-1p1                    (larger UMA)

        When no model is named, the default eSEN model is tried first and the
        local EquiformerV2-31M checkpoint is used as a fallback.  When a model
        is named explicitly, no other model is ever substituted: results are
        labelled with that model's protocol, so a failed load returns None.

    Args:
        model_name: Explicit model to load, or None for the default eSEN model
            with the legacy EquiformerV2 fallback.
        device: CPU or GPU device requested for execution.

    Returns:
        The requested calculator, or None when it cannot be loaded.
    """
    allow_fallback = model_name is None
    if model_name is None:
        model_name = _DEFAULT_MODEL_NAME
    _ensure_hf_token()

    # Try 1: fairchem v2 API (eSen/UMA — gated, needs HF_TOKEN)
    try:
        from fairchem.core.calculate.pretrained_mlip import get_predict_unit
        from fairchem.core import FAIRChemCalculator

        # fairchem v2 does not infer the caller's requested device from the ASE
        # calculator. Forward it explicitly so CPU requests and CUDA worker
        # masks cannot silently fall back to another device.
        # Recent Fairchem releases accept only ``cuda``/``cpu`` here. GPU
        # worker affinity is established before this call, so normalizing
        # ``cuda:N`` does not change which device the worker owns.
        fairchem_device = 'cuda' if device.startswith('cuda') else device
        # 'default' keeps the validated screening numerics; 'batch' is
        # fairchem's mode for heterogeneous batched evaluation and can be
        # benchmarked via HYDROGEN_FAIRCHEM_INFERENCE without code changes.
        inference = os.environ.get('HYDROGEN_FAIRCHEM_INFERENCE', 'default')
        unit = get_predict_unit(
            model_name, inference_settings=inference, device=fairchem_device)
        # Fairchem v2 requires the task explicitly.  eSEN OC25 checkpoints
        # expose ``oc25``; UMA surface checkpoints expose ``oc20``.
        task_name = 'oc25' if 'oc25' in model_name else 'oc20'
        calc = FAIRChemCalculator(predict_unit=unit, task_name=task_name)
        logger.info(f"Loaded {model_name} via fairchem v2")
        return calc
    except Exception as e:
        if not allow_fallback:
            logger.warning(
                f"Requested model {model_name} failed to load ({e}); "
                "not substituting another model"
            )
            return None
        logger.debug(f"fairchem v2 load failed: {e}")

    # Try 2: Local OC20 checkpoint (EquiformerV2-31M)
    local_ckpt = BASE_DIR / 'models' / 'eq2_31M_ec4_allmd.pt'
    if local_ckpt.exists():
        try:
            _verify_checkpoint(local_ckpt)
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

        _verify_checkpoint(local_ckpt)

        from fairchem.core import OCPCalculator

        calc = OCPCalculator(checkpoint_path=str(local_ckpt), cpu=(device == 'cpu'))
        logger.info("Loaded EquiformerV2 from downloaded checkpoint")
        return calc
    except Exception as e:
        logger.warning(f"OC20 calculator not available: {e}")

    return None


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 2: QUANTUM ESPRESSO (full DFT — highest fidelity)
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
# COMMAND-LINE AVAILABILITY CHECK
# ═══════════════════════════════════════════════════════════════════════════════





if __name__ == '__main__':
    print("=== Surface Calculator Fidelity Tiers ===")
    ocp = get_ocp_calculator()
    if ocp:
        print("Tier 1: eSEN/UMA surface GNN — AVAILABLE")
    else:
        print("Tier 1: eSEN/UMA surface GNN — NOT AVAILABLE")

    qe = get_qe_calculator(Atoms('Pt'))
    if qe:
        print("Tier 2: Quantum ESPRESSO — AVAILABLE")
    else:
        print("Tier 2: Quantum ESPRESSO — NEEDS PSEUDOPOTENTIALS")
