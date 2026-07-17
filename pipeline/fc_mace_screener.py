#!/usr/bin/env python3
"""
Multi-GPU MACE-MP-0 Screening for ORR Fuel Cell Cathode Catalysts.

For each catalyst candidate from the 25.3B design space, computes:
  1. ΔG_OH* — hydroxyl adsorption free energy
  2. ΔG_O*  — oxygen adsorption free energy
  3. ΔG_OOH* — peroxyl adsorption free energy
  4. ORR overpotential (4-electron CHE method)
  5. Fenton susceptibility index (radical degradation risk)
  6. Dissolution stability estimate

Uses the SAME design space and structure generators as the methane
pyrolysis screener, but evaluates ORR-specific descriptors.

Parallelized across all available GPUs with 2 workers per device.
"""

import os
import sys
import time
import random
import numpy as np
import torch
import torch.multiprocessing as mp
from typing import List, Tuple, Dict, Optional

from ase import Atoms, Atom
from ase.optimize import BFGS

from pipeline.utils import (
    BASE_DIR, FUEL_CELL_DIR, setup_logger, print_banner,
    save_screening_db, orr_overpotential, abundance_cost_penalty,
)
from pipeline.mace_screener import generate_structure

logger = setup_logger('fc_mace_screener', 'fuel_cell/fc_mace_screening.log')


# ═══════════════════════════════════════════════════════════════════════════════
# ORR REFERENCE ENERGIES
# ═══════════════════════════════════════════════════════════════════════════════

# Zero-point energy and entropy corrections for ORR intermediates (eV)
# From standard DFT+CHE literature (Nørskov et al., J. Phys. Chem. B, 2004)
ZPE_CORRECTIONS = {
    'OH': 0.35,   # ZPE(OH*) - 0.5*ZPE(H2O)
    'O': 0.05,    # ZPE(O*)
    'OOH': 0.40,  # ZPE(OOH*)
}
TS_CORRECTIONS = {
    'OH': -0.07,
    'O': 0.00,
    'OOH': -0.10,
}

# Elements with known Fenton reactivity (radical generation in acid)
FENTON_RISK = {
    'Fe': 3, 'Cu': 2, 'Co': 1, 'Mn': 1, 'Cr': 1,
    'V': 1, 'Ti': 0, 'Ni': 0, 'Zn': 0, 'Mo': 0,
}


def compute_water_ref(calc) -> float:
    """Compute H₂O reference energy."""
    from ase.build import molecule
    h2o = molecule('H2O')
    h2o.set_cell([10, 10, 10])
    h2o.center()
    h2o.pbc = True
    h2o.calc = calc
    BFGS(h2o, logfile=None).run(fmax=0.05)
    return h2o.get_potential_energy()


def compute_h2_ref(calc) -> float:
    """Compute H₂ reference energy."""
    from ase.build import molecule
    h2 = molecule('H2')
    h2.set_cell([10, 10, 10])
    h2.center()
    h2.pbc = True
    h2.calc = calc
    BFGS(h2, logfile=None).run(fmax=0.05)
    return h2.get_potential_energy()


# ═══════════════════════════════════════════════════════════════════════════════
# ORR EVALUATION
# ═══════════════════════════════════════════════════════════════════════════════

def evaluate_orr_candidate(genome: tuple, calc, e_h2o: float, e_h2: float) -> dict:
    """
    Evaluate a single catalyst candidate for ORR activity using MACE-MP-0.

    Computes adsorption free energies for OH*, O*, OOH* intermediates
    and derives the theoretical ORR overpotential via the CHE method.
    """
    mat_class = genome[0]
    result = {
        'genome': str(genome),
        'material_class': mat_class,
        'valid': False,
    }

    try:
        structure, active_idx, _ = generate_structure(genome)

        # 1. Relax clean surface
        structure.calc = calc
        BFGS(structure, logfile=None).run(fmax=0.08, steps=150)
        e_clean = structure.get_potential_energy()
        result['e_clean'] = e_clean

        # Active site position
        if len(active_idx) > 0 and active_idx[0] < len(structure):
            ads_base = structure[active_idx[0]].position.copy()
        else:
            ads_base = structure.positions.mean(axis=0)

        # 2. OH* adsorption
        slab_oh = structure.copy()
        oh_pos = ads_base + np.array([0.0, 0.0, 1.9])
        slab_oh.append(Atom('O', position=oh_pos))
        slab_oh.append(Atom('H', position=oh_pos + np.array([0.0, 0.0, 0.97])))
        slab_oh.calc = calc
        BFGS(slab_oh, logfile=None).run(fmax=0.08, steps=100)
        e_oh = slab_oh.get_potential_energy()
        # ΔG_OH* = E(slab+OH) - E(slab) - (E(H2O) - 0.5*E(H2)) + ZPE + TS
        dG_OH = (e_oh - e_clean) - (e_h2o - 0.5 * e_h2) + ZPE_CORRECTIONS['OH'] + TS_CORRECTIONS['OH']
        result['dG_OH_eV'] = float(dG_OH)

        # 3. O* adsorption
        slab_o = structure.copy()
        o_pos = ads_base + np.array([0.0, 0.0, 1.7])
        slab_o.append(Atom('O', position=o_pos))
        slab_o.calc = calc
        BFGS(slab_o, logfile=None).run(fmax=0.08, steps=100)
        e_o = slab_o.get_potential_energy()
        # ΔG_O* = E(slab+O) - E(slab) - (E(H2O) - E(H2)) + ZPE + TS
        dG_O = (e_o - e_clean) - (e_h2o - e_h2) + ZPE_CORRECTIONS['O'] + TS_CORRECTIONS['O']
        result['dG_O_eV'] = float(dG_O)

        # 4. OOH* adsorption
        slab_ooh = structure.copy()
        o1_pos = ads_base + np.array([0.0, 0.0, 1.9])
        o2_pos = o1_pos + np.array([1.2, 0.0, 0.6])
        h_pos = o2_pos + np.array([0.0, 0.0, 0.97])
        slab_ooh.append(Atom('O', position=o1_pos))
        slab_ooh.append(Atom('O', position=o2_pos))
        slab_ooh.append(Atom('H', position=h_pos))
        slab_ooh.calc = calc
        BFGS(slab_ooh, logfile=None).run(fmax=0.08, steps=100)
        e_ooh = slab_ooh.get_potential_energy()
        # ΔG_OOH* = E(slab+OOH) - E(slab) - (2*E(H2O) - 1.5*E(H2)) + ZPE + TS
        dG_OOH = (e_ooh - e_clean) - (2 * e_h2o - 1.5 * e_h2) + ZPE_CORRECTIONS['OOH'] + TS_CORRECTIONS['OOH']
        result['dG_OOH_eV'] = float(dG_OOH)

        # 5. ORR overpotential (4e⁻ CHE method)
        eta, rds = orr_overpotential(dG_OH, dG_O, dG_OOH)
        result['orr_overpotential_V'] = float(eta)
        result['rate_determining_step'] = rds

        # 6. Fenton susceptibility
        elements = _extract_elements(genome)
        fenton_score = sum(FENTON_RISK.get(e, 0) for e in elements)
        result['fenton_stability'] = max(0, 10 - fenton_score)

        # 7. Cost
        result['cost_penalty'] = abundance_cost_penalty(elements)

        # 8. Stability estimate (binding strength of active metal)
        result['binding_strength'] = float(abs(dG_OH) + abs(dG_O))

        result['valid'] = True

    except Exception as e:
        result['error'] = str(e)[:200]

    return result


def _extract_elements(genome: tuple) -> List[str]:
    """Extract metallic elements from genome for cost/Fenton scoring."""
    mat_class = genome[0]
    elements = []
    if mat_class == 'MoltenMetal':
        elements.append(genome[1])
        if genome[2] != 'None': elements.append(genome[2])
    elif mat_class == 'SolidCatalyst':
        elements.append(genome[1])
        for d in genome[5]: elements.append(d)
    elif mat_class == 'SAC':
        elements.append(genome[1])
    elif mat_class == 'DAC':
        elements.extend([genome[1], genome[2]])
    elif mat_class in ('MOF', 'COF'):
        if genome[1] != 'None': elements.append(genome[1])
    elif mat_class == 'Perovskite':
        elements.extend([genome[1], genome[2]])
        if genome[3] != 'None': elements.append(genome[3])
    elif mat_class == 'MetalHydride':
        elements.append(genome[1])
        if genome[3] != 'None': elements.append(genome[3])
    elif mat_class == 'MAXPhase':
        elements.extend([genome[1], genome[2]])
        if genome[5] != 'None': elements.append(genome[5])
    elif mat_class == 'HEA':
        elements.extend(list(genome[1]))
    return [e for e in elements if e != 'None']


# ═══════════════════════════════════════════════════════════════════════════════
# MULTI-GPU WORKER
# ═══════════════════════════════════════════════════════════════════════════════

def orr_worker(worker_id: int, gpu_id: int, task_queue: mp.Queue,
               result_queue: mp.Queue):
    """Worker process: loads MACE on assigned GPU, evaluates ORR candidates."""
    try:
        from mace.calculators import mace_mp
        calc = mace_mp(model="medium", device=f"cuda:{gpu_id}")

        e_h2o = compute_water_ref(calc)
        e_h2 = compute_h2_ref(calc)

        while True:
            item = task_queue.get()
            if item is None:
                break
            idx, genome = item
            try:
                result = evaluate_orr_candidate(genome, calc, e_h2o, e_h2)
                result['worker_id'] = worker_id
                result['gpu_id'] = gpu_id
                result_queue.put((idx, result))
            except Exception as e:
                result_queue.put((idx, {
                    'genome': str(genome),
                    'material_class': genome[0],
                    'valid': False,
                    'error': str(e)[:200],
                }))
    except Exception as e:
        logger.error(f"ORR Worker {worker_id} failed: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN ORR SCREENING ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

def run_orr_screening(genomes: List[tuple], db_filename: str = "fc_mace_screening.csv",
                      workers_per_gpu: int = 2) -> 'pd.DataFrame':
    """
    Run parallel MACE ORR screening on a list of catalyst genomes.

    Same interface as mace_screener.run_screening but evaluates
    OH*, O*, OOH* for fuel cell cathode applications.
    """
    import pandas as pd

    mp.set_start_method('spawn', force=True)

    print_banner("MACE-MP-0 ORR CATHODE CATALYST SCREENING")
    logger.info(f"ORR screening {len(genomes)} candidates...")

    device_count = torch.cuda.device_count()
    num_workers = device_count * workers_per_gpu
    logger.info(f"Using {device_count} GPU(s), {num_workers} parallel workers")

    task_queue = mp.Queue()
    result_queue = mp.Queue()

    for idx, genome in enumerate(genomes):
        task_queue.put((idx, genome))
    for _ in range(num_workers):
        task_queue.put(None)

    workers = []
    for w_id in range(num_workers):
        gpu_id = w_id % device_count
        p = mp.Process(target=orr_worker, args=(w_id, gpu_id, task_queue, result_queue))
        p.start()
        workers.append(p)

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
                f"({rate:.1f} cand/sec, {n_valid} valid, {elapsed:.0f}s)"
            )

    for p in workers:
        p.join(timeout=30)

    df = pd.DataFrame(results)
    path = save_screening_db(df, db_filename, subdir="fuel_cell")
    logger.info(f"ORR screening complete. {len(df)} results saved to {path}")

    valid_df = df[df['valid'] == True]
    if len(valid_df) > 0:
        logger.info(f"  Valid: {len(valid_df)}/{len(df)}")
        logger.info(f"  Best overpotential: {valid_df['orr_overpotential_V'].min():.4f} V")
        logger.info(f"  ΔG_OH range: [{valid_df['dG_OH_eV'].min():.3f}, {valid_df['dG_OH_eV'].max():.3f}] eV")

    return df


if __name__ == '__main__':
    from pipeline.catalyst_spaces import generate_population
    pop = generate_population(20)
    df = run_orr_screening(pop, db_filename="test_orr_screening.csv", workers_per_gpu=2)
    print(f"\nORR screening complete: {len(df)} candidates")
