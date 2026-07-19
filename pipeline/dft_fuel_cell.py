#!/usr/bin/env python3
"""
DFT Validation for Fuel Cell ORR Catalysts.

Computes adsorption free energies for ORR intermediates (O*, OH*, OOH*)
on champion cathode catalysts using the Computational Hydrogen Electrode (CHE).

Generates QE input files for each intermediate adsorbed on the active site,
then computes the theoretical ORR overpotential.
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
from typing import Dict, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.utils import (
    DFT_DIR, QE_PSEUDO_DIR, Ry_to_eV, E_ORR_eq,
    ZPE_H2, TS_H2, G_H2_correction,
    setup_logger, save_json, orr_overpotential,
)
from pipeline.dft_validator import (
    generate_slab_scf_input, parse_total_energy, parse_convergence, PW_X,
    PSEUDO_MAP, ATOMIC_MASSES_QE,
)
from pipeline.surface_screener import generate_porphyrin_cluster

logger = setup_logger('dft_fuel_cell', 'dft/dft_fuel_cell.log')


# ═══════════════════════════════════════════════════════════════════════════════
# FREE ENERGY CORRECTIONS (ZPE + entropy at 298 K, 1 bar)
# ═══════════════════════════════════════════════════════════════════════════════
# From Nørskov et al., J. Phys. Chem. B 108, 17886 (2004)

ZPE_CORRECTIONS_eV = {
    'OH*': 0.35,
    'O*': 0.05,
    'OOH*': 0.40,
    'H2O_g': 0.56,
    'H2_g': 0.27,
    'O2_g': 0.10,
}

TS_CORRECTIONS_eV = {
    'OH*': 0.07,
    'O*': 0.00,
    'OOH*': 0.10,
    'H2O_g': 0.67,
    'H2_g': 0.40,
    'O2_g': 0.63,
}


def compute_dG_adsorbate(E_slab_ads: float, E_slab_clean: float,
                          E_ref: float, ads_type: str) -> float:
    """
    Compute adsorption free energy using CHE method.
    
    ΔG = ΔE + ΔZPE − TΔS
    
    Args:
        E_slab_ads: DFT energy of slab + adsorbate (eV)
        E_slab_clean: DFT energy of clean slab (eV)
        E_ref: Reference molecule energy (eV)
        ads_type: 'OH*', 'O*', or 'OOH*'
    
    Returns: Adsorption free energy (eV)
    """
    dE = E_slab_ads - E_slab_clean - E_ref
    dZPE = ZPE_CORRECTIONS_eV.get(ads_type, 0.0)
    dTS = TS_CORRECTIONS_eV.get(ads_type, 0.0)
    return dE + dZPE - dTS


def validate_orr_catalyst(catalyst_name: str, genome: tuple,
                           run_dft: bool = True) -> Dict:
    """
    Full ORR catalyst validation workflow.
    
    1. Clean slab/cluster SCF
    2. OH* adsorption SCF
    3. O* adsorption SCF  
    4. OOH* adsorption SCF
    5. Compute free energy diagram
    6. Compute overpotential
    """
    from pipeline.dft_validator import validate_catalyst

    logger.info(f"ORR DFT validation: {catalyst_name}")
    calc_dir = DFT_DIR / f"fc_{catalyst_name}"
    calc_dir.mkdir(parents=True, exist_ok=True)

    result = {
        'catalyst_name': catalyst_name,
        'genome': str(genome),
        'material_class': genome[0],
    }

    mat_class = genome[0]

    # Generate the base cluster/slab structure
    if mat_class in ('SAC', 'DAC'):
        metal = genome[1]
        coord = genome[2] if mat_class == 'SAC' else genome[3]
        cluster = generate_porphyrin_cluster(metal, coord)
    elif mat_class in ('MOF', 'COF'):
        metal = genome[1]
        coord = genome[3]
        cluster = generate_porphyrin_cluster(metal, coord)
    else:
        # Alloy catalyst — build from genome
        from pipeline.surface_screener import generate_structure
        cluster, _, _ = generate_structure(genome)

    base_elements = [a.symbol for a in cluster if a.symbol != 'X']
    base_positions = [tuple(a.position) for a in cluster if a.symbol != 'X']
    cell = 15.0
    cell_params = [[cell, 0, 0], [0, cell, 0], [0, 0, cell]]

    # ── 1. Clean slab SCF ───────────────────────────────────────────────────
    clean_input = generate_slab_scf_input(
        base_elements, base_positions, cell_params,
        calc_name=f"fc_{catalyst_name}_clean", ecutwfc=40.0, kpoints=(1, 1, 1)
    )
    clean_in = calc_dir / f"{catalyst_name}_clean.in"
    clean_out = calc_dir / f"{catalyst_name}_clean.out"
    with open(clean_in, 'w') as f:
        f.write(clean_input)

    if run_dft:
        _run_pw(clean_in, clean_out, calc_dir)

    E_clean = parse_total_energy(str(clean_out))
    result['E_clean_Ry'] = E_clean
    if E_clean:
        E_clean_eV = E_clean * Ry_to_eV
        result['E_clean_eV'] = E_clean_eV
    else:
        E_clean_eV = None

    # ── 2-4. Adsorbate calculations ─────────────────────────────────────────
    ads_configs = {
        'OH': {'elements': ['O', 'H'], 'offsets': [np.array([0, 0, 1.8]), np.array([0, 0, 2.8])]},
        'O': {'elements': ['O'], 'offsets': [np.array([0, 0, 1.6])]},
        'OOH': {'elements': ['O', 'O', 'H'], 'offsets': [np.array([0, 0, 1.8]), np.array([1.2, 0, 2.4]), np.array([1.2, 0, 3.3])]},
    }

    # Active site position (first metal atom or center)
    active_pos = np.array(base_positions[0]) if base_positions else np.array([cell/2, cell/2, cell/2])

    for ads_name, ads_info in ads_configs.items():
        ads_elements = list(base_elements)
        ads_positions = list(base_positions)

        for elem, offset in zip(ads_info['elements'], ads_info['offsets']):
            ads_elements.append(elem)
            ads_positions.append(tuple(active_pos + offset))

        ads_input = generate_slab_scf_input(
            ads_elements, ads_positions, cell_params,
            calc_name=f"fc_{catalyst_name}_{ads_name}", ecutwfc=40.0, kpoints=(1, 1, 1)
        )

        ads_in = calc_dir / f"{catalyst_name}_{ads_name}.in"
        ads_out = calc_dir / f"{catalyst_name}_{ads_name}.out"
        with open(ads_in, 'w') as f:
            f.write(ads_input)

        if run_dft:
            _run_pw(ads_in, ads_out, calc_dir)

        E_ads = parse_total_energy(str(ads_out))
        result[f'E_{ads_name}_Ry'] = E_ads
        if E_ads:
            result[f'E_{ads_name}_eV'] = E_ads * Ry_to_eV

    # ── 4b. Gas-phase references (H₂O and H₂) ─────────────────────────────────
    h2_elements = ['H', 'H']
    h2_positions = [(cell/2, cell/2, cell/2 - 0.37), (cell/2, cell/2, cell/2 + 0.37)]
    h2_input = generate_slab_scf_input(
        h2_elements, h2_positions, cell_params,
        calc_name=f"fc_{catalyst_name}_h2", ecutwfc=40.0, kpoints=(1, 1, 1)
    )
    h2_in = calc_dir / f"{catalyst_name}_h2.in"
    h2_out = calc_dir / f"{catalyst_name}_h2.out"
    with open(h2_in, 'w') as f:
        f.write(h2_input)

    h2o_elements = ['O', 'H', 'H']
    h2o_positions = [
        (cell/2, cell/2, cell/2),
        (cell/2 + 0.757, cell/2 + 0.586, cell/2),
        (cell/2 - 0.757, cell/2 + 0.586, cell/2)
    ]
    h2o_input = generate_slab_scf_input(
        h2o_elements, h2o_positions, cell_params,
        calc_name=f"fc_{catalyst_name}_h2o", ecutwfc=40.0, kpoints=(1, 1, 1)
    )
    h2o_in = calc_dir / f"{catalyst_name}_h2o.in"
    h2o_out = calc_dir / f"{catalyst_name}_h2o.out"
    with open(h2o_in, 'w') as f:
        f.write(h2o_input)

    if run_dft:
        logger.info(f"  Running gas-phase reference H₂...")
        _run_pw(h2_in, h2_out, calc_dir)
        logger.info(f"  Running gas-phase reference H₂O...")
        _run_pw(h2o_in, h2o_out, calc_dir)

    E_h2_Ry = parse_total_energy(str(h2_out))
    E_h2o_Ry = parse_total_energy(str(h2o_out))

    if E_h2_Ry and E_h2o_Ry:
        E_h2_eV = E_h2_Ry * Ry_to_eV
        E_h2o_eV = E_h2o_Ry * Ry_to_eV
        result['E_H2_eV'] = E_h2_eV
        result['E_H2O_eV'] = E_h2o_eV
    else:
        # Fallbacks for non-run or missing outputs
        E_h2_eV = -31.8
        E_h2o_eV = -432.6
        result['E_H2_eV_fallback'] = E_h2_eV
        result['E_H2O_eV_fallback'] = E_h2o_eV

    # ── 5. Free energy diagram ──────────────────────────────────────────────
    if E_clean_eV and all(result.get(f'E_{a}_eV') for a in ['OH', 'O', 'OOH']):
        # Reference energies: H₂O(g) and H₂(g)
        # Using standard CHE: dG_OH = E(slab+OH) - E(slab) - (E_H2O - 0.5*E_H2) + corrections
        E_OH = result['E_OH_eV']
        E_O = result['E_O_eV']
        E_OOH = result['E_OOH_eV']

        # Adsorption free energies (relative to H₂O and H₂ references)
        dG_OH = (E_OH - E_clean_eV) - (E_h2o_eV - 0.5 * E_h2_eV) + ZPE_CORRECTIONS_eV['OH*'] - TS_CORRECTIONS_eV['OH*']
        dG_O = (E_O - E_clean_eV) - (E_h2o_eV - E_h2_eV) + ZPE_CORRECTIONS_eV['O*'] - TS_CORRECTIONS_eV['O*']
        dG_OOH = (E_OOH - E_clean_eV) - (2 * E_h2o_eV - 1.5 * E_h2_eV) + ZPE_CORRECTIONS_eV['OOH*'] - TS_CORRECTIONS_eV['OOH*']

        result['dG_OH_eV'] = float(dG_OH)
        result['dG_O_eV'] = float(dG_O)
        result['dG_OOH_eV'] = float(dG_OOH)

        eta, rds = orr_overpotential(dG_OH, dG_O, dG_OOH)
        result['orr_overpotential_V'] = float(eta)
        result['rate_determining_step'] = rds
        result['limiting_potential_V'] = float(E_ORR_eq - eta)

        logger.info(f"  ORR overpotential: η = {eta:.4f} V (RDS: {rds})")
    else:
        logger.warning("  Cannot compute overpotential: missing DFT energies")

    required_outputs = [clean_out, h2_out, h2o_out] + [
        calc_dir / f"{catalyst_name}_{name}.out" for name in ('OH', 'O', 'OOH')]
    result['converged'] = bool(run_dft and all(
        parse_convergence(str(path)) for path in required_outputs))
    result['evidence_level'] = 'converged_dft' if result['converged'] else 'incomplete'
    if not result['converged']:
        # Fallback molecular energies may aid input-generation smoke tests but
        # can never establish an ORR champion.
        result.pop('orr_overpotential_V', None)
        result.pop('limiting_potential_V', None)

    save_json(result, f"fc_{catalyst_name}_orr.json", subdir="dft")
    return result


def _run_pw(input_file: Path, output_file: Path, cwd: Path):
    """Execute pw.x calculation."""
    import subprocess
    try:
        subprocess.run(
            f"{PW_X} < {input_file} > {output_file}",
            shell=True, cwd=str(cwd), timeout=3600,
            capture_output=True, text=True,
        )
    except subprocess.TimeoutExpired:
        logger.warning(f"  DFT timed out: {input_file}")
    except Exception as e:
        logger.error(f"  DFT error: {e}")


if __name__ == '__main__':
    # Test: generate ORR inputs for a Fe-N4 SAC
    genome = ('SAC', 'Fe', 'N4', 'N-graphene')
    result = validate_orr_catalyst("FeN4_SAC", genome, run_dft=False)
    print(json.dumps(result, indent=2, default=str))
