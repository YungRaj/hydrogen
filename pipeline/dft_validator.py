#!/usr/bin/env python3
"""
Quantum ESPRESSO DFT Validation Pipeline for Methane Pyrolysis Catalysts.

For champion catalysts from the MACE screening / genetic algorithm:
  1. Bulk structure optimization (vc-relax)
  2. Slab generation with correct Miller indices
  3. Adsorption energy calculations (H*, CH₃*, C*)
  4. NEB transition state for CH₄ → CH₃* + H*
  5. Electronic structure analysis (PDOS, d-band center)
  
Generates QE input files and submits calculations via pw.x.
"""

import os
import sys
import subprocess
import shutil
import json
import re
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.utils import (
    BASE_DIR, DFT_DIR, QE_PSEUDO_DIR, Ry_to_eV,
    setup_logger, print_banner, save_json,
)

logger = setup_logger('dft_validator', 'dft/dft_validation.log')

PW_X = os.environ.get("PW_X", shutil.which("pw.x") or "pw.x")


# ═══════════════════════════════════════════════════════════════════════════════
# PSEUDOPOTENTIAL MAPPING
# ═══════════════════════════════════════════════════════════════════════════════

PSEUDO_MAP = {
    'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
    'C': 'C.pbe-n-rrkjus_psl.1.0.0.UPF',
    'N': 'N.pbe-n-rrkjus_psl.1.0.0.UPF',
    'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
    'B': 'B.pbe-n-rrkjus_psl.1.0.0.UPF',
    'S': 'S.pbe-n-rrkjus_psl.1.0.0.UPF',
    'P': 'P.pbe-n-rrkjus_psl.1.0.0.UPF',
    'F': 'F.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Na': 'Na.pbe-spnl-rrkjus_psl.1.0.0.UPF',
    'Mg': 'Mg.pbe-spnl-rrkjus_psl.1.0.0.UPF',
    'Al': 'Al.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Ti': 'Ti.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'V': 'V.pbe-spnl-rrkjus_psl.1.0.0.UPF',
    'Cr': 'Cr.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Mn': 'Mn.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Co': 'Co.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Ni': 'Ni.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Cu': 'Cu.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'Zn': 'Zn.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'Ga': 'Ga.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'Ge': 'Ge.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'Mo': 'Mo.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Ru': 'Ru.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Rh': 'Rh.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Pd': 'Pd.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Ag': 'Ag.pbe-n-rrkjus_psl.1.0.0.UPF',
    'In': 'In.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'Sn': 'Sn.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'Sb': 'Sb.pbe-n-rrkjus_psl.1.0.0.UPF',
    'W': 'W.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Pt': 'Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF',
    'Au': 'Au.pbe-n-rrkjus_psl.1.0.0.UPF',
    'Pb': 'Pb.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'Bi': 'Bi.pbe-dn-rrkjus_psl.1.0.0.UPF',
    'La': 'La.pbe-spfn-rrkjus_psl.1.0.0.UPF',
    'Ce': 'Ce.pbe-spdn-rrkjus_psl.1.0.0.UPF',
    'Zr': 'Zr.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Y': 'Y.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Nb': 'Nb.pbe-spn-rrkjus_psl.1.0.0.UPF',
    'Te': 'Te.pbe-dn-rrkjus_psl.1.0.0.UPF',
}

# Atomic masses for QE input
ATOMIC_MASSES_QE = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'B': 10.81,
    'S': 32.06, 'P': 30.974, 'F': 18.998, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.086, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
    'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546,
    'Zn': 65.38, 'Ga': 69.723, 'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.906,
    'Pd': 106.42, 'Ag': 107.868, 'In': 114.818, 'Sn': 118.710,
    'Sb': 121.760, 'W': 183.84, 'Pt': 195.084, 'Au': 196.967,
    'Pb': 207.2, 'Bi': 208.980, 'La': 138.905, 'Ce': 140.116,
    'Zr': 91.224, 'Y': 88.906, 'Nb': 92.906, 'Te': 127.60,
    'Ge': 72.63,
}


# ═══════════════════════════════════════════════════════════════════════════════
# QE INPUT GENERATORS
# ═══════════════════════════════════════════════════════════════════════════════

def generate_bulk_scf_input(elements: List[str], positions_frac: List[Tuple],
                             celldm_bohr: float, ibrav: int = 2,
                             calc_name: str = "bulk_scf",
                             ecutwfc: float = 50.0,
                             ecutrho: float = 400.0,
                             kpoints: Tuple = (6, 6, 6)) -> str:
    """Generate a QE scf input for a bulk crystal."""
    nat = len(positions_frac)
    ntyp = len(set(elements))

    from pipeline.qe_workflows import verify_sssp, SSSP_DIR
    sssp = verify_sssp(elements)
    if not sssp['valid']:
        raise RuntimeError(f"SSSP verification failed: {sssp['errors']}")
    ecutwfc = max(float(ecutwfc), sssp['ecutwfc_Ry'])
    ecutrho = max(float(ecutrho), sssp['ecutrho_Ry'])
    species_block = ""
    unique_elements = sorted(set(elements))
    for elem in unique_elements:
        pseudo = sssp['records'][elem]['filename']
        mass = ATOMIC_MASSES_QE.get(elem, 50.0)
        species_block += f"  {elem}  {mass:.3f}  {pseudo}\n"

    atoms_block = ""
    for elem, pos in zip(elements, positions_frac):
        atoms_block += f"  {elem}  {pos[0]:.10f}  {pos[1]:.10f}  {pos[2]:.10f}\n"

    input_text = f"""&CONTROL
  calculation = 'scf'
  prefix = '{calc_name}'
  outdir = './tmp'
  pseudo_dir = '{SSSP_DIR}'
  tprnfor = .true.
  tstress = .true.
  verbosity = 'high'
/
&SYSTEM
  ibrav = {ibrav}
  celldm(1) = {celldm_bohr}
  nat = {nat}
  ntyp = {ntyp}
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.02
/
&ELECTRONS
  mixing_beta = 0.3
  conv_thr = 1.0d-8
  electron_maxstep = 200
/
ATOMIC_SPECIES
{species_block}
ATOMIC_POSITIONS {{crystal}}
{atoms_block}
K_POINTS {{automatic}}
  {kpoints[0]} {kpoints[1]} {kpoints[2]}  1 1 1
"""
    return input_text


def generate_slab_scf_input(elements: List[str], positions_ang: List[Tuple],
                             cell_params: List[List[float]],
                             calc_name: str = "slab_scf",
                             ecutwfc: float = 50.0,
                             kpoints: Tuple = (4, 4, 1)) -> str:
    """Generate QE scf input for a surface slab (ibrav=0)."""
    nat = len(positions_ang)
    ntyp = len(set(elements))

    from pipeline.qe_workflows import verify_sssp, SSSP_DIR
    sssp = verify_sssp(elements)
    if not sssp['valid']:
        raise RuntimeError(f"SSSP verification failed: {sssp['errors']}")
    ecutwfc = max(float(ecutwfc), sssp['ecutwfc_Ry'])
    ecutrho = sssp['ecutrho_Ry']
    species_block = ""
    unique_elements = sorted(set(elements))
    for elem in unique_elements:
        pseudo = sssp['records'][elem]['filename']
        mass = ATOMIC_MASSES_QE.get(elem, 50.0)
        species_block += f"  {elem}  {mass:.3f}  {pseudo}\n"

    atoms_block = ""
    for elem, pos in zip(elements, positions_ang):
        atoms_block += f"  {elem}  {pos[0]:.10f}  {pos[1]:.10f}  {pos[2]:.10f}\n"

    cell_block = ""
    for row in cell_params:
        cell_block += f"  {row[0]:.10f}  {row[1]:.10f}  {row[2]:.10f}\n"

    input_text = f"""&CONTROL
  calculation = 'relax'
  prefix = '{calc_name}'
  outdir = './tmp'
  pseudo_dir = '{SSSP_DIR}'
  tprnfor = .true.
  forc_conv_thr = 1.0d-3
/
&SYSTEM
  ibrav = 0
  nat = {nat}
  ntyp = {ntyp}
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.02
  nspin = 2
/
&ELECTRONS
  mixing_beta = 0.3
  conv_thr = 1.0d-7
  electron_maxstep = 200
/
&IONS
  ion_dynamics = 'bfgs'
/
CELL_PARAMETERS {{angstrom}}
{cell_block}
ATOMIC_SPECIES
{species_block}
ATOMIC_POSITIONS {{angstrom}}
{atoms_block}
K_POINTS {{automatic}}
  {kpoints[0]} {kpoints[1]} {kpoints[2]}  1 1 1
"""
    return input_text


# ═══════════════════════════════════════════════════════════════════════════════
# DFT OUTPUT PARSING
# ═══════════════════════════════════════════════════════════════════════════════

def parse_total_energy(output_file: str) -> Optional[float]:
    """Parse total energy from QE output file (Ry)."""
    if not os.path.exists(output_file):
        return None
    with open(output_file, 'r') as f:
        content = f.read()
    matches = re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", content)
    if matches:
        return float(matches[-1])
    # Fallback: non-converged energy
    matches = re.findall(r"total energy\s+=\s+(-?\d+\.\d+)\s+Ry", content)
    return float(matches[-1]) if matches else None


def parse_forces(output_file: str) -> Optional[float]:
    """Parse maximum force from QE output (Ry/bohr)."""
    if not os.path.exists(output_file):
        return None
    with open(output_file, 'r') as f:
        content = f.read()
    matches = re.findall(r"Total force\s+=\s+(\d+\.\d+)", content)
    return float(matches[-1]) if matches else None


def parse_convergence(output_file: str) -> bool:
    """Require clean QE termination and electronic convergence."""
    if not os.path.exists(output_file):
        return False
    with open(output_file, 'r') as f:
        content = f.read()
    lower = content.lower()
    fatal = ('error in routine', 'convergence not achieved', 'stopping ...')
    return ('job done' in lower and 'convergence has been achieved' in lower and
            not any(token in lower for token in fatal))


# ═══════════════════════════════════════════════════════════════════════════════
# DFT VALIDATION WORKFLOW
# ═══════════════════════════════════════════════════════════════════════════════

def validate_catalyst(catalyst_name: str, genome: tuple,
                       run_dft: bool = True) -> Dict:
    """
    Full DFT validation workflow for a champion catalyst.
    
    Steps:
      1. Generate bulk structure input
      2. Run SCF (if run_dft=True)
      3. Parse results
      
    Returns dict with DFT validation results.
    """
    print_banner(f"DFT VALIDATION: {catalyst_name}")

    calc_dir = DFT_DIR / catalyst_name
    calc_dir.mkdir(parents=True, exist_ok=True)

    mat_class = genome[0]
    result = {
        'catalyst_name': catalyst_name,
        'material_class': mat_class,
        'genome': str(genome),
    }

    # Generate input based on material class
    if mat_class in ('MoltenMetal', 'SolidCatalyst'):
        input_text = _generate_alloy_input(genome, catalyst_name)
    elif mat_class in ('SAC', 'DAC', 'MOF', 'COF'):
        input_text = _generate_cluster_input(genome, catalyst_name)
    else:
        # All other classes (Perovskite, HEA, MAXPhase, MetalHydride,
        # Spinel, MXene, SAA, MetalFreeCarbon): build structure from genome
        # and use the alloy slab DFT pathway
        from pipeline.surface_screener import generate_structure
        try:
            structure, _, _ = generate_structure(genome)
            elements = [a.symbol for a in structure if a.symbol != 'X']
            positions = [tuple(a.position) for a in structure if a.symbol != 'X']
            cell_size = 15.0
            cell_params = [[cell_size, 0, 0], [0, cell_size, 0], [0, 0, cell_size]]
            from pipeline.dft_fuel_cell import generate_slab_scf_input
            input_text = generate_slab_scf_input(
                elements, positions, cell_params,
                calc_name=catalyst_name, ecutwfc=40.0, kpoints=(1, 1, 1)
            )
        except Exception as e:
            logger.warning(f"Structure generation failed for {mat_class}: {e}")
            result['error'] = str(e)[:200]
            return result

    # Write input file
    input_file = calc_dir / f"{catalyst_name}.in"
    output_file = calc_dir / f"{catalyst_name}.out"
    with open(input_file, 'w') as f:
        f.write(input_text)
    logger.info(f"  Written QE input: {input_file}")

    # Run DFT
    if run_dft:
        logger.info(f"  Running pw.x for {catalyst_name}...")
        try:
            proc = subprocess.run(
                f"{PW_X} < {input_file} > {output_file}",
                shell=True, cwd=str(calc_dir),
                timeout=3600, capture_output=True, text=True,
            )
            result['returncode'] = proc.returncode
        except subprocess.TimeoutExpired:
            logger.warning(f"  DFT calculation timed out for {catalyst_name}")
            result['returncode'] = -1
            result['error'] = 'timeout'

    # Parse results
    if output_file.exists():
        energy = parse_total_energy(str(output_file))
        converged = parse_convergence(str(output_file))
        max_force = parse_forces(str(output_file))

        result['dft_energy_Ry'] = energy
        result['dft_energy_eV'] = energy * Ry_to_eV if energy else None
        result['converged'] = converged
        result['max_force_Ry_bohr'] = max_force

        if converged and energy:
            logger.info(f"  ✓ Converged: E = {energy:.6f} Ry ({energy * Ry_to_eV:.4f} eV)")
        else:
            logger.warning(f"  ✗ Did not converge or no energy found")
    else:
        logger.warning(f"  No output file found: {output_file}")

    save_json(result, f"{catalyst_name}_dft.json", subdir="dft")
    return result


def _generate_alloy_input(genome: tuple, name: str) -> str:
    """Generate QE input for an alloy bulk structure."""
    mat_class = genome[0]
    if mat_class == 'MoltenMetal':
        host = genome[1]
        promoter = genome[2]
        elements = [host, host, host]
        if promoter != 'None':
            elements[0] = promoter
        # Simple L1₂ structure
        positions = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5)]
    else:
        metal = genome[1]
        elements = [metal, metal, metal, metal]
        positions = [
            (0.0, 0.0, 0.0), (0.5, 0.5, 0.0),
            (0.5, 0.0, 0.5), (0.0, 0.5, 0.5),
        ]
        # Apply dopant substitution
        if genome[5]:
            elements[0] = genome[5][0]

    # Estimate lattice parameter (Bohr) from element
    lattice_params_bohr = {
        'Ni': 6.65, 'Fe': 5.42, 'Co': 6.68, 'Cu': 6.82, 'Al': 7.65,
        'Sn': 12.30, 'Bi': 9.54, 'In': 8.73, 'Ga': 5.72, 'Pb': 9.35,
        'Sb': 8.63, 'Mo': 5.95, 'W': 5.98, 'Mn': 16.85, 'Pt': 7.42,
        'Pd': 7.35, 'Ag': 7.73, 'Au': 7.71, 'Ti': 5.58, 'V': 5.72,
        'Cr': 5.45, 'Zn': 5.03, 'Zr': 6.10, 'Y': 6.89, 'La': 7.10,
    }
    primary = elements[0] if elements[0] in lattice_params_bohr else genome[1]
    celldm = lattice_params_bohr.get(primary, 7.0)

    return generate_bulk_scf_input(
        elements, positions, celldm, ibrav=2,
        calc_name=name, ecutwfc=50.0, kpoints=(6, 6, 6)
    )


def _generate_cluster_input(genome: tuple, name: str) -> str:
    """Generate QE input for a molecular cluster (SAC/MOF)."""
    from pipeline.surface_screener import generate_porphyrin_cluster

    mat_class = genome[0]
    if mat_class == 'SAC':
        metal, coord = genome[1], genome[2]
    elif mat_class == 'DAC':
        metal, coord = genome[1], genome[3]
    elif mat_class in ('MOF', 'COF'):
        metal, coord = genome[1], genome[3]
    else:
        metal, coord = 'Ni', 'N4'

    cluster = generate_porphyrin_cluster(metal, coord)

    elements = [atom.symbol for atom in cluster if atom.symbol != 'X']
    positions = [tuple(atom.position) for atom in cluster if atom.symbol != 'X']

    cell_size = 15.0
    cell_params = [
        [cell_size, 0.0, 0.0],
        [0.0, cell_size, 0.0],
        [0.0, 0.0, cell_size],
    ]

    return generate_slab_scf_input(
        elements, positions, cell_params,
        calc_name=name, ecutwfc=40.0, kpoints=(1, 1, 1)
    )


if __name__ == '__main__':
    # Test: generate input for a Ni-Bi molten metal catalyst
    genome = ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000)
    result = validate_catalyst("NiBi_10pct", genome, run_dft=False)
    print(json.dumps(result, indent=2, default=str))
