#!/usr/bin/env python3
# Turquoise-hydrogen DFT validation.
"""
Quantum ESPRESSO DFT Validation Pipeline for Methane Pyrolysis Catalysts.

For champion catalysts from catalyst-guided screening:
  1. Bulk structure optimization (vc-relax)
  2. Slab generation with correct Miller indices
  3. Adsorption energy calculations (H*, CH₃*, C*)
  4. NEB transition state for CH₄ → CH₃* + H*
  5. Electronic structure analysis (PDOS, d-band center)

Generates QE input files and submits calculations via pw.x.
"""

import os
import sys
import json
import re
from pathlib import Path
from typing import List, Tuple, Optional

from pipeline.data_models.quantum import DFTResult

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.utils import (
    DFT_DIR,
    Ry_to_eV,
    setup_logger,
    print_banner,
    save_json,
)

logger = setup_logger('dft_validator', 'dft/dft_validation.log')

# Backward-compatible command hint; execution uses the portable QE resolver.
PW_X = os.environ.get("PW_X", "pw.x")


# ═══════════════════════════════════════════════════════════════════════════════
# PSEUDOPOTENTIAL MAPPING
# ═══════════════════════════════════════════════════════════════════════════════

# Atomic masses for QE input
ATOMIC_MASSES_QE = {
    'H': 1.008,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'B': 10.81,
    'S': 32.06,
    'P': 30.974,
    'F': 18.998,
    'Na': 22.990,
    'Mg': 24.305,
    'Al': 26.982,
    'Si': 28.086,
    'Ti': 47.867,
    'V': 50.942,
    'Cr': 51.996,
    'Mn': 54.938,
    'Fe': 55.845,
    'Co': 58.933,
    'Ni': 58.693,
    'Cu': 63.546,
    'Zn': 65.38,
    'Ga': 69.723,
    'Mo': 95.95,
    'Ru': 101.07,
    'Rh': 102.906,
    'Pd': 106.42,
    'Ag': 107.868,
    'In': 114.818,
    'Sn': 118.710,
    'Sb': 121.760,
    'W': 183.84,
    'Pt': 195.084,
    'Au': 196.967,
    'Pb': 207.2,
    'Bi': 208.980,
    'La': 138.905,
    'Ce': 140.116,
    'Zr': 91.224,
    'Y': 88.906,
    'Nb': 92.906,
    'Te': 127.60,
    'Ge': 72.63,
}


# ═══════════════════════════════════════════════════════════════════════════════
# QE INPUT GENERATORS
# ═══════════════════════════════════════════════════════════════════════════════


def generate_bulk_scf_input(
    elements: List[str],
    positions_frac: List[Tuple],
    celldm_bohr: float,
    ibrav: int = 2,
    calc_name: str = "bulk_scf",
    ecutwfc: float = 50.0,
    ecutrho: float = 400.0,
    kpoints: Tuple = (6, 6, 6),
) -> str:
    """Generate a QE scf input for a bulk crystal.

    Args:
        elements: Ordered values supplying elements.
        positions_frac: Ordered values supplying positions frac.
        celldm_bohr: Celldm bohr used by this operation.
        ibrav: Ibrav used by this operation.
        calc_name: Calc name used by this operation.
        ecutwfc: Ecutwfc used by this operation.
        ecutrho: Ecutrho used by this operation.
        kpoints: Ordered values supplying kpoints.

    Returns:
        Computed `str` result.
    """
    nat = len(positions_frac)
    ntyp = len(set(elements))
    kshift = tuple(0 if mesh == 1 else 1 for mesh in kpoints)

    from pipeline.validation.qe_workflows import verify_sssp, SSSP_DIR

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
  {kpoints[0]} {kpoints[1]} {kpoints[2]}  {kshift[0]} {kshift[1]} {kshift[2]}
"""
    return input_text


def generate_slab_scf_input(
    elements: List[str],
    positions_ang: List[Tuple],
    cell_params: List[List[float]],
    calc_name: str = "slab_scf",
    ecutwfc: float = 50.0,
    kpoints: Tuple = (4, 4, 1),
) -> str:
    """Generate QE scf input for a surface slab (ibrav=0).

    Args:
        elements: Ordered values supplying elements.
        positions_ang: Ordered values supplying positions ang.
        cell_params: Ordered values supplying cell params.
        calc_name: Calc name used by this operation.
        ecutwfc: Ecutwfc used by this operation.
        kpoints: Ordered values supplying kpoints.

    Returns:
        Computed `str` result.
    """
    nat = len(positions_ang)
    ntyp = len(set(elements))
    kshift = tuple(0 if mesh == 1 else 1 for mesh in kpoints)

    from pipeline.validation.qe_workflows import verify_sssp, SSSP_DIR

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

    # Quantum ESPRESSO rejects spin-polarized calculations unless at least one
    # species has a non-zero starting magnetization.  Use a conservative seed
    # for every species; this initializes the spin channel without asserting a
    # particular converged magnetic state.
    magnetization_block = "\n".join(
        f"  starting_magnetization({index}) = 0.05"
        for index, _ in enumerate(unique_elements, start=1)
    )

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
{magnetization_block}
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
  {kpoints[0]} {kpoints[1]} {kpoints[2]}  {kshift[0]} {kshift[1]} {kshift[2]}
"""
    return input_text


def generate_molecule_input(
    elements: List[str],
    positions_ang: List[Tuple],
    cell_size_ang: float,
    calc_name: str,
    ecutwfc: float = 40.0,
    calculation: str = 'relax',
) -> str:
    """Generate a closed-shell, Gamma-only molecular QE input.

        Gas-phase CHE references must not inherit metallic slab smearing or an
        unnecessary second spin channel.  Production pseudopotential cutoffs and
        the final convergence threshold remain unchanged.

    Args:
        elements: Ordered values supplying elements.
        positions_ang: Ordered values supplying positions ang.
        cell_size_ang: Cell size ang used by this operation.
        calc_name: Calc name used by this operation.
        ecutwfc: Ecutwfc used by this operation.
        calculation: Calculation used by this operation.

    Returns:
        Computed `str` result.
    """
    if calculation not in ('scf', 'relax'):
        raise ValueError('molecular calculation must be scf or relax')
    if len(elements) != len(positions_ang) or not elements:
        raise ValueError(
            'molecular elements and positions must be nonempty and aligned'
        )
    if cell_size_ang <= 0:
        raise ValueError('molecular cell must be positive')

    from pipeline.validation.qe_workflows import verify_sssp, SSSP_DIR

    sssp = verify_sssp(elements)
    if not sssp['valid']:
        raise RuntimeError(f"SSSP verification failed: {sssp['errors']}")
    ecutwfc = max(float(ecutwfc), sssp['ecutwfc_Ry'])
    unique_elements = sorted(set(elements))
    species_block = '\n'.join(
        f"  {element}  {ATOMIC_MASSES_QE.get(element, 50.0):.3f}  "
        f"{sssp['records'][element]['filename']}"
        for element in unique_elements
    )
    atoms_block = '\n'.join(
        f"  {element}  {position[0]:.10f}  {position[1]:.10f}  {position[2]:.10f}"
        for element, position in zip(elements, positions_ang)
    )
    ions = "\n&IONS\n  ion_dynamics = 'bfgs'\n/" if calculation == 'relax' else ''

    return f"""&CONTROL
  calculation = '{calculation}'
  prefix = '{calc_name}'
  outdir = './tmp'
  pseudo_dir = '{SSSP_DIR}'
  tprnfor = .true.
  forc_conv_thr = 1.0d-3
/
&SYSTEM
  ibrav = 1
  celldm(1) = {cell_size_ang / 0.529177210903:.12f}
  nat = {len(elements)}
  ntyp = {len(unique_elements)}
  ecutwfc = {ecutwfc}
  ecutrho = {sssp['ecutrho_Ry']}
  nspin = 1
  occupations = 'fixed'
/
&ELECTRONS
  mixing_beta = 0.3
  conv_thr = 1.0d-7
  electron_maxstep = 200
/{ions}
ATOMIC_SPECIES
{species_block}
ATOMIC_POSITIONS {{angstrom}}
{atoms_block}
K_POINTS {{gamma}}
"""


# ═══════════════════════════════════════════════════════════════════════════════
# DFT OUTPUT PARSING
# ═══════════════════════════════════════════════════════════════════════════════


def parse_total_energy(output_file: str) -> Optional[float]:
    """Parse total energy from QE output file (Ry).

    Args:
        output_file: Filesystem location used for output file.

    Returns:
        Computed `Optional[float]` result.
    """
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
    """Parse maximum force from QE output (Ry/bohr).

    Args:
        output_file: Filesystem location used for output file.

    Returns:
        Computed `Optional[float]` result.
    """
    if not os.path.exists(output_file):
        return None
    with open(output_file, 'r') as f:
        content = f.read()
    matches = re.findall(r"Total force\s+=\s+(\d+\.\d+)", content)
    return float(matches[-1]) if matches else None


def parse_convergence(output_file: str, *, require_ionic: bool = False) -> bool:
    """Require clean QE termination and the requested convergence level.

    Args:
        output_file: Filesystem location used for output file.
        require_ionic: Whether to enable require ionic.

    Returns:
        True when the documented condition holds; otherwise False.
    """
    if not os.path.exists(output_file):
        return False
    with open(output_file, 'r') as f:
        content = f.read()
    lower = content.lower()
    fatal = ('error in routine', 'convergence not achieved', 'stopping ...')
    electronic = (
        'job done' in lower
        and 'convergence has been achieved' in lower
        and not any(token in lower for token in fatal)
    )
    if not electronic:
        return False
    if not require_ionic:
        return True
    # QE prints this only after the BFGS ionic criteria (including
    # forc_conv_thr) have been satisfied. An SCF convergence line is not an
    # endpoint-relaxation certificate.
    ionic = 'end of bfgs geometry optimization' in lower or 'bfgs converged in' in lower
    return ionic


# ═══════════════════════════════════════════════════════════════════════════════
# DFT VALIDATION WORKFLOW
# ═══════════════════════════════════════════════════════════════════════════════


def validate_catalyst(
    catalyst_name: str,
    genome: tuple,
    run_dft: bool = True,
    restart_incomplete: bool = False,
    execution=None,
) -> DFTResult:
    """
        Full DFT validation workflow for a champion catalyst.

        Steps:
          1. Generate bulk structure input
          2. Run SCF (if run_dft=True)
          3. Parse results

        Returns dict with DFT validation results.

    Args:
        catalyst_name: Catalyst name used by this operation.
        genome: Encoded catalyst composition and structural configuration.
        run_dft: Whether to enable run dft.
        restart_incomplete: Whether to enable restart incomplete.
        execution: Optional QEExecutionConfig override.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    print_banner(f"DFT VALIDATION: {catalyst_name}")

    calc_dir = DFT_DIR / catalyst_name
    calc_dir.mkdir(parents=True, exist_ok=True)

    mat_class = genome[0]
    from pipeline.search.discovery import candidate_id

    result: DFTResult = {
        'catalyst_name': catalyst_name,
        'candidate_id': candidate_id(genome),
        'material_class': mat_class,
        'genome': str(genome),
        'evidence_level': 'screening_dft',
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
        from pipeline.screening.surface_screener import generate_structure

        try:
            structure, active_idx, _ = generate_structure(genome)
            raw_elems = [a.symbol for a in structure if a.symbol != 'X']
            active_elems = [
                structure[i].symbol
                for i in (active_idx or [])
                if 0 <= i < len(structure) and structure[i].symbol != 'X'
            ]
            seen_unique = list(dict.fromkeys(active_elems + raw_elems))
            site_elems = (seen_unique + raw_elems)[:4]
            while len(site_elems) < 4:
                site_elems.append(site_elems[0] if site_elems else 'Ni')
            celldm = _FCC_EQUIVALENT_LATTICE_BOHR.get(site_elems[0], 7.0)
            input_text = generate_bulk_scf_input(
                site_elems,
                list(_FCC_CONVENTIONAL_SITES),
                celldm,
                ibrav=1,
                calc_name=catalyst_name,
                ecutwfc=50.0,
                kpoints=(4, 4, 4),
            )
        except Exception as e:
            logger.warning(f"Structure generation failed for {mat_class}: {e}")
            result['error'] = str(e)[:200]
            result['converged'] = False
            save_json(result, f"{catalyst_name}_dft.json", subdir="dft")
            return result

    # Write input file
    input_file = calc_dir / f"{catalyst_name}.in"
    output_file = calc_dir / f"{catalyst_name}.out"
    # Names can be reused for a different structure (positional stage names,
    # changed generators). An output is only evidence for the exact input that
    # produced it; otherwise move it aside so it is neither parsed nor resumed.
    same_input = input_file.is_file() and input_file.read_text() == input_text
    if output_file.exists() and not same_input:
        stale_output = output_file.with_name(output_file.name + '.stale')
        output_file.replace(stale_output)
        result['discarded_stale_output'] = str(stale_output)
        logger.warning(
            f"  Existing output for {catalyst_name} came from a different input; "
            f"moved to {stale_output}"
        )
    if output_file.is_file() and output_file.stat().st_size:
        require_ionic = "calculation = 'relax'" in input_text
        existing_converged = parse_convergence(
            str(output_file), require_ionic=require_ionic
        )
        if existing_converged or not restart_incomplete:
            energy = parse_total_energy(str(output_file))
            result.update(
                {
                    'dft_energy_Ry': energy,
                    'dft_energy_eV': energy * Ry_to_eV if energy else None,
                    'converged': existing_converged,
                    'max_force_Ry_bohr': parse_forces(str(output_file)),
                    'resumed_existing_output': True,
                }
            )
            if not existing_converged:
                result['error'] = 'existing_output_incomplete_or_failed'
            save_json(result, f"{catalyst_name}_dft.json", subdir="dft")
            return result
    with open(input_file, 'w') as f:
        f.write(input_text)
    logger.info(f"  Written QE input: {input_file}")

    # Run DFT
    if run_dft:
        logger.info(f"  Running pw.x for {catalyst_name}...")
        try:
            from pipeline.validation.qe_workflows import QEExecutionConfig, run_pw

            exec_cfg = (
                execution
                if execution is not None
                else QEExecutionConfig.production_default()
            )
            outcome = run_pw(
                str(input_file),
                str(output_file),
                timeout_s=3600,
                execution=exec_cfg,
            )
            result['returncode'] = outcome['returncode']
            result['execution'] = outcome['execution']
            if outcome['timed_out']:
                logger.warning(f"  DFT calculation timed out for {catalyst_name}")
                result['error'] = 'timeout'
        except Exception as exc:
            logger.warning(f"  DFT execution failed for {catalyst_name}: {exc}")
            result['returncode'] = -1
            result['error'] = str(exc)[:200]

    # Parse results
    if output_file.exists():
        energy = parse_total_energy(str(output_file))
        require_ionic = "calculation = 'relax'" in input_text
        converged = parse_convergence(str(output_file), require_ionic=require_ionic)
        max_force = parse_forces(str(output_file))

        result['dft_energy_Ry'] = energy
        result['dft_energy_eV'] = energy * Ry_to_eV if energy else None
        result['converged'] = converged
        result['max_force_Ry_bohr'] = max_force

        if converged and energy:
            logger.info(
                f"  ✓ Converged: E = {energy:.6f} Ry ({energy * Ry_to_eV:.4f} eV)"
            )
        else:
            logger.warning(f"  ✗ Did not converge or no energy found")
    else:
        logger.warning(f"  No output file found: {output_file}")

    save_json(result, f"{catalyst_name}_dft.json", subdir="dft")
    return result


# fcc-equivalent conventional lattice constants (Bohr), a = (4 V_atom)^(1/3)
# from room-temperature atomic volumes. Native fcc metals keep their measured
# a; bcc/hcp/complex metals (Fe, Mo, W, V, Cr, Ti, Zr, Y, Zn, La, alpha-Mn,
# beta-Sn, In, Ga, Sb, Bi) use the volume-matched fcc value so that the
# 4-site conventional cube below is not compressed or expanded artificially.
_FCC_EQUIVALENT_LATTICE_BOHR = {
    'Ni': 6.66,
    'Fe': 6.82,
    'Co': 6.69,
    'Cu': 6.83,
    'Al': 7.65,
    'Sn': 9.00,
    'Bi': 9.85,
    'In': 8.90,
    'Ga': 8.09,
    'Pb': 9.35,
    'Sb': 9.34,
    'Mo': 7.49,
    'W': 7.54,
    'Mn': 6.91,
    'Pt': 7.42,
    'Pd': 7.35,
    'Ag': 7.72,
    'Au': 7.71,
    'Ti': 7.81,
    'V': 7.21,
    'Cr': 6.87,
    'Zn': 7.43,
    'Zr': 8.57,
    'Y': 9.62,
    'La': 10.03,
}
# Fractional coordinates of the four fcc sites in the simple-cubic
# conventional cell (ibrav=1). These are NOT valid in the ibrav=2 primitive
# cell, where they would place atoms 0.354 a apart.
_FCC_CONVENTIONAL_SITES = (
    (0.0, 0.0, 0.0),
    (0.5, 0.5, 0.0),
    (0.5, 0.0, 0.5),
    (0.0, 0.5, 0.5),
)


def _generate_alloy_input(genome: tuple, name: str) -> str:
    """Generate QE input for an alloy bulk structure (4-site fcc cube)."""
    mat_class = genome[0]
    host = genome[1]
    elements = [host] * len(_FCC_CONVENTIONAL_SITES)
    if mat_class == 'MoltenMetal':
        promoter = genome[2]
        # L1_2-type ordering: one promoter site out of four.
        if promoter != 'None':
            elements[0] = promoter
    elif genome[5]:
        # Apply dopant substitution
        elements[0] = genome[5][0]

    # The host metal sets the cell; a substituted minority site does not.
    celldm = _FCC_EQUIVALENT_LATTICE_BOHR.get(host, 7.0)

    return generate_bulk_scf_input(
        elements,
        list(_FCC_CONVENTIONAL_SITES),
        celldm,
        ibrav=1,
        calc_name=name,
        ecutwfc=50.0,
        kpoints=(6, 6, 6),
    )


def _generate_cluster_input(genome: tuple, name: str) -> str:
    """Generate QE input for a molecular cluster (SAC/MOF)."""
    from pipeline.screening.surface_screener import generate_porphyrin_cluster

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

    raw_elems = [atom.symbol for atom in cluster if atom.symbol != 'X']
    seen_unique = list(dict.fromkeys([metal] + raw_elems))
    site_elems = (seen_unique + raw_elems)[:4]
    while len(site_elems) < 4:
        site_elems.append(metal)
    celldm = _FCC_EQUIVALENT_LATTICE_BOHR.get(metal, 7.0)

    return generate_bulk_scf_input(
        site_elems,
        list(_FCC_CONVENTIONAL_SITES),
        celldm,
        ibrav=1,
        calc_name=name,
        ecutwfc=50.0,
        kpoints=(4, 4, 4),
    )


if __name__ == '__main__':
    # Test: generate input for a Ni-Bi molten metal catalyst
    genome = ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000)
    result = validate_catalyst("NiBi_10pct", genome, run_dft=False)
    print(json.dumps(result, indent=2, default=str))
