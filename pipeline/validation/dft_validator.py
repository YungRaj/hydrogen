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
    'Sc': 44.956,
    'Ba': 137.327,
    'Pr': 140.908,
    'Sm': 150.36,
    'Gd': 157.25,
    'Hf': 178.49,
    'Ta': 180.948,
    'Re': 186.207,
    'Ir': 192.217,
}


# ═══════════════════════════════════════════════════════════════════════════════
# QE INPUT GENERATORS
# ═══════════════════════════════════════════════════════════════════════════════


def _prepare_qe_species_and_atoms(
    elements: List[str],
    positions: List[Tuple],
    ecutwfc: float,
    *,
    trailing_newline: bool = True,
) -> tuple[dict, str, float, list[str], str, str]:
    """Verify SSSP pseudopotentials and format QE ATOMIC_SPECIES and ATOMIC_POSITIONS blocks."""
    from pipeline.validation.qe_workflows import SSSP_DIR, verify_sssp

    sssp = verify_sssp(elements)
    if not sssp['valid']:
        raise RuntimeError(f"SSSP verification failed: {sssp['errors']}")
    ecutwfc_resolved = max(float(ecutwfc), sssp['ecutwfc_Ry'])
    unique_elements = sorted(set(elements))
    suffix = '\n' if trailing_newline else ''
    species_block = (
        '\n'.join(
            f"  {elem}  {ATOMIC_MASSES_QE.get(elem, 50.0):.3f}  "
            f"{sssp['records'][elem]['filename']}"
            for elem in unique_elements
        )
        + suffix
    )
    atoms_block = (
        '\n'.join(
            f"  {elem}  {pos[0]:.10f}  {pos[1]:.10f}  {pos[2]:.10f}"
            for elem, pos in zip(elements, positions)
        )
        + suffix
    )
    return (
        sssp,
        str(SSSP_DIR),
        ecutwfc_resolved,
        unique_elements,
        species_block,
        atoms_block,
    )


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
    kshift = tuple(0 if mesh == 1 else 1 for mesh in kpoints)
    (
        sssp,
        sssp_dir,
        ecutwfc,
        unique_elements,
        species_block,
        atoms_block,
    ) = _prepare_qe_species_and_atoms(elements, positions_frac, ecutwfc)
    ntyp = len(unique_elements)
    ecutrho = max(float(ecutrho), sssp['ecutrho_Ry'])

    # smearing = 'mv' with degauss = 0.02 Ry (~0.272 eV): Marzari-Vanderbilt
    # cold smearing (Phys. Rev. Lett. 82, 3296, 1999) for metallic Fermi surfaces.
    return f"""&CONTROL
  calculation = 'scf'
  prefix = '{calc_name}'
  outdir = './tmp'
  pseudo_dir = '{sssp_dir}'
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
    kshift = tuple(0 if mesh == 1 else 1 for mesh in kpoints)
    (
        sssp,
        sssp_dir,
        ecutwfc,
        unique_elements,
        species_block,
        atoms_block,
    ) = _prepare_qe_species_and_atoms(elements, positions_ang, ecutwfc)
    ntyp = len(unique_elements)
    ecutrho = sssp['ecutrho_Ry']

    # Quantum ESPRESSO rejects spin-polarized calculations unless at least one
    # species has a non-zero starting magnetization.  Use a conservative seed
    # for every species; this initializes the spin channel without asserting a
    # particular converged magnetic state.
    magnetization_block = "\n".join(
        f"  starting_magnetization({index}) = 0.05"
        for index, _ in enumerate(unique_elements, start=1)
    )
    cell_block = "".join(
        f"  {row[0]:.10f}  {row[1]:.10f}  {row[2]:.10f}\n" for row in cell_params
    )

    return f"""&CONTROL
  calculation = 'relax'
  prefix = '{calc_name}'
  outdir = './tmp'
  pseudo_dir = '{sssp_dir}'
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

    (
        sssp,
        sssp_dir,
        ecutwfc,
        unique_elements,
        species_block,
        atoms_block,
    ) = _prepare_qe_species_and_atoms(
        elements, positions_ang, ecutwfc, trailing_newline=False
    )
    ions = "\n&IONS\n  ion_dynamics = 'bfgs'\n/" if calculation == 'relax' else ''

    # 0.529177210903 Å/Bohr: CODATA Bohr radius a_0 converting cubic cell_size_ang to celldm(1) in Bohr
    return f"""&CONTROL
  calculation = '{calculation}'
  prefix = '{calc_name}'
  outdir = './tmp'
  pseudo_dir = '{sssp_dir}'
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
            input_text = _generate_motif_scf_input(
                active_elems + raw_elems, catalyst_name
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
            neb_output_file = calc_dir / f"{catalyst_name}.neb.out"
            if neb_output_file.is_file() and neb_output_file.stat().st_size:
                result.update(parse_champion_neb_barrier(str(neb_output_file)))
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

    neb_output_file = calc_dir / f"{catalyst_name}.neb.out"
    if neb_output_file.is_file() and neb_output_file.stat().st_size:
        result.update(parse_champion_neb_barrier(str(neb_output_file)))

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
    'Ru': 7.16,
    'Rh': 7.19,
    'Ir': 7.26,
    'Re': 7.35,
    'Ta': 7.85,
    'Nb': 7.85,
    'Ge': 8.16,
    'Hf': 8.51,
    'Sc': 8.80,
    'Te': 9.60,
    'Gd': 9.65,
    'Sm': 9.69,
    'Ce': 9.75,
    'Pr': 9.78,
    'Ba': 10.50,
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


def _generate_motif_scf_input(elements: List[str], name: str) -> str:
    """Build a 4-site conventional cell SCF input from ordered motif elements."""
    seen_unique = list(dict.fromkeys(elements))
    site_elems = (seen_unique + elements)[:4]
    while len(site_elems) < 4:
        site_elems.append(site_elems[0] if site_elems else 'Ni')
    celldm = _FCC_EQUIVALENT_LATTICE_BOHR.get(site_elems[0], 7.0)
    return generate_bulk_scf_input(
        site_elems,
        list(_FCC_CONVENTIONAL_SITES),
        celldm,
        ibrav=1,
        calc_name=name,
        ecutwfc=50.0,
        kpoints=(2, 2, 2),
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
        kpoints=(2, 2, 2),
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
    return _generate_motif_scf_input([metal] + raw_elems, name)


def _extract_genome_motif_elements(genome: tuple) -> List[str]:
    """Extract ordered active/host elements from a catalyst genome for slab/NEB models."""
    mat_class = str(genome[0])
    if mat_class == 'MoltenMetal':
        host = str(genome[1])
        promoter = str(genome[2])
        elems = [promoter, host] if promoter and promoter != 'None' else [host]
    elif mat_class == 'SolidCatalyst':
        host = str(genome[1])
        dopants = [str(d) for d in (genome[5] or []) if d and str(d) != 'None']
        elems = dopants + [host]
    elif mat_class in ('SAC', 'DAC', 'MOF', 'COF'):
        from pipeline.screening.surface_screener import generate_porphyrin_cluster

        if mat_class == 'SAC':
            metal, coord = str(genome[1]), str(genome[2])
        else:
            metal, coord = str(genome[1]), str(genome[3])
        cluster = generate_porphyrin_cluster(metal, coord)
        elems = [metal] + [atom.symbol for atom in cluster if atom.symbol != 'X']
    else:
        from pipeline.screening.surface_screener import generate_structure

        structure, active_idx, _ = generate_structure(genome)
        raw_elems = [a.symbol for a in structure if a.symbol != 'X']
        active_elems = [
            structure[i].symbol
            for i in (active_idx or [])
            if 0 <= i < len(structure) and structure[i].symbol != 'X'
        ]
        elems = active_elems + raw_elems
    seen = list(dict.fromkeys(e for e in elems if e and e != 'X'))
    return seen if seen else ['Ni']


def build_champion_neb_slab(genome: tuple):
    """Build a periodic 2x2x2 surface slab and active site index for CH4 -> CH3* + H* NEB.

    Args:
        genome: Encoded catalyst composition and structural configuration.

    Returns:
        Tuple `(slab, active_index)` where `slab` is an ASE `Atoms` surface model
        with bottom-layer `FixAtoms` constraints and `active_index` is the top-layer
        catalytic site index.
    """
    from ase.build import fcc111
    from ase.constraints import FixAtoms

    motif_elems = _extract_genome_motif_elements(genome)
    host_elem = (
        str(genome[1])
        if str(genome[0]) in ('MoltenMetal', 'SolidCatalyst')
        and str(genome[1]) in _FCC_EQUIVALENT_LATTICE_BOHR
        else motif_elems[0]
    )
    celldm_bohr = _FCC_EQUIVALENT_LATTICE_BOHR.get(host_elem, 7.0)
    a_ang = celldm_bohr * 0.529177210903
    slab = fcc111(host_elem, size=(2, 2, 2), a=a_ang, vacuum=6.0, periodic=True)
    top_indices = [i for i, atom in enumerate(slab) if atom.tag == 1]
    bottom_indices = [i for i, atom in enumerate(slab) if atom.tag > 1]
    if not top_indices:
        top_indices = list(range(len(slab) // 2, len(slab)))
        bottom_indices = list(range(0, len(slab) // 2))
    # Populate top surface sites with the candidate's active motif elements, placing
    # the primary active site at the final top-layer index.
    ordered_top = list(reversed(top_indices))
    for k, idx in enumerate(ordered_top):
        slab[idx].symbol = motif_elems[k % len(motif_elems)]
    if bottom_indices:
        slab.set_constraint(FixAtoms(indices=bottom_indices))
    active_index = ordered_top[0]
    return slab, active_index


def generate_champion_neb_input(
    catalyst_name: str,
    genome: tuple,
    *,
    n_images: int = 7,
    path_config=None,
    output_dir: Optional[Path] = None,
) -> dict:
    """Generate a candidate-specific Quantum ESPRESSO `neb.x` input for CH4 -> CH3* + H*.

    Args:
        catalyst_name: Identifier for the champion catalyst.
        genome: Encoded catalyst composition and structural configuration.
        n_images: Number of NEB images along the dissociation pathway (>= 5).
        path_config: Optional `NEBPathConfig` override.
        output_dir: Optional directory override; defaults to `DFT_DIR / catalyst_name`.

    Returns:
        Dictionary containing the written `.neb.in` path, image count, active site index,
        motif elements, and verified SSSP metadata.
    """
    from pipeline.validation.neb_path import NEBPathConfig, methane_dissociation_images
    from pipeline.validation.qe_workflows import write_qe_neb_input

    calc_dir = Path(output_dir) if output_dir is not None else (DFT_DIR / catalyst_name)
    calc_dir.mkdir(parents=True, exist_ok=True)
    slab, active_index = build_champion_neb_slab(genome)
    images = methane_dissociation_images(
        slab, active_index=active_index, n_images=n_images
    )
    cfg = path_config if path_config is not None else NEBPathConfig()
    neb_input_path = calc_dir / f"{catalyst_name}.neb.in"
    meta = write_qe_neb_input(
        images,
        str(neb_input_path),
        prefix=f"{catalyst_name}_neb",
        path_config=cfg,
    )
    return {
        **meta,
        'catalyst_name': catalyst_name,
        'active_index': active_index,
        'elements': sorted(set(images[0].get_chemical_symbols())),
    }


def parse_champion_neb_barrier(neb_output_file: str) -> dict:
    """Parse a QE `neb.x` output file and classify candidate kinetics status.

    Only a cleanly converged `neb.x` output (`JOB DONE` and
    `neb: convergence achieved`) with a finite positive forward barrier
    promotes `kinetics_status` to `'dft_neb_complete'`.

    Args:
        neb_output_file: Path to the QE `neb.x` output file.

    Returns:
        Dictionary with `neb_converged`, `neb_forward_barrier_eV`,
        `neb_reverse_barrier_eV`, `neb_path_metric`, and `kinetics_status`.
    """
    from pipeline.validation.qe_workflows import parse_neb_result

    if not os.path.exists(neb_output_file):
        return {
            'neb_converged': False,
            'neb_forward_barrier_eV': None,
            'neb_reverse_barrier_eV': None,
            'neb_path_metric': None,
            'kinetics_status': 'qe_dft_calibrated_bep',
        }
    parsed = parse_neb_result(str(neb_output_file))
    conv = bool(parsed.get('converged', False))
    fwd = parsed.get('forward_barrier_eV')
    rev = parsed.get('reverse_barrier_eV')
    metric = parsed.get('path_metric')
    valid_barrier = (
        conv
        and isinstance(fwd, (int, float))
        and float(fwd) > 0.0
    )
    return {
        'neb_converged': conv,
        'neb_forward_barrier_eV': float(fwd) if isinstance(fwd, (int, float)) else None,
        'neb_reverse_barrier_eV': float(rev) if isinstance(rev, (int, float)) else None,
        'neb_path_metric': float(metric) if isinstance(metric, (int, float)) else None,
        'kinetics_status': (
            'dft_neb_complete' if valid_barrier else 'qe_dft_calibrated_bep'
        ),
    }


if __name__ == '__main__':
    # Test: generate input for a Ni-Bi molten metal catalyst
    genome = ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000)
    result = validate_catalyst("NiBi_10pct", genome, run_dft=False)
    print(json.dumps(result, indent=2, default=str))

