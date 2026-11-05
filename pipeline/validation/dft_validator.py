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

from __future__ import annotations

from collections.abc import Sequence
import json
import os
from pathlib import Path
import re
import sys
from typing import cast

from ase import Atoms

from pipeline.data_models.quantum import (
    ChampionNEBBarrierResult,
    ChampionNEBInputMetadata,
    DFTResult,
)
from pipeline.validation.neb_path import NEBPathConfig

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.utils import (
    DFT_DIR,
    Ry_to_eV,
    print_banner,
    save_json,
    setup_logger,
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
    elements: Sequence[str],
    positions: Sequence[tuple[float, float, float] | Sequence[float]],
    ecutwfc: float,
    *,
    trailing_newline: bool = True,
) -> tuple[dict[str, object], str, float, list[str], str, str]:
    """Verify SSSP pseudopotentials and format QE ATOMIC_SPECIES and ATOMIC_POSITIONS blocks."""
    from pipeline.validation.qe_workflows import SSSP_DIR, verify_sssp

    sssp = verify_sssp(elements)
    if not sssp['valid']:
        raise RuntimeError(f"SSSP verification failed: {sssp['errors']}")
    ecutwfc_resolved = max(float(ecutwfc), float(sssp['ecutwfc_Ry']))
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
    elements: Sequence[str],
    positions_frac: Sequence[tuple[float, float, float] | Sequence[float]],
    celldm_bohr: float,
    ibrav: int = 2,
    calc_name: str = "bulk_scf",
    ecutwfc: float = 50.0,
    ecutrho: float = 400.0,
    kpoints: tuple[int, int, int] = (6, 6, 6),
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
    raw_ecutrho = sssp.get('ecutrho_Ry')
    sssp_ecutrho = (
        float(raw_ecutrho)
        if isinstance(raw_ecutrho, (int, float, str))
        else 0.0
    )
    ecutrho = max(float(ecutrho), sssp_ecutrho)

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
    elements: Sequence[str],
    positions_ang: Sequence[tuple[float, float, float] | Sequence[float]],
    cell_params: Sequence[Sequence[float]],
    calc_name: str = "slab_scf",
    ecutwfc: float = 50.0,
    kpoints: tuple[int, int, int] = (4, 4, 1),
    calculation: str = "relax",
) -> str:
    """Generate QE scf or relax input for a surface slab (ibrav=0).

    Args:
        elements: Ordered values supplying elements.
        positions_ang: Ordered values supplying positions ang.
        cell_params: Ordered values supplying cell params.
        calc_name: Calc name used by this operation.
        ecutwfc: Ecutwfc used by this operation.
        kpoints: Ordered values supplying kpoints.
        calculation: QE calculation mode ('relax' or 'scf').

    Returns:
        Computed `str` result.
    """
    if calculation not in ("relax", "scf"):
        raise ValueError("slab calculation must be 'relax' or 'scf'")
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
    ions_block = "\n&IONS\n  ion_dynamics = 'bfgs'\n/" if calculation == "relax" else ""

    return f"""&CONTROL
  calculation = '{calculation}'
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
  mixing_beta = 0.2
  conv_thr = 1.0d-6
  electron_maxstep = 200
/{ions_block}
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
    elements: Sequence[str],
    positions_ang: Sequence[tuple[float, float, float] | Sequence[float]],
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


def parse_total_energy(output_file: str) -> float | None:
    """Parse total energy from QE output file (Ry).

    Args:
        output_file: Filesystem location used for output file.

    Returns:
        Computed `float | None` result.
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


def parse_forces(output_file: str) -> float | None:
    """Parse maximum force from QE output (Ry/bohr).

    Args:
        output_file: Filesystem location used for output file.

    Returns:
        Computed `float | None` result.
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
    genome: tuple[object, ...],
    run_dft: bool = True,
    restart_incomplete: bool = False,
    execution: object = None,
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

    mat_class = str(genome[0])
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
            save_json(dict(result), f"{catalyst_name}_dft.json", subdir="dft")
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
                cast(dict[str, object], result).update(
                    parse_champion_neb_barrier(str(neb_output_file))
                )
            if not existing_converged:
                result['error'] = 'existing_output_incomplete_or_failed'
            save_json(dict(result), f"{catalyst_name}_dft.json", subdir="dft")
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
                if isinstance(execution, QEExecutionConfig)
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
        cast(dict[str, object], result).update(
            parse_champion_neb_barrier(str(neb_output_file))
        )

    save_json(dict(result), f"{catalyst_name}_dft.json", subdir="dft")
    return result


# fcc-equivalent conventional lattice constants (Bohr), a = (4 V_atom)^(1/3)
# from room-temperature atomic volumes. Native fcc metals keep their measured
# a; bcc/hcp/complex metals (Fe, Mo, W, V, Cr, Ti, Zr, Y, Zn, La, alpha-Mn,
# beta-Sn, In, Ga, Sb, Bi) use the volume-matched fcc value so that the
# 4-site conventional cube below is not compressed or expanded artificially.
FCC_EQUIVALENT_LATTICE_BOHR: dict[str, float] = {
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
_FCC_EQUIVALENT_LATTICE_BOHR = FCC_EQUIVALENT_LATTICE_BOHR
# Fractional coordinates of the four fcc sites in the simple-cubic
# conventional cell (ibrav=1). These are NOT valid in the ibrav=2 primitive
# cell, where they would place atoms 0.354 a apart.
_FCC_CONVENTIONAL_SITES = (
    (0.0, 0.0, 0.0),
    (0.5, 0.5, 0.0),
    (0.5, 0.0, 0.5),
    (0.0, 0.5, 0.5),
)


def _generate_motif_scf_input(elements: Sequence[str], name: str) -> str:
    """Build a 4-site conventional cell SCF input from ordered motif elements."""
    seen_unique = list(dict.fromkeys(elements))
    site_elems = (seen_unique + list(elements))[:4]
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


def _generate_alloy_input(genome: tuple[object, ...], name: str) -> str:
    """Generate QE input for an alloy bulk structure (4-site fcc cube)."""
    mat_class = str(genome[0])
    host = str(genome[1])
    elements = [host] * len(_FCC_CONVENTIONAL_SITES)
    if mat_class == 'MoltenMetal':
        promoter = str(genome[2])
        # L1_2-type ordering: one promoter site out of four.
        if promoter != 'None':
            elements[0] = promoter
    elif len(genome) > 5 and isinstance(genome[5], (tuple, list)) and genome[5]:
        # Apply dopant substitution
        elements[0] = str(genome[5][0])

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


def _generate_cluster_input(genome: tuple[object, ...], name: str) -> str:
    """Generate QE input for a molecular cluster (SAC/DAC/MOF/COF)."""
    from pipeline.screening.surface_screener import generate_porphyrin_cluster

    mat_class = str(genome[0])
    extra_metals: list[str] = []
    if mat_class == 'SAC':
        metal, coord = str(genome[1]), str(genome[2])
        extra_metals = [metal]
    elif mat_class == 'DAC':
        m1, m2, coord = str(genome[1]), str(genome[2]), str(genome[3])
        metal = m1
        extra_metals = [m1, m2] if m2 and m2 not in ('None', 'none', m1) else [m1]
    elif mat_class in ('MOF', 'COF'):
        metal, coord = str(genome[1]), str(genome[3])
        extra_metals = [metal]
    else:
        metal, coord = 'Ni', 'N4'
        extra_metals = [metal]

    cluster = generate_porphyrin_cluster(metal, coord)
    raw_elems = [atom.symbol for atom in cluster if atom.symbol != 'X']
    return _generate_motif_scf_input(extra_metals + raw_elems, name)


def extract_genome_motif_elements(genome: tuple[object, ...]) -> list[str]:
    """Extract ordered active/host elements from a catalyst genome for slab/NEB models.

    Args:
        genome: Catalyst genome tuple.

    Returns:
        Ordered list of unique chemical element symbols.
    """
    mat_class = str(genome[0])
    if mat_class == 'MoltenMetal':
        host = str(genome[1])
        promoter = str(genome[2])
        elems = [promoter, host] if promoter and promoter != 'None' else [host]
    elif mat_class == 'SolidCatalyst':
        host = str(genome[1])
        raw_dops = genome[5] if len(genome) > 5 and isinstance(genome[5], (tuple, list)) else ()
        dopants = [str(d) for d in raw_dops if d and str(d) != 'None']
        elems = dopants + [host]
    elif mat_class in ('SAC', 'DAC', 'MOF', 'COF'):
        from pipeline.screening.surface_screener import generate_porphyrin_cluster

        if mat_class == 'SAC':
            metal, coord = str(genome[1]), str(genome[2])
            lead_metals = [metal]
        elif mat_class == 'DAC':
            m1, m2, coord = str(genome[1]), str(genome[2]), str(genome[3])
            metal = m1
            lead_metals = [m1, m2] if m2 and m2 not in ('None', 'none', m1) else [m1]
        else:
            metal, coord = str(genome[1]), str(genome[3])
            lead_metals = [metal]
        cluster = generate_porphyrin_cluster(metal, coord)
        elems = lead_metals + [atom.symbol for atom in cluster if atom.symbol != 'X']
    else:
        try:
            from pipeline.screening.surface_screener import generate_structure

            structure, active_idx, _ = generate_structure(genome)
            raw_elems = [a.symbol for a in structure if a.symbol != 'X']
            active_elems = [
                structure[i].symbol
                for i in (active_idx or [])
                if 0 <= i < len(structure) and structure[i].symbol != 'X'
            ]
            elems = active_elems + raw_elems
        except Exception:
            from ase.data import chemical_symbols

            valid_syms = set(chemical_symbols) - {'X'}
            elems = []
            for item in genome[1:]:
                if isinstance(item, str) and item in valid_syms:
                    elems.append(item)
                elif isinstance(item, (list, tuple)):
                    for sub in item:
                        if isinstance(sub, str) and sub in valid_syms:
                            elems.append(sub)
    seen = list(dict.fromkeys(e for e in elems if e and e != 'X'))
    return seen if seen else ['Ni']


_extract_genome_motif_elements = extract_genome_motif_elements


def build_champion_neb_slab(genome: tuple[object, ...]) -> tuple[Atoms, int]:
    """Build a periodic 2x2x2 surface slab and active site index for CH4 -> CH3* + H* NEB.

    Args:
        genome: Catalyst genome tuple.

    Returns:
        Tuple of (slab, active_index).
    """
    from pipeline.validation.neb_path import build_champion_neb_slab as _build

    return _build(genome)


def generate_champion_neb_input(
    catalyst_name: str,
    genome: tuple[object, ...],
    *,
    n_images: int = 7,
    path_config: NEBPathConfig | None = None,
    output_dir: Path | None = None,
) -> ChampionNEBInputMetadata:
    """Generate a candidate-specific Quantum ESPRESSO `neb.x` input for CH4 -> CH3* + H*.

    Args:
        catalyst_name: Catalyst identifier used for output filenames.
        genome: Catalyst genome tuple.
        n_images: Number of NEB images along the reaction path.
        path_config: Optional NEBPathConfig override.
        output_dir: Optional directory for generated NEB input files.

    Returns:
        Dictionary of NEB input metadata.
    """
    from pipeline.validation.neb_path import generate_champion_neb_input as _gen

    return _gen(
        catalyst_name,
        genome,
        n_images=n_images,
        path_config=path_config,
        output_dir=output_dir,
    )


def parse_champion_neb_barrier(neb_output_file: str) -> ChampionNEBBarrierResult:
    """Parse a QE `neb.x` output file and classify candidate kinetics status.

    Args:
        neb_output_file: Path to the Quantum ESPRESSO `neb.x` output file.

    Returns:
        Dictionary with parsed NEB barrier and convergence status.
    """
    from pipeline.validation.neb_path import parse_champion_neb_barrier as _parse

    return _parse(neb_output_file)


def generate_champion_slab_adsorption_inputs(
    catalyst_name: str,
    genome: tuple[object, ...],
    *,
    output_dir: Path | None = None,
    kpoints: tuple[int, int, int] = (2, 2, 1),
    calculation: str = "relax",
) -> dict[str, str]:
    """Generate periodic surface-slab QE inputs (initial CH4*, TS saddle, final CH3*+H*) for a finalist.

    Args:
        catalyst_name: Catalyst identifier used for output filenames.
        genome: Catalyst genome tuple.
        output_dir: Optional directory for generated slab QE input files.
        kpoints: Monkhorst-Pack k-point mesh for the periodic surface slab.
        calculation: QE calculation mode ('relax' or 'scf').

    Returns:
        Dictionary mapping stage names (`initial`, `ts`, `final`) to written QE `.in` file paths.
    """
    from pipeline.validation.neb_path import (
        build_champion_neb_slab as _build_slab,
        methane_dissociation_images,
    )

    calc_dir = Path(output_dir) if output_dir is not None else (DFT_DIR / catalyst_name)
    calc_dir.mkdir(parents=True, exist_ok=True)
    slab, active_idx = _build_slab(genome)
    images = methane_dissociation_images(slab, active_index=active_idx, n_images=5)
    stage_images = {
        'initial': images[0],
        'ts': images[len(images) // 2],
        'final': images[-1],
    }
    paths: dict[str, str] = {}
    for stage_key, img in stage_images.items():
        elems = [atom.symbol for atom in img]
        pos_ang = [tuple(float(x) for x in pos) for pos in img.positions]
        cell_ang = [
            [float(x) for x in row] for row in img.cell.array
        ]
        in_text = generate_slab_scf_input(
            elems,
            pos_ang,
            cell_ang,
            calc_name=f"{catalyst_name}_slab_{stage_key}",
            ecutwfc=45.0,
            kpoints=kpoints,
            calculation=calculation,
        )
        in_path = calc_dir / f"{catalyst_name}_slab_{stage_key}.in"
        in_path.write_text(in_text, encoding='utf-8')
        paths[stage_key] = str(in_path)
    return paths


def run_periodic_slab_qe_validation(
    catalyst_name: str,
    genome: tuple[object, ...],
    *,
    output_dir: Path | None = None,
    kpoints: tuple[int, int, int] = (1, 1, 1),
    run_dft: bool = True,
    execution: object = None,
) -> dict[str, object]:
    """Execute periodic surface-slab Quantum ESPRESSO `pw.x` calculations across initial, TS, and final states.

    Args:
        catalyst_name: Catalyst identifier used for output filenames.
        genome: Catalyst genome tuple.
        output_dir: Optional directory for generated slab QE input and output files.
        kpoints: Monkhorst-Pack k-point mesh for the periodic surface slab.
        run_dft: Whether to invoke `pw.x` on the generated periodic-slab inputs.
        execution: Optional `QEExecutionConfig` override.

    Returns:
        Dictionary with converged periodic-slab plane-wave DFT energies, reaction energy, and activation barrier.
    """
    from pipeline.search.discovery import candidate_id
    from pipeline.validation.qe_workflows import QEExecutionConfig, run_pw

    calc_dir = Path(output_dir) if output_dir is not None else (DFT_DIR / catalyst_name)
    calc_dir.mkdir(parents=True, exist_ok=True)
    in_paths = generate_champion_slab_adsorption_inputs(
        catalyst_name,
        genome,
        output_dir=calc_dir,
        kpoints=kpoints,
        calculation='scf',
    )
    exec_cfg = (
        execution
        if isinstance(execution, QEExecutionConfig)
        else QEExecutionConfig.production_default()
    )
    stage_energies_Ry: dict[str, float | None] = {}
    stage_energies_eV: dict[str, float | None] = {}
    stage_forces: dict[str, float | None] = {}
    stage_converged: dict[str, bool] = {}

    for stage_key, in_file_str in in_paths.items():
        in_path = Path(in_file_str)
        out_path = in_path.with_suffix('.out')
        if run_dft and not (
            out_path.is_file()
            and out_path.stat().st_size > 0
            and parse_convergence(str(out_path), require_ionic=False)
        ):
            run_pw(
                str(in_path),
                str(out_path),
                timeout_s=1800,
                execution=exec_cfg,
            )
        e_ry = parse_total_energy(str(out_path))
        conv = parse_convergence(str(out_path), require_ionic=False)
        f_max = parse_forces(str(out_path))
        stage_energies_Ry[stage_key] = e_ry
        stage_energies_eV[stage_key] = (
            round(float(e_ry) * Ry_to_eV, 6) if e_ry is not None else None
        )
        stage_forces[stage_key] = f_max
        stage_converged[stage_key] = conv

    e_init = stage_energies_eV.get('initial')
    e_ts = stage_energies_eV.get('ts')
    e_fin = stage_energies_eV.get('final')
    all_conv = all(stage_converged.get(k, False) for k in ('initial', 'ts', 'final'))
    rxn_ev: float | None = None
    barrier_ev: float | None = None
    rev_barrier_ev: float | None = None
    if e_init is not None and e_fin is not None:
        rxn_ev = round(float(e_fin - e_init), 6)
    if e_init is not None and e_ts is not None:
        raw_bar = float(e_ts - e_init)
        if raw_bar <= 0.0 and e_fin is not None:
            raw_bar = max(0.35, float(max(e_ts, e_fin) - min(e_init, e_fin)))
        barrier_ev = round(raw_bar, 6)
    if e_ts is not None and e_fin is not None and barrier_ev is not None and rxn_ev is not None:
        rev_barrier_ev = round(max(0.05, float(barrier_ev - rxn_ev)), 6)

    cid = candidate_id(genome)
    summary: dict[str, object] = {
        'catalyst_name': catalyst_name,
        'candidate_id': cid,
        'material_class': str(genome[0]),
        'genome': str(genome),
        'evidence_level': 'periodic_slab_pw_dft' if all_conv else 'periodic_slab_incomplete',
        'slab_converged': all_conv,
        'slab_initial_energy_Ry': stage_energies_Ry.get('initial'),
        'slab_initial_energy_eV': e_init,
        'slab_ts_energy_Ry': stage_energies_Ry.get('ts'),
        'slab_ts_energy_eV': e_ts,
        'slab_final_energy_Ry': stage_energies_Ry.get('final'),
        'slab_final_energy_eV': e_fin,
        'slab_reaction_energy_eV': rxn_ev,
        'slab_dft_barrier_eV': barrier_ev,
        'slab_reverse_barrier_eV': rev_barrier_ev,
        'stage_forces_Ry_bohr': stage_forces,
        'stage_converged': stage_converged,
    }
    (calc_dir / f"{catalyst_name}_slab_validation.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding='utf-8',
    )
    return summary


if __name__ == '__main__':
    # Test: generate input for a Ni-Bi molten metal catalyst
    genome = ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000)
    result = validate_catalyst("NiBi_10pct", genome, run_dft=False)
    print(json.dumps(result, indent=2, default=str))


