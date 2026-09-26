#!/usr/bin/env python3
"""Run reproducible Quantum ESPRESSO bulk equation-of-state benchmarks."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from pipeline.validation.qe_workflows import (
    QEExecutionConfig,
    SSSP_DIR,
    build_qe_command,
    resolve_qe_executable,
    verify_sssp,
)


RY_TO_EV = 13.605693122994
BOHR_PER_ANGSTROM = 1.889726125458
REFERENCE_LATTICES = {"Fe": 2.8664, "Ni": 3.5238}
STRUCTURES = {"Fe": "bcc", "Ni": "fcc"}
IBRAV = {"Fe": 3, "Ni": 2}
REFERENCE_IDS = {
    "Fe": "nist-fe-bcc-lattice-293k",
    "Ni": "nist-ni-fcc-lattice-293k",
}


def _digest(path: Path) -> str:
    value = hashlib.sha256()
    value.update(path.read_bytes())
    return value.hexdigest()


def qe_bulk_input(element: str, lattice_A: float, prefix: str) -> str:
    """Build one spin-polarized primitive-cell QE input.

    Args:
        element: Fe or Ni element symbol.
        lattice_A: Conventional cubic lattice parameter in ångströms.
        prefix: Unique Quantum ESPRESSO scratch prefix.

    Returns:
        Complete pw.x SCF input using the verified SSSP recommendation.
    """
    if element not in REFERENCE_LATTICES:
        raise ValueError(f"unsupported EOS reference element: {element}")
    verified = verify_sssp([element])
    if not verified["valid"]:
        raise RuntimeError(f"SSSP verification failed: {verified['errors']}")
    pseudo = verified["records"][element]["filename"]
    celldm = lattice_A * BOHR_PER_ANGSTROM
    return f"""&CONTROL
 calculation='scf', prefix='{prefix}', pseudo_dir='{SSSP_DIR}', outdir='./tmp',
 disk_io='none', verbosity='high'
/
&SYSTEM
 ibrav={IBRAV[element]}, celldm(1)={celldm:.12f}, nat=1, ntyp=1,
 ecutwfc={verified['ecutwfc_Ry']}, ecutrho={verified['ecutrho_Ry']},
 occupations='smearing', smearing='mv', degauss=0.02, nspin=2,
 starting_magnetization(1)=0.6
/
&ELECTRONS
 conv_thr=1.0d-10, mixing_beta=0.3, electron_maxstep=200
/
ATOMIC_SPECIES
{element} 1.0 {pseudo}
ATOMIC_POSITIONS crystal
{element} 0.0 0.0 0.0
K_POINTS automatic
12 12 12 1 1 1
"""


def parse_qe_total_energy(path: str | Path) -> float:
    """Extract the final converged total energy from a raw pw.x output.

    Args:
        path: Completed Quantum ESPRESSO text output.

    Returns:
        Final total energy in electronvolts.
    """
    text = Path(path).read_text(errors="replace")
    if "JOB DONE." not in text:
        raise RuntimeError(f"Quantum ESPRESSO output did not complete: {path}")
    energies = re.findall(r"!\s+total energy\s+=\s+([-+0-9.Ee]+)\s+Ry", text)
    if not energies:
        raise RuntimeError(f"no converged total energy found: {path}")
    return float(energies[-1]) * RY_TO_EV


def fit_equilibrium_lattice(lattices_A: list[float], energies_eV: list[float]) -> float:
    """Fit a local quadratic EOS minimum after checking it is interpolated.

    Args:
        lattices_A: Ordered conventional lattice parameters in ångströms.
        energies_eV: Corresponding converged primitive-cell energies.

    Returns:
        Interpolated equilibrium lattice parameter in ångströms.
    """
    if len(lattices_A) != len(energies_eV) or len(lattices_A) < 5:
        raise ValueError("EOS fit requires at least five paired points")
    coefficients = np.polyfit(np.asarray(lattices_A), np.asarray(energies_eV), 2)
    if coefficients[0] <= 0:
        raise RuntimeError("EOS fit has no convex energy minimum")
    minimum = float(-coefficients[1] / (2.0 * coefficients[0]))
    if not min(lattices_A) < minimum < max(lattices_A):
        raise RuntimeError("EOS minimum is outside the sampled interval")
    return minimum


def run_reference_eos(element: str, output_dir: str | Path) -> dict:
    """Execute and preserve a seven-point QE equation-of-state benchmark.

    Args:
        element: Fe or Ni reference material.
        output_dir: Directory receiving raw inputs, outputs, and summary JSON.

    Returns:
        Checksum-bound summary containing the independently fitted observable.
    """
    root = Path(output_dir).expanduser().resolve() / element.lower()
    root.mkdir(parents=True, exist_ok=True)
    pw = resolve_qe_executable("pw.x")
    execution = QEExecutionConfig.production_default()
    center = REFERENCE_LATTICES[element]
    lattices = [center * factor for factor in (0.96, 0.9733333333, 0.9866666667,
                                                1.0, 1.0133333333, 1.0266666667,
                                                1.04)]
    energies: list[float] = []
    artifacts = []
    for index, lattice in enumerate(lattices):
        input_path = root / f"eos_{index:02d}.in"
        output_path = root / f"eos_{index:02d}.out"
        input_path.write_text(qe_bulk_input(element, lattice, f"{element.lower()}_{index}"))
        command = build_qe_command(pw, str(input_path), execution)
        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = str(execution.omp_threads)
        environment.setdefault("OPENBLAS_NUM_THREADS", "1")
        with output_path.open("wb") as stdout:
            completed = subprocess.run(
                command, stdout=stdout, stderr=subprocess.STDOUT, env=environment,
                cwd=root, timeout=1800, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"pw.x failed for {element} point {index}")
        energy = parse_qe_total_energy(output_path)
        energies.append(energy)
        artifacts.append({
            "lattice_A": lattice,
            "energy_eV": energy,
            "input": input_path.name,
            "input_sha256": _digest(input_path),
            "output": output_path.name,
            "output_sha256": _digest(output_path),
            "command": command,
        })
    equilibrium = fit_equilibrium_lattice(lattices, energies)
    summary = {
        "schema_version": 1,
        "solver": "Quantum ESPRESSO pw.x",
        "solver_family": "quantum_espresso",
        "element": element,
        "structure": STRUCTURES[element],
        "protocol": "bulk-eos-zero-pressure-v1",
        "equilibrium_lattice_parameter_A": equilibrium,
        "points": artifacts,
        "converged": True,
        "execution": {
            "mpi_ranks": execution.mpi_ranks,
            "omp_threads": execution.omp_threads,
            "kpoint_pools": execution.kpoint_pools,
        },
    }
    summary_path = root / "eos_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return {**summary, "summary_path": str(summary_path),
            "summary_sha256": _digest(summary_path)}


def collect_reference_observations(output_dir: str | Path) -> Path:
    """Build validator input from previously completed QE EOS summaries.

    Args:
        output_dir: Directory containing the Fe and Ni EOS subdirectories.

    Returns:
        Path to the checksum-bound observation manifest.
    """
    root = Path(output_dir).expanduser().resolve()
    observations = []
    for element in sorted(REFERENCE_LATTICES):
        summary_path = root / element.lower() / "eos_summary.json"
        summary = json.loads(summary_path.read_text())
        if summary.get("converged") is not True or len(summary.get("points", [])) != 7:
            raise RuntimeError(f"incomplete EOS summary: {summary_path}")
        observations.append({
            "reference_id": REFERENCE_IDS[element],
            "material_id": f"{element}/{STRUCTURES[element]}/bulk",
            "observable": "equilibrium_lattice_parameter",
            "value": summary["equilibrium_lattice_parameter_A"],
            "unit": "angstrom",
            "solver_family": "quantum_espresso",
            "protocol": "bulk-eos-zero-pressure-v1",
            "conditions": {
                "temperature_K": 293,
                "crystal_structure": STRUCTURES[element],
            },
            "calculation_conditions": {
                "model_temperature_K": 0,
                "electronic_smearing": "Marzari-Vanderbilt",
                "electronic_smearing_Ry": 0.02,
                "kpoint_mesh": [12, 12, 12],
                "spin_polarized": True,
                "exchange_correlation": "PBE",
            },
            "converged": True,
            "mock": False,
            "candidate_specific": True,
            "benchmarked": False,
            "artifact_path": str(summary_path.relative_to(root)),
            "artifact_sha256": _digest(summary_path),
        })
    manifest = root / "reference_observations.json"
    manifest.write_text(json.dumps({
        "schema_version": 1,
        "observations": observations,
    }, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    """Run selected Fe/Ni reference calculations from the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--element", action="append", choices=sorted(REFERENCE_LATTICES))
    parser.add_argument("--output-dir", default="results/validation/qe_reference_eos")
    parser.add_argument("--collect-existing", action="store_true",
                        help="only collect completed summaries; run no calculations")
    arguments = parser.parse_args()
    if arguments.collect_existing:
        print(collect_reference_observations(arguments.output_dir))
        return
    for element in arguments.element or sorted(REFERENCE_LATTICES):
        result = run_reference_eos(element, arguments.output_dir)
        print(json.dumps({
            "element": element,
            "equilibrium_lattice_parameter_A": result["equilibrium_lattice_parameter_A"],
            "summary_path": result["summary_path"],
        }))
    print(collect_reference_observations(arguments.output_dir))


if __name__ == "__main__":
    main()
