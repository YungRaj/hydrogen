#!/usr/bin/env python3
"""Prepare the literature Ni(111) first-C-H-cleavage QE/NEB benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import write as ase_write

from pipeline.validation.production_workflow import prepare_methane_neb
from pipeline.validation.qe_workflows import verify_sssp


REFERENCE_ID = "bengaard-ni111-ch4-barrier"
REFERENCE_BARRIER_EV = 1.05
REFERENCE_DOI = "https://doi.org/10.1006/jcat.2002.3579"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def build_ni111_methane_endpoints(
    lattice_A: float = 3.5238,
) -> tuple[Atoms, Atoms]:
    """Construct matched molecular and dissociated Ni(111) endpoint guesses.

    The structures are reproducible starting geometries, not converged evidence.
    A 2x2 three-layer slab is used, the lower two layers are fixed, and the top
    layer plus adsorbates are relaxed by the production workflow.

    Args:
        lattice_A: Conventional fcc Ni lattice parameter in ångströms.

    Returns:
        Initial CH4*/Ni(111) and final CH3*+H*/Ni(111) ASE structures with
        identical atom ordering, cells, periodicity, and constraints.
    """
    if lattice_A <= 0:
        raise ValueError("Ni lattice parameter must be positive")
    slab = fcc111("Ni", size=(2, 2, 3), a=lattice_A, vacuum=10.0,
                  orthogonal=False, periodic=True)
    z_values = np.asarray(slab.positions[:, 2])
    unique_layers = sorted({round(float(value), 6) for value in z_values})
    fixed = [index for index, value in enumerate(z_values)
             if round(float(value), 6) in set(unique_layers[:2])]
    slab.set_constraint(FixAtoms(indices=fixed))

    top_indices = [index for index, value in enumerate(z_values)
                   if round(float(value), 6) == unique_layers[-1]]
    anchor = slab.positions[top_indices[1]].copy()
    top_z = float(anchor[2])
    bond = 1.09
    radial = bond * np.sqrt(8.0 / 9.0)
    upper_z = bond / 3.0
    initial_vectors = [
        np.array([radial * np.cos(angle), radial * np.sin(angle), upper_z])
        for angle in (0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0)
    ] + [np.array([0.0, 0.0, -bond])]
    # The literature barrier is referenced to gas-phase methane and two clean
    # surface sites.  Plain PBE does not necessarily produce a bound molecular
    # CH4 precursor on Ni(111), so placing CH4 in the vacuum region avoids
    # inventing an adsorbed minimum that drifts during endpoint relaxation.
    initial_c = np.array([anchor[0], anchor[1], top_z + 5.5])
    initial_adsorbate = Atoms(
        "CH4", positions=[initial_c] +
        [(initial_c + vector).tolist() for vector in initial_vectors])

    final_c = np.array([anchor[0], anchor[1], top_z + 2.0])
    ch_bond = 1.08
    ch_radial = ch_bond * np.sqrt(8.0 / 9.0)
    ch_upper_z = ch_bond / 3.0
    final_vectors = [
        np.array([ch_radial * np.cos(angle), ch_radial * np.sin(angle), ch_upper_z])
        for angle in (0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0)
    ]
    dissociated_h = np.array([
        anchor[0] + lattice_A / np.sqrt(2.0) * 0.60,
        anchor[1],
        top_z + 1.05,
    ])
    final_adsorbate = Atoms(
        "CH4", positions=[final_c] +
        [(final_c + vector).tolist() for vector in final_vectors] +
        [dissociated_h.tolist()])

    initial = slab.copy() + initial_adsorbate
    final = slab.copy() + final_adsorbate
    initial.set_constraint(FixAtoms(indices=fixed))
    final.set_constraint(FixAtoms(indices=fixed))
    return initial, final


def prepare_ni111_methane_benchmark(output_dir: str | Path) -> dict:
    """Persist endpoint sources and prepare production QE relaxation inputs.

    Args:
        output_dir: Isolated benchmark directory receiving structures and QE
            inputs; it is not part of candidate discovery state.

    Returns:
        Preparation status and checksum-bound literature provenance.
    """
    verified = verify_sssp(["Ni", "C", "H"])
    if not verified["valid"]:
        raise RuntimeError(f"SSSP verification failed: {verified['errors']}")
    root = Path(output_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    initial, final = build_ni111_methane_endpoints()
    initial_path = root / "initial.traj"
    final_path = root / "final.traj"
    ase_write(initial_path, initial)
    ase_write(final_path, final)
    status = prepare_methane_neb(root, "ni111_ch4", initial, final)
    manifest = {
        "schema_version": 1,
        "benchmark_only": True,
        "candidate_selection_authority": False,
        "reference_id": REFERENCE_ID,
        "reference_barrier_eV": REFERENCE_BARRIER_EV,
        "reference_kind": "published_computation",
        "reference_doi": REFERENCE_DOI,
        "reaction": "CH4(g) + 2* -> CH3* + H* on Ni(111)",
        "surface_cell": "2x2 Ni(111), three layers, lower two fixed",
        "initial_structure": initial_path.name,
        "initial_sha256": _sha256(initial_path),
        "final_structure": final_path.name,
        "final_sha256": _sha256(final_path),
        "atom_count": len(initial),
        "fixed_atom_indices": [
            int(index) for index in initial.constraints[0].get_indices()],
        "sssp": verified,
        "acceptance": {
            "absolute_tolerance_eV": 0.15,
            "requires_relaxed_endpoints": True,
            "requires_converged_neb": True,
            "requires_one_imaginary_mode": True,
        },
        "limitations": [
            "finite 2x2 three-layer slab; production convergence must test slab and cell size",
            "gas-reference and adsorbed-product endpoints require QE relaxation and chemical inspection",
            "published target is a DFT value, not a direct experimental barrier",
        ],
    }
    manifest_path = root / "benchmark_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {
        "manifest": str(manifest_path),
        "manifest_sha256": _sha256(manifest_path),
        "status": status,
    }


def main() -> None:
    """Prepare the isolated Ni(111) reference benchmark from the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", default="results/validation/ni111_methane_reference")
    arguments = parser.parse_args()
    print(json.dumps(
        prepare_ni111_methane_benchmark(arguments.output_dir), indent=2))


if __name__ == "__main__":
    main()
