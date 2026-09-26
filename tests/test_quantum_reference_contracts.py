#!/usr/bin/env python3
"""Portable contracts for honest solver-to-literature validation."""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.data_models.reference_validation import ComputedObservation
from pipeline.validation.reference_benchmarks import (
    compare_with_reference,
    computed_observation,
    load_reference_manifest,
    sha256_file,
)
from pipeline.validation.qe_reference_eos import (
    fit_equilibrium_lattice,
    parse_qe_total_energy,
    qe_bulk_input,
)
from pipeline.validation.ni111_methane_reference import (
    REFERENCE_BARRIER_EV,
    build_ni111_methane_endpoints,
    prepare_ni111_methane_benchmark,
)


MANIFEST = Path(__file__).parent / "fixtures" / "quantum_reference_manifest.json"
FAKE_DIGEST = "a" * 64


def _observation(reference, value=None, **updates):
    values = dict(
        reference_id=reference.reference_id,
        material_id=reference.material_id,
        observable=reference.observable,
        value=reference.value if value is None else value,
        unit=reference.unit,
        solver_family=reference.solver_family,
        protocol=reference.protocol,
        conditions=reference.conditions,
        calculation_conditions={"model_temperature_K": 0},
        converged=True,
        mock=False,
        candidate_specific=True,
        benchmarked=True,
        artifact_sha256=FAKE_DIGEST,
    )
    values.update(updates)
    return ComputedObservation(**values)


def _raises(fragment, function, *args):
    try:
        function(*args)
    except ValueError as error:
        assert fragment in str(error), str(error)
    else:
        raise AssertionError(f"expected ValueError containing {fragment!r}")


def _raises_runtime(fragment, function, *args):
    try:
        function(*args)
    except RuntimeError as error:
        assert fragment in str(error), str(error)
    else:
        raise AssertionError(f"expected RuntimeError containing {fragment!r}")


def test_manifest_preserves_reference_kind_and_conditions():
    references = load_reference_manifest(MANIFEST)
    assert set(references) == {
        "nist-fe-bcc-lattice-293k", "nist-ni-fcc-lattice-293k",
        "bengaard-ni111-ch4-barrier",
    }
    assert references["nist-fe-bcc-lattice-293k"].reference_kind == "experiment"
    assert references["bengaard-ni111-ch4-barrier"].reference_kind == "published_computation"
    assert references["nist-ni-fcc-lattice-293k"].conditions["temperature_K"] == 293


def test_qe_observable_passes_and_fails_declared_tolerance():
    reference = load_reference_manifest(MANIFEST)["nist-ni-fcc-lattice-293k"]
    assert compare_with_reference(reference, _observation(reference, 3.55)).passed
    failure = compare_with_reference(reference, _observation(reference, 3.70))
    assert not failure.passed
    assert failure.absolute_error > failure.tolerance


def test_comparison_rejects_wrong_material_observable_units_and_protocol():
    reference = load_reference_manifest(MANIFEST)["nist-fe-bcc-lattice-293k"]
    for field, value in (
        ("material_id", "Ni/fcc/bulk"),
        ("observable", "total_energy"),
        ("unit", "Ry"),
        ("protocol", "single-point-energy-v1"),
    ):
        _raises(field, compare_with_reference, reference,
                replace(_observation(reference), **{field: value}))
    _raises("conditions", compare_with_reference, reference,
            replace(_observation(reference), conditions={"temperature_K": 0}))


def test_incomplete_mock_and_toy_calculations_fail_closed():
    reference = load_reference_manifest(MANIFEST)["bengaard-ni111-ch4-barrier"]
    _raises("unconverged", compare_with_reference, reference,
            replace(_observation(reference), converged=False))
    _raises("mock", compare_with_reference, reference,
            replace(_observation(reference), mock=True))
    _raises("toy or generic", compare_with_reference, reference,
            replace(_observation(reference), candidate_specific=False))


def test_vqe_numerical_success_is_not_physical_validation():
    """Exact agreement on a toy Hamiltonian cannot be compared to experiment."""
    reference = replace(
        load_reference_manifest(MANIFEST)["nist-ni-fcc-lattice-293k"],
        solver_family="cudaq_vqe",
    )
    toy = _observation(
        reference, 3.5238, candidate_specific=False, benchmarked=True)
    _raises("toy or generic", compare_with_reference, reference, toy)
    physical = replace(toy, candidate_specific=True, benchmarked=False)
    _raises("exact-solver", compare_with_reference, reference, physical)


def test_artifact_records_require_traceable_sha256():
    record = {
        "reference_id": "nist-fe-bcc-lattice-293k",
        "material_id": "Fe/bcc/bulk",
        "observable": "equilibrium_lattice_parameter",
        "value": 2.87,
        "unit": "angstrom",
        "solver_family": "quantum_espresso",
        "protocol": "bulk-eos-zero-pressure-v1",
        "conditions": {"temperature_K": 293, "crystal_structure": "bcc"},
        "calculation_conditions": {"model_temperature_K": 0},
        "converged": True,
        "mock": False,
        "candidate_specific": True,
        "benchmarked": False,
        "artifact_sha256": "not-a-digest",
    }
    _raises("SHA-256", computed_observation, record)
    with tempfile.TemporaryDirectory() as temporary:
        artifact = Path(temporary) / "qe.out"
        artifact.write_text("JOB DONE.\n")
        record["artifact_sha256"] = sha256_file(artifact)
        assert computed_observation(record).artifact_sha256 == sha256_file(artifact)


def test_manifest_schema_is_fail_closed():
    content = json.loads(MANIFEST.read_text())
    del content["references"][0]["source_url"]
    with tempfile.TemporaryDirectory() as temporary:
        invalid = Path(temporary) / "invalid.json"
        invalid.write_text(json.dumps(content))
        _raises("source_url", load_reference_manifest, invalid)


def test_qe_reference_runner_parses_raw_output_and_fits_eos():
    lattices = [2.7, 2.8, 2.9, 3.0, 3.1]
    energies = [(lattice - 2.91) ** 2 - 8.0 for lattice in lattices]
    assert abs(fit_equilibrium_lattice(lattices, energies) - 2.91) < 1e-10
    with tempfile.TemporaryDirectory() as temporary:
        output = Path(temporary) / "pw.out"
        output.write_text("!    total energy = -10.000000 Ry\nJOB DONE.\n")
        assert abs(parse_qe_total_energy(output) + 136.05693122994) < 1e-10
        output.write_text("!    total energy = -10.000000 Ry\n")
        _raises_runtime("did not complete", parse_qe_total_energy, output)


def test_qe_reference_inputs_preserve_structure_and_convergence_protocol():
    for element, ibrav in (("Fe", 3), ("Ni", 2)):
        text = qe_bulk_input(element, 3.0, f"{element.lower()}_contract")
        assert f"ibrav={ibrav}" in text
        assert "12 12 12 1 1 1" in text
        assert "conv_thr=1.0d-10" in text
        assert "nspin=2" in text


def test_ni111_methane_reference_geometry_is_matched_and_constrained():
    initial, final = build_ni111_methane_endpoints()
    assert len(initial) == len(final) == 17
    assert initial.get_chemical_symbols() == final.get_chemical_symbols()
    assert initial.get_chemical_symbols()[-5:] == ["C", "H", "H", "H", "H"]
    assert np.allclose(initial.cell.array, final.cell.array)
    assert initial.pbc.all() and final.pbc.all()
    initial_fixed = set(initial.constraints[0].get_indices())
    final_fixed = set(final.constraints[0].get_indices())
    assert initial_fixed == final_fixed and len(initial_fixed) == 8
    assert not initial_fixed.intersection(range(12, 17))
    initial_bond = np.linalg.norm(initial.positions[-1] - initial.positions[-5])
    final_bond = np.linalg.norm(final.positions[-1] - final.positions[-5])
    assert abs(initial_bond - 1.09) < 1e-10
    assert final_bond > 1.5


def test_ni111_benchmark_is_isolated_and_prepares_real_qe_inputs():
    with tempfile.TemporaryDirectory() as temporary:
        result = prepare_ni111_methane_benchmark(temporary)
        manifest = json.loads(Path(result["manifest"]).read_text())
        assert manifest["benchmark_only"] is True
        assert manifest["candidate_selection_authority"] is False
        assert manifest["reference_barrier_eV"] == REFERENCE_BARRIER_EV
        assert manifest["reference_kind"] == "published_computation"
        assert len(manifest["fixed_atom_indices"]) == 8
        for name in ("initial.relax.in", "final.relax.in"):
            text = (Path(temporary) / name).read_text()
            assert "calculation='relax'" in text
            assert sum(line.startswith("Ni ") and line.endswith(" 0 0 0")
                       for line in text.splitlines()) == 8
            assert "forc_conv_thr=1.0d-3" in text


def main():
    tests = (
        test_manifest_preserves_reference_kind_and_conditions,
        test_qe_observable_passes_and_fails_declared_tolerance,
        test_comparison_rejects_wrong_material_observable_units_and_protocol,
        test_incomplete_mock_and_toy_calculations_fail_closed,
        test_vqe_numerical_success_is_not_physical_validation,
        test_artifact_records_require_traceable_sha256,
        test_manifest_schema_is_fail_closed,
        test_qe_reference_runner_parses_raw_output_and_fits_eos,
        test_qe_reference_inputs_preserve_structure_and_convergence_protocol,
        test_ni111_methane_reference_geometry_is_matched_and_constrained,
        test_ni111_benchmark_is_isolated_and_prepares_real_qe_inputs,
    )
    for test in tests:
        test()
        print("PASS", test.__name__)
    print(f"{len(tests)}/{len(tests)} quantum-reference contracts passed")


if __name__ == "__main__":
    main()
