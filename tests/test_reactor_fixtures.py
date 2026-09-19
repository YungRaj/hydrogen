#!/usr/bin/env python3
"""End-to-end representative fixtures for every reactor execution path."""

from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.process.multiphysics_contract import load_validated_artifact
from pipeline.process.multiphysics_runner import run_backend
from pipeline.process.reactor_models import ReactorConfig, simulate_reactor
from tests.fixtures.reactor_cases import FixtureSolver, MODES, write_case


def _run_external_fixture(root: Path, reactor_type: str,
                          temperature_K: float) -> tuple[dict, dict]:
    """Run one external fixture through production artifact and reactor APIs."""
    candidate = f"fixture-{reactor_type.lower()}"
    mode = MODES[reactor_type]
    case, results = root / f"case-{reactor_type}", root / "artifacts"
    case.mkdir()
    write_case(case, reactor_type, candidate, temperature_K)
    model = case / "fixture-model.py"
    model.write_text("# deterministic test fixture; not a physical solver\n")
    target = run_backend(
        mode=mode, reactor_type=reactor_type, candidate_id=candidate,
        temperature_K=temperature_K, case_dir=case, results_dir=results,
        model_source="test-fixture:representative-reactor-v1",
        fenics_model=model if reactor_type in {"NTEC", "Electrochemical"} else None,
        execution=FixtureSolver(reactor_type, candidate, temperature_K).services())
    artifact = json.loads(target.read_text())
    loaded = load_validated_artifact(
        results, candidate, mode, reactor_type, temperature_K)
    assert loaded["valid"] is True
    mechanism = ""
    catalyst_name = candidate
    if reactor_type in {"Fluidized", "MMBCR"}:
        from pipeline.process.reactor_mechanisms import (
            CandidateKinetics, write_full_mechanism)
        catalyst_name = candidate.replace("-", "_")
        kinetics = CandidateKinetics.from_screening_row({
            "E_act": .72, "dE_H": -.3, "dE_CH3": -.55, "dE_C": -1.1,
            "screening_protocol": "test-fixture:screening-v1"},
            candidate_id=candidate)
        with patch("pipeline.process.reactor_mechanisms.MECHANISMS_DIR",
                   root / "mechanisms"):
            mechanism = str(write_full_mechanism(catalyst_name, kinetics=kinetics))
    with patch("pipeline.process.reactor_models.save_json"):
        result = simulate_reactor(ReactorConfig(
            reactor_type=reactor_type, pathway_mode=mode, candidate_id=candidate,
            catalyst_name=catalyst_name, mechanism_file=mechanism,
            T_inlet_K=temperature_K,
            material_class=("MoltenMetal" if reactor_type == "MMBCR" else
                            "SAC" if reactor_type == "Fluidized" else None),
            multiphysics_results_dir=str(results)))
    return artifact, result


def test_representative_external_reactor_paths() -> None:
    """Validate handoffs and numerical consumption for four external modes."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        matrix = {}
        for reactor_type, temperature_K in (
                ("Fluidized", 900.), ("MMBCR", 1000.),
                ("NTEC", 300.), ("Electrochemical", 300.)):
            artifact, result = _run_external_fixture(
                root, reactor_type, temperature_K)
            assert artifact["complete"] is True
            assert artifact["provenance"]["model_source"].startswith("test-fixture:")
            assert artifact["convergence"]["mesh_independent"] is True
            assert artifact["convergence"]["conservation_satisfied"] is True
            assert artifact["model_validation"]["passed"] is True
            assert result["status"] == "complete"
            if reactor_type in {"NTEC", "Electrochemical"}:
                assert result["valid"] is True
            assert result["can_exclude_candidate"] is False
            assert 0 <= result["CH4_conversion"] <= 1
            matrix[reactor_type] = result["CH4_conversion"]
        assert set(matrix) == {"Fluidized", "MMBCR", "NTEC", "Electrochemical"}


def test_representative_cantera_pfr_path() -> None:
    """Run the real Cantera PFR entry point with candidate screening kinetics."""
    from pipeline.process.reactor_mechanisms import CandidateKinetics, write_full_mechanism

    candidate = "fixture_pfr"
    kinetics = CandidateKinetics.from_screening_row({
        "E_act": .72, "dE_H": -.3, "dE_CH3": -.55, "dE_C": -1.1,
        "screening_protocol": "test-fixture:screening-v1"}, candidate_id=candidate)
    with tempfile.TemporaryDirectory() as tmp, \
            patch("pipeline.process.reactor_mechanisms.MECHANISMS_DIR", Path(tmp)), \
            patch("pipeline.process.reactor_models.save_json"):
        mechanism = write_full_mechanism(candidate, kinetics=kinetics)
        result = simulate_reactor(ReactorConfig(
            reactor_type="PFR", pathway_mode="thermocatalytic_pfr",
            candidate_id=candidate, catalyst_name=candidate,
            material_class="SAC", mechanism_file=str(mechanism),
            T_inlet_K=900., max_residence_time_s=.05))
    assert result["status"] == "complete" and result.get("mock", False) is False
    assert math.isfinite(result["CH4_conversion"])
    assert 0 <= result["CH4_conversion"] <= 1
    assert result["can_exclude_candidate"] is False
    assert "incomplete_candidate_kinetics" in result["reactor_evidence_limitations"]


def main() -> None:
    """Run this file without requiring pytest."""
    tests = [value for name, value in sorted(globals().items())
             if name.startswith("test_") and callable(value)]
    for test in tests:
        test()
        print("PASS", test.__name__)
    print(f"{len(tests)}/{len(tests)} representative reactor fixtures passed")


if __name__ == "__main__":
    main()
