#!/usr/bin/env python3
"""Focused contracts for PFR, MMBCR, carbon, kinetics, and provenance merges."""

import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.process.pathway_modes import (
    reactor_applicability, reactor_types_for_mode, validate_mode_reactors)
from pipeline.process.reactor_mechanisms import CandidateKinetics
from pipeline.process.reactor_models import ReactorConfig, _kinetics_evidence


def _raises(message, exception_type=Exception):
    class Raises:
        def __enter__(self):
            return self

        def __exit__(self, kind, value, traceback):
            assert kind is not None and issubclass(kind, exception_type)
            assert message in str(value), (message, str(value))
            return True
    return Raises()


def _validated_kinetics(candidate_id="merge-candidate"):
    row = {
        "E_act": 0.91, "dE_H": -0.20, "dE_CH3": -0.35, "dE_C": -0.60,
        "screening_protocol": "merge-contract-v1",
    }
    resolved = {
        "methane_activation_eV": 0.80,
        "ch3_dehydrogenation_eV": 0.70,
        "ch2_dehydrogenation_eV": 0.75,
        "ch_dehydrogenation_eV": 0.65,
        "h2_desorption_eV": 0.55,
        "carbon_transfer_eV": 0.45,
    }
    validation = {
        "candidate_id": candidate_id, "complete": True,
        "evidence_level": "converged_dft_neb_frequency",
        "resolved_kinetics_eV": resolved,
    }
    return CandidateKinetics.from_screening_row(
        row, candidate_id=candidate_id, validation=validation)


def test_reactor_modes_cannot_cross_wire_pfr_and_mmbcr():
    assert reactor_types_for_mode("thermocatalytic_pfr") == ("PFR",)
    assert reactor_types_for_mode("mmbcr") == ("MMBCR",)
    validate_mode_reactors("thermocatalytic_pfr", ["PFR"])
    validate_mode_reactors("mmbcr", ["MMBCR"])
    with _raises("requires reactors", ValueError):
        validate_mode_reactors("mmbcr", ["PFR"])
    assert reactor_applicability("PFR", "SAC") == (True, None)
    assert reactor_applicability("MMBCR", "MoltenMetal") == (True, None)
    assert reactor_applicability("MMBCR", "SAC")[0] is False


def test_candidate_specific_kinetics_identity_and_provenance_are_complete():
    kinetics = _validated_kinetics()
    resolved = kinetics.resolved()
    assert resolved["candidate_id"] == "merge-candidate"
    assert resolved["quantitative_status"] == "candidate_specific"
    barrier_names = {
        "methane_activation_eV", "ch3_dehydrogenation_eV",
        "ch2_dehydrogenation_eV", "ch_dehydrogenation_eV",
        "h2_desorption_eV", "carbon_transfer_eV",
    }
    assert all(
        resolved["provenance"][name] ==
        "candidate_specific:converged_dft_neb_frequency"
        for name in barrier_names)
    with _raises("candidate_id mismatch", ValueError):
        CandidateKinetics.from_screening_row(
            {"E_act": 0.8}, candidate_id="wrong",
            validation={"candidate_id": "other", "complete": True,
                        "evidence_level": "converged_dft_neb_frequency",
                        "resolved_kinetics_eV": {}})


def test_generated_mechanism_is_balanced_and_matches_carbon_provenance():
    import cantera as ct
    import pipeline.process.reactor_mechanisms as mechanisms

    with tempfile.TemporaryDirectory() as tmp, patch.object(
            mechanisms, "MECHANISMS_DIR", Path(tmp)):
        path = mechanisms.write_full_mechanism(
            "merge_candidate", kinetics=_validated_kinetics())
        metadata = json.loads(path.with_suffix(".kinetics.json").read_text())
        gas = ct.Solution(str(path), "gas")
        surface = ct.Interface(str(path), "merge_candidate_surface", [gas])

        assert metadata["catalyst_name"] == "merge_candidate"
        assert metadata["inputs"]["candidate_id"] == "merge-candidate"
        assert metadata["inputs"]["quantitative_status"] == "candidate_specific"
        # Cantera rejects elementally unbalanced reactions while loading. Check
        # both reaction collections are nonempty so that validation was real.
        assert gas.n_reactions > 0
        assert surface.n_reactions > 0

        carbon_is_gas = "C_graphite" in gas.species_names
        declared = metadata["carbon_phase_model"]
        if carbon_is_gas:
            assert declared == "legacy_gas_tracer"
        else:
            assert declared != "legacy_gas_tracer"

        evidence = _kinetics_evidence(ReactorConfig(
            catalyst_name="merge_candidate", mechanism_file=str(path)))
        assert evidence["kinetics_status"] == "candidate_specific"
        assert evidence["carbon_phase_model"] == declared
        assert evidence["can_exclude_candidate"] is (not carbon_is_gas)


TESTS = (
    test_reactor_modes_cannot_cross_wire_pfr_and_mmbcr,
    test_candidate_specific_kinetics_identity_and_provenance_are_complete,
    test_generated_mechanism_is_balanced_and_matches_carbon_provenance,
)


def main():
    for test in TESTS:
        test()
        print("PASS", test.__name__)
    print(f"{len(TESTS)}/{len(TESTS)} reactor merge contracts passed")


if __name__ == "__main__":
    main()
