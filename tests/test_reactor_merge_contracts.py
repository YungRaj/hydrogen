#!/usr/bin/env python3
"""Focused contracts for PFR, MMBCR, carbon, kinetics, and provenance merges."""

import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.reactors.modes import (
    reactor_applicability,
    reactor_types_for_mode,
    validate_mode_reactors,
)
from pipeline.reactors.mechanisms import CandidateKinetics
from pipeline.reactors.models import ReactorConfig, _kinetics_evidence


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
        "E_act": 0.91,
        "dE_H": -0.20,
        "dE_CH3": -0.35,
        "dE_C": -0.60,
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
        "candidate_id": candidate_id,
        "complete": True,
        "evidence_level": "converged_dft_neb_frequency",
        "resolved_kinetics_eV": resolved,
    }
    return CandidateKinetics.from_screening_row(
        row, candidate_id=candidate_id, validation=validation
    )


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
        "methane_activation_eV",
        "ch3_dehydrogenation_eV",
        "ch2_dehydrogenation_eV",
        "ch_dehydrogenation_eV",
        "h2_desorption_eV",
        "carbon_transfer_eV",
    }
    assert all(
        resolved["provenance"][name] == "candidate_specific:converged_dft_neb_frequency"
        for name in barrier_names
    )
    with _raises("candidate_id mismatch", ValueError):
        CandidateKinetics.from_screening_row(
            {"E_act": 0.8},
            candidate_id="wrong",
            validation={
                "candidate_id": "other",
                "complete": True,
                "evidence_level": "converged_dft_neb_frequency",
                "resolved_kinetics_eV": {},
            },
        )


def test_generated_mechanism_is_balanced_and_matches_carbon_provenance():
    import cantera as ct
    import pipeline.reactors.mechanisms as mechanisms

    with (
        tempfile.TemporaryDirectory() as tmp,
        patch.object(mechanisms, "MECHANISMS_DIR", Path(tmp)),
    ):
        path = mechanisms.write_full_mechanism(
            "merge_candidate", kinetics=_validated_kinetics()
        )
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

        evidence = _kinetics_evidence(
            ReactorConfig(catalyst_name="merge_candidate", mechanism_file=str(path))
        )
        assert evidence["kinetics_status"] == "candidate_specific"
        assert evidence["carbon_phase_model"] == declared
        assert evidence["can_exclude_candidate"] is (not carbon_is_gas)


def test_kinetics_fidelity_ladder_bep_mlip_neb_and_dft():
    bep_row = {
        "E_act": 0.95,
        "dE_H": -0.25,
        "dE_CH3": -0.40,
        "dE_C": -0.70,
        "screening_protocol": "ladder-v1",
    }
    bep_kin = CandidateKinetics.from_screening_row(
        bep_row, candidate_id="cand-ladder"
    ).resolved()
    assert bep_kin["quantitative_status"] == "screening_complete"
    assert bep_kin["provenance"]["methane_activation_eV"] == "screening:ladder-v1"
    assert bep_kin["provenance"]["ch2_dehydrogenation_eV"] == "template_default"

    neb_row = dict(bep_row, E_act_neb=0.84, ch3_dehydrogenation_eV=0.72)
    neb_row_kin = CandidateKinetics.from_screening_row(
        neb_row, candidate_id="cand-ladder"
    ).resolved()
    assert neb_row_kin["methane_activation_eV"] == 0.84
    assert (
        neb_row_kin["provenance"]["methane_activation_eV"]
        == "mlip_neb:ladder-v1:E_act_neb"
    )
    assert neb_row_kin["ch3_dehydrogenation_eV"] == 0.72
    assert (
        neb_row_kin["provenance"]["ch3_dehydrogenation_eV"]
        == "mlip_neb:ladder-v1:ch3_dehydrogenation_eV"
    )

    mlip_validation = {
        "candidate_id": "cand-ladder",
        "complete": True,
        "evidence_level": "converged_mlip_neb",
        "resolved_kinetics_eV": {
            "methane_activation_eV": 0.82,
            "ch3_dehydrogenation_eV": 0.71,
            "ch2_dehydrogenation_eV": 0.76,
            "ch_dehydrogenation_eV": 0.66,
            "h2_desorption_eV": 0.56,
            "carbon_transfer_eV": 0.46,
        },
    }
    mlip_kin = CandidateKinetics.from_screening_row(
        bep_row, candidate_id="cand-ladder", validation=mlip_validation
    ).resolved()
    assert mlip_kin["quantitative_status"] == "mlip_neb_complete"
    assert (
        mlip_kin["provenance"]["ch2_dehydrogenation_eV"]
        == "mlip_neb:converged_mlip_neb"
    )


TESTS = (
    test_reactor_modes_cannot_cross_wire_pfr_and_mmbcr,
    test_candidate_specific_kinetics_identity_and_provenance_are_complete,
    test_generated_mechanism_is_balanced_and_matches_carbon_provenance,
    test_kinetics_fidelity_ladder_bep_mlip_neb_and_dft,
)


def main():
    for test in TESTS:
        test()
        print("PASS", test.__name__)
    print(f"{len(TESTS)}/{len(TESTS)} reactor merge contracts passed")


if __name__ == "__main__":
    main()
