"""Tests for multi-objective vectors, uncertainty LCB scoring, and Pareto selection."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

from pipeline.screening.objectives import (
    ORR_OBJECTIVES,
    TURQUOISE_OBJECTIVES,
    active_objective_indices,
    objective_specs_for_application,
    pareto_ranks,
    predict_objectives_with_uncertainty,
    resolve_objective_specs,
    robust_objective_matrix,
)
from pipeline.screening.small_data_ranker import (
    RankerDiagnostics,
    fit_tree_ranker,
    orr_tree_objectives,
    turquoise_tree_objectives,
)
from pipeline.screening.stage_selection import select_for_reactor
from pipeline.search.adaptive_validation import (
    AdvancementCriterion,
    advancement_scorecard,
    allocate_validation_batch,
    experimental_slate,
    verify_portfolio_artifact,
    write_portfolio_artifact,
)
from pipeline.search.branch_search import BranchConfig, run_branch_and_bound, verify_branch_coverage
from pipeline.search.indexed_space import deterministic_tree_probes


class TestMultiObjectiveUncertaintyAndPareto(unittest.TestCase):
    """Verify 4-column physical objective vectors, LCB scoring, and Pareto selection."""

    def _make_turquoise_frame(self, count: int = 24) -> tuple[list[tuple], pd.DataFrame]:
        genomes = deterministic_tree_probes(count)
        idx = np.arange(len(genomes), dtype=float)
        de_h = -2.4 - 0.08 * idx
        de_c = -4.2 - 0.15 * idx
        coking = de_c - 2.0 * de_h
        frame = pd.DataFrame({
            "genome": [repr(g) for g in genomes],
            "valid": True,
            "E_act": np.linspace(0.65, 1.85, len(genomes)),
            "E_act_censored": False,
            "pyrolysis_viable": True,
            "dE_H": de_h,
            "dE_C": de_c,
            "coking_index": coking,
        })
        return genomes, frame

    def _make_orr_frame(self, count: int = 24) -> tuple[list[tuple], pd.DataFrame]:
        genomes = deterministic_tree_probes(count)
        idx = np.arange(len(genomes), dtype=float)
        dg_oh = 0.55 + 0.04 * idx
        dg_o = 1.20 + 0.07 * idx
        frame = pd.DataFrame({
            "genome": [repr(g) for g in genomes],
            "valid": True,
            "orr_overpotential_V": np.linspace(0.32, 1.15, len(genomes)),
            "dG_OH_eV": dg_oh,
            "dG_O_eV": dg_o,
            "binding_strength": np.abs(dg_oh) + np.abs(dg_o),
        })
        return genomes, frame

    def test_turquoise_and_orr_objectives_have_four_nonconstant_columns(self):
        self.assertEqual(len(TURQUOISE_OBJECTIVES), 4)
        self.assertEqual(len(ORR_OBJECTIVES), 4)
        self.assertEqual(len(objective_specs_for_application("turquoise_hydrogen")), 4)
        self.assertEqual(len(objective_specs_for_application("fuel_cell_orr")), 4)

        turq_genomes, turq_df = self._make_turquoise_frame(24)
        turq_ranker = fit_tree_ranker(turq_df, "turquoise_hydrogen")
        self.assertIn("coking_susceptibility_eV", turq_ranker.secondary_models)
        self.assertIn("desorption_overbinding_penalty_eV", turq_ranker.secondary_models)

        turq_objs, turq_stds = turquoise_tree_objectives(
            turq_genomes, turq_ranker, return_uncertainty=True
        )
        self.assertEqual(turq_objs.shape, (len(turq_genomes), 4))
        self.assertEqual(turq_stds.shape, (len(turq_genomes), 4))
        for col in range(4):
            self.assertGreater(
                float(np.std(turq_objs[:, col])),
                0.0,
                f"Turquoise col {col} is constant",
            )
        self.assertGreater(float(np.mean(turq_stds[:, 0])), 0.0)

        orr_genomes, orr_df = self._make_orr_frame(24)
        orr_ranker = fit_tree_ranker(orr_df, "fuel_cell_orr")
        self.assertIn("oxide_overbinding_penalty_eV", orr_ranker.secondary_models)

        orr_objs, orr_stds = orr_tree_objectives(
            orr_genomes, orr_ranker, return_uncertainty=True
        )
        self.assertEqual(orr_objs.shape, (len(orr_genomes), 4))
        self.assertEqual(orr_stds.shape, (len(orr_genomes), 4))
        for col in range(4):
            self.assertGreater(
                float(np.std(orr_objs[:, col])),
                0.0,
                f"ORR col {col} is constant",
            )
        self.assertGreater(float(np.mean(orr_stds[:, 0])), 0.0)

    def test_unvalidated_ranker_marks_secondary_objectives_inactive(self):
        genomes = deterministic_tree_probes(8)
        diagnostics = RankerDiagnostics(
            sample_count=20,
            held_out_group_count=5,
            validation_strategy="held_out_material_class_cv",
            spearman_correlation=-0.1,
            model_mae=2.0,
            baseline_mae=1.0,
            ranking_validated=False,
        )
        quality = np.linspace(0.8, 0.1, len(genomes))
        ranker = SimpleNamespace(
            application="turquoise_hydrogen",
            diagnostics=diagnostics,
            predict=lambda items, uncertainty=True: (
                quality[: len(items)],
                np.full(len(items), 0.05) if uncertainty else np.zeros(len(items)),
            ),
        )
        specs = resolve_objective_specs("turquoise_hydrogen", ranker)
        self.assertEqual([s.active for s in specs], [True, False, False, False])
        self.assertEqual(active_objective_indices("turquoise_hydrogen", ranker), (0,))

        objs, stds = predict_objectives_with_uncertainty(genomes, ranker)
        for col in range(4):
            self.assertGreater(float(np.std(objs[:, col])), 0.0)
        ranks = pareto_ranks(
            objs,
            active_indices=active_objective_indices("turquoise_hydrogen", ranker),
        )
        primary_order = np.argsort(objs[:, 0], kind="mergesort")
        self.assertEqual(int(ranks[primary_order[0]]), 0)

    def test_robust_lcb_monotone_in_kappa(self):
        means = np.array([[0.8, -0.5, 0.2, 0.1], [1.1, -0.2, 0.0, 0.4]], dtype=float)
        stds = np.array([[0.15, 0.05, 0.02, 0.0], [0.05, 0.10, 0.01, 0.0]], dtype=float)
        m0 = robust_objective_matrix(means, stds, kappa=0.0)
        m1 = robust_objective_matrix(means, stds, kappa=0.5)
        m2 = robust_objective_matrix(means, stds, kappa=1.5)
        np.testing.assert_allclose(m0, means)
        self.assertTrue(np.all(m1 >= m0 - 1e-12))
        self.assertTrue(np.all(m2 >= m1 - 1e-12))

    def test_pareto_selection_and_validation_batch_preserve_class_quotas(self):
        genomes = [
            ("MoltenMetal", "Ni", "Bi", 10.0, 1000.0),
            ("MoltenMetal", "Cu", "Sn", 15.0, 1050.0),
            ("SolidCatalyst", "Ni", "Al2O3", "111", 0.0, (), 0, 0),
            ("SolidCatalyst", "Fe", "MgO", "100", 0.0, (), 0, 0),
        ]
        df = pd.DataFrame([
            {
                "genome": repr(genomes[0]),
                "material_class": "MoltenMetal",
                "valid": True,
                "E_act": 0.65,
                "coking_index": 0.9,
            },
            {
                "genome": repr(genomes[1]),
                "material_class": "MoltenMetal",
                "valid": True,
                "E_act": 1.10,
                "coking_index": -0.4,
            },
            {
                "genome": repr(genomes[2]),
                "material_class": "SolidCatalyst",
                "valid": True,
                "E_act": 0.70,
                "coking_index": 1.1,
            },
            {
                "genome": repr(genomes[3]),
                "material_class": "SolidCatalyst",
                "valid": True,
                "E_act": 1.25,
                "coking_index": -0.2,
            },
        ])
        objectives = np.array([
            [0.65, -0.9, 0.1, 0.2],
            [1.10, 0.4, 0.5, 0.6],
            [0.70, -1.1, 0.0, 0.1],
            [1.25, 0.2, 0.8, 0.9],
        ], dtype=float)
        ranks = pareto_ranks(objectives)
        self.assertEqual(int(ranks[0]), 0)
        self.assertEqual(int(ranks[2]), 0)
        self.assertGreater(int(ranks[1]), 0)
        self.assertGreater(int(ranks[3]), 0)

        selected_reactor = select_for_reactor(
            df,
            n_select=2,
            primary="E_act",
            min_per_class=1,
            objectives=["E_act", "-coking_index"],
        )
        self.assertEqual(
            set(selected_reactor["genome"]),
            {repr(genomes[0]), repr(genomes[2])},
        )

        with tempfile.TemporaryDirectory() as tmp:
            db_path = str(Path(tmp) / "validation.sqlite")
            batch_indices = allocate_validation_batch(
                genomes,
                objectives,
                n_select=2,
                database=db_path,
                application="turquoise_hydrogen",
                min_per_class=1,
                multi_objective=True,
            )
            self.assertEqual(set(batch_indices), {0, 2})

            slate_indices = experimental_slate(
                genomes,
                objectives,
                n_select=2,
                multi_objective=True,
            )
            self.assertEqual(set(slate_indices), {0, 2})

    def test_portfolio_and_coverage_certificate_include_objective_specs(self):
        turq_genomes, turq_df = self._make_turquoise_frame(24)
        ranker = fit_tree_ranker(turq_df, "turquoise_hydrogen")
        objs, _ = turquoise_tree_objectives(
            turq_genomes, ranker, return_uncertainty=True
        )
        criteria = (AdvancementCriterion("E_act", 1.5, "min"),)
        annotated = advancement_scorecard(turq_df, criteria)

        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "portfolio.json"
            write_portfolio_artifact(
                out_path,
                application="turquoise_hydrogen",
                protocol_id="pyrolysis_v1",
                criteria=criteria,
                candidates=turq_genomes,
                objectives=objs,
                evidence=annotated,
                indices=[0, 1, 2],
            )
            verification = verify_portfolio_artifact(out_path)
            self.assertTrue(verification["valid"])
            loaded = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertIn("objective_specs", loaded)
            self.assertEqual(len(loaded["objective_specs"]), 4)

            db_path = str(Path(tmp) / "branch.sqlite")
            cfg = BranchConfig(
                application="turquoise_hydrogen",
                database=db_path,
                material_classes=("MoltenMetal",),
                leaf_size=5000,
                max_leaves=1,
                probe_count=4,
            )
            run_branch_and_bound(
                cfg,
                lambda batch: turquoise_tree_objectives(batch, ranker),
            )
            cert = verify_branch_coverage(
                db_path,
                "turquoise_hydrogen",
                material_classes=("MoltenMetal",),
            )
            self.assertIn("objective_specs", cert)
            self.assertEqual(len(cert["objective_specs"]), 4)

    def test_validation_feedback_calibrates_tree_ranker_class_bias(self):
        from pipeline.search.adaptive_validation import record_validation

        turq_genomes, turq_df = self._make_turquoise_frame(24)
        baseline_ranker = fit_tree_ranker(turq_df, "turquoise_hydrogen")
        base_pred, _ = baseline_ranker.predict([turq_genomes[0]])

        with tempfile.TemporaryDirectory() as tmp:
            db_path = str(Path(tmp) / "feedback.sqlite")
            for idx in range(4):
                record_validation(
                    db_path,
                    "turquoise_hydrogen",
                    turq_genomes[idx],
                    predicted=0.80,
                    observed=1.20,
                    fidelity="dft",
                    productive=False,
                    provenance={"source_id": f"dft:{idx}"},
                )
            calibrated_ranker = fit_tree_ranker(
                turq_df,
                "turquoise_hydrogen",
                validation_database=db_path,
            )
            cal_pred, _ = calibrated_ranker.predict([turq_genomes[0]])
            self.assertGreater(float(cal_pred[0]), float(base_pred[0]))

    def test_calibration_artifact_seeds_db_and_calibrates_ranker(self):
        from pipeline.search.adaptive_validation import class_calibration
        from pipeline.search.discovery import candidate_id, discovery_region

        turq_genomes, turq_df = self._make_turquoise_frame(24)
        baseline_ranker = fit_tree_ranker(turq_df, "turquoise_hydrogen")
        base_pred, _ = baseline_ranker.predict([turq_genomes[0]])

        with tempfile.TemporaryDirectory() as tmp:
            db_path = str(Path(tmp) / "round2.sqlite")
            artifact_path = Path(tmp) / "active_learning_calibration.json"
            observations = [
                {
                    "application": "turquoise_hydrogen",
                    "candidate_id": candidate_id(turq_genomes[idx]),
                    "material_class": str(turq_genomes[idx][0]),
                    "region": "|".join(discovery_region(turq_genomes[idx])),
                    "genome": repr(turq_genomes[idx]),
                    "predicted_eV": 0.75,
                    "observed_eV": 1.25,
                    "fidelity": "screening_dft",
                    "productive": False,
                }
                for idx in range(4)
            ]
            artifact_path.write_text(
                json.dumps(
                    {
                        "n_validation_samples": len(observations),
                        "class_bias_eV": {},
                        "observations": observations,
                    }
                ),
                encoding="utf-8",
            )
            calibrated_ranker = fit_tree_ranker(
                turq_df,
                "turquoise_hydrogen",
                validation_database=db_path,
                calibration_artifact=str(artifact_path),
            )
            cal_pred, _ = calibrated_ranker.predict([turq_genomes[0]])
            self.assertGreater(float(cal_pred[0]), float(base_pred[0]))
            self.assertTrue(calibrated_ranker.region_bias)
            seeded_classes = class_calibration(
                db_path, "turquoise_hydrogen", fidelity="screening_dft"
            )
            self.assertIn(str(turq_genomes[0][0]), seeded_classes)

    def test_interfacial_carbon_wired_into_ranker_and_stage_selection(self):
        from pipeline.screening.small_data_ranker import _extract_secondary_targets
        from pipeline.screening.stage_selection import (
            annotate_interfacial_pyrolysis,
            select_for_pyrolysis_reactor,
        )

        # MoltenMetal row with NaN coking_index uses |c_height_above_surface_A|
        molten_row = {
            "material_class": "MoltenMetal",
            "coking_index": float("nan"),
            "coking_index_raw_slab": 4.10,
            "c_height_above_surface_A": -0.185,
            "dE_H": -0.01,
            "dE_C": 3.59,
        }
        sec = _extract_secondary_targets(molten_row, "turquoise_hydrogen")
        self.assertAlmostEqual(sec["coking_susceptibility_eV"], 0.185, places=6)

        # Stage selection annotates and Pareto-ranks by interfacial_pyrolysis_score_eV
        frame = pd.DataFrame([
            {
                "genome": "('MoltenMetal', 'Cu', 'Bi', 0.45, 1000)",
                "material_class": "MoltenMetal",
                "valid": True,
                "E_act": 0.80,
                "E_act_raw": 0.80,
                "E_act_censored": False,
                "pyrolysis_viable": True,
                "c_height_above_surface_A": -0.15,
            },
            {
                "genome": "('MoltenMetal', 'Fe', 'Ga', 0.20, 1000)",
                "material_class": "MoltenMetal",
                "valid": True,
                "E_act": 0.78,
                "E_act_raw": 0.78,
                "E_act_censored": False,
                "pyrolysis_viable": True,
                "c_height_above_surface_A": -2.40,
            },
        ])
        annotated = annotate_interfacial_pyrolysis(frame)
        self.assertIn("interfacial_carbon_displacement_A", annotated.columns)
        self.assertIn("interfacial_pyrolysis_score_eV", annotated.columns)
        self.assertAlmostEqual(
            float(annotated.loc[0, "interfacial_pyrolysis_score_eV"]), 0.95, places=6
        )
        selected = select_for_pyrolysis_reactor(frame, n_select=2, min_per_class=1)
        self.assertEqual(len(selected), 2)
        self.assertIn("interfacial_pyrolysis_score_eV", selected.columns)

    def test_multipathway_champion_scorecard_aggregates_all_six_pathways(self):
        from pipeline.economics.tea import build_multipathway_champion_scorecard

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "screening").mkdir(parents=True)
            (root / "reactor").mkdir(parents=True)
            (root / "fuel_cell").mkdir(parents=True)

            pd.DataFrame([
                {
                    "candidate_id": "cat_hea_1",
                    "genome": "('HEA', ('Fe', 'Zn', 'Nb', 'Ir', 'Sc'), 'hcp', '111', 800)",
                    "material_class": "HEA",
                    "valid": True,
                    "E_act": 0.3168,
                    "E_act_raw": 0.3168,
                    "E_act_censored": False,
                    "pyrolysis_viable": True,
                    "c_height_above_surface_A": -0.142,
                    "coking_index": 3.95,
                },
                {
                    "candidate_id": "cat_molten_1",
                    "genome": "('MoltenMetal', 'Te', 'Ti', 0.1, 750)",
                    "material_class": "MoltenMetal",
                    "valid": True,
                    "E_act": 0.8186,
                    "E_act_raw": 0.8186,
                    "E_act_censored": False,
                    "pyrolysis_viable": True,
                    "c_height_above_surface_A": -0.1855,
                    "coking_index": float("nan"),
                },
            ]).to_csv(root / "screening" / "branch_ranker_evidence.csv", index=False)

            reactor_fixtures = [
                ("PFR_cat_hea_1_1100K.json", {
                    "catalyst_name": "cat_hea_1",
                    "reactor_type": "PFR",
                    "status": "complete",
                    "T_K": 1100.0,
                    "CH4_conversion": 0.5348,
                    "single_pass_CH4_conversion": 0.5348,
                    "H2_selectivity": 0.998,
                    "surface_loaded": True,
                    "exceeds_equilibrium": False,
                    "carbon_balance_ok": True,
                    "catalyst_E_act_eV": 0.3168,
                    "catalyst_dE_H_eV": -0.55,
                    "kinetics_status": "mlip_neb_complete",
                }),
                ("Fluidized_cat_hea_1_1100K.json", {
                    "catalyst_name": "cat_hea_1",
                    "reactor_type": "Fluidized",
                    "status": "complete",
                    "T_K": 1100.0,
                    "CH4_conversion": 0.3412,
                    "single_pass_CH4_conversion": 0.3412,
                    "H2_selectivity": 0.995,
                    "surface_loaded": True,
                    "exceeds_equilibrium": False,
                    "carbon_balance_ok": True,
                    "catalyst_E_act_eV": 0.3168,
                    "catalyst_dE_H_eV": -0.55,
                    "kinetics_status": "mlip_neb_complete",
                }),
                ("MMBCR_cat_molten_1_1300K.json", {
                    "catalyst_name": "cat_molten_1",
                    "reactor_type": "MMBCR",
                    "status": "complete",
                    "T_K": 1300.0,
                    "CH4_conversion": 0.9461,
                    "single_pass_CH4_conversion": 0.9461,
                    "H2_selectivity": 1.0,
                    "exceeds_equilibrium": False,
                    "carbon_balance_ok": True,
                }),
                ("NTEC_cat_hea_1_473K.json", {
                    "catalyst_name": "cat_hea_1",
                    "reactor_type": "NTEC",
                    "status": "complete",
                    "T_K": 473.0,
                    "CH4_conversion": 0.5920,
                    "single_pass_CH4_conversion": 0.5920,
                    "H2_selectivity": 0.990,
                    "exceeds_equilibrium": False,
                    "carbon_balance_ok": True,
                }),
                ("Electrochemical_cat_hea_1_823K.json", {
                    "catalyst_name": "cat_hea_1",
                    "reactor_type": "Electrochemical",
                    "status": "complete",
                    "T_K": 823.0,
                    "CH4_conversion": 0.7204,
                    "single_pass_CH4_conversion": 0.7204,
                    "H2_selectivity": 0.985,
                    "exceeds_equilibrium": False,
                    "carbon_balance_ok": True,
                }),
            ]
            for fname, payload in reactor_fixtures:
                (root / "reactor" / fname).write_text(
                    json.dumps(payload), encoding="utf-8"
                )
            (root / "fuel_cell" / "pemfc_cat1_AquivionE87.json").write_text(
                json.dumps({
                    "cathode_catalyst": "008ae0ff99146112d2f0",
                    "membrane": "AquivionE87",
                    "peak_power_W_cm2": 1.385,
                    "peak_voltage_V": 0.66,
                    "peak_current_A_cm2": 2.10,
                    "orr_overpotential_V": 0.285,
                    "efficiency_at_peak": 0.536,
                }),
                encoding="utf-8",
            )

            scorecard = build_multipathway_champion_scorecard(str(root))
            self.assertEqual(
                scorecard["evidence_level"], "screening_multipathway_scorecard"
            )
            self.assertEqual(len(scorecard["discovery_champions"]), 2)
            self.assertAlmostEqual(
                float(scorecard["discovery_champions"][0]["interfacial_pyrolysis_score_eV"]),
                0.4588,
                places=4,
            )
            leaderboard = scorecard["pathway_leaderboard"]
            for mode_name in (
                "thermocatalytic_pfr",
                "thermocatalytic_fluidized",
                "thermocatalytic",
                "mmbcr",
                "ntec",
                "electrochemical",
            ):
                self.assertIn(mode_name, leaderboard)
                self.assertEqual(leaderboard[mode_name]["status"], "complete")
                self.assertGreater(
                    float(leaderboard[mode_name]["single_pass_CH4_conversion"]), 0.0
                )
                self.assertIsNotNone(leaderboard[mode_name]["tea_base_h2_usd_kg"])
            pemfc = scorecard["pemfc_stack_summary"]
            self.assertEqual(pemfc["status"], "complete")
            self.assertGreater(float(pemfc["stack_net_power_kW"]), 0.0)
            self.assertEqual(len(scorecard["coupled_turquoise_pemfc_tea"]), 6)


if __name__ == "__main__":
    unittest.main()

