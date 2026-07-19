#!/usr/bin/env python3
"""Lock or analyze a prospective all-class computational pilot."""

from __future__ import annotations

import argparse
import ast
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.indexed_space import deterministic_tree_probes
from pipeline.pilot_benchmark import (_identity, _ridge_predict, default_specs,
                                      expert_score, load_legacy_outcomes)


ROOT = Path("results/pilot")
MANIFEST = ROOT / "prospective_manifest.json"


def prepare(pool_size: int = 42) -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    candidates = deterministic_tree_probes(pool_size + 14)[14:14 + pool_size]
    records = []
    for spec in default_specs():
        legacy = load_legacy_outcomes(spec)
        pool = pd.DataFrame({"parsed_genome": candidates,
                             "candidate_id": [_identity(g) for g in candidates]})
        predicted = _ridge_predict(legacy, pool, spec.outcome)
        expert = np.array([expert_score(g, spec.application) for g in candidates])
        budget = max(1, math.ceil(pool_size * spec.selection_fraction))
        learned_order = np.argsort(predicted, kind="stable")
        expert_order = np.lexsort((pool.candidate_id.to_numpy(), expert))
        records.append({
            "application": spec.application, "outcome": spec.outcome,
            "budget": budget, "pool": [repr(g) for g in candidates],
            "candidate_ids": pool.candidate_id.tolist(),
            "learned_predictions": predicted.tolist(), "expert_scores": expert.tolist(),
            "learned_selected_ids": pool.candidate_id.iloc[learned_order[:budget]].tolist(),
            "expert_selected_ids": pool.candidate_id.iloc[expert_order[:budget]].tolist(),
            "random_seed": 20260719, "random_trials": 20_000,
        })
    payload = {"schema_version": 1, "locked_before_outcomes": True,
               "created_utc": datetime.now(timezone.utc).isoformat(), "pilots": records}
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Locked {pool_size} candidates and both selection slates in {MANIFEST}")


def analyze() -> None:
    payload = json.loads(MANIFEST.read_text())
    replay_path = ROOT / "production_selector_replay.json"
    replay = json.loads(replay_path.read_text()) if replay_path.exists() else None
    paths = {"turquoise_hydrogen": Path("results/screening/pilot/prospective_pyrolysis.csv"),
             "fuel_cell_orr": Path("results/fuel_cell/pilot/prospective_orr.csv")}
    results = []
    for pilot in payload["pilots"]:
        frame = pd.read_csv(paths[pilot["application"]])
        frame = frame[frame.valid.eq(True)].copy()
        frame[pilot["outcome"]] = pd.to_numeric(frame[pilot["outcome"]], errors="coerce")
        frame = frame[np.isfinite(frame[pilot["outcome"]])]
        frame["candidate_id"] = [_identity(ast.literal_eval(g)) for g in frame.genome]
        values = dict(zip(frame.candidate_id, frame[pilot["outcome"]]))
        candidate_ids = pilot["candidate_ids"]
        available = [cid for cid in candidate_ids if cid in values]
        ordered = np.array([values[cid] for cid in available], float)
        cutoff = float(np.quantile(ordered, .20))
        predicted_by_id = dict(zip(candidate_ids, pilot["learned_predictions"]))
        predicted_available = np.array([predicted_by_id[cid] for cid in available], float)
        errors = predicted_available - ordered
        rank_correlation = float(pd.Series(predicted_available).corr(
            pd.Series(ordered), method="spearman"))
        def stats(ids):
            vals = np.array([values.get(cid, np.inf) for cid in ids], float)
            finite = vals[np.isfinite(vals)]
            return {"selected": len(vals), "valid_selected": len(finite),
                    "hits": int(np.sum(vals <= cutoff)),
                    "hit_rate": float(np.mean(vals <= cutoff)),
                    "best": float(finite.min()) if len(finite) else None,
                    "mean_valid": float(finite.mean()) if len(finite) else None}
        learned = stats(pilot["learned_selected_ids"])
        expert = stats(pilot["expert_selected_ids"])
        rng = np.random.default_rng(pilot["random_seed"])
        random_rates = []
        for _ in range(pilot["random_trials"]):
            ids = rng.choice(candidate_ids, size=pilot["budget"], replace=False)
            random_rates.append(stats(ids)["hit_rate"])
        random_rates = np.array(random_rates)
        unique, counts = np.unique(ordered, return_counts=True)
        def details(ids):
            return [{"candidate_id": cid, "genome": pilot["pool"][candidate_ids.index(cid)],
                     "predicted": float(predicted_by_id[cid]),
                     "outcome": float(values[cid]) if cid in values else None,
                     "hit": bool(cid in values and values[cid] <= cutoff)} for cid in ids]
        actual_top = sorted(available, key=lambda cid: values[cid])[:pilot["budget"]]
        result = {"application": pilot["application"], "evaluated": len(available),
                        "validity_rate": len(available) / len(pilot["candidate_ids"]),
                        "outcome": pilot["outcome"], "hit_cutoff": cutoff,
                        "unique_outcomes": len(unique), "largest_tie": int(counts.max()),
                        "prediction_mae": float(np.mean(np.abs(errors))),
                        "prediction_rmse": float(np.sqrt(np.mean(errors ** 2))),
                        "spearman_rank_correlation": rank_correlation,
                        "learned": learned, "expert": expert,
                        "random_mean_hit_rate": float(random_rates.mean()),
                        "random_95pct": [float(np.quantile(random_rates, .025)),
                                         float(np.quantile(random_rates, .975))],
                        "learned_enrichment_vs_random": learned["hit_rate"] / random_rates.mean(),
                        "beats_random_95pct": bool(learned["hit_rate"] > np.quantile(random_rates, .975)),
                        "beats_expert": bool(learned["hit_rate"] > expert["hit_rate"]),
                        "learned_slate": details(pilot["learned_selected_ids"]),
                        "expert_slate": details(pilot["expert_selected_ids"]),
                        "actual_top": details(actual_top)}
        if replay:
            replay_result = next(x for x in replay["results"]
                                 if x["application"] == pilot["application"])
            eligible_ids = replay_result["eligible_candidate_ids"]
            eligible_valid = [cid for cid in eligible_ids if cid in values]
            eligible_cutoff = float(np.quantile([values[cid] for cid in eligible_valid], .20))
            production_scores = np.array([replay_result["primary_scores"][cid]
                                          for cid in eligible_valid], float)
            production_outcomes = np.array([values[cid] for cid in eligible_valid], float)
            production_rank_correlation = float(pd.Series(production_scores).corr(
                pd.Series(production_outcomes), method="spearman"))
            cutoff_original = cutoff
            cutoff = eligible_cutoff
            production = stats(replay_result["production_selected_ids"])
            expert_order = sorted(eligible_ids, key=lambda cid: (
                pilot["expert_scores"][candidate_ids.index(cid)], cid))[:pilot["budget"]]
            eligible_expert = stats(expert_order)
            rng_replay = np.random.default_rng(pilot["random_seed"])
            replay_random = np.array([
                stats(rng_replay.choice(eligible_ids, pilot["budget"], replace=False))["hit_rate"]
                for _ in range(pilot["random_trials"])
            ])
            cutoff = cutoff_original
            result["production_replay"] = {
                "trained_only_on_pre_pilot_evidence": True,
                "eligible_pool_size": len(eligible_ids), "eligible_hit_cutoff": eligible_cutoff,
                "spearman_rank_correlation": production_rank_correlation,
                "production": production, "eligible_expert": eligible_expert,
                "random_mean_hit_rate": float(replay_random.mean()),
                "random_95pct": [float(np.quantile(replay_random, .025)),
                                  float(np.quantile(replay_random, .975))],
                "enrichment_vs_random": production["hit_rate"] / replay_random.mean(),
                "beats_random_95pct": bool(production["hit_rate"] > np.quantile(replay_random, .975)),
                "beats_expert": bool(production["hit_rate"] > eligible_expert["hit_rate"]),
                "production_slate": details(replay_result["production_selected_ids"]),
            }
        results.append(result)
    out = ROOT / "prospective_analysis.json"
    out.write_text(json.dumps({"manifest": str(MANIFEST), "results": results}, indent=2) + "\n")
    print(json.dumps(results, indent=2)); print(f"Analysis: {out}")


def evaluate(application: str) -> None:
    payload = json.loads(MANIFEST.read_text())
    pilot = next(x for x in payload["pilots"] if x["application"] == application)
    pool = [ast.literal_eval(raw) for raw in pilot["pool"]]
    if application == "turquoise_hydrogen":
        from pipeline.surface_screener import run_screening
        frame = run_screening(pool, db_filename="pilot/prospective_pyrolysis.csv",
                              workers_per_gpu=1)
    else:
        from pipeline.fc_screener import run_orr_screening
        frame = run_orr_screening(pool, db_filename="pilot/prospective_orr.csv",
                                  workers_per_gpu=1)
    print(f"{application}: {len(frame)} evaluated, {int(frame.valid.eq(True).sum())} valid")


def replay_production() -> None:
    """Replay production selectors trained only on pre-pilot evidence."""
    import torch
    from pipeline.catalyst_spaces import encode_population
    from pipeline.genetic_optimizer import (_train_ensemble_from_db,
                                            compute_objectives_surrogate)
    from pipeline.fc_genetic_optimizer import (_train_orr_ensemble_from_db,
                                               compute_orr_objectives_surrogate)
    from pipeline.application_scope import pemfc_cathode_scope
    torch.manual_seed(20260719); np.random.seed(20260719)
    payload = json.loads(MANIFEST.read_text())
    outputs = []
    for spec in default_specs():
        pilot = next(x for x in payload["pilots"] if x["application"] == spec.application)
        pool = [ast.literal_eval(raw) for raw in pilot["pool"]]
        legacy = load_legacy_outcomes(spec)
        if spec.application == "turquoise_hydrogen":
            model = _train_ensemble_from_db(legacy, "cuda", n_models=3)
            objectives = compute_objectives_surrogate(pool, model, "cuda")
            eligible = list(range(len(pool)))
        else:
            model = _train_orr_ensemble_from_db(legacy, "cuda", n_models=3)
            objectives = compute_orr_objectives_surrogate(pool, model, "cuda")
            eligible = [i for i, g in enumerate(pool)
                        if pemfc_cathode_scope(g)["status"] == "candidate"]
        order = sorted(eligible, key=lambda i: (float(objectives[i, 0]), pilot["candidate_ids"][i]))
        selected = order[:min(pilot["budget"], len(order))]
        outputs.append({"application": spec.application, "seed": 20260719,
                        "training_candidate_ids": legacy.candidate_id.tolist(),
                        "eligible_pool_size": len(eligible),
                        "eligible_candidate_ids": [pilot["candidate_ids"][i] for i in eligible],
                        "production_selected_ids": [pilot["candidate_ids"][i] for i in selected],
                        "primary_scores": {pilot["candidate_ids"][i]: float(objectives[i, 0])
                                           for i in range(len(pool))}})
        del model
        torch.cuda.empty_cache()
    path = ROOT / "production_selector_replay.json"
    path.write_text(json.dumps({"trained_only_on_pre_pilot_evidence": True,
                                "results": outputs}, indent=2) + "\n")
    print(f"Production replay: {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare", "analyze", "evaluate-pyrolysis",
                                           "evaluate-orr", "replay-production"))
    parser.add_argument("--pool-size", type=int, default=42)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.pool_size)
    elif args.action == "evaluate-pyrolysis":
        evaluate("turquoise_hydrogen")
    elif args.action == "evaluate-orr":
        evaluate("fuel_cell_orr")
    elif args.action == "replay-production":
        replay_production()
    else:
        analyze()
