#!/usr/bin/env python3
"""Second-round locked pilot for the improved divide-and-conquer selector."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.common.application_scope import pemfc_cathode_scope
from pipeline.search.branch_search import _probe_indices
from pipeline.common.catalyst_spaces import encode_population
from pipeline.search.discovery import candidate_id
from pipeline.search.indexed_space import CLASS_ORDER, CLASS_SIZES, candidate_at_class
from pipeline.evidence.pilot_benchmark import default_specs, load_legacy_outcomes


ROUND = os.environ.get("PILOT_ROUND", "v2")
ROOT = Path("results/pilot") / f"divide_conquer_{ROUND}"
MANIFEST = ROOT / "manifest.json"


def _pool(per_class: int = 4) -> list[tuple]:
    """Fresh low-discrepancy candidates, excluding the first pilot pool."""
    prior = set()
    old = Path("results/pilot/prospective_manifest.json")
    if old.exists():
        payload = json.loads(old.read_text())
        prior.update(candidate_id(ast.literal_eval(raw)) for raw in payload["pilots"][0]["pool"])
    for manifest_path in sorted(Path("results/pilot").glob("divide_conquer_v*/manifest.json")):
        if manifest_path.parent == ROOT:
            continue
        payload = json.loads(manifest_path.read_text())
        prior.update(candidate_id(ast.literal_eval(raw))
                     for raw in payload["records"][0]["pool"])
    pool = []
    for material_class in CLASS_ORDER:
        # Draw extra deterministic points so exclusion never shrinks a class.
        indices = _probe_indices(0, CLASS_SIZES[material_class], per_class * 5 + 7)
        choices = []
        for index in indices:
            genome = candidate_at_class(material_class, index)
            if candidate_id(genome) not in prior:
                choices.append(genome)
            if len(choices) == per_class:
                break
        if len(choices) != per_class:
            raise RuntimeError(f"could not create fresh pool for {material_class}")
        pool.extend(choices)
    return pool


def _training_frame(application: str) -> pd.DataFrame:
    spec = next(x for x in default_specs() if x.application == application)
    legacy = load_legacy_outcomes(spec).drop(columns=["parsed_genome", "candidate_id", "replicates"])
    extra_path = (Path("results/screening/pilot/prospective_pyrolysis.csv")
                  if application == "turquoise_hydrogen" else
                  Path("results/fuel_cell/pilot/prospective_orr.csv"))
    extra = pd.read_csv(extra_path)
    frames = [legacy, extra]
    if ROUND != "v2":
        subdir = "screening" if application == "turquoise_hydrogen" else "fuel_cell"
        suffix = "pyrolysis" if application == "turquoise_hydrogen" else "orr"
        for prior_round in sorted(Path(f"results/{subdir}/pilot").glob(
                f"divide_conquer_v*_{suffix}.csv")):
            if f"divide_conquer_{ROUND}_" not in prior_round.name:
                frames.append(pd.read_csv(prior_round))
    combined = pd.concat(frames, ignore_index=True, sort=False)
    combined = combined[combined.valid.eq(True)].copy()
    combined["_cid"] = [candidate_id(ast.literal_eval(g) if isinstance(g, str) else tuple(g))
                        for g in combined.genome]
    return combined.drop_duplicates("_cid", keep="last").drop(columns="_cid")


def prepare(per_class: int = 4, extra_slots: int = 6) -> None:
    import torch
    from pipeline.search.adaptive_validation import allocate_validation_batch
    from pipeline.screening.genetic_optimizer import (_train_ensemble_from_db,
                                            compute_objectives_surrogate)
    from pipeline.screening.fc_genetic_optimizer import (_train_orr_ensemble_from_db,
                                               compute_orr_objectives_surrogate,
                                               ORRSurrogateEnsemble)
    from pipeline.screening.surrogate_model import SurrogateEnsemble, predict_ensemble

    ROOT.mkdir(parents=True, exist_ok=True)
    pool = _pool(per_class)
    records = []
    for application in ("turquoise_hydrogen", "fuel_cell_orr"):
        torch.manual_seed(20260720); np.random.seed(20260720)
        training = _training_frame(application)
        if ROUND != "v2":
            from sklearn.ensemble import ExtraTreesRegressor
            outcome = "E_act" if application == "turquoise_hydrogen" else "orr_overpotential_V"
            eligible = (list(range(len(pool))) if application == "turquoise_hydrogen" else
                        [i for i, genome in enumerate(pool)
                         if pemfc_cathode_scope(genome)["status"] == "candidate"])
            training = training[np.isfinite(pd.to_numeric(training[outcome], errors="coerce"))]
            x_train = encode_population([
                ast.literal_eval(g) if isinstance(g, str) else tuple(g) for g in training.genome])
            y_train = training[outcome].to_numpy(float)
            leaf = 1 if application == "turquoise_hydrogen" else 3
            x_pool = encode_population(pool)
            if application == "fuel_cell_orr" and ROUND not in ("v2", "v3"):
                adsorption = ["dG_OH_eV", "dG_O_eV", "dG_OOH_eV"]
                finite = training[adsorption].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
                model = ExtraTreesRegressor(n_estimators=1024, min_samples_leaf=1,
                                            max_features=1.0, random_state=20260721, n_jobs=-1)
                model.fit(x_train[finite.to_numpy()], training.loc[finite, adsorption].to_numpy(float))
                def eta(prediction):
                    d_oh, d_o, d_ooh = prediction.T
                    return 1.23 + np.maximum.reduce([d_ooh - 4.92, d_o - d_ooh,
                                                     d_oh - d_o, -d_oh])
                primary = eta(model.predict(x_pool))
                members = np.column_stack([eta(tree.predict(x_pool)) for tree in model.estimators_])
            else:
                model = ExtraTreesRegressor(n_estimators=1024, min_samples_leaf=leaf,
                                            max_features=1.0, random_state=20260721, n_jobs=-1)
                model.fit(x_train, y_train)
                primary = model.predict(x_pool)
                members = np.column_stack([tree.predict(x_pool) for tree in model.estimators_])
            uncertainty = members.std(axis=1)
            objectives = np.column_stack([primary, np.zeros((len(pool), 3))])
        elif application == "turquoise_hydrogen":
            eligible = list(range(len(pool)))
            model = _train_ensemble_from_db(training, "cuda", n_models=3)
            objectives = compute_objectives_surrogate(pool, model, "cuda")
            uncertainty = (predict_ensemble(model, encode_population(pool), "cuda")["E_act_std"]
                           if isinstance(model, SurrogateEnsemble) else np.zeros(len(pool)))
        else:
            eligible = [i for i, genome in enumerate(pool)
                        if pemfc_cathode_scope(genome)["status"] == "candidate"]
            model = _train_orr_ensemble_from_db(training, "cuda", n_models=3)
            objectives = compute_orr_objectives_surrogate(pool, model, "cuda")
            uncertainty = np.zeros(len(pool))
            if isinstance(model, ORRSurrogateEnsemble):
                x = torch.FloatTensor(encode_population(pool)).to("cuda")
                members = []
                for member in model.models:
                    member.eval()
                    with torch.no_grad():
                        _, eta, _ = member(x)
                    members.append(eta.cpu().numpy().ravel())
                uncertainty = np.column_stack(members).std(axis=1)

        candidates = [pool[i] for i in eligible]
        classes = {g[0] for g in candidates}
        budget = min(len(candidates), len(classes) + extra_slots)
        if ROUND == "v2":
            eligible_objectives = objectives[eligible]
            eligible_uncertainty = uncertainty[eligible]
            local_selected = allocate_validation_batch(
                candidates, eligible_objectives, budget,
                str(ROOT / "allocation.sqlite"), application,
                min_per_class=1, uncertainties=eligible_uncertainty)
            selected = [eligible[i] for i in local_selected]
            coverage_selected = selected
        else:
            # Discovery and calibration are distinct ledgers. Discovery is the
            # reproducible priority ranking tested against random; calibration
            # retains one uncertain representative per chemistry class.
            selected = sorted(eligible, key=lambda i: (float(objectives[i, 0]), candidate_id(pool[i])))[:budget]
            coverage_selected = []
            for material_class in sorted(classes):
                indices = [i for i in eligible if pool[i][0] == material_class]
                coverage_selected.append(max(indices, key=lambda i: (uncertainty[i], candidate_id(pool[i]))))
        ids = [candidate_id(g) for g in pool]
        records.append({
            "application": application, "training_rows": len(training),
            "training_digest": hashlib.sha256("\n".join(sorted(
                candidate_id(ast.literal_eval(g) if isinstance(g, str) else tuple(g))
                for g in training.genome)).encode()).hexdigest(),
            "pool": [repr(g) for g in pool], "candidate_ids": ids,
            "eligible_ids": [ids[i] for i in eligible], "budget": budget,
            "selected_ids": [ids[i] for i in selected],
            "coverage_validation_ids": [ids[i] for i in coverage_selected],
            "primary_scores": {ids[i]: float(objectives[i, 0]) for i in eligible},
            "uncertainties": {ids[i]: float(uncertainty[i]) for i in eligible},
            "random_seed": 20260720, "random_trials": 50_000,
        })
        del model
        torch.cuda.empty_cache()
    selector = ("class-floor + improvement + uncertainty + calibration" if ROUND == "v2" else
                "deterministic small-data tree ranking; separate class calibration slate")
    payload = {"schema_version": 1, "locked_before_outcomes": True,
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "selector": selector,
               "records": records}
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Locked divide-and-conquer pilot: {MANIFEST}")


def evaluate(application: str) -> None:
    payload = json.loads(MANIFEST.read_text())
    record = next(x for x in payload["records"] if x["application"] == application)
    allowed = set(record["eligible_ids"])
    pool = [ast.literal_eval(raw) for raw, cid in zip(record["pool"], record["candidate_ids"])
            if cid in allowed]
    if application == "turquoise_hydrogen":
        from pipeline.screening.surface_screener import run_screening
        frame = run_screening(pool, db_filename=f"pilot/divide_conquer_{ROUND}_pyrolysis.csv",
                              workers_per_gpu=1)
    else:
        from pipeline.screening.fc_screener import run_orr_screening
        frame = run_orr_screening(pool, db_filename=f"pilot/divide_conquer_{ROUND}_orr.csv",
                                  workers_per_gpu=1)
    print(f"{application}: evaluated={len(frame)}, valid={int(frame.valid.eq(True).sum())}")


def analyze() -> None:
    payload = json.loads(MANIFEST.read_text())
    paths = {"turquoise_hydrogen": Path(f"results/screening/pilot/divide_conquer_{ROUND}_pyrolysis.csv"),
             "fuel_cell_orr": Path(f"results/fuel_cell/pilot/divide_conquer_{ROUND}_orr.csv")}
    results = []
    combined_selected = combined_hits = 0
    random_trial_hits = None
    for record in payload["records"]:
        frame = pd.read_csv(paths[record["application"]])
        outcome = "E_act" if record["application"] == "turquoise_hydrogen" else "orr_overpotential_V"
        values = {}
        for _, row in frame.iterrows():
            if bool(row.valid) and np.isfinite(float(row[outcome])):
                genome = ast.literal_eval(row.genome)
                values[candidate_id(genome)] = float(row[outcome])
        eligible = record["eligible_ids"]
        valid_values = np.array([values[cid] for cid in eligible if cid in values])
        cutoff = float(np.quantile(valid_values, .20))
        def hits(ids):
            return sum(cid in values and values[cid] <= cutoff for cid in ids)
        selected_hits = hits(record["selected_ids"])
        rng = np.random.default_rng(record["random_seed"])
        random_hits = np.array([hits(rng.choice(eligible, record["budget"], replace=False))
                                for _ in range(record["random_trials"])])
        # Policy-matched random control: the same one-per-class coverage floor,
        # with remaining slots sampled uniformly. This isolates whether scoring
        # adds value beyond the coverage policy itself.
        by_class = {}
        genome_by_id = dict(zip(record["candidate_ids"], map(ast.literal_eval, record["pool"])))
        for cid in eligible:
            by_class.setdefault(genome_by_id[cid][0], []).append(cid)
        matched_random_hits = []
        for _ in range(record["random_trials"]):
            chosen = [rng.choice(by_class[material_class]) for material_class in sorted(by_class)]
            remaining = [cid for cid in eligible if cid not in chosen]
            chosen.extend(rng.choice(remaining, record["budget"] - len(chosen), replace=False))
            matched_random_hits.append(hits(chosen))
        matched_random_hits = np.asarray(matched_random_hits)
        scored_valid = [cid for cid in eligible if cid in values]
        score_values = np.array([record["primary_scores"][cid] for cid in scored_valid])
        outcome_values = np.array([values[cid] for cid in scored_valid])
        rank_correlation = float(pd.Series(score_values).corr(
            pd.Series(outcome_values), method="spearman"))
        primary_only = sorted(eligible, key=lambda cid: (record["primary_scores"][cid], cid))[
            :record["budget"]]
        primary_hits = hits(primary_only)
        unique, counts = np.unique(valid_values, return_counts=True)
        hit_rate = selected_hits / record["budget"]
        random_rate = float(random_hits.mean() / record["budget"])
        results.append({"application": record["application"], "pool": len(eligible),
                        "valid": len(values), "budget": record["budget"],
                        "hit_cutoff": cutoff, "hits": selected_hits, "hit_rate": hit_rate,
                        "random_mean_hits": float(random_hits.mean()),
                        "random_mean_hit_rate": random_rate,
                        "random_hits_95pct": [float(np.quantile(random_hits, .025)),
                                              float(np.quantile(random_hits, .975))],
                        "coverage_matched_random_mean_hits": float(matched_random_hits.mean()),
                        "coverage_matched_random_95pct": [
                            float(np.quantile(matched_random_hits, .025)),
                            float(np.quantile(matched_random_hits, .975))],
                        "enrichment_vs_coverage_matched_random":
                            selected_hits / matched_random_hits.mean(),
                        "beats_coverage_matched_random_95pct":
                            bool(selected_hits > np.quantile(matched_random_hits, .975)),
                        "enrichment_vs_random": hit_rate / random_rate,
                        "beats_random_95pct": bool(selected_hits > np.quantile(random_hits, .975)),
                        "spearman_rank_correlation": rank_correlation,
                        "unique_outcomes": len(unique), "largest_tie": int(counts.max()),
                        "primary_only_diagnostic_hits": primary_hits,
                        "primary_only_diagnostic_enrichment":
                            (primary_hits / record["budget"]) / random_rate})
        combined_selected += record["budget"]; combined_hits += selected_hits
        random_trial_hits = random_hits if random_trial_hits is None else random_trial_hits + random_hits
    combined = {"selected": combined_selected, "hits": combined_hits,
                "hit_rate": combined_hits / combined_selected,
                "random_mean_hits": float(random_trial_hits.mean()),
                "random_hits_95pct": [float(np.quantile(random_trial_hits, .025)),
                                      float(np.quantile(random_trial_hits, .975))],
                "enrichment_vs_random": combined_hits / random_trial_hits.mean(),
                "beats_random_95pct": bool(combined_hits > np.quantile(random_trial_hits, .975))}
    output = {"results": results, "combined": combined}
    path = ROOT / "analysis.json"; path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2)); print(f"Analysis: {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare", "evaluate-pyrolysis", "evaluate-orr", "analyze"))
    parser.add_argument("--per-class", type=int, default=4)
    parser.add_argument("--extra-slots", type=int, default=6)
    args = parser.parse_args()
    if args.action == "prepare": prepare(args.per_class, args.extra_slots)
    elif args.action == "evaluate-pyrolysis": evaluate("turquoise_hydrogen")
    elif args.action == "evaluate-orr": evaluate("fuel_cell_orr")
    else: analyze()
