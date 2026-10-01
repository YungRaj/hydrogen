"""Chronological benchmark for catalyst discovery acquisition policies.

Each pilot round is predicted from outcomes available before that round.  This
prevents a completed candidate from influencing its own selection score and
makes the comparison a search-policy test rather than an in-sample fit metric.
"""

from __future__ import annotations

import ast
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from pipeline.evidence.pilot_benchmark import default_specs, load_legacy_outcomes
from pipeline.search.design_space import encode_population
from pipeline.search.discovery import candidate_id
from pipeline.screening.small_data_ranker import (
    CLASS_BIAS_PRIOR_ROWS,
    CLASS_BIAS_WEIGHT,
)


@dataclass(frozen=True)
class PolicyResult:
    """Aggregate equal-budget discovery result for one acquisition policy."""

    hits: int
    selected: int
    hit_rate: float
    development_hits: int
    frozen_hits: int


@dataclass(frozen=True)
class BenchmarkSpec:
    """Inputs and deterministic settings for the chronological comparison."""

    root: str = "results"
    rounds: tuple[str, ...] = ("v2", "v3", "v4", "v5", "v6", "v7")
    random_trials: int = 50_000
    random_seed: int = 20260930


def _genome(raw: object) -> tuple:
    value = ast.literal_eval(raw) if isinstance(raw, str) else tuple(raw)
    if not isinstance(value, tuple) or not value:
        raise ValueError("invalid catalyst genome")
    return value


def _clean(frame: pd.DataFrame, outcome: str) -> pd.DataFrame:
    cleaned = frame.copy()
    cleaned["_genome"] = [_genome(raw) for raw in cleaned["genome"]]
    cleaned["_candidate_id"] = [candidate_id(g) for g in cleaned["_genome"]]
    cleaned[outcome] = pd.to_numeric(cleaned[outcome], errors="coerce")
    return cleaned.drop_duplicates("_candidate_id", keep="last")


def _initial_training(root: Path, application: str, outcome: str) -> pd.DataFrame:
    spec = next(spec for spec in default_specs() if spec.application == application)
    legacy = load_legacy_outcomes(spec).drop(
        columns=["parsed_genome", "candidate_id", "replicates"]
    )
    path = (
        root / "screening/pilot/prospective_pyrolysis.csv"
        if application == "turquoise_hydrogen"
        else root / "fuel_cell/pilot/prospective_orr.csv"
    )
    return _clean(pd.concat([legacy, pd.read_csv(path)], ignore_index=True), outcome)


def _round_path(root: Path, round_name: str, application: str) -> Path:
    if application == "turquoise_hydrogen":
        return root / f"screening/pilot/divide_conquer_{round_name}_pyrolysis.csv"
    return root / f"fuel_cell/pilot/divide_conquer_{round_name}_orr.csv"


def _coverage_select(
    scores: np.ndarray, genomes: Sequence[tuple], budget: int
) -> np.ndarray:
    """Apply the production class floor, then spend remaining slots by score."""
    by_class: dict[str, list[int]] = {}
    for index, genome in enumerate(genomes):
        by_class.setdefault(str(genome[0]), []).append(index)
    selected = [
        max(indices, key=lambda i: (float(scores[i]), candidate_id(genomes[i])))
        for _, indices in sorted(by_class.items())
    ]
    selected.sort(key=lambda i: (-float(scores[i]), candidate_id(genomes[i])))
    chosen = set(selected)
    remainder = sorted(
        (i for i in range(len(genomes)) if i not in chosen),
        key=lambda i: (-float(scores[i]), candidate_id(genomes[i])),
    )
    return np.asarray((selected + remainder)[:budget], dtype=int)


def run_search_policy_benchmark(spec: BenchmarkSpec = BenchmarkSpec()) -> dict:
    """Compare catalyst quality with random, uncertainty, and validity ranking.

    Args:
        spec: Input locations, locked rounds, and deterministic random settings.

    Returns:
        Chronological per-application and combined policy metrics.
    """
    from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor
    from sklearn.model_selection import KFold

    root = Path(spec.root)
    rng = np.random.default_rng(spec.random_seed)
    applications: list[dict] = []
    combined_hits = {name: 0 for name in ("catalyst", "uncertainty", "validity")}
    combined_selected = 0
    combined_random: np.ndarray | None = None

    for application, outcome in (
        ("turquoise_hydrogen", "E_act"),
        ("fuel_cell_orr", "orr_overpotential_V"),
    ):
        training = _initial_training(root, application, outcome)
        policy_hits = {name: 0 for name in combined_hits}
        development_hits = {name: 0 for name in combined_hits}
        frozen_hits = {name: 0 for name in combined_hits}
        selected_count = 0
        random_hits = np.zeros(spec.random_trials, dtype=int)
        rounds: list[dict] = []
        for round_name in spec.rounds:
            manifest = json.loads(
                (root / f"pilot/divide_conquer_{round_name}/manifest.json").read_text()
            )
            record = next(
                row for row in manifest["records"] if row["application"] == application
            )
            observed = _clean(
                pd.read_csv(_round_path(root, round_name, application)), outcome
            )
            rows = {
                str(row["_candidate_id"]): row for _, row in observed.iterrows()
            }
            ids = [value for value in record["eligible_ids"] if value in rows]
            genomes = [rows[value]["_genome"] for value in ids]
            outcomes = np.asarray([rows[value][outcome] for value in ids], dtype=float)
            valid = np.asarray(
                [bool(rows[value]["valid"]) for value in ids], dtype=bool
            ) & np.isfinite(outcomes)

            train_valid = training["valid"].eq(True).to_numpy() & np.isfinite(
                training[outcome].to_numpy(float)
            )
            train_x = encode_population(training["_genome"].tolist())
            test_x = encode_population(genomes)
            ranker = ExtraTreesRegressor(
                n_estimators=256,
                min_samples_leaf=2,
                max_features=1.0,
                random_state=20260721,
                n_jobs=-1,
            ).fit(train_x[train_valid], training.loc[train_valid, outcome])
            members = np.vstack([tree.predict(test_x) for tree in ranker.estimators_])
            mean = members.mean(axis=0)
            uncertainty = members.std(axis=0)
            train_y = training.loc[train_valid, outcome].to_numpy(float)
            valid_x = train_x[train_valid]
            out_of_fold = np.empty(len(train_y), dtype=float)
            folds = KFold(n_splits=5, shuffle=True, random_state=77)
            for train_indices, test_indices in folds.split(valid_x):
                calibration = ExtraTreesRegressor(
                    n_estimators=64,
                    min_samples_leaf=2,
                    max_features=1.0,
                    random_state=21,
                    n_jobs=-1,
                ).fit(valid_x[train_indices], train_y[train_indices])
                out_of_fold[test_indices] = calibration.predict(valid_x[test_indices])
            residual = train_y - out_of_fold
            global_bias = float(residual.mean())
            training_classes = np.asarray(
                [str(genome[0]) for genome in training.loc[train_valid, "_genome"]]
            )
            class_bias: dict[str, float] = {}
            for material_class in np.unique(training_classes):
                mask = training_classes == material_class
                class_bias[material_class] = float(
                    (
                        residual[mask].sum()
                        + CLASS_BIAS_PRIOR_ROWS * global_bias
                    )
                    / (mask.sum() + CLASS_BIAS_PRIOR_ROWS)
                )
            mean += CLASS_BIAS_WEIGHT * np.asarray(
                [class_bias.get(str(genome[0]), global_bias) for genome in genomes]
            )

            validity = ExtraTreesClassifier(
                n_estimators=256,
                min_samples_leaf=2,
                max_features=1.0,
                class_weight="balanced",
                random_state=20260722,
                n_jobs=-1,
            ).fit(train_x, training["valid"].eq(True).to_numpy(dtype=int))
            valid_class = list(validity.classes_).index(1)
            validity_probability = validity.predict_proba(test_x)[:, valid_class]

            cutoff = float(np.quantile(outcomes[valid], 0.20))
            hits = valid & (outcomes <= cutoff)
            budget = min(int(record["budget"]), len(genomes))
            scores = {
                "catalyst": -mean,
                "uncertainty": uncertainty,
                "validity": validity_probability,
            }
            round_hits: dict[str, int] = {}
            for name, score in scores.items():
                count = int(hits[_coverage_select(score, genomes, budget)].sum())
                policy_hits[name] += count
                combined_hits[name] += count
                target = frozen_hits if round_name in ("v6", "v7") else development_hits
                target[name] += count
                round_hits[name] = count
            for trial in range(spec.random_trials):
                random_hits[trial] += int(
                    hits[rng.choice(len(genomes), budget, replace=False)].sum()
                )
            selected_count += budget
            combined_selected += budget
            rounds.append(
                {
                    "round": round_name,
                    "eligible": len(genomes),
                    "valid": int(valid.sum()),
                    "budget": budget,
                    "hit_cutoff": cutoff,
                    "hits": round_hits,
                }
            )
            training = _clean(
                pd.concat([training, observed], ignore_index=True, sort=False), outcome
            )

        policies = {
            name: asdict(
                PolicyResult(
                    hits,
                    selected_count,
                    hits / selected_count,
                    development_hits[name],
                    frozen_hits[name],
                )
            )
            for name, hits in policy_hits.items()
        }
        random_upper = float(np.quantile(random_hits, 0.975))
        applications.append(
            {
                "application": application,
                "rounds": rounds,
                "selected": selected_count,
                "policies": policies,
                "random": {
                    "trials": spec.random_trials,
                    "mean_hits": float(random_hits.mean()),
                    "hits_95pct": [
                        float(np.quantile(random_hits, 0.025)),
                        random_upper,
                    ],
                },
                "catalyst_beats_all_baselines": bool(
                    policy_hits["catalyst"] > random_upper
                    and policy_hits["catalyst"] > policy_hits["uncertainty"]
                    and policy_hits["catalyst"] > policy_hits["validity"]
                ),
            }
        )
        combined_random = (
            random_hits if combined_random is None else combined_random + random_hits
        )

    assert combined_random is not None
    combined_upper = float(np.quantile(combined_random, 0.975))
    return {
        "schema_version": 1,
        "method": "chronological finished-candidate walk-forward",
        "hit_definition": "valid candidate in best 20% of its locked round",
        "applications": applications,
        "combined": {
            "selected": combined_selected,
            "hits": combined_hits,
            "random_mean_hits": float(combined_random.mean()),
            "random_hits_95pct": [
                float(np.quantile(combined_random, 0.025)),
                combined_upper,
            ],
            "catalyst_beats_all_baselines": bool(
                combined_hits["catalyst"] > combined_upper
                and combined_hits["catalyst"] > combined_hits["uncertainty"]
                and combined_hits["catalyst"] > combined_hits["validity"]
            ),
        },
    }


def write_search_policy_benchmark(result: dict, output: str | Path) -> Path:
    """Write a benchmark result as stable, reviewable JSON.

    Args:
        result: Completed benchmark metrics to preserve.
        output: Destination JSON path.

    Returns:
        Path to the written evidence artifact.
    """
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return path
