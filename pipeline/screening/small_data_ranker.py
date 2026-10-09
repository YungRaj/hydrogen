"""Deterministic screening rankers for small catalyst datasets."""

from __future__ import annotations

import ast
from dataclasses import dataclass

import numpy as np

from pipeline.search.design_space import encode_population

TREE_ENSEMBLE_SIZE = 256
MIN_TRAINING_ROWS = 20
MIN_RANK_CORRELATION = 0.20
CLASS_BIAS_WEIGHT = 0.25
CLASS_BIAS_PRIOR_ROWS = 5
INVALIDITY_PENALTY = 0.75
ORR_UNCERTAINTY_BONUS = 0.25
CLASS_SUCCESS_WEIGHT = 0.75
CLASS_SUCCESS_PRIOR_ROWS = 5
CENSORED_E_ACT_PENALTY_EV = 5.0


def screening_metric_eligibility(frame, application: str) -> np.ndarray:
    """Return rows whose metric may guide discovery or count as a hit.

    Pyrolysis values clipped at the model floor are censored observations, not
    evidence of a low barrier.  Non-viable pyrolysis structures likewise must
    not train or reward the search policy.  Missing gate columns fail closed.

    Args:
        frame: Screening records containing validity and application metrics.
        application: ``turquoise_hydrogen`` or ``fuel_cell_orr``.

    Returns:
        Boolean eligibility mask aligned with ``frame``.
    """
    outcome = {
        "turquoise_hydrogen": "E_act",
        "fuel_cell_orr": "orr_overpotential_V",
    }.get(application)
    if outcome is None:
        raise ValueError(f"unknown application {application}")
    if outcome not in frame.columns or "valid" not in frame.columns:
        return np.zeros(len(frame), dtype=bool)
    import pandas as pd

    values = pd.to_numeric(frame[outcome], errors="coerce").to_numpy(dtype=float)
    eligible = frame["valid"].eq(True).to_numpy() & np.isfinite(values)
    if application == "turquoise_hydrogen":
        if (
            "E_act_censored" not in frame.columns
            or "pyrolysis_viable" not in frame.columns
        ):
            return np.zeros(len(frame), dtype=bool)
        eligible &= frame["E_act_censored"].eq(False).to_numpy()
        eligible &= frame["pyrolysis_viable"].eq(True).to_numpy()
    return eligible


def training_metric_eligibility(frame, application: str) -> np.ndarray:
    """Return rows permitted to train the continuous metric ranker.

    Censored pyrolysis barriers remain useful negative evidence: they identify
    the over-binding region, but their numerical floor is not a measured target.
    Non-viable phases remain excluded because selection cannot act on them.

    Args:
        frame: Screening records containing validity and application metrics.
        application: ``turquoise_hydrogen`` or ``fuel_cell_orr``.

    Returns:
        Boolean training mask aligned with ``frame``.
    """
    outcome = {
        "turquoise_hydrogen": "E_act",
        "fuel_cell_orr": "orr_overpotential_V",
    }.get(application)
    if outcome is None:
        raise ValueError(f"unknown application {application}")
    if outcome not in frame.columns or "valid" not in frame.columns:
        return np.zeros(len(frame), dtype=bool)
    import pandas as pd

    values = pd.to_numeric(frame[outcome], errors="coerce").to_numpy(dtype=float)
    eligible = frame["valid"].eq(True).to_numpy() & np.isfinite(values)
    if application == "turquoise_hydrogen":
        if "pyrolysis_viable" not in frame.columns:
            return np.zeros(len(frame), dtype=bool)
        eligible &= frame["pyrolysis_viable"].eq(True).to_numpy()
    return eligible


@dataclass(frozen=True)
class RankerDiagnostics:
    """Out-of-fold evidence deciding whether predictions may rank candidates."""

    sample_count: int
    held_out_group_count: int
    validation_strategy: str
    spearman_correlation: float
    model_mae: float
    baseline_mae: float
    ranking_validated: bool

    @property
    def acquisition_mode(self) -> str:
        """Return the acquisition behavior authorized by this evidence.

        Returns:
            ``validated_quality`` or conservative catalyst-quality ranking.
        """
        if self.ranking_validated:
            return "validated_quality"
        return "catalyst_quality"


def valid_training_row_count(frame, application: str) -> int:
    """Count rows that can actually train the requested ranker.

    Args:
        frame: Tabular candidate or result records.
        application: Scientific objective, such as pyrolysis or ORR.

    Returns:
        Computed `int` value in the units documented above.
    """
    eligible = training_metric_eligibility(frame, application)
    count = 0
    for position, (_, row) in enumerate(frame.iterrows()):
        if not eligible[position]:
            continue
        try:
            (
                ast.literal_eval(row["genome"])
                if isinstance(row["genome"], str)
                else tuple(row["genome"])
            )
            count += 1
        except (ValueError, TypeError, SyntaxError, KeyError):
            continue
    return count


def merge_compatible_evidence(current, ledger, protocol_id: str):
    """Merge prior same-protocol observations without inventing compatibility.

    Args:
        current: Current used by this operation.
        ledger: Ledger used by this operation.
        protocol_id: Protocol id used by this operation.

    Returns:
        Computed result described above.
    """
    import pandas as pd

    frames = [current]
    if ledger is not None and len(ledger) and 'screening_protocol' in ledger.columns:
        compatible = ledger[ledger['screening_protocol'] == protocol_id]
        if len(compatible):
            frames.insert(0, compatible)
    merged = pd.concat(frames, ignore_index=True)
    if 'genome' in merged.columns:
        merged = merged.drop_duplicates('genome', keep='last').reset_index(drop=True)
    return merged


@dataclass
class TreeRanker:
    """Wrap a fitted small-data ranker and its uncertainty ensemble.

    Attributes:
        application: Configured application value.
        model: Configured model value.
        target_columns: Configured target columns value.
    """

    application: str
    model: object
    target_columns: tuple[str, ...]
    diagnostics: RankerDiagnostics
    class_bias: dict[str, float]

    def _primary(self, raw) -> np.ndarray:
        return _primary_target(raw, self.application)

    def predict(
        self, genomes, *, uncertainty: bool = True
    ) -> tuple[np.ndarray, np.ndarray]:
        # sklearn validates and converts X on every individual DecisionTree
        # call.  The encoder already guarantees a finite dense matrix, so do
        # that conversion once instead of hundreds of times per scan batch.
        """Predict outputs and uncertainty for one feature vector.

        Args:
            genomes: Encoded catalyst candidates to score.
            uncertainty: Whether prediction uncertainty should be returned.

        Returns:
            Predicted outputs paired with model uncertainty.
        """
        x = np.asarray(encode_population(genomes), dtype=np.float32, order="C")
        total = np.zeros(len(x), dtype=float)
        total_sq = np.zeros(len(x), dtype=float) if uncertainty else None
        for tree in self.model.estimators_:
            member = self._primary(tree.predict(x, check_input=False))
            total += member
            if total_sq is not None:
                total_sq += member * member
        count = len(self.model.estimators_)
        mean = total / count
        if total_sq is None:
            uncertainty_values = np.zeros(len(x), dtype=float)
        else:
            variance = np.maximum(0.0, total_sq / count - mean * mean)
            uncertainty_values = np.sqrt(variance)
        mean += CLASS_BIAS_WEIGHT * np.asarray(
            [self.class_bias.get(str(genome[0]), 0.0) for genome in genomes]
        )
        # Class calibration shifts every tree equally and therefore cannot
        # change ensemble variance. Streaming moments avoid a batch-by-tree
        # allocation in every shard.
        return mean, uncertainty_values


def fit_tree_ranker(
    frame, application: str, random_state: int = 20260721
) -> TreeRanker:
    """Fit the application-specific form validated by prospective pilots.

    Args:
        frame: Tabular candidate or result records.
        application: Scientific objective, such as pyrolysis or ORR.
        random_state: Random state used by this operation.

    Returns:
        Computed `TreeRanker` result.
    """
    from scipy.stats import spearmanr
    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.model_selection import GroupKFold, KFold

    if application == "turquoise_hydrogen":
        columns = ("E_act",)
    elif application == "fuel_cell_orr":
        columns = ("orr_overpotential_V",)
    else:
        raise ValueError(f"unknown application {application}")
    rows, genomes = [], []
    eligible = training_metric_eligibility(frame, application)
    for position, (_, row) in enumerate(frame.iterrows()):
        if not eligible[position]:
            continue
        try:
            values = [float(row[column]) for column in columns]
            if (
                application == "turquoise_hydrogen"
                and row.get("E_act_censored", False) == True
            ):
                values[0] = CENSORED_E_ACT_PENALTY_EV
            genome = (
                ast.literal_eval(row["genome"])
                if isinstance(row["genome"], str)
                else tuple(row["genome"])
            )
            rows.append(values)
            genomes.append(genome)
        except (ValueError, TypeError, SyntaxError, KeyError):
            continue
    if len(rows) < MIN_TRAINING_ROWS:
        raise ValueError(
            f"tree ranker requires at least {MIN_TRAINING_ROWS} valid rows; got {len(rows)}"
        )
    y = np.asarray(rows, float)
    y = y[:, 0]
    encoded = encode_population(genomes)
    observed_primary = _primary_target(y, application)
    # Hold out whole material classes. Fine-grained discovery regions are often
    # unique in a small calibration set and would therefore provide no stronger
    # separation than ordinary row-wise folds.
    groups = np.asarray([str(genome[0]) for genome in genomes])
    group_count = len(np.unique(groups))
    out_of_fold = np.empty_like(y)
    baseline_predictions = np.empty(len(observed_primary), dtype=float)
    if group_count >= 2:
        folds = GroupKFold(n_splits=min(5, group_count))
        for train, test in folds.split(encoded, y, groups):
            validation_model = ExtraTreesRegressor(
                n_estimators=TREE_ENSEMBLE_SIZE,
                min_samples_leaf=1,
                max_features=1.0,
                random_state=random_state,
                n_jobs=-1,
            )
            validation_model.fit(encoded[train], y[train])
            out_of_fold[test] = validation_model.predict(encoded[test])
            baseline_predictions[test] = np.median(observed_primary[train])
        predicted_primary = _primary_target(out_of_fold, application)
    else:
        predicted_primary = np.full(len(observed_primary), np.median(observed_primary))
        baseline_predictions.fill(np.median(observed_primary))
    correlation = float(spearmanr(observed_primary, predicted_primary).statistic)
    if not np.isfinite(correlation):
        correlation = 0.0
    model_mae = float(np.mean(np.abs(observed_primary - predicted_primary)))
    baseline_mae = float(np.mean(np.abs(observed_primary - baseline_predictions)))
    diagnostics = RankerDiagnostics(
        sample_count=len(rows),
        held_out_group_count=group_count,
        validation_strategy="held_out_material_class_cv",
        spearman_correlation=correlation,
        model_mae=model_mae,
        baseline_mae=baseline_mae,
        ranking_validated=(
            group_count >= 2
            and correlation >= MIN_RANK_CORRELATION
            and model_mae < baseline_mae
        ),
    )
    # Correct only persistent within-class residuals. Five pseudo-observations
    # at the global residual shrink sparse classes toward no correction.
    calibration_predictions = np.empty_like(y)
    calibration_folds = KFold(n_splits=5, shuffle=True, random_state=77)
    for train, test in calibration_folds.split(encoded):
        calibration_model = ExtraTreesRegressor(
            n_estimators=64,
            min_samples_leaf=2,
            max_features=1.0,
            random_state=21,
            n_jobs=-1,
        )
        calibration_model.fit(encoded[train], y[train])
        calibration_predictions[test] = calibration_model.predict(encoded[test])
    residual = observed_primary - _primary_target(
        calibration_predictions, application
    )
    global_bias = float(np.mean(residual))
    class_bias = {}
    for material_class in np.unique(groups):
        mask = groups == material_class
        class_bias[str(material_class)] = float(
            (residual[mask].sum() + CLASS_BIAS_PRIOR_ROWS * global_bias)
            / (mask.sum() + CLASS_BIAS_PRIOR_ROWS)
        )
    # With tens to hundreds of calibration rows, 1024 trees add repeated
    # inference work without useful independent evidence.  The fixed 256-tree
    # ensemble retains deterministic uncertainty and ranking while keeping a
    # production-scale indexed scan tractable.
    model = ExtraTreesRegressor(
        n_estimators=TREE_ENSEMBLE_SIZE,
        min_samples_leaf=1,
        max_features=1.0,
        random_state=random_state,
        n_jobs=-1,
    )
    model.fit(encoded, y)
    return TreeRanker(application, model, columns, diagnostics, class_bias)


def _primary_target(raw, application: str) -> np.ndarray:
    """Convert fitted targets to the scalar quantity used for candidate ranking."""
    values = np.asarray(raw, float)
    if values.ndim == 1 or values.shape[1] == 1:
        return values.reshape(-1)
    d_oh, d_o, d_ooh = values.T
    # Keep the continuous CHE value: clipping destroys rank information.
    return 1.23 + np.maximum.reduce([d_ooh - 4.92, d_o - d_ooh, d_oh - d_o, -d_oh])


def turquoise_tree_objectives(genomes, ranker: TreeRanker) -> np.ndarray:
    """Extract turquoise-hydrogen ranking targets from a screening frame.

    Args:
        genomes: Encoded catalyst candidates to score.
        ranker: Fitted tree ranker used for objective prediction.

    Returns:
        A `np.ndarray` containing the turquoise tree objectives result.
    """
    from pipeline.utils import abundance_cost_penalty
    from pipeline.screening.genetic_optimizer import _extract_elements_from_genome

    primary, _ = ranker.predict(genomes, uncertainty=False)
    costs = (
        [
            -abundance_cost_penalty(_extract_elements_from_genome(genome))
            for genome in genomes
        ]
        if ranker.diagnostics.ranking_validated
        else np.zeros(len(genomes))
    )
    return np.column_stack(
        [primary, np.zeros(len(genomes)), np.zeros(len(genomes)), costs]
    )


def orr_tree_objectives(genomes, ranker: TreeRanker) -> np.ndarray:
    """Extract ORR ranking targets from a screening frame.

    Args:
        genomes: Encoded catalyst candidates to score.
        ranker: Fitted tree ranker used for objective prediction.

    Returns:
        A `np.ndarray` containing the orr tree objectives result.
    """
    from pipeline.screening.fc_genetic_optimizer import (
        _cost_from_genome,
        _fenton_from_genome,
    )
    from pipeline.search.scope import pemfc_cathode_scope

    primary, _ = ranker.predict(genomes, uncertainty=False)
    objectives = np.column_stack(
        [
            primary,
            (
                [-_fenton_from_genome(g) for g in genomes]
                if ranker.diagnostics.ranking_validated
                else np.zeros(len(genomes))
            ),
            (
                [_cost_from_genome(g) for g in genomes]
                if ranker.diagnostics.ranking_validated
                else np.zeros(len(genomes))
            ),
            np.zeros(len(genomes)),
        ]
    )
    for i, genome in enumerate(genomes):
        if pemfc_cathode_scope(genome)["status"] != "candidate":
            objectives[i] = [5.0, 0.0, 100.0, 0.0]
    return objectives


def orr_catalyst_acquisition(
    predicted_overpotential: np.ndarray,
    uncertainty: np.ndarray,
    validity_probability: np.ndarray,
) -> np.ndarray:
    """Score ORR candidates by quality, validity, and bounded exploration.

    Higher values are preferred. Coefficients are fixed in volts and preserve
    every historical v2-v8 hit count while improving the first unified batch.

    Args:
        predicted_overpotential: Predicted ORR overpotential in volts.
        uncertainty: Ensemble standard deviation in volts.
        validity_probability: Probability that screening will produce a valid result.

    Returns:
        One finite acquisition score per candidate; higher values rank first.
    """
    predicted = np.asarray(predicted_overpotential, dtype=float)
    spread = np.asarray(uncertainty, dtype=float)
    validity = np.asarray(validity_probability, dtype=float)
    if predicted.shape != spread.shape or predicted.shape != validity.shape:
        raise ValueError("ORR acquisition inputs must have identical shapes")
    if np.any(~np.isfinite(predicted)) or np.any(~np.isfinite(spread)):
        raise ValueError("ORR acquisition inputs must be finite")
    if np.any((validity < 0.0) | (validity > 1.0)):
        raise ValueError("ORR validity probabilities must be in [0, 1]")
    return (
        -predicted
        - INVALIDITY_PENALTY * (1.0 - validity)
        + ORR_UNCERTAINTY_BONUS * spread
    )


def turquoise_catalyst_acquisition(
    predicted_activation_energy: np.ndarray,
    validity_probability: np.ndarray,
) -> np.ndarray:
    """Score low-barrier candidates while penalizing failed evaluations.

    Args:
        predicted_activation_energy: Predicted activation energy in electronvolts.
        validity_probability: Probability that screening produces a valid result.

    Returns:
        One finite acquisition score per candidate; higher values rank first.
    """
    predicted = np.asarray(predicted_activation_energy, dtype=float)
    validity = np.asarray(validity_probability, dtype=float)
    if predicted.shape != validity.shape:
        raise ValueError("turquoise acquisition inputs must have identical shapes")
    if np.any(~np.isfinite(predicted)):
        raise ValueError("turquoise activation predictions must be finite")
    if np.any((validity < 0.0) | (validity > 1.0)):
        raise ValueError("turquoise validity probabilities must be in [0, 1]")
    return -predicted - INVALIDITY_PENALTY * (1.0 - validity)


def class_success_probability(
    frame, genomes, application: str
) -> np.ndarray:
    """Estimate smoothed catalyst-hit probability for each material class.

    Args:
        frame: Completed screening records available before selection.
        genomes: Candidate genomes receiving class feedback.
        application: ``turquoise_hydrogen`` or ``fuel_cell_orr``.

    Returns:
        One probability per candidate, estimated only from completed outcomes.
    """
    outcome = {
        "turquoise_hydrogen": "E_act",
        "fuel_cell_orr": "orr_overpotential_V",
    }.get(application)
    if outcome is None:
        raise ValueError(f"unknown application {application}")
    values = np.asarray(frame[outcome], dtype=float)
    valid = screening_metric_eligibility(frame, application)
    if not np.any(valid):
        raise ValueError("class success feedback requires valid outcomes")
    threshold = float(np.quantile(values[valid], 0.20))
    hits = valid & (np.nan_to_num(values, nan=np.inf) <= threshold)
    classes = np.asarray(
        [
            str(ast.literal_eval(raw)[0]) if isinstance(raw, str) else str(raw[0])
            for raw in frame["genome"]
        ]
    )
    global_rate = float(
        (hits.sum() + 2.0) / (len(hits) + 2 * CLASS_SUCCESS_PRIOR_ROWS)
    )
    rates = {}
    for material_class in np.unique(classes):
        mask = classes == material_class
        rates[material_class] = float(
            (hits[mask].sum() + CLASS_SUCCESS_PRIOR_ROWS * global_rate)
            / (mask.sum() + CLASS_SUCCESS_PRIOR_ROWS)
        )
    return np.asarray(
        [rates.get(str(genome[0]), global_rate) for genome in genomes], dtype=float
    )


def guide_branch_objectives(
    objectives: np.ndarray, frame, genomes, application: str
) -> np.ndarray:
    """Apply completed-candidate class feedback to branch scheduling only.

    Args:
        objectives: Physical candidate objectives with lower primary values preferred.
        frame: Completed screening records available before branch scheduling.
        genomes: Candidate genomes represented by the objective rows.
        application: ``turquoise_hydrogen`` or ``fuel_cell_orr``.

    Returns:
        A copied objective matrix whose primary column includes branch feedback.
    """
    guided = np.asarray(objectives, dtype=float).copy()
    if guided.ndim != 2 or len(guided) != len(genomes) or guided.shape[1] == 0:
        raise ValueError("branch objectives must have shape (N, M) with M >= 1")
    guided[:, 0] -= CLASS_SUCCESS_WEIGHT * class_success_probability(
        frame, genomes, application
    )
    return guided
