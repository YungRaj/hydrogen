"""Typed data exchanged by full-physics and transport-surrogate components."""

from __future__ import annotations

from typing import Literal, TypedDict


SurrogateDecision = Literal['surrogate_closure', 'full_physics_required']
ClosureSource = Literal[
    'validated_full_physics', 'calibrated_transport_surrogate',
    'full_physics_required']


class _TransportPredictionOptional(TypedDict, total=False):
    predictions: dict[str, float]
    uncertainty_1sigma: dict[str, float]
    validation_rmse: dict[str, float]


class TransportPrediction(_TransportPredictionOptional):
    """Fail-closed prediction returned by one calibrated transport model."""

    usable: bool
    decision: SurrogateDecision
    reason: str | None
    candidate_exclusion_authorized: Literal[False]


class TransportModelDocument(TypedDict):
    """JSON representation of a fitted and holdout-validated surrogate."""

    schema_version: int
    pathway_mode: str
    reactor_type: str
    feature_names: list[str]
    target_names: list[str]
    center: list[float]
    scale: list[float]
    feature_min: list[float]
    feature_max: list[float]
    coefficients: list[list[list[float]]]
    validation_rmse: list[float]
    uncertainty_limits: list[float]
    training_case_ids: list[str]
    validation_case_ids: list[str]
    scope: Literal['transport_closure_only']
    candidate_exclusion_authorized: Literal[False]


class _ReactorClosureOptional(TypedDict, total=False):
    reason: str
    outputs: dict[str, float]
    uncertainty_1sigma: dict[str, float]
    validation_rmse: dict[str, float]
    model_sha256: str
    training_case_ids: list[str]
    validation_case_ids: list[str]
    surrogate_prediction: TransportPrediction
    full_physics: dict[str, object]
    evidence: dict[str, object]


class ReactorClosure(_ReactorClosureOptional):
    """Transport closure selected for a reduced Cantera reactor calculation."""

    available: bool
    source: ClosureSource
    candidate_exclusion_authorized: Literal[False]


class TransportTrainingResult(TypedDict):
    """Published surrogate identity and its disjoint validation evidence."""

    status: Literal['published']
    pathway_mode: str
    reactor_type: str
    model_sha256: str
    manifest_path: str
    training_case_ids: list[str]
    validation_case_ids: list[str]
    validation_rmse: dict[str, float]
    candidate_exclusion_authorized: Literal[False]
