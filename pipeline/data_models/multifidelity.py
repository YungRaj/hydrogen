"""Typed plans, queries, decisions, and summaries for multi-fidelity campaigns."""

from __future__ import annotations

from typing import Literal, TypedDict

from pipeline.data_models.transport import TransportPrediction


CasePartition = Literal['training', 'validation']


class RepresentativeCase(TypedDict):
    """Solver-independent point in a representative physical-case design."""

    case_id: str
    pathway_mode: str
    reactor_type: str
    features: dict[str, float]
    design_role: Literal['regime_anchor', 'stratified_coverage']


class _DesignedCaseOptional(TypedDict, total=False):
    pathway_mode: str
    reactor_type: str
    temperature_K: float
    case_dir: str
    model_source: str
    fenics_model: str
    timeout_s: int
    max_coupling_iterations: int


class DesignedCase(_DesignedCaseOptional):
    """Preassigned full-physics execution plan consumed by a campaign."""

    case_id: str
    partition: CasePartition


class _ScreeningQueryOptional(TypedDict, total=False):
    expected_improvement: float
    uncertainty: float
    calibration_error: float
    disagreement: float
    unproductive_history: float


class ScreeningQuery(_ScreeningQueryOptional):
    """One cheap inference request and its adaptive-allocation signals."""

    query_id: str
    region: str
    features: dict[str, float]


class ScreeningDecision(TransportPrediction):
    """Traceable surrogate/referral decision for one screening query."""

    query_id: str
    region: str
    input_sha256: str


class FullPhysicsSelection(TypedDict):
    """Budgeted referral preserving regional coverage."""

    case_id: str
    region: str
    priority: float
    selection_reason: Literal['fixed_regional_coverage', 'adaptive_priority']
    candidate_exclusion_authorized: Literal[False]


class MultiFidelityIterationResult(TypedDict):
    """Complete lineage-ready outcome of one active-learning iteration."""

    pathway_mode: str
    reactor_type: str
    model_sha256: str
    executed_case_ids: list[str]
    prior_case_count: int
    training_case_ids: list[str]
    validation_case_ids: list[str]
    screening_decisions: list[ScreeningDecision]
    accepted_surrogate_count: int
    full_physics_required_count: int
    scheduled_referrals: list[FullPhysicsSelection]
    candidate_exclusion_authorized: Literal[False]
    lineage_event_sha256: str | None
