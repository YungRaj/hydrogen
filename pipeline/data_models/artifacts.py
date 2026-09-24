"""Typed representations of validated JSON artifacts and solver handoffs."""

from __future__ import annotations

from typing import TypedDict

from pipeline.reactors.modes import PathwayModeName, ReactorTypeName


class FeedDocument(TypedDict):
    """Species composition and its provenance source."""

    composition: dict[str, float]
    source: str


class KineticsDocument(TypedDict):
    """Declared source of candidate-specific reaction kinetics."""

    source: str


class _CalibrationOptional(TypedDict, total=False):
    paired_control: bool


class CalibrationDocument(_CalibrationOptional):
    """Disjoint calibration and holdout definition for a physical case."""

    training_ids: list[str]
    validation_ids: list[str]
    source: str
    metric: str
    acceptance_threshold: float


class _PhysicalCaseOptional(TypedDict, total=False):
    electrolyte_phase: str
    electrolyte: dict[str, str]


class PhysicalCaseDocument(_PhysicalCaseOptional):
    """Validated, unit-bearing external reactor input document."""

    schema_version: int
    template: bool
    candidate_id: str
    pathway_mode: PathwayModeName
    reactor_type: ReactorTypeName
    geometry: dict[str, float]
    operating: dict[str, float]
    properties: dict[str, float]
    models: dict[str, str]
    feed: FeedDocument
    kinetics: KineticsDocument
    parameter_sources: dict[str, str]
    calibration: CalibrationDocument


class PhysicalCaseSummary(TypedDict):
    """Compact provenance retained in a completed solver artifact."""

    schema_version: int
    kinetics_source: str
    feed_source: str
    calibration_source: str
    calibration_count: int
    holdout_validation_count: int
    disjoint_holdout: bool


class SurrogateInputDocument(TypedDict):
    """Numeric physical-case snapshot used to train transport surrogates."""

    schema_version: int
    reactor_type: ReactorTypeName
    units_in_field_names: bool
    values: dict[str, float]
