"""Named evidence, selection, campaign-status, and lineage documents."""

from __future__ import annotations

from typing import Any, Literal, TypedDict


CandidateDisposition = Literal[
    'hard_excluded', 'validation_required', 'quantitative_screening']


class CampaignStatus(TypedDict):
    """Fail-closed summary of campaign acceptance criteria."""

    ready: bool
    criteria: dict[str, bool]
    missing: list[str]
    evidence_manifest_valid: bool
    evidence_errors: list[str]


class LedgerEvent(TypedDict):
    """One hash-linked, JSON-compatible multi-fidelity campaign event."""

    sequence: int
    event_type: str
    previous_event_sha256: str | None
    payload: dict[str, Any]
    event_sha256: str


class LedgerDocument(TypedDict):
    """Versioned persisted campaign event chain."""

    schema_version: int
    events: list[LedgerEvent]
