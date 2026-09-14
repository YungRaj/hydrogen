"""Atomic, hash-chained lineage for multi-fidelity campaign iterations."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Mapping


LEDGER_SCHEMA_VERSION = 1


def _digest(value: Mapping) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(',', ':')).encode()
    return hashlib.sha256(payload).hexdigest()


class CampaignLedger:
    """Maintain an atomic, tamper-evident campaign event chain.
    """
    def __init__(self, path: str | Path):
        self.path = Path(path)

    def read_verified(self) -> list[dict]:
        """Read the campaign ledger and verify every hash-chain link.

        Returns:
            A list of ledger events after schema and hash-chain verification.
        """
        if not self.path.is_file():
            return []
        try:
            value = json.loads(self.path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError('campaign ledger is unreadable') from exc
        if value.get('schema_version') != LEDGER_SCHEMA_VERSION:
            raise ValueError('campaign ledger schema is unsupported')
        events = value.get('events')
        if not isinstance(events, list):
            raise ValueError('campaign ledger events are invalid')
        previous = None
        for index, event in enumerate(events):
            if not isinstance(event, dict):
                raise ValueError('campaign ledger event is invalid')
            claimed = event.get('event_sha256')
            body = {key: val for key, val in event.items()
                    if key != 'event_sha256'}
            if (body.get('sequence') != index or
                    body.get('previous_event_sha256') != previous or
                    claimed != _digest(body)):
                raise ValueError('campaign ledger hash chain is invalid')
            previous = claimed
        return events

    def append(self, event_type: str, payload: Mapping) -> dict:
        """Append an event and atomically persist the updated hash chain.

        Args:
            event_type: Stable category assigned to the new ledger event.
            payload: JSON-compatible evidence stored in the event.

        Returns:
            The newly appended event, including its sequence and SHA-256 digest.
        """
        if not event_type or not isinstance(payload, Mapping):
            raise ValueError('ledger event type and payload are required')
        events = self.read_verified()
        body = {
            'sequence': len(events),
            'event_type': event_type,
            'previous_event_sha256': (
                events[-1]['event_sha256'] if events else None),
            'payload': dict(payload),
        }
        event = {**body, 'event_sha256': _digest(body)}
        document = {'schema_version': LEDGER_SCHEMA_VERSION,
                    'events': [*events, event]}
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(self.path.suffix + '.tmp')
        temporary.write_text(json.dumps(document, indent=2, sort_keys=True) + '\n')
        temporary.replace(self.path)
        return event
