"""Shared file-integrity helpers."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


def file_sha256(path: str | Path) -> str:
    """Return the SHA-256 digest for ``path``.

    Args:
        path: File whose bytes should be hashed.

    Returns:
        Lowercase hexadecimal SHA-256 digest.
    """
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_json_sha256(payload: Mapping[str, Any]) -> str:
    """Return the SHA-256 digest of a mapping serialized as canonical sorted JSON.

    Args:
        payload: JSON-serializable mapping to hash.

    Returns:
        Lowercase hexadecimal SHA-256 digest of canonical JSON bytes.
    """
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
    return hashlib.sha256(encoded).hexdigest()
