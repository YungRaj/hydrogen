"""Shared file-integrity helpers."""

from __future__ import annotations

import hashlib
from pathlib import Path


def file_sha256(path: str | Path) -> str:
    """Return the SHA-256 digest for ``path``.

    Args:
        path: File whose bytes should be hashed.

    Returns:
        Lowercase hexadecimal SHA-256 digest.
    """
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
