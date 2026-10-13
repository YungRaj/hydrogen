"""Small atomic JSON persistence primitive."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any, Mapping


def write_atomic_json(path: str | Path, payload: Mapping[str, Any]) -> None:
    """Write sorted, indented JSON and replace the destination atomically.

    Args:
        path: Destination JSON path.
        payload: JSON-serializable mapping to write.
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f'.{destination.name}.', dir=destination.parent
    )
    try:
        with os.fdopen(descriptor, 'w') as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
