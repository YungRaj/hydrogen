"""Portable, checksum-verified registry for transport-closure surrogates."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

from pipeline.process.multifidelity_surrogate import TransportSurrogate


REGISTRY_SCHEMA_VERSION = 1
_SAFE_NAME = re.compile(r'^[A-Za-z0-9_.-]+$')


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


class TransportModelRegistry:
    """Publish and retrieve models without assuming a machine directory layout."""

    def __init__(self, root: str | Path):
        self.root = Path(root).expanduser().resolve()

    def _stem(self, mode: str, reactor: str) -> str:
        if not _SAFE_NAME.fullmatch(mode) or not _SAFE_NAME.fullmatch(reactor):
            raise ValueError('mode and reactor must be safe path components')
        return f'{mode}__{reactor}'

    def publish(self, model: TransportSurrogate) -> Path:
        """Atomically publish a model and a manifest tied to its exact bytes."""
        self.root.mkdir(parents=True, exist_ok=True)
        stem = self._stem(model.pathway_mode, model.reactor_type)
        target = self.root / f'{stem}.model.json'
        temporary = self.root / f'{stem}.model.json.tmp'
        temporary.write_text(
            json.dumps(model.to_dict(), indent=2, sort_keys=True) + '\n')
        temporary.replace(target)
        manifest = {
            'schema_version': REGISTRY_SCHEMA_VERSION,
            'pathway_mode': model.pathway_mode,
            'reactor_type': model.reactor_type,
            'model_file': target.name,
            'model_file_sha256': _sha256(target),
            'model_content_sha256': model.sha256(),
            'training_case_ids': list(model.training_case_ids),
            'validation_case_ids': list(model.validation_case_ids),
            'candidate_exclusion_authorized': False,
        }
        manifest_target = self.root / f'{stem}.manifest.json'
        manifest_tmp = self.root / f'{stem}.manifest.json.tmp'
        manifest_tmp.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + '\n')
        manifest_tmp.replace(manifest_target)
        return manifest_target

    def load(self, mode: str, reactor: str) -> TransportSurrogate | None:
        """Return a verified identity-matched model, or None when unpublished."""
        stem = self._stem(mode, reactor)
        manifest_path = self.root / f'{stem}.manifest.json'
        if not manifest_path.is_file():
            return None
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError('transport model manifest is unreadable') from exc
        if (manifest.get('schema_version') != REGISTRY_SCHEMA_VERSION or
                manifest.get('pathway_mode') != mode or
                manifest.get('reactor_type') != reactor or
                manifest.get('candidate_exclusion_authorized') is not False):
            raise ValueError('transport model manifest identity is invalid')
        model_name = manifest.get('model_file')
        if not isinstance(model_name, str) or Path(model_name).name != model_name:
            raise ValueError('transport model manifest path is unsafe')
        model_path = self.root / model_name
        if (not model_path.is_file() or
                _sha256(model_path) != manifest.get('model_file_sha256')):
            raise ValueError('transport model file checksum mismatch')
        model = TransportSurrogate.load(model_path)
        if (model.pathway_mode != mode or model.reactor_type != reactor or
                model.sha256() != manifest.get('model_content_sha256') or
                list(model.training_case_ids) != manifest.get('training_case_ids') or
                list(model.validation_case_ids) != manifest.get('validation_case_ids')):
            raise ValueError('transport model content does not match manifest')
        return model
