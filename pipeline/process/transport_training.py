"""Fail-closed workflow from validated full-physics cases to a published model."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping

from pipeline.process.multifidelity_surrogate import (
    fit_transport_surrogate, record_from_artifact)
from pipeline.process.multiphysics_contract import load_validated_artifact
from pipeline.process.transport_model_registry import TransportModelRegistry


@dataclass(frozen=True)
class CaseArtifactReference:
    """An explicit solver artifact and its preassigned data partition."""

    candidate_id: str
    pathway_mode: str
    reactor_type: str
    temperature_K: float
    partition: str

    def validate(self) -> None:
        """Reject invalid configuration before it reaches scientific execution.
        """
        if self.partition not in {'training', 'validation'}:
            raise ValueError('case partition must be training or validation')
        if not self.candidate_id or not self.pathway_mode or not self.reactor_type:
            raise ValueError('case artifact identity fields are required')


def train_and_publish_transport_model(
        references: Iterable[CaseArtifactReference], *, results_dir: str | Path,
        model_registry: TransportModelRegistry, targets: Iterable[str],
        validation_rmse_limits: Mapping[str, float], ensemble_size: int = 8,
        ridge: float = 1e-8, random_seed: int = 0,
        artifact_loader: Callable = load_validated_artifact,
        model_fitter: Callable = fit_transport_surrogate) -> dict:
    """Validate every declared label, fit a blind holdout, and publish atomically.

        Partitions are supplied before fitting rather than selected opportunistically
        after inspecting errors. Any invalid artifact aborts the complete operation;
        the workflow never trains on a silently reduced collection.

    Args:
        references: References used by this operation.
        results_dir: Directory containing or receiving calculation results.
        model_registry: Model registry used by this operation.
        targets: Output properties learned or evaluated by the model.
        validation_rmse_limits: Mapping supplying validation rmse limits.
        ensemble_size: Number of ensemble size to use.
        ridge: Ridge used by this operation.
        random_seed: Seed controlling deterministic sampling or fitting.
        artifact_loader: Injected callable used to perform artifact loader.
        model_fitter: Injected callable used to perform model fitter.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    refs = list(references)
    if not refs:
        raise ValueError('at least one artifact reference is required')
    identity = set()
    records = {'training': [], 'validation': []}
    for reference in refs:
        if not isinstance(reference, CaseArtifactReference):
            raise TypeError('references must be CaseArtifactReference values')
        reference.validate()
        key = (reference.candidate_id, reference.pathway_mode,
               reference.reactor_type, float(reference.temperature_K))
        if key in identity:
            raise ValueError('duplicate physical case reference')
        identity.add(key)
        loaded = artifact_loader(
            results_dir, reference.candidate_id, reference.pathway_mode,
            reference.reactor_type, reference.temperature_K)
        if loaded.get('valid') is not True:
            raise ValueError(
                f'invalid full-physics artifact for {reference.candidate_id}: '
                f"{loaded.get('reason', 'validation_failed')}")
        record = record_from_artifact(loaded['artifact'])
        if ((record.pathway_mode, record.reactor_type) !=
                (reference.pathway_mode, reference.reactor_type)):
            raise ValueError('validated artifact identity differs from its reference')
        records[reference.partition].append(record)
    model = model_fitter(
        records['training'], records['validation'], targets=tuple(targets),
        validation_rmse_limits=validation_rmse_limits,
        ensemble_size=ensemble_size, ridge=ridge, random_seed=random_seed)
    manifest = model_registry.publish(model)
    return {
        'status': 'published',
        'pathway_mode': model.pathway_mode,
        'reactor_type': model.reactor_type,
        'model_sha256': model.sha256(),
        'manifest_path': str(manifest),
        'training_case_ids': list(model.training_case_ids),
        'validation_case_ids': list(model.validation_case_ids),
        'validation_rmse': dict(zip(
            model.target_names, model.validation_rmse.tolist())),
        'candidate_exclusion_authorized': False,
    }
