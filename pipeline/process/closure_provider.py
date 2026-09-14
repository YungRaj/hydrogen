"""Select full-physics or calibrated-surrogate reactor transport closures."""

from __future__ import annotations

import math
from typing import Mapping


def resolve_reactor_closure(*, full_physics: Mapping, surrogate,
                            features: Mapping[str, float] | None,
                            pathway_mode: str, reactor_type: str,
                            temperature_K: float) -> dict:
    """Prefer validated physics, otherwise admit only a safe surrogate closure.

    Args:
        full_physics: Mapping supplying full physics.
        surrogate: Surrogate used by this operation.
        features: Numerical model features.
        pathway_mode: Configured methane-conversion pathway.
        reactor_type: Physical reactor implementation identifier.
        temperature_K: Absolute temperature in kelvin.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    if full_physics.get('valid') is True:
        return {
            'available': True, 'source': 'validated_full_physics',
            'outputs': dict(full_physics['artifact']['outputs']),
            'evidence': full_physics,
            'candidate_exclusion_authorized': False,
        }
    if surrogate is None:
        return {
            'available': False, 'source': 'full_physics_required',
            'reason': 'transport_surrogate_unavailable',
            'full_physics': dict(full_physics),
            'candidate_exclusion_authorized': False,
        }
    if (surrogate.pathway_mode != pathway_mode or
            surrogate.reactor_type != reactor_type):
        return {
            'available': False, 'source': 'full_physics_required',
            'reason': 'transport_surrogate_identity_mismatch',
            'full_physics': dict(full_physics),
            'candidate_exclusion_authorized': False,
        }
    if not features:
        return {
            'available': False, 'source': 'full_physics_required',
            'reason': 'transport_closure_features_missing',
            'full_physics': dict(full_physics),
            'candidate_exclusion_authorized': False,
        }
    try:
        feature_temperature = float(features['operating.temperature_K'])
        temperature_matches = math.isclose(
            feature_temperature, float(temperature_K), abs_tol=1e-6)
    except (KeyError, TypeError, ValueError):
        temperature_matches = False
    if not temperature_matches:
        return {
            'available': False, 'source': 'full_physics_required',
            'reason': 'transport_closure_temperature_mismatch',
            'full_physics': dict(full_physics),
            'candidate_exclusion_authorized': False,
        }
    prediction = surrogate.predict(features)
    if prediction.get('usable') is not True:
        return {
            'available': False, 'source': 'full_physics_required',
            'reason': prediction.get('reason', 'transport_surrogate_rejected'),
            'surrogate_prediction': prediction,
            'full_physics': dict(full_physics),
            'candidate_exclusion_authorized': False,
        }
    return {
        'available': True, 'source': 'calibrated_transport_surrogate',
        'outputs': dict(prediction['predictions']),
        'uncertainty_1sigma': dict(prediction['uncertainty_1sigma']),
        'validation_rmse': dict(prediction['validation_rmse']),
        'model_sha256': surrogate.sha256(),
        'training_case_ids': list(surrogate.training_case_ids),
        'validation_case_ids': list(surrogate.validation_case_ids),
        'candidate_exclusion_authorized': False,
    }
