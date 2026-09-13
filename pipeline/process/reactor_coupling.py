"""Validated adapter from multiphysics evidence to reduced reactor inputs."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable


_CLOSURE_FIELDS = {
    'Fluidized': {
        'gas_velocity_m_s': 'gas_velocity_m_s',
        'u_mf_m_s': 'u_mf_m_s',
        'bubble_fraction': 'fluidized_bubble_fraction',
    },
    'MMBCR': {
        'gas_velocity_m_s': 'gas_velocity_m_s',
        'gas_holdup_fraction': 'gas_holdup_fraction',
        'bubble_diameter_mm': 'bubble_diameter_mm',
    },
    'NTEC': {},
    'Electrochemical': {},
}


@dataclass(frozen=True)
class ReactorCouplingServices:
    """Replaceable evidence-loading, compatibility, and closure operations."""

    load_artifact: Callable
    validate_compatibility: Callable
    couple_evidence: Callable
    load_surrogate: Callable | None = None


def validate_evidence_compatibility(config, loaded: dict) -> dict:
    """Apply mode-specific compatibility checks after generic validation."""
    if loaded.get('valid') and config.reactor_type == 'Electrochemical':
        from pipeline.process.electrochemical_model import (
            conditions_from_environment)
        requested = conditions_from_environment()
        artifact_phase = loaded['artifact'].get('electrolyte_phase')
        if (requested.electrolyte_phase is not None and
                requested.electrolyte_phase.lower() != artifact_phase):
            return {
                'valid': False, 'reason': 'electrolyte_phase_mismatch',
                'requested': requested.electrolyte_phase.lower(),
                'artifact': artifact_phase, 'path': loaded['path']}
    return loaded


def default_reactor_coupling_services() -> ReactorCouplingServices:
    from pipeline.process.multiphysics_contract import load_validated_artifact
    return ReactorCouplingServices(
        load_artifact=load_validated_artifact,
        validate_compatibility=validate_evidence_compatibility,
        couple_evidence=couple_multiphysics_evidence)


def couple_multiphysics_evidence(config, loaded: dict):
    """Attach identity-matched evidence and map established closure fields.

    ``loaded`` must be the result of the strict artifact loader. The identity is
    rechecked here so this adapter remains safe when used independently.
    """
    if loaded.get('valid') is not True:
        raise ValueError('cannot couple invalid multiphysics evidence')
    artifact = loaded.get('artifact', {})
    expected = {
        'candidate_id': str(config.candidate_id),
        'pathway_mode': config.pathway_mode,
        'reactor_type': config.reactor_type,
    }
    failures = [
        name for name, value in expected.items()
        if artifact.get(name) != value]
    try:
        if not math.isclose(float(artifact.get('temperature_K')),
                            float(config.T_inlet_K), abs_tol=1e-6):
            failures.append('temperature_K')
    except (TypeError, ValueError):
        failures.append('temperature_K')
    if failures:
        raise ValueError(
            'multiphysics evidence identity mismatch: ' +
            ', '.join(sorted(set(failures))))
    fields = _CLOSURE_FIELDS.get(config.reactor_type)
    if fields is None:
        raise ValueError(f'no coupling contract for {config.reactor_type}')
    _apply_closure_outputs(config, artifact.get('outputs', {}))
    config.multiphysics_artifact = loaded
    config.reactor_closure_evidence = {
        'source': 'validated_full_physics',
        'candidate_exclusion_authorized': False}
    return config


def _apply_closure_outputs(config, outputs: dict) -> None:
    """Map established closure scalars without assigning evidence identity."""
    fields = _CLOSURE_FIELDS.get(config.reactor_type)
    if fields is None:
        raise ValueError(f'no coupling contract for {config.reactor_type}')
    for output_name, config_name in fields.items():
        try:
            value = float(outputs[output_name])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f'invalid coupling output: {output_name}') from exc
        if not math.isfinite(value):
            raise ValueError(f'invalid coupling output: {output_name}')
        setattr(config, config_name, value)


def couple_surrogate_closure(config, decision: dict):
    """Apply an explicitly calibrated surrogate closure without relabeling it."""
    if (decision.get('available') is not True or
            decision.get('source') != 'calibrated_transport_surrogate' or
            decision.get('candidate_exclusion_authorized') is not False):
        raise ValueError('cannot couple an unaccepted surrogate closure')
    _apply_closure_outputs(config, decision.get('outputs', {}))
    config.reactor_closure_evidence = decision
    return config
