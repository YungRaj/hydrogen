"""NTEC and electrochemical reactor result adapters.

These modes consume validated multiphysics outputs.  They deliberately return
non-excluding validation requirements when those outputs are unavailable.
"""

from typing import Protocol

from pipeline.data_models.reactors import ReactorResult


class ElectrochemicalReactorConfig(Protocol):
    """Configuration fields required by electrochemical pathway adapters."""

    multiphysics_artifact: dict | None
    catalyst_name: str
    candidate_id: str
    T_inlet_K: float


def simulate_ntec_pathway(config: ElectrochemicalReactorConfig) -> ReactorResult:
    """Build an NTEC result from validated coupled-physics evidence.

    Args:
        config: Reactor configuration containing identity and optional evidence.

    Returns:
        A complete result or a non-excluding validation requirement.
    """
    from pipeline.electrochemistry.ntec import (
        conditions_from_environment, ntec_assistance)

    assistance = ntec_assistance(conditions_from_environment())
    if config.multiphysics_artifact:
        evidence = config.multiphysics_artifact
        artifact = evidence['artifact']
        outputs = artifact['outputs']
        return {
            'status': 'complete', 'valid': True,
            'reactor_type': 'NTEC', 'pathway_mode': 'ntec',
            'catalyst_name': config.catalyst_name,
            'candidate_id': config.candidate_id,
            'T_K': config.T_inlet_K,
            'CH4_conversion': float(outputs['CH4_conversion']),
            'H2_selectivity': float(outputs['H2_selectivity']),
            'solid_C_selectivity': float(outputs['solid_C_selectivity']),
            'specific_energy_kWh_kg_H2': float(
                outputs['specific_energy_kWh_kg_H2']),
            'ntec_assistance': assistance,
            'multiphysics_evidence': evidence,
            'reactor_evidence_tier': 'calibrated_multiphysics_screening',
            'can_exclude_candidate': False,
        }
    return {
        'status': 'validation_required', 'valid': False,
        'reactor_type': 'NTEC', 'pathway_mode': 'ntec',
        'catalyst_name': config.catalyst_name, 'T_K': config.T_inlet_K,
        'ntec_assistance': assistance,
        'reactor_evidence_tier': 'pathway_model_pending',
        'can_exclude_candidate': False,
        'limitations': [
            'candidate_specific_ntec_pathway_kinetics_required',
            'liquid_solid_hydrodynamic_model_required',
        ],
    }


def simulate_electrochemical_pathway(
        config: ElectrochemicalReactorConfig) -> ReactorResult:
    """Build an electrochemical result without inventing conversion evidence.

    Args:
        config: Reactor configuration containing identity and optional evidence.

    Returns:
        A complete result or a non-excluding validation requirement.
    """
    from pipeline.electrochemistry.model import (
        conditions_from_environment, electrochemical_evidence)

    evidence = electrochemical_evidence(conditions_from_environment())
    phase = evidence['conditions'].get('electrolyte_phase')
    if config.multiphysics_artifact:
        solver_evidence = config.multiphysics_artifact
        artifact = solver_evidence['artifact']
        outputs = artifact['outputs']
        phase = artifact.get('electrolyte_phase', phase)
        return {
            'status': 'complete', 'valid': True,
            'reactor_type': 'Electrochemical',
            'pathway_mode': 'electrochemical', 'electrolyte_phase': phase,
            'catalyst_name': config.catalyst_name,
            'candidate_id': config.candidate_id, 'T_K': config.T_inlet_K,
            **{name: float(outputs[name]) for name in (
                'CH4_conversion', 'H2_selectivity',
                'faradaic_efficiency_H2', 'current_density_A_cm2',
                'cell_voltage_V', 'electrical_power_density_W_cm2')},
            'electrochemical_evidence': evidence,
            'multiphysics_evidence': solver_evidence,
            'reactor_evidence_tier': 'mechanistic_multiphysics_screening',
            'can_exclude_candidate': False,
        }
    return {
        'status': 'validation_required', 'valid': False,
        'reactor_type': 'Electrochemical', 'pathway_mode': 'electrochemical',
        'electrolyte_phase': phase, 'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K, 'electrochemical_evidence': evidence,
        'reactor_evidence_tier': 'pathway_model_pending',
        'can_exclude_candidate': False,
        'limitations': [
            'candidate_specific_electrochemical_kinetics_required',
            f'{phase or "unspecified"}_electrolyte_transport_model_required',
        ],
    }
