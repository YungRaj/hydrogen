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
    electrochemical_descriptor_eV: float | None
    electrolyte_conductivity_S_m: float | None
    electrochemical_gap_m: float | None
    electrochemical_area_m2: float | None
    ntec_shear_rate_s: float | None
    ntec_particle_diameter_m: float | None
    ntec_number_density_m3: float | None
    ntec_effective_mass_kg: float | None
    ntec_adhesion_energy_J: float | None
    native_electrochemical_parameters: dict | None
    native_ntec_parameters: dict | None


def _pathway_result(
    config: ElectrochemicalReactorConfig,
    reactor_type: str,
    pathway_mode: str,
    *,
    complete: bool,
) -> ReactorResult:
    """Return the shared identity and disposition for a coupled pathway."""
    return {
        'status': 'complete' if complete else 'validation_required',
        'valid': complete,
        'reactor_type': reactor_type,
        'pathway_mode': pathway_mode,
        'catalyst_name': config.catalyst_name,
        'candidate_id': config.candidate_id,
        'T_K': config.T_inlet_K,
        'can_exclude_candidate': False,
    }


def simulate_ntec_pathway(config: ElectrochemicalReactorConfig) -> ReactorResult:
    """Build an NTEC result from validated coupled-physics evidence.

    Args:
        config: Reactor configuration containing identity and optional evidence.

    Returns:
        A complete result or a non-excluding validation requirement.
    """
    from pipeline.electrochemistry.ntec import (
        conditions_from_environment,
        ntec_assistance,
    )

    assistance = ntec_assistance(conditions_from_environment())
    if config.multiphysics_artifact:
        evidence = config.multiphysics_artifact
        artifact = evidence['artifact']
        outputs = artifact['outputs']
        return {
            **_pathway_result(config, 'NTEC', 'ntec', complete=True),
            'CH4_conversion': float(outputs['CH4_conversion']),
            'H2_selectivity': float(outputs['H2_selectivity']),
            'solid_C_selectivity': float(outputs['solid_C_selectivity']),
            'specific_energy_kWh_kg_H2': float(outputs['specific_energy_kWh_kg_H2']),
            'ntec_assistance': assistance,
            'multiphysics_evidence': evidence,
            'reactor_evidence_tier': 'calibrated_multiphysics_screening',
        }
    if config.native_ntec_parameters:
        from pipeline.electrochemistry.native_models import ntec_screening_point
        diagnostic = ntec_screening_point(config.native_ntec_parameters)
        diagnostic.update(status='validation_required', valid=False)
        return {
            **_pathway_result(config, 'NTEC', 'ntec', complete=False),
            **diagnostic,
            'limitations': ['native_ntec_diagnostic_not_kinetically_coupled'],
        }
    return {
        **_pathway_result(config, 'NTEC', 'ntec', complete=False),
        'ntec_assistance': assistance,
        'reactor_evidence_tier': 'pathway_model_pending',
        'limitations': [
            'candidate_specific_ntec_pathway_kinetics_required',
            'liquid_solid_hydrodynamic_model_required',
        ],
    }


def simulate_electrochemical_pathway(
    config: ElectrochemicalReactorConfig,
) -> ReactorResult:
    """Build an electrochemical result without inventing conversion evidence.

    Args:
        config: Reactor configuration containing identity and optional evidence.

    Returns:
        A complete result or a non-excluding validation requirement.
    """
    from pipeline.electrochemistry.model import (
        conditions_from_environment,
        electrochemical_evidence,
    )

    evidence = electrochemical_evidence(conditions_from_environment())
    phase = evidence['conditions'].get('electrolyte_phase')
    if config.multiphysics_artifact:
        solver_evidence = config.multiphysics_artifact
        artifact = solver_evidence['artifact']
        outputs = artifact['outputs']
        phase = artifact.get('electrolyte_phase', phase)
        return {
            **_pathway_result(
                config, 'Electrochemical', 'electrochemical', complete=True
            ),
            'electrolyte_phase': phase,
            **{
                name: float(outputs[name])
                for name in (
                    'CH4_conversion',
                    'H2_selectivity',
                    'faradaic_efficiency_H2',
                    'current_density_A_cm2',
                    'cell_voltage_V',
                    'electrical_power_density_W_cm2',
                )
            },
            'electrochemical_evidence': evidence,
            'multiphysics_evidence': solver_evidence,
            'reactor_evidence_tier': 'mechanistic_multiphysics_screening',
        }
    if config.native_electrochemical_parameters:
        from pipeline.electrochemistry.native_models import electrochemical_screening_point
        diagnostic = electrochemical_screening_point(
            config.native_electrochemical_parameters
        )
        diagnostic.update(status='validation_required', valid=False)
        return {
            **_pathway_result(config, 'Electrochemical', 'electrochemical', complete=False),
            **diagnostic,
            'limitations': ['native_electrochemical_diagnostic_not_feed_coupled'],
        }
    return {
        **_pathway_result(
            config, 'Electrochemical', 'electrochemical', complete=False
        ),
        'electrolyte_phase': phase,
        'electrochemical_evidence': evidence,
        'reactor_evidence_tier': 'pathway_model_pending',
        'limitations': [
            'candidate_specific_electrochemical_kinetics_required',
            f'{phase or "unspecified"}_electrolyte_transport_model_required',
        ],
    }
