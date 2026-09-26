#!/usr/bin/env python3
# Turquoise-hydrogen reactor process models.
"""
Cantera Reactor-Scale Simulation Models for Methane Pyrolysis.

Three Cantera-backed reactor archetypes with distinct carbon-handling physics:
  A. MMBCR — melt ODE to tabulated X_eq with bubble-area flotation. Residence
     time is column height over the Mendelson bubble rise velocity; gas holdup
     is derived, not an input. Carbon leaves the bubble; no site lattice.
  B. PFR — one shared Langmuir surface marched through stages (θ_C is
     time-on-stream, not axial); optional discrete non-oxidative regen.
  C. Fluidized — emulsion pass with in-step C_s removal (circulating) or
     batch regen, then bubble bypass mixed on molar flows via the Ar tracer.

NTEC and electrochemical pathways consume validated external multiphysics
artifacts. Cantera supplies chemistry within those coupled calculations, but
does not itself solve their mechanical/electrical/charge-transport physics.

Solid carbon is never a gas-phase species. Surface C_s blocks sites on solid
paths until removed by a named policy or a gated off-site channel (B6).
Oxidative regen requires co2_permitted. Γ is a monolayer and is locked (B1).

Reactor routing is owned by ``pipeline.reactors.modes``: a reactor
type must belong to the selected pathway mode, and MMBCR is applicable only
to MoltenMetal candidates.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

try:
    import cantera as ct
    HAS_CANTERA = True
except ImportError:
    HAS_CANTERA = False

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.utils import (
    REACTOR_DIR,
    setup_logger, print_banner, save_json,
)
from pipeline.reactors.modes import (
    DEFAULT_MODE, REACTOR_MODELS, reactor_applicability,
    reactor_types_for_mode, validate_mode_reactors)
from pipeline.data_models.reactors import ReactorResult

from pipeline.reactors.reactor_core import *  # noqa: F403

logger = setup_logger('reactor_models', 'reactor/reactor_simulation.log')

def simulate_mmbcr(config: ReactorConfig) -> ReactorResult:
    """Run the molten-metal bubble-column implementation.

    Args:
        config: Validated reactor configuration.

    Returns:
        MMBCR screening result with evidence metadata.
    """
    from pipeline.reactors.mmbcr import simulate_mmbcr as implementation
    return implementation(
        config, cantera_available=HAS_CANTERA,
        mock_result=_mock_reactor_result)


def simulate_pfr(config: ReactorConfig) -> ReactorResult:
    """Run the thermocatalytic plug-flow implementation.

    Args:
        config: Validated reactor configuration.

    Returns:
        Thermocatalytic PFR screening result with evidence metadata.
    """
    from pipeline.reactors.pfr import simulate_pfr as implementation
    return implementation(
        config, cantera_available=HAS_CANTERA,
        mock_result=_mock_reactor_result)


def simulate_fluidized_bed(config: ReactorConfig) -> ReactorResult:
    """Run the thermocatalytic fluidized-bed implementation.

    Args:
        config: Validated reactor configuration.

    Returns:
        Fluidized-bed screening result with evidence metadata.
    """
    from pipeline.reactors.fluidized_bed import (
        simulate_fluidized_bed as implementation)
    return implementation(
        config, cantera_available=HAS_CANTERA,
        mock_result=_mock_reactor_result)


# ═══════════════════════════════════════════════════════════════════════════════
# MOCK RESULTS (for testing without Cantera)
# ═══════════════════════════════════════════════════════════════════════════════

def _mock_reactor_result(config: ReactorConfig,
                         reactor_type: str) -> ReactorResult:
    logger.warning(f"Cantera not available. Generating mock {reactor_type} results.")
    _validate_carbon_policy(config)
    E_act = config.catalyst_E_act_eV
    k_B_eV = 8.617e-5
    k = 1e13 * np.exp(-E_act / (k_B_eV * config.T_inlet_K))
    extra: Dict = {}
    if reactor_type == 'MMBCR':
        hydro = _mmbcr_hydrodynamics(config)
        tau = hydro['residence_time_s']
        extra = hydro
    elif reactor_type == 'PFR':
        tau = config.bed_length_m * config.bed_porosity / max(config.gas_velocity_m_s, 0.1)
    else:
        tau = config.bed_height_m / max(config.gas_velocity_m_s, 0.05)
    conversion = float(np.clip(1.0 - np.exp(-k * tau * 1e-12), 0.01, 0.99))
    return {
        'reactor_type': reactor_type,
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'catalyst_E_act_eV': E_act,
        'residence_time_s': tau,
        **extra,
        'WHSV_h-1': reciprocal_residence_h(tau),
        'CH4_conversion': conversion,
        'single_pass_CH4_conversion': conversion,
        'conversion_basis': (
            'melt_ode_to_Xeq' if reactor_type == 'MMBCR' else 'argon_tracer'),
        'thermal_mode': (
            'isothermal_by_construction' if reactor_type == 'MMBCR'
            else 'isothermal_energy_disabled'),
        'exit_T_K': config.T_inlet_K,
        **kinetics_fields(config),
        'H2_atom_balance': 0.95,
        'H2_selectivity': 0.95,
        'solid_C_selectivity': 0.90,
        'exit_x_H2': conversion * 0.95 * 2.0 / (1.0 + conversion * 0.95),
        'mock': True,
        **_kinetics_evidence(config),
        **_policy_metadata(config),
        **_load_status_fields(config, None, None),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED SIMULATION INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

from pipeline.reactors.electrochemical_models import (
    simulate_electrochemical_pathway, simulate_ntec_pathway)


def simulate_reactor(config: ReactorConfig,
                     coupling_services=None) -> ReactorResult:
    """Run the appropriate reactor simulation based on config.reactor_type.

    Args:
        config: Configuration controlling this operation.
        coupling_services: Coupling services used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    simulators = {
        'MMBCR': simulate_mmbcr,
        'PFR': simulate_pfr,
        'Fluidized': simulate_fluidized_bed,
        'NTEC': simulate_ntec_pathway,
        'Electrochemical': simulate_electrochemical_pathway,
    }
    if config.reactor_type not in simulators:
        raise ValueError(f'unsupported reactor type: {config.reactor_type}')
    selected_mode = reactor_types_for_mode(config.pathway_mode)
    if config.reactor_type not in selected_mode:
        raise ValueError(
            f'reactor {config.reactor_type!r} is not valid for pathway mode '
            f'{config.pathway_mode!r}; expected {selected_mode}')
    applicable, reason = reactor_applicability(
        config.reactor_type, config.material_class)
    spec = REACTOR_MODELS[config.reactor_type]
    if not applicable:
        result = {
            'status': 'not_applicable', 'valid': False,
            'reactor_type': config.reactor_type,
            'pathway_mode': config.pathway_mode,
            'catalyst_name': config.catalyst_name,
            'T_K': config.T_inlet_K,
            'material_class': config.material_class,
            'reason': reason, 'can_exclude_candidate': False,
            'reactor_evidence_tier': 'not_applicable_non_excluding',
            'bed_or_interface': spec.bed_or_interface,
            'cantera_reactor_model': spec.cantera_model,
            'reaction_domain': spec.reaction_domain,
            'reactor_model_fidelity': spec.fidelity,
        }
    else:
        from pipeline.simulation.result_contract import EXTERNAL_SOLVERS
        if config.reactor_type in EXTERNAL_SOLVERS:
            from pipeline.simulation.reactor_handoff import (
                default_reactor_coupling_services)
            coupling_services = (coupling_services or
                                 default_reactor_coupling_services())
            loaded = coupling_services.load_artifact(
                config.multiphysics_results_dir, config.candidate_id,
                config.pathway_mode, config.reactor_type, config.T_inlet_K)
            loaded = coupling_services.validate_compatibility(config, loaded)
            if not loaded['valid']:
                from pipeline.transport.closure_provider import (
                    resolve_reactor_closure)
                surrogate = (coupling_services.load_surrogate(
                    config.pathway_mode, config.reactor_type)
                    if coupling_services.load_surrogate else None)
                closure = resolve_reactor_closure(
                    full_physics=loaded, surrogate=surrogate,
                    features=config.closure_features,
                    pathway_mode=config.pathway_mode,
                    reactor_type=config.reactor_type,
                    temperature_K=config.T_inlet_K)
                if closure['available']:
                    from pipeline.simulation.reactor_handoff import (
                        couple_surrogate_closure)
                    couple_surrogate_closure(config, closure)
                else:
                    loaded = {**loaded, 'closure_resolution': closure}
            if (not loaded['valid'] and not config.reactor_closure_evidence
                    and config.reactor_type in ANALYTICAL_CLOSURE_REACTORS):
                # No validated artifact or calibrated surrogate. Fluidized and
                # MMBCR carry an explicit analytical bubbly-flow closure
                # (Mendelson u_b / derived eps_g; clipped (u0-umf)/u0) so the
                # screening run can proceed with its provenance labelled.
                # It can never exclude a candidate. NTEC / Electrochemical
                # have no such closure and stay validation_required.
                config.reactor_closure_evidence = {
                    'source': ANALYTICAL_CLOSURE_SOURCE,
                    'candidate_exclusion_authorized': False,
                    'basis': ANALYTICAL_CLOSURE_BASIS[config.reactor_type],
                    'multiphysics_evidence': loaded,
                }
            if not loaded['valid'] and not config.reactor_closure_evidence:
                result = {
                    'status': 'validation_required', 'valid': False,
                    'reactor_type': config.reactor_type,
                    'pathway_mode': config.pathway_mode,
                    'catalyst_name': config.catalyst_name,
                    'T_K': config.T_inlet_K,
                    'candidate_id': config.candidate_id,
                    'material_class': config.material_class,
                    'multiphysics_evidence': loaded,
                    'can_exclude_candidate': False,
                    'reactor_evidence_tier': 'required_solver_evidence_missing',
                    'bed_or_interface': spec.bed_or_interface,
                    'cantera_reactor_model': spec.cantera_model,
                    'reaction_domain': spec.reaction_domain,
                    'reactor_model_fidelity': spec.fidelity,
                }
                REACTOR_DIR.mkdir(parents=True, exist_ok=True)
                fname = (f"{config.reactor_type}_{config.catalyst_name}_"
                         f"{int(config.T_inlet_K)}K.json")
                save_json(result, fname, subdir='reactor')
                return result
            if loaded['valid']:
                coupling_services.couple_evidence(config, loaded)
        _validate_reactor_config(config)
        _validate_carbon_policy(config)
        result = simulators[config.reactor_type](config)
        result.setdefault('status', 'complete')
        result.update({
            'pathway_mode': config.pathway_mode,
            'material_class': config.material_class,
            'bed_or_interface': spec.bed_or_interface,
            'cantera_reactor_model': spec.cantera_model,
            'reaction_domain': spec.reaction_domain,
            'reactor_model_fidelity': spec.fidelity,
            'multiphysics_evidence': config.multiphysics_artifact,
            'reactor_closure_evidence': config.reactor_closure_evidence,
        })
        if spec.fidelity != 'validated_predictive':
            limitations = list(result.get('reactor_evidence_limitations', []))
            limitation = f'{config.reactor_type.lower()}_screening_model_not_validated'
            if limitation not in limitations:
                limitations.append(limitation)
            result['reactor_evidence_limitations'] = limitations
            result['can_exclude_candidate'] = False

    # Save result
    REACTOR_DIR.mkdir(parents=True, exist_ok=True)
    fname = f"{config.reactor_type}_{config.catalyst_name}_{int(config.T_inlet_K)}K.json"
    save_json(result, fname, subdir='reactor')

    return result


def run_reactor_sweep(catalyst_name: str, mechanism_file: str,
                      temperatures: List[float] = None,
                      reactor_types: List[str] = None,
                      catalyst_E_act_eV: float = 0.8,
                      pathway_mode: str = DEFAULT_MODE,
                      material_class: Optional[str] = None,
                      candidate_id: str = 'unknown',
                      multiphysics_results_dir: Optional[str] = None,
                      coupling_services=None,
                      catalyst_dE_H_eV: float = 0.0,
                      reactor_config_kwargs: Optional[Dict] = None
                      ) -> list[ReactorResult]:
    """Sweep operating conditions across the reactors owned by a pathway.

    Args:
        catalyst_name: Human-readable catalyst identifier used in output names.
        mechanism_file: Candidate-specific Cantera mechanism path.
        temperatures: Absolute operating temperatures in kelvin; defaults to the
            standard four-point screening sweep.
        reactor_types: Ordered reactor implementations; defaults to those owned
            by ``pathway_mode``.
        catalyst_E_act_eV: Activation barrier in eV. Used by the MMBCR melt
            k_if and by the explicitly labeled mock path.
        pathway_mode: Methane-conversion pathway controlling valid reactors.
        material_class: Catalyst class used for physical-bed compatibility.
        candidate_id: Stable identity required by multiphysics evidence.
        multiphysics_results_dir: Root containing validated external artifacts.
        coupling_services: Optional injected evidence-loading/coupling adapters.
        catalyst_dE_H_eV: Screening H adsorption energy, carried into results.
        reactor_config_kwargs: Extra ``ReactorConfig`` fields (carbon policy,
            solids inventory, melt properties) applied to every condition.

    Returns:
        One result record per requested reactor and temperature. Individual
        condition failures remain non-excluding records rather than aborting the
        complete sweep.
    """
    if temperatures is None:
        temperatures = [773.15, 900.0, 1100.0, 1300.0]
    if reactor_types is None:
        # The default family is conventional thermocatalysis. MMBCR, NTEC, and
        # electrochemical physics must be selected explicitly.
        reactor_types = list(reactor_types_for_mode(pathway_mode))
    validate_mode_reactors(pathway_mode, reactor_types)
    extra = dict(reactor_config_kwargs or {})

    results = []
    for rt in reactor_types:
        for T in temperatures:
            config = ReactorConfig(
                T_inlet_K=T,
                reactor_type=rt,
                mechanism_file=str(mechanism_file),
                catalyst_name=catalyst_name,
                catalyst_E_act_eV=catalyst_E_act_eV,
                catalyst_dE_H_eV=catalyst_dE_H_eV,
                pathway_mode=pathway_mode,
                material_class=material_class,
                candidate_id=candidate_id,
                multiphysics_results_dir=multiphysics_results_dir,
                **extra,
            )
            try:
                result = (simulate_reactor(config) if coupling_services is None
                          else simulate_reactor(
                              config, coupling_services=coupling_services))
            except Exception as exc:
                # A stiff condition must not discard the other temperatures or
                # reactor types. Preserve it as non-excluding failed evidence.
                result = {
                    'status': 'failed',
                    'valid': False,
                    'reactor_type': rt,
                    'pathway_mode': pathway_mode,
                    'material_class': material_class,
                    'temperature_K': float(T),
                    'catalyst': catalyst_name,
                    'error_type': type(exc).__name__,
                    'error': str(exc),
                    'can_exclude_candidate': False,
                    'reactor_evidence_tier': 'failed_simulation',
                }
                REACTOR_DIR.mkdir(parents=True, exist_ok=True)
                fname = f"{rt}_{catalyst_name}_{int(T)}K.json"
                save_json(result, fname, subdir='reactor')
                logger.error(
                    'Reactor condition failed without excluding candidate: '
                    f'{rt} {catalyst_name} {T} K: {exc}')
            results.append(result)

    return results


if __name__ == '__main__':
    print_banner("REACTOR SIMULATION TEST")

    # Test with mock data (no Cantera needed)
    config = ReactorConfig(
        T_inlet_K=1000.0,
        reactor_type='MMBCR',
        catalyst_name='NiBi_test',
        material_class='MoltenMetal',
        pathway_mode='mmbcr',
    )
    result = simulate_reactor(config)
    print(json.dumps(result, indent=2))
