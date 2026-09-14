"""Unit tests for independently usable multiphysics data/model components."""

import json
import math
import tempfile
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.process.multifidelity_surrogate import (
    PhysicsRecord, TransportSurrogate, fit_transport_surrogate,
    record_from_artifact)
from pipeline.process.physical_case import case_template, surrogate_inputs
from pipeline.stages.reactor import (
    ReactorStageServices, simulate_candidate, summarize_reactor_sweep)
from pipeline.stages.orchestration import PipelineComponents, PipelineRuntime
from pipeline.stages.contracts import StageOutcome


@contextmanager
def _raises(message, exception_type=ValueError):
    try:
        yield
    except exception_type as exc:
        assert message in str(exc), (message, str(exc))
    else:
        raise AssertionError(
            f'expected {exception_type.__name__} containing {message!r}')


def _record(case_id, x, *, mode='mmbcr', reactor='MMBCR'):
    return PhysicsRecord(
        case_id=case_id, pathway_mode=mode, reactor_type=reactor,
        features={'operating.temperature_K': float(x),
                  'geometry.column_diameter_m': float(1 + x / 1000)},
        outputs={'gas_velocity_m_s': 0.25 + 0.002 * x,
                 'gas_holdup_fraction': 0.1 + 0.0001 * x,
                 'bubble_diameter_mm': 2.0 + 0.00001 * x * x})


def _fit():
    training = [_record(f'train-{x}', x) for x in range(100, 701, 30)]
    validation = [_record(f'holdout-{x}', x) for x in (175, 325, 475, 625)]
    return fit_transport_surrogate(
        training, validation,
        targets=('gas_velocity_m_s', 'gas_holdup_fraction',
                 'bubble_diameter_mm'),
        validation_rmse_limits={
            'gas_velocity_m_s': 1e-6, 'gas_holdup_fraction': 1e-6,
            'bubble_diameter_mm': 1e-6},
        ensemble_size=12, random_seed=7)


def test_surrogate_recovers_known_relationship_and_round_trips():
    model = _fit()
    features = _record('query', 400).features
    result = model.predict(features)
    assert result['usable'] is True
    assert result['candidate_exclusion_authorized'] is False
    assert abs(result['predictions']['gas_velocity_m_s'] - 1.05) < 1e-6
    with tempfile.TemporaryDirectory() as tmp:
        target = model.save(Path(tmp) / 'transport.json')
        restored = TransportSurrogate.load(target)
        assert restored.predict(features)['predictions'] == result['predictions']


def test_surrogate_escalates_outside_domain_and_on_schema_mismatch():
    model = _fit()
    outside = model.predict(_record('query', 900).features)
    assert outside['decision'] == 'full_physics_required'
    assert outside['reason'] == 'outside_calibrated_domain'
    missing = model.predict({'operating.temperature_K': 400.0})
    assert missing['decision'] == 'full_physics_required'


def test_surrogate_escalates_a_nonphysical_in_domain_closure():
    model = _fit()
    coefficients = model.coefficients.copy()
    coefficients[:, :, 2] = 0.0
    coefficients[:, 0, 2] = -1.0
    unsafe = replace(model, coefficients=coefficients)
    result = unsafe.predict(_record('query', 400).features)
    assert result['usable'] is False
    assert result['decision'] == 'full_physics_required'
    assert result['reason'] == 'nonphysical_closure'


def test_training_rejects_leakage_cross_mode_and_failed_holdout():
    training = [_record(f'case-{x}', x) for x in range(100, 701, 30)]
    with _raises('disjoint'):
        fit_transport_surrogate(
            training, [_record('case-100', 100), _record('holdout', 200)],
            targets=('gas_velocity_m_s',),
            validation_rmse_limits={'gas_velocity_m_s': 1.0})
    with _raises('exactly one pathway'):
        fit_transport_surrogate(
            training, [_record('h1', 200, mode='ntec'),
                       _record('h2', 300, mode='ntec')],
            targets=('gas_velocity_m_s',),
            validation_rmse_limits={'gas_velocity_m_s': 1.0})
    bad = [PhysicsRecord('h1', 'mmbcr', 'MMBCR', _record('x', 200).features,
                         {'gas_velocity_m_s': 99.0}),
           PhysicsRecord('h2', 'mmbcr', 'MMBCR', _record('x', 300).features,
                         {'gas_velocity_m_s': 99.0})]
    with _raises('holdout validation failed'):
        fit_transport_surrogate(
            training, bad, targets=('gas_velocity_m_s',),
            validation_rmse_limits={'gas_velocity_m_s': 0.01})


def test_overall_ntec_performance_is_not_mislabeled_as_transport_closure():
    def row(case_id, x):
        return PhysicsRecord(
            case_id, 'ntec', 'NTEC', {'operating.temperature_K': float(x)},
            {'CH4_conversion': 0.4, 'H2_selectivity': 0.8,
             'solid_C_selectivity': 0.8,
             'specific_energy_kWh_kg_H2': 10.0})
    with _raises('not established transport closures'):
        fit_transport_surrogate(
            [row(f't{x}', x) for x in range(4)],
            [row('v1', 5), row('v2', 6)], targets=('CH4_conversion',),
            validation_rmse_limits={'CH4_conversion': 0.1})


def test_artifact_record_requires_numeric_snapshot_and_convergence():
    artifact = {
        'complete': True, 'convergence': {'converged': True},
        'candidate_id': 'candidate', 'pathway_mode': 'mmbcr',
        'reactor_type': 'MMBCR',
        'provenance': {'input_sha256': 'a' * 64},
        'surrogate_inputs': {
            'schema_version': 1, 'reactor_type': 'MMBCR',
            'units_in_field_names': True, 'values': {'temperature_K': 900.0}},
        'outputs': {'gas_velocity_m_s': 1.0, 'gas_holdup_fraction': 0.2,
                    'bubble_diameter_mm': 3.0},
    }
    row = record_from_artifact(artifact)
    assert row.case_id == 'a' * 64
    del artifact['surrogate_inputs']
    with _raises('snapshot'):
        record_from_artifact(artifact)


def test_physical_case_snapshot_is_unit_bearing_and_normalizes_feed():
    case = case_template(candidate_id='c', mode='mmbcr', reactor_type='MMBCR',
                         temperature_K=900.0)
    case['feed']['composition'] = {'CH4': 3.0, 'Ar': 1.0}
    snapshot = surrogate_inputs(case)
    assert snapshot['units_in_field_names'] is True
    assert snapshot['values']['feed.mole_fraction.CH4'] == 0.75
    assert snapshot['values']['feed.mole_fraction.Ar'] == 0.25
    assert all(not isinstance(value, str) for value in snapshot['values'].values())


def test_unsafe_serialized_model_is_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        value = _fit().to_dict()
        value['candidate_exclusion_authorized'] = True
        target = Path(tmp) / 'unsafe.json'
        target.write_text(json.dumps(value))
        with _raises('unsafe'):
            TransportSurrogate.load(target)


def test_malformed_serialized_model_shape_and_split_are_rejected():
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / 'malformed.json'
        value = _fit().to_dict()
        value['center'] = [1.0]
        target.write_text(json.dumps(value))
        with _raises('invalid shape'):
            TransportSurrogate.load(target)
        value = _fit().to_dict()
        value['validation_case_ids'][0] = value['training_case_ids'][0]
        target.write_text(json.dumps(value))
        with _raises('split-disjoint'):
            TransportSurrogate.load(target)


def test_reactor_stage_runs_with_in_memory_components():
    calls = []

    class Mode:
        requires_specialized_validation = False

    def build(row, **kwargs):
        calls.append(('kinetics', kwargs['candidate_id']))
        return {'barrier': row['E_act']}

    def write(name, **kwargs):
        calls.append(('mechanism', name, kwargs['kinetics']))
        return Path('/portable/generated.yaml')

    def sweep(name, mechanism, **kwargs):
        calls.append(('sweep', name, mechanism, kwargs))
        return [
            {'status': 'complete', 'CH4_conversion': 0.4,
             'can_exclude_candidate': False},
            {'status': 'complete', 'CH4_conversion': 0.6,
             'can_exclude_candidate': False},
        ]

    services = ReactorStageServices(
        resolve_mode=lambda name: Mode(), build_kinetics=build,
        write_mechanism=write, run_sweep=sweep)
    result = simulate_candidate(
        {'candidate_id': 'stable-id', 'E_act': 0.7,
         'material_class': 'SolidCatalyst'},
        'candidate', [800.0, 900.0], services=services)
    assert [call[0] for call in calls] == ['kinetics', 'mechanism', 'sweep']
    assert result['sweep_status'] == 'complete'
    assert result['best_condition']['CH4_conversion'] == 0.6
    assert result['can_exclude_candidate'] is False


def test_reactor_sweep_summary_preserves_non_excluding_partial_failure():
    summary = summarize_reactor_sweep([
        {'status': 'complete', 'CH4_conversion': 0.8,
         'can_exclude_candidate': True},
        {'status': 'failed', 'can_exclude_candidate': False},
    ])
    assert summary['sweep_status'] == 'partial'
    assert summary['completed_conditions'] == 1
    assert summary['failed_conditions'] == 1
    assert summary['can_exclude_candidate'] is False


def test_orchestrator_runtime_is_injectable_and_config_is_not_mutated():
    from pipeline.orchestrator import (
        PipelineConfig, normalized_pipeline_config, run_pipeline)

    original = PipelineConfig(quick_mode=True, initial_fairchem_samples=999)
    effective = normalized_pipeline_config(original)
    assert original.initial_fairchem_samples == 999
    assert effective.initial_fairchem_samples == 50

    saved, modes, banners = [], [], []
    ticks = iter((10.0, 14.0))
    runtime = PipelineRuntime(
        clock=lambda: next(ticks), banner=banners.append,
        load_state=lambda: {'existing': True},
        save_state=lambda state: saved.append(dict(state)),
        select_pathway_mode=modes.append)
    result = run_pipeline(original, start_phase=7, end_phase=6, runtime=runtime)
    assert result['existing'] is True and result['total_elapsed_s'] == 4.0
    assert modes == ['thermocatalytic']
    assert len(saved) == 1
    assert original.initial_fairchem_samples == 999


def test_orchestrator_stage_graph_is_replaceable():
    from pipeline.orchestrator import PipelineConfig, run_pipeline

    calls, saved = [], []

    def forbidden(*args, **kwargs):
        raise AssertionError('unselected component executed')

    def vqe(**kwargs):
        calls.append(kwargs)
        return StageOutcome(
            state={'n_vqe_runs': 1}, products={'vqe_results': [{'ok': True}]})

    components = PipelineComponents(
        discovery=forbidden, reactor_batch=forbidden, dft=forbidden,
        vqe=vqe, fuel_cell=forbidden, report=forbidden,
        load_candidates=forbidden)
    ticks = iter((0.0, 1.0, 3.0, 3.0, 4.0))
    runtime = PipelineRuntime(
        clock=lambda: next(ticks), banner=lambda value: None,
        load_state=dict, save_state=lambda state: saved.append(dict(state)),
        select_pathway_mode=lambda mode: None)
    state = run_pipeline(
        PipelineConfig(top_k_vqe=1, run_vqe=False), start_phase=4,
        end_phase=4, runtime=runtime, components=components)
    assert calls == [{'top_k': 1, 'execute_quantum': False}]
    assert state['phase4'] == {'n_vqe_runs': 1, 'elapsed_s': 2.0}
    assert len(saved) == 2


def test_all_six_stages_compose_in_memory_with_explicit_handoffs():
    import pandas as pd
    from pipeline.orchestrator import PipelineConfig, run_pipeline

    order, saved = [], []
    screening = pd.DataFrame([{
        'candidate_id': 'c', 'E_act': 0.5, 'material_class': 'SAC'}])

    def discovery(**kwargs):
        order.append('discovery')
        return StageOutcome(
            {'valid_count': 1}, {
                'design_space_sizes': {'TOTAL': 1}, 'pareto_genomes': ['g'],
                'screening_database': screening,
                'top_catalysts': screening, 'dft_candidates': screening})

    def reactor(candidates, **kwargs):
        order.append(('reactor', candidates is screening))
        return StageOutcome(
            {'n_simulations': 1}, {'reactor_results': [{'complete': True}]})

    def dft(candidates, **kwargs):
        order.append(('dft', candidates is screening))
        return StageOutcome(
            {'n_validated': 1}, {'dft_results': [{'converged': True}]})

    def vqe(**kwargs):
        order.append('vqe')
        return StageOutcome(
            {'n_vqe_runs': 1}, {'vqe_results': [{'converged': True}]})

    def fuel(**kwargs):
        order.append('fuel')
        return StageOutcome(
            {'n_cathodes_screened': 1}, {
                'cathode_database': screening, 'valid_cathodes': screening,
                'pemfc_results': [{}], 'stack_result': {}})

    def report(state):
        order.append(('report', tuple(sorted(state))))
        return StageOutcome(
            {'report_path': 'memory.md'}, {'report_path': Path('memory.md')})

    components = PipelineComponents(
        discovery=discovery, reactor_batch=reactor, dft=dft, vqe=vqe,
        fuel_cell=fuel, report=report,
        load_candidates=lambda *args, **kwargs: None)
    counter = iter(float(value) for value in range(100))
    runtime = PipelineRuntime(
        clock=lambda: next(counter), banner=lambda value: None,
        load_state=dict, save_state=lambda state: saved.append(dict(state)),
        select_pathway_mode=lambda mode: None)
    state = run_pipeline(
        PipelineConfig(), runtime=runtime, components=components)
    assert order[:5] == [
        'discovery', ('reactor', True), ('dft', True), 'vqe', 'fuel']
    assert order[5][0] == 'report'
    assert order[5][1] == ('phase1', 'phase2', 'phase3', 'phase4', 'phase5')
    assert all(f'phase{number}' in state for number in range(1, 7))
    assert len(saved) == 7


def test_vqe_stage_is_independent_and_target_is_explicit():
    from pipeline.stages.vqe import run_vqe_stage

    calls = []

    def validator(candidate, reaction, *, target):
        calls.append((candidate, reaction, target))
        return {'candidate': candidate, 'converged': True}

    outcome = run_vqe_stage(
        top_k=5, execute_quantum=True, validator=validator)
    assert outcome.state == {'n_vqe_runs': 3}
    assert len(outcome.products['vqe_results']) == 3
    assert calls == [
        ('champion_0', 'CH_split', 'nvidia'),
        ('champion_1', 'CH_split', 'nvidia'),
        ('champion_2', 'CH_split', 'nvidia')]
    assert run_vqe_stage(
        top_k=0, execute_quantum=False,
        validator=validator).products['vqe_results'] == []


def test_report_stage_is_independent_and_does_not_mutate_state():
    from pipeline.stages.report import run_report_stage

    state = {'phase1': {'valid_count': 4}}
    observed = []
    outcome = run_report_stage(
        state, generator=lambda value: observed.append(value) or Path('report.md'))
    assert outcome.state == {'report_path': 'report.md'}
    assert observed == [state] and observed[0] is not state
    assert 'phase6' not in state


def test_dft_stage_parses_genomes_and_isolates_candidate_failures():
    from pipeline.stages.dft import run_dft_stage

    calls = []

    def validator(name, genome, *, run_dft):
        calls.append((name, genome, run_dft))
        if genome[0] == 'bad':
            raise RuntimeError('injected failure')
        return {'candidate': name, 'converged': genome[0] == 'good'}

    outcome = run_dft_stage(
        [{'genome': "('good', 1)"}, {'genome': ('bad', 2)},
         {'genome': ('other', 3)}],
        top_k=3, execute_dft=False, validator=validator)
    assert outcome.state == {
        'n_validated': 2, 'n_converged': 1, 'n_failed': 1}
    assert len(outcome.products['failures']) == 1
    assert calls[0] == ('dft_cat_0', ('good', 1), False)
    assert calls[-1] == ('dft_cat_2', ('other', 3), False)


def test_discovery_stage_components_preserve_selection_products():
    import pandas as pd
    from pipeline.stages.discovery import DiscoveryServices, run_discovery_stage

    database = pd.DataFrame([
        {'valid': True, 'E_act': 0.5, 'coking_index': 0.8,
         'candidate_disposition': 'admit'},
        {'valid': False, 'E_act': 0.9, 'coking_index': 0.2,
         'candidate_disposition': 'rescue'},
    ])
    built = []
    services = DiscoveryServices(
        estimate_space=lambda: {'A': 10, 'TOTAL': 10},
        build_config=lambda **kwargs: built.append(kwargs) or kwargs,
        run_search=lambda config: (['pareto'], database),
        annotate_evidence=lambda frame, target: frame,
        select_reactor=lambda frame, count, target, **kwargs: frame.head(1),
        select_validation=lambda frame, count, target, **kwargs: frame.tail(1))
    outcome = run_discovery_stage(
        initial_samples=2, leaf_size=5, max_leaves=1,
        top_k_reactor=1, top_k_dft=1, services=services)
    assert built[0]['expected_space_size'] == 10
    assert outcome.state['valid_count'] == 1
    assert outcome.state['best_E_act'] == 0.5
    assert len(outcome.products['top_catalysts']) == 1
    assert len(outcome.products['dft_candidates']) == 1


def test_fuel_cell_stage_components_preserve_ranking_and_stack_inputs():
    import pandas as pd
    from pipeline.stages.fuel_cell import FuelCellServices, run_fuel_cell_stage

    cathodes = pd.DataFrame([
        {'valid': True, 'name': 'slow', 'orr_overpotential_V': 0.5,
         'material_class': 'SAC'},
        {'valid': True, 'name': 'fast', 'orr_overpotential_V': 0.2,
         'material_class': 'SAA'},
        {'valid': False, 'name': 'invalid', 'orr_overpotential_V': 0.1},
    ])
    swept, stack_configs = [], []

    def sweep(name, eta, **kwargs):
        swept.append(name)
        return [{'orr_overpotential_V': eta, 'efficiency_at_peak': 0.4,
                 'peak_power_W_cm2': 1.0, 'peak_voltage_V': 0.7,
                 'peak_current_A_cm2': 1.4}]

    def stack_config(**kwargs):
        stack_configs.append(kwargs)
        return kwargs

    services = FuelCellServices(
        screen_cathodes=lambda: cathodes, sweep_membranes=sweep,
        build_stack_config=stack_config,
        model_stack=lambda config: {'received': config})
    outcome = run_fuel_cell_stage(
        top_k_pemfc=1, stack_cells=300, services=services)
    assert swept == ['fast']
    assert stack_configs == [{
        'n_cells': 300, 'cell_voltage_V': 0.7,
        'current_density_A_cm2': 1.4}]
    assert outcome.state['n_valid'] == 2
    assert outcome.state['min_overpotential_V'] == 0.2


def test_reactor_batch_stage_uses_candidate_component_and_preserves_results():
    import pandas as pd
    from pipeline.stages.reactor_batch import (
        ReactorBatchServices, run_reactor_batch_stage)

    prepared, calls = [], []
    candidates = pd.DataFrame([{'candidate_id': 'one', 'E_act': 0.6}])

    def simulate(row, name, temperatures, reactors, **kwargs):
        calls.append((name, tuple(temperatures), tuple(reactors), kwargs))
        return {'sweep': [
            {'status': 'complete', 'CH4_conversion': 0.55},
            {'status': 'validation_required', 'CH4_conversion': 0.0}]}

    services = ReactorBatchServices(
        prepare_gas_mechanism=lambda: prepared.append(True),
        simulate_candidate=simulate,
        write_mock_mechanism=lambda *args, **kwargs: None,
        run_mock_sweep=lambda *args, **kwargs: [])
    outcome = run_reactor_batch_stage(
        candidates, temperatures=(800.0, 900.0), reactor_types=('PFR',),
        pathway_mode='thermocatalytic', multiphysics_results_dir='/portable',
        allow_mock_inputs=False, services=services)
    assert prepared == [True]
    assert calls[0][0] == 'cat_0'
    assert calls[0][3]['forbid_mock'] is True
    assert outcome.state == {'n_simulations': 2, 'best_conversion': 0.55}
    assert len(outcome.products['reactor_results']) == 2


def test_reactor_batch_mock_path_remains_explicitly_gated():
    from pipeline.stages.reactor_batch import (
        ReactorBatchServices, run_reactor_batch_stage)

    services = ReactorBatchServices(
        prepare_gas_mechanism=lambda: None,
        simulate_candidate=lambda *args, **kwargs: None,
        write_mock_mechanism=lambda *args, **kwargs: Path('mock.yaml'),
        run_mock_sweep=lambda *args, **kwargs: [])
    with _raises('mock catalyst fallback is disabled', RuntimeError):
        run_reactor_batch_stage(
            None, temperatures=(800.0,), reactor_types=('PFR',),
            pathway_mode='thermocatalytic',
            multiphysics_results_dir='/portable', allow_mock_inputs=False,
            services=services)


def test_artifact_store_is_independent_and_removes_rejected_evidence():
    from pipeline.process.artifact_store import persist_validated_artifact

    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / 'artifact.json'
        path_builder = lambda *args: target
        accepted = persist_validated_artifact(
            {'complete': True}, results_dir=tmp, candidate_id='c',
            mode='mmbcr', reactor_type='MMBCR', temperature_K=900,
            path_builder=path_builder,
            validator=lambda *args: {'valid': True})
        assert accepted == target
        assert json.loads(target.read_text()) == {'complete': True}
        with _raises('failed-check', RuntimeError):
            persist_validated_artifact(
                {'complete': False}, results_dir=tmp, candidate_id='c',
                mode='mmbcr', reactor_type='MMBCR', temperature_K=900,
                path_builder=path_builder,
                validator=lambda *args: {
                    'valid': False, 'reason': 'invalid',
                    'failed_checks': ['failed-check']})
        assert not target.exists()


def test_multiphysics_solver_discovery_is_injectable_before_case_access():
    from pipeline.process.multiphysics_runner import run_backend
    from pipeline.process.solver_execution import SolverExecutionServices

    execution = SolverExecutionServices(
        preflight=lambda mode: {
            'missing': ['openfoam'], 'solvers': {'openfoam': {}}},
        execute=lambda *args: (_ for _ in ()).throw(
            AssertionError('solver must not execute')),
        openfoam_version=lambda value: 'unused',
        fenics_command=lambda value: ([], 'unused'),
        cantera_version=lambda: 'unused')
    with _raises('required solver(s) unavailable', RuntimeError):
        run_backend(
            mode='mmbcr', reactor_type='MMBCR', candidate_id='c',
            temperature_K=900, case_dir='/does/not/exist', results_dir='/tmp',
            model_source='test', execution=execution)


def test_reactor_coupling_adapter_maps_only_validated_identity_matched_fields():
    from types import SimpleNamespace
    from pipeline.process.reactor_coupling import couple_multiphysics_evidence

    config = SimpleNamespace(
        candidate_id='candidate', pathway_mode='mmbcr', reactor_type='MMBCR',
        T_inlet_K=900.0, multiphysics_artifact=None,
        gas_velocity_m_s=0.0, gas_holdup_fraction=0.0,
        bubble_diameter_mm=0.0)
    loaded = {'valid': True, 'artifact': {
        'candidate_id': 'candidate', 'pathway_mode': 'mmbcr',
        'reactor_type': 'MMBCR', 'temperature_K': 900.0,
        'outputs': {'gas_velocity_m_s': 0.8,
                    'gas_holdup_fraction': 0.25,
                    'bubble_diameter_mm': 4.0}}}
    assert couple_multiphysics_evidence(config, loaded) is config
    assert config.gas_velocity_m_s == 0.8
    assert config.gas_holdup_fraction == 0.25
    assert config.bubble_diameter_mm == 4.0
    assert config.multiphysics_artifact is loaded
    mismatched = {'valid': True, 'artifact': {
        **loaded['artifact'], 'candidate_id': 'wrong'}}
    with _raises('identity mismatch'):
        couple_multiphysics_evidence(config, mismatched)


def test_simulate_reactor_accepts_in_injected_coupling_services():
    from pipeline.process.reactor_coupling import ReactorCouplingServices
    from pipeline.process.reactor_models import ReactorConfig, simulate_reactor

    calls = []

    def load(*args):
        calls.append(('load', args))
        return {'valid': False, 'reason': 'injected_missing'}

    services = ReactorCouplingServices(
        load_artifact=load,
        validate_compatibility=lambda config, value: value,
        couple_evidence=lambda *args: (_ for _ in ()).throw(
            AssertionError('invalid evidence must not couple')))
    result = simulate_reactor(ReactorConfig(
        reactor_type='MMBCR', pathway_mode='mmbcr',
        material_class='MoltenMetal', candidate_id='candidate',
        catalyst_name='candidate', T_inlet_K=900.0,
        multiphysics_results_dir='/portable'), coupling_services=services)
    assert calls and calls[0][0] == 'load'
    assert result['status'] == 'validation_required'
    assert result['can_exclude_candidate'] is False
    assert result['multiphysics_evidence']['reason'] == 'injected_missing'


def test_closure_provider_prefers_full_physics_and_rejects_temperature_mismatch():
    from pipeline.process.closure_provider import resolve_reactor_closure

    class MustNotRun:
        pathway_mode = 'mmbcr'
        reactor_type = 'MMBCR'
        def predict(self, features):
            raise AssertionError('full physics must take priority')

    full = {'valid': True, 'artifact': {'outputs': {'gas_velocity_m_s': 0.9}}}
    decision = resolve_reactor_closure(
        full_physics=full, surrogate=MustNotRun(), features={},
        pathway_mode='mmbcr', reactor_type='MMBCR', temperature_K=400.0)
    assert decision['source'] == 'validated_full_physics'
    assert decision['outputs']['gas_velocity_m_s'] == 0.9

    rejected = resolve_reactor_closure(
        full_physics={'valid': False, 'reason': 'missing'}, surrogate=_fit(),
        features=_record('query', 400).features,
        pathway_mode='mmbcr', reactor_type='MMBCR', temperature_K=401.0)
    assert rejected['available'] is False
    assert rejected['reason'] == 'transport_closure_temperature_mismatch'


def test_calibrated_surrogate_closure_flows_into_reduced_reactor_without_relabeling():
    from pipeline.process.reactor_coupling import ReactorCouplingServices
    from pipeline.process.reactor_models import ReactorConfig, simulate_reactor

    model = _fit()
    services = ReactorCouplingServices(
        load_artifact=lambda *args: {'valid': False, 'reason': 'missing'},
        validate_compatibility=lambda config, value: value,
        couple_evidence=lambda *args: (_ for _ in ()).throw(
            AssertionError('invalid full physics must not couple')),
        load_surrogate=lambda mode, reactor: model)
    config = ReactorConfig(
        reactor_type='MMBCR', pathway_mode='mmbcr',
        material_class='MoltenMetal', candidate_id='candidate',
        catalyst_name='candidate', T_inlet_K=400.0,
        closure_features=dict(_record('query', 400).features),
        multiphysics_results_dir='/portable')
    with patch('pipeline.process.reactor_models.HAS_CANTERA', False):
        result = simulate_reactor(config, coupling_services=services)
    evidence = result['reactor_closure_evidence']
    assert result['status'] == 'complete'
    assert result['multiphysics_evidence'] is None
    assert evidence['source'] == 'calibrated_transport_surrogate'
    assert evidence['candidate_exclusion_authorized'] is False
    assert len(evidence['model_sha256']) == 64
    assert evidence['training_case_ids'] == list(model.training_case_ids)


def test_transport_model_registry_round_trip_and_tamper_rejection():
    from pipeline.process.transport_model_registry import TransportModelRegistry

    with tempfile.TemporaryDirectory() as tmp:
        registry = TransportModelRegistry(tmp)
        model = _fit()
        assert registry.load('mmbcr', 'MMBCR') is None
        manifest = registry.publish(model)
        restored = registry.load('mmbcr', 'MMBCR')
        assert restored is not None and restored.sha256() == model.sha256()
        assert manifest.is_file()
        model_path = Path(tmp) / 'mmbcr__MMBCR.model.json'
        model_path.write_text(model_path.read_text() + ' ')
        with _raises('checksum mismatch'):
            registry.load('mmbcr', 'MMBCR')
        with _raises('safe path components'):
            registry.load('../escape', 'MMBCR')


def test_representative_case_design_is_deterministic_stratified_and_anchored():
    from pipeline.process.representative_cases import (
        ParameterRange, design_representative_cases)

    ranges = {
        'operating.temperature_K': ParameterRange(700.0, 1300.0),
        'properties.viscosity_Pa_s': ParameterRange(1e-5, 1e-2, 'log'),
    }
    anchor = {
        'operating.temperature_K': 900.0,
        'properties.viscosity_Pa_s': 1e-3,
    }
    first = design_representative_cases(
        pathway_mode='mmbcr', reactor_type='MMBCR', ranges=ranges,
        sample_count=8, anchors=[anchor], random_seed=11)
    second = design_representative_cases(
        pathway_mode='mmbcr', reactor_type='MMBCR', ranges=ranges,
        sample_count=8, anchors=[anchor], random_seed=11)
    assert first == second and len(first) == 9
    assert first[-1]['features'] == anchor
    assert len({row['case_id'] for row in first}) == len(first)
    sampled = first[:8]
    for name, limits in ranges.items():
        if limits.scale == 'log':
            values = [math.log(row['features'][name]) for row in sampled]
            low, high = math.log(limits.minimum), math.log(limits.maximum)
        else:
            values = [row['features'][name] for row in sampled]
            low, high = limits.minimum, limits.maximum
        bins = {
            min(7, int((value - low) / (high - low) * 8))
            for value in values}
        assert bins == set(range(8))


def test_representative_case_design_rejects_bad_anchors():
    from pipeline.process.representative_cases import (
        ParameterRange, design_representative_cases)

    with _raises('outside range'):
        design_representative_cases(
            pathway_mode='mmbcr', reactor_type='MMBCR',
            ranges={'operating.temperature_K': ParameterRange(700, 1300)},
            sample_count=4,
            anchors=[{'operating.temperature_K': 1400}])


def test_training_workflow_validates_partitions_before_atomic_publish():
    from pipeline.process.transport_model_registry import TransportModelRegistry
    from pipeline.process.transport_training import (
        CaseArtifactReference, train_and_publish_transport_model)

    records = [_record(f'case-{index}', 300 + index * 20)
               for index in range(6)]
    artifacts = {}
    references = []
    for index, record in enumerate(records):
        reference = CaseArtifactReference(
            candidate_id=f'candidate-{index}', pathway_mode='mmbcr',
            reactor_type='MMBCR', temperature_K=record.features[
                'operating.temperature_K'],
            partition='training' if index < 4 else 'validation')
        references.append(reference)
        artifacts[reference.candidate_id] = {
            'complete': True, 'pathway_mode': 'mmbcr', 'reactor_type': 'MMBCR',
            'convergence': {'converged': True},
            'surrogate_inputs': {
                'schema_version': 1, 'reactor_type': 'MMBCR',
                'units_in_field_names': True, 'values': dict(record.features)},
            'outputs': dict(record.outputs),
            'provenance': {'input_sha256': record.case_id},
        }

    def loader(root, candidate, mode, reactor, temperature):
        return {'valid': True, 'artifact': artifacts[candidate]}

    with tempfile.TemporaryDirectory() as tmp:
        registry = TransportModelRegistry(tmp)
        report = train_and_publish_transport_model(
            references, results_dir='/portable', model_registry=registry,
            targets=['gas_velocity_m_s'],
            validation_rmse_limits={'gas_velocity_m_s': 1.0},
            artifact_loader=loader, random_seed=2)
        assert report['status'] == 'published'
        assert report['candidate_exclusion_authorized'] is False
        assert registry.load('mmbcr', 'MMBCR') is not None


def test_training_workflow_never_publishes_a_partial_invalid_dataset():
    from pipeline.process.transport_model_registry import TransportModelRegistry
    from pipeline.process.transport_training import (
        CaseArtifactReference, train_and_publish_transport_model)

    references = [CaseArtifactReference(
        candidate_id=f'candidate-{index}', pathway_mode='mmbcr',
        reactor_type='MMBCR', temperature_K=400 + index,
        partition='training' if index < 4 else 'validation')
        for index in range(6)]

    def loader(root, candidate, mode, reactor, temperature):
        return {'valid': False, 'reason': 'not_converged'}

    with tempfile.TemporaryDirectory() as tmp:
        registry = TransportModelRegistry(tmp)
        with _raises('invalid full-physics artifact'):
            train_and_publish_transport_model(
                references, results_dir='/portable', model_registry=registry,
                targets=['gas_velocity_m_s'],
                validation_rmse_limits={'gas_velocity_m_s': 1.0},
                artifact_loader=loader)
        assert not list(Path(tmp).iterdir())


def test_full_physics_scheduler_preserves_coverage_then_uses_feedback():
    from pipeline.process.full_physics_scheduler import (
        FullPhysicsRequest, schedule_full_physics_cases)

    requests = [
        FullPhysicsRequest('a-low', 'class-a', unproductive_history=1.0),
        FullPhysicsRequest('a-high', 'class-a', expected_improvement=1.0),
        FullPhysicsRequest('b-only', 'class-b', unproductive_history=1.0),
        FullPhysicsRequest('c-disagree', 'class-c', disagreement=1.0),
        FullPhysicsRequest('c-low', 'class-c'),
    ]
    selected = schedule_full_physics_cases(
        requests, total_budget=4, minimum_per_region=1)
    assert {row['region'] for row in selected} == {
        'class-a', 'class-b', 'class-c'}
    assert any(row['case_id'] == 'b-only' for row in selected)
    assert any(row['case_id'] == 'c-disagree' for row in selected)
    assert all(row['candidate_exclusion_authorized'] is False
               for row in selected)
    assert sum(row['selection_reason'] == 'fixed_regional_coverage'
               for row in selected) == 3


def test_full_physics_scheduler_fails_when_coverage_budget_is_impossible():
    from pipeline.process.full_physics_scheduler import (
        FullPhysicsRequest, schedule_full_physics_cases)

    with _raises('fixed regional coverage'):
        schedule_full_physics_cases([
            FullPhysicsRequest('a', 'class-a'),
            FullPhysicsRequest('b', 'class-b')], total_budget=1)


def test_case_partitions_are_deterministic_preassigned_and_nonmutating():
    from pipeline.process.representative_cases import assign_case_partitions

    source = [{'case_id': f'case-{index}', 'features': {'x': index}}
              for index in range(8)]
    first = assign_case_partitions(source, validation_count=2, partition_seed='v1')
    second = assign_case_partitions(source, validation_count=2, partition_seed='v1')
    assert first == second
    assert sum(row['partition'] == 'validation' for row in first) == 2
    assert all('partition' not in row for row in source)


def test_campaign_ledger_detects_lineage_tampering():
    from pipeline.process.campaign_ledger import CampaignLedger

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'campaign.json'
        ledger = CampaignLedger(path)
        first = ledger.append('designed', {'case_ids': ['a', 'b']})
        second = ledger.append('trained', {'model_sha256': 'a' * 64})
        assert second['previous_event_sha256'] == first['event_sha256']
        assert len(ledger.read_verified()) == 2
        value = json.loads(path.read_text())
        value['events'][0]['payload']['case_ids'].append('injected')
        path.write_text(json.dumps(value))
        with _raises('hash chain'):
            ledger.read_verified()


def test_complete_multifidelity_iteration_executes_trains_screens_and_refers():
    from pipeline.process.campaign_ledger import CampaignLedger
    from pipeline.process.multifidelity_campaign import (
        MultiFidelityCampaignServices, run_multifidelity_iteration)
    from pipeline.process.transport_training import CaseArtifactReference

    model = _fit()
    plans = [{'case_id': f'case-{index}',
              'partition': 'training' if index < 4 else 'validation'}
             for index in range(6)]

    def execute(plan):
        return CaseArtifactReference(
            candidate_id=plan['case_id'], pathway_mode='mmbcr',
            reactor_type='MMBCR', temperature_K=400.0,
            partition=plan['partition'])

    def train(references):
        return {
            'status': 'published', 'pathway_mode': 'mmbcr',
            'reactor_type': 'MMBCR', 'model_sha256': model.sha256(),
            'training_case_ids': [row.candidate_id for row in references
                                  if row.partition == 'training'],
            'validation_case_ids': [row.candidate_id for row in references
                                    if row.partition == 'validation'],
        }

    queries = [
        {'query_id': 'inside', 'region': 'class-a',
         'features': dict(_record('q', 400).features)},
        {'query_id': 'outside-a', 'region': 'class-a',
         'features': dict(_record('q', 900).features), 'uncertainty': 0.7},
        {'query_id': 'outside-b', 'region': 'class-b',
         'features': dict(_record('q', 950).features), 'disagreement': 0.8},
    ]
    with tempfile.TemporaryDirectory() as tmp:
        ledger = CampaignLedger(Path(tmp) / 'ledger.json')
        services = MultiFidelityCampaignServices(
            execute_case=execute, train_and_publish=train,
            load_published_model=lambda mode, reactor: model,
            append_lineage=ledger.append)
        result = run_multifidelity_iteration(
            designed_cases=plans, screening_queries=queries, services=services,
            referral_budget=2, minimum_per_region=1)
        assert result['accepted_surrogate_count'] == 1
        assert result['full_physics_required_count'] == 2
        assert {row['region'] for row in result['scheduled_referrals']} == {
            'class-a', 'class-b'}
        assert result['candidate_exclusion_authorized'] is False
        assert len(result['lineage_event_sha256']) == 64
        events = ledger.read_verified()
        assert [event['event_type'] for event in events] == [
            'multifidelity_iteration_started',
            'multifidelity_iteration_complete']
        assert events[1]['previous_event_sha256'] == events[0]['event_sha256']


def test_multifidelity_iteration_retrains_with_prior_validated_cases():
    from pipeline.process.multifidelity_campaign import (
        MultiFidelityCampaignServices, run_multifidelity_iteration)
    from pipeline.process.transport_training import CaseArtifactReference

    model = _fit()
    prior = CaseArtifactReference(
        'prior', 'mmbcr', 'MMBCR', 350.0, 'training')
    observed = {}

    def execute(plan):
        return CaseArtifactReference(
            plan['case_id'], 'mmbcr', 'MMBCR', 400.0, plan['partition'])

    def train(references):
        observed['ids'] = [reference.candidate_id for reference in references]
        return {'status': 'published', 'pathway_mode': 'mmbcr',
                'reactor_type': 'MMBCR', 'model_sha256': model.sha256()}

    services = MultiFidelityCampaignServices(
        execute_case=execute, train_and_publish=train,
        load_published_model=lambda mode, reactor: model,
        append_lineage=lambda event_type, payload: {})
    result = run_multifidelity_iteration(
        designed_cases=[{'case_id': 'new', 'partition': 'training'}],
        screening_queries=[{
            'query_id': 'inside', 'region': 'class-a',
            'features': dict(_record('q', 400).features)}],
        services=services, referral_budget=1, prior_references=[prior])
    assert observed['ids'] == ['prior', 'new']
    assert result['prior_case_count'] == 1


def test_default_campaign_services_bind_portable_production_inputs():
    from pipeline.process.multifidelity_campaign import default_campaign_services

    observed = {}

    def solver_runner(**kwargs):
        observed['solver'] = kwargs

    def trainer(references, **kwargs):
        observed['training'] = (references, kwargs)
        return {'status': 'published'}

    with tempfile.TemporaryDirectory() as tmp:
        services = default_campaign_services(
            results_dir=Path(tmp) / 'results',
            registry_dir=Path(tmp) / 'models',
            ledger_path=Path(tmp) / 'ledger.json',
            targets=['gas_velocity_m_s'],
            validation_rmse_limits={'gas_velocity_m_s': 0.1},
            solver_runner=solver_runner, trainer=trainer)
        plan = {
            'candidate_id': 'candidate', 'pathway_mode': 'mmbcr',
            'reactor_type': 'MMBCR', 'temperature_K': 900,
            'case_dir': 'cases/mmbcr', 'model_source': 'reviewed-case-v1',
            'partition': 'training'}
        reference = services.execute_case(plan)
        services.train_and_publish([reference])
        assert observed['solver']['case_dir'] == 'cases/mmbcr'
        assert observed['solver']['results_dir'] == Path(tmp) / 'results'
        assert reference.partition == 'training'
        assert observed['training'][1]['model_registry'].root == (
            Path(tmp) / 'models').resolve()


TESTS = [value for name, value in sorted(globals().items())
         if name.startswith('test_') and callable(value)]


if __name__ == '__main__':
    failures = []
    for test in TESTS:
        try:
            test()
            print(f'PASS {test.__name__}')
        except Exception as exc:
            failures.append((test.__name__, exc))
            print(f'FAIL {test.__name__}: {type(exc).__name__}: {exc}')
    print(f'\n{len(TESTS)-len(failures)}/{len(TESTS)} modular contracts passed')
    raise SystemExit(1 if failures else 0)
