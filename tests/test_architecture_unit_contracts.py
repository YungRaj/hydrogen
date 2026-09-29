#!/usr/bin/env python3
"""Unit contracts for modular adapters, factories, and persistence boundaries."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from pipeline.stages.candidate_io import load_selected_candidates
from pipeline.stages.orchestration import MemoryStateStore


def _raises(message, exception_type=Exception):
    class Raises:
        def __enter__(self):
            return self

        def __exit__(self, kind, value, traceback):
            assert kind is not None and issubclass(kind, exception_type)
            assert message in str(value), (message, str(value))
            return True
    return Raises()


def test_memory_state_store_deeply_isolates_initial_load_save_and_history():
    initial = {'phase': {'values': [1]}}
    store = MemoryStateStore(initial)
    initial['phase']['values'].append(2)
    assert store.state == {
        'schema_version': 1, 'phase': {'values': [1]}}
    loaded = store.load()
    loaded['phase']['values'].append(3)
    assert store.state == {
        'schema_version': 1, 'phase': {'values': [1]}}
    store.save(loaded)
    loaded['phase']['values'].append(4)
    assert store.state == {
        'schema_version': 1, 'phase': {'values': [1, 3]}}
    store.state['phase']['values'].append(5)
    assert store.history == [{
        'schema_version': 1, 'phase': {'values': [1, 3]}}]


def test_candidate_file_adapter_preserves_distinct_reactor_and_validation_routes():
    rows = pd.DataFrame([
        {'candidate_id': 'complete', 'material_class': 'SAC', 'valid': True,
         'E_act': 0.5, 'E_act_censored': False},
        {'candidate_id': 'unresolved', 'material_class': 'MOF', 'valid': False,
         'E_act': float('nan'), 'E_act_censored': False},
        {'candidate_id': 'hazard', 'material_class': 'HEA', 'valid': False,
         'E_act': float('nan'), 'error': 'toxic/radioactive element'},
    ])
    with tempfile.TemporaryDirectory() as tmp:
        missing = Path(tmp) / 'missing.csv'
        assert load_selected_candidates(
            missing, top_k_reactor=2, top_k_dft=2) is None
        source = Path(tmp) / 'candidates.csv'
        rows.to_csv(source, index=False)
        restored = load_selected_candidates(
            source, top_k_reactor=2, top_k_dft=2)
    reactor_ids = set(restored['top_catalysts']['candidate_id'])
    validation_ids = set(restored['dft_candidates']['candidate_id'])
    assert reactor_ids == {'complete'}
    assert 'unresolved' in validation_ids
    assert 'hazard' not in validation_ids


def test_typed_contracts_preserve_runtime_compatibility_and_fail_early():
    """Typed boundaries must retain mappings while rejecting bad config."""
    from dataclasses import FrozenInstanceError

    from pipeline.data_models.core import (
        CandidateId, ConvergenceStatus, EvidenceLevel)
    from pipeline.orchestrator import PipelineConfig
    from pipeline.stages.contracts import StageOutcome, require_stage_outcome

    outcome = StageOutcome(
        state={'n_vqe_runs': 1}, products={'vqe_results': [{'ok': True}]})
    checked = require_stage_outcome(
        outcome, stage='vqe', required_products=('vqe_results',))
    assert checked.products['vqe_results'][0]['ok'] is True
    with _raises('cannot assign', FrozenInstanceError):
        outcome.state = {}
    assert CandidateId('candidate-1') == 'candidate-1'
    assert str(EvidenceLevel.DFT) == 'dft'
    assert ConvergenceStatus.CONVERGED.value == 'converged'
    with _raises('reactor temperatures must be positive kelvin', ValueError):
        PipelineConfig(reactor_temperatures=(0.0,))
    with _raises('requires reactors', ValueError):
        PipelineConfig(reactor_types=('MMBCR',))


def test_data_models_are_the_initial_strict_type_checking_boundary():
    """The checked boundary must be strict and expand deliberately over time."""
    config = json.loads((REPO_ROOT / 'pyrightconfig.json').read_text())
    assert config['typeCheckingMode'] == 'strict'
    assert config['pythonVersion'] == '3.10'
    assert 'pipeline/data_models' in config['include']
    assert 'pipeline/stages/contracts.py' in config['include']


def test_pipeline_state_and_scientific_result_models_trace_real_boundaries():
    """Named models must remain compatible with validated runtime records."""
    from pipeline.data_models.campaigns import validate_pipeline_state
    from pipeline.validation.qe_workflows import partial_hessian

    state = validate_pipeline_state({
        'phase1': {'total_evaluated': 4}, 'total_elapsed_s': 1})
    assert state['schema_version'] == 1
    assert state['phase1']['total_evaluated'] == 4
    assert state['total_elapsed_s'] == 1.0
    with _raises('phase2 pipeline state must be a mapping', ValueError):
        validate_pipeline_state({'phase2': []})
    with _raises('must be finite and nonnegative', ValueError):
        validate_pipeline_state({'total_elapsed_s': float('nan')})
    with _raises('newer than supported', ValueError):
        validate_pipeline_state({'schema_version': 999})
    with _raises('phase2 elapsed_s must be finite', ValueError):
        validate_pipeline_state({'phase2': {'elapsed_s': -1}})

    # A positive diagonal Hessian has no imaginary transition-state mode. The
    # named result makes that scientific distinction visible to every caller.
    import numpy as np
    plus = np.zeros((3, 1, 3))
    minus = np.zeros((3, 1, 3))
    result = partial_hessian(plus, minus, 0.01, np.asarray([1.0]))
    assert result['imaginary_count'] == 0
    assert result['valid_transition_state'] is False


def test_transport_and_fuel_cell_payloads_have_named_component_contracts():
    """High-value numerical handoffs must not regress to anonymous dictionaries."""
    from typing import get_type_hints

    from pipeline.data_models.fuel_cells import PEMFCResult, StackResult
    from pipeline.data_models.stages import FuelCellProducts
    from pipeline.data_models.transport import (
        ReactorClosure, TransportModelDocument, TransportPrediction,
        TransportTrainingResult)
    from pipeline.fuel_cell.pemfc import simulate_pemfc
    from pipeline.fuel_cell.stack import model_stack
    from pipeline.transport.closure_provider import resolve_reactor_closure
    from pipeline.transport.surrogate import TransportSurrogate
    from pipeline.transport.training import train_and_publish_transport_model

    assert get_type_hints(TransportSurrogate.predict)['return'] is TransportPrediction
    assert get_type_hints(TransportSurrogate.to_dict)['return'] is TransportModelDocument
    assert get_type_hints(resolve_reactor_closure)['return'] is ReactorClosure
    assert get_type_hints(train_and_publish_transport_model)['return'] is \
        TransportTrainingResult
    assert get_type_hints(simulate_pemfc)['return'] is PEMFCResult
    assert get_type_hints(model_stack)['return'] is StackResult
    assert get_type_hints(FuelCellProducts)['pemfc_results'] == list[PEMFCResult]

    assert {'usable', 'decision', 'reason',
            'candidate_exclusion_authorized'} <= TransportPrediction.__required_keys__
    assert {'predictions', 'uncertainty_1sigma',
            'validation_rmse'} <= TransportPrediction.__optional_keys__
    assert {'peak_power_W_cm2', 'peak_voltage_V',
            'current_density'} <= PEMFCResult.__required_keys__
    assert {'net_power_kW', 'system_efficiency',
            'cost_per_kW'} <= StackResult.__required_keys__


def test_discovery_reactor_and_solver_artifacts_have_named_contracts():
    """Selection, scorecard, and solver evidence retain traceable identities."""
    from typing import get_type_hints

    from pipeline.data_models.artifacts import (
        ArtifactValidationResult, ExternalSolverArtifact,
        InvalidArtifactResult, ValidArtifactResult)
    from pipeline.data_models.discovery import AdmissibilitySummary
    from pipeline.data_models.reactors import SolidsMetricRow, SolidsScorecard
    from pipeline.reactors.scorecard import build_solids_scorecard, metric_row
    from pipeline.search.scope import scope_pyrolysis_pool
    from pipeline.simulation.result_contract import load_validated_artifact

    assert get_type_hints(scope_pyrolysis_pool)['return'] == \
        tuple[get_type_hints(scope_pyrolysis_pool)['df'], AdmissibilitySummary]
    assert get_type_hints(metric_row)['return'] is SolidsMetricRow
    assert get_type_hints(build_solids_scorecard)['return'] is SolidsScorecard
    assert get_type_hints(load_validated_artifact)['return'] == \
        ArtifactValidationResult
    assert {'filter'} <= AdmissibilitySummary.__required_keys__
    assert {'judge_reason', 'headline', 'mmbcr_max_conversion'} <= \
        SolidsScorecard.__required_keys__
    assert {'candidate_id', 'outputs', 'provenance', 'surrogate_inputs'} <= \
        ExternalSolverArtifact.__required_keys__
    assert ValidArtifactResult.__required_keys__ == {'valid', 'artifact', 'path'}
    assert {'valid', 'reason'} <= InvalidArtifactResult.__required_keys__

    missing = load_validated_artifact(
        None, 'candidate', 'mmbcr', 'MMBCR', 1000.0)
    assert missing == {
        'valid': False,
        'reason': 'multiphysics_results_dir_not_configured',
    }


def test_multifidelity_plans_queries_and_decisions_are_explicit():
    """Active-learning control flow must expose stable, named record shapes."""
    from typing import get_type_hints

    from pipeline.campaigns.full_physics import schedule_full_physics_cases
    from pipeline.campaigns.multifidelity import run_multifidelity_iteration
    from pipeline.data_models.multifidelity import (
        DesignedCase, FullPhysicsSelection, MultiFidelityIterationResult,
        RepresentativeCase, ScreeningDecision, ScreeningQuery)
    from pipeline.transport.representative_cases import (
        assign_case_partitions, design_representative_cases)

    assert get_type_hints(schedule_full_physics_cases)['return'] == \
        list[FullPhysicsSelection]
    assert get_type_hints(run_multifidelity_iteration)['return'] is \
        MultiFidelityIterationResult
    assert get_type_hints(design_representative_cases)['return'] == \
        list[RepresentativeCase]
    assert get_type_hints(assign_case_partitions)['return'] == list[DesignedCase]
    assert {'case_id', 'partition'} <= DesignedCase.__required_keys__
    assert {'query_id', 'region', 'features'} == ScreeningQuery.__required_keys__
    assert {'query_id', 'region', 'input_sha256'} <= \
        ScreeningDecision.__required_keys__
    assert {'screening_decisions', 'scheduled_referrals',
            'lineage_event_sha256'} <= MultiFidelityIterationResult.__required_keys__


def test_all_default_service_factories_produce_callable_boundaries():
    from pipeline.simulation.reactor_handoff import default_reactor_coupling_services
    from pipeline.simulation.external_runner import default_solver_execution_services
    from pipeline.stages.discovery import default_discovery_services
    from pipeline.stages.fuel_cell import default_fuel_cell_services
    from pipeline.stages.orchestration import default_pipeline_components
    from pipeline.stages.reactor import default_reactor_services
    from pipeline.stages.reactor_batch import default_reactor_batch_services

    bundles = [
        default_pipeline_components(), default_discovery_services(),
        default_reactor_batch_services(), default_reactor_services(),
        default_fuel_cell_services(), default_reactor_coupling_services(),
        default_solver_execution_services(),
    ]
    for bundle in bundles:
        for name, value in vars(bundle).items():
            if (type(bundle).__name__, name) == (
                    'ReactorCouplingServices', 'load_surrogate'):
                assert value is None or callable(value)
            else:
                assert callable(value), (type(bundle).__name__, name, value)


def test_reactor_model_decomposition_preserves_the_public_boundary():
    """Focused reactor modules must retain the established import surface."""
    from pipeline.reactors import models
    from pipeline.reactors import reactor_core
    from pipeline.reactors.fluidized_bed import simulate_fluidized_bed
    from pipeline.reactors.mmbcr import simulate_mmbcr
    from pipeline.reactors.pfr import simulate_pfr

    assert models.ReactorConfig is reactor_core.ReactorConfig
    assert models.DEFAULT_SOLIDS_PARTICLE_MM == (
        reactor_core.DEFAULT_SOLIDS_PARTICLE_MM)
    for implementation in (
            simulate_mmbcr, simulate_pfr, simulate_fluidized_bed):
        assert callable(implementation)
    # The stable facade owns routing and sweep persistence; algorithms do not.
    assert models.simulate_reactor.__module__ == 'pipeline.reactors.models'
    assert models.run_reactor_sweep.__module__ == 'pipeline.reactors.models'


def test_every_pipeline_module_imports_without_starting_external_processes():
    script = r'''
import importlib, json, pkgutil, subprocess
def forbidden(*args, **kwargs):
    raise AssertionError('module import launched an external process')
subprocess.run = forbidden
subprocess.Popen = forbidden
import pipeline
names = [entry.name for entry in pkgutil.walk_packages(
    pipeline.__path__, pipeline.__name__ + '.')]
for name in names:
    importlib.import_module(name)
print(json.dumps({'count': len(names)}))
'''
    completed = subprocess.run(
        [sys.executable, '-c', script], capture_output=True, text=True,
        cwd=REPO_ROOT, timeout=60)
    assert completed.returncode == 0, completed.stderr
    result = json.loads(completed.stdout.strip().splitlines()[-1])
    assert result['count'] >= 80


def test_architecture_visual_atlas_is_linked_and_structurally_complete():
    root = REPO_ROOT
    atlas = (root / 'docs' / 'ARCHITECTURE_DIAGRAMS.md').read_text()
    assert atlas.count('```mermaid') == 9
    assert atlas.count('```') == 18
    for section in (
            'End-to-end catalyst discovery engine',
            'Divide-and-conquer traversal',
            'Multi-fidelity learning and referral loop',
            'Scientific software responsibilities',
            'Methane-conversion mode routing',
            'Evidence authority ladder',
            'Replaceable component architecture',
            'Provenance chain', 'Compute allocation and escalation'):
        assert f'## {section}' in atlas
    assert 'ARCHITECTURE_DIAGRAMS.md' in (root / 'README.md').read_text()
    rendered_locations = {
        root / 'README.md': 2,
        root / 'docs' / 'TECHNICAL_ARCHITECTURE.md': 1,
        root / 'docs' / 'MODULAR_MULTIFIDELITY_WORKFLOW.md': 1,
    }
    for document, minimum_diagrams in rendered_locations.items():
        text = document.read_text()
        assert 'ARCHITECTURE_DIAGRAMS.md' in text
        assert text.count('```mermaid') >= minimum_diagrams
        assert text.count('```') % 2 == 0


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
    print(f'\n{len(TESTS)-len(failures)}/{len(TESTS)} architecture unit contracts passed')
    raise SystemExit(1 if failures else 0)
