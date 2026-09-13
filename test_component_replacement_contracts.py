#!/usr/bin/env python3
"""Integration contracts proving pipeline components are independently replaceable."""

from dataclasses import replace
import inspect
from pathlib import Path

import pandas as pd

from pipeline.orchestrator import PipelineConfig, run_pipeline
from pipeline.stages.contracts import StageOutcome, require_stage_outcome
from pipeline.stages.orchestration import (
    MemoryStateStore, PipelineComponents, PipelineRuntime,
    default_pipeline_components)


def _raises(message, exception_type=Exception):
    class Raises:
        def __enter__(self):
            return self

        def __exit__(self, kind, value, traceback):
            assert kind is not None and issubclass(kind, exception_type)
            assert message in str(value), (message, str(value))
            return True
    return Raises()


def _runtime(store, calls):
    tick = iter(range(1000))
    return PipelineRuntime(
        clock=lambda: float(next(tick)),
        banner=lambda value: calls.append(('banner', value)),
        load_state=store.load, save_state=store.save,
        select_pathway_mode=lambda value: calls.append(('mode', value)))


def _never(name):
    def fail(*args, **kwargs):
        raise AssertionError(f'unselected component executed: {name}')
    return fail


def _base_components():
    return PipelineComponents(
        discovery=_never('discovery'), reactor_batch=_never('reactor_batch'),
        dft=_never('dft'), vqe=_never('vqe'), fuel_cell=_never('fuel_cell'),
        report=_never('report'), load_candidates=_never('load_candidates'))


def test_each_pipeline_phase_runs_with_only_its_replacement_and_contract_inputs():
    candidates = pd.DataFrame([{
        'candidate_id': 'candidate-1', 'genome': "('SAC', 'Fe', 'N4')",
        'E_act': 0.7, 'valid': True, 'material_class': 'SAC'}])
    phase_calls = []

    def discovery(**kwargs):
        phase_calls.append((1, kwargs))
        return StageOutcome({'ok': True}, {
            'design_space_sizes': {'TOTAL': 1}, 'pareto_genomes': [],
            'screening_database': candidates, 'top_catalysts': candidates,
            'dft_candidates': candidates})

    def loader(path, **kwargs):
        phase_calls.append(('load', Path(path).name, kwargs))
        return {'screening_database': candidates, 'top_catalysts': candidates,
                'dft_candidates': candidates}

    def reactor(rows, **kwargs):
        phase_calls.append((2, len(rows), kwargs))
        return StageOutcome({'n_simulations': 1}, {'reactor_results': [{}]})

    def dft(rows, **kwargs):
        phase_calls.append((3, len(rows), kwargs))
        return StageOutcome({'n_validated': 1}, {'dft_results': [{}]})

    def vqe(**kwargs):
        phase_calls.append((4, kwargs))
        return StageOutcome({'n_vqe_runs': 1}, {'vqe_results': [{}]})

    def fuel(**kwargs):
        phase_calls.append((5, kwargs))
        return StageOutcome({'n_valid': 0}, {
            'cathode_database': candidates, 'valid_cathodes': candidates,
            'pemfc_results': [], 'stack_result': {}})

    def report(state):
        phase_calls.append((6, dict(state)))
        return StageOutcome({'report_path': '/portable/report.md'},
                            {'report_path': Path('/portable/report.md')})

    replacements = {
        1: {'discovery': discovery},
        2: {'reactor_batch': reactor, 'load_candidates': loader},
        3: {'dft': dft, 'load_candidates': loader},
        4: {'vqe': vqe},
        5: {'fuel_cell': fuel},
        6: {'report': report},
    }
    for phase, changes in replacements.items():
        store, calls = MemoryStateStore(), []
        components = replace(_base_components(), **changes)
        state = run_pipeline(
            PipelineConfig(run_dft=False, run_vqe=False),
            start_phase=phase, end_phase=phase,
            runtime=_runtime(store, calls), components=components)
        assert f'phase{phase}' in state
        assert store.state == state
    assert any(call[0] == 'load' for call in phase_calls)
    assert {call[0] for call in phase_calls if isinstance(call[0], int)} == set(range(1, 7))


def test_replacement_contract_rejects_wrong_type_and_missing_products_early():
    with _raises('must return StageOutcome', TypeError):
        require_stage_outcome({}, stage='replacement')
    with _raises('omitted required products', ValueError):
        require_stage_outcome(
            StageOutcome({}, {}), stage='replacement',
            required_products=('scientific_result',))


def test_bad_stage_replacement_fails_at_named_boundary_before_persistence():
    store, calls = MemoryStateStore(), []
    components = replace(_base_components(), vqe=lambda **kwargs: {'bad': True})
    with _raises('vqe stage must return StageOutcome', TypeError):
        run_pipeline(
            PipelineConfig(run_vqe=False), start_phase=4, end_phase=4,
            runtime=_runtime(store, calls), components=components)
    assert store.history == []


def test_every_stage_boundary_rejects_invalid_replacement_output():
    candidates = pd.DataFrame([{
        'candidate_id': 'candidate-1', 'genome': "('SAC', 'Fe', 'N4')",
        'E_act': 0.7, 'valid': True, 'material_class': 'SAC'}])
    loader = lambda *args, **kwargs: {
        'screening_database': candidates, 'top_catalysts': candidates,
        'dft_candidates': candidates}
    cases = {
        1: ('discovery', {}), 2: ('reactor_batch', {'load_candidates': loader}),
        3: ('dft', {'load_candidates': loader}), 4: ('vqe', {}),
        5: ('fuel_cell', {}), 6: ('report', {}),
    }
    for phase, (stage, supporting) in cases.items():
        store, calls = MemoryStateStore(), []
        changes = {stage: lambda *args, **kwargs: {'invalid': True}, **supporting}
        with _raises(f'{stage} stage must return StageOutcome', TypeError):
            run_pipeline(
                PipelineConfig(run_dft=False, run_vqe=False),
                start_phase=phase, end_phase=phase,
                runtime=_runtime(store, calls),
                components=replace(_base_components(), **changes))
        assert store.history == []


def test_default_component_graph_conforms_to_orchestrator_call_signatures():
    components = default_pipeline_components()
    expected = {
        'discovery': {'initial_samples', 'leaf_size', 'max_leaves',
                      'top_k_reactor', 'top_k_dft'},
        'reactor_batch': {'candidates', 'temperatures', 'reactor_types',
                          'pathway_mode', 'multiphysics_results_dir',
                          'allow_mock_inputs'},
        'dft': {'candidates', 'top_k', 'execute_dft'},
        'vqe': {'top_k', 'execute_quantum'},
        'fuel_cell': {'top_k_pemfc', 'stack_cells'},
        'report': {'pipeline_state'},
        'load_candidates': {'path', 'top_k_reactor', 'top_k_dft'},
    }
    for name, required in expected.items():
        parameters = set(inspect.signature(getattr(components, name)).parameters)
        assert required.issubset(parameters), (name, required - parameters)


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
    print(f'\n{len(TESTS)-len(failures)}/{len(TESTS)} replacement contracts passed')
    raise SystemExit(1 if failures else 0)
