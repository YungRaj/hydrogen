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
    assert store.state == {'phase': {'values': [1]}}
    loaded = store.load()
    loaded['phase']['values'].append(3)
    assert store.state == {'phase': {'values': [1]}}
    store.save(loaded)
    loaded['phase']['values'].append(4)
    assert store.state == {'phase': {'values': [1, 3]}}
    store.state['phase']['values'].append(5)
    assert store.history == [{'phase': {'values': [1, 3]}}]


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
