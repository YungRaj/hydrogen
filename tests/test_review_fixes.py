"""Focused regression tests for the code-review fixes (search/screening/evidence).

Plain pytest functions. No GPU, fairchem model, torch or network access is
required: heavy dependencies are faked with monkeypatch or skipped through
pytest.importorskip. Imports live inside each test.
"""

# ======================================================================
# pipeline/screening review fixes
# ======================================================================

def _install_fake_modules(monkeypatch, modules):
    """Register lightweight stand-in modules (e.g. fairchem, torch)."""
    import sys
    import types

    for name, attrs in modules.items():
        module = types.ModuleType(name)
        for key, value in attrs.items():
            setattr(module, key, value)
        monkeypatch.setitem(sys.modules, name, module)


# ── #1 mock fallback rows are never valid ────────────────────────────────────
def test_cathode_mock_rows_are_invalid_placeholders():
    import pytest

    pytest.importorskip('numpy')
    pytest.importorskip('pandas')
    pytest.importorskip('ase')
    from pipeline.screening.fc_cathode_screener import (
        _mock_orr_result,
        generate_fc_catalyst_list,
    )

    for candidate in generate_fc_catalyst_list()[:5]:
        row = _mock_orr_result(candidate)
        assert row['mock'] is True
        assert row['valid'] is False
        assert row['candidate_disposition'] == 'validation_required'
        # #2: mock rows also carry the identity columns the stage reads.
        assert row['genome'] == str(candidate['genome'])
        assert row['material_class'] == candidate['genome'][0]


# ── #2 genome/material_class columns + #3 unconverged relaxation fails closed
def test_cathode_unconverged_relaxation_is_invalid_and_has_identity(monkeypatch):
    import pytest

    np = pytest.importorskip('numpy')
    pytest.importorskip('pandas')
    pytest.importorskip('ase')
    from ase import Atoms
    from ase.calculators.calculator import Calculator, all_changes

    from pipeline.screening import fc_cathode_screener, surface_screener

    class ConstantForce(Calculator):
        """E = -sum(x): a uniform +x force that no optimizer can relax."""

        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy',),
                      system_changes=all_changes):
            Calculator.calculate(self, atoms, properties, system_changes)
            positions = np.asarray(self.atoms.positions)
            forces = np.zeros_like(positions)
            forces[:, 0] = 1.0
            self.results = {
                'energy': float(-positions[:, 0].sum()),
                'forces': forces,
            }

    def fake_structure(genome, structure_seed=0):
        atoms = Atoms('Pt', positions=[[5.0, 5.0, 5.0]], cell=[10, 10, 10])
        return atoms, [0], genome[0]

    monkeypatch.setattr(surface_screener, 'generate_structure', fake_structure)
    candidate = fc_cathode_screener.generate_fc_catalyst_list()[0]
    result = fc_cathode_screener.screen_orr_candidate(
        candidate, ConstantForce(), e_h2o=-14.0, e_h2=-6.8
    )
    assert result['valid'] is False
    assert result['relax_clean_converged'] is False
    assert result['candidate_disposition'] == 'validation_required'
    assert 'orr_overpotential_V' not in result
    assert result['genome'] == str(candidate['genome'])
    assert result['material_class'] == candidate['genome'][0]


# ── #4 cache digest covers every module that shapes cached fields ───────────
def test_result_cache_digest_covers_descriptor_modules():
    import ast
    from pathlib import Path

    import pipeline.screening as screening_pkg

    source = Path(screening_pkg.__file__).with_name('result_cache.py').read_text()
    tree = ast.parse(source)
    modules = None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == '_SCIENTIFIC_MODULES'
            for t in node.targets
        ):
            modules = set(ast.literal_eval(node.value))
    assert modules is not None
    for required in (
        'pipeline.utils',
        'pipeline.screening.ood',
        'pipeline.search.scope',
        'pipeline.search.design_space',
        'pipeline.screening.relaxation',
        'pipeline.screening.protocols',
    ):
        assert required in modules


# ── #5 no silent model substitution for an explicitly requested model ───────
def _fake_fairchem(monkeypatch, tmp_path, calls):
    from pipeline.screening import surface_calculator

    def failing_unit(*args, **kwargs):
        raise RuntimeError('gated checkpoint unavailable')

    def fake_ocp(**kwargs):
        calls.append(kwargs)
        return 'equiformer-v2'

    _install_fake_modules(monkeypatch, {
        'fairchem': {},
        'fairchem.core': {
            'FAIRChemCalculator': object,
            'OCPCalculator': fake_ocp,
        },
        'fairchem.core.calculate': {},
        'fairchem.core.calculate.pretrained_mlip': {
            'get_predict_unit': failing_unit,
        },
    })
    checkpoint = tmp_path / 'models' / 'eq2_31M_ec4_allmd.pt'
    checkpoint.parent.mkdir(parents=True)
    checkpoint.write_bytes(b'weights')
    monkeypatch.setattr(surface_calculator, 'BASE_DIR', tmp_path)
    monkeypatch.delenv('FAIRCHEM_CHECKPOINT_SHA256', raising=False)
    monkeypatch.setenv('HF_TOKEN', 'test-token')
    return surface_calculator


def test_explicit_model_failure_returns_none_without_fallback(monkeypatch, tmp_path):
    import pytest

    pytest.importorskip('numpy')
    pytest.importorskip('ase')
    calls = []
    surface_calculator = _fake_fairchem(monkeypatch, tmp_path, calls)
    assert surface_calculator.get_ocp_calculator(
        'esen-sm-conserving-all-oc25', device='cpu'
    ) is None
    assert calls == []


def test_unnamed_model_keeps_legacy_fallback(monkeypatch, tmp_path):
    import pytest

    pytest.importorskip('numpy')
    pytest.importorskip('ase')
    calls = []
    surface_calculator = _fake_fairchem(monkeypatch, tmp_path, calls)
    assert surface_calculator.get_ocp_calculator(device='cpu') == 'equiformer-v2'
    assert len(calls) == 1


def test_gpu_worker_reports_fatal_with_pid_when_model_missing(monkeypatch):
    import os
    import queue
    import threading

    import pytest

    pytest.importorskip('numpy')
    pytest.importorskip('pandas')
    pytest.importorskip('ase')
    from pipeline.screening import gpu_executor, surface_calculator

    _install_fake_modules(monkeypatch, {
        'torch': {
            'set_num_threads': lambda n: None,
            'set_num_interop_threads': lambda n: None,
        },
    })
    monkeypatch.setattr(
        surface_calculator, 'get_ocp_calculator', lambda **kwargs: None
    )
    # run_gpu_worker mutates these; monkeypatch restores them afterwards.
    for name in ('CUDA_VISIBLE_DEVICES', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS',
                 'NUMEXPR_NUM_THREADS'):
        monkeypatch.setenv(name, os.environ.get(name, '1'))

    class Logger:
        def error(self, *args):
            pass

    statuses = queue.Queue()
    gpu_executor.run_gpu_worker(
        0, 0, 'GPU-fake', queue.Queue(), statuses, threading.Event(), 1,
        False, lambda calc: None, 'protocol', Logger(), 'test worker',
    )
    kind, worker_id, _, payload = statuses.get_nowait()
    assert kind == 'fatal' and worker_id == 0
    assert payload['pid'] == os.getpid()
    assert 'could not be loaded' in payload['error']


# ── #6 stale fatal from a replaced process must not kill the replacement ────
def test_supervisor_ignores_stale_fatal_after_exit_restart(tmp_path):
    import json
    import queue
    import threading

    from pipeline.screening.worker_supervisor import collect_results, emit

    class FakeProcess:
        def __init__(self, pid, alive):
            self.pid, self.alive, self.exitcode = pid, alive, None

        def is_alive(self):
            return self.alive

        def terminate(self):
            self.alive, self.exitcode = False, -15

        def join(self, timeout=None):
            pass

    statuses, tasks, stopped = queue.Queue(), queue.Queue(), threading.Event()
    tasks.put((0, ('SAC',)))
    # The original worker (pid 101) already exited; its 'fatal' event is
    # queued behind another message, so the liveness check sees the exit first.
    initial = FakeProcess(101, alive=False)
    initial.exitcode = 1
    emit(statuses, 'ready', 0)
    emit(statuses, 'fatal', 0, {'pid': 101, 'error': 'init failed'})

    def restart(worker_id):
        replacement = FakeProcess(202, alive=True)

        def complete():
            index, genome = tasks.get(timeout=2)
            emit(statuses, 'ready', worker_id)
            emit(statuses, 'started', worker_id, index)
            emit(statuses, 'result', worker_id,
                 (index, {'genome': repr(genome), 'valid': True}))

        threading.Thread(target=complete, daemon=True).start()
        return replacement

    results = collect_results(
        statuses, tasks, stopped, {0: initial}, restart, [('SAC',)],
        'contract', tmp_path / 'health.json',
        startup_timeout_s=5, heartbeat_timeout_s=5, result_timeout_s=5,
    )
    assert results[0]['valid'] is True
    payload = json.loads((tmp_path / 'health.json').read_text())
    kinds = [event['kind'] for event in payload['events']]
    assert kinds.count('worker_restarted') == 1
    assert 'stale_fatal_ignored' in kinds


def test_supervisor_fatal_payload_parsing_keeps_legacy_strings():
    from pipeline.screening.worker_supervisor import _fatal_details

    assert _fatal_details('boom') == ('boom', None)
    assert _fatal_details({'pid': 7, 'error': 'boom'}) == ('boom', 7)


# ── #7 batched inference isolates failures per request ──────────────────────
def _fake_batch_service(monkeypatch):
    import pytest

    np = pytest.importorskip('numpy')
    pytest.importorskip('ase')
    from pipeline.screening.batched_calculator import BatchedInferenceService

    _install_fake_modules(monkeypatch, {
        'fairchem': {},
        'fairchem.core': {},
        'fairchem.core.datasets': {},
        'fairchem.core.datasets.atomic_data': {
            'atomicdata_list_to_batch': lambda data: list(data),
        },
    })

    class Tensor:
        def __init__(self, array):
            self.array = np.asarray(array, dtype=float)

        def detach(self):
            return self

        def cpu(self):
            return self

        def numpy(self):
            return self.array

    class Predictor:
        def validate_atoms_data(self, atoms, task_name):
            if 'Xe' in atoms.get_chemical_symbols():
                raise ValueError('unsupported element')

        def predict(self, batch):
            # Any batch containing a 3-atom system fails inference.
            if any(len(atoms) == 3 for atoms in batch):
                raise RuntimeError('poisoned batch')
            energies = [float(len(atoms)) for atoms in batch]
            forces = np.vstack([np.zeros((len(atoms), 3)) for atoms in batch])
            return {'energy': Tensor(energies), 'forces': Tensor(forces)}

    class FakeCalculator:
        implemented_properties = ['energy', 'free_energy', 'forces']
        task_name = 'oc25'
        predictor = Predictor()

        def _check_atoms_pbc(self, atoms):
            pass

        def a2g(self, atoms):
            return atoms

    return BatchedInferenceService(FakeCalculator(), batch_wait_ms=0)


def test_batched_validation_failure_does_not_poison_batch(monkeypatch):
    service = _fake_batch_service(monkeypatch)
    from ase import Atoms

    from pipeline.screening.batched_calculator import _Request

    try:
        good = _Request(Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]]))
        bad = _Request(Atoms('Xe', positions=[[0, 0, 0]]))
        service._run_batch([good, bad])
        assert good.error is None and good.result['energy'] == 2.0
        assert good.result['forces'].shape == (2, 3)
        assert isinstance(bad.error, ValueError)
        assert good.event.is_set() and bad.event.is_set()
    finally:
        service.close()


def test_batched_inference_failure_is_retried_per_request(monkeypatch):
    service = _fake_batch_service(monkeypatch)
    from ase import Atoms

    from pipeline.screening.batched_calculator import _Request

    try:
        good = _Request(Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]]))
        poison = _Request(
            Atoms('H3', positions=[[0, 0, 0], [0, 0, 0.74], [0, 0, 1.5]])
        )
        service._run_batch([good, poison])
        assert good.error is None and good.result['energy'] == 2.0
        assert isinstance(poison.error, RuntimeError)
        assert poison.result is None
    finally:
        service.close()


# ── #8 relaxation error message when every attempt raised ───────────────────
def test_require_relaxation_reports_optimizer_error_not_format_crash():
    import pytest

    pytest.importorskip('numpy')
    pytest.importorskip('ase')
    from ase import Atoms
    from ase.calculators.calculator import Calculator, all_changes

    from pipeline.screening.relaxation import require_relaxation

    class Broken(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy',),
                      system_changes=all_changes):
            raise RuntimeError('device lost')

    atoms = Atoms('H', positions=[[0.2, 0.0, 0.0]])
    atoms.calc = Broken()
    result = {'valid': True}
    assert require_relaxation(result, atoms, 'probe', 0.05, 5) is False
    assert result['valid'] is False
    assert result['relax_probe_max_force_eV_A'] is None
    assert 'device lost' in result['error']
    assert result['error'].startswith('Unconverged probe relaxation')
    assert result['candidate_disposition'] == 'validation_required'


# ── #9 one Fenton table for every ORR screener ──────────────────────────────
def test_cathode_screener_uses_canonical_fenton_table():
    import inspect

    import pytest

    pytest.importorskip('numpy')
    pytest.importorskip('pandas')
    pytest.importorskip('ase')
    from pipeline.screening import fc_cathode_screener, fc_screener

    assert fc_cathode_screener.FENTON_RISK is fc_screener.FENTON_RISK
    source = inspect.getsource(fc_cathode_screener.screen_orr_candidate)
    assert 'FENTON_RISK.get' in source

# ======================================================================
# pipeline/evidence review fixes
# ======================================================================

def test_experimental_evidence_records_do_not_duplicate_manifest_identity(tmp_path, monkeypatch):
    import hashlib
    import json

    from pipeline.evidence import experimental_dataset
    from pipeline.evidence.manifest import verify_evidence_manifest

    dataset = tmp_path / 'dataset.json'
    dataset.write_text('{"placeholder": true}\n')
    digest = hashlib.sha256(dataset.read_bytes()).hexdigest()
    report = {
        'dataset_id': 'ds-1',
        'protocol_id': 'protocol-1',
        'dataset_sha256': digest,
        'groups': [
            {'candidate_id': 'cand-1', 'measurement_type': 'ntec_pair', 'split_complete': True},
            {'candidate_id': 'cand-1', 'measurement_type': 'pyrolysis_reactor', 'split_complete': True},
            {'candidate_id': 'cand-2', 'measurement_type': 'pyrolysis_reactor', 'split_complete': True},
        ],
    }
    monkeypatch.setattr(experimental_dataset, 'validate_experimental_dataset', lambda path: report)

    records = experimental_dataset.experimental_evidence_records(dataset)

    assert [r['candidate_id'] for r in records['measured_reactor']] == ['cand-1', 'cand-2']
    assert [r['candidate_id'] for r in records['measured_deactivation']] == ['cand-1', 'cand-2']
    assert [r['candidate_id'] for r in records['ntec_control_pair']] == ['cand-1']
    manifest = tmp_path / 'evidence_manifest.json'
    manifest.write_text(json.dumps({'schema_version': 2, 'records': records}))
    verification = verify_evidence_manifest(manifest)
    assert verification['valid'], verification['errors']
    assert verification['counts']['measured_reactor'] == 2


def test_run_log_short_row_is_rejected_with_reason_error(tmp_path):
    import hashlib
    import json

    from pipeline.evidence.experimental_dataset import validate_experimental_dataset

    run_log = tmp_path / 'run_log.csv'
    # The data row omits failure_reason entirely; csv.DictReader yields None.
    run_log.write_text(
        'measurement_id,experimental_unit_id,split,status,failure_reason\n'
        'm-1,u-1,calibration,failed\n'
    )
    dataset = tmp_path / 'dataset.json'
    dataset.write_text(json.dumps({
        'schema_version': 1,
        'dataset_id': 'ds-1',
        'protocol_id': 'protocol-1',
        'records': [{'measurement_id': 'm-1'}],
        'run_log_path': run_log.name,
        'run_log_sha256': hashlib.sha256(run_log.read_bytes()).hexdigest(),
        'analysis_plan': {'plan_id': 'p'},
    }))
    try:
        validate_experimental_dataset(dataset)
    except ValueError as exc:
        assert 'require a reason' in str(exc)
    else:
        raise AssertionError('short run-log row must be rejected')


def test_campaign_status_ready_requires_valid_manifest(tmp_path, monkeypatch):
    import json

    from pipeline.evidence import campaign_status

    for relative in (
        'screening/turquoise_hydrogen_coverage_certificate.json',
        'fuel_cell/coverage_certificate.json',
    ):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({
            'complete': True,
            'declared_encoded_population': campaign_status.TOTAL_SIZE,
        }))
    counts = {
        key: 1 for key in (
            'converged_dft', 'converged_orr_dft', 'measured_reactor',
            'measured_deactivation', 'measured_mea', 'measured_durability',
            'hydrogen_impurity_test', 'time_split_benchmark',
            'curated_prior_art_source', 'ntec_control_pair',
        )
    }

    def verification(valid):
        return lambda path: {
            'valid': valid,
            'errors': [] if valid else ['measured_mea[1]:checksum_mismatch'],
            'counts': counts,
        }

    monkeypatch.setattr(campaign_status, 'verify_evidence_manifest', verification(False))
    blocked = campaign_status.assess_campaign(str(tmp_path))
    assert all(blocked['criteria'].values())
    assert not blocked['ready']
    assert not blocked['evidence_manifest_valid']

    monkeypatch.setattr(campaign_status, 'verify_evidence_manifest', verification(True))
    assert campaign_status.assess_campaign(str(tmp_path))['ready']


def test_report_top_barriers_exclude_censored_and_non_viable(tmp_path, monkeypatch):
    import pandas as pd

    monkeypatch.chdir(tmp_path)
    from pipeline.evidence.report_generator import _screening_section

    frame = pd.DataFrame([
        {'genome': "('censored',)", 'valid': True, 'E_act': 0.01,
         'E_act_censored': True, 'pyrolysis_viable': True},
        {'genome': "('nonviable',)", 'valid': True, 'E_act': 0.20,
         'E_act_censored': False, 'pyrolysis_viable': False},
        {'genome': "('eligible',)", 'valid': True, 'E_act': 0.90,
         'E_act_censored': False, 'pyrolysis_viable': True},
    ])
    text = '\n'.join(_screening_section({'surface_screening_full': frame}))
    assert "('eligible',)" in text
    assert "('censored',)" not in text
    assert "('nonviable',)" not in text
    # Without the gate columns the legacy valid-only ranking is unchanged.
    legacy = frame.drop(columns=['E_act_censored', 'pyrolysis_viable'])
    legacy_text = '\n'.join(_screening_section({'surface_screening_full': legacy}))
    assert "('censored',)" in legacy_text


def test_time_split_recall_counts_unique_candidates_not_citations():
    from pipeline.evidence.novelty_benchmark import time_split_recovery

    recovered = ('SAC', 'Fe', 'N4', 'N-graphene', 'OH')
    missed = ('SAC', 'Ni', 'N4', 'N-graphene', 'none')
    training = [{
        'genome': ('SAC', 'Co', 'N4', 'N-graphene', 'none'),
        'publication_year': 2023, 'source_id': 'doi:train', 'citation': 'Train',
    }]
    held_out = [
        {'genome': recovered, 'publication_year': 2025, 'source_id': f'doi:{n}',
         'citation': f'Paper {n}'}
        for n in range(3)
    ] + [{
        'genome': missed, 'publication_year': 2025, 'source_id': 'doi:missed',
        'citation': 'Missed',
    }]
    result = time_split_recovery(
        [recovered], held_out, cutoff_year=2024, training_records=training, k=1
    )
    assert result['valid']
    assert result['n_held_out'] == 2
    assert result['exact_recovered'] == 1
    assert result['exact_recall_at_k'] == 0.5


def test_pilot_outcomes_empty_input_raises_value_error(tmp_path):
    import pandas as pd

    from pipeline.evidence.pilot_benchmark import PilotSpec, load_legacy_outcomes

    path = tmp_path / 'legacy.csv'
    pd.DataFrame([
        {'genome': repr(('SAC', 'Fe', 'N4', 'N-graphene', 'none')), 'valid': False, 'score': 0.4},
        {'genome': 'not a genome', 'valid': True, 'score': 0.5},
    ]).to_csv(path, index=False)
    try:
        load_legacy_outcomes(PilotSpec('fuel_cell_orr', (str(path),), 'score'))
    except ValueError as exc:
        assert 'no valid' in str(exc)
    else:
        raise AssertionError('empty pilot input must raise ValueError')


def test_pilot_outcomes_merge_canonically_equivalent_genomes(tmp_path):
    import pandas as pd

    from pipeline.evidence.pilot_benchmark import PilotSpec, load_legacy_outcomes
    from pipeline.search.discovery import candidate_id

    first = ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, ('Cu', 'Fe'), 1, 0)
    swapped = ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, ('Fe', 'Cu'), 1, 0)
    path = tmp_path / 'legacy.csv'
    pd.DataFrame([
        {'genome': repr(first), 'valid': True, 'score': 0.2},
        {'genome': repr(swapped), 'valid': True, 'score': 0.6},
    ]).to_csv(path, index=False)
    data = load_legacy_outcomes(PilotSpec('turquoise_hydrogen', (str(path),), 'score'))
    assert len(data) == 1
    assert abs(float(data.score.iloc[0]) - 0.4) < 1e-12
    assert int(data.replicates.iloc[0]) == 2
    assert data.candidate_id.iloc[0] == candidate_id(first)


def test_campaign_qualification_closes_archive_connections():
    import ast
    import inspect

    from pipeline.evidence import campaign_qualification

    source = inspect.getsource(campaign_qualification._search_resume_check)
    tree = ast.parse(source.lstrip())
    bare = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.withitem)
        and isinstance(node.context_expr, ast.Call)
        and ast.unparse(node.context_expr.func) == 'sqlite3.connect'
    ]
    assert not bare, 'sqlite3.connect used as a context manager does not close'
    assert 'closing(sqlite3.connect(' in source

# ======================================================================
# pipeline/search review fixes
# ======================================================================

def test_orr_validation_quota_skips_pemfc_out_of_scope_classes(tmp_path):
    """#11: MoltenMetal has no PEMFC cathode realization and gets no ORR slot."""
    import numpy as np
    from pipeline.search.adaptive_validation import allocate_validation_batch

    candidates = [
        ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000),
        ('SAC', 'Fe', 'N4', 'N-graphene', 'OH'),
    ]
    # MoltenMetal looks best, so it would win any unreserved slot.
    objectives = np.array([[0.1, 1.0], [0.9, 1.0]])
    db = str(tmp_path / 'quota.sqlite')
    # Before the fix both classes reserved a slot: quota 2 > budget 1 raised.
    selected = allocate_validation_batch(
        candidates, objectives, 1, db, 'fuel_cell_orr', min_per_class=1
    )
    assert [candidates[i][0] for i in selected] == ['SAC']
    # Non-PEMFC applications keep their previous quota semantics.
    both = allocate_validation_batch(
        candidates, objectives, 2, db, 'test', min_per_class=1
    )
    assert {candidates[i][0] for i in both} == {'MoltenMetal', 'SAC'}


def test_experimental_slate_region_tie_goes_to_smallest_candidate_id():
    """#15d: region champion tie-break matches the ascending-ID ranking sort."""
    import numpy as np
    from pipeline.search.adaptive_validation import experimental_slate
    from pipeline.search.discovery import candidate_id, discovery_region

    a = ('SAC', 'Fe', 'N4', 'N-graphene', 'OH')
    b = ('SAC', 'Fe', 'N4', 'N-graphene', 'Cl')
    assert discovery_region(a) == discovery_region(b)
    candidates = [a, b]
    objectives = np.array([[0.5, 1.0], [0.5, 1.0]])
    expected = min(range(2), key=lambda i: candidate_id(candidates[i]))
    assert experimental_slate(candidates, objectives, 1) == [expected]
    assert experimental_slate(candidates[::-1], objectives, 1) == [1 - expected]


def test_discovery_batch_region_tie_goes_to_smallest_candidate_id():
    """#15d: select_discovery_batch uses the same tie-break as its sort."""
    import numpy as np
    from pipeline.search.discovery import candidate_id, select_discovery_batch

    candidates = [
        ('SAC', 'Fe', 'N4', 'N-graphene', 'OH'),
        ('SAC', 'Fe', 'N4', 'N-graphene', 'Cl'),
    ]
    objectives = np.array([[0.5, 1.0], [0.5, 1.0]])
    expected = min(range(2), key=lambda i: candidate_id(candidates[i]))
    assert select_discovery_batch(candidates, objectives, 1) == [expected]
