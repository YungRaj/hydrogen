#!/usr/bin/env python3
"""
Hydrogen Pipeline — Comprehensive Test Suite

Run:  python tests/test_pipeline.py
Exit: 0 = all pass, 1 = failures

Tests every component that has broken before, plus integration across
all 14 material classes. If this passes, the campaign is safe to launch.
"""

import sys, os, time, traceback
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import torch

PASS = 0
FAIL = 0
ERRORS = []


def test(name, fn):
    """Run a test, print result, track pass/fail."""
    global PASS, FAIL
    try:
        fn()
        print(f"  ✅ {name}")
        PASS += 1
    except Exception as e:
        print(f"  ❌ {name}: {e}")
        ERRORS.append((name, traceback.format_exc()))
        FAIL += 1


# ═══════════════════════════════════════════════════════════════════════════════
# 1. DESIGN SPACE
# ═══════════════════════════════════════════════════════════════════════════════

def test_all_14_classes_generate():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome
    pop = [generate_random_genome() for _ in range(3000)]
    classes_seen = set(g[0] for g in pop)
    assert len(ALL_MATERIAL_CLASSES) == 14, f"Expected 14 classes, got {len(ALL_MATERIAL_CLASSES)}"
    assert classes_seen == set(ALL_MATERIAL_CLASSES), f"Missing: {set(ALL_MATERIAL_CLASSES) - classes_seen}"


def test_no_toxic_elements():
    from pipeline.common.catalyst_spaces import generate_random_genome
    from pipeline.common.utils import TOXIC_ELEMENTS
    pop = [generate_random_genome() for _ in range(2000)]
    for g in pop:
        for field in g[1:]:
            if isinstance(field, str) and field in TOXIC_ELEMENTS:
                raise AssertionError(f"Toxic element {field} in genome {g}")


def test_sac_has_axial_ligand():
    from pipeline.common.catalyst_spaces import generate_random_genome, SAC_AXIAL_LIGANDS
    sacs = [generate_random_genome('SAC') for _ in range(100)]
    assert all(len(g) == 5 for g in sacs), "SAC genome should have 5 fields"
    axials = set(g[4] for g in sacs)
    assert len(axials) > 3, f"Only {len(axials)} axial variants seen, expected >3"


def test_class_weights_sum_to_1():
    from pipeline.common.catalyst_spaces import CLASS_WEIGHTS
    total = sum(CLASS_WEIGHTS.values())
    assert abs(total - 1.0) < 0.01, f"Class weights sum to {total}, expected ~1.0"


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ENCODING
# ═══════════════════════════════════════════════════════════════════════════════

def test_encode_all_classes_no_nan():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome, encode_genome, FEATURE_DIM
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(50):
            g = generate_random_genome(cls)
            feat = encode_genome(g)
            assert feat.shape == (FEATURE_DIM,), f"{cls}: shape {feat.shape} != ({FEATURE_DIM},)"
            assert not np.any(np.isnan(feat)), f"{cls}: NaN in encoding"


def test_feature_dim_matches_components():
    from pipeline.common.catalyst_spaces import (
        FEATURE_DIM, N_CLASSES, N_METALS, N_SUPPORTS, N_FACETS,
        N_COORDS, N_DOPANTS, N_CONTINUOUS
    )
    expected = N_CLASSES + 2 * N_METALS + N_SUPPORTS + N_FACETS + N_COORDS + N_DOPANTS + N_CONTINUOUS
    assert FEATURE_DIM == expected, f"FEATURE_DIM={FEATURE_DIM} != computed {expected}"


# ═══════════════════════════════════════════════════════════════════════════════
# 3. SURROGATES
# ═══════════════════════════════════════════════════════════════════════════════

def test_ch4_surrogate_no_nan():
    from pipeline.common.catalyst_spaces import generate_random_genome, encode_genome, FEATURE_DIM
    from pipeline.screening.surrogate_model import CatalystSurrogate
    model = CatalystSurrogate(input_dim=FEATURE_DIM)
    model.eval()
    pop = [generate_random_genome() for _ in range(500)]
    X = torch.FloatTensor(np.array([encode_genome(g) for g in pop]))
    with torch.no_grad():
        out = model(X)
    for i, o in enumerate(out):
        assert not torch.any(torch.isnan(o)), f"CH4 surrogate head {i} has NaN"


def test_orr_surrogate_no_nan():
    from pipeline.common.catalyst_spaces import generate_random_genome, encode_genome, FEATURE_DIM
    from pipeline.screening.fc_genetic_optimizer import ORRCatalystSurrogate
    model = ORRCatalystSurrogate(input_dim=FEATURE_DIM)
    model.eval()
    pop = [generate_random_genome() for _ in range(500)]
    X = torch.FloatTensor(np.array([encode_genome(g) for g in pop]))
    with torch.no_grad():
        out = model(X)
    for i, o in enumerate(out):
        assert not torch.any(torch.isnan(o)), f"ORR surrogate head {i} has NaN"


# ═══════════════════════════════════════════════════════════════════════════════
# 4. NSGA-II
# ═══════════════════════════════════════════════════════════════════════════════

def test_nsga2_sorts_correctly():
    from pipeline.common.catalyst_spaces import generate_random_genome, FEATURE_DIM
    from pipeline.screening.fc_genetic_optimizer import (
        ORRCatalystSurrogate, compute_orr_objectives_surrogate, fast_non_dominated_sort
    )
    model = ORRCatalystSurrogate(input_dim=FEATURE_DIM)
    model.eval()
    pop = [generate_random_genome() for _ in range(500)]
    obj = compute_orr_objectives_surrogate(pop, model, 'cpu')
    fronts = fast_non_dominated_sort(obj)
    total = sum(len(f) for f in fronts)
    assert total == 500, f"NSGA-II lost genomes: {total} != 500"
    assert not np.any(np.isnan(obj)), "NaN in objectives"


def test_cost_and_fenton_ranges():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome
    from pipeline.screening.fc_genetic_optimizer import _cost_from_genome, _fenton_from_genome
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(20):
            g = generate_random_genome(cls)
            c = _cost_from_genome(g)
            f = _fenton_from_genome(g)
            assert 0 <= c <= 2.0, f"{cls}: cost {c} out of range"
            assert 0 <= f <= 10, f"{cls}: fenton {f} out of range"


def test_metalfreecarbon_zero_cost():
    from pipeline.common.catalyst_spaces import generate_random_genome
    from pipeline.screening.fc_genetic_optimizer import _cost_from_genome
    for _ in range(20):
        g = generate_random_genome('MetalFreeCarbon')
        c = _cost_from_genome(g)
        assert c == 0.0, f"MetalFreeCarbon cost should be 0, got {c}"


def test_pemfc_application_scope():
    from pipeline.common.application_scope import pemfc_cathode_scope
    assert pemfc_cathode_scope(('SAC', 'Fe'))['status'] == 'candidate'
    assert pemfc_cathode_scope(('MoltenMetal', 'Ga'))['status'] == 'out_of_scope'
    assert pemfc_cathode_scope(('MetalHydride', 'La'))['status'] == 'out_of_scope'


def test_novelty_time_split_benchmark():
    from pipeline.evidence.novelty_benchmark import time_split_recovery
    known = ('SAC', 'Fe', 'N4', 'N-graphene', 'OH')
    training = [{'genome': ('SAC', 'Co', 'N4', 'N-graphene', 'none'),
                 'publication_year': 2023, 'source_id': 'doi:train',
                 'citation': 'Training et al.'}]
    held_out = [{'genome': known, 'publication_year': 2025,
                 'source_id': 'doi:test', 'citation': 'Test et al.'}]
    result = time_split_recovery(
        [known], held_out, cutoff_year=2024,
        training_records=training, k=1)
    assert result['valid'] and result['exact_recall_at_k'] == 1.0
    malformed = time_split_recovery(
        [known], [{'genome': known}], cutoff_year=2024,
        training_records=training, k=1)
    assert not malformed['valid']
    leaked = time_split_recovery(
        [known], held_out, cutoff_year=2024,
        training_records=[dict(held_out[0], publication_year=2024)], k=1)
    assert not leaked['valid']


def test_pilot_benchmark_deduplicates_candidates():
    import tempfile
    from pathlib import Path
    import pandas as pd
    from pipeline.evidence.pilot_benchmark import PilotSpec, load_legacy_outcomes
    rows = [
        {'genome': repr(('SAC', 'Fe', 'N4', 'N-graphene', 'none')),
         'valid': True, 'score': value} for value in (0.4, 0.6)
    ]
    rows.append({'genome': repr(('SAC', 'Co', 'N4', 'N-graphene', 'none')),
                 'valid': True, 'score': 0.8})
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'legacy.csv'
        pd.DataFrame(rows).to_csv(path, index=False)
        data = load_legacy_outcomes(PilotSpec('fuel_cell_orr', (str(path),), 'score'))
    assert len(data) == 2
    assert sorted(data.score.tolist()) == [0.5, 0.8]
    assert data.candidate_id.nunique() == 2


def test_small_data_rankers_preserve_continuous_targets():
    import pandas as pd
    from pipeline.search.indexed_space import deterministic_tree_probes
    from pipeline.screening.small_data_ranker import (
        TREE_ENSEMBLE_SIZE, fit_tree_ranker, merge_compatible_evidence)
    genomes = deterministic_tree_probes(24)
    pyro = pd.DataFrame({'genome': [repr(g) for g in genomes], 'valid': True,
                         'E_act': np.linspace(0.1, 2.0, len(genomes))})
    ranker = fit_tree_ranker(pyro, 'turquoise_hydrogen')
    assert len(ranker.model.estimators_) == TREE_ENSEMBLE_SIZE == 256
    mean, uncertainty = ranker.predict(genomes[:5])
    assert np.all(np.isfinite(mean)) and np.all(uncertainty >= 0)
    fast_mean, omitted_uncertainty = ranker.predict(
        genomes[:5], uncertainty=False)
    assert np.allclose(mean, fast_mean) and np.all(omitted_uncertainty == 0)
    orr = pd.DataFrame({'genome': [repr(g) for g in genomes], 'valid': True,
                        'dG_OH_eV': np.linspace(-1, 1, len(genomes)),
                        'dG_O_eV': np.linspace(-2, 2, len(genomes)),
                        'dG_OOH_eV': np.linspace(3, 5, len(genomes))})
    ranker = fit_tree_ranker(orr, 'fuel_cell_orr')
    mean, _ = ranker.predict(genomes[:5])
    assert np.all(np.isfinite(mean)) and len(set(mean.tolist())) > 1

    current = pd.DataFrame({
        'genome': ['new', 'duplicate'],
        'screening_protocol': ['v2', 'v2'],
        'value': [1, 2],
    })
    ledger = pd.DataFrame({
        'genome': ['old', 'wrong', 'duplicate'],
        'screening_protocol': ['v2', 'v1', 'v2'],
        'value': [3, 4, 5],
    })
    merged = merge_compatible_evidence(current, ledger, 'v2')
    assert set(merged.genome) == {'old', 'new', 'duplicate'}
    assert merged.loc[merged.genome == 'duplicate', 'value'].iloc[0] == 2


def test_six_point_status_fails_closed():
    import tempfile
    from pipeline.evidence.campaign_status import assess_campaign
    with tempfile.TemporaryDirectory() as tmp:
        result = assess_campaign(tmp)
        assert not result['ready']
        assert 'calibrated_ntec' not in result['criteria']
        assert len(result['missing']) == 5
        ntec = assess_campaign(tmp, pyrolysis_mode='ntec')
        assert 'calibrated_ntec' in ntec['criteria']
        assert len(ntec['missing']) == 6


def test_adaptive_validation_policy():
    import tempfile
    from pathlib import Path
    from pipeline.search.adaptive_validation import (
        allocate_validation_batch, experimental_slate, record_validation,
        regional_calibration, priority_adjustment, persist_experimental_slate)
    candidates = [
        ('SAC', 'Fe', 'N4', 'N-graphene', 'OH'),
        ('SAC', 'Co', 'N4', 'N-graphene', 'OH'),
        ('MoltenMetal', 'Ga', 'Sn', 10, 1000),
        ('MoltenMetal', 'Bi', 'Ni', 10, 1000),
        ('MXene', 'Ti', 'C', 2, 'O', 'Fe'),
    ]
    objectives = np.array([[0.2, 1], [0.3, 1], [0.4, 1], [0.5, 1], [0.6, 1]])
    with tempfile.TemporaryDirectory() as tmp:
        db = str(Path(tmp) / 'adaptive.sqlite')
        selected = allocate_validation_batch(candidates, objectives, 3, db, 'test',
                                             min_per_class=1,
                                             uncertainties=[0, 0, 0, 1, 0])
        assert {candidates[i][0] for i in selected} == {'SAC', 'MoltenMetal', 'MXene'}
        # MetalHydride is enumerated but does not consume reserved validation slots.
        hydride_cands = candidates + [('MetalHydride', 'La', 'H2', 'None', 'None', 400)]
        hydride_obj = np.vstack([objectives, [[0.01, 1]]])
        selected_h = allocate_validation_batch(
            hydride_cands, hydride_obj, 3, db, 'test', min_per_class=1)
        assert 'MetalHydride' not in {hydride_cands[i][0] for i in selected_h}
        # Pyrolysis reserved slots skip every phase-unstable class (ADR 0001).
        pyro_cands = [
            ('SAC', 'Fe', 'N4', 'N-graphene', 'OH'),
            ('MOF', 'W', 'Triazolate', 'N2P2', 16.0),
            ('MoltenMetal', 'Bi', 'Ni', 10, 1000),
        ]
        pyro_obj = np.array([[0.9, 1], [0.01, 1], [0.8, 1]])
        selected_p = allocate_validation_batch(
            pyro_cands, pyro_obj, 2, db, 'turquoise_pyrolysis', min_per_class=1)
        assert {pyro_cands[i][0] for i in selected_p} == {'SAC', 'MoltenMetal'}
        try:
            allocate_validation_batch(candidates, objectives, 2, db, 'test', min_per_class=1)
        except ValueError as exc:
            assert 'cannot satisfy class quota' in str(exc)
        else:
            raise AssertionError('undersized class quota budget must fail closed')
        record_validation(db, 'test', candidates[0], 0.2, 1.2, 'dft', False,
                          {'source_id': 'calc:1'})
        stats = regional_calibration(db, 'test')
        region = '|'.join(__import__('pipeline.search.discovery', fromlist=['discovery_region']).discovery_region(candidates[0]))
        assert stats[region]['mae'] == 1.0 and stats[region]['productivity'] == 0.0
        assert priority_adjustment(db, 'test', [candidates[0]]) < 0.5
        # A new region in the same class inherits class-level disagreement
        # until it accumulates its own calibration evidence.
        assert priority_adjustment(db, 'test', [candidates[1]]) < 0.0
        slate = experimental_slate(candidates, objectives, 5)
        assert len(slate) == 5 and len(set(slate)) == 5
        persist_experimental_slate(db, 'test', candidates, objectives, slate)
        import sqlite3
        with sqlite3.connect(db) as conn:
            assert conn.execute("SELECT COUNT(*) FROM experimental_slate").fetchone()[0] == 5


def test_sssp_and_candidate_neb_workflow():
    import tempfile
    from pathlib import Path
    from ase.build import fcc111
    from pipeline.validation.qe_workflows import (verify_sssp, methane_dissociation_images,
                                       write_qe_neb_input)
    encoded_elements = "Ag Al Au B Ba Bi Br C Ca Ce Cl Co Cr Cs Cu Dy Er Eu F Fe Ga Gd Ge Hf I In Ir K La Li Mg Mn Mo N Na Nb Nd Ni O P Pb Pd Pr Pt Rb Re Rh Ru S Sb Sc Se Si Sm Sn Sr Ta Te Ti V W Y Yb Zn Zr".split()
    verified = verify_sssp(encoded_elements)
    assert verified['valid'], verified['errors']
    slab = fcc111('Ni', size=(2, 2, 3), vacuum=10.0)
    images = methane_dissociation_images(slab, active_index=len(slab)-1, n_images=5)
    with tempfile.TemporaryDirectory() as tmp:
        result = write_qe_neb_input(images, str(Path(tmp) / 'ni_ch4.neb.in'), 'ni_ch4')
        text = Path(result['path']).read_text()
        assert "CI_scheme='auto'" in text and 'nspin=2' in text
        assert result['n_images'] == 5

    from pipeline.validation.dft_validator import generate_slab_scf_input
    slab_input = generate_slab_scf_input(
        ['Pd', 'N'], [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)],
        [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 15.0]],
        kpoints=(1, 1, 1),
    )
    assert 'nspin = 2' in slab_input
    assert 'starting_magnetization(1) = 0.05' in slab_input
    assert 'starting_magnetization(2) = 0.05' in slab_input
    assert '1 1 1  0 0 0' in slab_input


def test_orr_multisite_and_corrections():
    from ase.build import fcc111
    from pipeline.validation.orr_workflows import (ORRCorrections, enumerate_surface_sites,
                                        apply_orr_corrections, select_lowest_site)
    slab = fcc111('Pt', size=(2, 2, 3), vacuum=8.0)
    sites = enumerate_surface_sites(slab)
    assert any(x['kind'] == 'atop' for x in sites)
    corrected = apply_orr_corrections(1.0, 2.0, 3.5,
        ORRCorrections(electrode_potential_V=0.8, source_id='protocol:test'))
    assert corrected['dG_OH_eV'] < 1.0
    best = select_lowest_site([{'site_id': 'b', 'converged': True, 'E': -1},
                               {'site_id': 'a', 'converged': True, 'E': -2}], 'E')
    assert best['site_id'] == 'a'


def test_production_qe_workflow_fails_closed_and_resumes():
    import tempfile
    from pathlib import Path
    from pipeline.validation.dft_fuel_cell import converged_energy
    from pipeline.validation.production_workflow import (
        ORR_STAGES, orr_campaign_status, qe_output_status)
    from pipeline.validation.dft_validator import PW_X, generate_molecule_input
    from pipeline.validation.qe_workflows import (
        QEExecutionConfig, build_qe_command)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        partial = root / 'candidate_clean.out'
        partial.write_text('! total energy = -10.000000 Ry\niteration # 2\n')
        assert qe_output_status(partial) == 'incomplete'
        assert converged_energy(partial) is None
        status = orr_campaign_status(root, 'candidate')
        assert not status['complete'] and status['next_stage'] == 'clean'
        complete = ('! total energy = -10.000000 Ry\n'
                    'convergence has been achieved\n'
                    'End of BFGS Geometry Optimization\nJOB DONE\n')
        for stage in ORR_STAGES:
            (root / f'candidate_{stage}.out').write_text(complete)
        status = orr_campaign_status(root, 'candidate')
        assert status['complete'] and status['orr_result_allowed']
        assert converged_energy(root / 'candidate_h2.out') == -10.0

        molecule = generate_molecule_input(
            ['H', 'H'], [(7.5, 7.5, 7.13), (7.5, 7.5, 7.87)],
            15.0, 'h2_test', calculation='scf')
        assert "nspin = 1" in molecule
        assert "occupations = 'fixed'" in molecule
        assert "K_POINTS {gamma}" in molecule
        assert 'smearing' not in molecule
        command = build_qe_command(
            PW_X, str(root / 'h2.in'),
            QEExecutionConfig(mpi_ranks=4, omp_threads=1, kpoint_pools=2))
        assert isinstance(command, list)
        assert command[-4:] == ['-nk', '2', '-in',
                                str((root / 'h2.in').resolve())]
        assert command[-5].endswith('/pw.x')
        try:
            QEExecutionConfig(mpi_ranks=4, image_groups=3).validate(neb=True)
        except ValueError:
            pass
        else:
            raise AssertionError('incompatible MPI/image layout must fail')


def test_validation_task_queue_is_candidate_keyed_and_resumable():
    import tempfile
    from pathlib import Path
    from pipeline.validation.task_queue import ValidationTaskQueue
    from pipeline.search.discovery import candidate_id
    genome = ('SAC', 'Fe', 'N4', 'N-graphene', 'OH')
    with tempfile.TemporaryDirectory() as tmp:
        queue = ValidationTaskQueue(Path(tmp) / 'tasks.sqlite')
        cid = queue.enqueue('turquoise_hydrogen', genome, 'screening_dft', 'v1')
        assert cid == candidate_id(genome)
        assert queue.claim('turquoise_hydrogen', cid, 'screening_dft', 'v1')
        assert not queue.claim('turquoise_hydrogen', cid, 'screening_dft', 'v1')
        queue.finish('turquoise_hydrogen', cid, 'screening_dft', 'v1',
                     True, result_path='result.json')
        assert queue.summary('turquoise_hydrogen', 'screening_dft') == {
            'converged': 1}
        queue.enqueue('turquoise_hydrogen', genome, 'screening_dft', 'v1')
        assert queue.summary('turquoise_hydrogen', 'screening_dft') == {
            'converged': 1}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. ELEMENT EXTRACTORS (4 copies must agree)
# ═══════════════════════════════════════════════════════════════════════════════

def test_element_extractors_consistent():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome
    from pipeline.screening.fc_genetic_optimizer import _extract_elements_from_genome as fc_extract
    from pipeline.screening.genetic_optimizer import _extract_elements_from_genome as methane_extract
    from pipeline.screening.surface_screener import _extract_elements as screener_extract
    from pipeline.screening.fc_screener import _extract_elements as fc_screener_extract
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(20):
            g = generate_random_genome(cls)
            e1 = sorted(fc_extract(g))
            e2 = sorted(methane_extract(g))
            e3 = sorted(screener_extract(g))
            e4 = sorted(fc_screener_extract(g))
            assert e1 == e2 == e3 == e4, (
                f"{cls}: extractors disagree\n"
                f"  FC GA:      {e1}\n"
                f"  Methane GA: {e2}\n"
                f"  Screener:   {e3}\n"
                f"  FC Scr:     {e4}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# 6. STRUCTURE GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

def test_structure_generation_all_classes():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome
    from pipeline.screening.surface_screener import generate_structure
    from scipy.spatial.distance import pdist
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(20):
            g = generate_random_genome(cls)
            atoms, idx, mc = generate_structure(g)
            assert len(atoms) > 0, f"{cls}: empty structure"
            assert atoms.pbc.any(), f"{cls}: PBC not set"
            assert mc == cls, f"{cls}: returned class {mc}"
            assert 'X' not in atoms.get_chemical_symbols(), f"{cls}: dummy atom"
            if len(atoms) > 1:
                assert pdist(atoms.positions).min() > 0.45, (
                    f"{cls}: overlapping atoms")


def test_structure_generation_is_deterministic_and_genome_sensitive():
    import hashlib
    import random
    from pipeline.screening.surface_screener import generate_structure

    def fingerprint(genome):
        atoms, _, _ = generate_structure(genome)
        payload = ('|'.join(atoms.get_chemical_symbols()).encode() +
                   np.asarray(atoms.positions).round(8).tobytes() +
                   np.asarray(atoms.cell).round(8).tobytes())
        return hashlib.sha256(payload).hexdigest()

    candidate = ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0,
                 ('Cu', 'Fe', 'Co'), 3, 1)
    random.seed(1); first = fingerprint(candidate)
    random.seed(999); second = fingerprint(candidate)
    assert first == second

    pairs = [
        (('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, (), 1, 0),
         ('SolidCatalyst', 'Ni', 'Graphene', 'fcc111', 0.0, (), 1, 0)),
        (('SAC', 'Fe', 'N4', 'N-graphene', 'none'),
         ('SAC', 'Fe', 'N4', 'BN_sheet', 'OH')),
        (('MOF', 'Fe', 'BDC', 'N4', 4.0),
         ('MOF', 'Fe', 'TCPP', 'N4', 20.0)),
        (('MAXPhase', 'Ti', 'Al', 'C', 1, 'None', 'basal_0001'),
         ('MAXPhase', 'Ti', 'Al', 'N', 3, 'None', 'edge_1120')),
        (('MXene', 'Ti', 'C', 1, 'O', 'None'),
         ('MXene', 'Ti', 'N', 3, 'F', 'None')),
        (('MetalFreeCarbon', 'pyridinic', 0.02, 'none', 'graphene', 'B'),
         ('MetalFreeCarbon', 'graphitic', 0.20, 'divacancy', 'CNT', 'P')),
    ]
    for left, right in pairs:
        assert fingerprint(left) != fingerprint(right)


# ═══════════════════════════════════════════════════════════════════════════════
# 7. PEMFC MODEL
# ═══════════════════════════════════════════════════════════════════════════════

def test_tafel_covers_all_classes():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.process.pemfc_model import TAFEL_SLOPE_BY_CLASS
    for cls in ALL_MATERIAL_CLASSES:
        assert cls in TAFEL_SLOPE_BY_CLASS, f"Missing Tafel slope for {cls}"


def test_pemfc_power_monotonic_with_eta():
    from pipeline.process.pemfc_model import simulate_pemfc, PEMFCConfig
    etas = [0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80]
    results = [simulate_pemfc(PEMFCConfig(orr_overpotential_V=e)) for e in etas]
    powers = [r['peak_power_W_cm2'] for r in results]
    for i in range(len(powers) - 1):
        assert powers[i] >= powers[i+1], f"Power not monotonic: eta={etas[i]}->{etas[i+1]}"


def test_tafel_slope_affects_power():
    from pipeline.process.pemfc_model import simulate_pemfc, PEMFCConfig
    r_low = simulate_pemfc(PEMFCConfig(orr_overpotential_V=0.35, orr_tafel_slope_mV_dec=65))
    r_high = simulate_pemfc(PEMFCConfig(orr_overpotential_V=0.35, orr_tafel_slope_mV_dec=130))
    assert r_low['peak_power_W_cm2'] > r_high['peak_power_W_cm2'], \
        f"Lower Tafel should give higher power: {r_low['peak_power_W_cm2']} vs {r_high['peak_power_W_cm2']}"


def test_pemfc_efficiency_in_range():
    from pipeline.process.pemfc_model import simulate_pemfc, PEMFCConfig
    r = simulate_pemfc(PEMFCConfig(orr_overpotential_V=0.35))
    eff = r['efficiency_at_peak']
    assert 0.10 < eff < 0.60, f"Efficiency {eff} out of physical range"


# ═══════════════════════════════════════════════════════════════════════════════
# 8. STACK MODEL
# ═══════════════════════════════════════════════════════════════════════════════

def test_stack_model():
    from pipeline.process.fuel_cell_stack import model_stack, StackConfig
    stack = model_stack(StackConfig(cell_voltage_V=0.65, current_density_A_cm2=1.5))
    assert stack['net_power_kW'] > 0, f"Negative net power: {stack['net_power_kW']}"
    assert 0 < stack['system_efficiency'] < 1, f"Efficiency out of range: {stack['system_efficiency']}"
    assert stack['cost_per_kW'] > 0, f"Negative cost: {stack['cost_per_kW']}"


def test_tea_is_scenario_labelled_and_unclamped():
    from pipeline.process.tea import estimate_scenario_range
    low_conversion = estimate_scenario_range(0.10)
    high_conversion = estimate_scenario_range(0.90)
    assert low_conversion['max_usd_kg'] > low_conversion['min_usd_kg']
    assert (low_conversion['estimates']['base']['h2_cost_usd_kg'] >
            high_conversion['estimates']['base']['h2_cost_usd_kg'])
    assert all(
        value['evidence_level'] == 'screening_scenario_not_measured_tea'
        for value in low_conversion['estimates'].values())


# ═══════════════════════════════════════════════════════════════════════════════
# 9. CHE / UTILS
# ═══════════════════════════════════════════════════════════════════════════════

def test_orr_overpotential_ideal():
    from pipeline.common.utils import orr_overpotential
    eta, rds = orr_overpotential(1.23, 2.46, 3.69)
    assert abs(eta) < 0.001, f"Ideal overpotential should be ~0, got {eta}"


def test_abundance_cost_penalty():
    from pipeline.common.utils import abundance_cost_penalty
    assert abundance_cost_penalty(['Fe']) == 0.0, "Fe should be zero cost"
    assert abundance_cost_penalty(['Ir']) == -2.0, "Ir should be max penalty"
    assert abundance_cost_penalty(['Fe', 'Ir']) < 0, "Geo mean should catch Ir"
    assert abundance_cost_penalty([]) == 0.0, "Empty should be zero"


# ═══════════════════════════════════════════════════════════════════════════════
# 10. REPORT GENERATOR
# ═══════════════════════════════════════════════════════════════════════════════

def test_report_campaign_path():
    from pipeline.evidence.report_generator import generate_full_report
    path = generate_full_report({
        'phase5_ga': {'total_evaluated': 1000},
        'phase5_stack': {'best_power_W_cm2': 1.35, 'best_efficiency': 0.31, 'min_overpotential_V': 0.32},
    })
    with open(path) as f:
        content = f.read()
    assert 'N/A' not in content.split('FC catalysts')[1].split('\n')[0], "FC row shows N/A in campaign mode"


def test_report_orchestrator_path():
    from pipeline.evidence.report_generator import generate_full_report
    path = generate_full_report({
        'phase5': {'n_cathodes_screened': 137, 'best_power_W_cm2': 1.2,
                   'best_efficiency': 0.28, 'min_overpotential_V': 0.35},
    })
    with open(path) as f:
        content = f.read()
    assert 'N/A' not in content.split('PEMFC power')[1].split('\n')[0], "PEMFC row shows N/A in orchestrator mode"


def test_report_empty_no_crash():
    from pipeline.evidence.report_generator import generate_full_report
    path = generate_full_report({})
    assert os.path.exists(path)


# ═══════════════════════════════════════════════════════════════════════════════
# 11. CROSSOVER & MUTATION
# ═══════════════════════════════════════════════════════════════════════════════

def test_crossover_preserves_class():
    from pipeline.common.catalyst_spaces import generate_random_genome, crossover
    for _ in range(100):
        p1 = generate_random_genome('SAC')
        p2 = generate_random_genome('SAC')
        child = crossover(p1, p2)
        assert child[0] == 'SAC', f"Crossover changed class to {child[0]}"


def test_mutation_preserves_class():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome, mutate
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(20):
            g = generate_random_genome(cls)
            m = mutate(g, rate=1.0)  # force mutation
            assert m[0] == cls, f"Mutation changed {cls} to {m[0]}"


# ═══════════════════════════════════════════════════════════════════════════════
# 12. OOD CONFIDENCE
# ═══════════════════════════════════════════════════════════════════════════════

def test_ood_high_confidence_metals():
    from pipeline.common.catalyst_spaces import generate_random_genome
    from pipeline.common.ood_detector import compute_model_confidence
    from pipeline.screening.fc_genetic_optimizer import _extract_elements_from_genome
    # Metal slabs should have high confidence (>0.7)
    for cls in ['SolidCatalyst', 'HEA', 'SAA']:
        for _ in range(10):
            g = generate_random_genome(cls)
            elements = _extract_elements_from_genome(g)
            conf = compute_model_confidence(g, elements)
            assert conf > 0.6, f"{cls}: confidence {conf:.2f} too low (expected >0.6)"


def test_ood_low_confidence_ood():
    from pipeline.common.catalyst_spaces import generate_random_genome
    from pipeline.common.ood_detector import compute_model_confidence
    from pipeline.screening.fc_genetic_optimizer import _extract_elements_from_genome
    # OOD classes should have low confidence (<0.5)
    for cls in ['MOF', 'COF', 'MetalFreeCarbon']:
        for _ in range(10):
            g = generate_random_genome(cls)
            elements = _extract_elements_from_genome(g)
            conf = compute_model_confidence(g, elements)
            assert conf < 0.5, f"{cls}: confidence {conf:.2f} too high (expected <0.5)"


def test_ood_penalty_scales_objectives():
    from pipeline.common.ood_detector import confidence_penalty
    # High confidence → penalty near 0.0 (no shift)
    assert abs(confidence_penalty(1.0) - 0.0) < 0.01, "conf=1.0 should give penalty=0.0"
    # Low confidence → penalty > 0.5 (significant shift)
    assert confidence_penalty(0.2) > 0.5, "conf=0.2 should give penalty>0.5"
    # Zero confidence → maximum penalty = 1.0
    assert abs(confidence_penalty(0.0) - 1.0) < 0.01, "conf=0.0 should give penalty=1.0"


def test_ood_nsga2_integration():
    from pipeline.common.catalyst_spaces import generate_random_genome, FEATURE_DIM
    from pipeline.screening.fc_genetic_optimizer import (
        ORRCatalystSurrogate, compute_orr_objectives_surrogate
    )
    model = ORRCatalystSurrogate(input_dim=FEATURE_DIM); model.eval()
    # Generate OOD and in-distribution populations
    in_dist = [generate_random_genome('SolidCatalyst') for _ in range(50)]
    ood = [generate_random_genome('MetalFreeCarbon') for _ in range(50)]
    obj_in = compute_orr_objectives_surrogate(in_dist, model, 'cpu')
    obj_ood = compute_orr_objectives_surrogate(ood, model, 'cpu')
    # OOD overpotentials (obj[:, 0]) should be inflated by penalty
    mean_in = obj_in[:, 0].mean()
    mean_ood = obj_ood[:, 0].mean()
    # OOD should have higher (worse) mean overpotential after penalty
    # (MetalFreeCarbon conf ~0.15 → penalty ~2.7×)
    assert mean_ood > mean_in, (
        f"OOD penalty not working: mean_in={mean_in:.3f}, mean_ood={mean_ood:.3f}"
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 13. EXHAUSTIVE COVERAGE — NO CLASS GETS PRUNED
# ═══════════════════════════════════════════════════════════════════════════════

def test_all_elements_in_abundance_table():
    from pipeline.common.catalyst_spaces import generate_random_genome
    from pipeline.screening.fc_genetic_optimizer import _extract_elements_from_genome
    from pipeline.common.utils import CRUSTAL_ABUNDANCE_PPM
    missing = set()
    for _ in range(3000):
        g = generate_random_genome()
        for e in _extract_elements_from_genome(g):
            if e not in CRUSTAL_ABUNDANCE_PPM:
                missing.add(e)
    assert not missing, f"Elements missing from CRUSTAL_ABUNDANCE: {sorted(missing)}"


def test_all_elements_in_price_table():
    from pipeline.common.catalyst_spaces import generate_random_genome
    from pipeline.screening.fc_genetic_optimizer import _extract_elements_from_genome
    from pipeline.common.utils import METAL_PRICE_USD_KG
    missing = set()
    for _ in range(3000):
        g = generate_random_genome()
        for e in _extract_elements_from_genome(g):
            if e not in METAL_PRICE_USD_KG:
                missing.add(e)
    assert not missing, f"Elements missing from METAL_PRICE_USD_KG: {sorted(missing)}"


def test_all_classes_viable_both_applications():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.common.utils import VALID_CLASSES_PYROLYSIS, VALID_CLASSES_FUEL_CELL
    for cls in ALL_MATERIAL_CLASSES:
        assert cls in VALID_CLASSES_PYROLYSIS, f"{cls} excluded from pyrolysis"
        assert cls in VALID_CLASSES_FUEL_CELL, f"{cls} excluded from fuel cell"


def test_bep_params_all_classes():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.common.utils import bep_activation_energy
    for cls in ALL_MATERIAL_CLASSES:
        # Should not fall back to default — each class must have specific params
        e1 = bep_activation_energy(0.5, material_class=cls)
        e_default = bep_activation_energy(0.5)  # no class = default
        # At least some classes should differ from default
        assert isinstance(e1, float) and e1 > 0


def test_ood_confidence_all_classes():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.common.ood_detector import CLASS_CONFIDENCE
    for cls in ALL_MATERIAL_CLASSES:
        assert cls in CLASS_CONFIDENCE, f"OOD CLASS_CONFIDENCE missing: {cls}"


def test_tafel_all_classes():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.process.pemfc_model import TAFEL_SLOPE_BY_CLASS
    for cls in ALL_MATERIAL_CLASSES:
        assert cls in TAFEL_SLOPE_BY_CLASS, f"TAFEL_SLOPE missing: {cls}"


def test_pyrolysis_mode_coking_bonus():
    import os
    from pipeline.screening.genetic_optimizer import compute_objectives_surrogate, GAConfig
    from pipeline.screening.surrogate_model import CatalystSurrogate

    # Mock surrogate and population
    model = CatalystSurrogate()
    # We want a genome with Ga/In/Sn/Bi (e.g. MoltenMetal with Sb/Ga) and one without (e.g. SolidCatalyst with Fe)
    pop = [
        ('MoltenMetal', 'Ga', 'None', 0.0, 1000),  # low-melting metal Ga
        ('SolidCatalyst', 'Fe', 'None', 'fcc111', 0.0, ['Fe'], 0, 0),  # non-liquid Fe
    ]

    from unittest.mock import patch
    with patch('pipeline.screening.genetic_optimizer.predict_batch') as mock_predict:
        mock_predict.return_value = {
            'valid_prob': np.array([1.0, 1.0]),
            'E_act': np.array([0.5, 0.6]),
            'coking_index': np.array([1.0, 2.0]),
            'segregation_energy': np.array([-0.1, -0.2]),
        }

        # NTEC without measured operating evidence must not create a bonus.
        os.environ['PYROLYSIS_MODE'] = 'ntec'
        os.environ.pop('NTEC_CONDITIONS_JSON', None)
        objs_unknown = compute_objectives_surrogate(pop, model, device='cpu')

        # Operating evidence alone is insufficient without a paired control.
        import json
        os.environ['NTEC_CONDITIONS_JSON'] = json.dumps({
            'shear_rate_s': 1e4, 'interfacial_field_V_m': 1e8,
            'mechanical_power_W_kg': 1e3, 'carbon_detachment_fraction': 1.0,
            'field_measurement_source': 'test:field-measurement',
        })
        objs_uncalibrated = compute_objectives_surrogate(pop, model, device='cpu')

        # Explicit paired NTEC/control measurements activate the bounded transfer.
        os.environ['NTEC_CONDITIONS_JSON'] = json.dumps({
            'shear_rate_s': 1e4, 'interfacial_field_V_m': 1e8,
            'mechanical_power_W_kg': 1e3, 'carbon_detachment_fraction': 1.0,
            'field_measurement_source': 'test:field-measurement',
            'paired_control_source': 'test:paired-control',
            'paired_control_count': 2,
            'measured_barrier_reduction_eV': 0.25,
            'measured_coking_delta_eV': 3.0,
        })
        objs_ntec = compute_objectives_surrogate(pop, model, device='cpu')

        # Test under thermocatalytic mode
        os.environ['PYROLYSIS_MODE'] = 'thermocatalytic'
        objs_thermo = compute_objectives_surrogate(pop, model, device='cpu')

    # Reset environment
    if 'PYROLYSIS_MODE' in os.environ:
        del os.environ['PYROLYSIS_MODE']
    os.environ.pop('NTEC_CONDITIONS_JSON', None)

    # For Ga molten metal candidate, NTEC coking index objective (index 1) should be lower (more negative = better coking resistance)
    # Since obj2 = -(coking_index + bonus), objs_ntec[0, 1] = objs_thermo[0, 1] - 3.0
    diff_ga = objs_ntec[0, 1] - objs_thermo[0, 1]
    assert np.isclose(diff_ga, -3.0), f"Liquid metal Ga coking bonus not applied correctly, got diff: {diff_ga}"
    assert np.isclose(objs_unknown[0, 1], objs_thermo[0, 1]), \
        "NTEC without measured inputs must receive zero bonus"
    assert np.isclose(objs_uncalibrated[0, 1], objs_thermo[0, 1]), \
        "NTEC operating inputs without a paired control must receive zero bonus"

    # For Fe catalyst, there should be no bonus, so diff should be 0.0
    diff_fe = objs_ntec[1, 1] - objs_thermo[1, 1]
    assert np.isclose(diff_fe, 0.0), f"Non-liquid metal Fe coking bonus applied incorrectly, got diff: {diff_fe}"


def test_cathode_sac_genome_5tuple():
    from pipeline.screening.fc_cathode_screener import generate_fc_catalyst_list
    candidates = generate_fc_catalyst_list()
    sacs = [c for c in candidates if c['type'] == 'SAC']
    assert len(sacs) > 0, "No SAC candidates generated"
    for c in sacs:
        assert len(c['genome']) == 5, f"SAC genome should be 5-tuple, got {len(c['genome'])}: {c['genome']}"


def test_deterministic_hierarchical_pool():
    from pipeline.common.catalyst_spaces import generate_hierarchical_htvs_pool
    # Test fallback behavior when model is None
    pool = generate_hierarchical_htvs_pool(pool_size=100, scorer=None)
    assert len(pool) == 100, f"Expected pool size 100, got {len(pool)}"
    # Check that genomes are generated and valid
    for g in pool:
        assert isinstance(g, tuple), "Genome must be a tuple"
        assert len(g) >= 2, "Genome must have at least class name and one parameter"


def test_hierarchical_rounds_cover_complementary_cells():
    from pipeline.common.catalyst_spaces import generate_hierarchical_htvs_pool
    from pipeline.search.discovery import candidate_id
    first = generate_hierarchical_htvs_pool(500, campaign_round=0)
    second = generate_hierarchical_htvs_pool(500, campaign_round=1)
    ids_first = {candidate_id(g) for g in first}
    ids_second = {candidate_id(g) for g in second}
    assert ids_first != ids_second, "Campaign rounds must not regenerate the same fixed-stride pool"
    assert len(ids_first) == len(first), "Round 0 contains duplicate canonical candidates"
    assert len(ids_second) == len(second), "Round 1 contains duplicate canonical candidates"


def test_discovery_batch_prioritizes_unseen_regions():
    from pipeline.search.discovery import discovery_region, select_discovery_batch
    candidates = [
        ('SAC', 'Fe', 'N4', 'N-graphene', 'none'),
        ('SAC', 'Co', 'N4', 'N-graphene', 'none'),
        ('SAC', 'Ni', 'N3C', 'N-CNT', 'OH'),
        ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000),
    ]
    objectives = np.array([[0.10, 0, 0, 0], [0.11, 0, 0, 0],
                           [0.30, 0, 0, 0], [0.25, 0, 0, 0]])
    selected = select_discovery_batch(candidates, objectives, 3, evaluated=[candidates[0]])
    regions = {discovery_region(candidates[i]) for i in selected}
    assert len(selected) == 3
    assert len(regions) == 3, "Discovery acquisition collapsed into an already-covered chemistry region"


def test_candidate_ids_are_canonical():
    from pipeline.search.discovery import candidate_id
    a = ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.01, ('B', 'N'), 2, 0)
    b = ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0100000001, ('N', 'B'), 2, 0)
    assert candidate_id(a) == candidate_id(b), "Equivalent dopant permutations need one candidate ID"

    molten_a = ('MoltenMetal', 'Ga', 'Fe', 0.0, 1000)
    molten_b = ('MoltenMetal', 'Ga', 'None', 0.0, 1000)
    assert candidate_id(molten_a) == candidate_id(molten_b)

    perovskite_a = ('Perovskite', 'La', 'Fe', 'Co', 0.0, 'none')
    perovskite_b = ('Perovskite', 'La', 'Fe', 'None', 0.0, 'none')
    assert candidate_id(perovskite_a) == candidate_id(perovskite_b)


def test_design_space_audit_preserves_all_sizable_classes():
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.common.design_space_provenance import validate_provenance
    from pipeline.evidence.design_space_audit import audit_design_space
    from pipeline.search.indexed_space import is_physically_admissible

    provenance = validate_provenance(ALL_MATERIAL_CLASSES)
    assert provenance['valid'] and provenance['classes'] == 14
    report = audit_design_space(sample_per_class=128)
    assert report['valid'], report['failures']
    assert report['classes_represented'] == 14
    assert set(report['classes']) == set(ALL_MATERIAL_CLASSES)
    assert report['minimum_raw_per_class'] >= 1000
    assert report['minimum_projected_admissible_per_class'] >= 1000
    assert report['canonical_total'] < report['raw_cartesian_total']

    assert not is_physically_admissible(
        ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0,
         ('N', 'B'), 2, 0))[0]
    assert not is_physically_admissible(
        ('MoltenMetal', 'Ga', 'Fe', 0.0, 1000))[0]


def test_discovery_metadata_is_persistable():
    import pandas as pd
    from pipeline.search.discovery import add_discovery_metadata
    genome = ('SAC', 'Fe', 'N4', 'N-graphene', 'OH')
    out = add_discovery_metadata(pd.DataFrame({'genome': [str(genome)]}))
    assert out.loc[0, 'candidate_id']
    assert out.loc[0, 'discovery_region'].startswith('SAC|')


def test_indexed_space_boundaries_and_classes():
    from pipeline.search.indexed_space import (CLASS_OFFSETS, CLASS_ORDER, CLASS_SIZES,
                                        TOTAL_SIZE, candidate_at, candidate_at_class)
    from pipeline.common.catalyst_spaces import estimate_design_space_size
    assert TOTAL_SIZE == estimate_design_space_size()['TOTAL']
    for cls in CLASS_ORDER:
        assert candidate_at(CLASS_OFFSETS[cls])[0] == cls
        assert candidate_at_class(cls, CLASS_SIZES[cls] - 1)[0] == cls


def test_indexed_worker_shards_are_disjoint():
    from pipeline.search.indexed_space import iter_shard
    a = {i for i, _ in iter_shard(0, 101, 0, 3)}
    b = {i for i, _ in iter_shard(0, 101, 1, 3)}
    c = {i for i, _ in iter_shard(0, 101, 2, 3)}
    assert not (a & b or a & c or b & c)
    assert a | b | c == set(range(101))


def test_streaming_scan_resumes_without_rescoring():
    import tempfile
    from pathlib import Path
    from pipeline.search.exhaustive_search import ScanConfig, run_streaming_scan
    calls = []
    def scorer(genomes):
        calls.append(len(genomes))
        return np.column_stack([np.arange(len(genomes), dtype=float),
                                np.zeros((len(genomes), 3))])
    with tempfile.TemporaryDirectory() as tmp:
        db = str(Path(tmp) / 'scan.sqlite')
        first = run_streaming_scan(ScanConfig('test', db, stop=40, batch_size=10, max_batches=2), scorer)
        second = run_streaming_scan(ScanConfig('test', db, stop=40, batch_size=10), scorer)
        assert first['processed_this_run'] == 20
        assert second['processed_this_run'] == 20
        assert second['complete']
        assert sum(calls) <= 40  # invalid candidates may be rejected before scoring


def test_sharded_scan_matches_serial_and_fails_closed():
    import sqlite3
    import tempfile
    from pathlib import Path
    from pipeline.search.exhaustive_search import (
        ScanConfig, run_sharded_scan, run_streaming_scan)
    from pipeline.common.catalyst_spaces import encode_population

    def scorer(genomes):
        encoded = encode_population(genomes)
        return np.column_stack([
            encoded[:, 0] + 0.01 * encoded[:, 1],
            encoded[:, 2], encoded[:, 3], encoded[:, 4]])

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        serial_db = root / 'serial.sqlite'
        sharded_db = root / 'sharded.sqlite'
        serial = run_streaming_scan(ScanConfig(
            'shard_equivalence', str(serial_db), stop=2000,
            batch_size=128, global_archive_size=100,
            state_id='equivalence'), scorer)
        sharded = run_sharded_scan(ScanConfig(
            'shard_equivalence', str(sharded_db), stop=2000,
            batch_size=128, global_archive_size=100,
            state_id='equivalence'), scorer, workers=4)
        assert serial['complete'] and sharded['complete']
        tables = ('region_champions', 'global_archive', 'objective_archive',
                  'regional_objective_champions')
        with sqlite3.connect(serial_db) as left, sqlite3.connect(sharded_db) as right:
            for table in tables:
                left_rows = left.execute(f'SELECT * FROM {table} ORDER BY 1,2,3').fetchall()
                right_rows = right.execute(f'SELECT * FROM {table} ORDER BY 1,2,3').fetchall()
                assert left_rows == right_rows, f'shard merge changed {table}'
            progress = right.execute(
                "SELECT next_index, processed FROM scan_progress "
                "WHERE application=? AND state_id=?",
                ('shard_equivalence', 'equivalence')).fetchone()
            assert progress == (2000, 2000)

        incomplete_db = root / 'incomplete.sqlite'
        partial = run_sharded_scan(ScanConfig(
            'shard_resume', str(incomplete_db), stop=400,
            batch_size=10, max_batches=1, state_id='resume'),
            scorer, workers=4)
        assert not partial['complete']
        assert not incomplete_db.exists(), 'partial shards must not create aggregate evidence'
        resumed = run_sharded_scan(ScanConfig(
            'shard_resume', str(incomplete_db), stop=400,
            batch_size=10, state_id='resume'), scorer, workers=4)
        assert resumed['complete']
        with sqlite3.connect(incomplete_db) as conn:
            progress = conn.execute(
                "SELECT next_index, processed FROM scan_progress "
                "WHERE application='shard_resume' AND state_id='resume'").fetchone()
        assert progress == (400, 400)

        tiny = run_sharded_scan(ScanConfig(
            'tiny_range', str(root / 'tiny.sqlite'), stop=3,
            batch_size=10, state_id='tiny'), scorer, workers=8)
        assert tiny['complete'] and tiny['workers'] == 3


def test_branch_search_resolves_without_surrogate_pruning():
    import tempfile
    from pathlib import Path
    from pipeline.search.branch_search import BranchConfig, run_branch_and_bound
    from pipeline.search.indexed_space import CLASS_SIZES
    def deliberately_bad_scorer(genomes):
        # A poor probe is not permission to remove its branch.
        return np.column_stack([np.full(len(genomes), 99.0),
                                np.zeros((len(genomes), 3))])
    with tempfile.TemporaryDirectory() as tmp:
        result = run_branch_and_bound(BranchConfig(
            application='branch_test', database=str(Path(tmp) / 'branch.sqlite'),
            leaf_size=CLASS_SIZES['MetalFreeCarbon'] + 1,
            scan_batch_size=512, max_leaves=1, scan_workers=2,
            material_classes=('MetalFreeCarbon',),
        ), deliberately_bad_scorer)
        assert result['complete']
        assert result['node_status_counts'].get('scanned') == 1
        assert result['node_status_counts'].get('pruned', 0) == 0


def test_branch_probes_are_deterministic_low_discrepancy():
    from pipeline.search.branch_search import _probe_indices
    first = _probe_indices(100, 10100, 9)
    second = _probe_indices(100, 10100, 9)
    assert first == second and len(first) == 9
    assert first[0] == 100 and first[-1] == 10099
    gaps = np.diff(first)
    assert len(set(gaps.tolist())) > 1, 'probes regressed to alias-prone uniform spacing'


def test_branch_finite_budget_preserves_class_floor():
    import tempfile
    from pathlib import Path
    from pipeline.search.branch_search import BranchConfig, run_branch_and_bound
    def biased_scorer(genomes):
        primary = np.array([0.0 if genome[0] == 'SAA' else 100.0 for genome in genomes])
        return np.column_stack([primary, np.zeros((len(genomes), 3))])
    with tempfile.TemporaryDirectory() as tmp:
        result = run_branch_and_bound(BranchConfig(
            application='finite_budget_classes', database=str(Path(tmp) / 'branch.sqlite'),
            leaf_size=1000, scan_batch_size=512, max_leaves=2,
            material_classes=('SAA', 'MXene'), min_resolved_leaves_per_class=1,
        ), biased_scorer)
        resolved = result['resolved_terminal_nodes_by_class']
        assert resolved['SAA'] >= 1 and resolved['MXene'] >= 1
        assert result['scheduling_decisions']['class_floor'] > 0


def test_branch_certificate_detects_incomplete_and_gaps():
    import sqlite3
    import tempfile
    from pathlib import Path
    from pipeline.search.branch_search import BranchConfig, run_branch_and_bound, verify_branch_coverage
    def scorer(genomes):
        return np.zeros((len(genomes), 4))
    with tempfile.TemporaryDirectory() as tmp:
        db = str(Path(tmp) / 'certificate.sqlite')
        result = run_branch_and_bound(BranchConfig(
            application='certificate_test', database=db, leaf_size=1000,
            scan_batch_size=256, max_leaves=1,
            material_classes=('MetalFreeCarbon',),
        ), scorer)
        cert = result['coverage_certificate']
        assert cert['gap_free'] and cert['overlap_free']
        assert not cert['complete'] and cert['unresolved_terminal_nodes'] > 0

        conn = sqlite3.connect(db)
        row = conn.execute("SELECT node_id FROM branch_nodes WHERE application=? "
                           "AND status!='expanded' LIMIT 1", ('certificate_test',)).fetchone()
        conn.execute("DELETE FROM branch_nodes WHERE application=? AND node_id=?",
                     ('certificate_test', row[0]))
        conn.commit(); conn.close()
        broken = verify_branch_coverage(db, 'certificate_test', ('MetalFreeCarbon',))
        assert not broken['complete']
        assert broken['errors'], "Deleted terminal interval must invalidate certificate"


def test_branch_rejects_population_mismatch():
    import tempfile
    from pathlib import Path
    from pipeline.search.branch_search import BranchConfig, run_branch_and_bound
    def scorer(genomes):
        return np.zeros((len(genomes), 4))
    with tempfile.TemporaryDirectory() as tmp:
        try:
            run_branch_and_bound(BranchConfig(
                application='mismatch', database=str(Path(tmp) / 'mismatch.sqlite'),
                material_classes=('SAA',), expected_population=25_300_000_000,
            ), scorer)
        except ValueError as exc:
            assert 'denominator mismatch' in str(exc)
        else:
            raise AssertionError('Population mismatch was silently accepted')


def test_tree_calibration_probes_cover_all_classes_deterministically():
    from pipeline.search.indexed_space import deterministic_tree_probes
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES
    first = deterministic_tree_probes(100)
    second = deterministic_tree_probes(100)
    assert first == second
    assert {g[0] for g in first} == set(ALL_MATERIAL_CLASSES)


def test_production_has_only_branch_candidate_search():
    from pathlib import Path
    source = (REPO_ROOT / 'run_production_campaign.py').read_text()
    forbidden = [
        'run_genetic_algorithm', 'run_fc_genetic_algorithm',
        '--exhaustive-scan', '--branch-search', '--pop', '--gens',
        'generate_population', 'generate_random_genome',
    ]
    present = [token for token in forbidden if token in source]
    assert not present, f"Legacy candidate-search paths remain in production: {present}"
    assert 'run_branch_discovery' in source
    assert 'run_fc_branch_discovery' in source
    assert 'QE executables are resolved at execution time' in source
    resolver = (REPO_ROOT / 'pipeline/common/executables.py').read_text()
    assert "env_var=variables.get(name), conda_env='qe-env'" in resolver
    assert ('/' + 'home/') not in source + resolver
    assert "'conda', 'run', '-n', 'quantum-env'" in source
    assert "result.get('mock')" in source


def test_readme_matches_branch_only_contract():
    from pathlib import Path
    readme = (REPO_ROOT / 'README.md').read_text()
    assert '21,092,645,031' in readme
    assert 'Deterministic Branch-and-Bound Discovery' in readme
    assert '--calibration-probes' in readme
    assert '--branch-leaf-size' in readme
    forbidden = ['--pop', '--gens', '--exhaustive-scan', '--branch-search',
                 '25.3-billion-configuration', '21.3-billion-configuration']
    present = [token for token in forbidden if token in readme]
    assert not present, f"README advertises retired search controls: {present}"


def test_readme_contains_no_machine_specific_paths():
    from pathlib import Path
    readme = (REPO_ROOT / 'README.md').read_text()
    forbidden = [
        '/' + 'home/', '/' + 'Users/', 'mini' + 'conda3',
        'ana' + 'conda3', '/' + 'opt/conda', '.gem' + 'ini/antigravity',
    ]
    present = [token for token in forbidden if token in readme]
    assert not present, f'Machine-specific paths remain in README: {present}'


def test_root_documentation_is_canonical():
    from pathlib import Path
    root = REPO_ROOT
    root_docs = {
        path.name for pattern in ('*.md', '*.rst') for path in root.glob(pattern)
    }
    assert root_docs <= {'README.md', 'CHANGELOG.md'}, root_docs
    assert (root / 'docs/FUEL_CELL.md').is_file()
    assert (root / 'docs/TURQUOISE_HYDROGEN.md').is_file()


def test_retired_ga_entry_points_are_blocked():
    from pipeline.screening.genetic_optimizer import run_genetic_algorithm
    from pipeline.screening.fc_genetic_optimizer import run_fc_genetic_algorithm, FCGAConfig
    for fn, args in ((run_genetic_algorithm, ()),
                     (run_fc_genetic_algorithm, (FCGAConfig(),))):
        try:
            fn(*args)
        except RuntimeError as exc:
            assert 'retired' in str(exc)
        else:
            raise AssertionError(f"{fn.__name__} still permits legacy search")


def test_industrial_viability_gates_fail_closed():
    from pipeline.validation.viability import (
        evaluate_turquoise, evaluate_fuel_cell, TurquoiseHydrogenBounds)
    assert evaluate_turquoise({})['status'] == 'unknown'
    good_h2 = evaluate_turquoise({
        'temperature_K': 1000, 'H2_selectivity': 0.98, 'CH4_conversion': 0.8,
        'deactivation_fraction_per_h': 0.005, 'coke_fraction': 0.02,
        'net_energy_kWh_kg_h2': 12.0, 'measured_reactor': 1})
    assert good_h2['status'] == 'pass'
    # H2 selectivity gate is off by default (lumped metric is not true selectivity).
    assert evaluate_turquoise({'H2_selectivity': 0.8, 'CH4_conversion': 0.8,
                               'temperature_K': 1000,
                               'deactivation_fraction_per_h': 0.005,
                               'coke_fraction': 0.02,
                               'net_energy_kWh_kg_h2': 12.0,
                               'measured_reactor': 1})['status'] == 'pass'
    gated = evaluate_turquoise(
        {'H2_selectivity': 0.8, 'CH4_conversion': 0.8, 'temperature_K': 1000,
         'deactivation_fraction_per_h': 0.005, 'coke_fraction': 0.02,
         'net_energy_kWh_kg_h2': 12.0, 'measured_reactor': 1},
        TurquoiseHydrogenBounds(enforce_h2_selectivity_gate=True))
    assert gated['status'] == 'fail'
    assert evaluate_turquoise({
        'temperature_K': 1000, 'CH4_conversion': 0.8,
        'deactivation_fraction_per_h': 0.005, 'coke_fraction': 0.02,
        'net_energy_kWh_kg_h2': 12.0, 'measured_reactor': 1,
        'co2_permitted': True})['status'] == 'fail'
    good_fc = evaluate_fuel_cell({
        'orr_overpotential_V': 0.3, 'peak_power_W_cm2': 1.2,
        'system_efficiency': 0.5, 'voltage_degradation_uV_h': 5,
        'measured_hours': 500, 'measured_mea': 1})
    assert good_fc['status'] == 'pass'
    assert evaluate_fuel_cell({'orr_overpotential_V': 0.6})['status'] == 'fail'


def test_coking_loss_masks_nan_targets():
    """NaN coking targets must not train the coking head or poison other heads."""
    import torch
    import torch.nn as nn
    from pipeline.common.catalyst_spaces import FEATURE_DIM
    from pipeline.screening.surrogate_model import _masked_mse, train_surrogate

    mse = nn.MSELoss()
    pred = torch.tensor([[1.0], [2.0], [3.0]])
    target = torch.tensor([[1.0], [float('nan')], [3.0]])
    valid = torch.tensor([True, True, True])
    loss = _masked_mse(pred, target, valid, mse)
    assert torch.isfinite(loss)
    assert torch.isclose(loss, torch.tensor(0.0))

    rng = np.random.default_rng(0)
    n = 16
    X = rng.standard_normal((n, FEATURE_DIM)).astype(np.float32)
    y_valid = np.ones(n, dtype=np.float32)
    y_de = rng.standard_normal(n).astype(np.float32)
    y_coking = rng.standard_normal(n).astype(np.float32)
    y_coking[:6] = np.nan
    y_seg = rng.standard_normal(n).astype(np.float32)
    y_e = np.abs(rng.standard_normal(n)).astype(np.float32) + 0.2
    model = train_surrogate(
        X, y_valid, y_de, y_coking, y_seg, y_e,
        epochs=2, batch_size=8, device='cpu')
    for p in model.parameters():
        assert torch.isfinite(p).all()


def test_slab_coking_scope_excludes_molten_metal():
    from pipeline.common.application_scope import slab_coking_index_scope
    assert slab_coking_index_scope(('MoltenMetal', 'Bi', 'Ni', 10.0, 1000))['status'] == 'out_of_scope'
    assert slab_coking_index_scope(('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, ('Fe',), 1, 0))['status'] == 'candidate'


def test_phase_stable_at_application_t_per_class():
    from pipeline.common.application_scope import (
        phase_stable_at_application_T, is_turquoise_pyrolysis_candidate,
        VALIDATION_QUOTA_EXEMPT_CLASSES, validation_quota_class_count)
    from pipeline.common.catalyst_spaces import ALL_MATERIAL_CLASSES
    assert phase_stable_at_application_T(('MetalHydride', 'La'))['status'] == 'out_of_scope'
    assert phase_stable_at_application_T(('MOF', 'W', 'Triazolate', 'N2P2', 16.0))['status'] == 'out_of_scope'
    assert phase_stable_at_application_T(('COF', 'W', 'Imide', 'P4', 40.0))['status'] == 'out_of_scope'
    assert phase_stable_at_application_T(('MXene', 'Ti', 'C', 2, 'O', 'Fe'))['status'] == 'out_of_scope'
    assert phase_stable_at_application_T(('Perovskite', 'La', 'Fe', 'None'))['status'] == 'out_of_scope'
    assert phase_stable_at_application_T(('MoltenMetal', 'Bi', 'Ni', 10.0, 1000))['status'] == 'candidate'
    assert is_turquoise_pyrolysis_candidate("('MOF', 'W', 'Triazolate', 'N2P2', 16.0)") is False
    assert VALIDATION_QUOTA_EXEMPT_CLASSES == frozenset({'MetalHydride'})
    assert validation_quota_class_count(ALL_MATERIAL_CLASSES) == len(ALL_MATERIAL_CLASSES) - 1
    assert validation_quota_class_count(
        ALL_MATERIAL_CLASSES, 'turquoise_pyrolysis') == len(ALL_MATERIAL_CLASSES) - 5


def test_turquoise_pyrolysis_select_excludes_metal_hydride():
    import pandas as pd
    from pipeline.common.application_scope import (
        is_turquoise_pyrolysis_candidate, select_turquoise_pyrolysis_candidates)
    hydride = "('MetalHydride', 'La', 'H2', 'None', 'None', 400)"
    melt = "('MoltenMetal', 'Bi', 'Ni', 10.0, 1000)"
    solid = "('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, ('Fe',), 1, 0)"
    mof = "('MOF', 'W', 'Triazolate', 'N2P2', 16.0)"
    assert is_turquoise_pyrolysis_candidate(hydride) is False
    assert is_turquoise_pyrolysis_candidate(melt) is True
    df = pd.DataFrame([
        {'genome': hydride, 'E_act': 0.05, 'valid': True},
        {'genome': mof, 'E_act': 0.02, 'valid': True},
        {'genome': melt, 'E_act': 0.80, 'valid': True},
        {'genome': solid, 'E_act': 0.90, 'valid': True},
        {'genome': 'not-a-genome', 'E_act': 0.01, 'valid': True},
    ])
    selected = select_turquoise_pyrolysis_candidates(df, top_k=2)
    genomes = set(selected['genome'])
    assert hydride not in genomes
    assert mof not in genomes
    assert 'not-a-genome' not in genomes
    assert list(selected['E_act']) == [0.80, 0.90]


def test_mechanism_has_condensed_graphite_not_gas_carbon():
    from pipeline.process.reactor_mechanisms import write_full_mechanism, write_gas_only_mechanism
    gas_path = write_gas_only_mechanism()
    full_path = write_full_mechanism('test_cat_scope', E_act_CH4=0.9)
    gas_txt = gas_path.read_text(encoding='utf-8')
    full_txt = full_path.read_text(encoding='utf-8')
    assert 'C_graphite' not in gas_txt and 'C_graphite' not in full_txt
    assert 'thermo: fixed-stoichiometry' in gas_txt
    assert 'C(gr)' in gas_txt and 'C(gr)' in full_txt
    assert 'C_s => C_graphite' not in full_txt
    assert 'C_s => C(gr)' not in full_txt
    assert 'name: C_s' in full_txt


def test_off_site_carbon_gated_to_nanoparticle_metals():
    import json
    from pipeline.process.reactor_mechanisms import (
        CandidateKinetics, off_site_carbon_allowed, write_full_mechanism)

    ni = ('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)
    pd = ('SolidCatalyst', 'Pd', 'SiO2', 'fcc111', 0.0, (), 1, 0)
    sac = ('SAC', 'Fe', 'N4', 'N-graphene', 'none')
    saa_host = ('SAA', 'Rh', 'Ni', '111', 1000)
    saa_guest = ('SAA', 'Ni', 'Sn', '111', 1000)
    hea = ('HEA', ('V', 'Fe', 'Ni', 'Ag', 'Hf', 'Al'), 'fcc', '111', 800)
    assert off_site_carbon_allowed(ni) is True
    assert off_site_carbon_allowed(pd) is False
    assert off_site_carbon_allowed(sac) is False
    assert off_site_carbon_allowed(saa_host) is True
    assert off_site_carbon_allowed(saa_guest) is False
    assert off_site_carbon_allowed(hea) is True

    ungated = write_full_mechanism('test_cat_scope', E_act_CH4=0.9)
    ungated_txt = ungated.read_text(encoding='utf-8')
    assert 'C_s => C(gr) + site' not in ungated_txt
    assert 'C_encap_s' not in ungated_txt

    ni_kin = CandidateKinetics(
        methane_activation_eV=0.65, material_class='SolidCatalyst', genome=ni)
    ni_path = write_full_mechanism('ni_np', kinetics=ni_kin)
    ni_txt = ni_path.read_text(encoding='utf-8')
    assert 'C_s => C(gr) + site' in ni_txt
    assert 'C_s => C_encap_s' in ni_txt
    assert 'adjacent-phases: [gas, graphite]' in ni_txt
    assert 'transport-to-edge' in ni_txt
    meta = json.loads(ni_path.with_suffix('.kinetics.json').read_text(encoding='utf-8'))
    assert meta['off_site_carbon'] is True
    assert meta['coking_index_mapped_to_off_site'] is False
    assert meta['off_site_channels']['C_gamma']['barrier_eV'] == 1.5
    assert abs(meta['off_site_channels']['C_delta']['barrier_eV'] - 1.53) < 1e-12

    sac_kin = CandidateKinetics(
        methane_activation_eV=0.43, material_class='SAC', genome=sac)
    sac_path = write_full_mechanism('test_sac_fe', kinetics=sac_kin)
    assert 'C_s => C(gr) + site' not in sac_path.read_text(encoding='utf-8')
    try:
        write_full_mechanism('test_sac_forced', kinetics=sac_kin, off_site_carbon=True)
    except ValueError as exc:
        assert 'Ni/Fe/Co' in str(exc)
    else:
        raise AssertionError('forcing off-site carbon on SAC must fail closed')


def test_surface_thermo_references_and_tst_prefactors():
    """Writer puts adsorption energies on Cantera's absolute scale and uses
    TST prefactors for bimolecular surface steps (B6-5 corrections)."""
    import json
    from pipeline.process.reactor_mechanisms import (
        EV_TO_J_MOL, H_F_CH3_RADICAL_J_MOL, H_F_CH4_J_MOL,
        MONOLAYER_SITE_DENSITY_MOL_CM2, CandidateKinetics, write_full_mechanism)

    ni = ('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)
    kin = CandidateKinetics(
        methane_activation_eV=1.0, h_adsorption_eV=-0.5, ch3_adsorption_eV=-1.95,
        c_adsorption_eV=1.3, encapsulation_crossover_coverage=0.5,
        material_class='SolidCatalyst', genome=ni)
    path = write_full_mechanism('test_ni_refs', kinetics=kin)
    meta = json.loads(path.with_suffix('.kinetics.json').read_text(encoding='utf-8'))
    h = meta['surface_enthalpies_J_mol']
    assert abs(h['H_s'] - (-0.5 * EV_TO_J_MOL)) < 1e-6
    assert abs(h['CH3_s'] - (-1.95 * EV_TO_J_MOL + H_F_CH3_RADICAL_J_MOL)) < 1e-6
    assert abs(h['C_s'] - (1.3 * EV_TO_J_MOL + H_F_CH4_J_MOL)) < 1e-6
    # ladder interpolation on the same scale
    assert abs(h['CH2_s'] - (h['CH3_s'] + (h['C_s'] - h['CH3_s']) / 3)) < 1e-6
    assert abs(h['CH_s'] - (h['CH3_s'] + 2 * (h['C_s'] - h['CH3_s']) / 3)) < 1e-6
    # bimolecular prefactor is k_TST / Gamma in cm^2/mol/s
    a_bimol = meta['surface_prefactors']['bimolecular_cm2_mol_s']
    assert abs(a_bimol - 1e13 / MONOLAYER_SITE_DENSITY_MOL_CM2) / a_bimol < 1e-12
    txt = path.read_text(encoding='utf-8')
    assert txt.count('A: 4e+21') == 3          # three dehydrogenation steps
    assert 'A: 2e+22' in txt                    # H2 desorption
    # C_delta is theta_C^2 via Cantera's power-law coverage dependency (m, not a)
    assert 'C_s: {a: 0.0, m: 1.0, E: 0.0}' in txt
    cd = meta['off_site_channels']['C_delta']
    assert cd['form'] == 'coverage_dependent_theta_C_squared'
    assert cd['crossover_coverage_theta_star'] == 0.5
    assert abs(cd['preexponential_1_s'] - 2e13) < 1e-3   # A_gamma / theta*
    assert meta['off_site_channels']['C_gamma']['preexponential_1_s'] == 1e13

    # the mechanism loads in Cantera and the C_delta rate really is second order
    import cantera as ct
    gas = ct.Solution(str(path), 'gas')
    graphite = ct.Solution(str(path), 'graphite')
    surf = ct.Interface(str(path), 'test_ni_refs_surface', [gas, graphite])
    surf.TP = 923.15, ct.one_atm
    i_d = [r.equation for r in surf.reactions()].index('C_s => C_encap_s')
    surf.coverages = {'site': 0.8, 'C_s': 0.2}
    r1 = surf.net_rates_of_progress[i_d]
    surf.coverages = {'site': 0.6, 'C_s': 0.4}
    r2 = surf.net_rates_of_progress[i_d]
    assert abs(r2 / r1 - 4.0) < 1e-6

    # a declared prefactor is honoured and labelled
    kin2 = CandidateKinetics(
        methane_activation_eV=1.0, carbon_transfer_prefactor_1_s=1e9,
        material_class='SolidCatalyst', genome=ni)
    p2 = write_full_mechanism('test_ni_ag', kinetics=kin2)
    m2 = json.loads(p2.with_suffix('.kinetics.json').read_text(encoding='utf-8'))
    assert m2['off_site_channels']['C_gamma']['preexponential_1_s'] == 1e9
    assert m2['off_site_channels']['C_delta']['preexponential_1_s'] == 2e9
    assert m2['inputs']['provenance']['carbon_transfer_prefactor_1_s'] == 'declared_not_measured'


def test_mechanical_regen_never_clears_encapsulating_carbon():
    import cantera as ct
    from pipeline.process.reactor_mechanisms import CandidateKinetics, write_full_mechanism
    from pipeline.process.reactor_models import _reset_surface_carbon

    ni = ('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)
    path = write_full_mechanism('test_ni_regen', kinetics=CandidateKinetics(
        methane_activation_eV=1.0, material_class='SolidCatalyst', genome=ni))
    gas = ct.Solution(str(path), 'gas')
    graphite = ct.Solution(str(path), 'graphite')
    surf = ct.Interface(str(path), 'test_ni_regen_surface', [gas, graphite])
    surf.coverages = {'site': 0.3, 'C_s': 0.4, 'H_s': 0.1, 'C_encap_s': 0.2}
    _reset_surface_carbon(surf)
    cov = dict(zip(surf.species_names, surf.coverages))
    assert abs(cov['C_encap_s'] - 0.2) < 1e-12
    assert abs(cov['site'] - 0.8) < 1e-12
    assert cov['C_s'] == 0.0 and cov['H_s'] == 0.0


def test_yaml_sweep_kinetics_keys_for_b6():
    import tempfile
    from pathlib import Path
    from pipeline.process.yaml_sweep import parse_sweep, _load_kinetics_row

    body = """
name: b6_keys
catalyst:
  name: ni_keys
  material_class: SolidCatalyst
  genome: "('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)"
  kinetics:
    E_act: 1.0
    dE_H: -0.5
    dE_CH3: -1.95
    dE_C: 1.3
    carbon_transfer_eV: 1.5
    carbon_encapsulation_eV: 1.53
    encapsulation_crossover_coverage: 0.5
    carbon_transfer_prefactor_1_s: 1.0e9
    provenance:
      E_act: Bengaard 2002
conditions:
  temperatures_K: [923.15]
  reactors: [PFR]
cells:
  - name: c
    catalyst_particle_mm: 0.13
    metal_loading: 0.5
    metal_dispersion: 0.3
"""
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / 'b6.yaml'
        p.write_text(body, encoding='utf-8')
        job = parse_sweep(p)
        assert job.kinetics['carbon_transfer_prefactor_1_s'] == 1e9
        assert job.kinetics['encapsulation_crossover_coverage'] == 0.5
        _, kin, cls, cid = _load_kinetics_row(job)
        assert cls == 'SolidCatalyst' and cid == 'ni_keys'
        assert kin.carbon_transfer_prefactor_1_s == 1e9
        assert kin.sources['methane_activation_eV'] == 'sweep_yaml: Bengaard 2002'
        assert kin.sources['carbon_transfer_prefactor_1_s'] == 'sweep_yaml'
        # unknown kinetics keys fail closed
        p.write_text(body.replace('carbon_transfer_prefactor_1_s: 1.0e9', 'k_bogus: 2'), encoding='utf-8')
        try:
            parse_sweep(p)
        except ValueError as exc:
            assert 'k_bogus' in str(exc)
        else:
            raise AssertionError('unknown kinetics key must be rejected')


def test_staged_sweep_preserves_coarse_and_proposes_roi():
    import tempfile
    from pathlib import Path
    from pipeline.process import staged_sweep
    old = staged_sweep.SWEEPS_DIR
    with tempfile.TemporaryDirectory() as tmp:
        staged_sweep.SWEEPS_DIR = Path(tmp)
        try:
            spec = staged_sweep.SweepSpec(
                name='demo',
                levers={'x': [1.0, 2.0, 4.0], 'y': [0.2, 1.0]},
                score_key='score',
                constraint_key='ok',
                keep_fraction_of_max=0.5,
                n_targeted_per_lever=4,
                hard_bounds={'x': (0.5, 8.0), 'y': (0.0, 1.0)},
            )
            records = [
                {'x': 1.0, 'y': 1.0, 'score': 0.1, 'ok': True},
                {'x': 2.0, 'y': 1.0, 'score': 0.4, 'ok': True},
                {'x': 4.0, 'y': 1.0, 'score': 1.0, 'ok': True},
                {'x': 4.0, 'y': 0.2, 'score': 0.2, 'ok': True},
                {'x': 8.0, 'y': 1.0, 'score': 2.0, 'ok': False},
            ]
            coarse = {'records': records, 'levels': spec.levers, 'n_grid_cells': 6}
            staged_sweep.write_stage('demo', 'coarse', coarse)
            try:
                staged_sweep.write_stage('demo', 'coarse', coarse)
            except FileExistsError:
                pass
            else:
                raise AssertionError('coarse sweep must be write-once')
            roi = staged_sweep.propose_roi(coarse, spec)
            # Peak feasible cell is x=4 (last coarse level). extend_edge
            # steps past that outer edge toward the hard wall; it does not
            # pull the lower bound inward.
            assert roi['bounds']['x'][0] <= 4.0
            assert roi['bounds']['x'][1] > 4.0
            staged_sweep.write_stage('demo', 'targeted', {
                'records': [{'score': 1.1}], 'roi': roi, 'n_grid_cells': roi['n_grid_cells']})
            bundle = staged_sweep.load_sweep('demo')
            assert bundle['coarse'] is not None
            assert bundle['targeted'] is not None
            assert bundle['coarse']['n_grid_cells'] == 6
            assert bundle['manifest']['coarse_preserved'] is True
        finally:
            staged_sweep.SWEEPS_DIR = old


def test_inventory_levers_preserve_baseline_area():
    from pipeline.process.reactor_models import (
        ReactorConfig, _validate_carbon_policy, active_sv, ergun_delta_p_pa,
        geometric_sv_pfr, inventory_grid_cells)
    from pipeline.process.reactor_models import (
        DEFAULT_METAL_DISPERSION, DEFAULT_METAL_LOADING,
        DEFAULT_SOLIDS_PARTICLE_MM)
    assert abs(DEFAULT_SOLIDS_PARTICLE_MM - 0.13) < 1e-15
    assert abs(DEFAULT_METAL_LOADING - 0.5) < 1e-15
    assert abs(DEFAULT_METAL_DISPERSION - 0.3) < 1e-15
    cfg = ReactorConfig(
        catalyst_particle_mm=2.0, metal_loading=1.0, metal_dispersion=1.0)
    _validate_carbon_policy(cfg)
    assert abs(geometric_sv_pfr(cfg) - 1800.0) < 1e-9
    assert abs(active_sv(geometric_sv_pfr(cfg), cfg) - 1800.0) < 1e-9
    small = ReactorConfig(
        catalyst_particle_mm=0.2, metal_loading=1.0, metal_dispersion=1.0)
    prod = ReactorConfig()
    assert abs(prod.catalyst_particle_mm - 0.13) < 1e-15
    assert abs(prod.metal_loading - 0.5) < 1e-15
    assert abs(prod.metal_dispersion - 0.3) < 1e-15
    assert geometric_sv_pfr(prod) > geometric_sv_pfr(cfg)
    assert abs(active_sv(geometric_sv_pfr(prod), prod)
               / geometric_sv_pfr(prod) - 0.15) < 1e-9
    assert abs(geometric_sv_pfr(small) / geometric_sv_pfr(cfg) - 10.0) < 1e-9
    half = ReactorConfig(
        catalyst_particle_mm=2.0, metal_loading=0.5, metal_dispersion=1.0)
    assert abs(active_sv(geometric_sv_pfr(half), half) - 900.0) < 1e-9
    assert ergun_delta_p_pa(small) > ergun_delta_p_pa(cfg)
    cells = inventory_grid_cells()
    assert len(cells) == 36
    assert cells[0] == {
        'catalyst_particle_mm': 2.0, 'metal_loading': 1.0, 'metal_dispersion': 1.0}
    from pipeline.process.reactor_models import inventory_roi_grid_cells
    roi = inventory_roi_grid_cells()
    assert len(roi) == 54
    assert roi[0]['catalyst_particle_mm'] == 0.25
    assert max(c['catalyst_particle_mm'] for c in roi) <= 0.25
    assert min(c['catalyst_particle_mm'] for c in roi) >= 0.08
    try:
        _validate_carbon_policy(ReactorConfig(metal_loading=1.5))
    except ValueError as exc:
        assert 'invent area' in str(exc)
    else:
        raise AssertionError('loading > 1 must fail closed')


def test_yaml_sweep_parses_headline_example():
    from pipeline.common.utils import BASE_DIR
    from pipeline.process.yaml_sweep import parse_sweep
    job = parse_sweep(BASE_DIR / 'sweeps' / 'headline_cat9_1300K.yaml')
    assert job.name == 'headline_cat9_1300K'
    assert job.catalyst_name == 'cat_9'
    assert job.screening_index == 9
    assert job.temperatures_K == [1300.0]
    assert job.reactor_types == ['PFR', 'Fluidized', 'MMBCR']
    assert job.policy['co2_permitted'] is False
    assert job.policy['fluidized_mode'] == 'circulating'
    assert job.policy['max_regen_cycles'] == 3
    assert job.policy['regen_mechanism'] == 'mechanical'
    assert len(job.cells) == 2
    frac, env = job.cells
    assert frac.name == 'fractional'
    assert abs(frac.catalyst_particle_mm - 0.13) < 1e-15
    assert abs(frac.metal_loading - 0.5) < 1e-15
    assert abs(frac.metal_dispersion - 0.3) < 1e-15
    assert env.name == 'envelope_1x1'
    assert abs(env.metal_loading * env.metal_dispersion - 1.0) < 1e-15
    template = (BASE_DIR / 'docs' / 'sweep-template.md').read_text(encoding='utf-8')
    assert 'python runsweep.py' in template
    assert 'headline_cat9_1300K.yaml' in template
    try:
        parse_sweep(BASE_DIR / 'sweeps' / 'headline_cat9_1300K.xml')
    except ValueError as exc:
        assert 'YAML' in str(exc)
    else:
        raise AssertionError('leftover XML specs must fail closed')


def test_yaml_sweep_rejects_invented_area():
    import tempfile
    from pathlib import Path
    from pipeline.process.yaml_sweep import parse_sweep
    spec = """
name: bad_loading
catalyst:
  name: explicit
  material_class: SolidCatalyst
  kinetics:
    E_act: 0.43
    dE_H: -0.90
conditions:
  temperatures_K: [1300]
  reactors: [PFR]
cells:
  - name: overloaded
    catalyst_particle_mm: 0.13
    metal_loading: 1.5
    metal_dispersion: 1.0
"""
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'bad.yaml'
        path.write_text(spec, encoding='utf-8')
        try:
            parse_sweep(path)
        except ValueError as exc:
            assert 'metal_loading' in str(exc)
        else:
            raise AssertionError('loading > 1 must fail closed')
        # Kinetics-only sweeps must declare the class: reactor applicability
        # and the B6 gate depend on it, and upstream refuses a reactor
        # without one.
        no_class = spec.replace('  material_class: SolidCatalyst\n', '')
        no_class = no_class.replace('metal_loading: 1.5', 'metal_loading: 0.5')
        path.write_text(no_class, encoding='utf-8')
        try:
            parse_sweep(path)
        except ValueError as exc:
            assert 'material_class' in str(exc)
        else:
            raise AssertionError('kinetics-only sweep without a class must fail closed')


def test_yaml_sweep_splits_reactors_by_pathway_mode():
    from pipeline.process.yaml_sweep import _reactors_by_mode
    groups = dict(_reactors_by_mode(['PFR', 'Fluidized', 'MMBCR']))
    assert groups == {
        'thermocatalytic_pfr': ['PFR'],
        'thermocatalytic_fluidized': ['Fluidized'],
        'mmbcr': ['MMBCR'],
    }


def test_solids_scorecard_judges_cat_9_not_h_parked():
    from pipeline.process.phase2_scorecard import build_solids_scorecard, is_h_parked
    h_parked = {
        'reactor_type': 'PFR', 'catalyst_name': 'cat_40', 'T_K': 1300.0,
        'CH4_conversion': 0.003, 'single_pass_CH4_conversion': 0.003,
        'catalyst_E_act_eV': 0.01, 'catalyst_dE_H_eV': -2.5,
        'active_sv_1_m': 4153.8, 'WHSV_h-1': 900.0,
        'ergun_delta_p_bar': 0.40, 'ergun_ok': True,
        'surface_loaded': True,
    }
    judge = {
        'reactor_type': 'PFR', 'catalyst_name': 'cat_9', 'T_K': 1300.0,
        'CH4_conversion': 0.0102, 'single_pass_CH4_conversion': 0.0102,
        'catalyst_E_act_eV': 0.43, 'catalyst_dE_H_eV': -0.90,
        'active_sv_1_m': 4153.8, 'WHSV_h-1': 900.0,
        'ergun_delta_p_bar': 0.40, 'ergun_ok': True,
        'surface_loaded': True,
    }
    fluid = dict(judge)
    fluid.update({
        'reactor_type': 'Fluidized',
        'single_pass_CH4_conversion': 0.0122,
        'CH4_conversion': 0.0122,
        'WHSV_h-1': 180.0,
    })
    melt = {
        'reactor_type': 'MMBCR', 'catalyst_name': 'cat_9', 'T_K': 1300.0,
        'CH4_conversion': 0.985, 'catalyst_E_act_eV': 0.43, 'catalyst_dE_H_eV': -0.90,
    }
    assert is_h_parked(h_parked)
    assert not is_h_parked(judge)
    card = build_solids_scorecard(
        [h_parked, judge, fluid, melt], judge_catalyst='cat_9')
    assert card['judge_catalyst'] == 'cat_9'
    assert card['judge_catalyst_requested'] == 'cat_9'
    assert abs(card['headline_solids_conversion'] - 0.0102) < 1e-9
    assert card['headline']['PFR']['WHSV_h-1'] == 900.0
    assert card['headline']['Fluidized']['single_pass_CH4_conversion'] == 0.0122
    assert abs(card['mmbcr_max_conversion'] - 0.985) < 1e-9
    assert card['solids_max_excluding_h_parked']['catalyst_name'] == 'cat_9'
    unnamed = build_solids_scorecard([h_parked, judge, fluid, melt])
    assert unnamed['judge_catalyst_requested'] is None
    assert unnamed['judge_catalyst'] is None
    assert unnamed['headline_catalyst'] == 'cat_9'


def test_solids_run_requires_loaded_surface():
    from pipeline.process.phase2_scorecard import is_solids_run, build_solids_scorecard
    base = {
        'reactor_type': 'PFR', 'catalyst_name': 't', 'T_K': 1300.0,
        'CH4_conversion': 1e-5, 'single_pass_CH4_conversion': 1e-5,
    }
    assert not is_solids_run(base)
    assert not is_solids_run({**base, 'surface_loaded': False})
    assert not is_solids_run({**base, 'surface_loaded': True, 'mock': True})
    assert not is_solids_run({**base, 'surface_loaded': True, 'gas_only': True})
    assert is_solids_run({**base, 'surface_loaded': True})
    blank = {**base, 'surface_loaded': False, 'gas_only': False}
    card = build_solids_scorecard([blank])
    assert card['n_solids_records'] == 0
    assert card['headline_solids_conversion'] is None


def test_mismatched_catalyst_name_fails_closed():
    try:
        import cantera  # noqa: F401
    except ImportError:
        return
    from pipeline.process.reactor_mechanisms import write_full_mechanism
    from pipeline.process.reactor_models import ReactorConfig, simulate_pfr
    path = write_full_mechanism('t_0_05', E_act_CH4=0.9)
    cfg = ReactorConfig(
        mechanism_file=str(path), catalyst_name='t',
        reactor_type='PFR', T_inlet_K=1000.0)
    try:
        simulate_pfr(cfg)
    except RuntimeError as exc:
        assert 't_surface' in str(exc)
        assert 'gas_only' in str(exc)
    else:
        raise AssertionError('name mismatch must not return a blank X')


def test_ch4_conversion_uses_argon_tracer_when_c2_present():
    from pipeline.process.equilibrium_check import (
        ch4_conversion_from_argon_tracer, ch4_conversion_from_mole_fractions)
    # Feed CH4:0.95 / Ar:0.05. Solid route p, C2 route q (each CH4 → 0.5 C2H6 + 0.5 H2).
    x_ch4_0, x_ar_0 = 0.95, 0.05
    p, q = 0.00, 0.10
    n_ch4_0, n_ar = 0.95, 0.05
    n_ch4 = n_ch4_0 * (1.0 - p - q)
    n_c2h6 = 0.5 * n_ch4_0 * q
    n_h2 = 2.0 * n_ch4_0 * p + 0.5 * n_ch4_0 * q
    n_tot = n_ch4 + n_c2h6 + n_h2 + n_ar
    true_x = p + q
    got = ch4_conversion_from_argon_tracer(
        n_ch4 / n_tot, n_ar / n_tot, x_ch4_0, x_ar_0)
    assert abs(got - true_x) < 1e-12, (got, true_x)
    wrong = ch4_conversion_from_mole_fractions(n_ch4 / n_tot, n_h2 / n_tot)
    assert wrong < true_x - 0.05
    # Pure CH4/H2/Ar still recovers X on the H2 formula.
    X = 0.5
    n_ch4, n_h2, n_ar = 0.475, 0.95, 0.05
    n_tot = n_ch4 + n_h2 + n_ar
    assert abs(ch4_conversion_from_mole_fractions(n_ch4 / n_tot, n_h2 / n_tot) - X) < 1e-12
    assert abs(ch4_conversion_from_argon_tracer(
        n_ch4 / n_tot, n_ar / n_tot, 0.95, 0.05) - X) < 1e-12


def test_mmbcr_k0_and_flotation_fail_closed():
    from pipeline.process.reactor_models import ReactorConfig, _validate_carbon_policy
    _validate_carbon_policy(ReactorConfig())
    try:
        _validate_carbon_policy(ReactorConfig(mmbcr_interfacial_k0_m_s=0.0))
    except ValueError as exc:
        assert 'guardrail' in str(exc)
    else:
        raise AssertionError('k0=0 must fail closed')
    try:
        _validate_carbon_policy(ReactorConfig(mmbcr_interfacial_k0_m_s=5.0))
    except ValueError as exc:
        assert 'guardrail' in str(exc)
    else:
        raise AssertionError('huge k0 must fail closed')
    try:
        _validate_carbon_policy(ReactorConfig(mmbcr_carbon_removal_rate_1_s=-1.0))
    except ValueError as exc:
        assert 'mmbcr_carbon_removal' in str(exc)
    else:
        raise AssertionError('negative flotation rate must fail closed')


def test_site_density_locked_to_monolayer():
    from pipeline.process.reactor_mechanisms import (
        MONOLAYER_SITE_DENSITY_MOL_CM2, write_full_mechanism)
    from pipeline.process.reactor_models import ReactorConfig, _validate_carbon_policy
    assert abs(MONOLAYER_SITE_DENSITY_MOL_CM2 - 2.5e-9) < 1e-15
    path = write_full_mechanism('test_gamma_lock', E_act_CH4=0.9)
    assert 'site-density: 2.500e-09 mol/cm^2' in path.read_text(encoding='utf-8')
    try:
        write_full_mechanism('test_gamma_lock', E_act_CH4=0.9, site_density=1e-6)
    except ValueError as exc:
        assert 'B1 locks' in str(exc)
    else:
        raise AssertionError('raised site_density must fail closed')
    try:
        _validate_carbon_policy(ReactorConfig(site_density_mol_cm2=1e-6))
    except ValueError as exc:
        assert 'B1 locks' in str(exc)
    else:
        raise AssertionError('ReactorConfig site_density override must fail closed')


def test_oxidative_regen_requires_co2_permitted():
    from pipeline.process.reactor_models import ReactorConfig, simulate_pfr
    cfg = ReactorConfig(
        reactor_type='PFR', catalyst_name='x', mechanism_file='',
        max_regen_cycles=1, regen_mechanism='oxidative', co2_permitted=False)
    try:
        simulate_pfr(cfg)
    except RuntimeError as exc:
        assert 'co2_permitted' in str(exc)
    else:
        raise AssertionError('oxidative regen should be blocked')



def test_prior_art_registry_tracks_exact_and_region_novelty():
    import tempfile
    from pathlib import Path
    from pipeline.evidence.prior_art import PriorArtRegistry
    known = ('SAC', 'Fe', 'N4', 'N-graphene', 'OH')
    related = ('SAC', 'Fe', 'N4', 'N-graphene', 'none')
    unseen = ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000)
    with tempfile.TemporaryDirectory() as tmp:
        registry = PriorArtRegistry(str(Path(tmp) / 'prior.sqlite'))
        registry.add(known, 'literature', 'doi:test')
        registry.add(known, 'patent', 'patent:test')
        assert registry.count() == 2
        assert len(registry.classify(known)['exact_records']) == 2
        assert registry.classify(known)['novelty_status'] == 'known'
        assert registry.classify(related)['novelty_status'] == 'region_known'
        assert registry.classify(unseen)['novelty_status'] == 'unseen'


def test_multiobjective_archive_preserves_conflicting_winners():
    import sqlite3
    import tempfile
    from pathlib import Path
    from pipeline.search.exhaustive_search import ScanConfig, run_streaming_scan
    from pipeline.search.indexed_space import CLASS_OFFSETS
    with tempfile.TemporaryDirectory() as tmp:
        db = str(Path(tmp) / 'multi.sqlite')
        start = CLASS_OFFSETS['SAC']
        def scorer(genomes):
            n = len(genomes)
            return np.column_stack([np.arange(n), np.arange(n)[::-1],
                                    np.zeros(n), np.ones(n)])
        run_streaming_scan(ScanConfig('multi', db, start=start, stop=start + 20,
                                      batch_size=20, global_archive_size=20,
                                      state_id='multi-test'), scorer)
        conn = sqlite3.connect(db)
        objectives = {r[0] for r in conn.execute(
            "SELECT DISTINCT objective_index FROM objective_archive WHERE application='multi'")}
        regions = {r[0] for r in conn.execute(
            "SELECT DISTINCT objective_index FROM regional_objective_champions WHERE application='multi'")}
        conn.close()
        assert objectives == {0, 1, 2, 3}
        assert regions == {0, 1, 2, 3}


def test_final_campaign_readiness_fails_closed():
    import json
    import tempfile
    from pathlib import Path
    from pipeline.search.indexed_space import TOTAL_SIZE
    from pipeline.evidence.prior_art import PriorArtRegistry
    from pipeline.evidence.readiness import campaign_readiness
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        cert = root / 'coverage.json'
        prior = root / 'prior.sqlite'
        result = campaign_readiness(str(cert), str(prior))
        assert not result['ready']
        assert 'coverage_certificate_missing' in result['failures']
        assert 'prior_art_registry_empty' in result['failures']
        cert.write_text(json.dumps({
            'declared_encoded_population': TOTAL_SIZE, 'complete': True}))
        PriorArtRegistry(str(prior)).add(
            ('SAC', 'Fe', 'N4', 'N-graphene', 'OH'), 'literature', 'doi:test')
        assert campaign_readiness(str(cert), str(prior))['ready']
        manifest = root / 'evidence.json'
        gated = campaign_readiness(str(cert), str(prior),
                                   evidence_manifest=str(manifest),
                                   application='turquoise_hydrogen')
        assert not gated['ready'] and 'evidence_manifest_missing' in gated['failures']
        import hashlib
        records = {}
        for key in ('converged_dft', 'measured_reactor',
                    'measured_deactivation', 'ntec_control_pair'):
            artifact = root / f'{key}.json'
            artifact.write_text(json.dumps({'kind': key, 'verified': True}))
            records[key] = [{
                'candidate_id': f'candidate:{key}',
                'source_path': artifact.name,
                'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest(),
                'protocol_id': 'test-protocol-v1',
                'status': 'converged' if key == 'converged_dft' else 'measured',
            }]
        manifest.write_text(json.dumps({
            'schema_version': 2, 'records': records}))
        assert campaign_readiness(str(cert), str(prior),
                                  evidence_manifest=str(manifest),
                                  application='turquoise_hydrogen')['ready']
        no_ntec = json.loads(manifest.read_text())
        no_ntec['records']['ntec_control_pair'] = []
        manifest.write_text(json.dumps(no_ntec))
        assert campaign_readiness(
            str(cert), str(prior), evidence_manifest=str(manifest),
            application='turquoise_hydrogen',
            pyrolysis_mode='thermocatalytic')['ready']
        assert not campaign_readiness(
            str(cert), str(prior), evidence_manifest=str(manifest),
            application='turquoise_hydrogen', pyrolysis_mode='ntec')['ready']


# ═══════════════════════════════════════════════════════════════════════════════
# RUN ALL
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    t0 = time.time()
    print("=" * 60)
    print("  HYDROGEN PIPELINE — TEST SUITE")
    print("=" * 60)

    print("\n── Design Space ──")
    test("14 classes generate", test_all_14_classes_generate)
    test("No toxic elements", test_no_toxic_elements)
    test("SAC axial ligands", test_sac_has_axial_ligand)
    test("Class weights sum to one", test_class_weights_sum_to_1)

    print("\n── Encoding ──")
    test("All classes encode (no NaN)", test_encode_all_classes_no_nan)
    test("FEATURE_DIM matches components", test_feature_dim_matches_components)

    print("\n── Surrogates ──")
    test("CH4 surrogate (no NaN)", test_ch4_surrogate_no_nan)
    test("ORR surrogate (no NaN)", test_orr_surrogate_no_nan)

    print("\n── Pareto Objectives ──")
    test("Pareto sorting is correct", test_nsga2_sorts_correctly)
    test("Cost & Fenton in range", test_cost_and_fenton_ranges)
    test("MetalFreeCarbon cost = 0", test_metalfreecarbon_zero_cost)
    test("PEMFC application scope", test_pemfc_application_scope)
    test("Novelty time-split benchmark", test_novelty_time_split_benchmark)
    test("Pilot benchmark candidate deduplication", test_pilot_benchmark_deduplicates_candidates)
    test("Small-data rankers preserve continuous targets", test_small_data_rankers_preserve_continuous_targets)
    test("Six-point status fails closed", test_six_point_status_fails_closed)
    test("Adaptive validation policy", test_adaptive_validation_policy)
    test("SSSP and candidate NEB workflow", test_sssp_and_candidate_neb_workflow)
    test("ORR multisite and corrections", test_orr_multisite_and_corrections)
    test("Production QE workflow fails closed and resumes",
         test_production_qe_workflow_fails_closed_and_resumes)
    test("Validation task queue is candidate-keyed and resumable",
         test_validation_task_queue_is_candidate_keyed_and_resumable)

    print("\n── Element Extractors ──")
    test("4 extractors consistent", test_element_extractors_consistent)

    print("\n── Structure Generation ──")
    test("All 14 classes build structures", test_structure_generation_all_classes)
    test("Structures are deterministic and genome-sensitive",
         test_structure_generation_is_deterministic_and_genome_sensitive)

    print("\n── PEMFC Model ──")
    test("Tafel covers all classes", test_tafel_covers_all_classes)
    test("Power monotonic with eta", test_pemfc_power_monotonic_with_eta)
    test("Tafel slope affects power", test_tafel_slope_affects_power)
    test("Efficiency in physical range", test_pemfc_efficiency_in_range)

    print("\n── Stack Model ──")
    test("Stack model", test_stack_model)
    test("TEA is scenario-labelled and unclamped",
         test_tea_is_scenario_labelled_and_unclamped)

    print("\n── CHE / Utils ──")
    test("ORR ideal overpotential = 0", test_orr_overpotential_ideal)
    test("Abundance cost penalty", test_abundance_cost_penalty)

    print("\n── Report Generator ──")
    test("Campaign path", test_report_campaign_path)
    test("Orchestrator path", test_report_orchestrator_path)
    test("Empty state (no crash)", test_report_empty_no_crash)

    print("\n── OOD Confidence ──")
    test("High confidence for metals", test_ood_high_confidence_metals)
    test("Low confidence for OOD classes", test_ood_low_confidence_ood)
    test("Penalty scales objectives", test_ood_penalty_scales_objectives)
    test("Confidence affects objectives", test_ood_nsga2_integration)

    print("\n── Exhaustive Coverage ──")
    test("All elements in abundance table", test_all_elements_in_abundance_table)
    test("All elements in price table", test_all_elements_in_price_table)
    test("All 14 classes viable both apps", test_all_classes_viable_both_applications)
    test("BEP params all 14 classes", test_bep_params_all_classes)
    test("OOD confidence all 14 classes", test_ood_confidence_all_classes)
    test("Tafel slope all 14 classes", test_tafel_all_classes)
    test("Cathode SAC genomes 5-tuple", test_cathode_sac_genome_5tuple)
    test("Pyrolysis mode coking bonus", test_pyrolysis_mode_coking_bonus)
    test("Discovery batch covers unseen regions", test_discovery_batch_prioritizes_unseen_regions)
    test("Canonical candidate IDs", test_candidate_ids_are_canonical)
    test("Crossover preserves class", test_crossover_preserves_class)
    test("Mutation preserves class", test_mutation_preserves_class)
    test("Deterministic hierarchical pool", test_deterministic_hierarchical_pool)
    test("Hierarchical rounds cover complementary cells",
         test_hierarchical_rounds_cover_complementary_cells)
    test("Design space remains sizable and provenance-backed",
         test_design_space_audit_preserves_all_sizable_classes)
    test("Discovery metadata persists", test_discovery_metadata_is_persistable)
    test("Indexed space boundaries", test_indexed_space_boundaries_and_classes)
    test("Indexed worker shards", test_indexed_worker_shards_are_disjoint)
    test("Streaming scan resumes", test_streaming_scan_resumes_without_rescoring)
    test("Sharded scan matches serial and fails closed",
         test_sharded_scan_matches_serial_and_fails_closed)
    test("Branch search never surrogate-prunes", test_branch_search_resolves_without_surrogate_pruning)
    test("Branch probes use low-discrepancy schedule", test_branch_probes_are_deterministic_low_discrepancy)
    test("Branch finite budget preserves class floor", test_branch_finite_budget_preserves_class_floor)
    test("Branch certificate detects gaps", test_branch_certificate_detects_incomplete_and_gaps)
    test("Branch rejects population mismatch", test_branch_rejects_population_mismatch)
    test("Tree probes deterministic across 14 classes", test_tree_calibration_probes_cover_all_classes_deterministically)
    test("Production search is branch-only", test_production_has_only_branch_candidate_search)
    test("README matches branch-only contract", test_readme_matches_branch_only_contract)
    test("README has no machine-specific paths",
         test_readme_contains_no_machine_specific_paths)
    test("Root documentation is canonical", test_root_documentation_is_canonical)
    test("Retired GA entry points are blocked", test_retired_ga_entry_points_are_blocked)
    test("Industrial viability gates fail closed", test_industrial_viability_gates_fail_closed)
    test("Phase stability per class", test_phase_stable_at_application_t_per_class)
    test("Pyrolysis select excludes unstable phases", test_turquoise_pyrolysis_select_excludes_metal_hydride)
    test("Slab coking excludes melts", test_slab_coking_scope_excludes_molten_metal)
    test("Mechanism uses condensed graphite", test_mechanism_has_condensed_graphite_not_gas_carbon)
    test("Off-site carbon gated to nanoparticle metals", test_off_site_carbon_gated_to_nanoparticle_metals)
    test("Surface thermo references and TST prefactors", test_surface_thermo_references_and_tst_prefactors)
    test("Mechanical regen never clears C_encap_s", test_mechanical_regen_never_clears_encapsulating_carbon)
    test("YAML sweep B6 kinetics keys", test_yaml_sweep_kinetics_keys_for_b6)
    test("Site density locked to monolayer", test_site_density_locked_to_monolayer)
    test("Inventory levers preserve baseline area", test_inventory_levers_preserve_baseline_area)
    test("YAML sweep parses headline example", test_yaml_sweep_parses_headline_example)
    test("YAML sweep rejects invented area", test_yaml_sweep_rejects_invented_area)
    test("YAML sweep splits reactors by pathway mode", test_yaml_sweep_splits_reactors_by_pathway_mode)
    test("Solids scorecard takes named judge as argument", test_solids_scorecard_judges_cat_9_not_h_parked)
    test("Solids run requires loaded surface", test_solids_run_requires_loaded_surface)
    test("Mismatched catalyst name fails closed", test_mismatched_catalyst_name_fails_closed)
    test("CH4 conversion uses Ar tracer", test_ch4_conversion_uses_argon_tracer_when_c2_present)
    test("MMBCR k0 and flotation fail closed", test_mmbcr_k0_and_flotation_fail_closed)
    test("Staged sweep keeps coarse and proposes ROI", test_staged_sweep_preserves_coarse_and_proposes_roi)
    test("Oxidative regen requires co2_permitted", test_oxidative_regen_requires_co2_permitted)
    test("Prior-art novelty states", test_prior_art_registry_tracks_exact_and_region_novelty)
    test("Multi-objective archive keeps conflicting winners", test_multiobjective_archive_preserves_conflicting_winners)
    test("Final campaign readiness fails closed", test_final_campaign_readiness_fails_closed)

    elapsed = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"  {PASS} passed, {FAIL} failed ({elapsed:.1f}s)")
    if FAIL == 0:
        print(f"  🟢 PIPELINE IS READY TO LAUNCH")
    else:
        print(f"  🔴 FIX {FAIL} FAILURE(S) BEFORE LAUNCHING")
        for name, tb in ERRORS:
            print(f"\n{'─' * 40}")
            print(f"FAILED: {name}")
            print(tb)
    print("=" * 60)

    sys.exit(0 if FAIL == 0 else 1)
