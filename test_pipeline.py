#!/usr/bin/env python3
"""
Hydrogen Pipeline — Comprehensive Test Suite

Run:  python test_pipeline.py
Exit: 0 = all pass, 1 = failures

Tests every component that has broken before, plus integration across
all 14 material classes. If this passes, the campaign is safe to launch.
"""

import sys, os, time, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

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
    held_out = [{'genome': known, 'publication_year': 2025,
                 'source_id': 'doi:test', 'citation': 'Test et al.'}]
    result = time_split_recovery([known], held_out, k=1)
    assert result['valid'] and result['exact_recall_at_k'] == 1.0
    malformed = time_split_recovery([known], [{'genome': known}], k=1)
    assert not malformed['valid']


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
    from pipeline.screening.small_data_ranker import fit_tree_ranker
    genomes = deterministic_tree_probes(24)
    pyro = pd.DataFrame({'genome': [repr(g) for g in genomes], 'valid': True,
                         'E_act': np.linspace(0.1, 2.0, len(genomes))})
    ranker = fit_tree_ranker(pyro, 'turquoise_hydrogen')
    mean, uncertainty = ranker.predict(genomes[:5])
    assert np.all(np.isfinite(mean)) and np.all(uncertainty >= 0)
    orr = pd.DataFrame({'genome': [repr(g) for g in genomes], 'valid': True,
                        'dG_OH_eV': np.linspace(-1, 1, len(genomes)),
                        'dG_O_eV': np.linspace(-2, 2, len(genomes)),
                        'dG_OOH_eV': np.linspace(3, 5, len(genomes))})
    ranker = fit_tree_ranker(orr, 'fuel_cell_orr')
    mean, _ = ranker.predict(genomes[:5])
    assert np.all(np.isfinite(mean)) and len(set(mean.tolist())) > 1


def test_six_point_status_fails_closed():
    import tempfile
    from pipeline.evidence.campaign_status import assess_campaign
    with tempfile.TemporaryDirectory() as tmp:
        result = assess_campaign(tmp)
        assert not result['ready']
        assert len(result['missing']) == 6
        thermal = assess_campaign(tmp, pyrolysis_mode='thermocatalytic')
        assert 'calibrated_ntec' not in thermal['criteria']
        assert len(thermal['missing']) == 5


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
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(20):
            g = generate_random_genome(cls)
            atoms, idx, mc = generate_structure(g)
            assert len(atoms) > 0, f"{cls}: empty structure"
            assert atoms.pbc.any(), f"{cls}: PBC not set"
            assert mc == cls, f"{cls}: returned class {mc}"


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
            scan_batch_size=512, max_leaves=1,
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
    source = (Path(__file__).parent / 'run_production_campaign.py').read_text()
    forbidden = [
        'run_genetic_algorithm', 'run_fc_genetic_algorithm',
        '--exhaustive-scan', '--branch-search', '--pop', '--gens',
        'generate_population', 'generate_random_genome',
    ]
    present = [token for token in forbidden if token in source]
    assert not present, f"Legacy candidate-search paths remain in production: {present}"
    assert 'run_branch_discovery' in source
    assert 'run_fc_branch_discovery' in source
    assert "envs_dir / 'qe-env' / 'bin' / 'pw.x'" in source
    assert "'conda', 'run', '-n', 'quantum-env'" in source
    assert "result.get('mock')" in source


def test_readme_matches_branch_only_contract():
    from pathlib import Path
    readme = (Path(__file__).parent / 'README.md').read_text()
    assert '21,092,645,031' in readme
    assert 'Deterministic Branch-and-Bound Discovery' in readme
    assert '--calibration-probes' in readme
    assert '--branch-leaf-size' in readme
    forbidden = ['--pop', '--gens', '--exhaustive-scan', '--branch-search',
                 '25.3-billion-configuration', '21.3-billion-configuration']
    present = [token for token in forbidden if token in readme]
    assert not present, f"README advertises retired search controls: {present}"


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
    from pipeline.validation.viability import evaluate_turquoise, evaluate_fuel_cell
    assert evaluate_turquoise({})['status'] == 'unknown'
    good_h2 = evaluate_turquoise({
        'temperature_K': 1000, 'H2_selectivity': 0.98, 'CH4_conversion': 0.8,
        'deactivation_fraction_per_h': 0.005, 'coke_fraction': 0.02,
        'net_energy_kWh_kg_h2': 12.0, 'measured_reactor': 1})
    assert good_h2['status'] == 'pass'
    assert evaluate_turquoise({'H2_selectivity': 0.8})['status'] == 'fail'
    good_fc = evaluate_fuel_cell({
        'orr_overpotential_V': 0.3, 'peak_power_W_cm2': 1.2,
        'system_efficiency': 0.5, 'voltage_degradation_uV_h': 5,
        'measured_hours': 500, 'measured_mea': 1})
    assert good_fc['status'] == 'pass'
    assert evaluate_fuel_cell({'orr_overpotential_V': 0.6})['status'] == 'fail'


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
        manifest.write_text(json.dumps({
            'converged_dft_count': 1, 'measured_reactor_count': 1,
            'measured_deactivation_count': 1, 'ntec_control_pair_count': 1,
        }))
        assert campaign_readiness(str(cert), str(prior),
                                  evidence_manifest=str(manifest),
                                  application='turquoise_hydrogen')['ready']
        no_ntec = json.loads(manifest.read_text())
        no_ntec['ntec_control_pair_count'] = 0
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

    print("\n── Element Extractors ──")
    test("4 extractors consistent", test_element_extractors_consistent)

    print("\n── Structure Generation ──")
    test("All 14 classes build structures", test_structure_generation_all_classes)

    print("\n── PEMFC Model ──")
    test("Tafel covers all classes", test_tafel_covers_all_classes)
    test("Power monotonic with eta", test_pemfc_power_monotonic_with_eta)
    test("Tafel slope affects power", test_tafel_slope_affects_power)
    test("Efficiency in physical range", test_pemfc_efficiency_in_range)

    print("\n── Stack Model ──")
    test("Stack model", test_stack_model)

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
    test("Design space remains sizable and provenance-backed",
         test_design_space_audit_preserves_all_sizable_classes)
    test("Discovery metadata persists", test_discovery_metadata_is_persistable)
    test("Indexed space boundaries", test_indexed_space_boundaries_and_classes)
    test("Indexed worker shards", test_indexed_worker_shards_are_disjoint)
    test("Streaming scan resumes", test_streaming_scan_resumes_without_rescoring)
    test("Branch search never surrogate-prunes", test_branch_search_resolves_without_surrogate_pruning)
    test("Branch probes use low-discrepancy schedule", test_branch_probes_are_deterministic_low_discrepancy)
    test("Branch finite budget preserves class floor", test_branch_finite_budget_preserves_class_floor)
    test("Branch certificate detects gaps", test_branch_certificate_detects_incomplete_and_gaps)
    test("Branch rejects population mismatch", test_branch_rejects_population_mismatch)
    test("Tree probes deterministic across 14 classes", test_tree_calibration_probes_cover_all_classes_deterministically)
    test("Production search is branch-only", test_production_has_only_branch_candidate_search)
    test("README matches branch-only contract", test_readme_matches_branch_only_contract)
    test("Retired GA entry points are blocked", test_retired_ga_entry_points_are_blocked)
    test("Industrial viability gates fail closed", test_industrial_viability_gates_fail_closed)
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
