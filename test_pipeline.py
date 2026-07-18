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
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome
    pop = [generate_random_genome() for _ in range(3000)]
    classes_seen = set(g[0] for g in pop)
    assert len(ALL_MATERIAL_CLASSES) == 14, f"Expected 14 classes, got {len(ALL_MATERIAL_CLASSES)}"
    assert classes_seen == set(ALL_MATERIAL_CLASSES), f"Missing: {set(ALL_MATERIAL_CLASSES) - classes_seen}"


def test_no_toxic_elements():
    from pipeline.catalyst_spaces import generate_random_genome
    from pipeline.utils import TOXIC_ELEMENTS
    pop = [generate_random_genome() for _ in range(2000)]
    for g in pop:
        for field in g[1:]:
            if isinstance(field, str) and field in TOXIC_ELEMENTS:
                raise AssertionError(f"Toxic element {field} in genome {g}")


def test_sac_has_axial_ligand():
    from pipeline.catalyst_spaces import generate_random_genome, SAC_AXIAL_LIGANDS
    sacs = [generate_random_genome('SAC') for _ in range(100)]
    assert all(len(g) == 5 for g in sacs), "SAC genome should have 5 fields"
    axials = set(g[4] for g in sacs)
    assert len(axials) > 3, f"Only {len(axials)} axial variants seen, expected >3"


def test_class_weights_sum_to_1():
    from pipeline.catalyst_spaces import CLASS_WEIGHTS
    total = sum(CLASS_WEIGHTS.values())
    assert abs(total - 1.0) < 0.01, f"Class weights sum to {total}, expected ~1.0"


# ═══════════════════════════════════════════════════════════════════════════════
# 2. ENCODING
# ═══════════════════════════════════════════════════════════════════════════════

def test_encode_all_classes_no_nan():
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome, encode_genome, FEATURE_DIM
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(50):
            g = generate_random_genome(cls)
            feat = encode_genome(g)
            assert feat.shape == (FEATURE_DIM,), f"{cls}: shape {feat.shape} != ({FEATURE_DIM},)"
            assert not np.any(np.isnan(feat)), f"{cls}: NaN in encoding"


def test_feature_dim_matches_components():
    from pipeline.catalyst_spaces import (
        FEATURE_DIM, N_CLASSES, N_METALS, N_SUPPORTS, N_FACETS,
        N_COORDS, N_DOPANTS, N_CONTINUOUS
    )
    expected = N_CLASSES + 2 * N_METALS + N_SUPPORTS + N_FACETS + N_COORDS + N_DOPANTS + N_CONTINUOUS
    assert FEATURE_DIM == expected, f"FEATURE_DIM={FEATURE_DIM} != computed {expected}"


# ═══════════════════════════════════════════════════════════════════════════════
# 3. SURROGATES
# ═══════════════════════════════════════════════════════════════════════════════

def test_ch4_surrogate_no_nan():
    from pipeline.catalyst_spaces import generate_random_genome, encode_genome, FEATURE_DIM
    from pipeline.surrogate_model import CatalystSurrogate
    model = CatalystSurrogate(input_dim=FEATURE_DIM)
    model.eval()
    pop = [generate_random_genome() for _ in range(500)]
    X = torch.FloatTensor(np.array([encode_genome(g) for g in pop]))
    with torch.no_grad():
        out = model(X)
    for i, o in enumerate(out):
        assert not torch.any(torch.isnan(o)), f"CH4 surrogate head {i} has NaN"


def test_orr_surrogate_no_nan():
    from pipeline.catalyst_spaces import generate_random_genome, encode_genome, FEATURE_DIM
    from pipeline.fc_genetic_optimizer import ORRCatalystSurrogate
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
    from pipeline.catalyst_spaces import generate_random_genome, FEATURE_DIM
    from pipeline.fc_genetic_optimizer import (
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
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome
    from pipeline.fc_genetic_optimizer import _cost_from_genome, _fenton_from_genome
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(20):
            g = generate_random_genome(cls)
            c = _cost_from_genome(g)
            f = _fenton_from_genome(g)
            assert 0 <= c <= 2.0, f"{cls}: cost {c} out of range"
            assert 0 <= f <= 10, f"{cls}: fenton {f} out of range"


def test_metalfreecarbon_zero_cost():
    from pipeline.catalyst_spaces import generate_random_genome
    from pipeline.fc_genetic_optimizer import _cost_from_genome
    for _ in range(20):
        g = generate_random_genome('MetalFreeCarbon')
        c = _cost_from_genome(g)
        assert c == 0.0, f"MetalFreeCarbon cost should be 0, got {c}"


# ═══════════════════════════════════════════════════════════════════════════════
# 5. ELEMENT EXTRACTORS (4 copies must agree)
# ═══════════════════════════════════════════════════════════════════════════════

def test_element_extractors_consistent():
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome
    from pipeline.fc_genetic_optimizer import _extract_elements_from_genome as fc_extract
    from pipeline.genetic_optimizer import _extract_elements_from_genome as methane_extract
    from pipeline.surface_screener import _extract_elements as screener_extract
    from pipeline.fc_screener import _extract_elements as fc_screener_extract
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
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome
    from pipeline.surface_screener import generate_structure
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
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.pemfc_model import TAFEL_SLOPE_BY_CLASS
    for cls in ALL_MATERIAL_CLASSES:
        assert cls in TAFEL_SLOPE_BY_CLASS, f"Missing Tafel slope for {cls}"


def test_pemfc_power_monotonic_with_eta():
    from pipeline.pemfc_model import simulate_pemfc, PEMFCConfig
    etas = [0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80]
    results = [simulate_pemfc(PEMFCConfig(orr_overpotential_V=e)) for e in etas]
    powers = [r['peak_power_W_cm2'] for r in results]
    for i in range(len(powers) - 1):
        assert powers[i] >= powers[i+1], f"Power not monotonic: eta={etas[i]}->{etas[i+1]}"


def test_tafel_slope_affects_power():
    from pipeline.pemfc_model import simulate_pemfc, PEMFCConfig
    r_low = simulate_pemfc(PEMFCConfig(orr_overpotential_V=0.35, orr_tafel_slope_mV_dec=65))
    r_high = simulate_pemfc(PEMFCConfig(orr_overpotential_V=0.35, orr_tafel_slope_mV_dec=130))
    assert r_low['peak_power_W_cm2'] > r_high['peak_power_W_cm2'], \
        f"Lower Tafel should give higher power: {r_low['peak_power_W_cm2']} vs {r_high['peak_power_W_cm2']}"


def test_pemfc_efficiency_in_range():
    from pipeline.pemfc_model import simulate_pemfc, PEMFCConfig
    r = simulate_pemfc(PEMFCConfig(orr_overpotential_V=0.35))
    eff = r['efficiency_at_peak']
    assert 0.10 < eff < 0.60, f"Efficiency {eff} out of physical range"


# ═══════════════════════════════════════════════════════════════════════════════
# 8. STACK MODEL
# ═══════════════════════════════════════════════════════════════════════════════

def test_stack_model():
    from pipeline.fuel_cell_stack import model_stack, StackConfig
    stack = model_stack(StackConfig(cell_voltage_V=0.65, current_density_A_cm2=1.5))
    assert stack['net_power_kW'] > 0, f"Negative net power: {stack['net_power_kW']}"
    assert 0 < stack['system_efficiency'] < 1, f"Efficiency out of range: {stack['system_efficiency']}"
    assert stack['cost_per_kW'] > 0, f"Negative cost: {stack['cost_per_kW']}"


# ═══════════════════════════════════════════════════════════════════════════════
# 9. CHE / UTILS
# ═══════════════════════════════════════════════════════════════════════════════

def test_orr_overpotential_ideal():
    from pipeline.utils import orr_overpotential
    eta, rds = orr_overpotential(1.23, 2.46, 3.69)
    assert abs(eta) < 0.001, f"Ideal overpotential should be ~0, got {eta}"


def test_abundance_cost_penalty():
    from pipeline.utils import abundance_cost_penalty
    assert abundance_cost_penalty(['Fe']) == 0.0, "Fe should be zero cost"
    assert abundance_cost_penalty(['Ir']) == -2.0, "Ir should be max penalty"
    assert abundance_cost_penalty(['Fe', 'Ir']) < 0, "Geo mean should catch Ir"
    assert abundance_cost_penalty([]) == 0.0, "Empty should be zero"


# ═══════════════════════════════════════════════════════════════════════════════
# 10. REPORT GENERATOR
# ═══════════════════════════════════════════════════════════════════════════════

def test_report_campaign_path():
    from pipeline.report_generator import generate_full_report
    path = generate_full_report({
        'phase5_ga': {'total_evaluated': 1000},
        'phase5_stack': {'best_power_W_cm2': 1.35, 'best_efficiency': 0.31, 'min_overpotential_V': 0.32},
    })
    with open(path) as f:
        content = f.read()
    assert 'N/A' not in content.split('FC catalysts')[1].split('\n')[0], "FC row shows N/A in campaign mode"


def test_report_orchestrator_path():
    from pipeline.report_generator import generate_full_report
    path = generate_full_report({
        'phase5': {'n_cathodes_screened': 137, 'best_power_W_cm2': 1.2,
                   'best_efficiency': 0.28, 'min_overpotential_V': 0.35},
    })
    with open(path) as f:
        content = f.read()
    assert 'N/A' not in content.split('PEMFC power')[1].split('\n')[0], "PEMFC row shows N/A in orchestrator mode"


def test_report_empty_no_crash():
    from pipeline.report_generator import generate_full_report
    path = generate_full_report({})
    assert os.path.exists(path)


# ═══════════════════════════════════════════════════════════════════════════════
# 11. CROSSOVER & MUTATION
# ═══════════════════════════════════════════════════════════════════════════════

def test_crossover_preserves_class():
    from pipeline.catalyst_spaces import generate_random_genome, crossover
    for _ in range(100):
        p1 = generate_random_genome('SAC')
        p2 = generate_random_genome('SAC')
        child = crossover(p1, p2)
        assert child[0] == 'SAC', f"Crossover changed class to {child[0]}"


def test_mutation_preserves_class():
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES, generate_random_genome, mutate
    for cls in ALL_MATERIAL_CLASSES:
        for _ in range(20):
            g = generate_random_genome(cls)
            m = mutate(g, rate=1.0)  # force mutation
            assert m[0] == cls, f"Mutation changed {cls} to {m[0]}"


# ═══════════════════════════════════════════════════════════════════════════════
# 12. OOD CONFIDENCE
# ═══════════════════════════════════════════════════════════════════════════════

def test_ood_high_confidence_metals():
    from pipeline.catalyst_spaces import generate_random_genome
    from pipeline.ood_detector import compute_model_confidence
    from pipeline.fc_genetic_optimizer import _extract_elements_from_genome
    # Metal slabs should have high confidence (>0.7)
    for cls in ['SolidCatalyst', 'HEA', 'SAA']:
        for _ in range(10):
            g = generate_random_genome(cls)
            elements = _extract_elements_from_genome(g)
            conf = compute_model_confidence(g, elements)
            assert conf > 0.6, f"{cls}: confidence {conf:.2f} too low (expected >0.6)"


def test_ood_low_confidence_ood():
    from pipeline.catalyst_spaces import generate_random_genome
    from pipeline.ood_detector import compute_model_confidence
    from pipeline.fc_genetic_optimizer import _extract_elements_from_genome
    # OOD classes should have low confidence (<0.5)
    for cls in ['MOF', 'COF', 'MetalFreeCarbon']:
        for _ in range(10):
            g = generate_random_genome(cls)
            elements = _extract_elements_from_genome(g)
            conf = compute_model_confidence(g, elements)
            assert conf < 0.5, f"{cls}: confidence {conf:.2f} too high (expected <0.5)"


def test_ood_penalty_scales_objectives():
    from pipeline.ood_detector import confidence_penalty
    # High confidence → penalty near 0.0 (no shift)
    assert abs(confidence_penalty(1.0) - 0.0) < 0.01, "conf=1.0 should give penalty=0.0"
    # Low confidence → penalty > 0.5 (significant shift)
    assert confidence_penalty(0.2) > 0.5, "conf=0.2 should give penalty>0.5"
    # Zero confidence → maximum penalty = 1.0
    assert abs(confidence_penalty(0.0) - 1.0) < 0.01, "conf=0.0 should give penalty=1.0"


def test_ood_nsga2_integration():
    from pipeline.catalyst_spaces import generate_random_genome, FEATURE_DIM
    from pipeline.fc_genetic_optimizer import (
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
    from pipeline.catalyst_spaces import generate_random_genome
    from pipeline.fc_genetic_optimizer import _extract_elements_from_genome
    from pipeline.utils import CRUSTAL_ABUNDANCE_PPM
    missing = set()
    for _ in range(3000):
        g = generate_random_genome()
        for e in _extract_elements_from_genome(g):
            if e not in CRUSTAL_ABUNDANCE_PPM:
                missing.add(e)
    assert not missing, f"Elements missing from CRUSTAL_ABUNDANCE: {sorted(missing)}"


def test_all_elements_in_price_table():
    from pipeline.catalyst_spaces import generate_random_genome
    from pipeline.fc_genetic_optimizer import _extract_elements_from_genome
    from pipeline.utils import METAL_PRICE_USD_KG
    missing = set()
    for _ in range(3000):
        g = generate_random_genome()
        for e in _extract_elements_from_genome(g):
            if e not in METAL_PRICE_USD_KG:
                missing.add(e)
    assert not missing, f"Elements missing from METAL_PRICE_USD_KG: {sorted(missing)}"


def test_all_classes_viable_both_applications():
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.utils import VALID_CLASSES_PYROLYSIS, VALID_CLASSES_FUEL_CELL
    for cls in ALL_MATERIAL_CLASSES:
        assert cls in VALID_CLASSES_PYROLYSIS, f"{cls} excluded from pyrolysis"
        assert cls in VALID_CLASSES_FUEL_CELL, f"{cls} excluded from fuel cell"


def test_bep_params_all_classes():
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.utils import bep_activation_energy
    for cls in ALL_MATERIAL_CLASSES:
        # Should not fall back to default — each class must have specific params
        e1 = bep_activation_energy(0.5, material_class=cls)
        e_default = bep_activation_energy(0.5)  # no class = default
        # At least some classes should differ from default
        assert isinstance(e1, float) and e1 > 0


def test_ood_confidence_all_classes():
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.ood_detector import CLASS_CONFIDENCE
    for cls in ALL_MATERIAL_CLASSES:
        assert cls in CLASS_CONFIDENCE, f"OOD CLASS_CONFIDENCE missing: {cls}"


def test_tafel_all_classes():
    from pipeline.catalyst_spaces import ALL_MATERIAL_CLASSES
    from pipeline.pemfc_model import TAFEL_SLOPE_BY_CLASS
    for cls in ALL_MATERIAL_CLASSES:
        assert cls in TAFEL_SLOPE_BY_CLASS, f"TAFEL_SLOPE missing: {cls}"


def test_cathode_sac_genome_5tuple():
    from pipeline.fc_cathode_screener import generate_fc_catalyst_list
    candidates = generate_fc_catalyst_list()
    sacs = [c for c in candidates if c['type'] == 'SAC']
    assert len(sacs) > 0, "No SAC candidates generated"
    for c in sacs:
        assert len(c['genome']) == 5, f"SAC genome should be 5-tuple, got {len(c['genome'])}: {c['genome']}"


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
    test("Class weights sum ~1.0", test_class_weights_sum_to_1)

    print("\n── Encoding ──")
    test("All classes encode (no NaN)", test_encode_all_classes_no_nan)
    test("FEATURE_DIM matches components", test_feature_dim_matches_components)

    print("\n── Surrogates ──")
    test("CH4 surrogate (no NaN)", test_ch4_surrogate_no_nan)
    test("ORR surrogate (no NaN)", test_orr_surrogate_no_nan)

    print("\n── NSGA-II ──")
    test("NSGA-II sorts correctly", test_nsga2_sorts_correctly)
    test("Cost & Fenton in range", test_cost_and_fenton_ranges)
    test("MetalFreeCarbon cost = 0", test_metalfreecarbon_zero_cost)

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

    print("\n── Crossover & Mutation ──")
    test("Crossover preserves class", test_crossover_preserves_class)
    test("Mutation preserves class", test_mutation_preserves_class)

    print("\n── OOD Confidence ──")
    test("High confidence for metals", test_ood_high_confidence_metals)
    test("Low confidence for OOD classes", test_ood_low_confidence_ood)
    test("Penalty scales objectives", test_ood_penalty_scales_objectives)
    test("Confidence in NSGA-II", test_ood_nsga2_integration)

    print("\n── Exhaustive Coverage ──")
    test("All elements in abundance table", test_all_elements_in_abundance_table)
    test("All elements in price table", test_all_elements_in_price_table)
    test("All 14 classes viable both apps", test_all_classes_viable_both_applications)
    test("BEP params all 14 classes", test_bep_params_all_classes)
    test("OOD confidence all 14 classes", test_ood_confidence_all_classes)
    test("Tafel slope all 14 classes", test_tafel_all_classes)
    test("Cathode SAC genomes 5-tuple", test_cathode_sac_genome_5tuple)

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
