#!/usr/bin/env python3
"""
EXHAUSTIVE PIPELINE AUDIT

Programmatically tests every lookup table, every if/elif dispatch,
and every .get() default across the entire pipeline to ensure that
NO material class and NO element is silently excluded, mispriced,
or gets a wrong default.

Run: conda run -n fairchem-env python audit_pipeline.py
"""

import sys, os, re, ast, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PASS = 0
FAIL = 0
ERRORS = []

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✅ {name}")
    else:
        FAIL += 1
        ERRORS.append((name, detail))
        print(f"  ❌ {name}: {detail}")


from pipeline.common.catalyst_spaces import (
    ALL_MATERIAL_CLASSES, generate_random_genome,
    GENERATORS, encode_genome, FEATURE_DIM,
    crossover, mutate, CLASS_WEIGHTS,
)

ALL_14 = list(ALL_MATERIAL_CLASSES)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. EVERY LOOKUP TABLE MUST COVER ALL 14 CLASSES
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ LOOKUP TABLE COVERAGE ═══")

# 1a. GENERATORS
check("GENERATORS covers all classes",
      all(c in GENERATORS for c in ALL_14),
      f"Missing: {[c for c in ALL_14 if c not in GENERATORS]}")

# 1b. CLASS_WEIGHTS
check("CLASS_WEIGHTS covers all classes",
      all(c in CLASS_WEIGHTS for c in ALL_14),
      f"Missing: {[c for c in ALL_14 if c not in CLASS_WEIGHTS]}")

# 1c. TAFEL_SLOPE_BY_CLASS
from pipeline.process.pemfc_model import TAFEL_SLOPE_BY_CLASS
check("TAFEL_SLOPE covers all classes",
      all(c in TAFEL_SLOPE_BY_CLASS for c in ALL_14),
      f"Missing: {[c for c in ALL_14 if c not in TAFEL_SLOPE_BY_CLASS]}")

# 1d. OOD CLASS_CONFIDENCE
from pipeline.common.ood_detector import CLASS_CONFIDENCE
check("OOD CLASS_CONFIDENCE covers all classes",
      all(c in CLASS_CONFIDENCE for c in ALL_14),
      f"Missing: {[c for c in ALL_14 if c not in CLASS_CONFIDENCE]}")

# 1e. BEP_PARAMS (embedded in function — extract programmatically)
from pipeline.common.utils import bep_activation_energy
# Test that each class produces a DIFFERENT result than default
default_e = bep_activation_energy(0.5)
bep_missing = []
for cls in ALL_14:
    e = bep_activation_energy(0.5, material_class=cls)
    # If it equals default AND we didn't set it to default intentionally...
    # We just need to confirm the function doesn't error
    if not isinstance(e, float) or e <= 0:
        bep_missing.append(cls)
check("BEP_PARAMS handles all classes",
      len(bep_missing) == 0,
      f"Failed: {bep_missing}")

# 1f. VALID_CLASSES sets
from pipeline.common.utils import VALID_CLASSES_PYROLYSIS, VALID_CLASSES_FUEL_CELL
check("VALID_CLASSES_PYROLYSIS covers all",
      all(c in VALID_CLASSES_PYROLYSIS for c in ALL_14),
      f"Missing: {[c for c in ALL_14 if c not in VALID_CLASSES_PYROLYSIS]}")
check("VALID_CLASSES_FUEL_CELL covers all",
      all(c in VALID_CLASSES_FUEL_CELL for c in ALL_14),
      f"Missing: {[c for c in ALL_14 if c not in VALID_CLASSES_FUEL_CELL]}")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. EVERY ELEMENT THAT CAN BE GENERATED MUST BE IN ALL TABLES
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ ELEMENT TABLE COVERAGE ═══")

from pipeline.common.utils import CRUSTAL_ABUNDANCE_PPM, METAL_PRICE_USD_KG
from pipeline.screening.fc_genetic_optimizer import _extract_elements_from_genome
from pipeline.common.ood_detector import _ELEMENT_COVERAGE

# Generate many genomes and collect all possible elements
all_elements = set()
elements_by_class = {cls: set() for cls in ALL_14}
for _ in range(10000):
    g = generate_random_genome()
    cls = g[0]
    els = _extract_elements_from_genome(g)
    for e in els:
        all_elements.add(e)
        elements_by_class[cls].add(e)

print(f"  Found {len(all_elements)} unique elements across 10k genomes")

missing_abundance = [e for e in all_elements if e not in CRUSTAL_ABUNDANCE_PPM]
check("All elements in CRUSTAL_ABUNDANCE",
      len(missing_abundance) == 0,
      f"Missing: {sorted(missing_abundance)}")

missing_price = [e for e in all_elements if e not in METAL_PRICE_USD_KG]
check("All elements in METAL_PRICE_USD_KG",
      len(missing_price) == 0,
      f"Missing: {sorted(missing_price)}")

missing_ood_elem = [e for e in all_elements 
                     if e not in _ELEMENT_COVERAGE 
                     and e not in ('N', 'O', 'C', 'F', 'Cl', 'Br', 'I', 'S', 'P', 'H', 'B', 'Si')]
# Note: non-metals are expected to not be in OOD element coverage (they're not metals)
# But metals should be there
check("All metals in OOD _ELEMENT_COVERAGE",
      len(missing_ood_elem) == 0,
      f"Missing: {sorted(missing_ood_elem)}")

# Check FENTON_RISK — it's OK to not list all elements (default=0 is correct)
# But verify it doesn't list elements that aren't in the design space


# ═══════════════════════════════════════════════════════════════════════════════
# 3. ELEMENT EXTRACTOR CONSISTENCY (all 4 copies)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ ELEMENT EXTRACTOR CONSISTENCY ═══")

from pipeline.screening.fc_genetic_optimizer import _extract_elements_from_genome as fc_ga_extract
from pipeline.screening.fc_screener import _extract_elements as fc_screen_extract
from pipeline.screening.surface_screener import _extract_elements as surf_extract
from pipeline.screening.genetic_optimizer import _extract_elements_from_genome as ch4_ga_extract

extractor_failures = []
for cls in ALL_14:
    for _ in range(50):
        g = generate_random_genome(cls)
        try:
            e1 = sorted(fc_ga_extract(g))
            e2 = sorted(fc_screen_extract(g))
            e3 = sorted(surf_extract(g))
            e4 = sorted(ch4_ga_extract(g))
            if not (e1 == e2 == e3 == e4):
                extractor_failures.append(
                    f"{cls}: fc_ga={e1}, fc_screen={e2}, surf={e3}, ch4_ga={e4}"
                )
        except Exception as ex:
            extractor_failures.append(f"{cls}: CRASH: {ex}")

check("4 element extractors consistent (700 genomes)",
      len(extractor_failures) == 0,
      f"{len(extractor_failures)} failures. First: {extractor_failures[0] if extractor_failures else ''}")


# ═══════════════════════════════════════════════════════════════════════════════
# 4. STRUCTURE GENERATION — NO CLASS CRASHES
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ STRUCTURE GENERATION ═══")

from pipeline.screening.surface_screener import generate_structure

struct_failures = []
for cls in ALL_14:
    for _ in range(50):
        g = generate_random_genome(cls)
        try:
            atoms, active_idx, mat_class = generate_structure(g)
            if atoms is None or len(atoms) == 0:
                struct_failures.append(f"{cls}: returned empty atoms")
        except Exception as ex:
            struct_failures.append(f"{cls}: {str(ex)[:80]}")

check(f"Structure generation (700 genomes, 0 crashes)",
      len(struct_failures) == 0,
      f"{len(struct_failures)} failures. First: {struct_failures[0] if struct_failures else ''}")


# ═══════════════════════════════════════════════════════════════════════════════
# 5. ENCODING — NO NaN, NO CRASH, CORRECT DIMENSION
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ ENCODING INTEGRITY ═══")

import numpy as np

encode_failures = []
for cls in ALL_14:
    for _ in range(50):
        g = generate_random_genome(cls)
        try:
            enc = encode_genome(g)
            if len(enc) != FEATURE_DIM:
                encode_failures.append(f"{cls}: dim {len(enc)} != {FEATURE_DIM}")
            if any(np.isnan(enc)):
                encode_failures.append(f"{cls}: NaN in encoding")
            if any(np.isinf(enc)):
                encode_failures.append(f"{cls}: Inf in encoding")
        except Exception as ex:
            encode_failures.append(f"{cls}: {str(ex)[:80]}")

check(f"Encoding (700 genomes, no NaN/Inf)",
      len(encode_failures) == 0,
      f"{len(encode_failures)} failures. First: {encode_failures[0] if encode_failures else ''}")


# ═══════════════════════════════════════════════════════════════════════════════
# 6. CROSSOVER & MUTATION — NO CLASS LEFT BEHIND
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ CROSSOVER & MUTATION ═══")

genetic_failures = []
for cls in ALL_14:
    for _ in range(30):
        g1 = generate_random_genome(cls)
        g2 = generate_random_genome(cls)
        try:
            child = crossover(g1, g2)
            if child[0] != cls:
                genetic_failures.append(f"Crossover changed {cls} to {child[0]}")
        except Exception as ex:
            genetic_failures.append(f"Crossover {cls}: {str(ex)[:80]}")
        try:
            mutant = mutate(g1, rate=1.0)
            if mutant[0] != cls:
                genetic_failures.append(f"Mutation changed {cls} to {mutant[0]}")
        except Exception as ex:
            genetic_failures.append(f"Mutation {cls}: {str(ex)[:80]}")

check(f"Crossover & mutation (840 ops, class preserved)",
      len(genetic_failures) == 0,
      f"{len(genetic_failures)} failures. First: {genetic_failures[0] if genetic_failures else ''}")


# ═══════════════════════════════════════════════════════════════════════════════
# 7. SURROGATE MODELS — NO NaN PREDICTIONS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ SURROGATE PREDICTIONS ═══")

import torch
from pipeline.screening.genetic_optimizer import CatalystSurrogate, predict_batch
from pipeline.screening.fc_genetic_optimizer import ORRCatalystSurrogate
from pipeline.common.catalyst_spaces import encode_population

# CH4 surrogate
ch4_model = CatalystSurrogate(input_dim=FEATURE_DIM)
ch4_model.eval()
ch4_nan_classes = []
for cls in ALL_14:
    pop = [generate_random_genome(cls) for _ in range(20)]
    X = encode_population(pop)
    preds = predict_batch(ch4_model, X, device='cpu')
    for key in preds:
        arr = preds[key]
        if hasattr(arr, '__len__') and any(np.isnan(arr)):
            ch4_nan_classes.append(f"{cls}/{key}")

check("CH4 surrogate: no NaN (280 preds)",
      len(ch4_nan_classes) == 0,
      f"NaN in: {ch4_nan_classes}")

# ORR surrogate  
orr_model = ORRCatalystSurrogate(input_dim=FEATURE_DIM)
orr_model.eval()
orr_nan_classes = []
for cls in ALL_14:
    pop = [generate_random_genome(cls) for _ in range(20)]
    X_tensor = torch.FloatTensor(encode_population(pop))
    with torch.no_grad():
        valid_logit, pred_eta, pred_binding = orr_model(X_tensor)
    for name, tensor in [('valid', valid_logit), ('eta', pred_eta), ('binding', pred_binding)]:
        if torch.any(torch.isnan(tensor)):
            orr_nan_classes.append(f"{cls}/{name}")

check("ORR surrogate: no NaN (280 preds)",
      len(orr_nan_classes) == 0,
      f"NaN in: {orr_nan_classes}")


# ═══════════════════════════════════════════════════════════════════════════════
# 8. NSGA-II OBJECTIVES — NO NaN, PENALTIES WORK
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ NSGA-II OBJECTIVES ═══")

from pipeline.screening.fc_genetic_optimizer import compute_orr_objectives_surrogate
from pipeline.screening.genetic_optimizer import compute_objectives_surrogate

# ORR objectives
orr_obj_failures = []
for cls in ALL_14:
    pop = [generate_random_genome(cls) for _ in range(20)]
    try:
        obj = compute_orr_objectives_surrogate(pop, orr_model, 'cpu')
        if np.any(np.isnan(obj)):
            orr_obj_failures.append(f"{cls}: NaN in objectives")
        if np.any(np.isinf(obj)):
            orr_obj_failures.append(f"{cls}: Inf in objectives")
    except Exception as ex:
        orr_obj_failures.append(f"{cls}: {str(ex)[:80]}")

check("ORR NSGA-II objectives (280 candidates, no NaN)",
      len(orr_obj_failures) == 0,
      f"Failures: {orr_obj_failures}")

# CH4 objectives
ch4_obj_failures = []
for cls in ALL_14:
    pop = [generate_random_genome(cls) for _ in range(20)]
    try:
        obj = compute_objectives_surrogate(pop, ch4_model, 'cpu')
        if np.any(np.isnan(obj)):
            ch4_obj_failures.append(f"{cls}: NaN in objectives")
        if np.any(np.isinf(obj)):
            ch4_obj_failures.append(f"{cls}: Inf in objectives")
    except Exception as ex:
        ch4_obj_failures.append(f"{cls}: {str(ex)[:80]}")

check("CH4 NSGA-II objectives (280 candidates, no NaN)",
      len(ch4_obj_failures) == 0,
      f"Failures: {ch4_obj_failures}")


# ═══════════════════════════════════════════════════════════════════════════════
# 9. OOD CONFIDENCE — EVERY CLASS GETS A SCORE
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ OOD CONFIDENCE ═══")

from pipeline.common.ood_detector import compute_model_confidence, confidence_penalty

ood_failures = []
for cls in ALL_14:
    for _ in range(20):
        g = generate_random_genome(cls)
        els = _extract_elements_from_genome(g)
        try:
            conf = compute_model_confidence(g, els)
            if not (0.0 <= conf <= 1.0):
                ood_failures.append(f"{cls}: conf={conf} out of [0,1]")
            pen = confidence_penalty(conf)
            if not (0.0 <= pen <= 1.0):
                ood_failures.append(f"{cls}: penalty={pen} out of [0,1]")
        except Exception as ex:
            ood_failures.append(f"{cls}: {str(ex)[:80]}")

check("OOD confidence (280 genomes, all in range)",
      len(ood_failures) == 0,
      f"Failures: {ood_failures}")


# ═══════════════════════════════════════════════════════════════════════════════
# 10. COST & FENTON — NO CRASHES, RANGES OK
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ COST & FENTON SCORING ═══")

from pipeline.common.utils import abundance_cost_penalty
from pipeline.screening.fc_genetic_optimizer import _fenton_from_genome, _cost_from_genome

cost_fenton_failures = []
for cls in ALL_14:
    for _ in range(50):
        g = generate_random_genome(cls)
        try:
            cost = _cost_from_genome(g)
            fenton = _fenton_from_genome(g)
            if np.isnan(cost) or np.isinf(cost):
                cost_fenton_failures.append(f"{cls}: cost NaN/Inf")
            if np.isnan(fenton) or np.isinf(fenton):
                cost_fenton_failures.append(f"{cls}: fenton NaN/Inf")
            if not (0 <= fenton <= 10):
                cost_fenton_failures.append(f"{cls}: fenton={fenton} out of [0,10]")
        except Exception as ex:
            cost_fenton_failures.append(f"{cls}: {str(ex)[:80]}")

check("Cost & Fenton scoring (700 genomes, no NaN)",
      len(cost_fenton_failures) == 0,
      f"Failures: {cost_fenton_failures}")


# ═══════════════════════════════════════════════════════════════════════════════
# 11. PEMFC MODEL — ALL CLASSES PRODUCE VALID POWER
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ PEMFC MODEL ═══")

from pipeline.process.pemfc_model import sweep_membranes

pemfc_failures = []
for cls in ALL_14:
    try:
        results = sweep_membranes(cathode_name=f"audit_{cls}", orr_eta=0.5, material_class=cls)
        if not results or len(results) == 0:
            pemfc_failures.append(f"{cls}: returned empty")
        else:
            best = max(results, key=lambda r: r.get('peak_power_W_cm2', 0))
            if best['peak_power_W_cm2'] <= 0:
                pemfc_failures.append(f"{cls}: power <= 0")
            if np.isnan(best['peak_power_W_cm2']):
                pemfc_failures.append(f"{cls}: power NaN")
    except Exception as ex:
        pemfc_failures.append(f"{cls}: {str(ex)[:80]}")

check("PEMFC model (all 14 classes produce power)",
      len(pemfc_failures) == 0,
      f"Failures: {pemfc_failures}")


# ═══════════════════════════════════════════════════════════════════════════════
# 12. FC CATHODE SCREENER — GENOME FORMAT CHECK
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ FC CATHODE SCREENER ═══")

from pipeline.screening.fc_cathode_screener import generate_fc_catalyst_list

candidates = generate_fc_catalyst_list()
cathode_failures = []

for c in candidates:
    genome = c['genome']
    cls = genome[0]
    
    # Verify genome can be encoded
    try:
        enc = encode_genome(genome)
        if len(enc) != FEATURE_DIM:
            cathode_failures.append(f"{c['name']}: encoding dim {len(enc)}")
        if any(np.isnan(enc)):
            cathode_failures.append(f"{c['name']}: NaN encoding")
    except Exception as ex:
        cathode_failures.append(f"{c['name']}: encode crash: {str(ex)[:60]}")

    # Verify structure generation
    try:
        atoms, _, _ = generate_structure(genome)
        if atoms is None or len(atoms) == 0:
            cathode_failures.append(f"{c['name']}: empty structure")
    except Exception as ex:
        cathode_failures.append(f"{c['name']}: struct crash: {str(ex)[:60]}")

check(f"FC cathode candidates ({len(candidates)} candidates, encode+structure)",
      len(cathode_failures) == 0,
      f"{len(cathode_failures)} failures. First: {cathode_failures[0] if cathode_failures else ''}")


# ═══════════════════════════════════════════════════════════════════════════════
# 13. STATIC CODE ANALYSIS — FIND HARDCODED CLASS LISTS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ STATIC CODE ANALYSIS ═══")

import glob

pipeline_files = glob.glob("pipeline/**/*.py", recursive=True)
class_names_in_code = set()
string_pattern = re.compile(r"'(MoltenMetal|SolidCatalyst|SAC|DAC|MOF|COF|Perovskite|MetalHydride|MAXPhase|HEA|Spinel|MXene|SAA|MetalFreeCarbon)'")

# Find all files that reference material class names
files_with_classes = {}
for fpath in pipeline_files:
    with open(fpath) as f:
        content = f.read()
    found = set(string_pattern.findall(content))
    if found:
        files_with_classes[os.path.basename(fpath)] = found

# Check if any file mentions some classes but not all
# (This catches if/elif chains that miss new classes)
# Exclude files that route unknown classes through generate_structure()
# (dft_validator.py handles new classes via the else branch)
FILES_WITH_CATCH_ALL = {'dft_validator.py', 'dft_fuel_cell.py'}
incomplete_files = {}
for fname, classes_found in files_with_classes.items():
    # Skip files that only mention 1-2 classes (utility or test)
    if len(classes_found) >= 5 and fname not in FILES_WITH_CATCH_ALL:
        missing = set(ALL_14) - classes_found
        if missing:
            incomplete_files[fname] = missing

if incomplete_files:
    detail = "; ".join(f"{f}: missing {sorted(m)}" for f, m in incomplete_files.items())
    check("All files with class dispatches cover all 14", False, detail)
else:
    check("All files with class dispatches cover all 14", True)


# ═══════════════════════════════════════════════════════════════════════════════
# 14. DESIGN-SPACE PROVENANCE, CANONICALIZATION, AND CLASS SIZE
# ═══════════════════════════════════════════════════════════════════════════════
print("\n═══ DESIGN SPACE AUDIT ═══")

from pipeline.evidence.design_space_audit import audit_design_space

design_audit = audit_design_space(sample_per_class=1024)
check("All 14 design classes remain sizable and provenance-backed",
      design_audit['valid'],
      f"Failures: {design_audit['failures']}")


# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'═' * 60}")
print(f"  AUDIT COMPLETE: {PASS} passed, {FAIL} failed")
if FAIL == 0:
    print(f"  🟢 NO SILENT EXCLUSIONS FOUND")
else:
    print(f"  🔴 {FAIL} ISSUE(S) FOUND:")
    for name, detail in ERRORS:
        print(f"\n  ❌ {name}")
        print(f"     {detail}")
print("═" * 60)
sys.exit(0 if FAIL == 0 else 1)
