# *De Novo* Catalyst Wet-Lab Validation Protocol

This protocol defines the experimental synthesis, characterization, reactor/cell testing, and checksum-bound data ingestion requirements needed to promote the four *de novo* computational champions from **predicted benchmarks** (`results/evidence/de_novo_predicted_benchmarks.json`) to **measured experimental evidence** in `results/experimental/literature_experimental_dataset.json` and `results/evidence_manifest.json`.

> **Claim boundary:** Until wet-lab measurements are completed and ingested with raw instrument CSVs and SHA-256 digests, all performance metrics for these four champions remain computational predictions (QE DFT + MLIP-NEB + Cantera/PEMFC).

---

## 1. Target *De Novo* Champions & Pathways

| Candidate ID | Material / Composition | Pathway(s) | Required `measurement_type` |
|---|---|---|---|
| `9b83e69315874c5abca7` | `Cu,Ru-doped Ni(fcc211)/CeZrO2` ($d = 10.0\text{ nm}$, $25.0\text{ wt\%}$) | PBR, PFR, Fluidized Bed | `pyrolysis_reactor` |
| `017a44c337ece52f029b` | `Ga85Ni15` Liquid Metal Alloy | MMBCR (Bubble Column) | `pyrolysis_reactor` |
| `eb64b3aa6d3913c21f5d` | `(Ta,Ti)2AlC1(basal_0001)` MAX Phase | NTEC, Electrochemical | `ntec_pair` |
| `dbc15731f8fdb0f99d43` | `W-doped NiCo2O4/N-CNT` Spinel Oxide | PEMFC ORR Cathode | `mea_performance`, `durability`, `hydrogen_impurity` |

---

## 2. Experimental Design & Split-Completeness Requirements

Every `(candidate_id, measurement_type)` group must satisfy the calibration/blinded-holdout schema enforced by `pipeline/evidence/experimental_dataset.py`:

1. **Minimum Independent Experimental Units:**
   - At least **3 independent calibration units** (`split = "calibration"`) and **3 independent holdout units** (`split = "holdout"`).
   - For `ntec_pair`, each unit consists of a matched **NTEC-treated** and **control** run under identical $T$, $P$, and inlet $\text{CH}_4$ flow (minimum 6 pairs total: 3 calibration pairs + 3 holdout pairs).
2. **Blinded Split Assignment:**
   - Split labels (`calibration` vs. `holdout`) must be assigned and timestamped (`split_assigned_at_utc`) prior to data analysis, with `split_assignment_source` recorded.
3. **Mass & Atom Balance Closure:**
   - Pyrolysis runs (`pyrolysis_reactor`, `ntec_pair`) must report `CH4_conversion_fraction`, `H2_selectivity_fraction`, `solid_carbon_yield_fraction`, `carbon_balance_closure_fraction`, `hydrogen_balance_closure_fraction`, `net_energy_kWh_kg_H2`, and `deactivation_fraction_per_h`.
   - Solid carbon product must be characterized (Raman $I_D/I_G$, TGA, SEM/TEM) to distinguish filamentous/graphitic carbon from encapsulating carbon ($C_\delta$).
4. **PEMFC Cathode Testing (`dbc15731f8fdb0f99d43`):**
   - `mea_performance`: Polarization curves ($j$–$V$) at $T_{\text{cell}} = 353\text{ K}$, recording peak power density, cell efficiency, and ohmic/mass-transport losses.
   - `durability`: Accelerated stress test (AST) voltage/power degradation rate ($\mu\text{V h}^{-1}$).
   - `hydrogen_impurity`: Performance retention under trace $\text{CO}$ / $\text{CH}_4$ impurity challenge simulating PSA-purified turquoise $\text{H}_2$.

---

## 3. Data Provenance & Ingestion Workflow

1. Save raw instrument time-series CSVs under `results/experimental/raw/` and log every attempted run (including failed or excluded runs) in `results/experimental/raw/run_log.csv`.
2. Append the structured records (with `raw_source_path` and `raw_source_sha256`) to `results/experimental/literature_experimental_dataset.json`.
3. Validate dataset integrity and regenerate evidence records:
   ```bash
   python3 -m pipeline.evidence.experimental_dataset results/experimental/literature_experimental_dataset.json
   ```
4. Update `results/evidence_manifest.json` and verify hash integrity with `verify_evidence_manifest('results/evidence_manifest.json')`.

---

## 4. Expanding Thermocatalytic $\Delta$-ML QE Calibration (`n > 11`)

In parallel with wet-lab synthesis, expanding the periodic-slab Quantum ESPRESSO (`pw.x`) $\Delta$-ML calibration set beyond the current `n = 11` turquoise-hydrogen runs (versus `n = 48` ORR GC-DFT runs) should prioritize:
- Additional transition-metal-promoted Ni(211) and MAX-phase basal-plane slabs spanning intermediate C–H activation barriers ($0.6\text{–}1.5\text{ eV}$).
- Explicit `pw.x` CI-NEB (`neb.x`) saddle-point convergence checks on the 3 thermocatalytic *de novo* champions (`9b83e69315874c5abca7`, `017a44c337ece52f029b`, `eb64b3aa6d3913c21f5d`), replacing any pending external-solver slots in `results/evidence/de_novo_discovery_validation.json`.
