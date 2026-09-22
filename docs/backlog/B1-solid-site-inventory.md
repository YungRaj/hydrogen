# B1 — Solid-catalyst site inventory without inventing Γ

**Working checklist** (sim after each step; judge solids X on `cat_9`)

- [x] **B1-1** Lock Γ at monolayer (`2.5×10⁻⁹ mol/cm²`). No inventory change. Baseline sim.
- [x] **B1-2** Coarse gridded sweep. Archive: `results/sweeps/inventory_b1/coarse.json` (write-once). Legacy copy: `results/reactor/inventory_sweep_b1_2.json`.
- [x] **B1-2b** Targeted ROI sweep from the coarse file via `propose_roi` (general staged-sweep harness, not a one-off domain). Archive: `results/sweeps/inventory_b1/targeted.json`. Legacy copy: `results/reactor/inventory_sweep_b1_2_roi.json`. Coarse remains readable after targeted.
- [x] **B1-3** Production default `d_p` = 0.13 mm (ROI; Ergun ~0.34 bar at 1000 K config / ~0.40 bar at 1300 K). Loading/dispersion still 1. Phase 2: PFR max 6.82%, fluidized 8.01% (`cat_9` 1300 K). MMBCR unchanged at 98.5%.
- [x] **B1-4** Production defaults `metal_loading=0.5`, `metal_dispersion=0.3` (supported-TCD area fraction; Alves 2021; Sánchez-Bastardo 2021; Gili 2024). Product `a` = 0.15 × `a_geom`. Phase 2: PFR `cat_9` 1300 K **1.02%**, fluidized **1.22%** (was 6.82% / 8.01% at 1×1). MMBCR unchanged at 98.5%.
- [x] **B1-5** Production results report single-pass X, `a`, WHSV, ΔP. Scorecard `results/reactor/phase2_solids_scorecard.json` judges solids on `cat_9` 1300 K (PFR 1.02%, a=4154 m⁻¹, WHSV=900 h⁻¹, ΔP=0.40 bar; fluidized 1.22%, a=3808 m⁻¹, WHSV=180 h⁻¹). H-parked 0.01 eV cats and MMBCR 98.5% are logged, not ranks. Do not start B2 yet.

**Acceptance**

- Keep `site-density` at a physical monolayer (~10¹⁹ m⁻² / 2.5×10⁻⁹ mol/cm²). Do not raise Γ to force Damköhler.
- Solids levers are particle S/V, loading, and dispersion only, inside B4 packed-bed limits.
- MMBCR does not use a site lattice (see B3).
- Document the chosen S/V / loading in `ReactorConfig` with a cite.

**Why**

A metal monolayer is already ~10¹⁹ atoms/m². Extra “sites” in the CMD literature come from dispersion and high Ni loading, then die to encapsulating carbon and plugging [1,2,3].

B1’s X ∝ a (corr ≈ 1) is the monolayer identity on a path that ends at `C_s`, not kinetic closure. Intra-pass site return is [B6](B6-off-site-carbon-nucleation.md). Do not raise Γ.

**Refs (same numbering as ADR 0001)**

1. Alves et al., *Renew. Sustain. Energy Rev.* **2021**, *137*, 110465.
2. Sánchez-Bastardo et al., *Ind. Eng. Chem. Res.* **2021**, *60*, 11855–11881.
3. Gili et al., *ChemCatChem* **2024**.
