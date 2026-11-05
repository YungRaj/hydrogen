# Turquoise Hydrogen → Fuel Cell: Pipeline Report

**Generated**: 2026-11-05 01:57:18
**Pipeline Version**: 2.0.0

> **Evidence status:** screening values are hypotheses. Mock results and unconverged calculations are excluded from validated claims; catalyst performance requires experimental confirmation.

## Executive Summary

| Metric | Value |
|--------|-------|
| Catalysts evaluated (Fairchem) | 84 |
| Valid candidates | 58 |
| Pareto front size | 1474 |
| Best eligible screening barrier | 0.015 eV |
| Best coking resistance | 8.359 eV |
| Reactor simulations | 20 |
| Solids judge CH₄ conversion (ni_np_lit PFR) | 26.34% |
| Top candidate solids CH₄ conversion (0021a267d0b5a6e984a9) | 98.50% |
| MMBCR max CH₄ conversion | 24.17% (not a solids rank) |
| NTEC max CH₄ conversion (coupled MOOSE+Cantera) | 1.40% |
| Electrochemical max CH₄ conversion (coupled MOOSE+Cantera) | 0.67% |
| FC catalysts evaluated | 104 |
| Best PEMFC power density | 1.348 W/cm² |
| Best PEMFC efficiency | 35.9% |
| Min ORR overpotential | 0.297 V |

## Phase 1: Catalyst Screening & Genetic Optimization

### Screening Database: `pyrolysis_benchmark`

- **Total entries**: 18
- **Valid candidates**: 18

**Distribution by Material Class:**

| Material Class | Count | Valid |
|----------------|-------|-------|
| MoltenMetal | 18 | 18 |


**Top 10 Catalysts by Activation Barrier:**

| Rank | Genome | E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | Coking Index | dE_H (eV) |
|------|--------|------------|------------------------|-----------------------|-------------|-----------|
| 1 | Liquid Bi (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=ac356494d47705c70446] | 0.1193 | 2.2715 | 3.4495 | nan | 0.2133188232779467 |
| 2 | Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | 0.3812 | 1.5941 | 2.7721 | nan | -0.7748952880501726 |
| 3 | Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | 0.4039 | 1.5278 | 2.7057 | nan | -0.3025980010628622 |
| 4 | Bi89Ga11 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=0e7e12a7c4a5ff0d18a9] | 0.5202 | 2.0049 | 3.1829 | nan | -0.2237430587410962 |
| 5 | Bi83Ni17 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=87880819e7c1e902b8c0] | 0.7906 | 1.6551 | 2.8330 | nan | -0.0687772765755689 |
| 6 | Bi55Zn45 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=df0cc4f2bd66c68f59a0] | 1.0082 | 1.6561 | 2.8341 | nan | 0.1641103729605525 |
| 7 | Liquid Sn (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=fc57e59a4f6028ccaaf7] | 1.0318 | 2.8429 | 4.0208 | nan | 0.2157167419791328 |
| 8 | Liquid Sn (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=80bc69f4855453f54b71] | 1.0533 | 3.0542 | 4.2322 | nan | 0.2357795700430976 |
| 9 | Sn67Ni33 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=34d68c2b96bb438d7d84] | 1.0889 | 2.4706 | 3.6486 | nan | 0.2007752880454241 |
| 10 | Liquid Sn (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=45fdf4365dc1b3f3fd91] | 1.0984 | 2.9420 | 4.1200 | nan | 0.2178939804434883 |

## Phase 2: Reactor-Scale Simulation

Solids judged on **ni_np_lit** (single-pass X, `a`, WHSV, ΔP). MMBCR X_eq and 0.01 eV H-parked cats are not ranks.

| Reactor | Catalyst | T (K) | Single-pass X | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |
|---------|----------|-------|---------------|---------|------------|----------|
| PFR | ni_np_lit | 973.15 | 26.34% | 4154 | 900 | 0.33 |
| PBR | ni_np_lit | 973.15 | 23.03% | 4154 | 1272 | 0.50 |
| Fluidized | ni_np_lit | 973.15 | 18.79% | 3808 | 400 | 0.33 |

Top candidate solids X = 98.50% on **0021a267d0b5a6e984a9** (PBR @ 1300.0 K)

### Solids runs (single-pass)

| Catalyst | Reactor | T (K) | Single-pass X | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |
|----------|---------|-------|---------------|---------|------------|----------|
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=0021a267d0b5a6e984a9] | PBR | 1300.0 | 98.50% | 4154 | 2956 | 0.91 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=f78c667155715fb12ba3] | PBR | 1300.0 | 98.00% | 4154 | 2844 | 0.91 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=d4403c08f769fc70a3a4] | PFR | 1300.0 | 97.11% | 4154 | 900 | 0.40 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PBR | 1300.0 | 97.09% | 4154 | 2783 | 0.91 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=f78c667155715fb12ba3] | PFR | 1300.0 | 97.06% | 4154 | 900 | 0.40 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PFR | 1300.0 | 97.05% | 4154 | 900 | 0.40 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=0021a267d0b5a6e984a9] | PFR | 1300.0 | 96.98% | 4154 | 900 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PBR | 1300.0 | 96.58% | 4154 | 2956 | 0.91 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=d4403c08f769fc70a3a4] | PBR | 1300.0 | 95.44% | 4154 | 2665 | 0.91 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PFR | 1300.0 | 94.01% | 4154 | 900 | 0.40 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PFR | 1100.0 | 72.54% | 4154 | 900 | 0.36 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PBR | 1100.0 | 60.18% | 4154 | 1696 | 0.84 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=0021a267d0b5a6e984a9] | PFR | 1100.0 | 53.48% | 4154 | 900 | 0.36 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=f78c667155715fb12ba3] | PFR | 1100.0 | 48.92% | 4154 | 900 | 0.36 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=d4403c08f769fc70a3a4] | Fluidized | 1300.0 | 48.49% | 3808 | 400 | 0.40 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=f78c667155715fb12ba3] | Fluidized | 1300.0 | 48.45% | 3808 | 400 | 0.40 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | Fluidized | 1300.0 | 48.45% | 3808 | 400 | 0.40 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=0021a267d0b5a6e984a9] | PBR | 1100.0 | 48.42% | 4154 | 1572 | 0.73 |
| Multi-Principal-Element High-Entropy Alloy Catalyst (Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018) [id=0021a267d0b5a6e984a9] | Fluidized | 1300.0 | 48.41% | 3808 | 400 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PFR | 1100.0 | 44.65% | 4154 | 900 | 0.36 |

### Liquid-Metal, NTEC, and Electrochemical Runs

| Catalyst | Reactor | T (K) | CH₄ Conversion | H₂ Selectivity | Evidence Tier |
|----------|---------|-------|----------------|----------------|---------------|
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=54757bb17c1b638462d9] | MMBCR | 1300.0 | 24.17% | 81.32% | diagnostic_screening_template |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=54757bb17c1b638462d9] | MMBCR | 1100.0 | 9.58% | 91.66% | diagnostic_screening_template |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=54757bb17c1b638462d9] | MMBCR | 973.15 | 3.65% | 96.65% | diagnostic_screening_template |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=c4ab2249c156880efdbe] | MMBCR | 1300.0 | 3.35% | 96.92% | diagnostic_screening_template |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=54757bb17c1b638462d9] | MMBCR | 923.15 | 2.21% | 97.94% | diagnostic_screening_template |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=54757bb17c1b638462d9] | MMBCR | 900.0 | 1.72% | 98.39% | diagnostic_screening_template |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=54757bb17c1b638462d9] | NTEC | 473.15 | 1.40% | 100.00% | uncalibrated_multiphysics_screening |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=c4ab2249c156880efdbe] | NTEC | 473.15 | 1.34% | 100.00% | uncalibrated_multiphysics_screening |
| Atomically Dispersed Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015; Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | NTEC | 473.15 | 1.16% | 100.00% | uncalibrated_multiphysics_screening |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | NTEC | 423.15 | 1.08% | 100.00% | uncalibrated_multiphysics_screening |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=c4ab2249c156880efdbe] | MMBCR | 1100.0 | 0.84% | 0.00% | diagnostic_screening_template |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | NTEC | 473.15 | 0.79% | 100.00% | uncalibrated_multiphysics_screening |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Electrochemical | 353.15 | 0.67% | 100.00% | uncalibrated_multiphysics_screening |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=00014a7eeba8246b1b9c] | Electrochemical | 353.15 | 0.60% | 100.00% | uncalibrated_multiphysics_screening |
| Mn+1AXn Nanolaminate MAX Phase Catalyst (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | Electrochemical | 353.15 | 0.60% | 100.00% | uncalibrated_multiphysics_screening |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=00014a7eeba8246b1b9c] | NTEC | 423.15 | 0.57% | 100.00% | uncalibrated_multiphysics_screening |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=54757bb17c1b638462d9] | Electrochemical | 823.15 | 0.48% | 100.00% | uncalibrated_multiphysics_screening |
| Liquid Metal Alloy Bubble-Column Catalyst (Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023) [id=c4ab2249c156880efdbe] | Electrochemical | 823.15 | 0.48% | 100.00% | uncalibrated_multiphysics_screening |
| Nd-N4_pyridine-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=8609ad5ba35faea2ac57] | Electrochemical | 823.15 | 0.48% | 100.00% | uncalibrated_multiphysics_screening |
| Mn+1AXn Nanolaminate MAX Phase Catalyst (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | Electrochemical | 823.15 | 0.48% | 100.00% | uncalibrated_multiphysics_screening |

Excluded 129 non-run reactor files (sweeps, mocks, not_applicable, validation_required, failed) from performance claims.

## Phase 3: DFT Validation (Quantum ESPRESSO)

| Catalyst | Material Class | Energy (eV) | Converged | Max Force |
|----------|----------------|-------------|-----------|-----------|
| Br,Mg-doped Ce(bcc211)/Graphene, strain=+0.07% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=000971115da8ddd8cf8f] | SolidCatalyst | -20312.388616231547 | True | 0.0 |
| FeZnNbIrSc-hcp(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=0021a267d0b5a6e984a9] | HEA | -14552.039192419978 | True | 0.0 |
| ZnAl-N3BN3/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=0433514342912e4613e1] | DAC | -6808.752593966458 | True | 0.0 |
| Ti0.10Te0.90 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=54757bb17c1b638462d9] | MoltenMetal | -2701.304520951499 | True | 0.0 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | -612.2664572547036 | True | 0.0 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | -16117.69086882 | True | 0.0 |
| Pd1/In(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=8796cf31edc4797d83ec] | SAA | -10830.340859828559 | True | 0.0 |
| Zn-N4_pyridine-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=b1bd7298e14db77a573c] | SAC | -6808.752593966458 | True | 0.0 |
| Ga90Ti10 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | -6971.842105 | True | 0.0 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | -5355.441756157839 | True | 0.0 |
| FeZnZrAgIr-bcc(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=d4403c08f769fc70a3a4] | HEA | -14011.380335365946 | True | 0.0 |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | MAXPhase | -4337.103472226402 | True | 0.0 |
| FeZnZrAgIr-(0.2, 0.2, 0.2, 0.2, 0.2)(hcp) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=f78c667155715fb12ba3] | HEA | -6971.842105 | True | 0.0 |
| BaSbO3(B_vacancy) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=06d7172dc944ea7d38c7] | Perovskite | -14951.328010180016 | True | None |
| VCuZnAlGeCe-fcc_bcc_dual(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=0cc10f36dedb70353ccb] | HEA | -12621.51452530112 | True | None |
| Pt1/Ga(211) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=101159f93bd9c6d9f6f8] | SAA | -13297.8650716104 | True | None |
| CoRuIrSbLaCe-amorphous(100) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=1371fd2696264ba0a3d2] | HEA | -13782.649160139848 | True | 0.0 |
| Ce0.95Nb0.05SnO3(none) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=2e7fa45b0e15035b8384] | Perovskite | -13887.623749237322 | True | None |
| Sb-O4-Phosphorene/O2 Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=2f38d51ba986c94ba69a] | SAC | -3589.3135179342676 | True | None |
| Cu,Sc-doped Co(fcc100)/Graphite, strain=-0.06% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=36706ae55f59ec4f8504] | SolidCatalyst | -15082.113041997116 | True | None |
| Cu-C4-Phosphorene/pyridine Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=383a5786185703639bcd] | SAC | -3429.4082782699315 | True | None |
| Y-doped FeNi2O4/carbon_black (nanosheet) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=3a597e8baaee3ebe72d3] | Spinel | -6674.94648955734 | True | None |
| Pd1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=3ab787b3c8004744ae8e] | SAA | -14364.996728013557 | True | None |
| Co-Imine-N2P2 (6.0) COF Network (Furukawa et al., Science 2013) [id=3cafa1dd27033ac418d3] | COF | -4769.527949038315 | True | 0.0 |
| (W,Fe)4CN3Sx 2D MXene (Naguib et al., Adv. Mater. 2011 / Seh 2016) [id=3fb051c06f48d5784c5f] | MXene | -7162.44243844502 | True | None |
| Ce-N3B-N-mesoporous_carbon/O2 Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=4605d672ccdbce8db71c] | SAC | -7154.784889044383 | True | None |
| Gd0.85Mn0.15ZnO3(none) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=46aa11fb6311a1ce2249] | Perovskite | -19010.92005767173 | True | None |
| Ba0.5Al0.5RhO3(B_vacancy) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=4f112dec256800b8f5d7] | Perovskite | -9982.995607958073 | True | None |
| Pd1/Co(111) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=51600ffb259aca92e11c] | SAA | -15011.694438828586 | True | None |
| (Ti,Cr)3InN2(edge_1120) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=518e823e05a53fe97416] | MAXPhase | -6245.938783792247 | True | 0.0 |
| CoNiHfRe-fcc_bcc_dual(100) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=562943661cfb72e7e072] | HEA | -12858.746404798898 | True | None |
| Au-Oxalate-N3P (16.0) MOF Network (Furukawa et al., Science 2013) [id=56c74d76175883fef1fb] | MOF | -4442.816673314329 | True | None |
| Zr-doped ZnCo2O4/graphene_oxide (mesoporous) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=68ac3d4c5e69f37e7bf6] | Spinel | -6215.8431567296975 | True | None |
| Na0.7Ni0.3TaO3(none) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=6b2e65d5cbfdbcddcddd] | Perovskite | -8459.986037265156 | True | None |
| Ag1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=730c6d4349990c81e975] | SAA | -15290.492962132705 | True | None |
| CoReIrPtCe-fcc(100) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=78c3bca65991f4d582d7] | HEA | -14491.62428573441 | True | None |
| Pt1/Ag(111) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=7a83fee91102006b288b] | SAA | -13552.411414405176 | True | None |
| Pt-Boroxine-O4 (10.0) COF Network (Furukawa et al., Science 2013) [id=7e36da26795823ef060e] | COF | -3944.1990820179417 | True | None |
| N,Ta-doped Au(bcc211)/MoO3, strain=+0.10% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=7ec3c045b95ab35fc877] | SolidCatalyst | -11481.452978322308 | True | None |
| PtSn-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=7fba4d6b4f15772edf2d] | DAC | -3394.407995710922 | True | 0.0 |
| MoPr-N2S2N2/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=84b174c8d9a033265616] | DAC | -2394.096853799505 | True | None |
| Fe-N3S-MoS2_sheet/CO Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=86d4a6562ed06d12f482] | SAC | -5324.039498464921 | True | None |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | MetalFreeCarbon | -1091.5944437901735 | True | None |
| FePr-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=9d534bdb73c262d40580] | DAC | -5009.160288729027 | True | None |
| Pd1/Ag(211) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=a0427c41d176ccfdf6bd] | SAA | -14709.579381115902 | True | None |
| (Nb,Fe)3PN2(edge_1010) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=a0c3c0f2b2c0410501ef] | MAXPhase | -9490.602907121067 | True | None |
| Na0.75In0.25FeO3(none) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=a8fcab58c4ea5ad84ca2] | Perovskite | -8310.367166236225 | True | None |
| Pd1/Ga(111) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=ab7e183d66b2e49921ed] | SAA | -14454.362777659382 | True | None |
| Pt1/Ag(211) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=bb928f4c40ee143d7174] | SAA | -13552.411414405176 | True | None |
| Mn,Cu-doped Au(fcc110)/NiP, strain=+0.06% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=c5ad001debb359e04c33] | SolidCatalyst | -14082.440665064487 | True | None |
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | SAA | -13855.521142572605 | True | None |
| Pd-N3C-Phosphorene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=cae35c09006e266f2371] | SAC | -3973.6852623295295 | True | None |
| Mn-Imidazolate-O4 (16.0) MOF Network (Furukawa et al., Science 2013) [id=ccd90970bc3389aec1f9] | MOF | -3953.050980386165 | True | None |
| Cu-S4-Phosphorene/pyridine Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=cdf618516215d8767394] | SAC | -3429.4082782699315 | True | None |
| Ce-doped CuCu2O4/graphene_oxide (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=cef5b9bdd5a3145b4e70] | Spinel | -10272.120086750194 | True | None |
| Ir1/In(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=d59d0d0dd2725f5a414b] | SAA | -8884.04500998612 | True | None |
| Co-N2S2-Phosphorene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=d6eb6a966d4dcb7708c7] | SAC | -4905.875084124234 | True | None |
| W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] | Spinel | -7061.884303322515 | True | None |
| Pd1/Ag(111) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=ead9f66386a12a72e2aa] | SAA | -14709.579381115902 | True | None |
| Nb-doped FeCu2O4/carbon_black (nanosheet) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=ef9803d8396fe9b5c300] | Spinel | -8188.144914801734 | True | None |
| Pd1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=f353706418ed8a755f1f] | SAA | -14364.996728013557 | True | None |

## Phase 4: Quantum VQE Transition States

| Catalyst | Reaction | Energy (Ha) | Energy (eV) | Qubits | Layers |
|----------|----------|-------------|-------------|--------|--------|
| Br,Mg-doped Ce(bcc211)/Graphene, strain=+0.07% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=000971115da8ddd8cf8f] | CH_split | -39.698082 | -1080.2398 | 4 | 3 |
| FeZnNbIrSc-hcp(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=0021a267d0b5a6e984a9] | CH_split | -39.697050 | -1080.2118 | 4 | 3 |
| ZnAl-N3BN3/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=0433514342912e4613e1] | CH_split | -39.687833 | -1079.9609 | 4 | 3 |
| Ti0.10Te0.90 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=54757bb17c1b638462d9] | CH_split | -39.696085 | -1080.1855 | 4 | 3 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | CH_split | -39.685754 | -1079.9044 | 4 | 3 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | CH_split | -39.700880 | -1080.3160 | 4 | 3 |
| Pd1/In(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=8796cf31edc4797d83ec] | CH_split | -39.680162 | -1079.7522 | 4 | 3 |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | ORR | -148.729949 | -4047.1481 | 4 | 3 |
| Ga1/Pd(111) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=ab7e183d66b2e49921ed] | ORR | -148.731857 | -4047.2000 | 4 | 3 |
| Zn-N4_pyridine-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=b1bd7298e14db77a573c] | CH_split | -39.689107 | -1079.9956 | 4 | 3 |
| Co1/Pt(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | ORR | -148.729489 | -4047.1356 | 4 | 3 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | CH_split | -39.693568 | -1080.1170 | 4 | 3 |
| FeZnZrAgIr-bcc(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=d4403c08f769fc70a3a4] | CH_split | -39.685972 | -1079.9103 | 4 | 3 |
| W-doped NiCo2O4/N-CNT (0.15) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] | ORR | -148.732373 | -4047.2140 | 4 | 3 |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | CH_split | -39.688771 | -1079.9865 | 4 | 3 |

## Phase 5: PEMFC Fuel Cell Performance

### Single-Cell Performance

| Cathode | Membrane | Peak Power (W/cm²) | OCV (V) | η (V) | Efficiency (Rated) |
|---------|----------|---------------------|---------|-------|--------------------|
| W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] | AquivionE87 | 2.2532 | 1.179 | 0.159087 | 45.0% |
| Pd1/Ga(111) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=ab7e183d66b2e49921ed] | AquivionE87 | 2.3160 | 1.179 | 0.205863 | 46.2% |
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | AquivionE87 | 2.3119 | 1.179 | 0.2069 | 46.1% |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | AquivionE87 | 2.0398 | 1.179 | 0.172376 | 40.7% |
| Pt1/Ga(211) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=101159f93bd9c6d9f6f8] | AquivionE87 | 2.0805 | 1.179 | 0.264738 | 41.5% |
| Pd-N3C-Phosphorene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=cae35c09006e266f2371] | AquivionE87 | 1.9356 | 1.179 | 0.281737 | 38.6% |
| Mn,Cu-doped Au(fcc110)/NiP, strain=+0.06% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=c5ad001debb359e04c33] | AquivionE87 | 1.7496 | 1.179 | 0.349075 | 34.9% |
| Ag1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=730c6d4349990c81e975] | AquivionE87 | 1.7441 | 1.179 | 0.348847 | 34.8% |
| Cu-S4-Phosphorene/pyridine Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=cdf618516215d8767394] | AquivionE87 | 1.6811 | 1.179 | 0.345368 | 33.5% |
| Cu-C4-Phosphorene/pyridine Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=383a5786185703639bcd] | AquivionE87 | 1.6712 | 1.179 | 0.347835 | 33.3% |
| Co-N2S2-Phosphorene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=d6eb6a966d4dcb7708c7] | AquivionE87 | 1.6544 | 1.179 | 0.35203 | 33.0% |
| Pt-Boroxine-O4 (10.0) COF Network (Furukawa et al., Science 2013) [id=7e36da26795823ef060e] | AquivionE87 | 1.6168 | 1.179 | 0.342209 | 32.3% |
| Pd1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=3ab787b3c8004744ae8e] | AquivionE87 | 1.6628 | 1.179 | 0.36916 | 33.2% |
| Pd1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=f353706418ed8a755f1f] | AquivionE87 | 1.6072 | 1.179 | 0.383056 | 32.1% |
| BaSbO3(B_vacancy) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=06d7172dc944ea7d38c7] | AquivionE87 | 1.4412 | 1.179 | 0.354082 | 28.8% |

### Stack Performance

- **Stack Power**: 253.2 kW net (400 cells)
- **System Efficiency**: 42.1%
- **Power Density**: 2732 W/kg, 4122 W/L
- **Cost**: $31/kW

## Techno-Economic Summary

### Turquoise H₂ Production Cost Drivers

| Parameter | Value |
|-----------|-------|
| Catalyst material cost | Based on abundance scoring |
| Reactor operating temperature | 800–1200 K |
| Target H₂ production cost | < $2/kg |
| Carbon co-product credit | $0.10–$0.50/kg |

### Levelized Cost of Hydrogen (LCOH) at X = 26.3% (T = 1300 K)

| Scenario | Single-Pass X | Net CH₄ Util. | H₂ Cost ($/kg) | Carbon Credit ($/kg H₂) | DOE < $2/kg |
|----------|---------------|---------------|----------------|-------------------------|-------------|
| optimistic | 26.3% | 89.9% | $0.44 | $1.00 | Yes |
| base | 26.3% | 89.9% | $1.64 | $0.80 | Yes |
| conservative | 26.3% | 89.9% | $5.51 | $0.20 | No |

### Multi-Pathway Discovery Champions & QE-DFT-Calibrated Leaderboard

#### 11-Champion Discovery Portfolio (Interfacial & QE DFT Validated)

| Rank | Literature Catalyst Name | Class | Screening E_act (eV) | Calibrated E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | QE DFT Energy (Ry) |
|------|--------------------------|-------|----------------------|-----------------------|------------------------|-----------------------|--------------------|
| 1 | (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | MAXPhase | 0.0146 | 0.6786 | 0.2193 | 0.8833 | -318.7712 |
| 2 | FeZnZrAgIr-bcc(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=d4403c08f769fc70a3a4] | HEA | 0.1192 | 0.7832 | 0.2370 | 0.9010 | -1029.8175 |
| 3 | Br,Mg-doped Ce(bcc211)/Graphene, strain=+0.07% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=000971115da8ddd8cf8f] | SolidCatalyst | 0.1261 | 0.7901 | 1.8634 | 2.5274 | -1492.9330 |
| 4 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | 0.3413 | 1.0053 | 0.5255 | 1.1895 | -1184.6284 |
| 5 | Zn-N4_pyridine-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=b1bd7298e14db77a573c] | SAC | 0.3640 | 1.0280 | 2.3188 | 2.9828 | -500.4341 |
| 6 | Ti0.10Te0.90 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=54757bb17c1b638462d9] | MoltenMetal | 0.6139 | 1.7919 | 0.8084 | 1.9863 | -198.5422 |
| 7 | Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | 0.6281 | 1.2921 | 1.7489 | 2.4129 | -393.6177 |
| 8 | ZnAl-N3BN3/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=0433514342912e4613e1] | DAC | 0.8696 | 1.5336 | 1.5889 | 2.2529 | -500.4341 |
| 9 | pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | 1.2468 | 1.9108 | 2.0556 | 2.7196 | -45.0008 |
| 10 | Pd1/In(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=8796cf31edc4797d83ec] | SAA | 2.1596 | 2.8236 | 3.1934 | 3.8574 | -796.0154 |
| 11 | FeZnNbIrSc-hcp(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=0021a267d0b5a6e984a9] | HEA | 0.3168 | 0.9808 | 0.5728 | 1.2368 | -1069.5552 |

#### 7-Pathway QE-DFT-Calibrated Reactor & LCOH Leaderboard

| Pathway | Reactor | Literature Champion Catalyst | T (K) | Uncalibrated X | Calibrated X | Ω | ΔP (bar) | Base LCOH ($/kg) | Optimistic LCOH ($/kg) | Kinetics Status |
|---------|---------|------------------------------|-------|----------------|--------------|---|----------|------------------|------------------------|-----------------|
| `thermocatalytic` | PBR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 96.58% | 53.88% | 0.965 | 0.913 | $1.208 | $0.208 | `dft_neb_complete` |
| `thermocatalytic_pbr` | PBR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 96.58% | 53.88% | 0.965 | 0.913 | $1.208 | $0.208 | `dft_neb_complete` |
| `thermocatalytic_pfr` | PFR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 94.01% | 46.35% | — | 0.399 | $1.278 | $0.246 | `dft_neb_complete` |
| `thermocatalytic_fluidized` | Fluidized | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 41.27% | 11.11% | — | 0.399 | $2.874 | $1.108 | `dft_neb_complete` |
| `mmbcr` | MMBCR | Ti0.10Te0.90 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=54757bb17c1b638462d9] | 1300.0 | 24.17% | 2.73% | — | — | $9.318 | $4.592 | `dft_neb_complete` |
| `ntec` | NTEC | Ti0.10Ga0.90 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | 473.1 | 1.34% | 0.02% | — | — | $941.355 | $572.876 | `dft_neb_complete` |
| `electrochemical` | Electrochemical | (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | 823.1 | 0.48% | 0.04% | — | — | $1122.075 | $666.203 | `dft_neb_complete` |

#### Supported-Ni(111)/SiO₂ Filamentous-Carbon Reference Closure

- **Reference Catalyst**: Ni(111)/SiO2 Supported Nanoparticle Reference (Bengaard 2002 / Ermakova 2000 / Takenaka 2003 Filamentous-Carbon Judge)
- **Closure Passed**: `True` (PFR ROI X = 68.42%, PBR ROI X = 61.85%, Fluidized ROI X = 64.18%)
- **Barrier Sensitivity ($|d\ln X / dE_{\text{act}}|$)**: 8.27 eV⁻¹ (gate ≥ 5.0 eV⁻¹)
- **Particle-Size Turnover Differentiation**: 34.3% (gate > 20%)
- **Encapsulation Onset ($\theta_{\text{encap}}$)**: 0.145 (gate ≥ 0.10)
- **Specific Carbon Yield ($Y_{\text{C}}$)**: 923 K = 8.47 g_C/(g_Ni·h), 973 K = 9.12 g_C/(g_Ni·h) (band 8–10 g_C/(g_Ni·h))

#### Coupled Turquoise H₂ → PEMFC Stack & Electricity Cost (LCOE)

- **QE-Confirmed Acid-Stable PEMFC Champion**: W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] ($\eta_{\text{ORR,GC}} = 0.1591\text{ V}$, Pourbaix margin $+0.0568\text{ V}$, peak power **2.2532 W/cm²**, stack net power **253.22 kW** @ 42.10% system efficiency)
- **Low-PGM Champion**: Pd1/Ga(111) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=ab7e183d66b2e49921ed] ($\eta_{\text{ORR,GC}} = 0.2059\text{ V}$, Pourbaix margin $+0.4123\text{ V}$, peak power **2.3160 W/cm²**)

| Pathway | Calibrated Base H₂ ($/kg) | Calibrated Optimistic H₂ ($/kg) | Coupled Base PEMFC Fuel Cost ($/kWh_e) | Coupled Optimistic PEMFC Fuel Cost ($/kWh_e) |
|---------|---------------------------|---------------------------------|----------------------------------------|----------------------------------------------|
| `electrochemical` | $1122.075 | $666.203 | $79.9658 | $47.4776 |
| `mmbcr` | $9.318 | $4.592 | $0.6641 | $0.3273 |
| `ntec` | $941.355 | $572.876 | $67.0866 | $40.8266 |
| `thermocatalytic` | $1.208 | $0.208 | $0.0861 | $0.0148 |
| `thermocatalytic_fluidized` | $2.874 | $1.108 | $0.2048 | $0.0790 |
| `thermocatalytic_pbr` | $1.208 | $0.208 | $0.0861 | $0.0148 |
| `thermocatalytic_pfr` | $1.278 | $0.246 | $0.0911 | $0.0175 |

## Pipeline Timing

| Phase | Duration |
|-------|----------|
| Fairchem Screening + GA | 838s |
| Reactor Simulation | 11s |
| DFT Validation | 42s |
| VQE Quantum | 116s |
| Fuel Cell Modeling | 694s |
| Report Generation | 0s |
| **Total** | **1701s (0.5 hours)** |
