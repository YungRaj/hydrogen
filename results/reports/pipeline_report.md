# Turquoise Hydrogen → Fuel Cell: Pipeline Report

**Generated**: 2026-11-07 05:52:38
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
| Solids judge CH₄ conversion (Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] PFR) | 26.34% |
| Top candidate solids CH₄ conversion (FeZnNbIrSc-hcp(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=0021a267d0b5a6e984a9]) | 98.50% |
| MMBCR max CH₄ conversion | 73.83% (not a solids rank) |
| NTEC max CH₄ conversion (coupled MOOSE+Cantera) | 1.34% |
| Electrochemical max CH₄ conversion (coupled MOOSE+Cantera) | 0.67% |
| FC catalysts evaluated | 104 |
| Best PEMFC power density | 1.348 W/cm² |
| Best PEMFC efficiency | 35.9% |
| Min ORR overpotential | 0.297 V |

## Phase 1: Catalyst Screening & Genetic Optimization

### Screening Database: `branch_ranker_evidence`

- **Total entries**: 90
- **Valid candidates**: 64

**Distribution by Material Class:**

| Material Class | Count | Valid |
|----------------|-------|-------|
| MoltenMetal | 6 | 6 |
| SolidCatalyst | 7 | 6 |
| SAC | 6 | 6 |
| DAC | 7 | 7 |
| MOF | 3 | 3 |
| COF | 3 | 2 |
| Perovskite | 3 | 3 |
| MetalHydride | 3 | 3 |
| MAXPhase | 7 | 4 |
| HEA | 25 | 8 |
| Spinel | 6 | 5 |
| MXene | 3 | 3 |
| SAA | 4 | 3 |
| MetalFreeCarbon | 6 | 4 |
| SupportedNanoparticle | 1 | 1 |


**Top 10 Catalysts by Activation Barrier:**

| Rank | Genome | E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | Coking Index | dE_H (eV) |
|------|--------|------------|------------------------|-----------------------|-------------|-----------|
| 1 | Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | 0.2753 | 0.5169 | 1.1809 | 2.479209050536163 | -0.9519897475838696 |
| 2 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 0.3413 | 0.5255 | 1.1895 | 0.7676210526315789 | -0.4912 |
| 3 | La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | 0.4655 | 0.7784 | 1.4424 | 2.131796821951859 | -0.6455789580941236 |
| 4 | FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | 0.4742 | 0.8961 | 1.5601 | 2.334771141409924 | -0.9486385360360464 |
| 5 | (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | 0.4989 | 0.6573 | 1.3213 | 2.9311808198690343 | -0.5749312415719068 |
| 6 | Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | 0.6281 | 1.7489 | 2.4129 | 1.9931449118678044 | -0.29857622185007 |
| 7 | CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | 0.7474 | 1.1337 | 1.7977 | 3.669448837637922 | -0.3556362167000806 |
| 8 | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | 0.8215 | 1.0165 | 2.1945 | nan | 0.0601983252881605 |
| 9 | Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | 0.8274 | 0.9524 | 1.6164 | 5.305162179504616 | 0.8368801762489819 |
| 10 | Nd-N4_pyridine-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=8609ad5ba35faea2ac57] | 0.9407 | 1.6838 | 2.3478 | 5.021425515515548 | -0.7594639325536012 |

## Phase 2: Reactor-Scale Simulation

Solids judged on **Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit]** (single-pass X, `a`, WHSV, ΔP). MMBCR X_eq and 0.01 eV H-parked cats are not ranks.

| Reactor | Catalyst | T (K) | Single-pass X | H₂ Selectivity | H₂ Yield | r_H₂ (mol/m³·s) | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |
|---------|----------|-------|---------------|----------------|----------|-----------------|---------|------------|----------|
| PFR | Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | 973.15 | 26.34% | 100.00% | 26.34% | — | 4154 | 900 | 0.33 |
| PBR | Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | 973.15 | 23.03% | 100.00% | 23.03% | — | 4154 | 1272 | 0.50 |
| Fluidized | Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | 973.15 | 18.79% | 100.00% | 18.79% | — | 3808 | 400 | 0.33 |

Top candidate solids X = 98.50% on **FeZnNbIrSc-hcp(111) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=0021a267d0b5a6e984a9]** (PBR @ 1300.0 K)

### Solids runs (single-pass)

| Catalyst | Reactor | T (K) | Single-pass X | H₂ Selectivity | H₂ Yield | r_H₂ (mol/m³·s) | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |
|----------|---------|-------|---------------|----------------|----------|-----------------|---------|------------|----------|
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PBR | 1300.0 | 97.09% | 99.87% | 96.97% | 5.34 | 4154 | 2783 | 0.91 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PFR | 1300.0 | 97.05% | 99.62% | 96.68% | 1.72 | 4154 | 900 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PBR | 1300.0 | 96.58% | 99.87% | 96.46% | 5.64 | 4154 | 2956 | 0.91 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PFR | 1300.0 | 94.01% | 99.54% | 93.57% | 1.67 | 4154 | 900 | 0.40 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PFR | 1100.0 | 72.54% | 100.00% | 72.54% | 1.53 | 4154 | 900 | 0.36 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PBR | 1100.0 | 60.18% | 100.00% | 60.18% | 2.39 | 4154 | 1696 | 0.84 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | Fluidized | 1300.0 | 48.45% | 97.75% | 47.36% | 0.37 | 3808 | 400 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PFR | 1100.0 | 44.65% | 99.99% | 44.64% | 0.94 | 4154 | 900 | 0.36 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | Fluidized | 1100.0 | 43.84% | 99.99% | 43.84% | 0.41 | 3808 | 400 | 0.36 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | Fluidized | 1300.0 | 41.27% | 97.53% | 40.25% | 0.32 | 3808 | 400 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PBR | 1100.0 | 40.05% | 100.00% | 40.05% | 1.47 | 4154 | 1572 | 0.73 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | Fluidized | 1100.0 | 27.35% | 99.97% | 27.35% | 0.26 | 3808 | 400 | 0.36 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PFR | 973.15 | 26.34% | 100.00% | 26.34% | 0.63 | 4154 | 900 | 0.33 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PBR | 973.15 | 23.03% | 100.00% | 23.03% | 0.77 | 4154 | 1272 | 0.50 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | Fluidized | 973.15 | 18.79% | 100.00% | 18.79% | 0.20 | 3808 | 400 | 0.33 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PFR | 923.15 | 14.56% | 100.00% | 14.56% | 0.37 | 4154 | 900 | 0.32 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PBR | 923.15 | 13.04% | 100.00% | 13.04% | 0.43 | 4154 | 1185 | 0.44 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | Fluidized | 923.15 | 10.70% | 100.00% | 10.70% | 0.12 | 3808 | 400 | 0.32 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PFR | 900.0 | 10.69% | 100.00% | 10.69% | 0.27 | 4154 | 900 | 0.32 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PBR | 900.0 | 9.66% | 100.00% | 9.66% | 0.32 | 4154 | 1155 | 0.42 |

### Liquid-Metal, NTEC, and Electrochemical Runs

| Catalyst | Reactor | T (K) | CH₄ Conversion | H₂ Selectivity | H₂ Yield | r_H₂ (mol/m³·s) | Evidence Tier |
|----------|---------|-------|----------------|----------------|----------|-----------------|---------------|
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | MMBCR | 1300.0 | 73.83% | 98.49% | 72.72% | 0.43 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | MMBCR | 1300.0 | 67.88% | 98.29% | 66.72% | 0.40 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | MMBCR | 1100.0 | 51.16% | 99.99% | 51.15% | 0.36 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | MMBCR | 1100.0 | 44.51% | 99.99% | 44.51% | 0.31 | diagnostic_screening_template |
| Bi89Ga11 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=0e7e12a7c4a5ff0d18a9] | MMBCR | 1300.0 | 40.58% | 97.23% | 39.46% | 0.23 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | MMBCR | 973.15 | 31.12% | 100.00% | 31.12% | 0.25 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | MMBCR | 973.15 | 25.69% | 100.00% | 25.69% | 0.20 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | MMBCR | 923.15 | 23.26% | 100.00% | 23.26% | 0.19 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | MMBCR | 900.0 | 20.04% | 100.00% | 20.04% | 0.17 | diagnostic_screening_template |
| Bi89Ga11 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=0e7e12a7c4a5ff0d18a9] | MMBCR | 1100.0 | 20.01% | 99.95% | 20.00% | 0.14 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | MMBCR | 923.15 | 18.73% | 100.00% | 18.73% | 0.16 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | MMBCR | 900.0 | 15.93% | 100.00% | 15.93% | 0.14 | diagnostic_screening_template |
| Bi83Ni17 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=87880819e7c1e902b8c0] | MMBCR | 1300.0 | 11.61% | 96.23% | 11.17% | 0.07 | diagnostic_screening_template |
| Bi89Ga11 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=0e7e12a7c4a5ff0d18a9] | MMBCR | 973.15 | 8.97% | 100.00% | 8.97% | 0.07 | diagnostic_screening_template |
| Bi89Ga11 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=0e7e12a7c4a5ff0d18a9] | MMBCR | 923.15 | 5.83% | 100.00% | 5.83% | 0.05 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | MMBCR | 773.15 | 5.62% | 100.00% | 5.62% | 0.06 | diagnostic_screening_template |
| Bi89Ga11 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=0e7e12a7c4a5ff0d18a9] | MMBCR | 900.0 | 4.69% | 100.00% | 4.69% | 0.04 | diagnostic_screening_template |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | MMBCR | 773.15 | 4.13% | 100.00% | 4.13% | 0.04 | diagnostic_screening_template |
| Bi83Ni17 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=87880819e7c1e902b8c0] | MMBCR | 1100.0 | 3.75% | 99.74% | 3.74% | 0.03 | diagnostic_screening_template |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MMBCR | 1300.0 | 3.35% | 95.94% | 3.21% | 0.02 | diagnostic_screening_template |

Excluded 324 non-run reactor files (sweeps, mocks, not_applicable, validation_required, failed) from performance claims.

## Phase 3: DFT Validation (Quantum ESPRESSO)

| Catalyst | Material Class | Energy (eV) | Converged | Max Force |
|----------|----------------|-------------|-----------|-----------|
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | -7684.891466 | True | 0.0 |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | -5835.70603 | True | 0.0 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | -612.2664572547036 | True | 0.0 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | -16117.69086882 | True | 0.0 |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | MoltenMetal | -24768.9406 | True | 0.000412 |
| Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | SAA | -11458.505617 | True | 0.0 |
| La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | SolidCatalyst | -17938.329641 | True | 0.0 |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | -6971.842105 | True | 0.0 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | -5355.441756157839 | True | 0.0 |
| FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | HEA | -15302.842685 | True | 0.0 |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | MAXPhase | -4337.103472226402 | True | 0.0 |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | -7053.083039 | True | 0.0 |
| Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | Spinel | -5671.253041 | True | 0.0 |

### Periodic Surface-Slab Quantum ESPRESSO Verification (`pw.x` Slab + Adsorbate & Saddle)

| Catalyst | Class | E_clean,slab (eV) | ΔE_rxn,slab(CH₄→CH₃*+H*) (eV) | E_act,slab,DFT (eV) | ΔE_C,slab (eV) | Converged |
|----------|-------|-------------------|-------------------------------|---------------------|----------------|-----------|
| Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | SAA | -35844.5773 | +0.4702 | 1.1809 | +1.7792 | True |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | -31555.4425 | +0.4764 | 1.1895 | +0.0676 | True |
| La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | SolidCatalyst | -37298.4091 | +0.6585 | 1.4424 | +1.4318 | True |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | -28490.5673 | +0.5713 | 1.3213 | +2.2312 | True |
| FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | HEA | -28677.0160 | +0.7433 | 1.5601 | +1.6348 | True |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | -28003.6274 | +2.2013 | 2.4129 | +1.2931 | True |
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | -13955.7626 | +1.3348 | 1.7977 | +2.9694 | True |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | -26244.2822 | +1.2000 | 2.1945 | -0.2000 | True |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | -10492.2452 | +0.7838 | 1.6164 | +4.6052 | True |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | -2006.5712 | +2.6332 | 2.7196 | +5.9712 | True |
| Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | Spinel | -28003.7011 | +3.0300 | 3.0013 | +1.9695 | True |

### 7-Step Elementary Dehydrogenation & Coking CI-NEB Network (Fairchem 3D + DFT Calibrated)

| Catalyst | Class | CH₄→CH₃* (eV) | CH₃*→CH₂* (eV) | CH₂*→CH* (eV) | CH*→C* (eV) | 2H*→H₂ (eV) | C* Diff. (eV) | C*+C*→C₂* (eV) | Steps Conv. |
|----------|-------|---------------|----------------|---------------|-------------|-------------|---------------|----------------|-------------|
| Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | SAA | 1.181 | 1.152 | 1.192 | 1.112 | 1.450 | 1.613 | 1.643 | 7/7 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | 1.190 | 1.018 | 1.058 | 0.978 | 1.450 | 1.408 | 1.438 | 7/7 |
| La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | SolidCatalyst | 1.442 | 1.385 | 1.425 | 1.345 | 1.450 | 1.572 | 1.602 | 7/7 |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | 1.321 | 1.330 | 1.370 | 1.290 | 1.450 | 1.640 | 1.670 | 7/7 |
| FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | HEA | 1.560 | 1.519 | 1.559 | 1.479 | 1.450 | 1.596 | 1.626 | 7/7 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | 2.413 | 2.344 | 2.384 | 2.304 | 1.450 | 1.555 | 1.585 | 7/7 |
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | 1.798 | 1.868 | 1.908 | 1.828 | 1.450 | 1.640 | 1.670 | 7/7 |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | 2.195 | 2.001 | 2.041 | 1.961 | 1.450 | 1.376 | 1.406 | 7/7 |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | 1.616 | 1.823 | 1.863 | 1.783 | 1.450 | 1.640 | 1.670 | 7/7 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | 2.720 | 3.040 | 3.080 | 3.000 | 1.450 | 1.640 | 1.670 | 7/7 |
| Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | Spinel | 3.001 | 2.989 | 3.029 | 2.949 | 1.450 | 1.636 | 1.666 | 7/7 |

### Finite-Temperature Langevin NVT Liquid-Alloy Molecular Dynamics (`LiquidMDProtocol`)

- **Catalyst**: Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] at $T = 1273.2\text{ K}$ (8 decorrelated NVT snapshots)
- **Mean Potential Energy**: -62.2034 ± 0.2095 eV
- **Active Promoter Vertical Offset ($\Delta z_{\text{promoter-host}}$)**: -0.2816 Å (Gibbs surface-tension sub-surface relaxation)

## Phase 4: Quantum VQE Transition States

| Catalyst | Reaction | Energy (Ha) | Energy (eV) | Qubits | Layers |
|----------|----------|-------------|-------------|--------|--------|
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | CH_split | -125.485697 | -3414.6398 | 8 | 3 |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | CH_split | -208.197442 | -5665.3410 | 8 | 3 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | CH_split | -94.606711 | -2574.3797 | 8 | 3 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | CH_split | -208.192487 | -5665.2062 | 8 | 3 |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | ORR | -205.049811 | -5579.6896 | 8 | 3 |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | CH_split | -208.193116 | -5665.2233 | 8 | 3 |
| Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | CH_split | -158.278058 | -4306.9654 | 8 | 3 |
| La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | CH_split | -208.194538 | -5665.2620 | 8 | 3 |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | CH_split | -97.537007 | -2654.1172 | 8 | 3 |
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | ORR | -268.549084 | -7307.5929 | 8 | 3 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | CH_split | -184.646811 | -5024.4957 | 8 | 3 |
| W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] | ORR | -295.287102 | -8035.1714 | 8 | 3 |
| FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | CH_split | -133.073573 | -3621.1164 | 8 | 3 |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | CH_split | -97.539149 | -2654.1755 | 8 | 3 |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | CH_split | -97.537646 | -2654.1346 | 8 | 3 |
| Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | CH_split | -184.646075 | -5024.4757 | 8 | 3 |

### Multi-Reference `CASSCF(4,4)` Active-Space Diagnostics (`lanl2dz+6-31g*`, 8 Qubits vs. Single-Reference DFT)

| Catalyst | Active TS Motif | ΔE_corr,MR (eV) | Leading CI |c₀|² | CCSD T₁ Diag. | E_act,DFT-PBE (eV) | E_act,VQE-CAS (eV) | ΔΔE‡_MR-DFT (eV) |
|----------|-----------------|-----------------|------------------|---------------|--------------------|--------------------|------------------|
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | `Cr-CH4_ts` | -4.1029 | 0.9750 | 0.1128 | 0.0800 | 1.4732 | +1.3932 |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | `Ni-CH4_ts` | -0.4010 | 0.9623 | 0.0827 | 0.0800 | 2.2182 | +2.1382 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | `NH-CH4_ts` | -3.6095 | 0.4954 | 0.0136 | 7.0071 | 9.1224 | +2.1153 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | `Ni-CH4_ts` | -0.4022 | 0.9622 | 0.0834 | 0.0800 | 2.3560 | +2.2760 |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | `NH-O2H2_ts` | -5.2236 | 0.4836 | 0.1096 | 0.9026 | 4.1749 | +3.2723 |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | `Ni-CH4_ts` | -0.4021 | 0.9622 | 0.0833 | 0.0800 | 1.9624 | +1.8824 |
| Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | `Pt-CH4_ts` | -0.4761 | 0.9728 | 0.0295 | 0.6459 | 2.1062 | +1.4602 |
| La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | `Ni-CH4_ts` | -0.4017 | 0.9622 | 0.0831 | 0.0800 | 2.1954 | +2.1154 |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | `Ti-CH4_ts` | -6.3843 | 0.4751 | 0.1663 | 2.4706 | 0.1740 | -2.2966 |
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | `Pt-O2H2_ts` | -5.6491 | 0.4468 | 0.0935 | 0.0800 | 0.1201 | +0.0401 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | `CoH-CH4_ts` | -1.8800 | 0.9609 | 0.0265 | 0.0800 | 0.7899 | +0.7099 |
| W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] | `CoH-O2H2_ts` | -3.2279 | 0.8056 | 0.0339 | 0.0800 | 0.0800 | +0.0000 |
| FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | `Ru-CH4_ts` | -3.0018 | 0.4309 | 0.0168 | 0.0800 | 0.6251 | +0.5451 |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | `Ti-CH4_ts` | -6.3827 | 0.4750 | 0.1651 | 1.0431 | 0.1153 | -0.9278 |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | `Ti-CH4_ts` | -6.3838 | 0.4750 | 0.1653 | 2.1998 | 0.1566 | -2.0432 |
| Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | `CoH-CH4_ts` | -1.8801 | 0.9609 | 0.0265 | 0.0800 | 0.8041 | +0.7241 |

## Phase 5: PEMFC Fuel Cell Performance

### Single-Cell Performance

| Cathode | Membrane | Peak Power (W/cm²) | OCV (V) | η (V) | Efficiency (Rated) |
|---------|----------|---------------------|---------|-------|--------------------|
| W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] | AquivionE87 | 2.2531 | 1.179 | 0.1591 | 45.0% |
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | AquivionE87 | 2.3119 | 1.179 | 0.2069 | 46.1% |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | AquivionE87 | 2.0398 | 1.179 | 0.172376 | 40.7% |
| Mn,Cu-doped Au(fcc110)/NiP, strain=+0.06% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=c5ad001debb359e04c33] | AquivionE87 | 1.7496 | 1.179 | 0.349075 | 34.9% |
| Ag1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=730c6d4349990c81e975] | AquivionE87 | 1.7441 | 1.179 | 0.348847 | 34.8% |
| Pt-Boroxine-O4 (10.0) COF Network (Furukawa et al., Science 2013) [id=7e36da26795823ef060e] | AquivionE87 | 1.6168 | 1.179 | 0.342209 | 32.3% |
| Pd1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=3ab787b3c8004744ae8e] | AquivionE87 | 1.6628 | 1.179 | 0.36916 | 33.2% |
| Pd1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=f353706418ed8a755f1f] | AquivionE87 | 1.6072 | 1.179 | 0.383056 | 32.1% |
| Fe_N_C | AquivionE87 | 1.3278 | 1.179 | 0.38 | 35.3% |
| Fe_N_C | Nafion_211 | 1.3278 | 1.179 | 0.38 | 35.3% |
| BaSbO3(B_vacancy) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=06d7172dc944ea7d38c7] | AquivionE87 | 1.4412 | 1.179 | 0.354082 | 28.8% |
| Pt_wet | Nafion_212 | 1.2283 | 1.179 | 0.35 | 32.7% |
| N,Ta-doped Au(bcc211)/MoO3, strain=+0.10% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=7ec3c045b95ab35fc877] | AquivionE87 | 1.4751 | 1.179 | 0.417691 | 29.4% |
| PtSn-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=7fba4d6b4f15772edf2d] | AquivionE87 | 1.4187 | 1.179 | 0.420581 | 28.3% |
| Fe_N_C | Nafion_212 | 1.1028 | 1.179 | 0.38 | 29.3% |

### Stack Performance

- **Stack Power**: 82.6 kW net (300 cells)
- **System Efficiency**: 48.9%
- **Power Density**: 1171 W/kg, 1766 W/L
- **Cost**: $65/kW

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

### Levelized Cost of Hydrogen (LCOH) at X = 26.3% (T = 973 K)

| Scenario | Single-Pass X | Net CH₄ Util. | H₂ Cost ($/kg) | Carbon Credit ($/kg H₂) | DOE < $2/kg |
|----------|---------------|---------------|----------------|-------------------------|-------------|
| optimistic | 26.3% | 89.9% | $0.38 | $1.00 | Yes |
| base | 26.3% | 89.9% | $1.50 | $0.80 | Yes |
| conservative | 26.3% | 89.9% | $5.12 | $0.20 | No |

### Multi-Pathway Discovery Champions & QE-DFT-Calibrated Leaderboard

#### 11-Champion Discovery Portfolio (Interfacial & QE DFT Validated)

| Rank | Literature Catalyst Name | Class | Screening E_act (eV) | Calibrated E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | QE DFT Energy (Ry) |
|------|--------------------------|-------|----------------------|-----------------------|------------------------|-----------------------|--------------------|
| 1 | Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | SAA | 0.2753 | 0.9393 | 0.5169 | 1.1809 | -842.1846 |
| 2 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | 0.3413 | 1.0053 | 0.5255 | 1.1895 | -1184.6284 |
| 3 | La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | SolidCatalyst | 0.4655 | 1.1295 | 0.7784 | 1.4424 | -1318.4429 |
| 4 | (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | 0.4989 | 1.1629 | 0.6573 | 1.3213 | -428.9165 |
| 5 | FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | HEA | 0.4742 | 1.1382 | 0.8961 | 1.5601 | -1124.7382 |
| 6 | Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | 0.6281 | 1.2921 | 1.7489 | 2.4129 | -393.6177 |
| 7 | CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | 0.7474 | 1.4114 | 1.1337 | 1.7977 | -564.8291 |
| 8 | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | 0.8215 | 1.9994 | 1.0165 | 2.1945 | -512.4211 |
| 9 | Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | 0.8274 | 1.4914 | 0.9524 | 1.6164 | -518.3920 |
| 10 | pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | 1.2468 | 1.9108 | 2.0556 | 2.7196 | -45.0008 |
| 11 | Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | Spinel | 1.2189 | 1.8829 | 2.3373 | 3.0013 | -416.8294 |

#### 7-Pathway QE-DFT-Calibrated Reactor & LCOH Leaderboard

| Pathway | Reactor | Literature Champion Catalyst | T (K) | Uncalibrated X | Calibrated X | Ω | ΔP (bar) | Base LCOH ($/kg) | Optimistic LCOH ($/kg) | Kinetics Status |
|---------|---------|------------------------------|-------|----------------|--------------|---|----------|------------------|------------------------|-----------------|
| `thermocatalytic` | PBR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 96.58% | 54.87% | 0.993 | 0.913 | $1.190 | $0.200 | `dft_neb_complete` |
| `thermocatalytic_pbr` | PBR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 96.58% | 54.87% | 0.993 | 0.913 | $1.190 | $0.200 | `dft_neb_complete` |
| `thermocatalytic_pfr` | PFR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 94.01% | 46.32% | — | 0.399 | $1.273 | $0.249 | `dft_neb_complete` |
| `thermocatalytic_fluidized` | Fluidized | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 41.27% | 11.10% | — | 0.399 | $2.884 | $1.143 | `dft_neb_complete` |
| `mmbcr` | MMBCR | Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | 1300.0 | 73.83% | 40.37% | — | — | $1.362 | $0.308 | `dft_neb_complete` |
| `ntec` | NTEC | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | 473.1 | 1.34% | 0.36% | — | — | $52.902 | $31.769 | `dft_neb_complete` |
| `electrochemical` | Electrochemical | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | 823.1 | 0.48% | 0.14% | — | — | $321.607 | $190.625 | `dft_neb_complete` |

#### 7-Pathway H₂/C₂ Selectivity, Yield, & Volumetric / Site-Normalized Production Rates

| Pathway | Reactor | H₂ Selectivity | Solid C Selectivity | C₂ Byproduct Selectivity | Calibrated H₂ Yield | Uncal. r_H₂ (mol/m³·s) | Cal. r_CH₄ (mol/m³·s) | Cal. r_H₂ (mol/m³·s) | Cal. r_H₂ (kg/m³·h) | Cal. TOF_H₂ (s⁻¹) |
|---------|---------|----------------|---------------------|--------------------------|---------------------|------------------------|-----------------------|----------------------|---------------------|-------------------|
| `thermocatalytic` | PBR | 99.87% | 99.75% | 0.25% | 54.80% | 5.644 | 1.605 | 3.206 | 23.267 | 30.8734 |
| `thermocatalytic_pbr` | PBR | 99.87% | 99.75% | 0.25% | 54.80% | 5.644 | 1.605 | 3.206 | 23.267 | 30.8734 |
| `thermocatalytic_pfr` | PFR | 99.54% | 99.10% | 0.90% | 46.11% | 1.667 | 0.412 | 0.821 | 5.960 | 7.9079 |
| `thermocatalytic_fluidized` | Fluidized | 97.53% | 95.18% | 4.82% | 10.83% | 0.319 | 0.044 | 0.086 | 0.622 | 0.9003 |
| `mmbcr` | MMBCR | 98.49% | 97.13% | 2.87% | 39.76% | 0.432 | 0.120 | 0.236 | 1.713 | 35.7290 |
| `ntec` | NTEC | 100.00% | 100.00% | 0.00% | 0.36% | 0.210 | 0.028 | 0.056 | 0.409 | 0.6264 |
| `electrochemical` | Electrochemical | 100.00% | 100.00% | 0.00% | 0.14% | 0.140 | 0.006 | 0.013 | 0.091 | 0.1400 |

#### Supported-Ni(111)/SiO₂ Filamentous-Carbon Reference Closure

- **Reference Catalyst**: Ni(111)/SiO2 Supported Nanoparticle Reference (Bengaard 2002 / Ermakova 2000 / Takenaka 2003 Filamentous-Carbon Judge)
- **Closure Passed**: `True` (PFR ROI X = 68.42%, PBR ROI X = 61.85%, Fluidized ROI X = 64.18%)
- **Barrier Sensitivity ($|d\ln X / dE_{\text{act}}|$)**: 8.27 eV⁻¹ (gate ≥ 5.0 eV⁻¹)
- **Particle-Size Turnover Differentiation**: 34.3% (gate > 20%)
- **Encapsulation Onset ($\theta_{\text{encap}}$)**: 0.145 (gate ≥ 0.10)
- **Specific Carbon Yield ($Y_{\text{C}}$)**: 923 K = 8.47 g_C/(g_Ni·h), 973 K = 9.12 g_C/(g_Ni·h) (band 8–10 g_C/(g_Ni·h))

#### Coupled Turquoise H₂ → PEMFC Stack & Electricity Cost (LCOE)

- **QE-Confirmed Acid-Stable PEMFC Champion**: W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] ($\eta_{\text{ORR,GC}} = 0.1591\text{ V}$, Pourbaix margin $+0.0568\text{ V}$, peak power **2.2532 W/cm²**, stack net power **253.22 kW** @ 42.10% system efficiency)
- **Low-PGM Champion**: Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] ($\eta_{\text{ORR,GC}} = 0.2069\text{ V}$, Pourbaix margin $+0.6816\text{ V}$, peak power **2.3119 W/cm²**)

| Pathway | Calibrated Base H₂ ($/kg) | Calibrated Optimistic H₂ ($/kg) | Coupled Base PEMFC Fuel Cost ($/kWh_e) | Coupled Optimistic PEMFC Fuel Cost ($/kWh_e) |
|---------|---------------------------|---------------------------------|----------------------------------------|----------------------------------------------|
| `thermocatalytic_pfr` | $1.273 | $0.249 | $0.0907 | $0.0177 |
| `thermocatalytic_pbr` | $1.190 | $0.200 | $0.0848 | $0.0143 |
| `thermocatalytic_fluidized` | $2.884 | $1.143 | $0.2055 | $0.0815 |
| `thermocatalytic` | $1.190 | $0.200 | $0.0848 | $0.0143 |
| `mmbcr` | $1.362 | $0.308 | $0.0971 | $0.0219 |
| `ntec` | $52.902 | $31.769 | $3.7701 | $2.2641 |
| `electrochemical` | $321.607 | $190.625 | $22.9197 | $13.5851 |

### De Novo Multi-Step Atomistic, Langevin MD & Periodic Slab-DFT Discovery Validation

#### Periodic Surface-Slab Quantum ESPRESSO (`pw.x` / `neb.x`) Validation

| Rank | Catalyst | Class | Clean Slab E (eV) | ΔE_rxn,slab (eV) | E_act,slab-DFT (eV) | ΔE_C,slab (eV) | Converged |
|------|----------|-------|-------------------|------------------|---------------------|----------------|-----------|
| 1 | Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | SAA | -35844.58 | +0.4702 | 1.1809 | +1.7792 | `True` |
| 2 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | -31555.44 | +0.4764 | 1.1895 | +0.0676 | `True` |
| 3 | La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | SolidCatalyst | -37298.41 | +0.6585 | 1.4424 | +1.4318 | `True` |
| 4 | (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | -28490.57 | +0.5713 | 1.3213 | +2.2312 | `True` |
| 5 | FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | HEA | -28677.02 | +0.7433 | 1.5601 | +1.6348 | `True` |
| 6 | Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | -28003.63 | +2.2013 | 2.4129 | +1.2931 | `True` |
| 7 | CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | -13955.76 | +1.3348 | 1.7977 | +2.9694 | `True` |
| 8 | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | -26244.28 | +1.2000 | 2.1945 | -0.2000 | `True` |
| 9 | Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | -10492.25 | +0.7838 | 1.6164 | +4.6052 | `True` |
| 10 | pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | -2006.57 | +2.6332 | 2.7196 | +5.9712 | `True` |
| 11 | Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | Spinel | -28003.70 | +3.0300 | 3.0013 | +1.9695 | `True` |

#### 7-Step Methane Dehydrogenation & Coking CI-NEB Network (`CH4* → C* + 4H*`)

| Rank | Catalyst | Class | CH₄→CH₃*+H* (eV) | CH₃*→CH₂*+H* (eV) | CH₂*→CH*+H* (eV) | CH*→C*+H* (eV) | 2H*→H₂ (eV) | C_surf→C_sub (eV) | C_graphene (eV) | Steps Converged |
|------|----------|-------|------------------|-------------------|------------------|----------------|-------------|-------------------|-----------------|-----------------|
| 1 | Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | SAA | 1.1809 | 1.1522 | 1.1922 | 1.1122 | 1.4500 | 1.6135 | 1.6435 | 7/7 |
| 2 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | 1.1895 | 1.0182 | 1.0582 | 0.9782 | 1.4500 | 1.4081 | 1.4381 | 7/7 |
| 3 | La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | SolidCatalyst | 1.4424 | 1.3848 | 1.4248 | 1.3448 | 1.4500 | 1.5718 | 1.6018 | 7/7 |
| 4 | (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | 1.3213 | 1.3303 | 1.3703 | 1.2903 | 1.4500 | 1.6400 | 1.6700 | 7/7 |
| 5 | FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | HEA | 1.5601 | 1.5194 | 1.5594 | 1.4794 | 1.4500 | 1.5962 | 1.6262 | 7/7 |
| 6 | Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | 2.4129 | 2.3437 | 2.3837 | 2.3037 | 1.4500 | 1.5552 | 1.5852 | 7/7 |
| 7 | CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | 1.7977 | 1.8682 | 1.9082 | 1.8282 | 1.4500 | 1.6400 | 1.6700 | 7/7 |
| 8 | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | 2.1945 | 2.0009 | 2.0409 | 1.9609 | 1.4500 | 1.3760 | 1.4060 | 7/7 |
| 9 | Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | 1.6164 | 1.8232 | 1.8632 | 1.7832 | 1.4500 | 1.6400 | 1.6700 | 7/7 |
| 10 | pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | 2.7196 | 3.0403 | 3.0803 | 3.0003 | 1.4500 | 1.6400 | 1.6700 | 7/7 |
| 11 | Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | Spinel | 3.0013 | 2.9885 | 3.0285 | 2.9485 | 1.4500 | 1.6363 | 1.6663 | 7/7 |

#### Ab Initio / MLIP Langevin NVT Liquid-Metal Alloy Trajectory (`MoltenMetal`)

- **Champion Melt**: `('MoltenMetal', 'Ga', 'Ti', 0.1, 1273.15)` (`id=c4ab2249c156880efdbe`)
- **Ensemble & Thermostat**: `NVT_Langevin` at $T = 1273.2\text{ K}$ (8 snapshots)
- **Mean Potential Energy**: -62.2034 ± 0.2095 eV
- **Gibbs Promoter Vertical Offset ($\Delta z_{\text{promoter}}$)**: -0.2816 Å

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
