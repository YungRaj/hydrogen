# Turquoise Hydrogen → Fuel Cell: Pipeline Report

**Generated**: 2026-11-07 22:39:13
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
| Solids judge CH₄ conversion (solids PFR) | 26.34% |
| Top candidate solids CH₄ conversion (Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7]) | 98.50% |
| MMBCR max CH₄ conversion | 77.39% (not a solids rank) |
| NTEC max CH₄ conversion (coupled MOOSE+Cantera) | 84.07% |
| Electrochemical max CH₄ conversion (coupled MOOSE+Cantera) | 84.08% |
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
| 1 | Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | 0.2753 | 0.5169 | 1.1809 | 2.4792 | -0.9519897475838696 |
| 2 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 0.3413 | 0.5255 | 1.1895 | 0.7676 | -0.4912 |
| 3 | La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | 0.4655 | 0.7784 | 1.4424 | 2.1318 | -0.6455789580941236 |
| 4 | FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | 0.4742 | 0.8961 | 1.5601 | 2.3348 | -0.9486385360360464 |
| 5 | (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | 0.4989 | 0.6573 | 1.3213 | 2.9312 | -0.5749312415719068 |
| 6 | Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | 0.6281 | 1.7489 | 2.4129 | 1.9931 | -0.29857622185007 |
| 7 | CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | 0.7474 | 1.1337 | 1.7977 | 3.6694 | -0.3556362167000806 |
| 8 | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | 0.8215 | 1.0165 | 2.1945 | — (liquid melt) | 0.0601983252881605 |
| 9 | Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | 0.8274 | 0.9524 | 1.6164 | 5.3052 | 0.8368801762489819 |
| 10 | Nd-N4_pyridine-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=8609ad5ba35faea2ac57] | 0.9407 | 1.6838 | 2.3478 | 5.0214 | -0.7594639325536012 |

## Phase 2: Reactor-Scale Simulation

### Packed-Bed (PBR), Ergun Pilot Residence-Time, & Regenerated-Cycle Diagnostics

- **Completed solids runs**: 378 (PBR: 126)
- **Mean PBR overall effectiveness ($\Omega$)**: 0.946
- **Mean PBR Thiele modulus ($\phi$)**: 0.280
- **Ergun $\Delta P$ compliance**: 373/378
- **Ergun & Froment pilot residence-time band ($\tau \in [0.1, 10]\text{ s}$)**: 241/252
- **Mean Muradov regenerated-cycle minus single-pass conversion ($\Delta X$)**: +0.00%

### Solids runs (single-pass)

| Catalyst | Reactor | T (K) | Single-pass X | H₂ Selectivity | H₂ Yield | r_H₂ (mol/m³·s) | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |
|----------|---------|-------|---------------|----------------|----------|-----------------|---------|------------|----------|
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 1300.0 | 98.50% | 100.00% | 98.50% | 6.60 | 4154 | 3388 | 3.42 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 1300.0 | 98.50% | 100.00% | 98.50% | 6.60 | 4154 | 3388 | 3.42 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PBR | 1300.0 | 97.09% | 99.87% | 96.97% | 5.34 | 4154 | 2783 | 0.91 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PFR | 1300.0 | 97.05% | 99.62% | 96.68% | 1.72 | 4154 | 900 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PBR | 1300.0 | 96.58% | 99.87% | 96.46% | 5.64 | 4154 | 2956 | 0.91 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 1173.15 | 95.63% | 100.00% | 95.63% | 6.62 | 4154 | 3156 | 2.84 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 1173.15 | 95.63% | 100.00% | 95.63% | 6.62 | 4154 | 3156 | 2.84 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PFR | 1300.0 | 94.01% | 99.54% | 93.57% | 1.67 | 4154 | 900 | 0.40 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PFR | 1300.0 | 92.37% | 100.00% | 92.37% | 1.65 | 4154 | 900 | 0.40 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PFR | 1300.0 | 92.37% | 100.00% | 92.37% | 1.65 | 4154 | 900 | 0.40 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 1073.15 | 90.42% | 100.00% | 90.42% | 6.27 | 4154 | 2894 | 2.29 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 1073.15 | 90.42% | 100.00% | 90.42% | 6.27 | 4154 | 2894 | 2.29 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 973.15 | 79.96% | 100.00% | 79.96% | 5.26 | 4154 | 2487 | 1.56 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 973.15 | 79.96% | 100.00% | 79.96% | 5.26 | 4154 | 2487 | 1.56 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PFR | 1100.0 | 72.54% | 100.00% | 72.54% | 1.53 | 4154 | 900 | 0.36 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 923.15 | 72.06% | 100.00% | 72.06% | 4.40 | 4154 | 2189 | 1.08 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 923.15 | 72.06% | 100.00% | 72.06% | 4.40 | 4154 | 2189 | 1.08 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 873.15 | 62.07% | 100.00% | 62.07% | 3.37 | 4154 | 1844 | 0.73 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | PBR | 873.15 | 62.07% | 100.00% | 62.07% | 3.37 | 4154 | 1844 | 0.73 |
| Ni(fcc111)/SiO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=ni_np_lit] | PBR | 1100.0 | 60.18% | 100.00% | 60.18% | 2.39 | 4154 | 1696 | 0.84 |

### Liquid-Metal, NTEC, and Electrochemical Runs

| Catalyst | Reactor | T (K) | CH₄ Conversion | H₂ Selectivity | H₂ Yield | r_H₂ (mol/m³·s) | Evidence Tier |
|----------|---------|-------|----------------|----------------|----------|-----------------|---------------|
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | Electrochemical | 823.15 | 84.08% | 100.00% | 84.08% | 37.84 | uncalibrated_multiphysics_screening |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | NTEC | 473.15 | 84.07% | 100.00% | 84.07% | 65.82 | uncalibrated_multiphysics_screening |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | NTEC | 423.15 | 84.04% | 100.00% | 84.04% | 73.58 | uncalibrated_multiphysics_screening |
| Pt1/Ni(fcc110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | NTEC | 473.15 | 83.24% | 100.00% | 83.24% | 65.18 | uncalibrated_multiphysics_screening |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | NTEC | 473.15 | 82.93% | 100.00% | 82.93% | 64.94 | uncalibrated_multiphysics_screening |
| Cu,Fe-doped Ni(fcc111)/MgAl2O4 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=768cb0559591faec06ba] | NTEC | 473.15 | 82.78% | 100.00% | 82.78% | 64.81 | uncalibrated_multiphysics_screening |
| Pt1/Ni(fcc110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | NTEC | 423.15 | 82.66% | 100.00% | 82.66% | 72.37 | uncalibrated_multiphysics_screening |
| La,Y-doped Ni(fcc111)/CeZrO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | NTEC | 473.15 | 82.60% | 100.00% | 82.60% | 64.67 | uncalibrated_multiphysics_screening |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | NTEC | 423.15 | 82.08% | 100.00% | 82.08% | 71.87 | uncalibrated_multiphysics_screening |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | NTEC | 473.15 | 82.02% | 100.00% | 82.02% | 64.22 | uncalibrated_multiphysics_screening |
| Cu,Fe-doped Ni(fcc111)/MgAl2O4 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=768cb0559591faec06ba] | NTEC | 423.15 | 81.79% | 100.00% | 81.79% | 71.61 | uncalibrated_multiphysics_screening |
| La,Y-doped Ni(fcc111)/CeZrO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | NTEC | 423.15 | 81.45% | 100.00% | 81.45% | 71.31 | uncalibrated_multiphysics_screening |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | NTEC | 423.15 | 80.32% | 100.00% | 80.32% | 70.32 | uncalibrated_multiphysics_screening |
| Pt1/Ni(fcc110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | Electrochemical | 823.15 | 79.06% | 100.00% | 79.06% | 35.58 | uncalibrated_multiphysics_screening |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | Electrochemical | 823.15 | 77.83% | 100.00% | 77.83% | 35.03 | uncalibrated_multiphysics_screening |
| Ga85Ni15 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=017a44c337ece52f029b] | MMBCR | 1300.0 | 77.39% | 98.73% | 76.41% | 0.45 | diagnostic_screening_template |
| Ga85Ni15 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=017a44c337ece52f029b] | MMBCR | 1300.0 | 77.39% | 98.73% | 76.41% | 0.45 | diagnostic_screening_template |
| Cu,Fe-doped Ni(fcc111)/MgAl2O4 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=768cb0559591faec06ba] | Electrochemical | 823.15 | 77.28% | 100.00% | 77.28% | 34.78 | uncalibrated_multiphysics_screening |
| La,Y-doped Ni(fcc111)/CeZrO2 Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | Electrochemical | 823.15 | 76.67% | 100.00% | 76.67% | 34.51 | uncalibrated_multiphysics_screening |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=f2f562d4c2602522a2e3] | Electrochemical | 823.15 | 74.92% | 100.00% | 74.92% | 33.72 | uncalibrated_multiphysics_screening |

Excluded 341 non-run reactor files (sweeps, mocks, not_applicable, validation_required, failed) from performance claims.

## Phase 3: DFT Validation (Quantum ESPRESSO)

| Catalyst | Material Class | Energy (eV) | Converged | Max Force |
|----------|----------------|-------------|-----------|-----------|
| Ga85Ni15 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=017a44c337ece52f029b] | MoltenMetal | -8820.8348 | True | 0.00042 |
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | -7684.891466 | True | 0.0 |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | -5835.70603 | True | 0.0 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | -612.2664572547036 | True | 0.0 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | -16117.69086882 | True | 0.0 |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | MoltenMetal | -24768.9406 | True | 0.000412 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | SupportedNanoparticle | -17616.8719 | True | 0.00038 |
| Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | SAA | -11458.505617 | True | 0.0 |
| La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | SolidCatalyst | -17938.329641 | True | 0.0 |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | -6971.842105 | True | 0.0 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | -5355.441756157839 | True | 0.0 |
| FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | HEA | -15302.842685 | True | 0.0 |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | MAXPhase | -5671.1403 | True | 0.00045 |
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
| Ga85Ni15 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=017a44c337ece52f029b] | CH_split | -208.194705 | -5665.2665 | 8 | 3 |
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | CH_split | -125.485697 | -3414.6398 | 8 | 3 |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | CH_split | -208.197442 | -5665.3410 | 8 | 3 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | CH_split | -94.606711 | -2574.3797 | 8 | 3 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | CH_split | -208.192487 | -5665.2062 | 8 | 3 |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | ORR | -205.049811 | -5579.6896 | 8 | 3 |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | CH_split | -208.193116 | -5665.2233 | 8 | 3 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | CH_split | -133.093474 | -3621.6579 | 8 | 3 |
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
| Ga85Ni15 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=017a44c337ece52f029b] | `Ni-CH4_ts` | -0.4017 | 0.9622 | 0.0831 | 0.0800 | 1.8885 | +1.8085 |
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | `Cr-CH4_ts` | -4.1029 | 0.9750 | 0.1128 | 0.0800 | 1.4732 | +1.3932 |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | `Ni-CH4_ts` | -0.4010 | 0.9623 | 0.0827 | 0.0800 | 2.2182 | +2.1382 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | `NH-CH4_ts` | -3.6095 | 0.4954 | 0.0136 | 7.0071 | 9.1224 | +2.1153 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | `Ni-CH4_ts` | -0.4022 | 0.9622 | 0.0834 | 0.0800 | 2.3560 | +2.2760 |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | `NH-O2H2_ts` | -5.2236 | 0.4836 | 0.1096 | 0.9026 | 4.1749 | +3.2723 |
| Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | `Ni-CH4_ts` | -0.4021 | 0.9622 | 0.0833 | 0.0800 | 1.9624 | +1.8824 |
| Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | `Ru-CH4_ts` | -1.0310 | 0.7510 | 0.0280 | 0.0800 | 1.8980 | +1.8180 |
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
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | AquivionE87 | 2.3119 | 1.179 | 0.2069 | 46.1% |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | AquivionE87 | 2.0398 | 1.179 | 0.172376 | 40.7% |
| b5e764d896815bfcfca7 | AquivionE87 | 2.0294 | 1.179 | 0.2583 | 40.5% |
| W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] | AquivionE87 | 1.2628 | 1.179 | 0.1591 | 25.2% |
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

### Stack Performance

- **Stack Power**: 82.6 kW net (300 cells)
- **System Efficiency**: 48.9%
- **Power Density**: 1171 W/kg, 1766 W/L
- **Cost**: $65/kW

- **Stack Power**: 226.4 kW net (400 cells)
- **System Efficiency**: 37.6%
- **Power Density**: 2442 W/kg, 3685 W/L
- **Cost**: $34/kW

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

| Rank | Literature Catalyst Name | Class | Novelty | Screening E_act (eV) | Calibrated E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | QE DFT Energy (Ry) |
|------|--------------------------|-------|---------|----------------------|-----------------------|------------------------|-----------------------|--------------------|
| 1 | Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | SAA | `literature_benchmark` | 0.2753 | 0.9393 | 0.5169 | 1.1809 | -842.1846 |
| 2 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | `literature_benchmark` | 0.3413 | 1.0053 | 0.5255 | 1.1895 | -1184.6284 |
| 3 | La,Y-doped Ni(fcc211)/CeZrO2, strain=+0.09% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=be45ea11e23945dbbda5] | SolidCatalyst | `literature_benchmark` | 0.4655 | 1.1295 | 0.7784 | 1.4424 | -1318.4429 |
| 4 | (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | `de_novo_novel` | 0.4989 | 1.1629 | 0.6573 | 1.3213 | -428.9165 |
| 5 | FeNbRuSiGe-fcc_bcc_dual(211) High-Entropy Alloy (Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface) [id=e40b9a32452d7b1e6f1f] | HEA | `de_novo_novel` | 0.4742 | 1.1382 | 0.8961 | 1.5601 | -1124.7382 |
| 6 | Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Spinel | `de_novo_novel` | 0.6281 | 1.2921 | 1.7489 | 2.4129 | -393.6177 |
| 7 | CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | `de_novo_novel` | 0.7474 | 1.4114 | 1.1337 | 1.7977 | -564.8291 |
| 8 | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | `literature_benchmark` | 0.8215 | 1.9994 | 1.0165 | 2.1945 | -512.4211 |
| 9 | Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | `de_novo_novel` | 0.8274 | 1.4914 | 0.9524 | 1.6164 | -518.3920 |
| 10 | pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | `de_novo_novel` | 1.2468 | 1.9108 | 2.0556 | 2.7196 | -45.0008 |
| 11 | Ce-doped NiCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=f7e94c08fb1675631eb9] | Spinel | `de_novo_novel` | 1.2189 | 1.8829 | 2.3373 | 3.0013 | -416.8294 |

#### 7-Pathway QE-DFT-Calibrated Reactor & LCOH Leaderboard

| Pathway | Reactor | Champion Catalyst | Novelty | T (K) | Uncalibrated X | Calibrated X | Ω | ΔP (bar) | Base LCOH ($/kg) | Optimistic LCOH ($/kg) | Kinetics Status |
|---------|---------|-------------------|---------|-------|----------------|--------------|---|----------|------------------|------------------------|-----------------|
| `thermocatalytic` | PBR | Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | `de_novo_novel` | 1300.0 | 98.50% | 76.28% | 0.989 | 3.420 | $1.070 | $0.134 | `dft_neb_complete` |
| `thermocatalytic_pbr` | PBR | Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | `de_novo_novel` | 1300.0 | 98.50% | 76.28% | 0.989 | 3.420 | $1.070 | $0.134 | `dft_neb_complete` |
| `thermocatalytic_pfr` | PFR | Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | `de_novo_novel` | 1300.0 | 92.37% | 56.31% | — | 0.399 | $1.177 | $0.192 | `dft_neb_complete` |
| `thermocatalytic_fluidized` | Fluidized | Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | `de_novo_novel` | 1300.0 | 48.12% | 19.04% | — | 0.399 | $1.979 | $0.626 | `dft_neb_complete` |
| `mmbcr` | MMBCR | Ga85Ni15 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=017a44c337ece52f029b] | `de_novo_novel` | 1300.0 | 77.39% | 47.49% | — | — | $1.272 | $0.257 | `dft_neb_complete` |
| `ntec` | NTEC | (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | `de_novo_novel` | 473.1 | 84.07% | 49.77% | — | — | $0.736 | $0.005 | `dft_neb_complete` |
| `electrochemical` | Electrochemical | (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | `de_novo_novel` | 823.1 | 84.08% | 61.73% | — | — | $1.686 | $0.557 | `dft_neb_complete` |

#### De Novo Novel Catalyst Champions vs. Literature Reference Benchmarks (7-Pathway Comparison)

| Pathway | Reactor | De Novo Novel Champion | De Novo Cal. X | De Novo Base LCOH ($/kg) | Literature Benchmark Reference | Benchmark Cal. X | Benchmark Base LCOH ($/kg) | ΔX_cal Lift | ΔLCOH Saving ($/kg) |
|---------|---------|------------------------|----------------|--------------------------|--------------------------------|------------------|----------------------------|-------------|---------------------|
| `thermocatalytic` | PBR | Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | 76.28% | $1.070 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 54.87% | $1.190 | **+21.41%** | **$-0.120** |
| `thermocatalytic_pbr` | PBR | Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | 76.28% | $1.070 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 54.87% | $1.190 | **+21.41%** | **$-0.120** |
| `thermocatalytic_pfr` | PFR | Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | 56.31% | $1.177 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 46.32% | $1.273 | **+9.99%** | **$-0.096** |
| `thermocatalytic_fluidized` | Fluidized | Cu,Ru-doped Ni(fcc211)/CeZrO2 (d=10.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=9b83e69315874c5abca7] | 19.04% | $1.979 | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 11.10% | $2.884 | **+7.94%** | **$-0.905** |
| `mmbcr` | MMBCR | Ga85Ni15 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=017a44c337ece52f029b] | 47.49% | $1.272 | Bi73Ni27 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=99cbb4a430ec9289b09d] | 40.37% | $1.362 | **+7.12%** | **$-0.090** |
| `ntec` | NTEC | (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | 49.77% | $0.736 | Pt1/Ni(fcc110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | 48.82% | $0.743 | **+0.95%** | **$-0.007** |
| `electrochemical` | Electrochemical | (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | 61.73% | $1.686 | Pt1/Ni(fcc110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | 55.83% | $1.745 | **+5.90%** | **$-0.059** |

#### 7-Pathway H₂/C₂ Selectivity, Yield, & Volumetric / Site-Normalized Production Rates

| Pathway | Reactor | H₂ Selectivity | Solid C Selectivity | C₂ Byproduct Selectivity | Calibrated H₂ Yield | Uncal. r_H₂ (mol/m³·s) | Cal. r_CH₄ (mol/m³·s) | Cal. r_H₂ (mol/m³·s) | Cal. r_H₂ (kg/m³·h) | Cal. TOF_H₂ (s⁻¹) |
|---------|---------|----------------|---------------------|--------------------------|---------------------|------------------------|-----------------------|----------------------|---------------------|-------------------|
| `thermocatalytic` | PBR | 100.00% | 100.00% | 0.00% | 76.28% | 6.604 | 2.557 | 5.114 | 37.114 | 49.2473 |
| `thermocatalytic_pbr` | PBR | 100.00% | 100.00% | 0.00% | 76.28% | 6.604 | 2.557 | 5.114 | 37.114 | 49.2473 |
| `thermocatalytic_pfr` | PFR | 100.00% | 100.00% | 0.00% | 56.31% | 1.667 | 0.501 | 1.003 | 7.279 | 9.6580 |
| `thermocatalytic_fluidized` | Fluidized | 100.00% | 100.00% | 0.00% | 19.04% | 0.691 | 0.137 | 0.273 | 1.982 | 5.7396 |
| `mmbcr` | MMBCR | 98.73% | 97.59% | 2.41% | 46.89% | 0.454 | 0.141 | 0.278 | 2.020 | 51.3444 |
| `ntec` | NTEC | 100.00% | 100.00% | 0.00% | 49.77% | 65.824 | 19.485 | 38.969 | 282.808 | 432.9944 |
| `electrochemical` | Electrochemical | 100.00% | 100.00% | 0.00% | 61.73% | 37.843 | 13.891 | 27.783 | 201.623 | 308.6959 |

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
| `thermocatalytic_pfr` | $1.177 | $0.192 | $0.0839 | $0.0137 |
| `thermocatalytic_pbr` | $1.070 | $0.134 | $0.0763 | $0.0095 |
| `thermocatalytic_fluidized` | $1.979 | $0.626 | $0.1410 | $0.0446 |
| `thermocatalytic` | $1.070 | $0.134 | $0.0763 | $0.0095 |
| `mmbcr` | $1.272 | $0.257 | $0.0907 | $0.0183 |
| `ntec` | $0.736 | $0.005 | $0.0525 | $0.0004 |
| `electrochemical` | $1.686 | $0.557 | $0.1202 | $0.0397 |

#### Integrated Reactor-PSA-PEMFC Well-to-Power Coupling & CO/C₂H₄ Impurity Tolerance

| Pathway | Reactor | Post-PSA H₂ Purity | CO Slip (ppm) | C₂H₄ Slip (ppm) | Anode HOR Retention | Stack Net Power (kW) | Well-to-H₂ η (LHV) | Well-to-Power η (LHV) | Total LCOE ($/kWh_e) |
|---------|---------|--------------------|---------------|-----------------|---------------------|----------------------|--------------------|-----------------------|----------------------|
| `thermocatalytic_pfr` | PFR | 99.9970% | 0.000 | 0.000 | 100.00% | 253.21 | 68.41% | 28.80% | $0.0848 |
| `thermocatalytic_pbr` | PBR | 99.9970% | 0.000 | 0.000 | 100.00% | 253.21 | 74.40% | 31.32% | $0.0772 |
| `thermocatalytic_fluidized` | Fluidized | 99.9970% | 0.000 | 0.000 | 100.00% | 253.21 | 38.32% | 16.13% | $0.1419 |
| `thermocatalytic` | PBR | 99.9970% | 0.000 | 0.000 | 100.00% | 253.21 | 74.40% | 31.32% | $0.0772 |
| `mmbcr` | MMBCR | 99.9967% | 0.000 | 3.198 | 96.23% | 252.12 | 62.61% | 26.25% | $0.0920 |
| `ntec` | NTEC | 99.9970% | 0.000 | 0.000 | 100.00% | 253.21 | 84.77% | 35.69% | $0.0534 |
| `electrochemical` | Electrochemical | 99.9970% | 0.000 | 0.000 | 100.00% | 253.21 | 56.47% | 23.77% | $0.1211 |

#### 18-Feedstock Multi-Hydrocarbon, Waste Polyolefin & Biomass Flexibility Matrix

| Feedstock | Name | Category | H/C | kg_feed/kg_H₂ | kg_C/kg_H₂ | ΔH°₂₉₈ (kJ/mol H₂) | PFR X | PFR Net LCOH ($/kg) | MMBCR X | MMBCR Gross LCOH ($/kg) | MMBCR Net LCOH ($/kg) | MMBCR Opt. LCOH ($/kg) | Post-PSA CO (ppm) | Coupled LCOE ($/kWh_e) |
|-----------|------|----------|-----|---------------|------------|--------------------|-------|---------------------|---------|-------------------------|-----------------------|------------------------|-------------------|------------------------|
| `CH4` | Methane (Pipeline Natural Gas) | `natural_gas` | 4.00 | 3.98 | 2.98 | 37.4 | 56.31% | $1.177 | 47.49% | $2.053 | $1.272 | $0.257 | 0.00 | $0.0920 |
| `ASSOC_GAS` | Associated / Flare Gas (Wet C1-C4 Blend) | `natural_gas` | 3.50 | 4.41 | 3.40 | 23.3 | 80.21% | $0.373 | 79.86% | $1.257 | $0.373 | $-0.375 | 0.00 | $0.0276 |
| `C2H6` | Ethane (Shale NGL) | `light_alkane` | 3.00 | 4.97 | 3.97 | 27.9 | 90.09% | $0.795 | 89.85% | $1.820 | $0.794 | $-0.156 | 0.00 | $0.0579 |
| `C3H8` | Propane (Wholesale LPG) | `light_alkane` | 2.67 | 5.47 | 4.47 | 26.0 | 96.70% | $1.010 | 97.92% | $2.108 | $1.004 | $-0.040 | 0.00 | $0.0736 |
| `C4H10` | n-Butane (Refinery C4 Raffinate / LPG) | `light_alkane` | 2.50 | 5.77 | 4.77 | 25.2 | 98.19% | $0.936 | 99.31% | $2.090 | $0.930 | $-0.121 | 0.00 | $0.0686 |
| `C5H12` | n-Pentane (Natural Gasoline / C5+ Condensate) | `liquid_hydrocarbon` | 2.40 | 5.96 | 4.96 | 24.5 | 98.79% | $0.960 | 99.68% | $2.153 | $0.956 | $-0.120 | 0.00 | $0.0707 |
| `C6H14` | n-Hexane (Straight-Run Light Naphtha) | `liquid_hydrocarbon` | 2.33 | 6.11 | 5.11 | 23.9 | 98.91% | $1.042 | 99.78% | $2.261 | $1.038 | $-0.074 | 0.00 | $0.0768 |
| `C6H12` | Cyclohexane (Naphthenic Ring Cut / LOHC) | `cycloalkane_olefin` | 2.00 | 6.96 | 5.96 | 20.6 | 99.50% | $0.980 | 99.87% | $2.403 | $0.977 | $-0.221 | 0.00 | $0.0727 |
| `C8H18` | n-Octane (Heavy Naphtha / Gasoline Surrogate) | `liquid_hydrocarbon` | 2.25 | 6.30 | 5.30 | 23.2 | 99.47% | $1.127 | 99.88% | $2.383 | $1.124 | $-0.029 | 0.00 | $0.0834 |
| `C12H26` | n-Dodecane (Kerosene / Middle Distillate Diesel) | `liquid_hydrocarbon` | 2.17 | 6.50 | 5.50 | 22.4 | 99.69% | $1.208 | 99.90% | $2.501 | $1.207 | $0.013 | 0.00 | $0.0898 |
| `C16H34` | n-Hexadecane / Cetane (Heavy Gas Oil / Slack Wax) | `liquid_hydrocarbon` | 2.12 | 6.61 | 5.61 | 22.0 | 99.77% | $0.943 | 99.90% | $2.253 | $0.942 | $-0.185 | 0.00 | $0.0704 |
| `C10H20_PE` | Waste Polyolefin Pyrolysis Oil (PE/PP -(CH2)-n) | `waste_plastic` | 2.00 | 6.96 | 5.96 | 16.5 | 99.66% | $-0.043 | 99.90% | $1.363 | $-0.044 | $-0.932 | 0.00 | $0.0017 |
| `BIOMETHANE` | Biomethane / RNG (Upgraded Anaerobic Biogas) | `biomass_biogenic` | 4.00 | 3.98 | 2.98 | 37.4 | 56.31% | $1.063 | 47.49% | $1.938 | $1.157 | $0.175 | 0.00 | $0.0837 |
| `BIOGAS` | Raw Anaerobic Digester / Landfill Biogas (CH4 + CO2) | `biomass_biogenic` | 3.00 | 7.62 | 1.99 | 45.3 | 76.43% | $1.043 | 67.44% | $1.607 | $1.085 | $0.273 | 50.00 | $0.1453 |
| `C18H34O2_LIPID` | Waste Cooking Oil / Tall Oil Bio-Lipid (C18H34O2) | `biomass_biogenic` | 1.89 | 8.24 | 5.61 | 20.5 | 99.68% | $0.466 | 99.90% | $1.776 | $0.465 | $-0.523 | 22.52 | $0.0449 |
| `C3H8O3_GLYCEROL` | Crude Biodiesel Glycerol Byproduct (C3H8O3) | `biomass_biogenic` | 2.67 | 9.44 | 1.49 | 34.5 | 99.89% | $0.970 | 99.83% | $1.339 | $0.970 | $0.313 | 50.00 | $0.1309 |
| `C6H10O2_BIOOIL` | Fast-Pyrolysis Bio-Oil (Hydro-Deoxygenated C6H10O2) | `biomass_biogenic` | 1.67 | 11.32 | 4.77 | 24.4 | 99.41% | $0.459 | 99.78% | $1.606 | $0.458 | $-0.453 | 50.00 | $0.0648 |
| `C6H10O4_BIOMASS` | Torrefied Lignocellulosic Biomass (Woody/Ag C6H10O4) | `biomass_biogenic` | 1.67 | 14.50 | 2.38 | 23.6 | 99.55% | $0.874 | 99.53% | $1.457 | $0.874 | $0.148 | 50.00 | $0.1190 |

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
