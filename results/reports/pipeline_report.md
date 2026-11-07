# Turquoise Hydrogen → Fuel Cell: Pipeline Report

**Generated**: 2026-11-07 01:01:26
**Pipeline Version**: 2.0.0

> **Evidence status:** screening values are hypotheses. Mock results and unconverged calculations are excluded from validated claims; catalyst performance requires experimental confirmation.

## Executive Summary

| Metric | Value |
|--------|-------|
| Catalysts evaluated (Fairchem) | 84 |
| Valid candidates | 58 |
| Pareto front size | 11 |
| Best eligible screening barrier | 0.275 eV |
| Best coking resistance | 6.6712 |
| Reactor simulations | 20 |
| Solids judge CH₄ conversion (ni_np_lit PFR) | 26.34% |
| Top candidate solids CH₄ conversion (Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba]) | 96.58% |
| MMBCR max CH₄ conversion | 3.35% (not a solids rank) |
| NTEC max CH₄ conversion (coupled MOOSE+Cantera) | 1.34% |
| Electrochemical max CH₄ conversion (coupled MOOSE+Cantera) | 0.67% |
| FC catalysts evaluated | 104 |
| Best PEMFC power density | 2.253 W/cm² |
| Best PEMFC efficiency | 45.0% |
| Min ORR overpotential | 0.159 V |

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

Solids judged on **ni_np_lit** (single-pass X, `a`, WHSV, ΔP). MMBCR X_eq and 0.01 eV H-parked cats are not ranks.

| Reactor | Catalyst | T (K) | Single-pass X | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |
|---------|----------|-------|---------------|---------|------------|----------|
| PFR | ni_np_lit | 973.15 | 26.34% | 4154 | 900 | 0.33 |
| PBR | ni_np_lit | 973.15 | 23.03% | 4154 | 1272 | 0.50 |
| Fluidized | ni_np_lit | 973.15 | 18.79% | 3808 | 400 | 0.33 |

Top candidate solids X = 96.58% on **Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba]** (PBR @ 1300.0 K)

MMBCR max X = 3.4% (X_eq by construction for large Da; not a catalyst rank (B5))

### Solids runs (single-pass)

| Catalyst | Reactor | T (K) | Single-pass X | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |
|----------|---------|-------|---------------|---------|------------|----------|
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PBR | 1300.0 | 97.09% | 4154 | 2783 | 0.91 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PFR | 1300.0 | 97.05% | 4154 | 900 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PBR | 1300.0 | 96.58% | 4154 | 2956 | 0.91 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PFR | 1300.0 | 94.01% | 4154 | 900 | 0.40 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PFR | 1100.0 | 72.54% | 4154 | 900 | 0.36 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PBR | 1100.0 | 60.18% | 4154 | 1696 | 0.84 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | Fluidized | 1300.0 | 48.45% | 3808 | 400 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PFR | 1100.0 | 44.65% | 4154 | 900 | 0.36 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | Fluidized | 1100.0 | 43.84% | 3808 | 400 | 0.36 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | Fluidized | 1300.0 | 41.27% | 3808 | 400 | 0.40 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | PBR | 1100.0 | 40.05% | 4154 | 1572 | 0.73 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | Fluidized | 1100.0 | 27.35% | 3808 | 400 | 0.36 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PFR | 973.15 | 26.34% | 4154 | 900 | 0.33 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PBR | 973.15 | 23.03% | 4154 | 1272 | 0.50 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | Fluidized | 973.15 | 18.79% | 3808 | 400 | 0.33 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PFR | 923.15 | 14.56% | 4154 | 900 | 0.32 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PBR | 923.15 | 13.04% | 4154 | 1185 | 0.44 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | Fluidized | 923.15 | 10.70% | 3808 | 400 | 0.32 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PFR | 900.0 | 10.69% | 4154 | 900 | 0.32 |
| Supported Extended-Surface Transition Metal Catalyst (Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021) [id=ni_np_lit] | PBR | 900.0 | 9.66% | 4154 | 1155 | 0.42 |

### Liquid-Metal, NTEC, and Electrochemical Runs

| Catalyst | Reactor | T (K) | CH₄ Conversion | H₂ Selectivity | Evidence Tier |
|----------|---------|-------|----------------|----------------|---------------|
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MMBCR | 1300.0 | 3.35% | 96.92% | diagnostic_screening_template |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | NTEC | 473.15 | 1.34% | 100.00% | uncalibrated_multiphysics_screening |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | NTEC | 473.15 | 1.16% | 100.00% | uncalibrated_multiphysics_screening |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | NTEC | 423.15 | 1.08% | 100.00% | uncalibrated_multiphysics_screening |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MMBCR | 1100.0 | 0.84% | 0.00% | diagnostic_screening_template |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | NTEC | 473.15 | 0.79% | 100.00% | uncalibrated_multiphysics_screening |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | Electrochemical | 353.15 | 0.67% | 100.00% | uncalibrated_multiphysics_screening |
| Cu,Eu-doped Pr(bcc110)/BaO, strain=+0.07% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=00014a7eeba8246b1b9c] | Electrochemical | 353.15 | 0.60% | 100.00% | uncalibrated_multiphysics_screening |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | Electrochemical | 353.15 | 0.60% | 100.00% | uncalibrated_multiphysics_screening |
| Cu,Eu-doped Pr(bcc110)/BaO, strain=+0.07% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=00014a7eeba8246b1b9c] | NTEC | 423.15 | 0.57% | 100.00% | uncalibrated_multiphysics_screening |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | Electrochemical | 823.15 | 0.48% | 100.00% | uncalibrated_multiphysics_screening |
| Nd-N4_pyridine-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=8609ad5ba35faea2ac57] | Electrochemical | 823.15 | 0.48% | 100.00% | uncalibrated_multiphysics_screening |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | Electrochemical | 823.15 | 0.48% | 100.00% | uncalibrated_multiphysics_screening |
| ZnAl-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=0eb155d20835aac7b291] | Electrochemical | 823.15 | 0.44% | 100.00% | uncalibrated_multiphysics_screening |
| ZnAl-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=21724d09776007929597] | Electrochemical | 823.15 | 0.44% | 100.00% | uncalibrated_multiphysics_screening |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | Electrochemical | 353.15 | 0.42% | 100.00% | uncalibrated_multiphysics_screening |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | Electrochemical | 353.15 | 0.40% | 100.00% | uncalibrated_multiphysics_screening |
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | NTEC | 423.15 | 0.35% | 100.00% | uncalibrated_multiphysics_screening |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | NTEC | 423.15 | 0.31% | 100.00% | uncalibrated_multiphysics_screening |
| ZnAl-N2S2N2/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=21a2125a111a954bd9a7] | Electrochemical | 823.15 | 0.31% | 100.00% | uncalibrated_multiphysics_screening |

Excluded 324 non-run reactor files (sweeps, mocks, not_applicable, validation_required, failed) from performance claims.

## Phase 3: DFT Validation (Quantum ESPRESSO)

| Catalyst | Material Class | Energy (eV) | Converged | Max Force |
|----------|----------------|-------------|-----------|-----------|
| CrY-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=25810913ec6523464e81] | DAC | -7684.891466 | True | 0.0 |
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | MAXPhase | -5835.70603 | True | 0.0 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | MetalFreeCarbon | -612.2664572547036 | True | 0.0 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | SupportedNanoparticle | -16117.69086882 | True | 0.0 |
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
| Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | -26244.2822 | +1.2000 | 2.1945 | -0.2000 | True |
| Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | -10492.2452 | +0.7838 | 1.6164 | +4.6052 | True |

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
- **Mean Potential Energy**: 8.1282 ± 0.2164 eV
- **Active Promoter Vertical Offset ($\Delta z_{\text{promoter-host}}$)**: -0.1412 Å (Gibbs surface-tension sub-surface relaxation)

## Phase 4: Quantum VQE Transition States

| Catalyst | Reaction | Energy (Ha) | Energy (eV) | Qubits | Layers |
|----------|----------|-------------|-------------|--------|--------|
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | CH_split | -208.197442 | -5665.3410 | 8 | 3 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | CH_split | -94.606711 | -2574.3797 | 8 | 3 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | CH_split | -208.192487 | -5665.2062 | 8 | 3 |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | ORR | -205.049811 | -5579.6896 | 8 | 3 |
| Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | CH_split | -158.278058 | -4306.9654 | 8 | 3 |
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | ORR | -268.549084 | -7307.5929 | 8 | 3 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | CH_split | -184.646811 | -5024.4957 | 8 | 3 |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | CH_split | -97.539149 | -2654.1755 | 8 | 3 |

### Multi-Reference `CASSCF(4,4)` Active-Space Diagnostics (`lanl2dz+6-31g*`, 8 Qubits vs. Single-Reference DFT)

| Catalyst | Active TS Motif | ΔE_corr,MR (eV) | Leading CI |c₀|² | CCSD T₁ Diag. | E_act,DFT-PBE (eV) | E_act,VQE-CAS (eV) | ΔΔE‡_MR-DFT (eV) |
|----------|-----------------|-----------------|------------------|---------------|--------------------|--------------------|------------------|
| (Co,Ni)3AlC2(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=3def15c30d9a9b799254] | `Ni-CH4_ts` | -0.4010 | 0.9623 | 0.0827 | 0.0800 | 2.2182 | +2.1382 |
| pyrrolic,0.04-doped none graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=5ed9c863d22330991eb6] | `NH-CH4_ts` | -3.6095 | 0.4954 | 0.0136 | 7.0071 | 9.1224 | +2.1153 |
| Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | `Ni-CH4_ts` | -0.4022 | 0.9622 | 0.0834 | 0.0800 | 2.3560 | +2.2760 |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | `NH-O2H2_ts` | -5.2236 | 0.4836 | 0.1096 | 0.9026 | 4.1749 | +3.2723 |
| Pt1/Ni(110) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=b6702ffe2074faa4fa50] | `Pt-CH4_ts` | -0.4761 | 0.9728 | 0.0295 | 0.6459 | 2.1062 | +1.4602 |
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | `Pt-O2H2_ts` | -5.6491 | 0.4468 | 0.0935 | 0.0800 | 0.1201 | +0.0401 |
| Ce-doped MgCo2O4/N-graphene (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=c9bf4b755ab0023fffba] | `CoH-CH4_ts` | -1.8800 | 0.9609 | 0.0265 | 0.0800 | 0.7899 | +0.7099 |
| (Ta,Ti)2AlC1(basal_0001) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000) [id=eb64b3aa6d3913c21f5d] | `Ti-CH4_ts` | -6.3827 | 0.4750 | 0.1651 | 1.0431 | 0.1153 | -0.9278 |

## Phase 5: PEMFC Fuel Cell Performance

### Single-Cell Performance

| Cathode | Membrane | Peak Power (W/cm²) | OCV (V) | η (V) | Efficiency (Rated) |
|---------|----------|---------------------|---------|-------|--------------------|
| W-doped NiCo2O4/N-CNT (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=dbc15731f8fdb0f99d43] | AquivionE87 | 2.2532 | 1.179 | 0.159087 | 45.0% |
| Pt1/Co(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=c77b56fa00cdc515e3b7] | AquivionE87 | 2.3119 | 1.179 | 0.2069 | 46.1% |
| oxidized,0.01-doped edge_armchair graphene Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009) [id=87d86cdaf36fda7c6831] | AquivionE87 | 2.0398 | 1.179 | 0.172376 | 40.7% |
| Mn,Cu-doped Au(fcc110)/NiP, strain=+0.06% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=c5ad001debb359e04c33] | AquivionE87 | 1.7496 | 1.179 | 0.349075 | 34.9% |
| Ag1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=730c6d4349990c81e975] | AquivionE87 | 1.7441 | 1.179 | 0.348847 | 34.8% |
| Pt-Boroxine-O4 (10.0) COF Network (Furukawa et al., Science 2013) [id=7e36da26795823ef060e] | AquivionE87 | 1.6168 | 1.179 | 0.342209 | 32.3% |
| Pd1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=3ab787b3c8004744ae8e] | AquivionE87 | 1.6628 | 1.179 | 0.36916 | 33.2% |
| Pd1/Au(100) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=f353706418ed8a755f1f] | AquivionE87 | 1.6072 | 1.179 | 0.383056 | 32.1% |
| BaSbO3(B_vacancy) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide) [id=06d7172dc944ea7d38c7] | AquivionE87 | 1.4412 | 1.179 | 0.354082 | 28.8% |
| N,Ta-doped Au(bcc211)/MoO3, strain=+0.10% Supported Solid Catalyst (Bengaard 2002 / Alves 2021) [id=7ec3c045b95ab35fc877] | AquivionE87 | 1.4751 | 1.179 | 0.417691 | 29.4% |
| PtSn-N8/Graphene Dual-Atom Catalyst (Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site) [id=7fba4d6b4f15772edf2d] | AquivionE87 | 1.4187 | 1.179 | 0.420581 | 28.3% |
| Co-Imine-N2P2 (6.0) COF Network (Furukawa et al., Science 2013) [id=3cafa1dd27033ac418d3] | AquivionE87 | 1.2629 | 1.179 | 0.430677 | 25.2% |
| Au-Oxalate-N3P (16.0) MOF Network (Furukawa et al., Science 2013) [id=56c74d76175883fef1fb] | AquivionE87 | 1.2549 | 1.179 | 0.440691 | 25.0% |
| Pt1/Ag(111) Single-Atom Alloy (Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018) [id=7a83fee91102006b288b] | AquivionE87 | 1.1631 | 1.179 | 0.494092 | 23.2% |
| Ce-doped CuCu2O4/graphene_oxide (nanoparticle) Spinel Oxide (Liang et al., Nat. Mater. 2011) [id=cef5b9bdd5a3145b4e70] | AquivionE87 | 0.7207 | 1.179 | 0.553482 | 18.5% |

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
| `thermocatalytic` | PBR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 96.58% | 54.87% | 0.993 | 0.913 | $1.188 | $0.198 | `dft_neb_complete` |
| `thermocatalytic_pbr` | PBR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 96.58% | 54.87% | 0.993 | 0.913 | $1.188 | $0.198 | `dft_neb_complete` |
| `thermocatalytic_pfr` | PFR | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 94.01% | 46.32% | — | 0.399 | $1.265 | $0.240 | `dft_neb_complete` |
| `thermocatalytic_fluidized` | Fluidized | Cu,Fe-doped Ni(fcc211)/MgAl2O4 (d=12.0 nm, 25.0 wt%) Supported Bimetallic Nanoparticle (Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy) [id=768cb0559591faec06ba] | 1300.0 | 41.27% | 11.10% | — | 0.399 | $2.845 | $1.095 | `dft_neb_complete` |
| `mmbcr` | MMBCR | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | 1300.0 | 3.35% | 1.21% | — | — | $19.830 | $10.287 | `dft_neb_complete` |
| `ntec` | NTEC | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | 473.1 | 1.34% | 0.36% | — | — | $52.902 | $31.769 | `dft_neb_complete` |
| `electrochemical` | Electrochemical | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | 823.1 | 0.48% | 0.14% | — | — | $321.607 | $190.625 | `dft_neb_complete` |

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
| `thermocatalytic_pfr` | $1.265 | $0.240 | $0.0902 | $0.0171 |
| `thermocatalytic_pbr` | $1.188 | $0.198 | $0.0847 | $0.0141 |
| `thermocatalytic_fluidized` | $2.845 | $1.095 | $0.2028 | $0.0780 |
| `thermocatalytic` | $1.188 | $0.198 | $0.0847 | $0.0141 |
| `mmbcr` | $19.830 | $10.287 | $1.4132 | $0.7331 |
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
| 8 | Ga99.9Ti0.1 Liquid Metal Alloy (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt) [id=c4ab2249c156880efdbe] | MoltenMetal | -26244.28 | +1.2000 | 2.1945 | -0.2000 | `True` |
| 9 | Ta-N4-N-graphene/none Single-Atom Catalyst (Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019) [id=f58366bf5dfcdc4578f3] | SAC | -10492.25 | +0.7838 | 1.6164 | +4.6052 | `True` |

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
- **Mean Potential Energy**: 8.1282 ± 0.2164 eV
- **Gibbs Promoter Vertical Offset ($\Delta z_{\text{promoter}}$)**: -0.1412 Å

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
