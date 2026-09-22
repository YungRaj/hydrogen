# ADR 0001 — Phase stability at application temperature

- **Status:** Accepted
- **Date:** 2026-09-06
- **Deciders:** pyrolysis screening / Phase 2 reactor handoff

## Context

Phase 2 methane pyrolysis in this repo is **773–1300 K**, 1 bar, flowing CH₄-rich feed. Genomes encode a **named phase** (hydride family, MOF linker/node, MXene termination, perovskite stoichiometry, …). Ranking that genome as a reactor catalyst is only valid if that phase still exists under those conditions.

A one-class exception (`MetalHydride` only) is the weakest possible position: it is the class that had been winning on `E_act`, so a selective filter looks like score-shopping. The same costume argument applies at least as strongly to MOF, COF, MXene, and perovskite.

## Decision

Admissibility for a turquoise-pyrolysis reactor is a **generalized** check:

`phase_stable_at_application_T(genome, T_K, application='turquoise_pyrolysis')`

A class is `out_of_scope` when the **encoded phase** is not the material present at application T. Residual metal, carbonized debris, TiC, or exsolved nanoparticles may be interesting — they must be encoded as the class they actually are (`SolidCatalyst`, `MoltenMetal`, `MetalFreeCarbon`, …).

MetalHydride is the **first** documented entry, not the only one.

The check is applied at **candidate selection** for reactors (and any later DFT slate drawn from that list). It does **not** change the indexed-space coverage certificate or the 21,092,645,031 denominator. Those count enumerated encodings, not application-admissible phases.

### Arguments that hold across the whole 773–1300 K band

Do not lead with “hundreds of degrees hotter.” That is true at 1300 K and arguable at 773 K for the most stable interstitial hydrides (ZrH₁.₅–₁.₆ is a TRIGA moderator near 600–750 °C). Across the **entire** sweep:

1. **P_H₂ ≈ 0 in flowing CH₄.** Van ’t Hoff equilibrium hydrogen content of a hydride collapses when the gas-phase H₂ partial pressure is near zero, regardless of how stable the hydride is under 1 bar H₂.
2. **Carburization.** A Zr or Hf (or similar) surface under CH₄ at 500–1000 °C forms carbide. The working solid is not the encoded hydride.

Temperature margin is supporting evidence at the top of the band only.

### Per-class basis (pyrolysis)

| Class | Status | Basis |
|---|---|---|
| **MetalHydride** | out_of_scope | Encoded hydride phase. Storage families unload well below the band under 1 bar H₂ (AB5 already ~1.7 bar at 298 K [13]; MgH₂ TDS ~410 °C [14]; TiH₂ TPD completes ≲ 650 °C [15,16]). At 773 K the ZrHₓ case is the honest weak point — close it with P_H₂ ≈ 0 and carburization, not ΔT. Residual metal ≠ hydride genome. |
| **MOF** | out_of_scope | Carboxylates / imidazolates lose crystallinity typically 300–400 °C; the most robust (MIL-53, UiO-66) are cited around 500 °C [17,18]. Above ~550 °C in inert gas they pyrolyze to carbonized MOF-derived composites [17]. At 900–1300 K under CH₄ the encoded framework is gone. Precursor, not catalyst. |
| **COF** | out_of_scope | Same precursor-not-catalyst case: organic frameworks do not survive 773–1300 K as the encoded COF [19]. What remains is a carbonaceous residue. |
| **MXene** | out_of_scope | −OH / −F / −O terminations that define Tₓ leave well below the top of the band; Ti₃C₂Tₓ converts toward TiC (and oxides if any O is present) [20,21]. Encoded MXene ≠ high-T carbide. |
| **Perovskite** | out_of_scope | Under reducing CH₄ at ~500–1000 °C, A-site-deficient perovskites **exsolve** metal nanoparticles [22,23]. The working catalyst is socketed metal + a reduced host, not the encoded ABO₃ genome. |
| MoltenMetal | candidate | Liquid metal / alloy is the intended Phase 2 phase (Upham-type melt [5]). |
| SolidCatalyst, SAC, DAC, SAA, HEA, Spinel, MAXPhase, MetalFreeCarbon | candidate | Encoded as high-T solids or carbons. Revisit if a specific subclass is later shown to reconstruct. |

PEMFC cathode scope is unchanged and separate: `MetalHydride` and `MoltenMetal` have no solid catalyst-layer realization.

## Consequences

### Candidate selection

`select_turquoise_pyrolysis_candidates()` uses this check. Unparseable genomes fail closed. Phase 1→2 and Phase-2-only CSV load both re-select from the full valid pool (no `dir()` bypass).

### Coverage certificate (unchanged denominator)

All 14 classes remain in the indexed space. Branch class-floor and the 21.1B population count still include MetalHydride (796,068 encodings) and the other out-of-scope classes. The filter does not rewrite identity.

### Validation / DFT budget

MetalHydride cannot produce a result under **either** application (PEMFC cathode already out of scope; pyrolysis now out of scope). It remains enumerated for coverage and still receives `--branch-class-floor` / champion-archive slots, but it is **exempt from `--min-validation-per-class` reserved slots** on both applications so eSen and DFT quota is not spent on a dead class. It may still appear in a remainder allocation only if a caller passes it in; reactor and PEMFC slates should not.

MOF, COF, MXene, and Perovskite still earn PEMFC reserved slots (they are in-scope there). They do **not** earn reserved pyrolysis validation slots, because the encoded phase cannot produce a pyrolysis reactor result. `--validation-batch` floors use quota-eligible classes for that application (pyrolysis: 14 − 5; PEMFC: 14 − 1), not 14.

### What this is not

- Not a slab-coking rule. `coking_index = ΔE_C − 2ΔE_H` stays defined for any class with a slab, including MetalHydride. That index is not the filament vs encapsulating branch (B6).
- Not a claim that residual metal, MOF-derived carbon, TiC, or exsolved Ni is inactive. Those need their own genomes.
- Not a grant of the Baker/Helveg filament path. SAC remaining a pyrolysis **candidate** means the encoded single-atom phase can exist at application T. It does not mean Cα can dissolve, diffuse, and nucleate graphite at a separate face [24–26, 33]. That physics is B6 and is nanoparticle-only.

## References

1. Alves, L.; Pereira, V.; Lagarteira, T.; Mendes, A. Catalytic methane decomposition to boost the energy transition: Scientific and technological advancements. *Renew. Sustain. Energy Rev.* **2021**, *137*, 110465. DOI: [10.1016/j.rser.2020.110465](https://doi.org/10.1016/j.rser.2020.110465).
2. Sánchez-Bastardo, N.; Schlögl, R.; Ruland, H. Methane Pyrolysis for Zero-Emission Hydrogen Production: A Potential Bridge Technology from Fossil Fuels to a Renewable and Sustainable Hydrogen Economy. *Ind. Eng. Chem. Res.* **2021**, *60*, 11855–11881. DOI: [10.1021/acs.iecr.1c01679](https://doi.org/10.1021/acs.iecr.1c01679).
3. Gili, A.; et al. Toward Rational Design of Nickel Catalysts for Thermocatalytic Decomposition of Methane. *ChemCatChem* **2024**. DOI: [10.1002/cctc.202301629](https://doi.org/10.1002/cctc.202301629).
4. Pinilla, J. L.; et al. Hydrogen production by methane cracking using Ni-supported catalysts in a fluidized bed. *Int. J. Hydrogen Energy* **2012**.
5. Upham, D. C.; Agarwal, V.; Khechfe, A.; Snodgrass, Z. R.; Gordon, M. J.; Metiu, H.; McFarland, E. W. Catalytic molten metals for the direct conversion of methane to hydrogen and separable carbon. *Science* **2017**, *358*, 917–921. DOI: [10.1126/science.aao5023](https://doi.org/10.1126/science.aao5023).
6. Abdollahi, M. R.; Nathan, G. J.; Jafarian, M. Preliminary Evaluation of Methods for Continuous Carbon Removal from a Molten Catalyst Bubbling Methane Pyrolysis Reactor. *Energies* **2024**, *17*, 290. DOI: [10.3390/en17020290](https://doi.org/10.3390/en17020290).
7. Chen, L.; et al. Ternary NiMo-Bi liquid alloy catalyst for efficient hydrogen production from methane pyrolysis. *Science* **2023**. DOI: [10.1126/science.adh8872](https://doi.org/10.1126/science.adh8872).
8. Meloni, E.; et al. Methane pyrolysis in molten media: The interplay of physical properties and catalytic activity on carbon and hydrogen production. *J. Anal. Appl. Pyrolysis* **2024**.
9. Le, B. T.; Ngo, S. I.; Lim, Y.-I.; Lee, U.-D. One-dimensional kinetic model with novel bubble size equation in molten-metal bubble column reactors for CH₄ pyrolysis. *AIChE J.* **2024**. DOI: [10.1002/aic.18540](https://doi.org/10.1002/aic.18540).
10. Rahimi, N.; et al. Solid carbon production and recovery from high temperature methane pyrolysis in bubble columns containing molten metals and molten salts. *Carbon* **2019**, *151*, 181–191. DOI: [10.1016/j.carbon.2019.05.041](https://doi.org/10.1016/j.carbon.2019.05.041).
11. Abbas, H. F.; Daud, W. M. A. W. Hydrogen production by thermocatalytic decomposition of methane using a fixed bed activated carbon in a pilot scale unit. *Int. J. Hydrogen Energy* **2010**.
12. Scheiblehner, D.; et al. Hydrogen Production by Methane Pyrolysis in Molten Cu-Ni-Sn Alloys. *Metals* **2023**, *13*, 1310. DOI: [10.3390/met13071310](https://doi.org/10.3390/met13071310).
13. Joubert, J.-M.; et al. LaNi₅ related AB₅ compounds: Structure, properties and applications. *J. Alloys Compd.* **2021**, *862*, 158163. DOI: [10.1016/j.jallcom.2020.158163](https://doi.org/10.1016/j.jallcom.2020.158163).
14. El-Eskandarany, M. S.; et al. Reversible Hydrogen Storage Using Nanocomposites. *Appl. Sci.* **2020**, *10*, 4618. DOI: [10.3390/app10134618](https://doi.org/10.3390/app10134618).
15. Liu, H.; et al. Phase transformations of titanium hydride in thermal desorption process with different heating rates. *Int. J. Hydrogen Energy* **2015**.
16. Matz, N.; et al. Titanium hydride elevated-temperature literature review. OSTI 1822588, 2021. https://www.osti.gov/servlets/purl/1822588
17. Hupp, J. T.; Farha, O. K.; et al. Carbon-efficient conversion of natural gas … via catalytic metal–organic framework (MOF) chemistry. *Energy Environ. Sci.* **2022**. DOI: [10.1039/d2ee01010k](https://doi.org/10.1039/d2ee01010k). (MOFs complementary to high-T catalysts; pyrolyze ≳ 550 °C to carbonized composites.)
18. Howarth, A. J.; Liu, Y.; Li, P.; Li, Z.; Wang, T. C.; Hupp, J. T.; Farha, O. K. Chemical, thermal and mechanical stabilities of metal–organic frameworks. *Chem. Soc. Rev.* **2016**, *45*, 5107–5134.
19. Kandambeth, S.; Dey, K.; Banerjee, R. Covalent Organic Frameworks: Chemistry beyond the Structure. *J. Am. Chem. Soc.* **2019**, *141*, 1807–1822.
20. Naguib, M.; Kurtoglu, M.; Presser, V.; Lu, J.; Niu, J.; Heon, M.; Hultman, L.; Gogotsi, Y.; Barsoum, M. W. Two-Dimensional Nanocrystals Produced by Exfoliation of Ti₃AlC₂. *Adv. Mater.* **2011**, *23*, 4248–4253.
21. Seredych, M.; et al. High-Temperature Behavior and Surface Chemistry of Ti₃C₂Tₓ MXenes. *Chem. Mater.* (thermal loss of Tₓ; conversion toward TiC).
22. Neagu, D.; Tsekouras, G.; Miller, D. N.; Ménard, H.; Irvine, J. T. S. In situ growth of nanoparticles through control of non-stoichiometry. *Nat. Chem.* **2013**, *5*, 916–923.
23. Kousi, K.; Neagu, D.; Bekris, L.; Papaioannou, E. I.; Metcalfe, I. S. Endogenous nanoparticles strain perovskite host lattice providing oxygen capacity and driving oxygen exchange and CH₄ conversion. *Angew. Chem. Int. Ed.* **2020**, *59*, 2510–2519.
24. Baker, R. T. K.; Barber, M. A.; Harris, P. S.; Feates, F. S.; Waite, R. J. Nucleation and growth of carbon deposits from the nickel catalyzed decomposition of acetylene. *J. Catal.* **1972**, *26*, 51–62. DOI: [10.1016/0021-9517(72)90032-2](https://doi.org/10.1016/0021-9517(72)90032-2).
25. Helveg, S.; López-Cartes, C.; Sehested, J.; Hansen, P. L.; Clausen, B. S.; Rostrup-Nielsen, J. R.; Abild-Pedersen, F.; Nørskov, J. K. Atomic-scale imaging of carbon nanofibre growth. *Nature* **2004**, *427*, 426–429. DOI: [10.1038/nature02278](https://doi.org/10.1038/nature02278).
26. Abild-Pedersen, F.; Nørskov, J. K.; Rostrup-Nielsen, J. R.; Sehested, J.; Helveg, S. Mechanisms for catalytic carbon nanofiber growth studied by *ab initio* density functional theory calculations. *Phys. Rev. B* **2006**, *73*, 115419. DOI: [10.1103/PhysRevB.73.115419](https://doi.org/10.1103/PhysRevB.73.115419).
27. Snoeck, J.-W.; Froment, G. F.; Fowles, M. Filamentous carbon formation and gasification: thermodynamics, driving force, nucleation, and steady-state growth. *J. Catal.* **1997**, *169*, 240–249. DOI: [10.1006/jcat.1997.1634](https://doi.org/10.1006/jcat.1997.1634).
28. Snoeck, J.-W.; Froment, G. F.; Fowles, M. Kinetic study of the carbon filament formation by methane cracking on a nickel catalyst. *J. Catal.* **1997**, *169*, 250–262. DOI: [10.1006/jcat.1997.1635](https://doi.org/10.1006/jcat.1997.1635).
29. Xu, M.; Lopez-Ruiz, J. A.; Kovarik, L.; et al. Structure sensitivity and its effect on methane turnover and carbon co-product selectivity in thermocatalytic decomposition of methane over supported Ni catalysts. *Appl. Catal. A* **2021**, *611*, 117967. DOI: [10.1016/j.apcata.2020.117967](https://doi.org/10.1016/j.apcata.2020.117967).
30. Gili, A.; Schlicker, L.; Bekheet, M. F.; Görke, O.; Kober, D.; Simon, U.; Littlewood, P.; Schomäcker, R.; Doran, A.; Gaissmaier, D.; Jacob, T.; Selve, S.; Gurlo, A. Revealing the mechanism of multiwalled carbon nanotube growth on supported nickel nanoparticles by in situ synchrotron X-ray diffraction, density functional theory, and molecular dynamics simulations. *ACS Catal.* **2019**, *9*, 6999–7011. DOI: [10.1021/acscatal.9b00733](https://doi.org/10.1021/acscatal.9b00733).
31. Ermakova, M. A.; Ermakov, D. Yu. Ni/SiO₂ and Fe/SiO₂ catalysts for production of hydrogen and filamentous carbon via methane decomposition. *Catal. Today* **2002**, *77*, 225–235. DOI: [10.1016/S0920-5861(02)00248-1](https://doi.org/10.1016/S0920-5861(02)00248-1).
32. Amin, A.; Epling, W. S.; Croiset, E. Reaction and deactivation rates of methane catalytic cracking over nickel. *Ind. Eng. Chem. Res.* **2011**, *50*, 12460–12470. DOI: [10.1021/ie201194z](https://doi.org/10.1021/ie201194z).
33. Akri, M.; et al. Atomically dispersed nickel as coke-resistant active sites for methane dry reforming. *Nat. Commun.* **2019**, *10*, 5181. DOI: [10.1038/s41467-019-12843-w](https://doi.org/10.1038/s41467-019-12843-w).
34. Takenaka, S.; Ogihara, H.; Yamanaka, I.; Otsuka, K. Decomposition of methane over supported-Ni catalysts: effects of the supports on the catalytic lifetime. *Appl. Catal. A* **2001**, *217*, 137–146. DOI: [10.1016/S0926-860X(01)00593-2](https://doi.org/10.1016/S0926-860X(01)00593-2).
35. Bengaard, H. S.; Nørskov, J. K.; Sehested, J.; Clausen, B. S.; Nielsen, L. P.; Molenbroek, A. M.; Rostrup-Nielsen, J. R. Steam reforming and graphite formation on Ni catalysts. *J. Catal.* **2002**, *209*, 365–384. DOI: [10.1006/jcat.2002.3579](https://doi.org/10.1006/jcat.2002.3579). Ni(111) terrace CH₄ TS 101 kJ/mol (1.05 eV, no ZPE); recorded judge `E_act` is 1.00 eV.
36. Watwe, R. M.; Bengaard, H. S.; Rostrup-Nielsen, J. R.; Dumesic, J. A.; Nørskov, J. K. Theoretical studies of stability and reactivity of CHₓ species on Ni(111). *J. Catal.* **2000**, *189*, 16–30. DOI: [10.1006/jcat.1999.2699](https://doi.org/10.1006/jcat.1999.2699). Ni(111) fcc CH₃\* vs the CH₃ radical (screener window; typical later PW91 −1.84 eV).
37. Greeley, J.; Mavrikakis, M. Surface and subsurface hydrogen: adsorption properties on transition metals and near-surface alloys. *J. Phys. Chem. B* **2005**, *109*, 3460–3471. DOI: [10.1021/jp046540q](https://doi.org/10.1021/jp046540q). Ni(111) TPD ≈ −0.50 eV vs ½ H₂; PW91 DFT −0.61 eV. Dedicated Ni(111) H study: Greeley, J.; Mavrikakis, M. *Surf. Sci.* **2003**, *540*, 215–229. DOI: [10.1016/S0039-6028(03)00790-8](https://doi.org/10.1016/S0039-6028(03)00790-8).
38. Lander, J. J.; Kern, H. E.; Beach, A. L. Solubility and diffusion coefficient of carbon in nickel. *J. Appl. Phys.* **1952**, *23*, 1305–1309. DOI: [10.1063/1.1702064](https://doi.org/10.1063/1.1702064). ln *D* = 0.909 − 20200/*T* (*D* in cm²/s) → *D*₀ = 2.48 cm²/s = 2.48×10⁻⁴ m²/s.
