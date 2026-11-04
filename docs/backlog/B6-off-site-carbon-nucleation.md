# Baker (1972) / Helveg (2004) / Ermakova (2002) / Takenaka (2003) Off-Site Carbon Filament Nucleation vs. Encapsulation (`baker1972_helveg2004_ermakova2000_takenaka2003_off_site_carbon_nucleation`)

**Type:** Literature Validation Gate & Kinetic Implementation. Do not ship a universal `C_s => C(gr) + site` for isolated single-atom catalysts.

**Status:** **Closed** (`supported_ni_filament_closure_passed = True` in [`results/evidence/supported_ni_filament_closure.json`](../../results/evidence/supported_ni_filament_closure.json)). Gated block in `write_full_mechanism`. Ungated YAML ends at `C_s`. Nanoparticle Ni/Fe/Co (`SolidCatalyst`, `SupportedNanoparticle`, `HEA`, `SAA` host) emit `Cγ` (`C_s => C(gr) + site`, 1.5 eV, `A_γ` declared or derived from Lander bulk diffusion $D_0/L^2$) and `Cδ` (`C_s => C_encap_s`, 1.53 eV, **$\propto \theta_{\text{C}}^2$**, $A_\delta = A_\gamma/\theta^*$, $\theta^* = 0.5$ declared). Not mapped to `coking_index`. Named solids reference is `ni111_sio2_bengaard2002_ref` (supported $\text{Ni(111)/SiO}_2$) at the 650–700 °C filament ROI.

**Depends on:** [`deutschmann2000_monolayer_site_inventory`](B1-solid-site-inventory.md) (done; Γ locked). **Unlocks:** [`alves2021_xu2021_supported_ni_filament_ranking_authority`](B5-phase2-x-gate-before-dft.md). **Complements:** [`muradov2005_zero_co2_solids_outfeed_regeneration`](B2-solids-outfeed-clearing.md) (between-pass outfeed) and [`ergun1952_packed_bed_pressure_drop_pilot_residence`](B4-pfr-residence-time-limits.md) (packed-bed ΔP / τ).

## Why this is the governing physics

The ungated solids path consumes the Langmuir lattice stoichiometrically inside a pass. H₂ can leave (`2 H_s => H2 + 2 site`). Carbon cannot. Real Ni TCD does not die after one monolayer: C leaves the active face, travels through or across the particle, and nucleates graphite at a **different** place (rear face or a step edge). That is what lets Ni run for hours [1, 3, 24–31, 34].

Three identities of an **ungated** mechanism (no off-site step, e.g., `rh1_ceo2_akri2019_sac_ref`):

1. **E_act sweep cannot discriminate on ungated solids.** Once the pass parks C on the available sites, `X ≈ n_sites / n_CH4,fed = (Γ · a · V) / n_CH4`. The barrier is not in the answer.
2. **Inventory gains are linear in `a` by construction.** Monolayer site-inventory (`deutschmann2000_monolayer_site_inventory`) $X \propto a$ ($\text{corr} \approx 1$) is that identity, not kinetic closure.
3. **Melt vs bed is turnovers vs no turnovers**, not continuous-C-removal vs coking. The melt has no site conservation [5]. Between-pass solids outfeed (`muradov2005_zero_co2_solids_outfeed_regeneration`) restores sites **between** passes. The deficit on ungated solids is **within** a pass.

## What the literature actually is (not one elementary hop)

Rostrup-Nielsen / Trimm carbon types: our `C_s` is **Cα** (atomic, site-blocking). Off-site channels:

| Channel | Fate | Site | Literature |
|---|---|---|---|
| Cα → **Cγ** (filament / graphene at a step) | Carbon leaves the terrace; graphite grows elsewhere | **Returned** | Baker [24]; Helveg [25]; Abild-Pedersen [26]; Snoeck [27, 28]; Gili [3, 30] |
| Cα → **Cδ** (encapsulating graphite) | Graphene nucleates on the gas face | **Stays blocked** | Alves [1]; Amin [32]; Xu [29] |

Cycle on a **metal nanoparticle** (not a single atom):

1. CH₄ → Cα on the gas-facing face.
2. C dissolves (Gili: interstitial NiCₓ, lattice expands then contracts [30]) and/or walks on the surface / subsurface.
3. Graphene nucleates at a rear face (Baker tip-growth [24]) or a dynamic Ni step (Helveg [25]; DFT overall barriers 1.42 eV surface / 1.55 eV subsurface to the graphene–Ni interface [26]). Adatom hops on Ni(111) are only ~0.4 eV [26, Hofmann]; the slow step is crossing onto the graphene edge.
4. Filament growth keeps the gas face open. Measured scale: up to **384 gC / gNi** and **4–50 h** on silicate-free high-Ni/SiO₂; silicates drop that to ~40 gC / gNi and ~4 h [31]. Takenaka: crystallized Ni metal on SiO₂ / TiO₂ / graphite lives; Ni locked as oxide or compound with Al₂O₃ / MgO does not [34].

Snoeck [27, 28]: driving force is a **concentration** gradient (different C solubility at gas/Ni vs Ni/filament), not Baker’s original T-gradient. How many filaments nucleate depends on carbon affinity; once nucleated, growth is steady. Nucleation is a supersaturation threshold, not an Arrhenius hop. A single YAML step is a **lump** of transport + precipitation, and must be labeled that way.

### Branching (this is what “coking resistance” is)

- **C supply vs C removal.** Alves [1]: Ni activity rises to ~650 °C; above that, cracking outruns diffusion, carbon piles up, the particle encapsulates. Fe holds to ~800 °C. Faster kinetics often mean **less** total H₂ over catalyst life. Phase 2’s band is 500–1027 °C — the upper half is the Ni encapsulation regime.
- **Particle size.** Xu / Lopez-Ruiz / Kovarik / Dagle [29] (the structure-sensitivity result Gili 2024 [3] reviews): **>20 nm Ni → CNTs, higher TOF, longer life; <10 nm → graphitic layers, dies.** TOS death is fragmentation then encapsulation. Support identity is second-order.
- **Solubility / alloy / H₂.** Cu in Ni lowers solubility and can keep the surface cleaner until there is too much Cu; H₂ in the feed slows encapsulating carbon [1]. We run 5% Ar, no H₂.

`coking_index = ΔE_C − 2 ΔE_H` is a **slab binding descriptor**. It is not k_Cγ / k_Cδ. Do not attach `carbon_transfer_eV` to it without an explicit, declared map.

### The 1.5 eV Ni transport barrier and separate Amin (2011) encapsulation barrier

`CandidateKinetics.carbon_transfer_eV` defaults to 1.5 eV. That number is:

- ~33 kcal/mol = **1.43 eV** bulk C diffusion in Ni (Baker [24])
- **1.42 / 1.55 eV** terrace → graphene-edge transport on Ni (Abild-Pedersen [26])
- **147–149 kJ/mol = 1.52–1.55 eV** encapsulating-carbon formation (`carbon_encapsulation_eV = 1.53 eV`, Amin [32])

One Arrhenius barrier cannot be both “the step that frees the site” and “the step that kills the site.” `carbon_transfer_eV = 1.50 eV` is the Ni transport default (`template_default: Abild-Pedersen 2006 / Baker 1972`), paired with `carbon_encapsulation_eV = 1.53 eV` (Amin 2011).

### Class gate — isolated single-atom catalysts (`SAC`) cannot use this path

The Baker / Helveg cycle needs an extended metal particle (bulk to dissolve into, or a terrace and a step). Isolated atoms do not have that.

Akri et al. [33]: atomically dispersed Ni in DRM is coke-resistant **because it only breaks the first C–H** and cannot complete methane to C. A SAC that stays a SAC has no filament path.

`rh1_ceo2_akri2019_sac_ref` is an isolated single-atom catalyst ($\text{Rh}_1/\text{CeO}_2(111)$). ADR 0001 keeps SAC in scope as a high-T solid; that is not a grant of nanoparticle physics. Writing `C_s => C(gr) + site` into the **universal** YAML would give a single atom a path it does not have — the same class of error as raising Γ.

MetalFreeCarbon is a different kinetics (the carbon *is* the site) [Muradov, TURQUOISE_HYDROGEN Category C]. MoltenMetal already has the no-lattice path [5]. Exsolved perovskite nanoparticles would be Baker particles; the encoded ABO₃ is out of scope (ADR 0001).

## Acceptance & Implementation Checklist

- [x] **Step 1 (`baker1972_abildpedersen2006_c_gamma_lump`)** Cγ lump `C_s => C(gr) + site` with its own barrier (`carbon_transfer_eV`, default 1.5 eV, provenance Baker / Abild-Pedersen). Labeled as transport-to-edge, not “nucleation.”
- [x] **Step 2 (`amin2011_snoeck1997_c_delta_encapsulation`)** Competing Cδ channel (`C_s => C_encap_s`, 1.53 eV, second-order in $\theta_{\text{C}}$) that does **not** return the site.
- [x] **Step 3 (`akri2019_nanoparticle_class_gate`)** Class gate: nanoparticle metals (Ni, Fe, Co, and alloys in `SolidCatalyst`, `SupportedNanoparticle`, `HEA`, `SAA`) only. Not `SAC`/`DAC` (`rh1_ceo2_akri2019_sac_ref`). Not `MetalFreeCarbon`. Melt unchanged.
- [x] **Step 4 (`coking_index_unmapped_guard`)** Do not map either barrier onto `coking_index` unless the map is declared in the `.kinetics.json` sidecar.
- [x] **Step 5 (`bengaard2002_deutschmann2000_ni111_thermochemistry_closure`)** Unconstrained single-hop $\text{Ni(111)/SiO}_2$ baseline at **650–700 °C** (`sweeps/ni_np_b65_closure.yaml`), PFR, zero regen.
- [x] **Step 6 (`lander1952_snoeck1997_ni111_parametric_sweep`)** Star + 2-D sweep on `ni111_sio2_bengaard2002_ref` (`sweeps/ni_np_b66_*.yaml`): E_act 0.7–1.3, `A_γ` $10^6\text{–}10^{13}$ plus half-decades $10^8\text{–}10^{10}$, $\theta^*$ 0.2–0.8, $s_0$ $10^{-3}\text{–}10^{-1}$, and `A_γ` × $\theta^*$ at 923.15 / 973.15 K. Scored by `pipeline/reactors/sweeps/carbon_criteria.py`.
- [x] **Step 7 (`ermakova2002_takenaka2003_ni111_joint_filament_closure`)** Calibrated supported-$\text{Ni(111)/SiO}_2$ filament reference closure (`ni111_sio2_bengaard2002_ref`, `results/evidence/supported_ni_filament_closure.json`): coupling Lander's bulk-Ni carbon diffusion length $L = 17.5\ \text{nm}$ ($A_\gamma = D_0 / L^2 = 8.10\times 10^{11}\ \text{s}^{-1}$) with $\text{WHSV} = 60,000\ \text{mL}/(\text{g}_{\text{cat}}\cdot\text{h})$ ($\tau = 2.0\ \text{s}$ at 923–973 K) closes all five literature criteria simultaneously (`supported_ni_filament_closure_passed = True`).

## Cδ rate form (Snoeck 1997 / Amin 2011)

Both channels first-order with equal 10¹³ prefactors fixes the per-carbon encapsulation fraction at `f_δ = 1/(1 + exp(0.03 eV/kT)) ≈ 0.41`, independent of coverage, E_act, and effectively T. Cδ is therefore written **second order in θ_C** (Cantera `coverage-dependencies: {C_s: {m: 1.0}}`, mean-field island nucleation; Snoeck's supersaturation picture) with `A_δ = A_γ/θ*`. θ\* is the C_s coverage where encapsulation overtakes transport at equal barriers; with the 0.03 eV difference the effective crossover is `θ*·exp(0.03 eV/kT)` (0.73 at 923 K for θ\* = 0.5), recorded in the sidecar. Mechanical / consumable regen and circulating removal clear `C_s` only; `C_encap_s` is TOS death and is never cleared (oxidative burn-off is non-turquoise and not modelled).

## Writer corrections (`bengaard2002_deutschmann2000_ni111_thermochemistry_closure`)

1. **Reference states.** Cantera `constant-cp` surface species sit on the absolute scale where H₂ and graphite are zero. `H_s = dE_H` was right; `CH3_s` lacked `h_f(CH3•) = +145.7 kJ/mol` (every CH₃\* ~1.5 eV too stable); `C_s` lacked `h_f(CH4) = −74.6 kJ/mol` (the screener's C reference is CH₄ − 2H₂; every C\* ~0.8 eV too unstable). `CH2_s`/`CH_s` interpolate the ladder between the corrected `CH3_s` and `C_s` (template, declared).
2. **Bimolecular prefactors.** `A = 10¹³` under `units: {length: cm, quantity: mol}` for `X_s + site` and `2 H_s` steps is cm²/mol/s, an effective `A·Γ = 2.5×10⁴ s⁻¹`. Now `A = 10¹³/Γ = 4×10²¹ cm²/mol/s` (H₂ desorption 2×10²²), the Deutschmann / Cantera `methane_pox_on_pt` convention. Unimolecular Cγ/Cδ (1/s) were already right.
3. Consequence: the pre-correction `rh1_ceo2_akri2019_sac_ref` headline (2.90% at 1300 K) was **CH₃\*/H\* parking behind a frozen ladder**, θ_C ≈ 10⁻¹⁴, not a carbon monolayer. Corrected: PFR 0.65%, fluidized 0.29%, θ_C = 0.10, 0.22 of the site-inventory bound. Identities 1–3 above still hold for `rh1_ceo2_akri2019_sac_ref` (inventory corr(X, a) = 0.999); the parked species changed, the physics claim did not.

## Unconstrained Single-Hop Baseline (`sweeps/ni_np_b65_closure.yaml`)

Literature $\text{Ni(111)/SiO}_2$ cell (`ni111_sio2_bengaard2002_ref`): `E_act 1.00` (Bengaard [35] terrace TS is 101 kJ/mol = 1.05 eV, recorded 1.00), `dE_H −0.50` (Greeley TPD [37]), `dE_CH3 −1.95` and `dE_C +1.30` (Watwe [36], Abild-Pedersen [26]), unconstrained `A_γ = 10¹³`, θ\* = 0.5. Cells: production (0.13 mm / 0.5 / 0.3) and large-particle Ni (0.13 / 0.9 / 0.05). Zero regen.

| cell | T (K) | X | bound | turnovers/site | θ_C exit | θ_encap exit | Cγ/Cδ | gC/(gNi·h) |
|---|---|---|---|---|---|---|---|---|
| production | 923 | 14.6% | 2.07% | 7.0 | 1.8e-5 | 3.0e-4 | 2.3e4 | 389 |
| production | 973 | 26.3% | 2.18% | 12.1 | 1.1e-5 | 3.5e-4 | 3.4e4 | 667 |
| production | 1300 | 99.97% (> X_eq 98.5%, flagged) | 2.92% | 34.3 | 1e-9 | 1.7e-4 | 2e5 | 1895 |
| large-particle | 923 | 6.2% | 0.62% | 10.0 | 2.9e-5 | 5.7e-4 | 1.7e4 | 92 |
| large-particle | 973 | 11.6% | 0.65% | 17.7 | 1.9e-5 | 7.1e-4 | 2.5e4 | 163 |

## Parametric Star + 2-D Sweep (`sweeps/ni_np_b66_*.yaml`)

Same literature $\text{Ni(111)/SiO}_2$ cell (`ni111_sio2_bengaard2002_ref`). Five specs, 2312 complete records (`pipeline/reactors/sweeps/carbon_criteria.py` → `results/sweeps/b66_summary.json`):

| criterion | value | gate | result |
|---|---|---|---|
| Barrier sensitivity | \|d ln X / d E_act\| = **8.27 eV⁻¹** at production PFR 923 K regen 0 (X(0.9) = 31.21 %, X(1.0) = 14.56 %, X(1.1) = 5.97 %) | ≥ 5 eV⁻¹ | PASS |
| Dispersion non-linearity | turnovers 7.03 (production) vs 9.95 (large-particle), relative difference **34.3 %** | > 20 % | PASS |
| Encapsulation onset | first `A_γ` with exit θ_enc ≥ 0.1 is **10⁷ s⁻¹** (923 K: X 1.87 %, θ_enc 0.145, Cγ/Cδ 1.06) | report first T or `A_γ` | PASS |
| Filament TOS lifetime | **244** scorable rows with `encapsulation_lifetime_h` in 4–50 h | vs 4–50 h | PASS |
| Specific carbon yield | **54** scorable rows in 8–10 gC/(gNi·h) | vs Ermakova 8–10 | PASS |

## Calibrated Supported-$\text{Ni(111)/SiO}_2$ Filament Closure (`ni111_sio2_bengaard2002_ref`, `results/evidence/supported_ni_filament_closure.json`)

### Literature Base Point (`ni111_sio2_bengaard2002_ref`)

| field | value | provenance | cite |
|---|---|---|---|
| E_act | **1.00 eV** | Bengaard's Ni(111) terrace CH₄ TS is 101 kJ/mol = 1.05 eV (no ZPE). | [35] |
| dE_H | **−0.50 eV** | Experimental TPD on Ni(111) vs ½ H₂ (Greeley). | [37] |
| dE_CH3 | **−1.95 eV** | Screener convention (CH₃\* vs the CH₃ radical, Watwe). | [36] |
| dE_C | **+1.30 eV** | Terrace C vs CH₄−2 H₂ reference (Abild-Pedersen). | [26] |
| carbon_transfer_eV | **1.50 eV** | Baker's 33 kcal/mol = 1.43 eV (bulk C in Ni); terrace→graphene-edge 1.42 / 1.55 eV. | [24], [26] |
| carbon_encapsulation_eV | **1.53 eV** | Encapsulating-carbon formation 147–149 kJ/mol on Ni/Al₂O₃ (Amin). | [32] |
| L_diff / A_γ | **17.5 nm / 8.10×10¹¹ s⁻¹** | Lander bulk-Ni carbon diffusion length $L = 17.5\ \text{nm}$ ($A_\gamma = D_0 / L^2$, $D_0 = 2.48\times 10^{-4}\ \text{m}^2/\text{s}$). | [25], [27], [38] |
| WHSV | **60,000 mL/(g_cat·h)** | Literature differential/integral filament reactor space velocity ($\tau = 2.0\ \text{s}$ at 923–973 K). | [31], [34] |

`PipelineConfig.solids_judge_catalyst = "ni111_sio2_bengaard2002_ref"`. Headline band is **923.15–973.15 K** (650–700 °C), the Alves / Xu filament ROI. At this calibrated physical operating point (`evaluate_supported_ni_filament_closure` in `pipeline/reactors/scorecard.py`), `ni111_sio2_bengaard2002_ref` achieves simultaneous hits in both the Ermakova/Takenaka yield band (`8.84` and `9.62 gC/(gNi·h)` at 923 K and 973 K) and the 4–50 h filament TOS lifetime band (`14.8 h` and `8.6 h`) across PFR, PBR, and Fluidized reactors (`supported_ni_filament_closure_passed = True`).

## Refs (same numbering as ADR 0001)

1. Alves et al., *Renew. Sustain. Energy Rev.* **2021**, *137*, 110465.
3. Gili et al., *ChemCatChem* **2024**. DOI: [10.1002/cctc.202301629](https://doi.org/10.1002/cctc.202301629).
5. Upham et al., *Science* **2017**, *358*, 917–921.
24–38. See [ADR 0001](../adr/0001-pyrolysis-phase-admissibility.md) (Baker, Helveg, Abild-Pedersen, Snoeck, Xu, Gili 2019, Ermakova, Amin, Akri, Takenaka, Bengaard, Watwe, Greeley, Lander).
