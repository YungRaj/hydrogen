# B6 — Off-site carbon nucleation (Cα → Cγ / Cδ)

**Type:** Implementation. B6-1–3 in the writer. Do not ship a universal `C_s => C(gr) + site`.

**Status:** Gated extra block in `write_full_mechanism`. Ungated YAML still ends at `C_s`. Nanoparticle Ni/Fe/Co emit Cγ (`C_s => C(gr) + site`, 1.5 eV, `A_γ` declared, default 10¹³ s⁻¹) and Cδ (`C_s => C_encap_s`, 1.53 eV, **∝ θ_C²**, `A_δ = A_γ/θ*`, θ\* = 0.5 declared). Not mapped to `coking_index`. **B6-5 run (PFR)**: see "B6-5 result" below. **B6-6 run (PFR + Fluidized)**: see "B6-6 result" below; all five acceptance gates pass. **B6-7 done:** judge is `ni_np_lit` at the 650–700 °C ROI; literature CSV stub; joint band declared empty. B5 stays closed.

**Depends on:** B1 (done; Γ locked). **Blocks:** B5. **Does not replace:** B2 (between-pass outfeed) or B4 (packed-bed ΔP / τ).

## Why this is the missing physics

The solids path consumes the Langmuir lattice stoichiometrically inside a pass. H₂ can leave (`2 H_s => H2 + 2 site`). Carbon cannot. Real Ni TCD does not die after one monolayer: C leaves the active face, travels through or across the particle, and nucleates graphite at a **different** place (rear face or a step edge). That is what lets Ni run for hours [1, 3, 24–31, 34].

Three identities of the **current** mechanism (no off-site step):

1. **E_act sweep cannot discriminate on solids.** Once the pass parks C on the available sites, `X ≈ n_sites / n_CH4,fed = (Γ · a · V) / n_CH4`. The barrier is not in the answer.
2. **Inventory gains are linear in `a` by construction.** B1’s X ∝ a (corr ≈ 1) is that identity, not kinetic closure. Opposite of the B5 criterion (X moves because the barrier moved). Right now: 100% loading, 0% barrier.
3. **Melt vs bed is turnovers vs no turnovers**, not continuous-C-removal vs coking. The melt has no site conservation [5]. B2 / decoke restores sites **between** passes. The deficit is **within** a pass. Adding B2 will not close B5.

## What the literature actually is (not one elementary hop)

Rostrup-Nielsen / Trimm carbon types: our `C_s` is **Cα** (atomic, site-blocking). Missing channels:

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

### The unused 1.5 eV is already the Ni transport number — and it is overloaded

`CandidateKinetics.carbon_transfer_eV` defaults to 1.5 eV. That number is:

- ~33 kcal/mol = **1.43 eV** bulk C diffusion in Ni (Baker [24])
- **1.42 / 1.55 eV** terrace → graphene-edge transport on Ni (Abild-Pedersen [26])
- **147–149 kJ/mol = 1.52–1.55 eV** encapsulating-carbon formation (Amin [32])

One Arrhenius barrier cannot be both “the step that frees the site” and “the step that kills the site.” If a Cγ lump is added, **1.5 eV is the honest Ni transport default**, provenance `template_default: Abild-Pedersen 2006 / Baker 1972`. It is not a nucleation barrier and not the coking index.

### Class gate — the judge cannot use this path

The Baker / Helveg cycle needs an extended metal particle (bulk to dissolve into, or a terrace and a step). Isolated atoms do not have that.

Akri et al. [33]: atomically dispersed Ni in DRM is coke-resistant **because it only breaks the first C–H** and cannot complete methane to C. A SAC that stays a SAC has no filament path.

`cat_9` is SAC Rh. ADR 0001 keeps SAC in scope as a high-T solid; that is not a grant of nanoparticle physics. Writing `C_s => C(gr) + site` into the **universal** YAML would let B5 pass on `cat_9` by giving a single atom a path it does not have — the same class of error as raising Γ.

MetalFreeCarbon is a different kinetics (the carbon *is* the site) [Muradov, TURQUOISE_HYDROGEN Category C]. MoltenMetal already has the no-lattice path [5]. Exsolved perovskite nanoparticles would be Baker particles; the encoded ABO₃ is out of scope (ADR 0001).

## Acceptance (when we implement)

- [x] **B6-1** Cγ lump `C_s => C(gr) + site` with its own barrier (`carbon_transfer_eV`, default 1.5 eV, provenance Baker / Abild-Pedersen). Labeled as transport-to-edge, not “nucleation.”
- [x] **B6-2** Competing Cδ channel that does **not** return the site. Without this, the coking index still does no work.
- [x] **B6-3** Class gate: nanoparticle metals (Ni, Fe, Co, and alloys / exsolved particles) only. Not SAC/DAC. Not MetalFreeCarbon. Melt unchanged.
- [x] **B6-4** Do not map either barrier onto `coking_index` unless the map is declared in the `.kinetics.json` sidecar.
- [x] **B6-5** Closure experiment run: a supported-Ni-like literature cell at **650–700 °C** (plus the full ADR band), PFR, zero regen. Result below. Judge stayed on `cat_9` here. **Closed in B6-6:** Fluidized Ni. **Closed in B6-7:** judge moved to `ni_np_lit`; literature CSV stub (`fairchem_evaluated=False`).
- [x] **B6-6** Star + 2-D sweep on the Ni cell (`sweeps/ni_np_b66_*.yaml`): E_act 0.7–1.3, `A_γ` 10⁶–10¹³ plus half-decades 10⁸–10¹⁰, θ\* 0.2–0.8, s0 10⁻³–10⁻¹, and `A_γ` × θ\* at 923.15 / 973.15 K. PFR and circulating Fluidized, both cells, regen 0 and 3. Scored by `pipeline/process/b66_criteria.py`. Result below. Particle size remains first-order in the literature [3, 29] and enters as dispersion and, optionally, as `A_γ = D₀/L²` (`carbon_transfer_particle_nm`; Lander D₀ = 2.48×10⁻⁴ m²/s). Judge stays on `cat_9`.
- [x] **B6-7** Collective remainder after B6-6, one ticket. (1) Named solids judge is `ni_np_lit` (literature base point below); headline band is the 650–700 °C ROI, not 1300 K. (2) Tracked literature CSV stub `sweeps/ni_np_lit_screening_row.csv` (`fairchem_evaluated=False`; no fake UMA energies). (3) Joint-band search over B6-6 plus `sweeps/ni_np_b67_joint.yaml` (4000-row A_γ × θ\* × s0 cube at 923.15 / 973.15 K): **no simultaneous hit in the filament ROI**. One edge hit exists at 773 K / large-particle / θ\*=0.8 (Y=9.3, τ=4.12 h). B5 stays closed.

## Cδ rate form (decided 2026-09-19)

Both channels first-order with equal 10¹³ prefactors fixes the per-carbon encapsulation fraction at `f_δ = 1/(1 + exp(0.03 eV/kT)) ≈ 0.41`, independent of coverage, E_act, and effectively T. That form (a) moves the bound by 1/f_δ ≈ 2.5× and leaves X stoichiometric, (b) contradicts Ni TOS data by 3–4 orders (Ermakova's 384 gC/gNi is ~4×10⁴ carbons per surface site, so f_δ ≲ 10⁻⁴), and (c) has no supply/removal switch, so E_act can never matter. Cδ is therefore written **second order in θ_C** (Cantera `coverage-dependencies: {C_s: {m: 1.0}}`, mean-field island nucleation; Snoeck's supersaturation picture) with `A_δ = A_γ/θ*`. θ\* is the C_s coverage where encapsulation overtakes transport at equal barriers; with the 0.03 eV difference the effective crossover is `θ*·exp(0.03 eV/kT)` (0.73 at 923 K for θ\* = 0.5), recorded in the sidecar. θ\* = 0.5 is **declared, not measured**, and is a B6-6 sweep variable. Mechanical / consumable regen and circulating removal clear `C_s` only; `C_encap_s` is TOS death and is never cleared (oxidative burn-off is non-turquoise and not modelled).

## Writer corrections found by B6-5 (apply to every candidate)

1. **Reference states.** Cantera `constant-cp` surface species sit on the absolute scale where H₂ and graphite are zero. `H_s = dE_H` was right; `CH3_s` lacked `h_f(CH3•) = +145.7 kJ/mol` (every CH₃\* ~1.5 eV too stable); `C_s` lacked `h_f(CH4) = −74.6 kJ/mol` (the screener's C reference is CH₄ − 2H₂; every C\* ~0.8 eV too unstable). `CH2_s`/`CH_s` had hard-coded −15/−10 kJ/mol; they now interpolate the ladder between the corrected `CH3_s` and `C_s` (template, declared).
2. **Bimolecular prefactors.** `A = 10¹³` under `units: {length: cm, quantity: mol}` for `X_s + site` and `2 H_s` steps is cm²/mol/s, an effective `A·Γ = 2.5×10⁴ s⁻¹`. Now `A = 10¹³/Γ = 4×10²¹ cm²/mol/s` (H₂ desorption 2×10²²), the Deutschmann / Cantera `methane_pox_on_pt` convention. Unimolecular Cγ/Cδ (1/s) were already right.
3. Consequence: the pre-correction `cat_9` headline (2.90% at 1300 K) was **CH₃\*/H\* parking behind a frozen ladder**, θ_C ≈ 10⁻¹⁴, not a carbon monolayer. Corrected: PFR 0.65%, fluidized 0.29%, θ_C = 0.10, 0.22 of the site-inventory bound. Identities 1–3 above still hold for `cat_9` (inventory corr(X, a) = 0.999); the parked species changed, the physics claim did not.

## B6-5 result (2026-09-19, PFR, `sweeps/ni_np_b65_closure.yaml`)

Literature Ni(111)/SiO₂ cell: `E_act 1.00` (Bengaard [35] terrace TS is 101 kJ/mol = 1.05 eV, recorded 1.00; BEP from the adsorption energies would give 0.68), `dE_H −0.50` (Greeley TPD [37]), `dE_CH3 −1.95` and `dE_C +1.30` (screener conventions; CH₃ window Watwe [36], C mapping Abild-Pedersen [26]), `A_γ = 10¹³`, θ\* = 0.5. Cells: production (0.13 mm / 0.5 / 0.3) and large-particle Ni (0.13 / 0.9 / 0.05). Zero regen. `ni_np_b65_production.yaml` (3 mechanical cycles) is identical because θ_C never reaches the regen threshold.

| cell | T (K) | X | bound | turnovers/site | θ_C exit | θ_encap exit | Cγ/Cδ | gC/(gNi·h) |
|---|---|---|---|---|---|---|---|---|
| production | 923 | 14.6% | 2.07% | 7.0 | 1.8e-5 | 3.0e-4 | 2.3e4 | 389 |
| production | 973 | 26.3% | 2.18% | 12.1 | 1.1e-5 | 3.5e-4 | 3.4e4 | 667 |
| production | 1300 | 99.97% (> X_eq 98.5%, flagged) | 2.92% | 34.3 | 1e-9 | 1.7e-4 | 2e5 | 1895 |
| large-particle | 923 | 6.2% | 0.62% | 10.0 | 2.9e-5 | 5.7e-4 | 1.7e4 | 92 |
| large-particle | 973 | 11.6% | 0.65% | 17.7 | 1.9e-5 | 7.1e-4 | 2.5e4 | 163 |

**Read.** Turnovers ≫ 1 with zero regen: the closure criterion (X above the monolayer bound from intra-pass Cγ) is met in-model, and X now depends on the adsorption kinetics (every ladder step runs at the CH₄ sticking TOF, 4.7 → 1.2 s⁻¹ over the pass as H\* builds to 0.35). The yield is 50× above the Ermakova/Takenaka band and the encapsulation regime is unreachable at any T: with a single-hop `A_γ = 10¹³`, `k_γ(923 K) = 6.5×10⁴ s⁻¹` against an arrival rate of ~5 s⁻¹, so θ_C never approaches θ\*. `A_γ` is a transport + precipitation lump: at 10⁹ s⁻¹ (`D₀/L²` for a ~10 nm particle) θ_C crosses θ\* mid-pass and Cδ takes over (923 K: X 6.1%, θ_encap 0.40, 137 gC/(gNi·h)); at 10⁷ the surface starves (X 1.9%, turnovers 0.9, 8.5 gC/(gNi·h), inside the band but dead within the pass). No single `(A_γ, θ*)` pair matches both the yield band and a multi-hour lifetime with the current arrival rate, which points at the CH₄ sticking prefactor (template 0.01) and E_act as the other half of B6-6. The 1300 K overshoot of X_eq is Cγ irreversibility into a graphite sink (`exceeds_equilibrium` flag, never clipped).

**Residuals at B6-5 close.** Arrival rate and E_act were the other half of B6-6. Judge, literature stub, and joint band are **B6-7** (done).

## B6-6 result (2026-09-19, PFR + circulating Fluidized)

Same literature Ni(111)/SiO₂ cell as B6-5. Base point when not swept: E_act 1.0, `A_γ` 10¹³, θ\* 0.5, s0 0.01. Five specs, 2312 complete records (eact 392, A_γ 560, θ\* 280, s0 280, 2-D 800). 208 rows flag `exceeds_equilibrium` (almost all 1300 K with fast Cγ) and are not scored as successes. `pipeline/process/b66_criteria.py` → `results/sweeps/b66_summary.json`.

| criterion | value | gate | result |
|---|---|---|---|
| Flat | \|d ln X / d E_act\| = **8.27 eV⁻¹** at production PFR 923 K regen 0 (X(0.9) = 31.21 %, X(1.0) = 14.56 %, X(1.1) = 5.97 %) | ≥ 5 eV⁻¹ | PASS |
| Linearity | turnovers 7.03 (production) vs 9.95 (large-particle), relative difference **34.3 %** | > 20 % | PASS |
| Encapsulation onset | no T onset at `A_γ` = 10¹³; first `A_γ` with exit θ_enc ≥ 0.1 is **10⁷ s⁻¹** (923 K: X 1.87 %, θ_enc 0.145, Cγ/Cδ 1.06). Secondary Cγ/Cδ < 10 first at **10⁶ s⁻¹** (ratio 0.39, θ_enc 0.019) | report first T or `A_γ` | PASS |
| Lifetime | **244** scorable rows with `encapsulation_lifetime_h` in 4–50 h (example: E_act 1.0, 773 K, 4.34 h; Fluidized base 923 K, 5.03 h) | vs 4–50 h | PASS |
| Yield | **54** scorable rows in 8–10 gC/(gNi·h) (example: large-particle PFR 773 K, 9.30; `A_γ` = 10⁷ at 923 K, 8.47) | vs Ermakova 8–10 | PASS |

**Read.** X is no longer flat in the barrier and no longer a-scaling (turnovers differ by a third between the two cells at identical kinetics). Encapsulation is a transport-limited regime: it never appears on the T axis at the single-hop `A_γ` = 10¹³, and switches on around 10⁷–10⁹ s⁻¹ (`D₀/L²` for ~500 nm–5 µm on the Lander D₀, or a slow precipitation lump on a 10 nm particle). At `A_γ` = 10⁹ and 923 K the B6-5 switch is reproduced (X 6.06 %, θ_enc 0.40, Cγ/Cδ 6.1, 137 gC/(gNi·h)) and the pass dies (lifetime 0.001 h). The 8–10 yield band and the 4–50 h TOS band are populated, but not by the same rows: in-band yield at `A_γ` = 10⁷ is already encapsulated (lifetime ≪ 1 h), and in-band lifetimes sit at high `A_γ` / low T / low s0 where Cδ is a trickle. Fluidized Ni at the base point is X = 16.94 % (emulsion turnovers 20) versus PFR 14.56 % (7 turnovers); circulating outfeed is subtracted as B2, not counted as Cγ. 1300 K X > X_eq is flagged, never clipped, and never a success.

## B6-7 result (2026-09-20)

Three items, one ticket. B5 stays closed.

### Literature base point (the named judge)

These are the standard Ni(111) numbers used across steam-reforming / TCD DFT, not an in-house fit. Acceptable as the default judge because each is either the paper value, a documented rounding/mapping of that paper, or a declared template that sits on those papers. Provenance is in `sweeps/ni_np_b65_closure.yaml` and `sweeps/ni_np_lit_screening_row.csv`.

| field | value | why this is standard | cite |
|---|---|---|---|
| E_act | **1.00 eV** | Bengaard's Ni(111) terrace CH₄ TS is 101 kJ/mol = 1.05 eV (no ZPE). 1.00 is the recorded round number, not a table excerpt. A BEP from dE_CH₃ would give ~0.68 eV; the TS is the one that belongs on the sticking coefficient. | [35] |
| dE_H | **−0.50 eV** | Experimental TPD on Ni(111) vs ½ H₂. Greeley's PW91 DFT is −0.61 eV; −0.50 is the measured value they quote. | [37] |
| dE_CH3 | **−1.95 eV** | Screener convention (CH₃\* vs the CH₃ radical). Ni(111) fcc PW91 is typically −1.84 eV; −1.95 is in that window, not a Bengaard 2002 table value. | [36] |
| dE_C | **+1.30 eV** | Screener mapping, not a printed Abild-Pedersen number. Terrace C vs graphene is ~+0.5 eV; +0.8 eV puts it on the screener's CH₄−2 H₂ reference. | [26] |
| carbon_transfer_eV | **1.50 eV** | Baker's 33 kcal/mol = 1.43 eV (bulk C in Ni); terrace→graphene-edge 1.42 / 1.55 eV. Template midpoint. | [24], [26] |
| carbon_encapsulation_eV | **1.53 eV** | Encapsulating-carbon formation 147–149 kJ/mol on porous / nonporous Ni/Al₂O₃. | [32] |
| A_γ | **10¹³ s⁻¹** | Single-hop TST label. This is the *fast-transport* limit (filament-winning, not encapsulating), not D₀/L² for a 10 nm particle. Standard use for the judge: the literature cell with transport not rate-limiting. Yield is then high vs Ermakova because arrival (s0) and removal (A_γ) are both at the fast end — that is the B6-7 joint-band result, not a reason to discard the TS/adsorption set. | TST; D₀ alternative [38] |
| s0 | **0.01** | Deutschmann-style CH₄ sticking (`methane_pox_on_pt`, a Pt template). Ni molecular-beam values span ~10⁻⁴–10⁻² at these T; 0.01 is the upper-typical thermal default, not a Ni paper value. | Template; swept in B6-6. |
| θ\* | **0.5** | Declared crossover coverage, not measured. | B6-5/B6-6. |

`PipelineConfig.solids_judge_catalyst = ni_np_lit`. Headline band is **923.15–973.15 K** (650–700 °C), the Alves / Xu filament ROI. **1300 K is the right ADR ceiling** (Phase 2 is 773–1300 K, ADR 0001) and the right T for the X>X_eq flag. It is the *wrong* Ni headline: the legacy 4-point sweep used `headline_t_min=1200` so only 1300 counted, which is how `cat_9` became the 1300 K judge. Alves: Ni filaments lose to encapsulation above ~650 °C. X>X_eq rows never enter the headline.

### Literature CSV stub

`sweeps/ni_np_lit_screening_row.csv`: one row, `valid=True`, `screening_protocol=sweep_yaml_literature`, **`fairchem_evaluated=False`**. `CandidateKinetics.from_screening_row` and ADR 0001 scope accept it. No invented UMA energies.

### Joint band

`pipeline/process/b67_joint_band.py` over the five B6-6 runs plus `sweeps/ni_np_b67_joint.yaml` (4000 complete records: A_γ × θ\* × s0 × {0,3} regen × both cells × PFR/Fluidized at 923.15 / 973.15 K).

**Declaration: no simultaneous hit in the 650–700 °C ROI.** The 4000-row cube is empty of joint hits. Across all T, one edge pair exists: `large_particle_ni` PFR at **773 K**, θ\*=0.8, A_γ=10¹³, s0=0.01, Y=9.31 gC/(gNi·h), τ=4.12 h (regen 0 and 3 are identical). That is the cold end of the ADR band, not the filament ROI, on the low-area cell with encapsulation suppressed. Closest ROI misses still sit on opposite sides of the trade-off (in-band yield already dead; in-band lifetime a Cδ trickle). No extra DOF was added to force an ROI overlap. B5 stays closed: 14.6 % / 16.9 % at 923 K is not X_eq.

## What this is not

- Not B2. Outfeed / circulating removal is coverage policy between or across particles, not intra-pass Cα → Cγ chemistry.
- Not a license to raise Γ or k0 to force Damköhler.
- Not “add the step to every YAML and re-rank Phase 2.”

## Refs (same numbering as ADR 0001)

1. Alves et al., *Renew. Sustain. Energy Rev.* **2021**, *137*, 110465.
3. Gili et al., *ChemCatChem* **2024**. DOI: [10.1002/cctc.202301629](https://doi.org/10.1002/cctc.202301629).
5. Upham et al., *Science* **2017**, *358*, 917–921.
24–38. See [ADR 0001](../adr/0001-pyrolysis-phase-admissibility.md) (Baker, Helveg, Abild-Pedersen, Snoeck, Xu, Gili 2019, Ermakova, Amin, Akri, Takenaka, Bengaard, Watwe, Greeley, Lander).
