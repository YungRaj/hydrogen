# B6 — Off-site carbon nucleation (Cα → Cγ / Cδ)

**Type:** Implementation. B6-1–3 in the writer. Do not ship a universal `C_s => C(gr) + site`.

**Status:** Gated extra block in `write_full_mechanism`. Ungated YAML still ends at `C_s`. Nanoparticle Ni/Fe/Co emit Cγ (`C_s => C(gr) + site`, 1.5 eV) and Cδ (`C_s => C_encap_s`, 1.53 eV). Not mapped to `coking_index`. B6-5 / B5 still open.

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
- [ ] **B6-5** B5 judge moves off `cat_9` at 1300 K. Closure experiment: a supported-Ni-like genome at **650–700 °C**, where filaments actually win [1, 29]. Alves says Ni at 1300 K encapsulates — a fast off-site step that drives SAC Rh to X_eq there would contradict [1].
- [ ] **B6-6** After B6-1–3: E_act sweep on that Ni-like solids cell is no longer flat, and X is no longer linear in `a` alone (function of `k(E_act)·k_Cγ·a·τ` and the Cγ/Cδ branch). Particle size remains first-order in the literature [3, 29] and is **not** in the genome except as dispersion; a single k for every genome will still rank inventory. Document that residual.

## What this is not

- Not B2. Outfeed / circulating removal is coverage policy between or across particles, not intra-pass Cα → Cγ chemistry.
- Not a license to raise Γ or k0 to force Damköhler.
- Not “add the step to every YAML and re-rank Phase 2.”

## Refs (same numbering as ADR 0001)

1. Alves et al., *Renew. Sustain. Energy Rev.* **2021**, *137*, 110465.
3. Gili et al., *ChemCatChem* **2024**. DOI: [10.1002/cctc.202301629](https://doi.org/10.1002/cctc.202301629).
5. Upham et al., *Science* **2017**, *358*, 917–921.
24–34. See [ADR 0001](../adr/0001-pyrolysis-phase-admissibility.md) (Baker, Helveg, Abild-Pedersen, Snoeck, Xu, Gili 2019, Ermakova, Amin, Akri, Takenaka).
