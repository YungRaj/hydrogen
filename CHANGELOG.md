# Changelog

## [Unreleased]

### 2026-09-07 — B6-1–4: gated Cγ / Cδ on nanoparticle Ni/Fe/Co

`write_full_mechanism` emits `C_s => C(gr) + site` (1.5 eV, transport-to-edge) and `C_s => C_encap_s` (1.53 eV, encapsulating) only when the genome is SolidCatalyst / HEA / SAA-host Ni, Fe, or Co. SAC/DAC/`cat_9`/ungated names still end at `C_s`. Sidecar records the channels and `coking_index_mapped_to_off_site: false`. B5 judge is still `cat_9` at 1300 K (B6-5 open).

### 2026-09-06 — Surface/graphite load fail-closed

`_load_gas_and_surface` no longer swallows a missing Langmuir surface or graphite phase. A `catalyst_name` that does not match the YAML (or a stale mechanism) raises unless `ReactorConfig.gas_only=True`. Results record `surface_loaded` / `graphite_loaded`; `is_solids_run` requires `surface_loaded is True`, so a blank X cannot enter the scorecard. Remaining Phase 2 cleanups (PFR nonlocal closure, unused imports, hand-rolled test harness, `_ch4_extent` import cycle) are listed in the README, not done.

### 2026-09-06 — B6 documented: no intra-pass turnovers on solids

Surface YAML still ends at `C_s`. Not implemented. README, ADR 0001 refs [24]–[34], and [`docs/backlog/B6-off-site-carbon-nucleation.md`](docs/backlog/B6-off-site-carbon-nucleation.md) now state the three identities (E_act sweep is structurally flat; X ∝ a; melt vs bed is turnovers vs no turnovers) and the literature plan: Cα → Cγ lump plus a Cδ competitor, nanoparticle metals only, B5 judge off `cat_9` at 1300 K. B2 will not close B5. `carbon_transfer_eV = 1.5` remains unused (Baker / Abild-Pedersen Ni transport).

### 2026-09-06 — Solids X mole balance, wired melt flotation

PFR/fluidized `CH4_conversion` is the Ar-tracer ratio `1 - (x_CH4/x_Ar)/(x_CH4,0/x_Ar,0)`. The H2 mole-balance is exact only for CH4/H2/inert (equilibrium check). The solids YAML has an active C2 chain, so that denominator deflates bed X in the low-Da regime. `1 - x_CH4/x_CH4,0` remains `CH4_mole_fraction_drop` only.

PFR/fluidized disable the Cantera energy equation (`thermal_mode=isothermal_energy_disabled`) so they match isothermal MMBCR instead of running adiabatic and self-quenching. Exit T is reported. Melt wall-loss / interface sag is still a separate duty account, not “no energy balance.”

`mmbcr_carbon_removal_rate_1_s` is a flotation frequency. Production default is `None` = unconstrained (`η = 1`); the knob is wired but off. `0` fouls the interface. Detachment ablation labels PFR as not wired (discrete regen only). Circulating fluidized removal runs during integrate substeps (B2 slice).

Melt `k0 ≤ 1` m/s is a **guardrail** against inventing Da, not a physical bound like loading ≤ 1. Scorecard `judge_catalyst` is set only when a named judge is present; otherwise `headline_catalyst` is the best non-H-parked row. PFR `theta_C_axial` is labeled time-on-stream. MMBCR no longer reports fake C2 / `solid_C_selectivity=1`.

### 2026-09-06 — Phase 2 carbon model, phase admissibility, MMBCR rate form

Shipped working-tree changes. Durable rule: [ADR 0001](docs/adr/0001-pyrolysis-phase-admissibility.md). Open work: issues B1–B5.

**Admissibility.** `phase_stable_at_application_T` rejects encoded phases that do not exist at pyrolysis T. First entry MetalHydride; also MOF, COF, MXene, Perovskite. Filter is **candidate selection only** — coverage certificate / 21.1B denominator / 14-class enumeration unchanged. MetalHydride is exempt from reserved validation quotas on both applications. MOF/COF/MXene/Perovskite are exempt from reserved **pyrolysis** validation slots (they still earn PEMFC quota).

**Carbon / mechanisms.** Gas-phase `C_graphite` removed. Condensed `C(gr)` (`fixed-stoichiometry`). Surface path ends at `C_s`. `co2_permitted=False`. Mechanical PFR outfeed default `max_regen_cycles=3`. H₂ selectivity gate off by default.

**Coking / surrogate.** MoltenMetal `coking_index` is explicit NaN. Missing key fails closed. Surrogate coking head uses masked MSE (no `fillna(0)`).

**MMBCR.** Rate is `k(E_act,T) × a_bubble × (X_eq − X)` with flotation (no Langmuir lattice). This **cannot exceed X_eq** and **reaches X_eq for large k·a·τ by construction**. The 0.1 eV → X = 0.985 result is a sanity check on k0 and area, not kinetic closure. **Closure: sanity-checked on MMBCR; pending on PFR and fluidized** (B5). Cross-reactor comparison now systematically favors MMBCR; within MMBCR ranking saturates once Da is large.

**Coking descriptor vs hydrides.** MetalHydride remains in `SLAB_COKING_INDEX_CLASSES` (slab exists). That is not pyrolysis admissibility.
