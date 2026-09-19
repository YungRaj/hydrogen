# Changelog

## Unreleased

### Added

- Six explicit methane-conversion pathway modes, with thermocatalytic operation
  as the default.
- Portable OpenFOAM, FEniCSx, and Cantera solver discovery and execution.
- Strict external-result artifacts gated by identity, backend versions,
  convergence, mesh independence, conservation, physical bounds, provenance,
  and held-out model error.
- Mode-specific `hydrogen_case.json` validation for fluidized-bed, MMBCR, NTEC,
  and electrochemical calculations, including per-parameter sources and
  disjoint calibration/validation records.
- Reproducible `openfoam-env` and `fenicsx-env` definitions.
- Aqueous/molten electrochemical configuration and paired-control NTEC evidence
  boundaries.
- Checksum- and geometry-bound candidate Hamiltonian construction from standard
  FCIDUMP integrals using PySCF and OpenFermion.
- Pre-write validation for curated, checksum-bound prior-art manifests with an
  explicit blinded time split.
- A portable CUDA-Q/PySCF/OpenFermion environment definition.
- A checksum-bound experimental dataset contract for reactor, paired NTEC,
  MEA, durability, and hydrogen-impurity calibration/holdout evidence, including
  preregistered sample minima and experimental-unit leakage prevention.
- Strict OpenFOAM-to-FEniCSx field handoffs and Cantera coupling receipts,
  including mechanism/log/rate/history hashes and two-way residual convergence.
- Runner-owned, numbered OpenFOAM/FEniCSx outer iterations whose observed
  feedback states must match the final coupling receipt.
- Pristine-case enforcement and external FEniCSx script hashing for reproducible
  multiphysics input identities.
- An `analytical_hydrodynamic_closure` tier for Fluidized and MMBCR: when no
  validated OpenFOAM artifact or calibrated surrogate exists, the reduced
  Cantera model runs on a labelled bubbly-flow closure (Mendelson bubble rise
  / derived holdup; clipped `(u0-umf)/u0` bubble fraction) with
  `candidate_exclusion_authorized=False`. Artifact and surrogate closures
  take precedence when present. NTEC / electrochemical stay
  `validation_required`.

### Changed

- Packed-bed PFR geometry distinguishes total catalyst-bed area from void gas
  residence volume.
- Fluidized screening uses a reacting emulsion plus conservative bubble bypass,
  mixed on molar flows via the Ar tracer, populated by validated OpenFOAM
  hydrodynamics or the labelled analytical closure.
- MMBCR screening uses physical column flow, bubble diameter, and interfacial
  area. Residence time is `τ = H / u_b` with `u_b` from Mendelson; gas holdup
  is derived (`u_sup / u_b`) and fails closed above 0.3 (churn-turbulent).
  When a validated artifact supplies holdup, `τ = ε_g · H / u_sup`.
- Missing specialized physics remains visible and non-excluding instead of
  silently falling back to an unrelated reactor.
- Production and orchestrator CLIs propagate pathway and multiphysics artifact
  configuration end to end.
- Candidate slates (reactor and DFT) are drawn from the ADR 0001 admissibility
  pool in both the live discovery stage and the phase-2-only CSV restart. The
  DFT rescue route for invalid rows is preserved.
- Phase 2 headline is the solids scorecard (single-pass X of a named or best
  non-H-parked solids judge); MMBCR X_eq is reported separately and never ranks.

### Validation

- 74 pipeline tests, 42 scientific contracts, and 24 exclusion-audit checks
  pass locally (six contracts require QE pseudopotentials, `fairchem`, or
  OpenFOAM on the host).
- Real Cantera chemistry loading, OpenFOAM executable startup, and a FEniCSx
  finite-element solve were exercised successfully.
- External modes now include safe case-template generation, batch readiness,
  independently computed raw holdout error, quantitative mesh refinement,
  recomputed elemental/energy/charge balances, phase consistency, physical
  output identities, backend-specific logs, and an explicit NTEC handoff.

### Known scientific boundaries

- PFR is operational for qualified screening, not experimentally validated
  prediction; incomplete elementary kinetics remain explicitly labeled.
- Fluidized, MMBCR, NTEC, and electrochemical modes require real case inputs,
  candidate mechanisms, calibration data, and held-out validation before any
  result can exclude a candidate.
- Condensed `C(gr)` replaces the ideal-gas carbon tracer. Solids surface paths
  end at `C_s` unless the B6 nanoparticle gate admits Cγ / Cδ (Ni, Fe, Co).
  B5 (Phase 2 X before DFT) remains gated on B6-5/B6-6.

### Phase 2 log

Dated working-tree entries for the carbon-model / reactor-physics track.
Durable rule: [ADR 0001](docs/adr/0001-pyrolysis-phase-admissibility.md).
Open work: [`docs/backlog/`](docs/backlog/README.md) B1–B6.

#### 2026-09-19 — Sweep callers on pathway modes; emulsion area basis

`yaml_sweep`, `inventory_sweep`, and `eact_sensitivity` route each reactor
through its pathway mode (`SINGLE_REACTOR_MODE`), carry `material_class` /
`candidate_id` from the screening row (or `catalyst.material_class` +
`genome` for kinetics-only specs), and keep `not_applicable` /
`validation_required` / `failed` rows as records with a `status` column.
`phase2_scorecard` admits only `status == complete`. Fluidized emulsion area is
now per emulsion volume with gas volume `ε_mf = 0.45` (same basis as PFR;
the old pairing undercounted area by `1/ε_mf`). `inventory_sweep` builds the
probe / control mechanisms from `CandidateKinetics.from_screening_row`, so
the "cat_9 probe" is the same catalyst as the headline sweep (the template
path had default adsorbate enthalpies and gave ~100× lower X under the same
name). Headline `cat_9` 1300 K production cell: PFR 2.90 %, fluidized 1.19 %
(emulsion 2.38 %), MMBCR `not_applicable` (SAC). `run_eact_sweep` takes a
`kinetics=` record for B6-6.

#### 2026-09-19 — Merge onto upstream pathway-mode architecture

Upstream is the base. `reactor_models` keeps his `pathway_modes` routing,
`ReactorConfig` fields, `_validate_reactor_config`, NTEC / electrochemical
simulators, and artifact / surrogate coupling; our Γ lock, inventory, carbon
policies, fail-closed `[gas, graphite]` loading, melt ODE to X_eq, shared-surface
PFR, and emulsion fluidized bed are grafted on. Three physics adoptions: PFR
`A = sv · V_bed` (his area basis; corrects a 2.5× undercount), fluidized bubble
bypass on Ar-tracer molar flows, MMBCR `τ = H / u_b` (Mendelson). MMBCR on a
non-`MoltenMetal` class is `not_applicable`. `carbon_encapsulation_eV` counts
toward `quantitative_status` only when the B6 gate admits the candidate.
Headline MMBCR at 1300 K moves from X_eq-saturated (98.5 %) to Da-limited
(22.9 % in a 1.5 m column); column height is now the residence-time lever.

#### 2026-09-07 — B6-1–4: gated Cγ / Cδ on nanoparticle Ni/Fe/Co

`write_full_mechanism` emits `C_s => C(gr) + site` (1.5 eV, transport-to-edge) and `C_s => C_encap_s` (1.53 eV, encapsulating) only when the genome is SolidCatalyst / HEA / SAA-host Ni, Fe, or Co. SAC/DAC/`cat_9`/ungated names still end at `C_s`. Sidecar records the channels and `coking_index_mapped_to_off_site: false`. B5 judge is still `cat_9` at 1300 K (B6-5 open).

#### 2026-09-06 — Surface/graphite load fail-closed

`_load_gas_and_surface` no longer swallows a missing Langmuir surface or graphite phase. A `catalyst_name` that does not match the YAML (or a stale mechanism) raises unless `ReactorConfig.gas_only=True`. Results record `surface_loaded` / `graphite_loaded`; `is_solids_run` requires `surface_loaded is True`, so a blank X cannot enter the scorecard. Remaining Phase 2 cleanups (PFR nonlocal closure, unused imports, hand-rolled test harness, `_ch4_extent` import cycle) are listed in the README, not done.

#### 2026-09-06 — B6 documented: no intra-pass turnovers on solids

Surface YAML still ends at `C_s`. Not implemented. README, ADR 0001 refs [24]–[34], and [`docs/backlog/B6-off-site-carbon-nucleation.md`](docs/backlog/B6-off-site-carbon-nucleation.md) now state the three identities (E_act sweep is structurally flat; X ∝ a; melt vs bed is turnovers vs no turnovers) and the literature plan: Cα → Cγ lump plus a Cδ competitor, nanoparticle metals only, B5 judge off `cat_9` at 1300 K. B2 will not close B5. `carbon_transfer_eV = 1.5` remains unused (Baker / Abild-Pedersen Ni transport).

#### 2026-09-06 — Solids X mole balance, wired melt flotation

PFR/fluidized `CH4_conversion` is the Ar-tracer ratio `1 - (x_CH4/x_Ar)/(x_CH4,0/x_Ar,0)`. The H2 mole-balance is exact only for CH4/H2/inert (equilibrium check). The solids YAML has an active C2 chain, so that denominator deflates bed X in the low-Da regime. `1 - x_CH4/x_CH4,0` remains `CH4_mole_fraction_drop` only.

PFR/fluidized disable the Cantera energy equation (`thermal_mode=isothermal_energy_disabled`) so they match isothermal MMBCR instead of running adiabatic and self-quenching. Exit T is reported. Melt wall-loss / interface sag is still a separate duty account, not "no energy balance."

`mmbcr_carbon_removal_rate_1_s` is a flotation frequency. Production default is `None` = unconstrained (`η = 1`); the knob is wired but off. `0` fouls the interface. Detachment ablation labels PFR as not wired (discrete regen only). Circulating fluidized removal runs during integrate substeps (B2 slice).

Melt `k0 ≤ 1` m/s is a **guardrail** against inventing Da, not a physical bound like loading ≤ 1. Scorecard `judge_catalyst` is set only when a named judge is present; otherwise `headline_catalyst` is the best non-H-parked row. PFR `theta_C_axial` is labeled time-on-stream. MMBCR no longer reports fake C2 / `solid_C_selectivity=1`.

#### 2026-09-06 — Phase 2 carbon model, phase admissibility, MMBCR rate form

**Admissibility.** `phase_stable_at_application_T` rejects encoded phases that do not exist at pyrolysis T. First entry MetalHydride; also MOF, COF, MXene, Perovskite. Filter is **candidate selection only** — coverage certificate / 21.1B denominator / 14-class enumeration unchanged. MetalHydride is exempt from reserved validation quotas on both applications. MOF/COF/MXene/Perovskite are exempt from reserved **pyrolysis** validation slots (they still earn PEMFC quota).

**Carbon / mechanisms.** Gas-phase `C_graphite` removed. Condensed `C(gr)` (`fixed-stoichiometry`). Surface path ends at `C_s`. `co2_permitted=False`. Mechanical PFR outfeed default `max_regen_cycles=3`. H₂ selectivity gate off by default.

**Coking / surrogate.** MoltenMetal `coking_index` is explicit NaN. Missing key fails closed. Surrogate coking head uses masked MSE (no `fillna(0)`).

**MMBCR.** Rate is `k(E_act,T) × a_bubble × (X_eq − X)` with flotation (no Langmuir lattice). This **cannot exceed X_eq** and **reaches X_eq for large k·a·τ by construction**. The 0.1 eV → X = 0.985 result is a sanity check on k0 and area, not kinetic closure. **Closure: sanity-checked on MMBCR; pending on PFR and fluidized** (B5). Cross-reactor comparison now systematically favors MMBCR; within MMBCR ranking saturates once Da is large.

**Coking descriptor vs hydrides.** MetalHydride remains in `SLAB_COKING_INDEX_CLASSES` (slab exists). That is not pyrolysis admissibility.
