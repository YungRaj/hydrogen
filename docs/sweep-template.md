# Reactor-cell sweep YAML

Define a sweep in YAML, then run:

```powershell
conda activate cp2k-env
$env:PYTHONUTF8="1"
python runsweep.py sweeps/headline_cat9_1300K.yaml
```

`runsweep.py` loads the file, writes one `CandidateKinetics` mechanism, and
runs every cell at every listed temperature and reactor. That product is
the whole job: **cells × temperatures × reactors**. Results go to
`results/sweeps/<name>/run.json` (plus a copy of the input YAML). Specs under
`sweeps/` are git-tracked; run products under `results/sweeps/` are not.

This file is **not** a Cantera mechanism. Mechanism YAML lives under
`mechanisms/` (`phases:`, `species:`, `reactions:`). A sweep document uses
`name`, `catalyst`, `conditions`, and `cells`.

The stored example is [`sweeps/headline_cat9_1300K.yaml`](../sweeps/headline_cat9_1300K.yaml)
— `cat_9` at 1300 K, production fractional cell vs 1×1 envelope.

There is no `test` key and no pass/fail threshold in the sweep file.

## Reactor archetypes

`conditions.reactors` is a list of these three names. They are different carbon-handling pictures, not three views of the same bed. You can run any non-empty subset.

### PFR — packed solid bed

A packed column of catalyst particles. Feed flows through a 0.5 m bed at 0.05 m/s. The model is **one shared surface** marched through a chain of CSTRs: that is **time-on-stream**, not a true axial profile of fresh particles.

Carbon stays as site-blocking `C_s` until a named regen policy clears it. If coverage hits the threshold, the bed can do up to `max_regen_cycles` discrete mechanical clears and run again.

Cell inventory (`d_p`, loading, dispersion) sets active area `a = a_geom × loading × dispersion`. The table reports *a*, WHSV, and Ergun ΔP. This is the solids judge path.

### Fluidized — fluidized solid bed

The same solid-surface physics (`C_s` blocks sites; same inventory levers; Ergun ΔP), but the particles are fluidized (0.8 m bed, emulsion residence from `u0` / `u_mf`).

Two modes, chosen by `fluidized_mode`:

- **`circulating`** (default) — carbon is removed continuously during the integrate (a circulating-bed stand-in). No discrete regen loop.
- **`batch_regen`** — the bed parks and cokes; then the same discrete regen loop as the PFR (mechanical clear, up to `max_regen_cycles`).

`fluidized_mode` does nothing on PFR or MMBCR.

### MMBCR — molten-metal bubble column

A 1.5 m melt column. Gas rises as bubbles. There is **no site lattice** and no `C_s` inventory. Carbon leaves by flotation / transport, not by clearing a packed-bed surface.

The ODE is `dX/dt = k(E_act, T) · a_bubble · (X_eq − X)`, with `a_bubble = 6 / d_bubble` and residence `τ = H / u_b` (Mendelson bubble rise; ~7.3 s for the default column). For a large `k · a · τ` the conversion **is** `X_eq` by construction; the default column is Da-limited (~23% at 1300 K on a 0.43 eV melt). That number is not a catalyst rank and does not use the cell’s loading or dispersion. The table leaves *a*, WHSV, and ΔP blank.

**MMBCR only runs for `MoltenMetal` catalysts.** For any solids class (SAC, SolidCatalyst, HEA, …) the cell is reported as `not_applicable` with the reason, not dropped and not zero. `cat_9` is a SAC, so its MMBCR cells are `not_applicable`.

Use MMBCR as the melt contrast, not as a third solids score.

### Pathway modes and hydrodynamic closure

Upstream routes each reactor through a pathway mode (`thermocatalytic_pfr`, `thermocatalytic_fluidized`, `mmbcr`) and refuses a run that mixes solids and melt modes. A sweep file may still list all three reactors; the runner splits them into one call per mode.

Fluidized and MMBCR hydrodynamics are external-solver quantities upstream. When a validated OpenFOAM artifact or calibrated surrogate is present it is used. Without one the run proceeds on the labelled `analytical_hydrodynamic_closure` (Mendelson `u_b` / derived holdup for the melt; clipped `(u0−umf)/u0` bubble bypass for the bed) and records `closure_source`. Such a run can never exclude a candidate.

## What to set vs leave alone

**Set these.** They are why the file exists.

| Field | Why you set it |
|---|---|
| `name` | Output folder name under `results/sweeps/`. |
| `catalyst.name` | Cantera surface suffix (`mechanism_<name>.yaml`). |
| `catalyst.screening` **or** `catalyst.kinetics` | Where the barriers / adsorption energies come from. |
| `catalyst.material_class` (+ `genome`) | Required with `kinetics`. Reactor applicability (MMBCR = `MoltenMetal` only) and the B6 nanoparticle gate depend on it. Taken from the CSV row with `screening`. |
| `conditions.temperatures_K` | Which T points to run. |
| `conditions.reactors` | Which of the three reactor archetypes to run (see above). |
| `cells[]` inventory | Particle size, loading, and dispersion for that cell. |

**Copy the template defaults and leave them.** They are the turquoise production policy. Changing them is a different claim, not a more complete sweep.

| Field | Default | Why it is the default |
|---|---|---|
| `co2_permitted` | `false` | Turquoise pyrolysis does not burn carbon to CO₂. |
| `regen_mechanism` | `mechanical` | Solids: produce → mechanical outfeed/clear → return. No CO₂ chemistry. |
| `fluidized_mode` | `circulating` | Continuous carbon removal, not a parked batch-regen bed. |
| `max_regen_cycles` | `3` | Cap on discrete PFR / batch-fluidized regen loops. |

Omitting `policy` entirely uses those same four defaults. You do not need the block unless you are changing one of them on purpose.

**Do not put these in the sweep file.** They are locked in the reactor / mechanism code. The parser has no keys for them.

| Locked quantity | Value / rule |
|---|---|
| Site density Γ | `2.5e-9 mol/cm²` (monolayer). Extra sites come from particle S/V, loading, and dispersion only. |
| Loading / dispersion ceiling | Each ≤ 1. Do not invent BET. |
| Pressure, feed | 1 bar, `CH4:0.95, Ar:0.05`. |
| Bed geometry / velocity | `ReactorConfig` defaults (0.5 m PFR bed, 0.05 m/s, …). |
| MMBCR interfacial k₀ / flotation | Melt-side defaults. Not a DFT barrier. |
| `C_s => C(gr) + site` | Not a sweep key. Writer emits Cγ/Cδ only for nanoparticle Ni/Fe/Co. SAC/`cat_9` YAML still ends at `C_s`. |

## Options and ranges

### Catalyst

Use **one** source. Screening wins for the row; explicit `kinetics` is for a file that has no CSV.

| Option | Allowed values | Notes |
|---|---|---|
| `screening.csv` | Path relative to the repo root | Typical: `results/screening/ga_full_database.csv`. |
| `screening.index` | Integer ≥ 0 | 0-based pandas index after `read_csv`. `cat_9` is `index: 9`. |
| `kinetics.E_act` | Finite **eV > 0** | Methane activation barrier. Required if there is no screening row. |
| `kinetics.dE_H`, `.dE_CH3`, `.dE_C` | Finite eV, optional | Adsorption energies → surface enthalpies. **Not** barriers. Do not treat `\|dE_H\|` as H₂ desorption. |
| `material_class` | One of the 14 classes, or `MoltenMetal` | **Required with `kinetics`**; ignored with `screening` (the row wins). Decides which reactors apply and whether the B6 Cγ/Cδ channels are written. |
| `genome` | Genome tuple as a string, optional | e.g. `"('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)"`. Only the B6 gate reads it (metal must be Ni/Fe/Co on an extended particle). Omit for SAC / melts. |

### Conditions

| Option | Allowed values | Notes |
|---|---|---|
| `temperatures_K` | Number or list of K | Application band is **773.15–1300 K** (ADR 0001). The usual grid is `[773.15, 900, 1100, 1300]`. The parser does not reject numbers outside the band; those runs are off-contract. |
| `reactors` | `PFR`, `Fluidized`, `MMBCR` | Any non-empty subset. Unknown names fail closed. See **Reactor archetypes**. |

### Policy

| Option | Allowed values | Set this? |
|---|---|---|
| `co2_permitted` | `true` / `false` | Leave `false`. `true` is only for oxidative-regen tests. |
| `regen_mechanism` | `mechanical`, `consumable`, `oxidative` | Leave `mechanical`. `consumable` also clears `C_s` without CO₂. `oxidative` **requires** `co2_permitted: true` and is not a turquoise claim. |
| `fluidized_mode` | `circulating`, `batch_regen` | Leave `circulating`. `batch_regen` is the parked-bed contrast. Ignored by PFR and MMBCR. |
| `max_regen_cycles` | Integer. `0` disables discrete regen | Leave `3` unless you are probing the regen cap. |

### Cells

Each entry in `cells` is one inventory point. `a = a_geom(d_p) × loading × dispersion`.

| Field | Range | Production default | Why |
|---|---|---|---|
| `catalyst_particle_mm` | `> 0` | `0.13` | In-band Ergun (~0.40 bar) on this 0.5 m / 0.05 m/s bed. ~15× geometric *a* vs the old 2 mm pellet. Last Ergun-legal envelope cell with margin is ~0.10 mm (0.67 bar); **~0.08 mm is the 1 bar wall**. Cells above 1 bar still run and are marked. |
| `metal_loading` | `(0, 1]` | `0.5` | Area fraction of geometric pellet that is metal. **Not wt%.** 0.5 × 0.3 is the supported-TCD proxy (Alves 2021 / Sánchez-Bastardo 2021), not a bulk-metal pellet. |
| `metal_dispersion` | `(0, 1]` | `0.3` | Fraction of that metal that is surface-available (Gili 2024: accessible metal dies to encapsulation). Product *a* = **0.15 × a_geom**. |
| `name` | Non-empty string | — | Label in the results table. |

A 1×1 envelope cell (`loading: 1`, `dispersion: 1`) is the geometric upper bound at the same `d_p`. Envelope / fractional = `1 / (0.5 × 0.3)` = **6.67×** by identity. That ratio is not a catalyst finding.

B1 coarse archive levels (if you are repeating that grid, not inventing new ones): particles `2.0, 0.5, 0.2, 0.1` mm; loadings `1.0, 0.5, 0.2`; dispersions `1.0, 0.3, 0.1`. ROI refine: particles `0.25 … 0.08` mm; loadings `1.0, 0.7, 0.5`; dispersions `1.0, 0.5, 0.3`.

## Template

```yaml
name: my_sweep
description: Optional note. Printed with the results table.

catalyst:
  name: cat_9
  screening:
    csv: results/screening/ga_full_database.csv
    index: 9
  # Or, if there is no CSV (material_class is then required):
  # material_class: SolidCatalyst
  # genome: "('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)"
  # kinetics:
  #   E_act: 0.43
  #   dE_H: -0.90
  #   dE_CH3: -2.80
  #   dE_C: 2.52

conditions:
  temperatures_K: [773.15, 1300]
  reactors: [PFR, Fluidized, MMBCR]

policy:
  co2_permitted: false
  fluidized_mode: circulating
  max_regen_cycles: 3
  regen_mechanism: mechanical

cells:
  - name: fractional
    catalyst_particle_mm: 0.13
    metal_loading: 0.5
    metal_dispersion: 0.3
  - name: envelope_1x1
    catalyst_particle_mm: 0.13
    metal_loading: 1.0
    metal_dispersion: 1.0
```

## Output columns

| Field | Meaning |
|---|---|
| `cell` | `cells[].name` |
| `X` | Ar-tracer CH₄ conversion. Fluidized: after bubble bypass, `(1−δ)·X_emulsion`; `emulsion_CH4_conversion` is in the JSON. Blank when the cell did not run. |
| `a_1/m` | Active solids area; blank for MMBCR |
| `WHSV` | 1/τ in h⁻¹ (historical name) |
| `dP_bar` | Ergun ΔP; blank for MMBCR |
| `status` | `complete`, `not_applicable (reason)`, `validation_required`, or `failed (error)`. Only `complete` rows are catalyst evidence. |

The JSON also records `pathway_mode`, `material_class`, `closure_source`, `residence_time_s`, and `can_exclude_candidate` per row.
