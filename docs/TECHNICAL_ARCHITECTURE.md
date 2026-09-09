# Technical Architecture and Scientific Workflow

This document is the implementation-level guide to the turquoise-hydrogen and
fuel-cell catalyst discovery repository. It explains what each stage does, why
it exists, which library performs each calculation, where that integration
lives in the source tree, what data crosses each stage boundary, and what the
result is scientifically allowed to claim.

The most important distinction is this:

> The software is a search, screening, prioritization, simulation, and evidence
> management system. Enumeration, surrogate ranking, machine-learned atomistic
> energies, template microkinetics, toy VQE, and unconverged calculations are
> not treated as proof of a state-of-the-art catalyst. A final claim requires
> converged candidate-specific calculations, prior-art comparison, and measured
> reactor or MEA evidence.

## 1. Mission and scope

The project searches a very large, heterogeneous catalyst design space for two
connected applications:

1. Methane pyrolysis for turquoise hydrogen, in either ordinary
   thermocatalytic, MMBCR, NTEC, or broader electrochemical operation.
2. Oxygen-reduction cathodes and PEM fuel-cell systems that consume the
   resulting hydrogen.

The pipeline is designed to find promising and unusual candidates without
silently collapsing the search to familiar examples such as Fe, Pt, or Fe-N4.
It therefore combines deterministic coverage, class-preserving allocation,
surrogate-guided prioritization, atomistic screening, higher-fidelity
validation, process simulation, and evidence gates.

The repository does **not** model the entire physical chain at one uniform
level of theory. It is a multi-fidelity workflow. Cheap calculations traverse
the large space; more expensive calculations resolve a much smaller set; and
physical measurements remain necessary for discovery and performance claims.

## 2. Architecture at a glance

```text
14 catalyst classes / 21.1B raw encoded configurations
                         |
                         v
 deterministic index + design-space audit
                         |
                         v
 persistent divide-and-conquer branch search
  |                      |                         |
  | hard chemistry       | surrogate priorities  | coverage certificate
  | constraints only     | (never hard pruning)  | + persistent SQLite state
  v                      v                         v
 class/region champions + calibration probes + validation allocation
                         |
              +----------+-----------+
              |                      |
              v                      v
 turquoise-H2 atomistic screen       ORR atomistic screen
 Meta eSen / fairchem / ASE           Meta eSen / fairchem / ASE + CHE
              |                      |
              v                      v
 evidence-aware stage selection       DFT-resolution slate
              |
       +------+------+---------------------+
       |             |                     |
       v             v                     v
 Cantera          Quantum ESPRESSO       CUDA-Q VQE
 reactor models   DFT / NEB / ORR        optional model-Hamiltonian check
       |             |                     |
       v             v                     v
 conversion and   converged energies,   variational diagnostics;
 selectivity       paths, barriers,      not catalyst validation unless the
 diagnostics       frequencies           Hamiltonian becomes candidate-specific
       |
       +-------------------+------------------+
                           v
                  PEMFC cell + stack models
                           |
                           v
       prior art + novelty benchmark + hashed evidence manifest
                           |
                           v
        report, experimental slate, and fail-closed readiness result
```

The arrows represent information flow, not automatic upgrading of evidence.
For example, a surrogate barrier can parameterize a diagnostic Cantera
mechanism, but the resulting conversion is not a measured conversion and does
not become one merely because it is farther downstream.

## 3. Main entry points

### `run_production_campaign.py`

`--results-dir` gives each campaign an independent state, evidence, log, scan,
and report root. Generated Cantera mechanisms default to that root's
`mechanisms/` directory (or an explicit `--mechanisms-dir`), preventing
candidate names such as `catalyst_0` from colliding across runs. Ranker
calibration also distinguishes requested probes from valid converged rows: a
bounded deterministic refill evaluates new tree probes until the minimum
training contract is met or records explicit exhaustion.

This is the primary production-oriented campaign driver. Use it when the goal
is persistent, resumable traversal of the indexed population with explicit
resource controls and readiness checks.

It performs the following work:

- validates the raw and canonical design-space definitions;
- records run provenance, visible GPUs, command-line arguments, and Git state;
- runs deterministic branch discovery for turquoise-hydrogen candidates;
- sends quantitative candidates to diagnostic reactor simulation;
- creates a class-preserving DFT resolution slate;
- runs or queues DFT according to the chosen options;
- runs the corresponding fuel-cell branch discovery and ORR workflow;
- models PEMFC and stack performance for eligible candidates;
- writes reports and campaign readiness results;
- in `--final-campaign` mode, fails closed if required evidence is absent.

Important controls include `--hours`, `--branch-max-leaves`,
`--branch-leaf-size`, `--branch-probes`, `--scan-workers`,
`--validation-batch`, `--min-validation-per-class`, `--qe-mpi-ranks`,
`--qe-omp-threads`, and `--qe-max-concurrent`.

### `pipeline/orchestrator.py`

This is the compact six-phase end-to-end orchestrator. Its phases are:

1. deterministic branch-and-bound discovery;
2. Cantera reactor simulation;
3. Quantum ESPRESSO DFT validation;
4. CUDA-Q VQE diagnostics;
5. fuel-cell cathode screening, PEMFC modeling, and stack modeling;
6. report generation.

It supports phase-by-phase execution and a reduced `--quick` mode. Quick mode
is for integration testing and workflow inspection; it is not a completed
population search or a scientific validation campaign.

### `run_validation_campaign.py`

This is the safe inspection/resume interface for candidate-specific Quantum
ESPRESSO work. Without `--advance`, it reports which methane NEB or ORR stages
are missing, incomplete, failed, or converged. With `--advance`, it runs only
eligible next stages and preserves completed outputs. Nonempty incomplete ORR
outputs are left untouched unless `--restart-incomplete` is explicitly given,
because another process may still own them.

### `run_divide_conquer_pilot.py`

This constructs and evaluates the reproducible pilot benchmark used to compare
the divide-and-conquer method with random sampling and a deterministic expert
heuristic. The pilot uses legacy computational outcomes; its evidence level is
therefore explicitly lower than a blinded experimental benchmark. It measures
search-selection performance, not physical catalyst superiority.

### `audit_pipeline.py` and the test programs

`audit_pipeline.py` checks repository and scientific-policy contracts.
`test_pipeline.py`, `test_scientific_contracts.py`,
`test_gpu_affinity_contract.py`, and `test_vqe_solver_contract.py` cover general
behavior, equations/evidence rules, GPU assignment, and VQE solver invariants.

### `live_dashboard.py`

The dashboard reads logs, CSV counts, and NVIDIA telemetry to display campaign
progress. It is observational. It does not schedule work or alter evidence.

## 4. Source-tree ownership

| Directory or file | Responsibility |
|---|---|
| `pipeline/common/` | Paths, constants, executable discovery, catalyst definitions, provenance, applicability, confidence/OOD logic |
| `pipeline/search/` | Indexed population, deterministic branch traversal, persistent scans, coverage, adaptive validation, diverse discovery batches |
| `pipeline/screening/` | eSen/fairchem structure construction and energy evaluation, relaxation, surrogate models, GPU worker runtime, stage admission |
| `pipeline/process/` | Candidate mechanisms, Cantera reactor models, NTEC transfer model, PEMFC model, stack model, techno-economics |
| `pipeline/stages/` | Explicit boundaries that translate one phase's record into the next phase's inputs |
| `pipeline/validation/` | Quantum ESPRESSO input/execution/parsing, ORR corrections, production validation sequencing, task queue, viability, CUDA-Q VQE |
| `pipeline/evidence/` | Coverage/readiness, prior art, manifests, novelty and pilot benchmarks, status and reports |
| `quantum_espresso/` | Repository-owned SSSP manifest/data used to validate pseudopotentials and cutoffs |
| `mechanisms/` | Generated candidate-specific Cantera YAML and kinetics sidecars |
| `results/` | Generated screening databases, certificates, calculations, reports, logs, and run state |

The dependency direction is intentionally mostly one-way: common definitions
support search and screening; selected screening records feed process and
validation stages; evidence modules inspect artifacts rather than silently
changing scientific results.

## 5. The 21.1-billion-candidate design space

### 5.1 What the number means

`pipeline/common/catalyst_spaces.py::estimate_design_space_size` computes the
raw Cartesian population. The current declared total is **21,092,645,031 raw
encoded configurations** across 14 classes. Some encodings represent the same
canonical chemistry—for example, ordering a dopant pair in two ways—so the
design-space audit also reports **10,815,793,768 canonical identities**.

The raw total is the coverage denominator. Invalid or redundant Cartesian
entries are recorded with explicit rejection reasons; they do not disappear
from the denominator. This prevents a search from claiming exhaustive coverage
after silently skipping inconvenient parts of the space.

The 14 classes are:

1. `MoltenMetal`
2. `SolidCatalyst`
3. `SAC`
4. `DAC`
5. `MOF`
6. `COF`
7. `Perovskite`
8. `MetalHydride`
9. `MAXPhase`
10. `HEA`
11. `Spinel`
12. `MXene`
13. `SAA`
14. `MetalFreeCarbon`

The class genes include combinations of active metals, hosts, promoters,
supports, facets, strain, dopants, vacancies, coordination environments,
framework linkers, pore/cavity types, crystal families, surface terminations,
loadings, and temperatures. The exact axes live next to each class definition
in `pipeline/common/catalyst_spaces.py`; their deterministic inverse mapping
lives in `pipeline/search/indexed_space.py::candidate_at_class`.

### 5.2 O(1) indexed addressing

The software never materializes a Python list containing 21 billion tuples.
`candidate_at(global_index)` maps an integer in `[0, TOTAL_SIZE)` to a catalyst
in time proportional to the number of genes. Per-class offsets identify the
class, and mixed-radix decoding or combination unranking reconstructs the
genome.

This matters because deterministic integer addressing provides:

- reproducible candidates independent of process timing;
- disjoint worker shards through `iter_shard`;
- constant-size resume cursors;
- exact interval coverage certificates;
- no memory requirement proportional to 21 billion;
- stable candidate identities across similar campaign runs.

### 5.3 Raw admissibility versus model uncertainty

`is_physically_admissible` applies conservative, auditable constraints such as
invalid promoter fractions, duplicate ordered encodings, impossible dopant
settings, and application scope. Its policy version is stored as
`canonical-chemistry-v1`.

Model uncertainty is not a physical inadmissibility proof. Low eSen confidence,
surrogate disagreement, failed relaxation, absent pseudopotentials, or an
incomplete downstream calculation therefore must not be turned into a hard
chemical rejection.

## 6. Deterministic divide-and-conquer search

### 6.1 Why branch the space

A full high-fidelity calculation on every configuration is infeasible. The
branch search divides each class's integer interval into a binary tree. Probe
scores estimate which unresolved branches are promising or uncertain, while
leaf intervals are streamed through persistent scanners.

The implementation is in:

- `pipeline/search/branch_search.py` for the interval tree and scheduler;
- `pipeline/search/exhaustive_search.py` for leaf scanning and archives;
- `pipeline/screening/genetic_optimizer.py::run_branch_discovery` for the
  turquoise-hydrogen integration;
- `pipeline/screening/fc_genetic_optimizer.py::run_fc_branch_discovery` for the
  fuel-cell integration.

The filenames retain earlier “genetic optimizer” terminology for compatibility,
but the production entry points use branch discovery. Legacy GA entry points
are not the production search policy.

### 6.2 Branch lifecycle

Each branch node has an application, material class, half-open index interval,
depth, priority, probe statistics, status, and reason. A node is:

- **pending** while unresolved;
- **expanded** after it is divided into two child intervals;
- **scanned** after its terminal leaf cursor reaches the end;
- **pruned** only when every member of a sufficiently small interval is
  exhaustively proven to fail hard admissibility constraints.

Surrogate predictions can order pending nodes. They cannot mark a node pruned.
The source makes this invariant explicit: a hard prune is a proof over every
member, never a surrogate guess.

### 6.3 Coverage and exploration

The scheduler has three modes:

- `class_floor`: resolve a configured minimum number of leaves in every class;
- `balanced_exploration`: periodically choose the least-covered class;
- `priority`: exploit the best current branch priority.

This preserves coverage while spending additional work where expected value is
larger. Repeatedly unproductive regions receive lower priority, not zero
coverage. When a campaign resumes, a bounded number of stale pending priorities
can be refreshed using the latest calibration evidence.

### 6.4 Persistent scanning and archives

`run_sharded_scan` divides a leaf among independent CPU scanner processes.
SQLite stores scan cursors and bounded regional, global, and objective archives.
Worker databases are merged deterministically. This permits interruption and
resume without redoing the whole space.

These scans use cheap deterministic features or trained rankers. They are a way
to decide which candidates deserve expensive calculations; they are not a
substitute for those calculations.

### 6.5 Coverage certificates

`verify_branch_coverage` checks that terminal intervals form an exact partition
of every declared class. It detects gaps, overlaps, invalid bounds, unfinished
scan cursors, and prunes without a reproducible hard proof. It records a digest
and the unresolved encoded population.

A complete certificate proves address-space accounting under the recorded
policy. It does **not** prove that predictions match nature, that every geometry
was relaxed, or that every candidate was evaluated with DFT.

## 7. Feedback, calibration, and validation allocation

`pipeline/search/adaptive_validation.py` implements the coverage-feedback loop.
It stores paired predictions and higher-fidelity observations in SQLite with a
candidate ID, material class, chemistry region, fidelity, absolute error,
productivity flag, timestamp, and mandatory source provenance.

Calibration is tracked at two resolutions:

- exact discovery region when sufficient local evidence exists;
- material class as a fallback for unseen regions.

`allocate_validation_batch` first reserves `min_per_class` slots for every
represented class. Remaining slots use a deterministic combined score:

```text
0.30 quality
+ 0.20 expected improvement
+ 0.25 predictive uncertainty
+ 0.15 observed calibration error
+ 0.10 observed productivity
```

This realizes six important policies:

- every material class receives a fixed validation budget;
- extra calculations follow expected improvement and uncertainty;
- calibration error is tracked separately by chemistry region;
- disagreement with DFT or experiments increases validation priority;
- persistently unproductive regions lose priority but retain coverage;
- `experimental_slate` preserves one regional champion before selecting
  repeated candidates from the same region.

Candidate IDs and region definitions are centralized in
`pipeline/search/discovery.py`. Stable IDs make comparisons and feedback usable
across runs.

## 8. Atomistic screening with Meta eSen, fairchem, ASE, and PyTorch

### 8.1 What each library does

**ASE (Atomic Simulation Environment)** supplies atom containers, periodic
cells, slab/cluster construction, adsorbate placement, optimizers, calculator
interfaces, and Quantum ESPRESSO I/O. It is the common geometry/data interface
between ML screening and DFT validation.

**fairchem** supplies Meta's Open Catalyst model loading and ASE-compatible
inference calculator. The default model is
`esen-sm-conserving-all-oc25`.

**eSen** is the learned interatomic model used to estimate energies and forces
for clean and adsorbate-bearing candidate structures. It makes atomistic
screening of a validation batch much cheaper than running DFT for every member.
Its output remains ML screening evidence.

**PyTorch** executes the eSen/fairchem model on CUDA GPUs and also implements
the neural surrogate ensembles used in some ranking paths.

### 8.2 Where the integration lives

- `pipeline/screening/surface_calculator.py::get_ocp_calculator` loads the
  fairchem eSen model.
- `pipeline/screening/surface_screener.py::generate_structure` maps each
  catalyst genome to an ASE structure.
- `pipeline/screening/surface_screener.py::evaluate_candidate` runs the
  turquoise-hydrogen adsorbate protocol.
- `pipeline/screening/fc_screener.py::evaluate_orr_candidate` runs the ORR
  adsorbate protocol.
- `pipeline/screening/relaxation.py` owns relaxation records and fail-closed
  convergence handling.
- `pipeline/screening/batched_calculator.py` multiplexes requests through one
  loaded model per service.
- `pipeline/screening/gpu_executor.py` assigns worker processes to visible
  GPUs.
- `pipeline/screening/worker_supervisor.py` owns leases, heartbeats, requeues,
  health records, and abnormal-worker detection.

### 8.3 Structure representations

`generate_structure` has class-specific builders rather than mapping every
genome to a generic metal slab. Examples include alloy slabs, porphyrin-like
single/dual-atom environments, perovskite and spinel slabs, hydride slabs,
MAX/MXene representations, high-entropy alloys, and supported or carbon-based
sites.

These structures are screening representations. A selected candidate still
needs candidate-specific crystallography, defect/site enumeration, cell-size
checks, and converged DFT before a quantitative discovery claim.

### 8.4 Relaxation protocol

Both screening applications use immutable protocol identifiers from
`pipeline/screening/protocols.py`. Current IDs are:

- `esen-sm-conserving-all-oc25:relax-v3:pyrolysis-v2`
- `esen-sm-conserving-all-oc25:relax-v3:orr-che-v2`

Reference calculations use a target force of 0.05 eV/A with 200 steps. Clean
surfaces use 0.08 eV/A with 150 steps. Adsorbates use 0.08 eV/A with 100 steps.
`relax_with_record` attempts BFGS and can recover with FIRE. Records include
initial/final geometry digests, final maximum force, attempts, and failure
reasons.

An optimizer returning is not enough. The final force and geometry checks must
pass. A failed or censored relaxation becomes `validation_required`; it is not
evidence that the chemistry is poor.

### 8.5 Turquoise-hydrogen descriptors

The screening workflow evaluates a clean structure and adsorbed H, CH3, and C
states relative to reference calculations. It derives adsorption energies,
the methane reaction-energy split, a BEP-estimated activation barrier,
segregation/binding diagnostics where available, and a coking index.

The relevant physical helper functions are in `pipeline/common/utils.py`:

```text
k = A exp(-Ea / (kB T))                 Arrhenius rate, Ea in eV
k_TST = (kB T / h) exp(DeltaS / kB)     transition-state prefactor
Ea = max(E_floor, alpha DeltaE + beta)  bounded BEP estimate
```

BEP is a screening relationship, not a candidate-specific transition-state
calculation. A candidate's converged NEB barrier supersedes it.

### 8.6 ORR descriptors and the computational hydrogen electrode

The ORR screen evaluates clean, OH*, O*, and OOH* states plus H2/H2O reference
energies. The computational hydrogen electrode constructs free-energy steps for
the four-electron oxygen-reduction pathway and calculates the limiting step and
overpotential.

`pipeline/common/utils.py::orr_overpotential` performs the core screening
calculation. `pipeline/validation/orr_workflows.py` provides site enumeration,
explicit ORR corrections, and lowest-energy-site selection for higher-fidelity
work.

The screen also carries cost, elemental safety, applicability, Fenton-risk, and
confidence information. A finite ML overpotential is a ranking descriptor, not
a measured half-cell or MEA result.

### 8.7 Confidence and out-of-distribution behavior

`pipeline/common/ood_detector.py` combines elemental training coverage and,
when available, disagreement between independent calculators. Confidence can
change priority and request DFT. It does not grant permission to erase unfamiliar
chemistry. This is central to finding candidates outside familiar regions.

### 8.8 GPU execution and fault tolerance

`execution_layout` maps workers over all visible CUDA devices. Each worker is
given a numerical GPU ID and UUID, sets device affinity, and uses the shared
worker loop for either pyrolysis or ORR screening. Batched inference amortizes
model loading and improves GPU occupancy.

The supervisor leases tasks before execution, consumes heartbeats, requeues
work from dead or stale workers, records health, and writes results atomically.
Increasing workers per GPU can improve feeding and overlap until memory,
structure relaxation latency, or model throughput becomes limiting. The best
setting is machine-specific; the architecture does not assume a fixed GPU
count or model path.

## 9. Evidence-aware stage boundaries

`pipeline/screening/stage_selection.py` defines three dispositions:

- `quantitative_screening`: valid, finite, uncensored primary descriptor;
- `validation_required`: unresolved, failed, censored, uncertain, or explicitly
  marked for higher fidelity;
- `hard_excluded`: supported hard hazard such as toxic/radioactive chemistry.

`select_for_reactor` accepts only quantitative screening records because
Cantera requires finite kinetic inputs. `select_for_validation` admits both
quantitative and unresolved records, excludes only hard hazards, gives
unresolved candidates first access, and preserves class quotas.

This separation prevents a technical failure in a cheap stage from becoming a
parasitic chemical screen-out. It also prevents missing numerical inputs from
being silently invented for a quantitative reactor ranking.

## 10. Candidate-specific mechanism translation

### 10.1 Boundary module

`pipeline/stages/reactor.py::simulate_candidate` is the explicit boundary from
a selected screening row to Cantera. It:

1. builds `CandidateKinetics` from the row;
2. writes a candidate mechanism and provenance sidecar;
3. runs the requested temperature/reactor sweep;
4. rejects mock output when production forbids it;
5. returns the mechanism path, sweep, and best diagnostic condition.

### 10.2 What is candidate-specific today

`pipeline/process/reactor_mechanisms.py::CandidateKinetics` consumes:

- `E_act` as the candidate's methane activation estimate;
- `dE_H` as H adsorption thermochemistry;
- `dE_CH3` as CH3 adsorption thermochemistry;
- `dE_C` as surface-carbon adsorption thermochemistry;
- candidate ID and screening protocol provenance.

Adsorption energies modify surface thermochemistry. They are not silently
reinterpreted as elementary activation barriers.

### 10.3 What is still templated

When candidate-specific values are absent, CH3, CH2, and CH dehydrogenation,
H2 desorption, and carbon-transfer barriers use declared defaults. The sidecar
labels each value as screening-derived, candidate-specific, or
`template_default`. Any default makes the mechanism
`screening_template_incomplete`.

This makes Cantera more candidate-aware than a single universal mechanism, but
it is not a completed ab initio microkinetic model. Production NEB or measured
kinetics must replace template elementary barriers.

## 11. Cantera reactor simulation

### 11.1 Why Cantera is used

Quantum calculations estimate local energies and barriers. They do not by
themselves predict conversion through a finite reactor with residence time,
transport approximations, temperature, pressure, feed composition, and
surface/gas kinetics. **Cantera** provides thermodynamic phases, reaction-rate
evaluation, reactor-network integration, and surface chemistry for that scale.

The integration lives in:

- `pipeline/process/reactor_mechanisms.py` for Cantera YAML generation;
- `pipeline/process/reactor_models.py` for loading phases and simulating reactor
  configurations;
- `pipeline/stages/reactor.py` for the screening-to-reactor boundary.

### 11.2 Mechanism contents

Each generated YAML includes an ideal gas phase, an ideal surface phase,
surface intermediates, gas reactions, and candidate-parameterized surface
reactions. A `.kinetics.json` file beside the YAML records source provenance,
quantitative completeness, and carbon-phase status.

The surface sequence represents methane adsorption/dehydrogenation, hydrogen
desorption, and carbon transfer. Reverse reactions and thermochemistry are used
where the current mechanism declares them.

### 11.3 Reactor models

`ReactorConfig` supplies temperature, pressure, inlet composition, flow and
geometry parameters. `simulate_reactor` dispatches to:

- `simulate_mmbcr`: a molten-metal bubble-column representation using a CSTR
  cascade, bubble surface-to-volume ratio, and axial conversion profile;
- `simulate_pfr`: a staged plug-flow representation;
- `simulate_fluidized_bed`: a two-phase fluidized-bed proxy;
- `simulate_ntec_pathway`: consumes a validated OpenFOAM + FEniCSx + Cantera
  artifact for the liquid-solid NTEC pathway;
- `simulate_electrochemical_pathway`: consumes a validated FEniCSx + Cantera
  artifact whose configuration selects aqueous or molten electrolyte physics.

`run_reactor_sweep` runs combinations of temperature and reactor type and writes
JSON results under `results/reactor/`. Outputs include methane conversion,
hydrogen and carbon selectivity diagnostics, residence/axial information, and
the mechanism evidence metadata. Before dispatch,
`validate_mode_reactors` requires the exact mode-to-reactor mapping and
`reactor_applicability` checks the candidate material phase. Wrong pairings are
recorded as non-excluding `not_applicable` evidence and are never sent into
Cantera.

The thermal models have deliberately bounded interpretations:

| Reactor | Cantera construction | Geometry/phase contract | Limitation |
|---|---|---|---|
| Packed-bed PFR | Sequential isothermal ideal-gas reactors followed as a Lagrangian parcel, each with `ReactorSurface` | Solid catalyst; surface area from total packed volume, gas residence from void volume | Staged 1-D approximation; no axial dispersion, pressure drop, heat/mass-transfer limitation, or pellet diffusion |
| Fluidized bed | Isothermal reacting emulsion parcel plus explicit unreacted bubble bypass mixing | Particulate solid catalyst; superficial velocity must exceed minimum fluidization velocity | Conservative two-phase screening closure; no interphase exchange correlation, population balance, attrition, or CFD |
| MMBCR | Isothermal steady `IdealGasReactor` CSTRs in series with `ReactorSurface` | `MoltenMetal` candidates only; column area sets flow, gas holdup sets residence, bubble diameter and holdup set interfacial area; gas-liquid kinetics currently use a declared Cantera ideal-surface proxy | Idealized interface/bubbles/CSTRs; no liquid-phase activity model, coalescence, breakup, circulation, mass-transfer coefficient, or CFD |

Cantera integrates the declared gas and heterogeneous surface reaction network
inside those idealized thermal control volumes. Cantera also supports declared
electrochemical interface reactions, phase electric potentials, and associated
rate expressions. That capability is a kinetics component, not a complete
aqueous or molten electrochemical reactor: the methane-specific charge-transfer
mechanism, ionic/electronic transport closure, potential/current boundary
conditions, and reactor geometry must still be supplied and coupled. Cantera
does not supply NTEC contact electrification/mechanical coupling or detailed
multiphase hydrodynamics. Those equations belong to mode-specific
OpenFOAM/FEniCSx cases. The repository executes those cases, validates their
artifacts, and never falls back to a thermal bed; the scientific case inputs
must still be candidate-specific and calibrated.

See Cantera's official documentation for
[electrochemical interface reactions](https://cantera.org/stable/yaml/reactions.html#electrochemical)
and [reactor-network scope](https://cantera.org/stable/reference/reactors/).

The shared reactor stage generates a candidate thermal Cantera YAML only for
the three thermal dispatch targets. NTEC and electrochemical execution leaves
`mechanism_file` null and instead requires a mode-owned artifact. A valid
artifact completes the condition; a missing, mismatched, unconverged,
non-conservative, or non-mesh-independent artifact returns
`validation_required`. `failed`, `not_applicable`, and `validation_required`
conditions are counted separately.

External cases pass through three independently implemented gates:

1. `physical_case.py` verifies routing, identity, positive unit-bearing inputs,
   geometry relationships, non-placeholder parameter/model sources, feed and
   kinetics provenance, electrolyte phase, and disjoint calibration/holdout
   identifiers before solver launch. Its CLI creates deliberately non-runnable
   templates and validates completed cases.
2. `multiphysics_runner.py` rejects stale outputs, hashes pristine inputs plus
   any external FEniCSx script, writes backend-specific logs, captures solver
   versions, validates the OpenFOAM identity/unit/mesh/checksum-bound handoff,
   verifies checksummed Cantera mechanism/log/rate exchange and converged
   two-way coupling history, owns the numbered residual-based outer loop, and
   independently scores raw holdout records.
3. `multiphysics_contract.py` recomputes mesh change from at least three refined
   meshes, recomputes required mass/carbon/hydrogen/energy/charge residuals
   from inlet/outlet budgets, and checks reactor-specific identities such as
   fluidization velocity and electrochemical power. Self-attested Boolean
   success fields cannot substitute for these raw values.

`multiphysics_prepare.py` applies preflight across a candidate manifest and can
create missing non-runnable skeletons. It never invents physical values or
marks a template ready.

### 11.4 Current carbon-model limitation

The current mechanism still lists `C_graphite` inside the ideal-gas phase as a
legacy bookkeeping tracer. Therefore Cantera attaches ideal-gas chemical
potential behavior—including composition/partial-pressure dependence—to that
species. Deposited carbon should instead be represented by an appropriate
condensed fixed-stoichiometry phase and/or surface C* transfer to bulk graphite.

This limitation is deliberately exposed in evidence metadata as
`legacy_gas_carbon_tracer`. It is not fixed in this repository version because
that change is being handled separately. Until then, carbon selectivity and
coking conclusions from the reactor mechanism are diagnostic, not validated
physical evidence.

### 11.5 Fail-closed interpretation

`_kinetics_evidence` reports:

- `candidate_specific_kinetics` only when no template or carbon-model
  limitation remains;
- otherwise `diagnostic_screening_template`;
- `can_exclude_candidate = false` for incomplete kinetics or the legacy carbon
  representation.

Thus Cantera can help rank conditions and identify sensitivity targets, but the
current incomplete mechanism cannot eliminate a catalyst or count as measured
reactor validation.

## 12. Methane-conversion pathway modes

`PYROLYSIS_MODE` selects one of six explicitly routed pathway modes:
`thermocatalytic` (the default PFR + fluidized-bed family),
`thermocatalytic_pfr`, `thermocatalytic_fluidized`, `mmbcr`, `ntec`, or
`electrochemical`. Aqueous versus molten electrolyte is an electrochemical
condition, not a separate orchestration mode. NTEC and electrochemical routes
consume strict external-solver artifacts and return non-excluding
`validation_required` evidence when those artifacts are absent or invalid;
they never fall back silently to a thermal reactor.

Thermocatalytic screening uses the unassisted atomistic descriptors. NTEC uses
the same candidate base calculation and optionally applies a bounded transfer
from measured operating conditions and paired NTEC/control calibration.

`pipeline/process/ntec_model.py` requires:

- shear rate;
- interfacial electric field;
- mechanical power per catalyst mass;
- carbon-detachment fraction;
- a field-measurement source;
- paired-control source and count;
- measured barrier reduction;
- measured coking change.

If any required value or provenance is absent, the returned status is `unknown`
and both the barrier reduction and coking bonus are exactly zero. Invalid values
also produce zero assistance. When calibrated, the effect is bounded by the
measured delta and the degree to which the proposed operating conditions are
supported.

This is still a transfer model for a new catalyst. It does not replace a paired
NTEC/control experiment on that catalyst. `--ntec-conditions-json` supplies the
calibration record to the production driver.

## 13. Quantum ESPRESSO DFT and NEB validation

### 13.1 Why Quantum ESPRESSO is used

**Quantum ESPRESSO (QE)** supplies plane-wave, periodic density-functional
theory suitable for catalyst slabs, adsorbates, spin-polarized materials, and
nudged elastic band calculations. It is used to replace ML/BEP estimates for a
small, selected set with candidate-specific electronic-structure calculations.

QE is an external executable suite, not a Python package embedded in this
repository:

- `pw.x` performs SCF and ionic relaxation jobs;
- `neb.x` optimizes reaction paths and climbing images.

ASE is used to construct or read atomistic geometries, while repository code
writes QE inputs, launches executables, and validates outputs.

### 13.2 Portable executable resolution

`pipeline/common/executables.py::resolve_qe_executable` resolves tools in this
order:

1. explicit environment override (`PW_X` or `NEB_X`);
2. the current `PATH`;
3. a PATH-resolved `conda` executable querying the documented `qe-env`;
4. a clear failure if no usable executable exists.

There are no required user-home or absolute installation paths. MPI is similarly
resolved rather than assumed at a machine-specific location.

### 13.3 MPI and resource configuration

`QEExecutionConfig` records MPI ranks, OpenMP threads, k-point pools, and NEB
image groups. `build_qe_command` constructs the argv list rather than a shell
string. The production driver also limits the number of concurrent candidate
jobs and rejects a requested ranks × threads × concurrency product larger than
the visible CPU count.

MPI divides QE work among processes. K-point pools distribute independent
k-points, while NEB image groups can distribute path images. OpenMP divides work
inside each MPI rank. More ranks are not always faster; FFT grid sizes, memory
bandwidth, k-point count, image count, communication, and CPU affinity determine
the optimum. The configuration and execution record make this benchmarkable.

### 13.4 Pseudopotential and cutoff verification

`pipeline/validation/qe_workflows.py::verify_sssp` reads the repository's SSSP
manifest, checks required elements, filenames, checksums, and recommended wave
function/charge-density cutoffs. Input generation fails if verification fails.

The candidate space contains many more elements than a minimal local
pseudopotential collection may cover. Missing coverage is a validation-resource
gap, not evidence against the candidate. Production campaigns must acquire and
verify the needed pseudopotentials before those candidates can complete DFT.

### 13.5 General DFT screening validation

`pipeline/validation/dft_validator.py` generates bulk, slab, and molecule QE
inputs and parses total energy, forces, and convergence. `validate_catalyst`
records whether a run was executed, converged, or merely prepared. Its current
evidence label is `screening_dft`; a generated input or partial output is not a
converged scientific result.

### 13.6 Candidate-specific methane NEB sequence

`methane_dissociation_images` builds an initial path from adsorbed methane toward
CH3* + H*. `write_qe_relax_input` prepares force-converged endpoint relaxations.
`write_qe_neb_input` prepares a spin-polarized multi-image QE NEB job with
climbing-image selection. `run_neb` records command, resources, duration,
return code, and timeout state. `parse_neb_result` accepts a barrier only when
QE reports both normal completion and NEB convergence.

`pipeline/validation/production_workflow.py` now supports both a single NEB
directory and a manifest-driven elementary-kinetics campaign. The declared
fields are methane activation, three subsequent dehydrogenations, hydrogen
desorption, and carbon transport. Each manifest entry supplies explicit,
candidate-specific initial and final structures; the workflow deliberately
does not invent generic reaction geometries for arbitrary materials.

For each supplied elementary step, `advance_methane_neb`:

1. prepares or resumes initial and final `pw.x` relaxations;
2. requires both ionic relaxations to converge;
3. extracts their final geometries;
4. creates an IDPP-interpolated climbing-image path and runs `neb.x`;
5. accepts a barrier only from a converged path;
6. runs prepared central finite-difference force jobs for the proposed
   transition state;
7. builds a mass-weighted Hessian and requires one significant imaginary mode.

`prepare_pyrolysis_campaign` reads the explicit structure manifest and writes
the endpoint and optional frequency-force inputs. `advance_pyrolysis_campaign`
runs/resumes all supplied steps. `pyrolysis_campaign_status` lists unresolved
kinetic fields and emits `converged_dft_neb_frequency` only when every declared
field is resolved. `CandidateKinetics.from_screening_row` will replace Cantera
templates only when that evidence is complete and its candidate ID matches.
On a later production run, `--kinetics-validation-dir` lets Phase 2 locate
`<candidate_id>/pyrolysis_validation.json` and inject a matching completed
campaign into the Cantera mechanism. Missing or incomplete validation leaves
the declared screening templates and diagnostic evidence tier in place.

The proposed transition-state geometry must currently be supplied in the
manifest. Automatically extracting the highest NEB image would not eliminate
the need to inspect/refine that geometry; the explicit input prevents an
uninspected path image from silently becoming frequency evidence.

The path is not complete merely because one or more images printed energies.
Every required electronic calculation, both endpoints, the path, and the
transition-state frequency check must satisfy their contracts.

### 13.7 Transition-state frequency validation

`partial_hessian` constructs a symmetrized, mass-weighted partial Hessian from
central finite differences of forces. It converts eigenvalues to signed
wavenumbers and requires exactly one significant imaginary frequency for a
valid first-order transition state. This is distinct from VQE and from NEB path
convergence: NEB locates a path maximum; frequencies test the local character
of the proposed transition state.

### 13.8 Candidate-specific ORR DFT sequence

`pipeline/validation/dft_fuel_cell.py` and
`pipeline/validation/production_workflow.py::run_orr_sequence` manage clean,
OH*, O*, OOH*, H2, and H2O calculations. `orr_campaign_status` permits a final
ORR result only when every required output is converged. Corrected adsorption
free energies and overpotential are then computed with explicit correction
metadata.

The sequence is resumable. Completed stages are preserved; missing stages can
advance; incomplete nonempty output is not overwritten by default.

For production refinement, `build_orr_validation_plan` deterministically
enumerates atop, bridge, and hollow trial sites across declared coverages and
OH/O/OOH adsorbates. It labels every task as requiring an explicitly realized
structure: a coverage label alone is never presented as a physical supercell.
`evaluate_orr_ensemble` admits only site/coverage pathways where all three
adsorbates converged, applies one or more provenance-bearing solvation,
potential, pH, and temperature correction models, and reports the spread across
all cases as model uncertainty. Missing expected cases suppress the headline
overpotential and keep the ensemble evidence `incomplete`.

### 13.9 Persistent validation task queue

`pipeline/validation/task_queue.py::ValidationTaskQueue` stores pending,
running, and completed candidate jobs in SQLite. The production driver uses it
to coordinate multiple candidate calculations while respecting
`--qe-max-concurrent`. This separates candidate-level parallelism from MPI
parallelism inside each QE calculation.

## 14. CUDA-Q, CUDA Quantum, and VQE

### 14.1 Why a quantum-computing stage exists

Catalyst chemistry is an interacting-electron problem. Conventional DFT is the
practical high-fidelity method used elsewhere in this repository, but some
small, strongly correlated active spaces can be difficult for a single
mean-field-like electronic description. A future quantum workflow could take a
carefully selected active space around a bond-breaking event—such as C-H
cleavage during methane activation or O-O cleavage in OOH*—and estimate its
correlated electronic energy.

That is why the project contains a VQE path. Its intended long-term role is
**targeted refinement of a very small number of finalists**, not traversal of
the 21.1-billion-member population. Even an ideal quantum solver would sit near
the expensive end of the fidelity funnel:

```text
branch search -> eSen screen -> DFT geometry/path -> active-space construction
     -> VQE correlated energies -> compare with classical reference
     -> refine an energy difference or uncertainty for a finalist
```

VQE would be run separately for the reactant, product, or transition-state
geometries required by a defined quantity. An isolated VQE ground-state energy
is not itself an activation barrier. A barrier requires a consistently defined
energy difference, and kinetics additionally require the relevant thermal,
entropic, and dynamical treatment.

The built-in path establishes only the solver portion of that chain. The
candidate-specific boundary accepts externally generated standard FCIDUMP
integrals through `pipeline/validation/candidate_hamiltonian.py`. PySCF restores
the spatial integral tensor and OpenFermion performs the spin-orbital expansion
and Jordan-Wigner transform. A checksum-bound sidecar must identify the exact
geometry, electronic-structure protocol, basis, active orbitals/electrons,
charge, multiplicity, frozen orbitals, and integral source. The repository does
not infer these scientific choices from a plane-wave QE output.

### 14.2 Terminology and division of responsibilities

**CUDA-Q** is NVIDIA's hybrid quantum programming platform. “CUDA Quantum” is a
name used for the same platform, not a separate second solver in this pipeline.

The three relevant layers are different:

| Layer | Responsibility in or beneath this project |
|---|---|
| **VQE** | The hybrid quantum-classical algorithm: prepare a parameterized quantum state, measure the Hamiltonian expectation value, and use a classical optimizer to lower that energy. |
| **CUDA-Q / CUDA Quantum** | The programming and execution layer: define the quantum kernel, Pauli-spin Hamiltonian, target backend, expectation-value evaluation, and optimizer call through one Python API. |
| **cuQuantum** | NVIDIA's lower-level GPU libraries for accelerated state-vector/tensor-network quantum simulation. A CUDA-Q simulator may use this acceleration internally; this repository does not import or call cuQuantum directly. |

This means the repository should not say that CUDA-Q and CUDA Quantum are two
independent solvers, or that it explicitly schedules cuQuantum kernels. It asks
CUDA-Q to execute a quantum kernel on the selected target. Backend implementation
and any cuQuantum acceleration are owned by the installed NVIDIA software
stack.

The integration lives entirely in
`pipeline/validation/vqe_transition_state.py`.

### 14.3 Why CUDA-Q is used instead of a custom state-vector implementation

CUDA-Q supplies the plumbing needed to express and run a hybrid algorithm
without coupling the scientific workflow to one simulator implementation. In
this module it provides:

- `cudaq.kernel` and `cudaq.qvector` for defining the parameterized circuit;
- quantum gates (`x`, `ry`, `rz`, and `cx`) for state preparation and the
  hardware-efficient ansatz;
- `cudaq.spin` operators for translating Pauli words into an executable
  Hamiltonian;
- `cudaq.set_target` for choosing the NVIDIA GPU simulator or CPU backend;
- `cudaq.vqe` for the repeated expectation-value/optimization loop;
- `cudaq.optimizers.COBYLA` for the classical parameter update.

The architectural benefit is backend portability: the surrounding pipeline can
retain the same Hamiltonian, ansatz, result schema, and evidence checks while a
CUDA-Q target changes from a local CPU simulator to an NVIDIA GPU simulator or,
after appropriate validation, another supported execution target. The present
production driver deliberately asks for `target='nvidia'` so compatible GPU
hardware accelerates repeated circuit simulation and Hamiltonian expectation
evaluation.

CUDA-Q does not discover the catalyst, generate the geometry, select the active
orbitals, or decide whether an energy is chemically meaningful. Those remain
scientific responsibilities of the upstream workflow.

### 14.4 What cuQuantum contributes

A classical simulation of `n` qubits generally stores a state vector with
`2^n` complex amplitudes. Applying gates and evaluating many Pauli expectation
values therefore becomes increasingly expensive as the active space grows.
NVIDIA's cuQuantum libraries are designed to accelerate those underlying
linear-algebra and tensor-network operations on NVIDIA GPUs.

For this repository, cuQuantum's purpose is consequently **execution
acceleration beneath a compatible CUDA-Q simulation target**. It can reduce the
wall time of the many circuit evaluations made during VQE and use GPU memory and
throughput more effectively than a simple Python/CPU simulator. It does not:

- improve the physical accuracy of an incorrect Hamiltonian;
- turn a representative Hamiltonian into a candidate-specific one;
- remove the exponential scaling of an exact state-vector simulation;
- replace DFT geometry optimization, NEB, frequency analysis, Cantera, or
  experiment;
- provide a separate result consumed directly by another pipeline stage.

The code intentionally has no `import cuquantum`. Whether the installed CUDA-Q
`nvidia` target uses a particular cuQuantum backend is a deployment detail and
must be verified from that environment when reporting performance. The only
claim supported by this repository itself is that it requests the CUDA-Q
NVIDIA target.

### 14.5 What VQE does step by step

VQE minimizes the Rayleigh quotient

```text
E(theta) = <psi(theta) | H | psi(theta)>
E(theta) >= E0
```

where `H` is the qubit Hamiltonian, `|psi(theta)>` is the parameterized ansatz,
and `E0` is the exact ground-state energy in the represented space. The
variational inequality is an important correctness invariant: a noiseless VQE
energy below the exact ground-state energy indicates a solver, convention, or
Hamiltonian error.

`run_vqe` implements the following sequence:

1. It imports `cudaq` and `cudaq.spin`. If CUDA-Q is unavailable, it emits an
   explicitly marked mock result rather than pretending the solver ran.
2. It maps `target='default'` to the CUDA-Q `qpp-cpu` simulator and otherwise
   selects the requested target, normally `nvidia`.
3. It converts each `(coefficient, Pauli_word)` pair into a CUDA-Q spin operator
   and sums the terms into `H`.
4. It allocates four qubits for the current model and initializes a half-filled
   computational-basis reference with X gates.
5. For each of three layers it applies parameterized Ry and Rz rotations and a
   ring of CNOT entanglers. Four qubits × two rotations × three layers gives 24
   variational parameters.
6. It initializes those parameters and asks CUDA-Q's COBYLA optimizer to
   minimize the measured expectation value, with at most 3,000 iterations.
7. It independently constructs the complete small Hamiltonian matrix with
   NumPy and diagonalizes it with `eigvalsh`.
8. It records the VQE energy, exact energy, variational gap, target, parameters,
   and benchmark outcome.

COBYLA is used because it is derivative-free: the optimizer needs only energy
evaluations, which matches how a quantum backend exposes expectation values.
It is not guaranteed to find the global optimum; ansatz expressivity,
initialization, optimizer settings, sampling noise, and barren plateaus can all
affect convergence.

### 14.6 Hamiltonians implemented today

The module defines small Pauli Hamiltonians for a model C-H splitting or O-O
cleavage active space:

- `build_ch_splitting_hamiltonian` represents a conceptual sigma(C-H),
  sigma*(C-H), and metal-d interaction model;
- `build_orr_hamiltonian` represents a conceptual sigma(O-O), pi*(O-O), and
  metal-O interaction model.

Both return four-qubit lists of weighted Pauli strings. The identity term
contains a representative core offset, while Z, X, and Y products encode model
orbital energies, couplings, exchange, and interaction terms. These constants
exercise the complete Hamiltonian-to-circuit-to-optimizer path, but they were
not generated from the named catalyst passed to `validate_transition_state`.
The catalyst name currently labels the output artifact; it does not change the
Hamiltonian coefficients.

For small Hamiltonians, `exact_ground_energy` constructs the full matrix and
classically diagonalizes it. The VQE result must respect the variational lower
bound and can be checked against a 1.6e-3 Hartree chemical-accuracy tolerance.

### 14.7 What purpose this stage serves today

In the current repository, the quantum stage serves four concrete engineering
and research-preparation purposes:

1. **Integration validation.** It proves that the isolated `quantum-env`,
   CUDA-Q import, target selection, kernel compilation, GPU execution, optimizer
   contract, JSON serialization, and campaign handoff can operate together.
2. **Numerical contract testing.** Classical exact diagonalization checks
   Hermiticity, the variational lower bound, and the energy gap for the small
   model, catching Pauli-word, qubit-order, or solver-interface errors.
3. **Resource-path validation.** It exercises the NVIDIA quantum-simulation
   route separately from PyTorch/eSen, allowing GPU/runtime problems to be found
   before candidate-specific quantum calculations become expensive.
4. **Architecture reservation.** It defines the result and evidence fields that
   a future candidate-specific active-space generator must satisfy, without
   granting current model calculations scientific authority they do not have.

The stage currently serves **no direct candidate-selection purpose**. It does
not refine `E_act`, ORR overpotential, Cantera kinetics, PEMFC power, or branch
priority. In fact, `run_production_campaign.py` rejects its result as campaign
validation unless both `catalyst_specific_hamiltonian` and `benchmarked` are
true. The current implementation sets
`catalyst_specific_hamiltonian = false`, so the fail-closed production path will
skip/fail this phase rather than promote the output.

### 14.8 Scientific limitation and required upgrade path

The built-in Hamiltonian coefficients remain representative constants rather
than candidate-specific integrals. A successful built-in GPU run is therefore
labeled `toy_hamiltonian`, and a run without CUDA-Q is labeled `mock`. A sourced
FCIDUMP can instead be labeled `candidate_specific_integrals`; its VQE result is
only `candidate_specific_vqe_benchmarked` after the solver passes the classical
reference tolerance. None of these isolated energies is an activation barrier,
and none replaces DFT/NEB.

VQE becomes scientifically relevant to a candidate only after an upstream
electronic-structure workflow supplies:

1. converged candidate and reaction geometries from DFT/NEB;
2. a reproducible orbital basis and active-space selection rule;
3. candidate-specific one- and two-electron integrals;
4. charge, electron count, spin/multiplicity, frozen-core choices, and reference
   energy;
5. a documented fermion-to-qubit mapping and any symmetry tapering;
6. separate, consistent calculations whose energy difference represents the
   desired reactant, transition-state, or ORR quantity;
7. comparison against a trusted classical calculation for every tractable
   instance;
8. an uncertainty and convergence protocol appropriate to simulator or hardware
   noise.

Only then should the result set `catalyst_specific_hamiltonian = true`, enter
the adaptive-validation ledger, or influence a finalist. Until that upgrade,
this stage validates software and solver behavior and prepares a future quantum
research path; the actual catalyst evidence comes from eSen screening,
converged Quantum ESPRESSO calculations, and physical measurements.

## 15. Fuel-cell cathode, cell, and stack workflow

### 15.1 Cathode screening

Fuel-cell branch discovery uses the same indexed classes and coverage policy,
subject to `pipeline/common/application_scope.py::pemfc_cathode_scope`.
Chemistries that are not meaningful as a solid PEMFC cathode are marked out of
scope with an explicit reason instead of being mis-scored.

`pipeline/screening/fc_screener.py` performs class-spanning eSen/CHE screening.
`pipeline/screening/fc_cathode_screener.py` provides the direct cathode and
membrane workflow used by the compact orchestrator. Candidate records contain
adsorption/free-energy descriptors, overpotential, confidence, safety, cost,
and validation requirements.

### 15.2 PEMFC model

`pipeline/process/pemfc_model.py` implements a one-dimensional, lumped
through-MEA electrochemical model. For each current density it computes:

```text
Vcell = Enernst - eta_cathode - eta_anode - eta_ohmic - eta_transport
power density = current density * Vcell
```

Components include:

- temperature/pressure-adjusted Nernst voltage;
- Tafel cathode activation loss derived from the ORR overpotential;
- Butler-Volmer/asinh anode activation loss;
- membrane, electronic, contact, and catalyst-layer resistance;
- GDL oxygen-diffusion limiting current;
- CO, sulfur, and hydrogen-purity effects on the anode;
- optional measured voltage-degradation input;
- class-specific Tafel-slope assumptions;
- membrane property sweeps.

Outputs include the full polarization and power curves, peak and rated power,
efficiency, limiting current, resistance, impurity factor, and evidence labels.
The result is `modeled` and `requires_mea_validation = true`.

### 15.3 Stack model

`pipeline/process/fuel_cell_stack.py` scales a selected cell operating point to
an N-cell stack. It estimates gross and net power, compressor/pump/blower/control
parasitics, heat rejection, radiator area, hydrogen consumption, stack/system
efficiency, mass, volume, gravimetric/volumetric power, catalyst and membrane
cost, total cost, and cost per kW.

These are engineering screening estimates. Stack design, water management,
thermal transients, mechanical compression, controls, manufacturability, and
durability require dedicated models and physical validation.

### 15.4 End-to-end hydrogen quality

The process and fuel-cell stages are conceptually connected through hydrogen
purity and impurity inputs, not through a complete dynamic plant flowsheet.
Reactor conversion/selectivity and hydrogen cleanup specifications must
eventually be converted into measured or validated inlet compositions for the
PEMFC model. The current impurity model supports sensitivity analysis and makes
the required testing explicit.

## 16. Techno-economics

`pipeline/process/tea.py` converts modeled conversion and operating assumptions
into hydrogen-cost scenarios. It explicitly labels results
`screening_scenario_not_measured_tea` or `screening_sensitivity_range`.

These estimates are useful for rejecting obviously unattractive operating
regimes and identifying which measurement most affects economics. They are not
a bankable process design, vendor quote, or independently reviewed TEA.

## 17. Prior art, novelty, and claim readiness

### 17.1 Prior-art registry

`pipeline/evidence/prior_art.py::PriorArtRegistry` stores canonical candidate
identities, source type/ID, citation, evidence level, and publication year in
SQLite. It can import curated CSV sources and annotate screening frames.

An empty registry cannot establish novelty. Candidate uniqueness inside this
generated space also cannot establish novelty; novelty is relative to public
knowledge and earlier disclosure.

### 17.2 Time-split novelty benchmark

`pipeline/evidence/novelty_benchmark.py::time_split_recovery` freezes a
historical cutoff and tests whether rankings recover discoveries held out from
the earlier training/prior-art snapshot. Manifest hashing makes the input split
auditable. A successful blinded time split measures prospective ranking value
more credibly than re-ranking examples already known to the system.

### 17.3 Pilot benchmark

`pipeline/evidence/pilot_benchmark.py` compares divide-and-conquer, seeded
random sampling, and deterministic expert-style scoring over compatible legacy
outcomes. It reports enrichment and uncertainty statistics. Because the
outcomes are legacy computational screening values, this is methodological
evidence rather than a prospective discovery result.

### 17.4 Hash-verified evidence manifest

`pipeline/evidence/manifest.py` requires each evidence record to name a
candidate, source artifact, SHA-256 checksum, protocol ID, and exact status.
Recognized categories include:

- converged DFT and ORR DFT;
- paired NTEC/control measurements;
- measured reactor performance and deactivation;
- measured MEA and durability performance;
- hydrogen impurity tests;
- valid time-split benchmark evidence;
- curated prior-art sources.

Changing a source artifact invalidates its checksum. Merely entering a count in
JSON is insufficient.

### 17.5 Final readiness

`pipeline/evidence/readiness.py::campaign_readiness` checks the correct
population denominator, complete coverage when requested, nonempty prior-art
registry, and verified evidence categories.

For turquoise hydrogen it requires converged DFT, measured reactor performance,
and measured deactivation; NTEC mode additionally requires a measured paired
control. For fuel cells it requires converged ORR DFT, measured MEA performance,
durability, impurity testing, a valid time-split benchmark, and curated prior
art.

`--final-campaign` activates these fail-closed gates. The pipeline may still
produce useful intermediate artifacts when not ready, but must not call them a
completed discovery.

## 18. Evidence ladder and allowed language

| Level | Typical artifact | What it supports | What it does not support |
|---|---|---|---|
| Enumeration | candidate ID, coverage record | Candidate exists in declared space | Activity or novelty |
| Surrogate ranking | tree/neural score, uncertainty | Search priority | Physical performance |
| ML atomistic screening | eSen energies, forces, CHE/BEP descriptors | Shortlisting and DFT allocation | Converged DFT or experiment |
| Diagnostic process model | incomplete Cantera mechanism, PEMFC model | Sensitivity and system-level prioritization | Measured conversion, power, durability |
| Converged candidate DFT | verified QE outputs | Electronic/ionic energies under stated model | Reactor or MEA performance |
| Converged NEB/frequencies | barrier and one imaginary mode | Computed pathway evidence | Measured kinetics or durability |
| Candidate-specific complete kinetics | all elementary parameters with provenance | Quantitative microkinetic simulation | Reactor validation by itself |
| Measurement | reactor, deactivation, NTEC/control, MEA, impurity, durability | Performance under recorded protocol | Universal performance outside tested domain |
| Novelty validation | curated prior art + blinded time split | Evidence of novelty/search value | Patentability or commercial freedom to operate |

Always preserve the qualifier. “Predicted low barrier” is accurate for a BEP or
ML result. “Converged DFT barrier” requires a completed path. “High methane
conversion” without “modeled” requires reactor measurements.

## 19. Runtime environments and dependency boundaries

The repository documents multiple environments because the scientific stacks
have different binary, CUDA, and compiler requirements:

| Logical environment | Main purpose | Key software |
|---|---|---|
| `deepmd-env` | GPU atomistic screening | PyTorch, fairchem/eSen, ASE, pandas |
| `cp2k-env` | Reactor simulation (historical name) | Cantera 3.x |
| `qe-env` | Periodic DFT and NEB | Quantum ESPRESSO, MPI, ASE helpers |
| `quantum-env` | Quantum workflow | CUDA-Q and a compatible NVIDIA stack |
| `battery-env` | Fuel-cell/data utilities | NumPy, SciPy, pandas, ASE/pymatgen as needed |
| `openfoam-env` | Fluidized/MMBCR/NTEC hydrodynamics | OpenFOAM `multiphaseEulerFoam` |
| `fenicsx-env` | NTEC/electrochemical continuum transport | FEniCSx/dolfinx, MPI, Cantera coupling |

The environment names are conventions, not absolute paths. `run_in_env` and
the executable resolver locate `conda` from `PATH`. Users may instead provide
the required executable in `PATH` or set documented overrides.

`environment.yml` provides a useful core Python environment;
`environment-openfoam.yml` and `environment-fenicsx.yml` provide reproducible
external multiphysics environments. The runner discovers these environments,
`PATH`, or documented executable overrides without assuming a home directory.
It accepts solver output only after artifact identity, backend version,
convergence, mesh-independence, conservation, required-output, and provenance
validation. `requirements.txt` identifies phase-specific dependencies and
comments on external executables. GPU driver, CUDA, MPI, Quantum ESPRESSO,
pseudopotential, and CUDA-Q compatibility must still be verified on the target
machine.

## 20. Artifacts and provenance

All base paths are derived from the repository location in
`pipeline/common/utils.py`; no user-home layout is assumed.

Typical generated artifacts are:

```text
results/
  design_space_audit.json
  prior_art.sqlite
  screening/
    indexed_scan.sqlite
    turquoise_hydrogen_coverage_certificate.json
    branch_ranker_evidence.csv
    ga_full_database.csv
  reactor/
    *.json
  dft/
    validation_tasks.sqlite
    *.in / *.out / *.execution.json / *_dft.json
  vqe/
    vqe_*.json
  fuel_cell/
    indexed_scan.sqlite
    coverage_certificate.json
    cathode_screening.csv
    pemfc_*.json
    stack_*.json
  reports/
    pipeline_state.json
    campaign_readiness.json

mechanisms/
  *.yaml
  *.kinetics.json
```

Exact files depend on the chosen entry point and completed stages. Provenance
includes CLI arguments, Git SHA, timestamp, GPU identities, environment
snapshot, candidate IDs, protocol IDs, execution resources, and artifact
hashes where the evidence manifest requires them.

## 21. Failure and restart semantics

The pipeline distinguishes operational failure from scientific outcome:

- a worker crash requeues leased screening work;
- a time-limited branch run leaves pending intervals and an incomplete
  certificate, then resumes from SQLite;
- an unconverged relaxation becomes a validation request;
- a missing pseudopotential blocks that calculation, not that chemistry;
- a nonzero or incomplete QE output is not parsed as converged evidence;
- an existing incomplete output is preserved by default;
- a missing Cantera or CUDA-Q installation may produce mock output only where
  the caller explicitly permits it;
- production reactor execution forbids mock output;
- a mock or toy result carries an explicit evidence label;
- final-campaign mode refuses to pass with missing evidence.

This behavior is as important as the numerical formulas: it prevents equipment
or software limitations from masquerading as negative catalyst results.

## 22. Running and interpreting the workflow

### 22.1 Integration check

```bash
python -m pipeline.orchestrator --quick --no-dft --no-vqe
```

This exercises reduced workflow plumbing. It does not complete search coverage
or high-fidelity validation.

### 22.2 Staged production search

```bash
python run_production_campaign.py \
  --hours 4 \
  --branch-max-leaves 8 \
  --scan-workers 8 \
  --validation-batch 500 \
  --min-validation-per-class 2 \
  --mode thermocatalytic \
  --no-dft \
  --no-vqe
```

This advances persistent search state for a bounded time and can be rerun. The
resulting coverage certificate tells the reader whether traversal is complete.
Removing `--no-dft` authorizes the configured validation execution; it does not
guarantee convergence.

### 22.3 NTEC run

First populate a sourced solver case as specified in
[`PHYSICAL_CASES.md`](PHYSICAL_CASES.md), then generate its validated artifact:

```bash
python -m pipeline.process.multiphysics_runner \
  --mode ntec \
  --reactor-type NTEC \
  --candidate-id CANONICAL_ID \
  --temperature-K 300 \
  --case-dir cases/CANONICAL_ID/ntec \
  --fenics-model cases/CANONICAL_ID/ntec/model.py \
  --results-dir results/multiphysics \
  --model-source "immutable case revision or DOI"
```

The case must emit converged, conservative, mesh-independent outputs and
held-out validation metadata. Then run the campaign against that artifact root:

```bash
python run_production_campaign.py \
  --mode ntec \
  --ntec-conditions-json measured_ntec_conditions.json \
  --multiphysics-results-dir results/multiphysics \
  --hours 4
```

The JSON must contain measured operating conditions and paired-control
provenance. Without those measurements NTEC screening bonuses remain zero;
without the validated artifact the reactor condition remains non-excluding
`validation_required`.

### 22.4 Inspect or advance QE validation

```bash
python run_validation_campaign.py --pyro-dir results/dft/example

python run_validation_campaign.py \
  --pyro-dir results/dft/example \
  --advance \
  --mpi-ranks 4 \
  --omp-threads 1 \
  --neb-image-groups 1
```

Use the status-only form before restarting work. Do not use
`--restart-incomplete` while another scheduler or process may own the output.

A complete elementary-step campaign is described by JSON:

```json
{
  "candidate_id": "stable-candidate-id",
  "campaign_dir": "candidate_calculations",
  "steps": {
    "methane_activation": {
      "initial": "structures/ch4_adsorbed.traj",
      "final": "structures/ch3_h_adsorbed.traj",
      "transition_state": "structures/ch4_ts.traj",
      "n_images": 7,
      "displacement_A": 0.01
    }
  }
}
```

Prepare and then advance it with:

```bash
python run_validation_campaign.py --pyro-manifest campaign.json --prepare
python run_validation_campaign.py --pyro-manifest campaign.json --advance
```

Omitted elementary steps remain unresolved. They are never filled with a
successful result from a different step.

Generate an ORR site/coverage plan and evaluate completed free-energy records
with provenance-bearing correction models using:

```bash
python run_validation_campaign.py \
  --orr-structure relaxed_surface.traj \
  --orr-coverages 0.25,0.5,1.0 \
  --orr-plan-output results/dft/orr_plan.json

python run_validation_campaign.py \
  --orr-ensemble-results results/dft/orr_site_results.json \
  --orr-corrections-json orr_corrections.json \
  --orr-expected-cases 18
```

The first command creates calculation tasks; candidate-specific supercells and
adsorbate structures must realize those tasks before QE execution. The second
command refuses a headline ORR overpotential unless every expected case is
present and converged.

### 22.5 Final readiness check

```bash
python run_production_campaign.py \
  --final-campaign \
  --evidence-manifest results/evidence_manifest.json \
  --prior-art-db results/prior_art.sqlite
```

This should fail until search coverage, prior art, calculations, and physical
evidence are actually present. Failure is expected and useful during campaign
development.

## 23. How to extend the system safely

### Add a catalyst class

1. Add its axes, generator, validation, encoding, and size calculation in
   `pipeline/common/catalyst_spaces.py`.
2. Add deterministic inverse mapping in `candidate_at_class`.
3. Add canonical-count logic and provenance in the design-space audit.
4. Add a physically meaningful ASE screening representation.
5. Define application scope and extraction of constituent elements.
6. Update class-preserving tests, raw denominator, and expected certificates.
7. Never change the population denominator without a versioned migration of
   persistent scan state.

### Add a screening descriptor

1. Define its physical meaning, units, reference state, and protocol version.
2. Record convergence/censoring and provenance alongside the scalar value.
3. Decide whether it is a priority feature, quantitative stage requirement, or
   validation-only observation.
4. Add equation and limiting-case tests.
5. Do not let a missing optional descriptor become a hidden hard rejection.

### Replace a template kinetic barrier

1. Add the candidate-specific field to `CandidateKinetics`.
2. Populate it only from a source with an explicit protocol and candidate ID.
3. Preserve units and convert eV per event to J/mol exactly once.
4. Update the YAML reaction and kinetics sidecar.
5. Remove the relevant `template_default` limitation only after the new value
   is present and validated.

### Add a new high-fidelity backend

1. Keep executable discovery portable.
2. write inputs and execution metadata before accepting outputs;
3. distinguish prepared, running, incomplete, failed, and converged states;
4. require candidate identity and geometry hashes;
5. connect results through `record_validation` so calibration can learn from
   disagreement;
6. make absent software a calculation blockage, not a candidate failure.

### Make VQE candidate-specific

1. Generate candidate-specific orbitals and one-/two-electron integrals.
2. Document active-space, charge, multiplicity, and frozen orbitals.
3. Map the fermionic Hamiltonian reproducibly and retain the constant energy.
4. check symmetries, electron number, Hermiticity, and tapering;
5. benchmark against exact diagonalization or a trusted classical method while
   the active space permits it;
6. propagate uncertainty and never substitute VQE energy for an activation
   barrier without a defined reactant/transition-state energy difference.

## 24. Known limitations and remaining campaign work

The architecture is suitable for a rigorous search and validation campaign,
but repository capability is not the same thing as completed evidence. The
following remain essential:

- finish all unresolved branch leaves for complete address-space coverage;
- populate and curate prior art rather than relying on an empty registry;
- ensure verified pseudopotential coverage for every finalist's elements;
- relax candidate-specific endpoints and converge all NEB images;
- establish converged methane barriers and transition-state frequencies;
- converge clean/OH/O/OOH/H2/H2O ORR calculations for finalists;
- replace every material kinetic template needed for quantitative Cantera use;
- replace the legacy gas carbon tracer with a validated surface/condensed
  carbon treatment in the separately owned mechanism change;
- collect conversion, selectivity, carbon yield, energy, and deactivation data;
- collect paired NTEC/control calibration and candidate measurements;
- collect MEA power, impurity tolerance, and durability data;
- run a successful blinded time-split novelty benchmark;
- attach all evidence through immutable, checksum-verified manifests.

Until those tasks are complete, the defensible output is a prioritized,
diverse, reproducible set of candidates and validation plans—not a declaration
that a state-of-the-art material has been discovered.

## 25. Glossary

- **ASE** — Atomic Simulation Environment; Python objects and tools for atomic
  structures, calculators, optimization, and simulation I/O.
- **BEP relation** — Brønsted–Evans–Polanyi relationship used to estimate an
  activation barrier from reaction energy during screening.
- **BFGS** — Quasi-Newton geometry optimizer; here used for atomic relaxation.
- **CHE** — Computational hydrogen electrode; converts adsorbate energies and
  corrections into electrochemical free-energy steps.
- **CSTR** — Continuously stirred-tank reactor; a cascade approximates axial
  progression in the bubble-column model.
- **CUDA-Q / CUDA Quantum** — NVIDIA hybrid quantum programming platform used
  here for an optional VQE solver workflow.
- **cuQuantum** — NVIDIA GPU quantum-simulation primitives that can support
  CUDA-Q targets; not called directly by this code.
- **DFT** — Density-functional theory; the periodic electronic-structure method
  run through Quantum ESPRESSO.
- **eSen** — Meta/fairchem learned interatomic model used for fast atomistic
  screening energies and forces.
- **FIRE** — Inertial relaxation optimizer used as recovery after BFGS.
- **MEA** — Membrane electrode assembly in a PEM fuel cell.
- **MPI** — Message Passing Interface; process-level parallelism used by QE.
- **NEB** — Nudged elastic band; a chain-of-images method for reaction paths and
  activation barriers.
- **NTEC** — The project's paired-control-calibrated non-thermal enhancement
  mode for methane pyrolysis; speculative benefits default to zero.
- **OOD** — Out of distribution; chemistry insufficiently represented by model
  training/support evidence.
- **ORR** — Oxygen reduction reaction at the fuel-cell cathode.
- **PEMFC** — Proton-exchange-membrane fuel cell.
- **QE** — Quantum ESPRESSO.
- **SSSP** — Standard Solid State Pseudopotentials library/manifest used to
  verify QE inputs and recommended cutoffs.
- **TST** — Transition-state theory.
- **VQE** — Variational quantum eigensolver.

## 26. One-sentence interpretation of every major tool

- **Python/NumPy/SciPy/pandas:** data movement, numerical formulas, statistics,
  optimization support, and tabular artifacts.
- **SQLite:** durable search intervals, archives, calibration observations,
  prior art, experimental slates, and validation task state.
- **ASE:** the shared representation and manipulation layer for atomistic
  structures.
- **Meta eSen through fairchem:** fast learned energies and forces used to
  decide where expensive calculations are worth spending.
- **PyTorch/CUDA:** GPU execution of learned screening and surrogate models.
- **Quantum ESPRESSO:** candidate-specific periodic DFT, relaxation, ORR, and
  NEB calculations that can upgrade screening estimates when converged.
- **MPI/OpenMP:** parallel execution inside QE calculations; independent task
  concurrency is managed separately.
- **Cantera:** reaction-network and reactor-scale propagation of declared
  kinetic/thermodynamic inputs into diagnostic conversion/selectivity behavior.
- **CUDA-Q/CUDA Quantum:** optional VQE execution and solver verification for a
  future candidate-specific quantum chemistry path; currently model-level.
- **Checksums and Git provenance:** evidence integrity and reproducibility, not
  physical validation by themselves.

The system's competitive advantage is therefore not a single solver. It is the
combination of exact population accounting, diverse deterministic traversal,
feedback-directed allocation, multi-fidelity physics, explicit uncertainty,
portable high-performance execution, and refusal to upgrade incomplete work
into a stronger claim than the evidence supports.
