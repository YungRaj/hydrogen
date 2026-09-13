# Architecture and workflow visual atlas

These Mermaid diagrams describe what the engine executes, which tools own which
physics, how evidence gains authority, and where components may be replaced.
They render directly on GitHub while remaining reviewable, versioned text.

Read each diagram from top to bottom unless its arrows run left to right. Solid
arrows show execution or data flow; dashed arrows show feedback or a planned
extension. Blue denotes search or control, green denotes scientific execution,
amber denotes evidence, and red denotes a decision or rejection gate.

## End-to-end catalyst discovery engine

```mermaid
flowchart TB
    subgraph EXPLORE["1 · Explore the design space"]
        direction LR
        SPACE["21.1B indexed candidates<br/>across 14 material classes"]
        SEARCH["Coverage-guided<br/>branch-and-bound search"]
        SCREEN["Fast physical and<br/>machine-learned screening"]
        ARCHIVE["Diverse Pareto archive<br/>and regional champions"]
        SPACE --> SEARCH --> SCREEN --> ARCHIVE
    end

    subgraph VALIDATE["2 · Evaluate independent scientific objectives"]
        direction LR
        REACTOR["Methane conversion<br/>and reactor behavior"]
        ATOMISTIC["Quantum ESPRESSO<br/>relaxation, NEB, frequencies"]
        POWER["ORR, PEMFC<br/>power and durability"]
    end

    subgraph EVIDENCE["3 · Establish evidence and decide what runs next"]
        direction LR
        NOVELTY["Prior-art and<br/>novelty checks"]
        EXPERIMENT["Experimental calibration<br/>and blinded holdouts"]
        REPORT["Readiness report<br/>with claim gates"]
        NOVELTY --> REPORT
        EXPERIMENT --> REPORT
    end

    ARCHIVE --> REACTOR --> REPORT
    ARCHIVE --> ATOMISTIC --> REPORT
    ARCHIVE --> POWER --> REPORT
    ARCHIVE --> NOVELTY
    REPORT -. "uncertainty and disagreement guide the next round" .-> SEARCH

    classDef explore fill:#e8f1ff,stroke:#2563eb,color:#172554
    classDef validate fill:#e8f8ee,stroke:#15803d,color:#052e16
    classDef evidence fill:#fff7dc,stroke:#b45309,color:#451a03
    class SPACE,SEARCH,SCREEN,ARCHIVE explore
    class REACTOR,ATOMISTIC,POWER validate
    class NOVELTY,EXPERIMENT,REPORT evidence
```

The archive deliberately fans out. Methane conversion, ORR, novelty, and
experimental readiness are distinct objectives; no single score proves all of
them.

## Divide-and-conquer traversal

```mermaid
flowchart TB
    ROOT["Index the complete design space"]
    PARTITION["Partition by material class<br/>and chemical region"]
    PROBE["Probe every region with a<br/>deterministic low-discrepancy schedule"]
    PRIORITY{"Does this branch show<br/>promise or uncertainty?"}
    REFINE["Subdivide and evaluate sooner"]
    DEFER["Lower its priority<br/>but retain coverage"]
    FLOOR["Apply the fixed regional budget"]
    LEAF["Run resumable terminal-leaf scans"]
    RESULTS["Update regional champions<br/>and the Pareto archive"]
    CERTIFICATE["Issue a coverage certificate<br/>with visited intervals and gaps"]

    ROOT --> PARTITION --> PROBE --> PRIORITY
    PRIORITY -- "Yes" --> REFINE --> LEAF
    PRIORITY -- "Not yet" --> DEFER --> FLOOR --> LEAF
    LEAF --> RESULTS
    LEAF --> CERTIFICATE
    RESULTS -. "calibration feedback" .-> PROBE

    classDef structure fill:#eef2ff,stroke:#4338ca,color:#1e1b4b
    classDef decision fill:#fff1f2,stroke:#be123c,color:#4c0519
    classDef execution fill:#ecfdf5,stroke:#047857,color:#022c22
    classDef proof fill:#fffbeb,stroke:#b45309,color:#451a03
    class ROOT,PARTITION,PROBE structure
    class PRIORITY,REFINE,DEFER,FLOOR decision
    class LEAF execution
    class RESULTS,CERTIFICATE proof
```

Priority controls ordering and discretionary compute, not exclusion. Each
region retains a floor, and the certificate records what remains unvisited.

## Multi-fidelity learning and referral loop

```mermaid
flowchart LR
    subgraph CALIBRATE["A · Offline calibration"]
        direction TB
        DESIGN["Design representative cases<br/>including regime anchors"]
        SPLIT["Assign training and blind validation<br/>before execution"]
        FULL["Run the required<br/>full-physics solvers"]
        ARTIFACT{"Does the artifact pass identity,<br/>convergence and balance checks?"}
        TRAIN["Train one surrogate<br/>for one mode and reactor"]
        HOLDOUT{"Does blind error satisfy<br/>the declared RMSE limit?"}
        REGISTRY["Publish the model and<br/>checksum manifest atomically"]
        DESIGN --> SPLIT --> FULL --> ARTIFACT
        ARTIFACT -- "Yes" --> TRAIN --> HOLDOUT
        HOLDOUT -- "Yes" --> REGISTRY
    end

    subgraph INFERENCE["B · Online screening"]
        direction TB
        QUERY["Request a reduced-reactor closure"]
        DOMAIN{"Do identity, schema, domain,<br/>uncertainty and physics checks pass?"}
        CLOSURE["Use the calibrated closure<br/>in the reduced Cantera model"]
        QUERY --> DOMAIN
        DOMAIN -- "Yes" --> CLOSURE
    end

    subgraph REFINE["C · Evidence-directed refinement"]
        direction TB
        REFER["Require a full-physics calculation"]
        SCHEDULE["Reserve regional coverage,<br/>then apply adaptive priority"]
        LEDGER["Record inputs, decisions and<br/>model identity in the hash chain"]
        REFER --> SCHEDULE --> LEDGER
    end

    REGISTRY --> QUERY
    ARTIFACT -- "No" --> REFER
    HOLDOUT -- "No" --> REFER
    DOMAIN -- "No" --> REFER
    SCHEDULE -. "new validated cases" .-> FULL
    CLOSURE --> LEDGER
    REGISTRY --> LEDGER

    classDef calibration fill:#e8f8ee,stroke:#15803d,color:#052e16
    classDef inference fill:#e8f1ff,stroke:#2563eb,color:#172554
    classDef gate fill:#fff1f2,stroke:#be123c,color:#4c0519
    classDef record fill:#fff7dc,stroke:#b45309,color:#451a03
    class DESIGN,SPLIT,FULL,TRAIN,REGISTRY calibration
    class QUERY,CLOSURE inference
    class ARTIFACT,HOLDOUT,DOMAIN,REFER,SCHEDULE gate
    class LEDGER record
```

Full physics always takes precedence. A surrogate supplies only a calibrated
transport closure; it is never relabeled as multiphysics evidence and never
authorizes candidate exclusion.

## Scientific software responsibilities

```mermaid
flowchart TB
    INPUT["Candidate structure and operating case"]

    subgraph SCREENING["Fast screening · prioritizes what runs next"]
        direction LR
        META["Meta eSEN<br/>estimates structures, energies and forces"]
        SURROGATE["Analytical and learned models<br/>estimate rank and uncertainty"]
        META --> SURROGATE
    end

    subgraph ATOMISTIC["Atomistic validation · resolves candidate chemistry"]
        direction LR
        QE["Quantum ESPRESSO<br/>DFT relaxation, NEB and frequencies"]
        HAMILTONIAN["Candidate-specific<br/>active-space Hamiltonian"]
        CUDA["CUDA-Q and cuQuantum<br/>execute the VQE workflow"]
        QE --> HAMILTONIAN --> CUDA
    end

    subgraph REACTOR["Reactor physics · resolves chemistry and transport"]
        direction LR
        OPENFOAM["OpenFOAM<br/>flow, bubbles and multiphase mixing"]
        FENICS["FEniCSx<br/>species, heat, ions, charge and potential"]
        CANTERA["Cantera<br/>kinetics, thermochemistry and reduced reactors"]
        OPENFOAM -->|hydrodynamic closure| CANTERA
        OPENFOAM -->|NTEC handoff| FENICS
        CANTERA -->|reaction source terms| FENICS
        FENICS -->|updated transport state| CANTERA
    end

    DEVICE["PEMFC and stack models<br/>predict power and efficiency"]
    EVIDENCE["Evidence system<br/>records provenance, novelty and readiness"]

    INPUT --> META
    SURROGATE --> QE
    SURROGATE --> CANTERA
    CANTERA --> DEVICE
    QE --> EVIDENCE
    CUDA --> EVIDENCE
    OPENFOAM --> EVIDENCE
    FENICS --> EVIDENCE
    CANTERA --> EVIDENCE
    DEVICE --> EVIDENCE

    classDef input fill:#f8fafc,stroke:#475569,color:#0f172a
    classDef screening fill:#e0f2fe,stroke:#0369a1,color:#082f49
    classDef atomistic fill:#f3e8ff,stroke:#7e22ce,color:#3b0764
    classDef reactor fill:#e8f8ee,stroke:#15803d,color:#052e16
    classDef evidence fill:#fff7dc,stroke:#b45309,color:#451a03
    class INPUT input
    class META,SURROGATE screening
    class QE,HAMILTONIAN,CUDA atomistic
    class OPENFOAM,FENICS,CANTERA reactor
    class DEVICE,EVIDENCE evidence
```

Cantera owns chemistry and reduced reactor kinetics. OpenFOAM owns resolved
flow and multiphase hydrodynamics. FEniCSx owns coupled continuum transport.
Quantum ESPRESSO owns plane-wave DFT. CUDA-Q/cuQuantum execute the selected
quantum workflow. None is a universal simulator.

## Methane-conversion mode routing

```mermaid
flowchart TB
    REQUEST["Select a methane-conversion pathway"]
    MODE{"Which physics does<br/>the selected mode require?"}

    subgraph THERMAL["Thermocatalytic routes"]
        direction LR
        DEFAULT["thermocatalytic<br/>default route"]
        PFR_MODE["thermocatalytic_pfr"]
        FLUID_MODE["thermocatalytic_fluidized"]
        PFR["Plug-flow reactor<br/>fixed catalytic bed · Cantera"]
        FLUID["Fluidized bed<br/>gas-solid flow · OpenFOAM + Cantera"]
        DEFAULT --> PFR
        DEFAULT --> FLUID
        PFR_MODE --> PFR
        FLUID_MODE --> FLUID
    end

    subgraph SPECIALIZED["Specialized routes"]
        direction LR
        MMB_MODE["mmbcr"]
        NTEC_MODE["ntec"]
        ELECTRO_MODE["electrochemical"]
        MMB["Molten-metal bubble column<br/>OpenFOAM + Cantera"]
        NTEC["NTEC interface<br/>OpenFOAM + FEniCSx + Cantera"]
        ELECTRO["Electrode-electrolyte interface<br/>FEniCSx + Cantera"]
        MMB_MODE --> MMB
        NTEC_MODE --> NTEC
        ELECTRO_MODE --> ELECTRO
    end

    REQUEST --> MODE
    MODE --> DEFAULT
    MODE --> PFR_MODE
    MODE --> FLUID_MODE
    MODE --> MMB_MODE
    MODE --> NTEC_MODE
    MODE --> ELECTRO_MODE

    classDef control fill:#e8f1ff,stroke:#2563eb,color:#172554
    classDef thermal fill:#e8f8ee,stroke:#15803d,color:#052e16
    classDef specialized fill:#f3e8ff,stroke:#7e22ce,color:#3b0764
    class REQUEST,MODE control
    class DEFAULT,PFR_MODE,FLUID_MODE,PFR,FLUID thermal
    class MMB_MODE,NTEC_MODE,ELECTRO_MODE,MMB,NTEC,ELECTRO specialized
```

Mode selection chooses physics; it does not award performance. Material-bed
incompatibility produces a non-excluding, not-applicable result.

## Evidence authority ladder

```mermaid
flowchart TB
    subgraph PRIORITIZE["1 · Prioritization evidence"]
        direction LR
        ENUM["Enumerated candidate<br/>identity established"]
        SURROGATE["Surrogate estimate<br/>priority and uncertainty"]
        REDUCED["Reduced model<br/>diagnostic prediction"]
        ENUM --> SURROGATE --> REDUCED
    end

    NUMERICAL{"Are candidate-specific calculations<br/>complete and converged?"}

    subgraph VALIDATED["2 · Validated computational evidence"]
        direction LR
        ATOMISTIC["Converged atomistic result<br/>pathway or ORR evidence"]
        MULTIPHYSICS["Converged multiphysics artifact<br/>transport and coupling evidence"]
        ATOMISTIC --> MULTIPHYSICS
    end

    EMPIRICAL{"Do calibrated measurements<br/>confirm the prediction?"}

    subgraph CLAIM["3 · Defensible empirical claim"]
        direction LR
        EXPERIMENT["Calibrated experiment<br/>measured performance"]
        BLIND["Blinded holdout + prior-art check<br/>defensible discovery claim"]
        EXPERIMENT --> BLIND
    end

    REDUCED --> NUMERICAL
    NUMERICAL -- "Yes" --> ATOMISTIC
    NUMERICAL -- "No" --> SURROGATE
    MULTIPHYSICS --> EMPIRICAL
    EMPIRICAL -- "Yes" --> EXPERIMENT
    EMPIRICAL -- "No" --> SURROGATE

    classDef screening fill:#f1f5f9,stroke:#475569,color:#0f172a
    classDef validated fill:#e8f1ff,stroke:#2563eb,color:#172554
    classDef empirical fill:#e8f8ee,stroke:#15803d,color:#052e16
    classDef gate fill:#fff1f2,stroke:#be123c,color:#4c0519
    class ENUM,SURROGATE,REDUCED screening
    class ATOMISTIC,MULTIPHYSICS validated
    class EXPERIMENT,BLIND empirical
    class NUMERICAL,EMPIRICAL gate
```

Evidence advances only through completed gates. Enumeration, a favorable
surrogate, incomplete DFT, or a process exit code cannot become a discovery.

## Replaceable component architecture

```mermaid
flowchart TB
    subgraph CONTROL["Control plane"]
        direction LR
        ENTRY["CLI, notebook,<br/>scheduler or test"]
        COORDINATOR["Pipeline coordinator<br/>orders phases and handoffs"]
        RUNTIME["PipelineRuntime<br/>clock, state and mode"]
        COMPONENTS["PipelineComponents<br/>selected implementations"]
        ENTRY --> COORDINATOR
        RUNTIME -. "execution context" .-> COORDINATOR
        COMPONENTS -. "component graph" .-> COORDINATOR
    end

    subgraph BOUNDARY["Replaceable stage boundary"]
        direction LR
        STAGE["Selected stage<br/>Discovery · Reactor · DFT · VQE · Fuel cell · Report"]
        SERVICES["Stage-local service bundle<br/>replaceable scientific operations"]
        RESOURCES["Solvers, files,<br/>registries and GPUs"]
        RESOURCES --> SERVICES --> STAGE
    end

    OUTCOME["StageOutcome<br/>persistent state + downstream products"]
    CONTRACT{"Are all required products<br/>valid and complete?"}
    PERSIST["Persist state<br/>and start the next phase"]
    REJECT["Reject the stage result<br/>do not save partial state"]

    COORDINATOR --> STAGE --> OUTCOME --> CONTRACT
    CONTRACT -- "Yes" --> PERSIST --> COORDINATOR
    CONTRACT -- "No" --> REJECT
    COMPONENTS -. "chooses the stage implementation" .-> STAGE

    classDef control fill:#e8f1ff,stroke:#2563eb,color:#172554
    classDef stage fill:#e8f8ee,stroke:#15803d,color:#052e16
    classDef evidence fill:#fff7dc,stroke:#b45309,color:#451a03
    classDef gate fill:#fff1f2,stroke:#be123c,color:#4c0519
    class ENTRY,COORDINATOR,RUNTIME,COMPONENTS control
    class STAGE,SERVICES,RESOURCES stage
    class OUTCOME,PERSIST evidence
    class CONTRACT,REJECT gate
```

The replacement tests run every phase independently, make every unselected
stage raise if called, verify production signatures, and reject malformed
handoffs before persistence.

## Provenance chain

```mermaid
flowchart TB
    subgraph INPUTS["1 · Bind the requested calculation"]
        direction LR
        CONFIG["Campaign configuration<br/>seed, ranges and budgets"]
        CASE["hydrogen_case.json<br/>identity, schema and units"]
        HASH["Pristine case-tree<br/>and model hashes"]
        CONFIG --> CASE --> HASH
    end

    subgraph EXECUTION["2 · Record what actually ran"]
        direction LR
        RUN["Solver execution<br/>versions and iterations"]
        ARTIFACT["Validated artifact<br/>outputs, convergence and balances"]
        DECISION["Surrogate publication,<br/>closure or referral decision"]
        RUN --> ARTIFACT --> DECISION
    end

    subgraph PROOF["3 · Preserve an auditable conclusion"]
        direction LR
        LEDGER["SHA-256 event chain<br/>tamper-evident lineage"]
        CERTIFICATE["Coverage and readiness certificates<br/>scope of the supported claim"]
        LEDGER --> CERTIFICATE
    end

    HASH --> RUN
    CONFIG --> LEDGER
    ARTIFACT --> LEDGER
    DECISION --> LEDGER

    classDef input fill:#e8f1ff,stroke:#2563eb,color:#172554
    classDef execution fill:#f3e8ff,stroke:#7e22ce,color:#3b0764
    classDef proof fill:#fff7dc,stroke:#b45309,color:#451a03
    class CONFIG,CASE,HASH input
    class RUN,ARTIFACT,DECISION execution
    class LEDGER,CERTIFICATE proof
```

The chain answers: what was requested, what exact inputs and software ran,
which checks passed, and which later decisions consumed that evidence.

## Compute allocation and escalation

```mermaid
flowchart TB
    POOL["Pending candidates and cases"]

    subgraph POLICY["Implemented allocation policy"]
        direction LR
        COVERAGE["Reserve a fixed budget<br/>for every chemistry region"]
        PRIORITY["Rank remaining work by<br/>improvement, uncertainty and disagreement"]
        QUEUE["Build a fidelity-aware queue"]
        COVERAGE --> QUEUE
        PRIORITY --> QUEUE
    end

    subgraph FIDELITY["Execute only the fidelity the evidence requires"]
        direction LR
        FAST["Fast screening model"]
        REDUCED["Reduced reactor or<br/>ML atomistic model"]
        FULL["DFT, NEB, ORR or<br/>full multiphysics"]
    end

    RESULT["Validate the result<br/>and collect telemetry"]
    CALIBRATE["Update calibration error<br/>within its chemistry region"]
    COST["Planned extension<br/>predict CPU, GPU, memory and wall time"]

    POOL --> COVERAGE
    POOL --> PRIORITY
    QUEUE --> FAST --> RESULT
    QUEUE --> REDUCED --> RESULT
    QUEUE --> FULL --> RESULT
    RESULT --> CALIBRATE
    CALIBRATE -. "reallocate the next round" .-> PRIORITY
    CALIBRATE -. "preserve coverage" .-> COVERAGE
    COST -. "future cost-aware scheduling" .-> QUEUE

    classDef policy fill:#e8f1ff,stroke:#2563eb,color:#172554
    classDef execution fill:#e8f8ee,stroke:#15803d,color:#052e16
    classDef evidence fill:#fff7dc,stroke:#b45309,color:#451a03
    classDef planned fill:#f1f5f9,stroke:#64748b,color:#334155,stroke-dasharray: 5 5
    class POOL,COVERAGE,PRIORITY,QUEUE policy
    class FAST,REDUCED,FULL execution
    class RESULT,CALIBRATE evidence
    class COST planned
```

Coverage floors and scientific feedback are implemented. Resource-cost
prediction remains a recommended extension and is shown here as intended
architecture, not completed scientific authority.
