# Turquoise Hydrogen — Autonomous Multi-Scale Catalyst Discovery

A GPU-accelerated computational pipeline for autonomous catalyst discovery targeting **turquoise hydrogen production** (methane pyrolysis via NTEC) and **PEM fuel cell** energy conversion. Exhaustively traverses a **21.1-billion-configuration** encoded design space across 14 material classes using deterministic branch-and-bound, Meta's FAIR Chemistry equivariant graph neural networks, reactor-scale simulation, density functional theory, and variational quantum chemistry.

---

### 📖 References & Deep-Dives

* 🔬 **[Turquoise Hydrogen Reference Guide](TURQUOISE_HYDROGEN.md)**: Exhaustive literature review of thermocatalytic and nanotribo-mechano-electrochemical (NTEC) methane splitting.
* ⚡ **[Fuel Cell ORR & MEA Guide](FUEL_CELL.md)**: Comprehensive description of state-of-the-art catalysts, MEA designs, and large-scale PEMFC stack configurations.

---

## Table of Contents

- [Overview](#overview)
- [Where to Start](#where-to-start)
- [Architecture](#architecture)
- [Design Space](#design-space)
- [Simulation Software Stack](#simulation-software-stack)
- [Pipeline Modules](#pipeline-modules)
- [Physical Models](#physical-models)
- [Environment Setup](#environment-setup)
- [HuggingFace Token & Meta Models Setup](#huggingface-token--meta-models-setup)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Hardware Requirements](#hardware-requirements)
- [How It Works — Phase by Phase](#how-it-works--phase-by-phase)
- [Monitoring a Running Campaign](#monitoring-a-running-campaign)
- [Outputs & Results](#outputs--results)
- [Citation](#citation)
- [License](#license)

---

## Overview

The pipeline answers two questions end-to-end:

1. **Which catalyst best cracks methane into hydrogen + solid carbon** at high yield, low coking, and low cost?
2. **Which cathode catalyst + membrane + stack design** converts that hydrogen with the highest efficiency and least overvoltage, while maximizing output electrical power?

It does this autonomously across six phases:

```
Methane (CH₄) ──→ Phase 1-4: Catalyst Discovery ──→ H₂ + C(s)
                                                       │
                                                       ▼
                               Phase 5-6: Fuel Cell Optimization ──→ Electricity
```

## Where to Start

The repository has one production launcher and one current pilot launcher:

| Goal | Start here | Supporting code |
|------|------------|-----------------|
| Run or resume a discovery campaign | `run_production_campaign.py` | `pipeline/orchestrator.py`, `pipeline/branch_search.py` |
| Reproduce the locked divide-and-conquer pilot | `run_divide_conquer_pilot.py` | `pipeline/pilot_benchmark.py`, `pipeline/small_data_ranker.py` |
| Check scientific and implementation invariants | `test_pipeline.py`, `audit_pipeline.py` | readiness and claim gates under `pipeline/` |
| Monitor an active local run | `live_dashboard.py` | generated state under `results/` |

Inside `pipeline/`, the main code paths are grouped conceptually as follows:

- **Search:** `indexed_space.py`, `exhaustive_search.py`, `branch_search.py`,
  `discovery.py`, and `adaptive_validation.py`.
- **Candidate scoring:** `surface_screener.py`, `fc_screener.py`,
  `surrogate_model.py`, and `small_data_ranker.py`.
- **High-fidelity validation:** `qe_workflows.py`, `orr_workflows.py`,
  `dft_validator.py`, and `dft_fuel_cell.py`.
- **System models:** `reactor_models.py`, `ntec_model.py`, `pemfc_model.py`,
  and `fuel_cell_stack.py`.
- **Evidence and release gates:** `prior_art.py`, `novelty_benchmark.py`,
  `readiness.py`, `campaign_status.py`, and `report_generator.py`.

Generated outputs, downloaded model weights, pseudopotentials, mechanisms, and
Python caches are intentionally ignored. They are runtime assets, not source.

### Turquoise Hydrogen Regimes: NTEC vs. Thermocatalytic Pyrolysis

Methane splitting (pyrolysis) traditionally requires high temperatures due to the high activation barrier of the C-H bond. This pipeline supports dual-mode screening to optimize for both regimes over a shared sweep range of **500°C (773.15 K) to 1300 K** (using points `[773.15, 900.0, 1100.0, 1300.0] K`):

* **Nanotribo-Mechano-Electrochemical (NTEC) Pyrolysis (Default):**
  - **Mechanism:** Employs mechanical fluidization or shearing forces to create local triboelectric fields, facilitating C-H bond activation.
  - **Coking Resistance:** NTEC assistance is zero unless measured operating
    inputs and explicit paired NTEC/control effect measurements are supplied.
    The bounded transfer model in `pipeline/ntec_model.py` remains modeled
    evidence for a new catalyst, not candidate-specific validation.

* **Thermocatalytic Pyrolysis:**
  - **Mechanism:** Standard thermochemical activation where carbon splitting is driven purely by bulk temperature and traditional solid/alloy surface kinetics.
  - **Coking Resistance:** No mechanical coking bonuses are applied, focusing the optimization on high-temperature phase stability and traditional activation barriers.

---

## Architecture

```
Phase 1: SCREENING + OPTIMIZATION            Phase 2: REACTOR SIMULATION
┌─────────────────────────────────┐          ┌──────────────────────────┐
│  21.1B Design Space             │          │  Cantera 3.2             │
│  │                              │          │  ├─ MMBCR (bubble col.)  │
│  ▼                              │          │  ├─ PFR (plug flow)     │
│  eSen-SM (3 GPUs, 6 workers)    │──Top-K──→│  ├─ Fluidized bed        │
│  │                              │          │  └─ TST + BEP kinetics   │
│  ▼                              │          └──────────────────────────┘
│  Surrogate NN (PyTorch)         │                    │
│  │                              │                    ▼
│  ▼                              │          Phase 3: DFT VALIDATION
│  Branch-and-bound + archives    │          ┌──────────────────────────┐
└─────────────────────────────────┘          │  Quantum ESPRESSO pw.x   │
           │                                 │  ├─ Bulk SCF / vc-relax  │
           │                                 │  ├─ Slab relaxation      │
           │                                 │  └─ Adsorption energies  │
           │                                 └──────────────────────────┘
           │                                              │
           │                                              ▼
           │                                 Phase 4: QUANTUM CHEMISTRY
           │                                 ┌──────────────────────────┐
           │                                 │  CUDA-Q / cuQuantum      │
           │                                 │  VQE transition states   │
           │                                 │  ├─ C-H bond splitting   │
           │                                 │  └─ O-O bond activation  │
           │                                 └──────────────────────────┘
           │
           ▼
Phase 5: FUEL CELL                        Phase 6: REPORTING
┌─────────────────────────────────┐      ┌──────────────────────────┐
│  ORR cathode screening (eSen)   │      │  Auto-generated Markdown │
│  ├─ 137+ PGM-free candidates   │      │  + JSON pipeline state   │
│  ├─ Butler-Volmer kinetics      │──→   │  Pareto front analysis   │
│  ├─ 1D PEMFC polarization       │      │  Champion catalyst cards │
│  └─ N-cell stack + BOP + TEA    │      └──────────────────────────┘
└─────────────────────────────────┘
```

---

## Design Space

**21,092,645,031** (21.1 billion) encoded Cartesian configurations across 14 material classes:

This is the exact addressable denominator used by the indexed scanner. It is
not a claim of 21.1B symmetry-distinct physical structures: canonical IDs merge
representational duplicates, while conservative feasibility rules record invalid
Cartesian combinations as rejected rather than silently removing them.

| Class | Configs | Description |
|-------|--------:|-------------|
| **SolidCatalyst** | 20,890,448,640 | 35 active metals × 50 supports × 12 facets × 67 dopants with strain configurations |
| **HEA** | 200,344,320 | 4–6 component high-entropy alloys from 35 elements |
| **MetalHydride** | 796,068 | Alanates, borohydrides, amides, intermetallic AB₅/AB₂/AB with 13 additives |
| **Perovskite** | 449,280 | ABO₃ oxides — 16 A-site × 27 B-site with dopant fractions and defect types |
| **DAC** | 155,952 | Dual-atom metal pairs (37² combinations) × 12 coordination environments × 9 substrates |
| **MoltenMetal** | 134,640 | 13 low-melting hosts × 53 promoters × 16 concentrations × 15 temperatures |
| **MOF** | 112,047 | 35 metal nodes × 17 organic linkers × 13 cavity types × 13 pore sizes |
| **COF** | 81,120 | 36 metals × 12 covalent linkages (imine, triazine, boroxine, etc.) |
| **SAC** | 55,404 | 37 single-atom metals × 18 coordinations × 9 substrates × axial ligands |
| **MAXPhase** | 37,800 | M_{n+1}AX_n layered ternary ceramics (14 M × 13 A × 2 X × 3 n) |
| **Spinel** | 14,400 | AB₂O₄ spinel oxides (Ni, Co, Fe, Mn, Zn, Mg) × dopants × morphology × carbon supports |
| **MetalFreeCarbon** | 7,200 | Nitrogen-doped carbon (pyridinic, pyrrolic, graphitic) × defects × co-dopant (B, S, P, F) |
| **MXene** | 6,480 | M_{n+1}X_n carbides/nitrides (M elements × X elements × terminations × single-atom metals) |
| **SAA** | 1,680 | Single-atom alloys (dilute trace metals in molten metal host) × facets × loadings |

Each genome encodes into a **353-dimensional** feature vector for the surrogate neural network.

---

## Simulation Software Stack

### GPU-Accelerated Engines

| Software | Version | Role | Phase |
|----------|---------|------|-------|
| **Meta eSen-SM** | OC25 | Equivariant GNN surface-catalysis potential — slab relaxation, adsorption energies, barriers | 1, 5 |
| **PyTorch** | 2.11.0 | Multi-GPU GNN inference + surrogate NN training/prediction | 1, 5 |
| **CUDA-Q** | 0.12.0 | Variational Quantum Eigensolver on GPU quantum simulator | 4 |
| **cuQuantum** | 26.6.0 | Accelerated statevector simulation backend for CUDA-Q | 4 |
| **Cantera** | 3.2.0 | Chemical kinetics — reactor ODEs with custom YAML mechanisms | 2 |
| **Quantum ESPRESSO** | 7.x | Plane-wave DFT (pw.x) — SCF, relaxation, electronic structure | 3 |

### Scientific Libraries

| Library | Version | Role |
|---------|---------|------|
| **ASE** | 3.29.0 | Atomic structure generation (slabs, clusters, perovskites, hydrides) |
| **fairchem-core**| 2.x | Meta's FAIR Chemistry machine learning interatomic potentials framework |
| **NumPy** | 2.2.6 | Vectorized computation, objectives, feature encoding |
| **SciPy** | 1.15.2 | Electrode kinetics, Nernst equation, ODE integration |
| **Pandas** | 2.3.3 | Screening database I/O, population tracking |

### Custom Models (Pure Python/PyTorch)

| Module | Physics |
|--------|---------|
| `surrogate_model.py` | Multi-task NN predicting E_act, coking, validity (~1000× faster than eSen-SM) |
| `branch_search.py` | Persistent deterministic branch subdivision, priority, and coverage certification |
| `pemfc_model.py` | 1D through-MEA PEM fuel cell (Tafel + Ohmic + mass transport losses) |
| `fuel_cell_stack.py` | N-cell stack scaling with balance-of-plant and $/kW techno-economics |
| `reactor_mechanisms.py` | TST/BEP mechanism generator producing Cantera 3.x-compliant YAML |

---

## Pipeline Modules

| Module | Lines | Description |
|--------|------:|-------------|
| `catalyst_spaces.py` | 1116 | 14-class encoded design definitions and feature encoding |
| `surface_screener.py` | 814 | Multi-GPU parallelized Meta eSen-SM screening (slab relaxation, adsorption energies, coking index) |
| `fc_genetic_optimizer.py` | — | ORR surrogate objectives and branch-discovery orchestration |
| `genetic_optimizer.py` | — | Methane surrogate objectives and branch-discovery orchestration |
| `reactor_models.py` | 457 | Cantera reactor simulations — MMBCR, PFR, fluidized bed with custom surface kinetics |
| `dft_validator.py` | 433 | Quantum ESPRESSO input generation & parsing for champion catalysts |
| `orchestrator.py` | 396 | Core pipeline orchestrator managing phases and configurations |
| `fc_screener.py` | 374 | Meta eSen-SM-based ORR cathode screening (137+ PGM-free candidates) |
| `pemfc_model.py` | 326 | 1D PEMFC polarization model (OCV, Tafel, Ohmic, mass transport) |
| `reactor_mechanisms.py` | 302 | Cantera YAML mechanism generator with TST pre-exponentials |
| `fc_cathode_screener.py` | 297 | Generates, builds structures, and encodes cathode candidates |
| `report_generator.py` | 297 | Auto-generated Markdown + JSON pipeline reports |
| `dft_fuel_cell.py` | 278 | CHE-method ORR intermediate DFT validation |
| `vqe_transition_state.py` | 244 | CUDA-Q VQE transition states for C-H and O-O bond activation |
| `ood_detector.py` | 229 | Out-of-Distribution detector for training confidence scaling |
| `fuel_cell_stack.py` | 205 | Stack scaling, BOP parasitic loads, $/kW cost model |
| `surrogate_model.py` | 174 | Multi-task PyTorch neural network (valid/dE/E_act/coking heads) |
| `utils.py` | 511 | Constants, BEP correlations, Arrhenius rates, abundance costs, safety checks |

---

## Physical Models

### Methane Pyrolysis Descriptors

| Descriptor | Definition | Target |
|-----------|-----------|--------|
| **E_act** (activation barrier) | BEP correlation: `0.75 × ΔE_split + 0.95` eV | < 0.8 eV |
| **ΔE_H** (H* adsorption) | `E(slab+H) - E(slab) - 0.5×E(H₂)` | -0.3 to -0.5 eV |
| **ΔE_C** (C* adsorption) | `E(slab+C) - E(slab) - E(C)` | > -4.0 eV (resist coking) |
| **Coking index** | `ΔE_C - 2×ΔE_H` | Positive = resistant |
| **Segregation energy** | `E_clean - E_swapped` (dopant→surface preference) | Negative = stable |

Activation barriers at the numerical bounds (0.01 or 5.0 eV) are marked
`E_act_censored`, require DFT/NEB validation, and are excluded from champion
ranking whenever an uncensored candidate is available.

### Reactor Kinetics

| Model | Implementation |
|-------|---------------|
| Rate constants | Arrhenius: `k = A × exp(-E_act / k_B T)`, A from TST |
| Surface reactions | Cantera `ReactorSurface` with custom YAML mechanism |
| Solid carbon | Modeled as `C_graphite` gas-phase tracer species |
| Reactor types | MMBCR (molten metal bubble column), PFR, fluidized bed |

### Fuel Cell Models

| Model | Implementation |
|-------|---------------|
| ORR overpotential | Computational Hydrogen Electrode (4e⁻ pathway) |
| Cell voltage | `V = E_Nernst - η_act - η_ohm - η_mass` |
| Activation loss | Tafel: `η = (RT/αF) × ln(j/j₀)` |
| Ohmic loss | `η = j × (t_mem / σ_mem)` |
| Mass transport | `η = -c × ln(1 - j/j_L)` |
| Stack power | `P_net = n_cells × V × j × A_cell - P_BOP` |

The shared encoded space is not automatically a shared physical application
space. `MetalHydride` and `MoltenMetal` genomes are excluded from direct PEMFC
cathode ranking unless a future genome explicitly defines a solid, stable
catalyst-layer realization.

---

## Environment Setup

The pipeline requires **5 separate conda environments** due to incompatible dependency trees. Each environment serves specific phases.

### Prerequisites

- **OS**: Linux (Ubuntu 22.04+ recommended)
- **GPU**: NVIDIA GPU with CUDA 12+ (≥16 GB VRAM)
- **Conda**: Miniconda or Anaconda

### Environment 1: `fairchem-env` — Meta GNN + PyTorch (Phases 1, 5)

This environment is used for the high-throughput GNN screening and active learning loops.

```bash
conda create -n fairchem-env python=3.10 -y
conda activate fairchem-env
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121
pip install fairchem-core ase pandas scipy numpy
```

### Environment 2: `cp2k-env` — Cantera (Phase 2)

```bash
conda create -n cp2k-env python=3.12 -y
conda activate cp2k-env
conda install -c cantera cantera -y
pip install numpy scipy
```

### Environment 3: `qe-env` — Quantum ESPRESSO (Phase 3)

```bash
conda create -n qe-env python=3.10 -y
conda activate qe-env
conda install -c conda-forge qe -y
pip install numpy scipy
```

**Pseudopotentials** — download PBE RRKJUS PSL files:
```bash
mkdir -p quantum_espresso/pseudo && cd quantum_espresso/pseudo
# Download from https://www.quantum-espresso.org/pseudopotentials/
# Required elements: H, C, N, O, B, S, P, F, Na, Mg, Al, Si,
# Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, Mo, Ru, Rh, Pd,
# Ag, In, Sn, Sb, W, Pt, Au, Pb, Bi, La, Ce, Zr, Y, Nb, Te
```

### Environment 4: `quantum-env` — CUDA-Q (Phase 4)

```bash
conda create -n quantum-env python=3.10 -y
conda activate quantum-env
pip install cuda-quantum cuquantum-cu12 numpy scipy
```

### Environment 5: `battery-env` — Lightweight (Phase 6, utilities)

```bash
conda create -n battery-env python=3.10 -y
conda activate battery-env
pip install numpy scipy pandas
```

### Verify Installation

```bash
# Test all environments
conda run -n fairchem-env python -c "import torch, fairchem.core, ase; print('fairchem-env OK')"
conda run -n cp2k-env python -c "import cantera; print(f'cp2k-env OK: Cantera {cantera.__version__}')"
conda run -n qe-env bash -c "which pw.x && echo 'qe-env OK'"
conda run -n quantum-env python -c "import cudaq; print('quantum-env OK')"
conda run -n battery-env python -c "import numpy, scipy; print('battery-env OK')"
```

---

## HuggingFace Token & Meta Models Setup

The pipeline relies on Meta's FAIR Chemistry **eSen (EquiformerV2 Energy-Conserving)** model (`esen-sm-conserving-all-oc25`) for surface-catalysis energy evaluations. This model is hosted as a gated repository on HuggingFace Hub.

### Configuration Instructions

1. **Accept License Terms**: 
   Visit the model card on Hugging Face (e.g. [Meta FAIR Chemistry](https://huggingface.co/collections/facebook/fair-chemistry-671a556d11d0445a6c382218)) and request access to the gated checkpoints by agreeing to the academic use terms.

2. **Obtain HuggingFace Token**:
   Generate a **Read** access token from your HuggingFace account at: [huggingface.co/settings/tokens](https://huggingface.co/settings/tokens).

3. **Deploy Token to Workspace**:
   Create a file named `.hf_token` in the root of the project repository containing ONLY your token:
   ```bash
   echo "hf_your_token_here" > /home/ilhanraja/.gemini/antigravity/scratch/hydrogen/.hf_token
   chmod 600 /home/ilhanraja/.gemini/antigravity/scratch/hydrogen/.hf_token
   ```
   Alternatively, you can export it to your environment:
   ```bash
   export HF_TOKEN="hf_your_token_here"
   ```

---

## Usage

### Quick Test (5 min)

```bash
conda run -n fairchem-env python -m pipeline.orchestrator --quick --no-dft --no-vqe
```

Quick mode runs deterministic calibration and one resumable terminal branch leaf,
then exercises the downstream reactor, validation, and reporting path. It is a
smoke test and does not produce a `complete: true` 21.1B coverage certificate.

### Test Suite

```bash
conda run -n deepmd-env python test_pipeline.py
conda run -n deepmd-env python audit_pipeline.py
```

The active suite verifies indexed-space boundaries, disjoint shards, deterministic
tree probes across all 14 classes, branch resume, no surrogate-based pruning,
gap/overlap detection, population-denominator enforcement, coverage certificates,
blocked legacy GA entry points, and consistency between this README and the
branch-only production CLI.

### Production Campaign (48 hours)

To run a production-scale campaign, set OpenMP/MKL environment variables to prevent CPU thread over-subscription thrashing, and launch the campaign background script:

```bash
# Set OpenMP and MKL thread limits to 1 to avoid multiprocessing CPU bottlenecking
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Launch GPU-saturated campaign across all GPUs
nohup /home/ilhanraja/miniconda3/envs/fairchem-env/bin/python -u run_production_campaign.py \
  --calibration-probes 500 \
  --validation-batch 500 \
  --branch-leaf-size 1000000 \
  --prior-art-csv data/literature_registry.csv \
  --prior-art-csv data/patent_registry.csv \
  --hours 48 \
  --top-k 200 \
  > results/campaign_v6.log 2>&1 &
```

**Key parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--calibration-probes` | 500 | Deterministic binary-tree points used for initial surrogate evidence |
| `--validation-batch` | 500 | Global and regional champions sent to the atomistic model |
| `--min-validation-per-class` | 2 | Fixed validation quota reserved for every represented material class before adaptive allocation |
| `--branch-leaf-size` | 1,000,000 | Maximum indexed population per exhaustively streamed terminal leaf |
| `--branch-max-leaves` | 0 | Staged leaf limit; zero continues until the tree is complete |
| `--hours` | 0 | Campaign-wide wall-clock limit shared by both applications; zero is unlimited |
| `--top-k` | 200 | Top-K catalysts forwarded to reactor simulation |
| `--prior-art-db` | `results/prior_art.sqlite` | Persistent literature/patent/experimental identity registry |
| `--prior-art-csv` | — | Import a registry CSV (repeatable; requires a `genome` column) |
| `--final-campaign` | false | Exit nonzero unless both coverage certificates are complete and prior art is populated |
| `--evidence-manifest` | `results/evidence_manifest.json` | Counts of converged calculations and measured reactor/MEA/durability/NTEC-control evidence required in final mode |
| `--no-dft` | false | Skip Quantum ESPRESSO phase |
| `--no-vqe` | false | Skip CUDA-Q VQE phase |
| `--mode` | `ntec` | Pyrolysis mode: `ntec` (nanotriboelectric) or `thermocatalytic` |
| `--ntec-conditions-json` | — | Measured NTEC operating inputs plus paired-control effect calibration; incomplete inputs receive no numerical assistance |

### Pyrolysis Modes: NTEC vs. Thermocatalytic

The pipeline supports dual-mode screening of methane conversion mechanisms, toggled via the `--mode` flag. Both modes sweep the same temperature range from **500°C (773.15 K) to 1300 K** (`[773.15, 900.0, 1100.0, 1300.0] K`):

1. **NTEC Mode (Default):**
   * **Catalyst Physics:** Uses explicit NTEC operating inputs and paired-control
     effect measurements. Missing calibration yields zero assistance and
     `unknown` evidence status.

2. **Thermocatalytic Mode:**
   * **Catalyst Physics:** Standard thermal cracking without mechanical shear bonuses, prioritizing materials with high thermal stability and low activation energy.

### Single Phase Execution

```bash
# Run only Phase 2 (reactor simulation)
conda run -n cp2k-env python -m pipeline.orchestrator --phase 2

# Run Phases 1-3
conda run -n fairchem-env python -m pipeline.orchestrator --start 1 --end 3
```

### Standalone Module Testing

```bash
# Test design space
conda run -n battery-env python -m pipeline.catalyst_spaces

# Test eSen screening with deterministic tree probes
conda run -n fairchem-env python -m pipeline.surface_screener

# Test DFT input generation (no pw.x execution)
conda run -n battery-env python -m pipeline.dft_validator

# Test PEMFC model
conda run -n battery-env python -m pipeline.pemfc_model
```

---

## Project Structure

```
hydrogen/
├── README.md                      # This file
├── LICENSE                        # MIT license
├── requirements.txt               # Python dependencies
├── environment.yml                # Conda environment spec
├── run_production_campaign.py     # Production launcher (GPU-saturated)
├── .hf_token                      # HuggingFace token (chmod 600, gitignored)
│
├── pipeline/                      # Core pipeline package
│   ├── __init__.py
│   ├── catalyst_spaces.py         # 21.1B encoded design-space definitions
│   ├── indexed_space.py           # O(1) global candidate addressing + shards
│   ├── exhaustive_search.py       # Resumable bounded-memory population scan
│   ├── branch_search.py           # Persistent divide-and-conquer + certificate
│   ├── discovery.py               # Canonical IDs + novelty/coverage acquisition
│   ├── surface_screener.py        # Multi-GPU Meta eSen screening (methane pyrolysis)
│   ├── surrogate_model.py         # Multi-task PyTorch surrogate NN
│   ├── genetic_optimizer.py       # Methane objectives + branch orchestration
│   ├── reactor_mechanisms.py      # Cantera YAML mechanism generator
│   ├── reactor_models.py          # MMBCR / PFR / fluidized bed
│   ├── dft_validator.py           # Quantum ESPRESSO DFT validation
│   ├── dft_fuel_cell.py           # ORR intermediate DFT (CHE)
│   ├── vqe_transition_state.py    # CUDA-Q VQE transition states
│   ├── fc_screener.py             # Meta eSen screening (ORR fuel cell)
│   ├── pemfc_model.py             # 1D PEMFC polarization model
│   ├── fuel_cell_stack.py         # Stack scaling + TEA
│   ├── report_generator.py        # Auto-report generation
│   └── utils.py                   # Constants, helpers, I/O
│
├── quantum_espresso/              # QE pseudopotentials (gitignored)
│   └── pseudo/                    # .UPF files (download separately)
│
├── mechanisms/                    # Generated Cantera YAML (gitignored)
└── results/                       # Pipeline outputs (gitignored)
    ├── screening/                 # Branch database, certificates, GNN CSVs
    ├── reactor/                   # Cantera simulation results
    ├── dft/                       # QE input/output files
    ├── vqe/                       # VQE energetics (JSON)
    ├── fuel_cell/                 # Cathode screening + PEMFC curves
    └── reports/                   # Auto-generated pipeline report
```

---

## Hardware Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **GPU** | 1× NVIDIA GPU, 16 GB VRAM | 3× GPUs (e.g., RTX 5090 + 2× Blackwell) |
| **RAM** | 32 GB | 64+ GB |
| **Storage** | 50 GB | 200+ GB (for full DFT campaigns) |
| **CPU** | 8 cores | 32+ cores (for Cantera multi-process) |

---

## How It Works — Phase by Phase

### Methodology: coverage-guided scientific fuzzing

The discovery engine can be understood as **coverage-guided, physics-constrained
fuzzing for catalysts**. A candidate genome is the input; the turquoise-hydrogen
reactor or PEM fuel-cell pathway is the target; and progressively more expensive
models are feedback oracles. Unlike a software fuzzer, success is not a crash.
The objective is a physically admissible, synthesizable candidate with unusually
strong hydrogen-production or fuel-cell performance.

The feedback loop is deterministic and auditable:

1. Map every encoded catalyst to a stable global index and canonical identity.
2. Divide each material-class range recursively and probe it with a deterministic
   low-discrepancy sequence.
3. Reject only candidates that fail explicit hard physical constraints. Model
   uncertainty, unfamiliar chemistry, or a poor surrogate score never proves
   that a branch is empty.
4. Rank admissible branches with robust quantiles, predicted performance,
   uncertainty, and measured surrogate/high-fidelity disagreement.
5. Preserve global, per-objective, per-class, and chemistry-region champions so
   one familiar family cannot erase unusual candidates.
6. Send selected candidates through increasingly expensive eSen/OC25, DFT/NEB,
   reactor, ORR, MEA, and experimental checks.
7. Feed paired predictions and observations back into regional and class-level
   calibration, then resume from the persistent search tree.

Discovery and calibration use separate ledgers. The discovery slate exploits the
best reproducible ranking. A fixed class-coverage slate and additional
uncertainty/disagreement selections deliberately investigate poorly calibrated
regions. Coverage calculations are not relabeled as discovery hits, and failures
remain useful feedback.

Random sampling and transparent chemistry heuristics are evaluation controls,
not production search strategies. A locked pilot fixes its candidate pool,
training digest, selection slate, budget, random seed, and number of random trials
before hidden outcomes are calculated. Invalid selected candidates consume budget
and count as misses.

Prospective computational pilots currently support approximately **1.8×
enrichment over equal-budget random sampling** for the validated application-
specific rankers:

- Turquoise hydrogen: 9/20 top-quintile hits versus 5.00/20 random mean
  (1.80×; above the 95% random bound).
- Fuel-cell ORR: 7/18 hits versus 3.75/18 random mean
  (1.87×; above the 95% random bound).

These successes occurred in separate locked rounds. They demonstrate
computational enrichment, not a universal 1.8× guarantee or experimental catalyst
performance. Other rounds exposed model drift, establishing a release criterion:
a challenger should not replace an incumbent without new locked validation.
Persistent automatic model promotion is not yet treated as completed scientific
infrastructure. The reproducible pilot entry point is
`run_divide_conquer_pilot.py`; its locked result manifests and raw selections
are versioned under `results/pilot/`, `results/screening/pilot/`, and
`results/fuel_cell/pilot/`. Earlier exploratory launchers were removed after
their useful logic was incorporated into this runner and
`pipeline/pilot_benchmark.py`.

### Phase 1: Deterministic Branch-and-Bound Discovery

1. **Calibrate at deterministic tree probes** — recursively bisected probe points from all 14 class roots establish initial model evidence; these probes do not count as population coverage
2. **eSen-SM evaluation** — for each candidate, build an atomic slab or cluster, enforce periodic boundary conditions (`pbc=True`), relax with BFGS, compute H*/CH₃*/C* adsorption energies
3. **Train deterministic small-data rankers** — turquoise hydrogen ranks the
   continuous activation barrier; ORR predicts OH/O/OOH adsorption energies and
   derives an unclipped CHE overpotential so saturated labels cannot erase order
4. **Divide all class ranges recursively** — deterministic surrogate probes prioritize child branches but never authorize pruning
5. **Resolve every terminal branch** — a branch is either exhaustively streamed or hard-pruned only after every member fails conservative feasibility rules
6. **Retain global and regional champions** — unfamiliar chemistry regions remain represented even when familiar chemistry dominates the global scores
7. **Retain every objective's winners** — bounded global archives and per-region champions prevent a primary-objective ranking from discarding selectivity, stability, cost, or uncertainty extremes
8. **Output a coverage certificate** — exact terminal population, gap/overlap checks, scan cursors, pruning proofs, canonical candidate IDs, and application-specific champions

Expensive validation is allocated adaptively by `pipeline/adaptive_validation.py`.
Each represented material class receives a fixed quota first. Remaining slots
combine expected improvement, ensemble uncertainty, regional calibration error,
and observed productivity. Paired surrogate/Fairchem/DFT/experimental results are
stored by chemistry region in SQLite. Large disagreement moves a region earlier;
repeated low productivity moves it later but never prunes it. A separate
`experimental_slate` table preserves one champion per region before repeats.

#### What “novel” means

The discovery engine uses **campaign novelty**: a candidate has not previously been
evaluated under its canonical ID, or represents a chemistry region not yet covered
by the campaign. This maximizes the chance of finding unfamiliar viable chemistry
without pretending that model uncertainty is poor performance. The repository now
includes a versioned SQLite literature/patent/experimental registry keyed by the
same canonical candidate identity and reports `known`, `region_known`, or `unseen`.
Populate it with repeatable `--prior-art-csv` inputs before making external novelty
claims. `unseen` means absent from the supplied registry, not proof of worldwide
novelty; registry completeness and chemical identity resolution still require
curated external data.

External novelty claims additionally require a prospective/time-split recovery
benchmark (`pipeline.novelty_benchmark.time_split_recovery`): candidates reported
after the training cutoff are hidden, ranked blindly, and scored by exact and
chemistry-region recall at K. Missing publication year, source ID, or citation
invalidates the benchmark.

#### Six-point scientific status

`pipeline.campaign_status.assess_campaign()` reports one fail-closed status for
complete search, validated champions, calibrated NTEC, validated reactor,
validated PEMFC, and defensible novelty. A production run is not scientifically
ready until all six are true.

#### Coverage and exhaustiveness

The finite genome space is traversed through persistent binary subdivision and
indexed streaming. Random, stratified, genetic, and rotating-grid candidate
sampling are not production search strategies. Expensive atomistic calculations remain multi-fidelity:
the repository does not claim that all billions of candidates received DFT or
experimental validation. Coverage must be reported separately for generated,
surrogate-scored, GNN-validated, and DFT-validated candidates.

Industrial gates are fail-closed: missing measurements produce `unknown`, never a
pass or a pruning decision. The configurable defaults are 700–1300 K, at least
95% H2 selectivity and 70% methane conversion, at most 1%/h deactivation and 5%
coke for turquoise hydrogen; and at most 0.40 V ORR overpotential, at least
1.00 W/cm2 peak power and 40% system efficiency, and at most 10 uV/h voltage
degradation for fuel cells. These are campaign screening criteria, not universal
industrial standards.

#### Production divide-and-conquer search

Use deterministic hierarchical branch-and-bound to process the most promising,
uncertain, novel, and populous regions first while retaining exhaustive coverage:

```bash
python run_production_campaign.py \
  --branch-leaf-size 1000000 \
  --branch-probes 9 \
  --branch-class-floor 1 \
  --branch-exploration-interval 4
```

The tree begins with one root for each of the 14 material classes and recursively
bisects class-local indexed ranges. Deterministic low-discrepancy surrogate probes
establish a robust quantile priority for each child. A class floor gives every
chemistry family a finite-budget opportunity; every fourth resolved leaf then
returns to the least-covered family before exploitation resumes. On restart,
stale pending priorities are recalculated using accumulated regional or
class-level validation disagreement. Probe predictions **never authorize pruning**.
A branch is removed only when every encoded member has been checked against the
conservative hard-feasibility rules and all fail. Every other leaf is passed to
the exhaustive streaming scanner.

For staged campaigns, limit the number of leaves handled in one invocation:

```bash
python run_production_campaign.py --branch-max-leaves 100
```

Running the same command again resumes the persistent tree and each partially
processed leaf. SQLite records pending, expanded, hard-pruned, and fully scanned
nodes, including the unresolved encoded population. This behaves like binary
divide-and-conquer without making the unsafe monotonicity assumption required by
literal binary search.

For the final fail-closed check, run the same command with `--final-campaign`.
It writes `results/campaign_readiness.json` and exits nonzero if either application
lacks a complete, denominator-matched coverage certificate or the prior-art
registry is empty. A complete computational certificate still does not substitute
for reactor, stack, durability, synthesis, safety, or experimental validation.
Copy `evidence_manifest.example.json` to `results/evidence_manifest.json` and
update its counts only from traceable records. Final mode remains nonzero until
all six scientific criteria pass. Current four-qubit VQE Hamiltonians are labeled
toy models and are not accepted as catalyst evidence.

Each invocation also regenerates an application-specific coverage certificate:

- `results/screening/turquoise_hydrogen_coverage_certificate.json`
- `results/fuel_cell/coverage_certificate.json`

The certificate verifies that terminal intervals form a gap-free, non-overlapping
partition of all 14 indexed class ranges; every `scanned` leaf has a completed
resume cursor; and every `pruned` leaf has a rechecked all-members-fail hard-rule
proof. `complete` becomes true only when the terminal population equals exactly
**21,092,645,031** and no unresolved leaf remains. The production command defaults
to `--expected-space-size 21092645031` and stops on any denominator mismatch. In
particular, labeling the present repository population as 25.3B now produces an
error rather than a false exhaustive-coverage claim.

#### Structure Generation

The eSen screener builds physically realistic, periodic atomic structures for all 14 classes:

| Class | Structure Type | Method |
|-------|---------------|--------|
| SolidCatalyst | FCC/BCC/HCP slab (3×3×4) | ASE slab builders with explicit lattice constants for 60+ elements |
| MoltenMetal | FCC slab with promoter substitutions | Host slab + random dopant placement |
| SAC / DAC | Metal-porphyrin cluster (periodic) | Square/hexagonal N/S/O/P coordination |
| MOF / COF | Metal-cavity cluster (periodic) | Porphyrin-like with organic skeleton |
| Perovskite | ABO₃ 2×2×3 supercell | Simple cubic with A/B-site doping + O vacancies |
| MetalHydride | FCC slab + interstitial H | Metal surface with H at tetrahedral sites |
| MAXPhase | HCP slab with A-element substitution | M-layer slab with interstitial dopants |
| HEA | Random-substitution FCC slab | Host + 3-5 equimolar dopants |
| Spinel | AB₂O₄ spinel slab (2×2×1) | Spinel lattice creation with A/B-site placement |
| MXene | HCP slab with terminations | M-element slab layer + OH/O/F termination |
| SAA | Host slab with isolated trace metal | Single trace atom substituted at surface |
| MetalFreeCarbon| N-doped carbon structures | Graphene/CNT base with vacancy/nitrogen doping |

### Phase 2: Cantera Reactor Simulation

For each top-K catalyst from Phase 1:
1. Generate a Cantera YAML mechanism with TST-derived rate constants calibrated to the catalyst's E_act
2. Simulate three reactor types (MMBCR, PFR, fluidized bed) across the standardized 4 temperatures (500°C to 1300 K / 773.15–1300 K)
3. Record CH₄ conversion, H₂ selectivity, carbon yield, residence time

### Phase 3: DFT Validation (Quantum ESPRESSO)

For the top 10 champion catalysts:
1. Generate QE input files (SCF / relax) with proper pseudopotentials
2. Run `pw.x` for bulk optimization and slab relaxation
3. Parse converged total energies, forces, electronic structure

### Phase 4: VQE Quantum Chemistry (CUDA-Q)

For the top 3–5 champions:
1. Build molecular Hamiltonians for C-H and O-O bond-breaking transition states
2. Run VQE with hardware-efficient ansätze on the NVIDIA GPU quantum simulator
3. Extract refined activation barriers beyond DFT accuracy

### Phase 5: Fuel Cell Modeling

1. **Cathode screening** — evaluate 137+ PGM-free ORR catalysts using Meta's eSen-SM surface model, optimizing for highest efficiency and least overvoltage first and foremost, while maximizing output electrical power.
2. **PEMFC polarization** — 1D model: Nernst OCV → Tafel activation → Ohmic → mass transport
3. **Membrane sweep** — test Nafion 211/212, Gore-Select, Aquivion across operating conditions
4. **Stack scaling** — 300–400 cell stack with balance-of-plant, gravimetric/volumetric power density, $/kW, optimized for system efficiency

### Phase 6: Report Generation

Auto-generates a comprehensive Markdown report with:
- Champion catalyst cards (genome, E_act, coking index, cost)
- Reactor performance tables
- PEMFC polarization data
- Stack-level techno-economics

---

## Monitoring a Running Campaign

```bash
# GPU utilization
nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv -l 5

# Branch discovery progress
tail -f results/screening/genetic_optimizer.log

# eSen screening throughput
tail -f results/screening/surface_screening.log

# Campaign stdout
tail -f results/campaign_v6.log

# Check pipeline state
python -c "import json; print(json.dumps(json.load(open('pipeline_state.json')), indent=2))"
```

---

## Outputs & Results

After a campaign completes, key outputs include:

| File | Contents |
|------|----------|
| `results/screening/ga_full_database.csv` | Complete screening database (all eSen-evaluated candidates) |
| `results/screening/ga_surface_gen*.csv` | Per-round GNN validation results |
| `results/reports/pipeline_report.md` | Auto-generated comprehensive report |
| `pipeline_state.json` | Machine-readable pipeline state with timing and metrics |
| `results/reactor/*.json` | Cantera simulation results per catalyst |
| `results/dft/*/` | QE input/output files per catalyst |
| `results/fuel_cell/` | Cathode screening + PEMFC polarization data |

---

## Citation

```bibtex
@software{turquoise_h2_pipeline,
  title  = {Turquoise Hydrogen: Autonomous Multi-Scale Catalyst Discovery Pipeline},
  year   = {2026},
  url    = {https://github.com/YungRaj/hydrogen},
  note   = {21.1B design space, 14 material classes, 6-phase pipeline}
}
```

## License

MIT — see [LICENSE](LICENSE).
