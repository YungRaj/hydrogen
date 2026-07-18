# Turquoise Hydrogen — Autonomous Multi-Scale Catalyst Discovery

A GPU-accelerated computational pipeline for autonomous catalyst discovery targeting **turquoise hydrogen production** (methane pyrolysis via NTEC) and **PEM fuel cell** energy conversion. Explores a **25.3-billion-configuration** design space across 10 material classes using Meta's FAIR Chemistry equivariant graph neural networks, genetic optimization, reactor-scale simulation, density functional theory, and variational quantum chemistry.

---

## Table of Contents

- [Overview](#overview)
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

---

## Architecture

```
Phase 1: SCREENING + OPTIMIZATION            Phase 2: REACTOR SIMULATION
┌─────────────────────────────────┐          ┌──────────────────────────┐
│  25.3B Design Space             │          │  Cantera 3.2             │
│  │                              │          │  ├─ MMBCR (bubble col.)  │
│  ▼                              │          │  ├─ PFR (plug flow)     │
│  eSen-SM (3 GPUs, 6 workers)    │──Top-K──→│  ├─ Fluidized bed        │
│  │                              │          │  └─ TST + BEP kinetics   │
│  ▼                              │          └──────────────────────────┘
│  Surrogate NN (PyTorch)         │                    │
│  │                              │                    ▼
│  ▼                              │          Phase 3: DFT VALIDATION
│  NSGA-II (4-objective Pareto)   │          ┌──────────────────────────┐
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

**25,296,914,418** (25.3 billion) unique catalyst configurations across 10 material classes:

| Class | Configs | Description |
|-------|--------:|-------------|
| **SolidCatalyst** | 25,055,084,160 | 35 metals × 50 supports × 12 facets × 67 dopants (oxides, carbides, nitrides, zeolites, sulfides) |
| **HEA** | 240,018,240 | 4–6 component high-entropy alloys from 35 elements (C(35,4)+C(35,5)+C(35,6) combos) |
| **MetalHydride** | 796,068 | Alanates, borohydrides, amides, intermetallic AB₅/AB₂/AB with 13 additives |
| **Perovskite** | 449,280 | ABO₃ oxides — 16 A-site × 27 B-site with dopant fractions and defect types |
| **MoltenMetal** | 168,480 | 13 low-melting hosts × 53 promoters × 16 concentrations × 15 temperatures |
| **DAC** | 164,268 | 37² dual-atom metal pairs × 12 coordination environments × 9 substrates |
| **MOF** | 103,428 | 35 metal nodes × 17 organic linkers × 13 cavity types × 13 pore sizes |
| **COF** | 75,036 | 36 metals × 12 covalent linkages (imine, triazine, boroxine, etc.) |
| **MAXPhase** | 49,140 | M_{n+1}AX_n layered ternary carbides/nitrides (14 M × 13 A × 2 X × 3 n) |
| **SAC** | 6,318 | 37 single-atom metals × 18 N/S/O/P coordinations × 9 substrates |

Each genome encodes into a **324-dimensional** feature vector for the surrogate neural network.

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
| **NumPy** | 2.2.6 | Vectorized computation, NSGA-II, feature encoding |
| **SciPy** | 1.15.2 | Electrode kinetics, Nernst equation, ODE integration |
| **Pandas** | 2.3.3 | Screening database I/O, population tracking |

### Custom Models (Pure Python/PyTorch)

| Module | Physics |
|--------|---------|
| `surrogate_model.py` | Multi-task NN predicting E_act, coking, validity (~1000× faster than eSen-SM) |
| `genetic_optimizer.py` | NSGA-II with crowding distance — 4-objective Pareto optimization |
| `pemfc_model.py` | 1D through-MEA PEM fuel cell (Tafel + Ohmic + mass transport losses) |
| `fuel_cell_stack.py` | N-cell stack scaling with balance-of-plant and $/kW techno-economics |
| `reactor_mechanisms.py` | TST/BEP mechanism generator producing Cantera 3.x-compliant YAML |

---

## Pipeline Modules

| Module | Lines | Description |
|--------|------:|-------------|
| `catalyst_spaces.py` | 790 | 10-class design space definitions, genome generators, crossover/mutation, feature encoding |
| `surface_screener.py` | 662 | Multi-GPU parallelized Meta eSen-SM screening — slab generation, adsorption energies, coking index |
| `reactor_models.py` | 452 | Cantera reactor simulations — MMBCR, PFR, fluidized bed with custom surface kinetics |
| `genetic_optimizer.py` | 416 | NSGA-II GA — surrogate-accelerated 4-objective optimization |
| `dft_validator.py` | 416 | Quantum ESPRESSO input generation & parsing for champion catalysts |
| `utils.py` | 385 | Shared constants, BEP correlations, Arrhenius rates, I/O utilities |
| `reactor_mechanisms.py` | 302 | Cantera YAML mechanism generator with TST pre-exponentials |
| `report_generator.py` | 287 | Auto-generated Markdown + JSON pipeline reports |
| `fc_screener.py` | 285 | Meta eSen-SM-based ORR cathode screening (137+ PGM-free candidates) |
| `pemfc_model.py` | 275 | 1D PEMFC polarization curves — OCV, Tafel, Ohmic, mass transport |
| `vqe_transition_state.py` | 244 | CUDA-Q VQE for C-H and O-O transition state refinement |
| `dft_fuel_cell.py` | 232 | CHE-method ORR intermediate DFT validation |
| `fuel_cell_stack.py` | 205 | Stack scaling, BOP parasitic loads, $/kW cost model |
| `surrogate_model.py` | 174 | Multi-task PyTorch neural network (valid/dE/E_act/coking heads) |

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

This runs 50 generations with 100 population, validating the full pipeline end-to-end.

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
  --pop 1000 \
  --gens 3000 \
  --mace-batch 500 \
  --mace-per-round 500 \
  --mace-interval 5 \
  --hours 48 \
  --top-k 200 \
  > results/campaign_v6.log 2>&1 &
```

**Key parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pop` | 1000 | GA population size per generation |
| `--gens` | 3000 | Total GA generations |
| `--mace-batch` | 500 | Initial GNN screening batch size |
| `--mace-per-round` | 500 | GNN evaluations per validation round |
| `--mace-interval` | 5 | Generations between GNN validation rounds |
| `--hours` | 48 | Maximum wall-clock time |
| `--top-k` | 200 | Top-K catalysts forwarded to reactor simulation |
| `--no-dft` | false | Skip Quantum ESPRESSO phase |
| `--no-vqe` | false | Skip CUDA-Q VQE phase |

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

# Test eSen screening (20 random candidates)
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
│   ├── catalyst_spaces.py         # 25.3B design space definitions
│   ├── surface_screener.py        # Multi-GPU Meta eSen screening (methane pyrolysis)
│   ├── surrogate_model.py         # Multi-task PyTorch surrogate NN
│   ├── genetic_optimizer.py       # NSGA-II 4-objective GA
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
    ├── screening/                 # GA databases, GNN CSVs
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

### Phase 1: Meta eSen Screening + Genetic Optimization

1. **Generate initial population** — random genomes from all 10 material classes
2. **eSen-SM evaluation** — for each candidate, build an atomic slab or cluster, enforce periodic boundary conditions (`pbc=True`), relax with BFGS, compute H*/CH₃*/C* adsorption energies
3. **Train surrogate NN** — multi-task network learns to predict E_act, coking index, validity from the 324-dim genome encoding (~1000× faster than eSen-SM)
4. **NSGA-II loop** — evolve population via tournament selection, uniform crossover, class-aware mutation; evaluate with surrogate; periodically validate top candidates with full eSen-SM on GPU
5. **Class-diversity enforcement** — each generation guarantees ≥5% population from every material class, preventing any single class from dominating the front
6. **Output** — Pareto-optimal front of catalysts minimizing (E_act, -coking, segregation, cost)

#### Structure Generation

The eSen screener builds physically realistic, periodic atomic structures for all 10 classes:

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

### Phase 2: Cantera Reactor Simulation

For each top-K catalyst from Phase 1:
1. Generate a Cantera YAML mechanism with TST-derived rate constants calibrated to the catalyst's E_act
2. Simulate three reactor types (MMBCR, PFR, fluidized bed) across 5 temperatures (800–1200 K)
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

# GA evolution progress
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
  note   = {25.3B design space, 10 material classes, 6-phase pipeline}
}
```

## License

MIT — see [LICENSE](LICENSE).
