# Turquoise Hydrogen Pipeline

**Autonomous multi-scale catalyst discovery and system optimization for turquoise hydrogen production via methane pyrolysis (NTEC) and PEM fuel cell energy conversion.**

This pipeline integrates GPU-accelerated atomistic screening, reactor-scale microkinetic simulation, high-fidelity DFT validation, quantum chemistry, and techno-economic analysis into a closed-loop materials discovery workflow.

---

## Architecture

```
Phase 1                Phase 2              Phase 3           Phase 4
MACE-MP-0 Screening ─→ Cantera Reactor ──→ Quantum ESPRESSO ─→ CUDA-Q VQE
+ NSGA-II Genetic       MMBCR / PFR /       Bulk SCF / Slab     Transition State
  Algorithm             Fluidized Bed        Relaxation          Refinement
       │                                                              │
       │                                                              │
       ▼                                                              ▼
Phase 5                                                       Phase 6
PEMFC Cathode Screening ─→ Single-Cell ─→ Stack Scaling ─→ Report Generation
ORR Overpotential           Polarization     TEA               Markdown + JSON
```

## Modules

| Module | Description |
|--------|-------------|
| `pipeline/catalyst_spaces.py` | 107M-config chemical design space (alloys, SACs, DACs, MOFs, COFs) |
| `pipeline/mace_screener.py` | Multi-GPU MACE-MP-0 adsorption energy & barrier screening |
| `pipeline/surrogate_model.py` | Multi-task PyTorch surrogate for GA acceleration |
| `pipeline/genetic_optimizer.py` | NSGA-II multi-objective optimization (E_act, coking, cost) |
| `pipeline/reactor_mechanisms.py` | Cantera 3.x YAML mechanism generator (TST + BEP) |
| `pipeline/reactor_models.py` | MMBCR, PFR, and fluidized bed reactor simulators |
| `pipeline/dft_validator.py` | Quantum ESPRESSO structural relaxation & SCF |
| `pipeline/dft_fuel_cell.py` | ORR intermediate DFT (CHE method) |
| `pipeline/vqe_transition_state.py` | CUDA-Q VQE for C-H/O-O transition states |
| `pipeline/fuel_cell_cathode_screener.py` | MACE-based ORR cathode screening (137 candidates) |
| `pipeline/pemfc_model.py` | 1D through-MEA PEMFC polarization model |
| `pipeline/fuel_cell_stack.py` | N-cell stack scaling with BOP & techno-economics |
| `pipeline/orchestrator.py` | Master controller for all 6 phases |
| `pipeline/report_generator.py` | Auto-generated Markdown + JSON reports |

## Quick Start

### 1. Environment Setup

```bash
# Primary environment (MACE + PyTorch + ASE)
conda env create -f environment.yml
conda activate hydrogen-pipeline

# Additional environments for specific phases:
# Phase 2: conda install -c cantera cantera
# Phase 3: conda install -c conda-forge qe
# Phase 4: pip install cuda-quantum
```

### 2. Download Pseudopotentials (Phase 3 only)

```bash
mkdir -p quantum_espresso/pseudo
cd quantum_espresso/pseudo
# Download from https://www.quantum-espresso.org/pseudopotentials/
# Required: PBE RRKJUS PSL pseudopotentials for elements in PSEUDO_MAP
```

### 3. Run the Pipeline

```bash
# Quick test (50 generations, reduced population)
python -m pipeline.orchestrator --quick --no-dft --no-vqe

# Full production run
python -m pipeline.orchestrator

# Single phase
python -m pipeline.orchestrator --phase 2

# Phase range
python -m pipeline.orchestrator --start 1 --end 3
```

### 4. Results

Output is written to `results/`:
- `results/screening/` — GA population databases (CSV)
- `results/reactor/` — Cantera simulation results
- `results/dft/` — QE input/output files
- `results/vqe/` — VQE energetics (JSON)
- `results/fuel_cell/` — Cathode screening + PEMFC curves
- `results/reports/` — Auto-generated pipeline report

## Hardware Requirements

- **GPU**: NVIDIA GPU with ≥16 GB VRAM (Blackwell recommended)
- **RAM**: ≥32 GB system memory
- **Storage**: ≥50 GB for full screening campaigns

## Key Physical Models

| Model | Implementation |
|-------|---------------|
| **C-H activation barrier** | BEP correlation calibrated on TM surfaces |
| **Surface kinetics** | Transition State Theory (TST) pre-exponentials |
| **ORR overpotential** | Computational Hydrogen Electrode (4e⁻ pathway) |
| **Electrode kinetics** | Butler-Volmer equation |
| **Reactor mass transport** | Cantera 3.x IdealGasReactor + ReactorSurface |
| **PEMFC losses** | Nernst OCV + Tafel + Ohmic + mass transport |

## Citation

If you use this code in your research, please cite:

```bibtex
@software{turquoise_h2_pipeline,
  title = {Turquoise Hydrogen Pipeline: Autonomous Multi-Scale Catalyst Discovery},
  year = {2026},
  url = {https://github.com/YOUR_USERNAME/hydrogen}
}
```

## License

MIT — see [LICENSE](LICENSE).
