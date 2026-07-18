#!/usr/bin/env python3
"""
Master Pipeline Orchestrator.

Coordinates all 6 phases of the turquoise hydrogen → fuel cell pipeline:

  Phase 1: MACE Screening → Genetic Optimization → Pareto Front
  Phase 2: Cantera Reactor Simulation (top catalysts × 3 reactor types)
  Phase 3: DFT Validation (top catalysts)
  Phase 4: VQE Transition State (top 3 catalysts)
  Phase 5: Fuel Cell Cathode Screening → PEMFC Modeling
  Phase 6: Report Generation

Usage:
    python -m pipeline.orchestrator [--phase N] [--quick]
"""

import os
import sys
import time
import json
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.utils import (
    BASE_DIR, RESULTS_DIR, SCREENING_DIR, REACTOR_DIR, DFT_DIR, VQE_DIR,
    FUEL_CELL_DIR, REPORTS_DIR, MECHANISMS_DIR,
    ENV_MACE, ENV_BATTERY, ENV_CANTERA, ENV_QE, ENV_QUANTUM,
    setup_logger, print_banner, save_json, load_json,
    run_in_env,
)

logger = setup_logger('orchestrator', 'pipeline_orchestrator.log')


@dataclass
class PipelineConfig:
    """Pipeline-level configuration."""
    # Phase 1: Screening
    initial_mace_samples: int = 500       # Initial MACE evaluations
    ga_pop_size: int = 500                # GA population size
    ga_generations: int = 200             # GA generations
    top_k_reactor: int = 50              # Top K catalysts → reactor simulation
    top_k_dft: int = 10                  # Top K → DFT validation
    top_k_vqe: int = 3                   # Top K → VQE

    # Phase 2: Reactor
    reactor_temperatures: tuple = (800, 900, 1000, 1100, 1200)
    reactor_types: tuple = ('MMBCR', 'PFR', 'Fluidized')

    # Phase 5: Fuel Cell
    fc_top_k_pemfc: int = 20            # Top cathode catalysts → PEMFC model
    fc_stack_cells: int = 300            # Stack cells

    # Runtime
    run_dft: bool = True                 # Actually execute pw.x
    run_vqe: bool = True                 # Actually execute CUDA-Q
    quick_mode: bool = False             # Reduced parameters for testing
    seed: int = 42


def run_pipeline(config: PipelineConfig = PipelineConfig(),
                 start_phase: int = 1, end_phase: int = 6):
    """
    Execute the full multi-scale simulation pipeline.
    """
    t_total = time.time()
    print_banner("TURQUOISE HYDROGEN → FUEL CELL: MULTI-SCALE PIPELINE")
    logger.info(f"Starting pipeline: phases {start_phase}–{end_phase}")
    logger.info(f"Configuration: quick_mode={config.quick_mode}")

    if config.quick_mode:
        config.initial_mace_samples = 50
        config.ga_pop_size = 100
        config.ga_generations = 50
        config.top_k_reactor = 10
        config.top_k_dft = 3
        config.top_k_vqe = 1
        config.fc_top_k_pemfc = 5

    pipeline_state = load_json("pipeline_state.json") or {}

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 1: CATALYST SCREENING & GENETIC OPTIMIZATION
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 1 <= end_phase:
        print_banner("PHASE 1: META ESEN SCREENING & GENETIC OPTIMIZATION")
        t1 = time.time()

        from pipeline.catalyst_spaces import generate_population, estimate_design_space_size
        from pipeline.surface_screener import run_screening
        from pipeline.genetic_optimizer import run_genetic_algorithm, GAConfig

        # Report design space
        sizes = estimate_design_space_size()
        logger.info(f"Design space: {sizes['TOTAL']:,} total configurations")
        for cls, size in sizes.items():
            if cls != 'TOTAL':
                logger.info(f"  {cls}: {size:,}")

        # Run GA
        ga_config = GAConfig(
            pop_size=config.ga_pop_size,
            n_generations=config.ga_generations,
            initial_mace_samples=config.initial_mace_samples,
            mace_eval_interval=max(1, config.ga_generations // 4),
            surrogate_retrain_interval=max(1, config.ga_generations // 4),
            seed=config.seed,
        )
        pareto_genomes, screening_db = run_genetic_algorithm(ga_config)

        # Select top-K from Pareto front
        valid_db = screening_db[screening_db['valid'] == True].copy()
        if 'E_act' in valid_db.columns:
            top_catalysts = valid_db.nsmallest(config.top_k_reactor, 'E_act')
        else:
            top_catalysts = valid_db.head(config.top_k_reactor)

        pipeline_state['phase1'] = {
            'pareto_size': len(pareto_genomes),
            'total_evaluated': len(screening_db),
            'valid_count': len(valid_db),
            'top_catalysts_count': len(top_catalysts),
            'elapsed_s': time.time() - t1,
        }
        if len(valid_db) > 0 and 'E_act' in valid_db.columns:
            pipeline_state['phase1']['best_E_act'] = float(valid_db['E_act'].min())
            pipeline_state['phase1']['best_coking'] = float(valid_db['coking_index'].max())

        save_json(pipeline_state, "pipeline_state.json")
        logger.info(f"Phase 1 complete: {time.time()-t1:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 2: REACTOR-SCALE SIMULATION (CANTERA)
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 2 <= end_phase:
        print_banner("PHASE 2: CANTERA REACTOR SIMULATION")
        t2 = time.time()

        from pipeline.reactor_mechanisms import write_full_mechanism, write_gri30_subset
        from pipeline.reactor_models import run_reactor_sweep

        # Write gas-phase mechanism
        write_gri30_subset()

        # For each top catalyst, generate mechanism and run reactor sweep
        if 'top_catalysts' not in dir():
            # Load from previous phase
            import pandas as pd
            db_path = SCREENING_DIR / "ga_full_database.csv"
            if db_path.exists():
                valid_db = pd.read_csv(db_path)
                valid_db = valid_db[valid_db['valid'] == True]
                top_catalysts = valid_db.nsmallest(config.top_k_reactor, 'E_act')
            else:
                logger.warning("No screening database found. Using mock catalysts.")
                top_catalysts = None

        reactor_results = []
        if top_catalysts is not None:
            for idx, row in top_catalysts.iterrows():
                cat_name = f"cat_{idx}"
                E_act = row.get('E_act', 0.8)

                # Generate Cantera mechanism
                mech_path = write_full_mechanism(
                    cat_name, E_act_CH4=E_act,
                    E_act_H_desorb=max(0.3, abs(row.get('dE_H', -0.5))),
                )

                # Run reactor sweep
                results = run_reactor_sweep(
                    cat_name, str(mech_path),
                    temperatures=list(config.reactor_temperatures),
                    reactor_types=list(config.reactor_types),
                )
                reactor_results.extend(results)
        else:
            # Mock: run 3 test catalysts
            for name, e_act in [('NiBi_10', 0.85), ('FeC_supported', 0.65), ('CuSn_20', 1.1)]:
                mech_path = write_full_mechanism(name, E_act_CH4=e_act)
                results = run_reactor_sweep(
                    name, str(mech_path),
                    temperatures=list(config.reactor_temperatures),
                    reactor_types=list(config.reactor_types),
                )
                reactor_results.extend(results)

        pipeline_state['phase2'] = {
            'n_simulations': len(reactor_results),
            'elapsed_s': time.time() - t2,
        }
        if reactor_results:
            best_conv = max(r.get('CH4_conversion', 0) for r in reactor_results)
            pipeline_state['phase2']['best_conversion'] = best_conv

        save_json(pipeline_state, "pipeline_state.json")
        logger.info(f"Phase 2 complete: {len(reactor_results)} simulations, {time.time()-t2:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 3: DFT VALIDATION (QUANTUM ESPRESSO)
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 3 <= end_phase:
        print_banner("PHASE 3: DFT VALIDATION")
        t3 = time.time()

        from pipeline.dft_validator import validate_catalyst

        dft_results = []
        if 'top_catalysts' in dir() and top_catalysts is not None:
            top_dft = top_catalysts.head(config.top_k_dft)
            for idx, row in top_dft.iterrows():
                try:
                    genome = eval(row['genome'])
                    cat_name = f"dft_cat_{idx}"
                    result = validate_catalyst(cat_name, genome, run_dft=config.run_dft)
                    dft_results.append(result)
                except Exception as e:
                    logger.error(f"DFT failed for cat_{idx}: {e}")
        else:
            # Mock validation
            mock_genomes = [
                ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000),
                ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, ('Cu',), 1, 0),
                ('SAC', 'Fe', 'N4', 'N-graphene'),
            ]
            for i, genome in enumerate(mock_genomes[:config.top_k_dft]):
                result = validate_catalyst(f"dft_mock_{i}", genome, run_dft=config.run_dft)
                dft_results.append(result)

        pipeline_state['phase3'] = {
            'n_validated': len(dft_results),
            'n_converged': sum(1 for r in dft_results if r.get('converged', False)),
            'elapsed_s': time.time() - t3,
        }
        save_json(pipeline_state, "pipeline_state.json")
        logger.info(f"Phase 3 complete: {time.time()-t3:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 4: VQE TRANSITION STATE (CUDA-Q)
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 4 <= end_phase:
        print_banner("PHASE 4: CUDA-Q VQE TRANSITION STATE")
        t4 = time.time()

        from pipeline.vqe_transition_state import validate_transition_state

        vqe_results = []
        target = 'nvidia' if config.run_vqe else 'default'

        for i in range(min(config.top_k_vqe, 3)):
            result = validate_transition_state(f"champion_{i}", "CH_split", target=target)
            vqe_results.append(result)

        pipeline_state['phase4'] = {
            'n_vqe_runs': len(vqe_results),
            'elapsed_s': time.time() - t4,
        }
        save_json(pipeline_state, "pipeline_state.json")
        logger.info(f"Phase 4 complete: {time.time()-t4:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 5: FUEL CELL SCREENING & PEMFC MODELING
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 5 <= end_phase:
        print_banner("PHASE 5: FUEL CELL CATHODE SCREENING & PEMFC MODEL")
        t5 = time.time()

        from pipeline.fc_cathode_screener import run_cathode_screening, MEMBRANE_TYPES
        from pipeline.pemfc_model import PEMFCConfig, simulate_pemfc, sweep_membranes
        from pipeline.fuel_cell_stack import StackConfig, model_stack

        # Screen cathode catalysts
        cathode_df = run_cathode_screening()

        # Take top-K for PEMFC simulation
        valid_cathodes = cathode_df[cathode_df['valid'] == True].copy()
        if 'orr_overpotential_V' in valid_cathodes.columns:
            top_cathodes = valid_cathodes.nsmallest(config.fc_top_k_pemfc, 'orr_overpotential_V')
        else:
            top_cathodes = valid_cathodes.head(config.fc_top_k_pemfc)

        pemfc_results = []
        for _, row in top_cathodes.iterrows():
            cat_name = row['name']
            eta = row.get('orr_overpotential_V', 0.4)
            pgm = row.get('pgm_loading_mg_cm2', 0.0)

            # Sweep membranes
            mem_results = sweep_membranes(cat_name, eta)
            pemfc_results.extend(mem_results)

        # Stack model for the best catalyst + membrane combo
        if pemfc_results:
            # Optimize for highest efficiency and least overvoltage first and foremost,
            # while maximizing peak power output as much as possible.
            def fc_composite_score(r):
                eff = r.get('efficiency_at_rated', 0.0)
                power = r.get('peak_power_W_cm2', 0.0)
                eta = max(r.get('orr_overpotential_V', 0.4), 0.01)
                return (eff * power) / eta

            best_pemfc = max(pemfc_results, key=fc_composite_score)
            stack_config = StackConfig(
                n_cells=config.fc_stack_cells,
                cell_voltage_V=best_pemfc.get('peak_voltage_V', 0.65),
                current_density_A_cm2=best_pemfc.get('peak_current_A_cm2', 1.5),
            )
            stack_result = model_stack(stack_config)
        else:
            stack_result = {}

        pipeline_state['phase5'] = {
            'n_cathodes_screened': len(cathode_df),
            'n_valid': len(valid_cathodes),
            'n_pemfc_simulations': len(pemfc_results),
            'elapsed_s': time.time() - t5,
        }
        if pemfc_results:
            pipeline_state['phase5']['best_power_W_cm2'] = max(
                r.get('peak_power_W_cm2', 0) for r in pemfc_results
            )
            pipeline_state['phase5']['best_efficiency'] = max(
                r.get('efficiency_at_rated', 0) for r in pemfc_results
            )
            pipeline_state['phase5']['min_overpotential_V'] = min(
                r.get('orr_overpotential_V', 1.0) for r in pemfc_results
            )
        save_json(pipeline_state, "pipeline_state.json")
        logger.info(f"Phase 5 complete: {time.time()-t5:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 6: REPORT GENERATION
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 6 <= end_phase:
        print_banner("PHASE 6: REPORT GENERATION")
        t6 = time.time()

        from pipeline.report_generator import generate_full_report
        report_path = generate_full_report(pipeline_state)

        pipeline_state['phase6'] = {
            'report_path': str(report_path),
            'elapsed_s': time.time() - t6,
        }
        save_json(pipeline_state, "pipeline_state.json")

    # ═════════════════════════════════════════════════════════════════════════
    total_time = time.time() - t_total
    logger.info(f"\n{'='*70}")
    logger.info(f"  PIPELINE COMPLETE: {total_time:.0f}s ({total_time/3600:.1f} hours)")
    logger.info(f"{'='*70}")

    pipeline_state['total_elapsed_s'] = total_time
    save_json(pipeline_state, "pipeline_state.json")

    return pipeline_state


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Turquoise H₂ → Fuel Cell Pipeline')
    parser.add_argument('--phase', type=int, default=None, help='Run only this phase')
    parser.add_argument('--start', type=int, default=1, help='Start phase')
    parser.add_argument('--end', type=int, default=6, help='End phase')
    parser.add_argument('--quick', action='store_true', help='Quick mode (reduced parameters)')
    parser.add_argument('--no-dft', action='store_true', help='Skip DFT execution')
    parser.add_argument('--no-vqe', action='store_true', help='Skip VQE execution')
    args = parser.parse_args()

    config = PipelineConfig(
        quick_mode=args.quick,
        run_dft=not args.no_dft,
        run_vqe=not args.no_vqe,
    )

    if args.phase:
        run_pipeline(config, start_phase=args.phase, end_phase=args.phase)
    else:
        run_pipeline(config, start_phase=args.start, end_phase=args.end)
