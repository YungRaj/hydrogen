#!/usr/bin/env python3
"""
Master Pipeline Orchestrator.

Coordinates all 6 phases of the turquoise hydrogen → fuel cell pipeline:

  Phase 1: Deterministic Branch-and-Bound → Champion Archives
  Phase 2: Cantera Reactor Simulation (top catalysts × 3 reactor types)
  Phase 3: DFT Validation (top catalysts)
  Phase 4: VQE Transition State (top 3 catalysts)
  Phase 5: Fuel Cell Cathode Screening → PEMFC Modeling
  Phase 6: Report Generation

Usage:
    python -m pipeline.orchestrator [--phase N] [--quick]
"""

import sys
import argparse
from pathlib import Path
from typing import Optional
from dataclasses import dataclass, replace

sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.utils import (
    RESULTS_DIR, SCREENING_DIR, setup_logger,
)
from pipeline.reactors.modes import (
    DEFAULT_MODE, MODE_CHOICES, reactor_types_for_mode)
from pipeline.stages.orchestration import (
    PipelineComponents, PipelineRuntime, default_pipeline_components,
    default_pipeline_runtime)
from pipeline.stages.contracts import require_stage_outcome

logger = setup_logger('orchestrator', 'pipeline_orchestrator.log')


@dataclass
class PipelineConfig:
    """Pipeline-level configuration."""
    # Phase 1: Screening
    initial_fairchem_samples: int = 500       # Initial Fairchem evaluations
    branch_leaf_size: int = 1_000_000         # Exhaustive terminal range size
    branch_max_leaves: Optional[int] = None    # Staged execution; None = complete
    top_k_reactor: int = 50                  # Top K catalysts → reactor simulation
    top_k_dft: int = 14                  # At least one class champion → DFT
    top_k_vqe: int = 3                   # Top K → VQE

    # Phase 2: Reactor
    # Preserve the broad screening points and include the Ni reference band.
    reactor_temperatures: tuple = (773.15, 900.0, 923.15, 973.15, 1100.0, 1300.0)
    reactor_types: Optional[tuple] = None  # None derives routing from mode
    multiphysics_results_dir: Optional[str] = None
    # Named solids judge. B6-7: literature Ni cell (ni_np_lit), not cat_9.
    # None = best non-H-parked solids at the headline T band.
    solids_judge_catalyst: Optional[str] = 'ni_np_lit'
    # Ni judge headline is the 650–700 °C filament ROI. 1300 K is the ADR
    # ceiling, not this headline (X>X_eq at the hot end is a flag).
    solids_headline_t_min: float = 923.15
    solids_headline_t_max: Optional[float] = 973.15

    # Phase 5: Fuel Cell
    fc_top_k_pemfc: int = 20            # Top cathode catalysts → PEMFC model
    fc_stack_cells: int = 300            # Stack cells

    # Runtime
    run_dft: bool = True                 # Actually execute pw.x
    run_vqe: bool = True                 # Actually execute CUDA-Q
    quick_mode: bool = False             # Reduced parameters for testing
    pyrolysis_mode: str = DEFAULT_MODE
    allow_mock_inputs: bool = False       # Explicit test-only opt-in


def normalized_pipeline_config(config: PipelineConfig) -> PipelineConfig:
    """Apply runtime presets while preserving the caller's reactor conditions.

    Args:
        config: Configuration controlling this operation.

    Returns:
        Computed `PipelineConfig` result.
    """
    effective = replace(config)
    if effective.quick_mode:
        effective = replace(
            effective, initial_fairchem_samples=50, branch_leaf_size=10_000,
            branch_max_leaves=1, top_k_reactor=10, top_k_dft=3,
            top_k_vqe=1, fc_top_k_pemfc=5)
    return effective


def run_pipeline(config: PipelineConfig | None = None,
                 start_phase: int = 1, end_phase: int = 6,
                 runtime: PipelineRuntime | None = None,
                 components: PipelineComponents | None = None):
    """
        Execute the full multi-scale simulation pipeline.

    Args:
        config: Configuration controlling this operation.
        start_phase: Start phase used by this operation.
        end_phase: End phase used by this operation.
        runtime: Runtime used by this operation.
        components: Components used by this operation.

    Returns:
        Computed result described above.
    """
    config = normalized_pipeline_config(config or PipelineConfig())
    runtime = runtime or default_pipeline_runtime()
    components = components or default_pipeline_components()
    t_total = runtime.clock()
    runtime.banner("TURQUOISE HYDROGEN → FUEL CELL: MULTI-SCALE PIPELINE")
    logger.info(f"Starting pipeline: phases {start_phase}–{end_phase}")
    logger.info(f"Configuration: quick_mode={config.quick_mode}, pyrolysis_mode={config.pyrolysis_mode}")

    # Propagate pyrolysis mode to env
    runtime.select_pathway_mode(config.pyrolysis_mode)
    selected_reactors = (config.reactor_types if config.reactor_types is not None
                         else reactor_types_for_mode(config.pyrolysis_mode))
    multiphysics_results_dir = (config.multiphysics_results_dir or
                                str(RESULTS_DIR / 'multiphysics'))

    pipeline_state = runtime.load_state()

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 1: DETERMINISTIC BRANCH-AND-BOUND
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 1 <= end_phase:
        runtime.banner("PHASE 1: DETERMINISTIC BRANCH-AND-BOUND DISCOVERY")
        t1 = runtime.clock()

        outcome = require_stage_outcome(components.discovery(
            initial_samples=config.initial_fairchem_samples,
            leaf_size=config.branch_leaf_size,
            max_leaves=config.branch_max_leaves,
            top_k_reactor=config.top_k_reactor,
            top_k_dft=config.top_k_dft), stage='discovery', required_products=(
                'design_space_sizes', 'pareto_genomes', 'screening_database',
                'top_catalysts', 'dft_candidates'))
        sizes = outcome.products['design_space_sizes']
        logger.info(f"Design space: {sizes['TOTAL']:,} total configurations")
        for cls, size in sizes.items():
            if cls != 'TOTAL':
                logger.info(f"  {cls}: {size:,}")

        pareto_genomes = outcome.products['pareto_genomes']
        screening_db = outcome.products['screening_database']
        top_catalysts = outcome.products['top_catalysts']
        dft_candidates = outcome.products['dft_candidates']
        pipeline_state['phase1'] = {
            **outcome.state, 'elapsed_s': runtime.clock() - t1}

        runtime.save_state(pipeline_state)
        logger.info(f"Phase 1 complete: {runtime.clock()-t1:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 2: REACTOR-SCALE SIMULATION (CANTERA)
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 2 <= end_phase:
        runtime.banner("PHASE 2: PATHWAY-SPECIFIC REACTOR SIMULATION")
        t2 = runtime.clock()

        # For each top catalyst, generate mechanism and run reactor sweep
        if 'top_catalysts' not in dir():
            # Load from previous phase
            db_path = SCREENING_DIR / "ga_full_database.csv"
            restored = components.load_candidates(
                db_path, top_k_reactor=config.top_k_reactor,
                top_k_dft=config.top_k_dft)
            if restored is not None:
                screening_db = restored['screening_database']
                top_catalysts = restored['top_catalysts']
                dft_candidates = restored['dft_candidates']
            else:
                if not config.allow_mock_inputs:
                    raise RuntimeError(
                        'screening database is required; mock catalyst fallback is disabled')
                logger.warning("No screening database found. Using mock catalysts.")
                top_catalysts = None

        outcome = require_stage_outcome(components.reactor_batch(
            top_catalysts, temperatures=config.reactor_temperatures,
            reactor_types=selected_reactors,
            pathway_mode=config.pyrolysis_mode,
            multiphysics_results_dir=multiphysics_results_dir,
            allow_mock_inputs=config.allow_mock_inputs,
            judge_catalyst=config.solids_judge_catalyst,
            headline_t_min=config.solids_headline_t_min,
            headline_t_max=config.solids_headline_t_max), stage='reactor_batch',
            required_products=('reactor_results',))
        reactor_results = outcome.products['reactor_results']
        pipeline_state['phase2'] = {
            **outcome.state, 'elapsed_s': runtime.clock() - t2}

        runtime.save_state(pipeline_state)
        logger.info(f"Phase 2 complete: {len(reactor_results)} simulations, {runtime.clock()-t2:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 3: DFT VALIDATION (QUANTUM ESPRESSO)
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 3 <= end_phase:
        runtime.banner("PHASE 3: DFT VALIDATION")
        t3 = runtime.clock()

        if 'dft_candidates' not in dir():
            restored = components.load_candidates(
                SCREENING_DIR / "ga_full_database.csv",
                top_k_reactor=config.top_k_reactor,
                top_k_dft=config.top_k_dft)
            dft_candidates = (restored.get('dft_candidates')
                              if restored is not None else None)
        if dft_candidates is not None:
            outcome = components.dft(
                dft_candidates, top_k=config.top_k_dft,
                execute_dft=config.run_dft, error_sink=logger.error)
        else:
            # Mock validation
            if not config.allow_mock_inputs:
                raise RuntimeError(
                    'screening candidates are required; mock DFT fallback is disabled')
            mock_genomes = [
                ('MoltenMetal', 'Bi', 'Ni', 10.0, 1000),
                ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, ('Cu',), 1, 0),
                ('SAC', 'Fe', 'N4', 'N-graphene'),
            ]
            outcome = components.dft(
                mock_genomes, top_k=config.top_k_dft,
                execute_dft=config.run_dft, name_prefix='dft_mock',
                error_sink=logger.error)

        outcome = require_stage_outcome(
            outcome, stage='dft', required_products=('dft_results',))

        dft_results = outcome.products['dft_results']
        pipeline_state['phase3'] = {
            **outcome.state, 'elapsed_s': runtime.clock() - t3}
        runtime.save_state(pipeline_state)
        logger.info(f"Phase 3 complete: {runtime.clock()-t3:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 4: VQE TRANSITION STATE (CUDA-Q)
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 4 <= end_phase:
        runtime.banner("PHASE 4: CUDA-Q VQE TRANSITION STATE")
        t4 = runtime.clock()

        outcome = require_stage_outcome(components.vqe(
            top_k=config.top_k_vqe, execute_quantum=config.run_vqe),
            stage='vqe', required_products=('vqe_results',))
        vqe_results = outcome.products['vqe_results']
        pipeline_state['phase4'] = {
            **outcome.state, 'elapsed_s': runtime.clock() - t4}
        runtime.save_state(pipeline_state)
        logger.info(f"Phase 4 complete: {runtime.clock()-t4:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 5: FUEL CELL SCREENING & PEMFC MODELING
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 5 <= end_phase:
        runtime.banner("PHASE 5: FUEL CELL CATHODE SCREENING & PEMFC MODEL")
        t5 = runtime.clock()

        outcome = require_stage_outcome(components.fuel_cell(
            top_k_pemfc=config.fc_top_k_pemfc,
            stack_cells=config.fc_stack_cells), stage='fuel_cell',
            required_products=('cathode_database', 'valid_cathodes',
                               'pemfc_results', 'stack_result'))
        cathode_df = outcome.products['cathode_database']
        valid_cathodes = outcome.products['valid_cathodes']
        pemfc_results = outcome.products['pemfc_results']
        stack_result = outcome.products['stack_result']
        pipeline_state['phase5'] = {
            **outcome.state, 'elapsed_s': runtime.clock() - t5}
        runtime.save_state(pipeline_state)
        logger.info(f"Phase 5 complete: {runtime.clock()-t5:.0f}s")

    # ═════════════════════════════════════════════════════════════════════════
    # PHASE 6: REPORT GENERATION
    # ═════════════════════════════════════════════════════════════════════════
    if start_phase <= 6 <= end_phase:
        runtime.banner("PHASE 6: REPORT GENERATION")
        t6 = runtime.clock()

        outcome = require_stage_outcome(
            components.report(pipeline_state), stage='report',
            required_products=('report_path',))
        report_path = outcome.products['report_path']

        pipeline_state['phase6'] = {
            **outcome.state,
            'elapsed_s': runtime.clock() - t6,
        }
        runtime.save_state(pipeline_state)

    # ═════════════════════════════════════════════════════════════════════════
    total_time = runtime.clock() - t_total
    logger.info(f"\n{'='*70}")
    logger.info(f"  PIPELINE COMPLETE: {total_time:.0f}s ({total_time/3600:.1f} hours)")
    logger.info(f"{'='*70}")

    pipeline_state['total_elapsed_s'] = total_time
    runtime.save_state(pipeline_state)

    return pipeline_state


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Turquoise H₂ → Fuel Cell Pipeline')
    parser.add_argument('--phase', type=int, default=None, help='Run only this phase')
    parser.add_argument('--start', type=int, default=1, help='Start phase')
    parser.add_argument('--end', type=int, default=6, help='End phase')
    parser.add_argument('--quick', action='store_true', help='Quick mode (reduced parameters)')
    parser.add_argument('--no-dft', action='store_true', help='Skip DFT execution')
    parser.add_argument('--no-vqe', action='store_true', help='Skip VQE execution')
    parser.add_argument('--mode', choices=MODE_CHOICES, default=DEFAULT_MODE,
                        help=f'Methane-conversion pathway (default: {DEFAULT_MODE})')
    parser.add_argument('--multiphysics-results-dir', default=None)
    args = parser.parse_args()

    config = PipelineConfig(
        quick_mode=args.quick,
        run_dft=not args.no_dft,
        run_vqe=not args.no_vqe,
        pyrolysis_mode=args.mode,
        multiphysics_results_dir=args.multiphysics_results_dir,
    )

    if args.phase:
        run_pipeline(config, start_phase=args.phase, end_phase=args.phase)
    else:
        run_pipeline(config, start_phase=args.start, end_phase=args.end)
