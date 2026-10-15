#!/usr/bin/env python3
"""
GPU-SATURATED DETERMINISTIC DISCOVERY CAMPAIGN

The only production candidate-search strategy is persistent, exhaustive
branch-and-bound over all 14 classes in the 21.1B indexed encoded space.

``main`` reads as the campaign outline: configure, gate the design space,
record provenance, then run the six phases. Each phase records its summary
in ``pipeline_state`` and persists it. ``pipeline.utils`` reads the results
directory from the environment at import time, so every pipeline import
happens after ``configure_environment``.
"""

import os
import sys
import time
import json
import ast
import argparse
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from pipeline.reactors.modes import (
    DEFAULT_MODE, MODE_CHOICES, THERMAL_REACTOR_TEMPERATURES_K,
    reactor_temperatures_for_mode, reactor_types_for_mode)

REPO_DIR = Path(__file__).parent
REACTOR_TEMPERATURES_K = list(THERMAL_REACTOR_TEMPERATURES_K)


def _reactor_temperatures(mode: str) -> list[float]:
    """Return the shared mode temperature grid from pipeline.reactors.modes."""
    return list(reactor_temperatures_for_mode(mode))
MAX_REACTOR_CANDIDATES = 20
MAX_VQE_CANDIDATES = 5
PEMFC_CANDIDATES = 30


def _branch_config_kwargs(
    args: argparse.Namespace, deadline_epoch_s: float
) -> dict[str, object]:
    """Translate shared CLI settings into either branch-discovery config."""
    return {
        'initial_fairchem_samples': args.calibration_probes,
        'fairchem_eval_top_k': args.validation_batch,
        'exhaustive_batch_size': args.scan_batch_size,
        'branch_leaf_size': args.branch_leaf_size,
        'branch_probe_count': args.branch_probes,
        'branch_max_leaves': (
            None if args.branch_max_leaves == 0 else args.branch_max_leaves
        ),
        'expected_space_size': args.expected_space_size,
        'max_runtime_s': (
            None
            if args.hours == 0
            else max(0.0, deadline_epoch_s - time.time())
        ),
        'prior_art_db': args.prior_art_db,
        'min_validation_per_class': args.min_validation_per_class,
        'min_resolved_leaves_per_class': args.branch_class_floor,
        'branch_exploration_interval': args.branch_exploration_interval,
        'refresh_pending_priorities': args.branch_priority_refresh,
        'scan_workers': args.scan_workers,
        'workers_per_gpu': args.workers_per_gpu,
        'calibration_artifact': getattr(args, 'calibration_json', None),
    }


# ─── Configuration ──────────────────────────────────────────────────────────


def _build_parser() -> argparse.ArgumentParser:
    """Return the campaign command-line interface."""
    parser = argparse.ArgumentParser(
        description='GPU-Saturated Turquoise H₂ Catalyst Discovery'
    )
    parser.add_argument('--calibration-probes', type=int, default=500,
                        help='Deterministic binary-tree probes used to calibrate the surrogate')
    parser.add_argument('--validation-batch', type=int, default=500,
                        help='Branch/archive champions sent to the atomistic model')
    parser.add_argument('--min-validation-per-class', type=int, default=2,
                        help='Reserved atomistic validation quota per represented class')
    parser.add_argument('--hours', type=float, default=0,
                        help='Max wall-clock hours. 0 = unlimited (default: 0)')
    parser.add_argument('--top-k', type=int, default=200,
                        help='Top-K for reactor simulation (default: 200)')
    parser.add_argument('--no-dft', action='store_true')
    parser.add_argument('--no-vqe', action='store_true')
    parser.add_argument('--mode', choices=MODE_CHOICES, default=DEFAULT_MODE,
                        help=f'Methane-conversion pathway (default: {DEFAULT_MODE})')
    parser.add_argument('--ntec-conditions-json', default='',
                        help=('JSON with measured operating conditions and paired '
                              'NTEC/control effect calibration'))
    parser.add_argument('--electrochemical-conditions-json', default='',
                        help=('JSON measured electrochemical operating point; '
                              'electrolyte_phase is aqueous or molten'))
    parser.add_argument('--workers-per-gpu', type=int, default=32,
                        help='Candidate threads sharing each GPU in the batched '
                             'eSEN engine (relaxation CPU work overlaps inference)')
    parser.add_argument('--scan-batch-size', type=int, default=65536)
    parser.add_argument('--scan-workers', type=int, default=8,
                        help='Independent deterministic scanner shards per resolved leaf')
    parser.add_argument('--branch-leaf-size', type=int, default=1_000_000)
    parser.add_argument('--branch-probes', type=int, default=9)
    parser.add_argument('--branch-max-leaves', type=int, default=0,
                        help='Process at most this many leaves per run; 0 continues until complete')
    parser.add_argument('--branch-class-floor', type=int, default=1,
                        help='Resolve at least this many leaves per class before pure exploitation')
    parser.add_argument('--branch-exploration-interval', type=int, default=4,
                        help='After class floors, reserve every Nth leaf for the least-covered class')
    parser.add_argument('--branch-priority-refresh', type=int, default=10000,
                        help='Re-score up to this many pending nodes when resuming')
    parser.add_argument('--expected-space-size', type=int, default=21_092_657_271,
                        help='Fail if the indexed population denominator differs')
    parser.add_argument('--results-dir', default='results',
                        help='Isolated output root for this campaign (default: results)')
    parser.add_argument('--mechanisms-dir', default=None,
                        help='Generated Cantera mechanisms (default: <results-dir>/mechanisms)')
    parser.add_argument('--prior-art-db', default=None,
                        help='Prior-art SQLite database (default: <results-dir>/prior_art.sqlite)')
    parser.add_argument('--prior-art-csv', action='append', default=[],
                        help='CSV registry to import; repeat for multiple sources')
    parser.add_argument('--calibration-json', default=None,
                        help=('Active-learning DFT calibration JSON '
                              '(default: <results-dir>/dft/active_learning_calibration.json)'))
    parser.add_argument('--final-campaign', action='store_true',
                        help='Fail closed unless coverage and prior-art readiness requirements pass')
    parser.add_argument('--evidence-manifest', default=None,
                        help='Measured/validated evidence counts required by --final-campaign')
    parser.add_argument(
        '--kinetics-validation-dir', default=None,
        help=('Directory containing <candidate_id>/pyrolysis_validation.json; '
              'complete matching campaigns replace Cantera template barriers'))
    parser.add_argument(
        '--multiphysics-results-dir', default=None,
        help=('Required OpenFOAM/MOOSE artifacts (default: '
              '<results-dir>/multiphysics)'))
    parser.add_argument('--qe-mpi-ranks', type=int, default=4)
    parser.add_argument('--qe-omp-threads', type=int, default=1)
    default_qe_concurrency = max(1, min(4, (os.cpu_count() or 1) // 4))
    parser.add_argument('--qe-max-concurrent', type=int, default=default_qe_concurrency,
                        help='Maximum independent candidate DFT jobs launched together')
    return parser


def _resolve_paths(args: argparse.Namespace) -> Path:
    """Make every output path absolute and default it under the results root."""
    results_dir = Path(args.results_dir).expanduser().resolve()
    args.results_dir = str(results_dir)
    mechanisms_dir = (Path(args.mechanisms_dir).expanduser().resolve()
                      if args.mechanisms_dir else results_dir / 'mechanisms')
    args.mechanisms_dir = str(mechanisms_dir)
    if args.prior_art_db is None:
        args.prior_art_db = str(results_dir / 'prior_art.sqlite')
    if getattr(args, 'calibration_json', None) is None:
        args.calibration_json = str(
            results_dir / 'dft' / 'active_learning_calibration.json'
        )
    else:
        args.calibration_json = str(
            Path(args.calibration_json).expanduser().resolve()
        )
    if args.evidence_manifest is None:
        args.evidence_manifest = str(results_dir / 'evidence_manifest.json')
    if args.kinetics_validation_dir is None:
        args.kinetics_validation_dir = str(results_dir / 'dft' / 'pyrolysis_kinetics')
    if args.multiphysics_results_dir is None:
        args.multiphysics_results_dir = str(results_dir / 'multiphysics')
    else:
        args.multiphysics_results_dir = str(
            Path(args.multiphysics_results_dir).expanduser().resolve())
    os.environ['HYDROGEN_RESULTS_DIR'] = str(results_dir)
    os.environ['HYDROGEN_MECHANISMS_DIR'] = str(mechanisms_dir)
    return results_dir


def _validate_resources(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Reject budgets the branch search or the QE allocation cannot honour."""
    if args.scan_workers < 1 or args.workers_per_gpu < 1 or args.qe_mpi_ranks < 1 or \
            args.qe_omp_threads < 1 or args.qe_max_concurrent < 1:
        parser.error('scanner and QE resource dimensions must be positive')
    if args.calibration_probes < 20:
        parser.error('--calibration-probes must be at least 20 for the tree ranker')
    required_validation = 14 * args.min_validation_per_class
    if args.validation_batch < required_validation:
        parser.error(
            f'--validation-batch must be at least {required_validation} to reserve '
            f'{args.min_validation_per_class} candidate(s) across all 14 classes')
    requested_qe_cpus = (args.qe_mpi_ranks * args.qe_omp_threads *
                         args.qe_max_concurrent)
    available_cpus = os.cpu_count() or 1
    if requested_qe_cpus > available_cpus:
        parser.error(
            f'QE allocation requests {requested_qe_cpus} CPUs, '
            f'but only {available_cpus} are visible')


def _configure_environment(args: argparse.Namespace):
    """Export run settings, load prior art, and print the hardware banner.

    Returns:
        The prior-art registry, GPU count and GPU names.
    """
    # QE executables are resolved at execution time only from the explicit
    # PW_X/NEB_X and MPIEXEC set by the pinned Blackwell GPU activation file.
    os.environ['QE_MPI_RANKS'] = str(args.qe_mpi_ranks)
    os.environ['QE_OMP_THREADS'] = str(args.qe_omp_threads)

    from pipeline.evidence.prior_art import PriorArtRegistry
    prior_registry = PriorArtRegistry(args.prior_art_db)
    for prior_csv in args.prior_art_csv:
        prior_registry.import_csv(prior_csv)

    os.environ['PYROLYSIS_MODE'] = args.mode
    if args.ntec_conditions_json:
        os.environ['NTEC_CONDITIONS_JSON'] = args.ntec_conditions_json
    if args.electrochemical_conditions_json:
        os.environ['ELECTROCHEMICAL_CONDITIONS_JSON'] = \
            args.electrochemical_conditions_json
    import torch
    n_gpus = torch.cuda.device_count()
    gpu_names = [torch.cuda.get_device_name(i) for i in range(n_gpus)]
    gpu_mem = [torch.cuda.get_device_properties(i).total_memory / 1e9 for i in range(n_gpus)]

    print("=" * 80)
    print("  GPU-SATURATED CATALYST DISCOVERY CAMPAIGN v2")
    print("=" * 80)
    for i in range(n_gpus):
        print(f"  GPU[{i}] {gpu_names[i]} — {gpu_mem[i]:.1f} GB")
    print(f"  CPUs: {os.cpu_count()} cores")
    print(f"  Mode: {args.mode.upper()}")
    print("  Search: catalyst-guided, coverage-safe branch-and-bound")
    print(f"  Calibration probes: {args.calibration_probes:,}")
    print(f"  Validation batch: {args.validation_batch:,}")
    print(f"  Prior-art records: {prior_registry.count():,}")
    print("=" * 80)
    return prior_registry, n_gpus, gpu_names


def _load_hf_token() -> None:
    """Load the HuggingFace token needed for the gated eSEN surface model.

    Uses the environment, ``$HF_HOME/token`` or the repo ``.hf_token``; prompts
    only on an interactive terminal; and fails fast otherwise, because
    screening has no non-gated fallback model.
    """
    from pipeline.screening.surface_calculator import _ensure_hf_token

    if not _ensure_hf_token() and sys.stdin.isatty():
        print("\n  HuggingFace token needed for the eSEN OC25 surface model")
        print("    Get one at: https://huggingface.co/settings/tokens")
        token = input("  Enter HF token: ").strip()
        if token:
            os.environ['HF_TOKEN'] = token
    if not os.environ.get('HF_TOKEN'):
        raise SystemExit(
            'No HuggingFace token: set HF_TOKEN or run `huggingface-cli login` '
            '(writes $HF_HOME/token); eSEN screening cannot run without it')
    print("  ✓ HuggingFace token loaded (eSEN surface model enabled)")


def _check_design_space(args: argparse.Namespace) -> None:
    """Fail unless the indexed denominator and the canonical-identity audit pass."""
    from pipeline.search.design_space import estimate_design_space_size, ALL_MATERIAL_CLASSES
    from pipeline.utils import save_json

    sizes = estimate_design_space_size()
    if args.expected_space_size != sizes['TOTAL']:
        raise SystemExit(
            f"Population denominator mismatch: --expected-space-size={args.expected_space_size:,}, "
            f"but this repository indexes {sizes['TOTAL']:,}. Update the design-space "
            "definition or use the verified denominator; do not label it 25.3B."
        )
    from pipeline.evidence.design_space_audit import audit_design_space
    design_audit = audit_design_space(sample_per_class=2048)
    save_json(design_audit, 'design_space_audit.json', subdir='')
    if not design_audit['valid']:
        raise SystemExit(f"Design-space audit failed: {design_audit['failures']}")
    print(f"\n  Design Space: {sizes['TOTAL']:,} ({sizes['TOTAL']/1e9:.1f}B)")
    print(f"  Canonical identities: {design_audit['canonical_total']:,}")
    for cls in ALL_MATERIAL_CLASSES:
        print(f"    {cls:20s}: {sizes[cls]:>15,}")


def _get_git_sha():
    try:
        return subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=str(REPO_DIR), stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return 'unknown'


def _get_conda_env():
    try:
        return subprocess.check_output(
            ['conda', 'list', '--export'],
            stderr=subprocess.DEVNULL
        ).decode()[:2000]  # truncate for readability
    except Exception:
        return 'unavailable'


def _record_provenance(args: argparse.Namespace, n_gpus: int, gpu_names: list) -> dict:
    """Start ``pipeline_state.json`` with reproducibility metadata."""
    from pipeline.utils import save_json

    pipeline_state = {
        'provenance': {
            'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
            'git_sha': _get_git_sha(),
            'cli_args': vars(args),
            'n_gpus': n_gpus,
            'gpu_names': gpu_names,
            'conda_env_snapshot': _get_conda_env(),
        }
    }
    save_json(pipeline_state, "pipeline_state.json")
    return pipeline_state


# ─── Phases ─────────────────────────────────────────────────────────────────


def _run_phase1_discovery(args, pipeline_state: dict, t_deadline: float):
    """Branch-and-bound discovery, then admissibility-first slate selection.

    Returns:
        Reactor slate and DFT-resolution slate.
    """
    from pipeline.utils import print_banner, save_json

    print_banner("PHASE 1: DETERMINISTIC BRANCH-AND-BOUND + META ESEN-SM")
    t1 = time.time()

    from pipeline.screening.genetic_optimizer import run_branch_discovery, BranchDiscoveryConfig

    branch_config = BranchDiscoveryConfig(**_branch_config_kwargs(args, t_deadline))

    pareto_genomes, screening_db = run_branch_discovery(branch_config)

    from pipeline.search.scope import scope_pyrolysis_pool
    from pipeline.screening.stage_selection import (
        annotate_evidence, select_for_pyrolysis_reactor, select_for_validation)
    valid_db = screening_db[screening_db['valid'] == True].copy()
    evidence_db = annotate_evidence(screening_db, 'E_act')
    # Same split as run_discovery_stage: admissibility on the full table,
    # then each slate applies its own validity rule. Filtering valid==True
    # here drops unresolved rows the DFT rescue path is supposed to keep.
    pool, _ = scope_pyrolysis_pool(screening_db)
    top_catalysts = select_for_pyrolysis_reactor(
        pool, args.top_k, 'E_act',
        min_per_class=args.min_validation_per_class)
    from pipeline.screening.surface_screener import upgrade_candidates_with_mlip_neb
    top_catalysts = upgrade_candidates_with_mlip_neb(top_catalysts)
    dft_candidates = select_for_validation(
        pool, min(args.validation_batch, max(len(pool), 1)), 'E_act',
        min_per_class=args.min_validation_per_class)

    pipeline_state['phase1'] = {
        'pareto_size': len(pareto_genomes),
        'total_evaluated': len(screening_db),
        'valid_count': len(valid_db),
        'top_catalysts_count': len(top_catalysts),
        'dft_resolution_count': len(dft_candidates),
        'candidate_dispositions': evidence_db[
            'candidate_disposition'].value_counts().to_dict(),
        'elapsed_s': time.time() - t1,
        'search_strategy': 'deterministic_branch_and_bound',
    }
    # Censored barriers sit at the clip floor; only metric-eligible rows
    # (the same rule run_discovery_stage uses) may report a best barrier.
    from pipeline.screening.small_data_ranker import screening_metric_eligibility
    eligible_db = screening_db.loc[
        screening_metric_eligibility(screening_db, 'turquoise_hydrogen')]
    pipeline_state['phase1']['metric_eligible_count'] = len(eligible_db)
    if len(eligible_db) > 0:
        pipeline_state['phase1']['best_E_act'] = float(eligible_db['E_act'].min())
        if 'coking_index' in eligible_db.columns and eligible_db['coking_index'].notna().any():
            pipeline_state['phase1']['best_coking'] = f"{float(eligible_db['coking_index'].dropna().max()):.3f} eV"
    save_json(pipeline_state, "pipeline_state.json")
    print(f"\n  Phase 1: {time.time()-t1:.0f}s | Evaluated: {len(screening_db):,}")
    return top_catalysts, dft_candidates


def _load_kinetics_validation(args, candidate_key: str):
    """Return a complete matching pyrolysis validation campaign, if present."""
    validation_path = (Path(args.kinetics_validation_dir) / candidate_key /
                       'pyrolysis_validation.json')
    if validation_path.is_file():
        available = json.loads(validation_path.read_text())
        if available.get('complete'):
            return available
    return None


def _reactor_candidate_key(row, index: int) -> str:
    """Return the candidate id; an unsafe path component aborts the campaign."""
    candidate_key = str(row.get('candidate_id', f'catalyst_{index}'))
    if Path(candidate_key).name != candidate_key:
        raise ValueError('candidate_id is not a safe path component')
    return candidate_key


def _simulate_reactor_candidate(args, row, candidate_key: str, simulate_candidate,
                                sweep_records: list):
    """Run one candidate's reactor sweep and summarise its best condition."""
    e_act = row.get('E_act', 1.0)
    stage_result = simulate_candidate(
        row, candidate_key, _reactor_temperatures(args.mode),
        reactor_types=reactor_types_for_mode(args.mode),
        forbid_mock=True,
        kinetics_validation=_load_kinetics_validation(args, candidate_key),
        pathway_mode=args.mode,
        multiphysics_results_dir=args.multiphysics_results_dir)
    sweep = stage_result['sweep']
    sweep_records.extend(sweep)
    best_condition = stage_result['best_condition']
    return {
        'catalyst': candidate_key,
        'candidate_id': candidate_key,
        'E_act': e_act,
        'best_conversion': best_condition.get('CH4_conversion', 0),
        'n_conditions': len(sweep),
        'completed_conditions': stage_result['completed_conditions'],
        'failed_conditions': stage_result['failed_conditions'],
        'sweep_status': stage_result['sweep_status'],
        # Completed reactor records carry T_K; temperature_K is only set on
        # failed-condition records.
        'temperature_K': best_condition.get('T_K', best_condition.get('temperature_K')),
        **{k: best_condition.get(k) for k in (
            'H2_selectivity', 'CH4_conversion',
            'deactivation_fraction_per_h', 'coke_fraction',
            'kinetics_status', 'reactor_evidence_tier',
            'can_exclude_candidate')},
    }


def _tea_estimates(reactor_results: list) -> list:
    """Pyrolysis TEA ($/kg H₂) for candidates above 1% conversion."""
    from pipeline.economics.tea import estimate_scenario_range
    tea_results = []
    for r in reactor_results:
        conv = r.get('best_conversion', 0)
        if conv > 0.01:
            estimate = estimate_scenario_range(conv)
            tea_results.append({
                'catalyst': r['catalyst'],
                'h2_cost_usd_kg': estimate['estimates']['base'][
                    'h2_cost_usd_kg'],
                **estimate,
            })
    return tea_results


def _simulate_judge_reference(args, simulate_candidate, sweep_records: list) -> None:
    """Add the literature Ni judge's sweep to the scorecard records only.

    The reference is never a candidate: it is excluded from reactor_results,
    TEA and viability, exactly like the orchestrator reactor batch.
    """
    from pipeline.reactors.scorecard import NI_JUDGE_CATALYST
    from pipeline.stages.candidate_io import load_reactor_reference

    reactor_types = reactor_types_for_mode(args.mode)
    if not {'PFR', 'Fluidized'}.intersection(reactor_types):
        return
    if any(r.get('catalyst_name') == NI_JUDGE_CATALYST for r in sweep_records):
        return
    reference = load_reactor_reference(NI_JUDGE_CATALYST)
    if reference is None:
        return
    try:
        result = simulate_candidate(
            reference, NI_JUDGE_CATALYST, _reactor_temperatures(args.mode),
            reactor_types=reactor_types, forbid_mock=True,
            pathway_mode=args.mode,
            multiphysics_results_dir=args.multiphysics_results_dir)
        sweep_records.extend(result['sweep'])
    except Exception as e:
        print(f"    Judge reference error: {e}")


class _PrintLogger:
    def info(self, msg):
        print(f"  {msg}")


def _run_phase2_reactors(args, pipeline_state: dict, top_catalysts) -> None:
    """Pathway-specific reactor sweeps, solids scorecard, TEA and viability."""
    from pipeline.utils import print_banner, save_json

    print_banner("PHASE 2: PATHWAY-SPECIFIC REACTOR SIMULATION")
    t2 = time.time()
    try:
        from pipeline.stages.reactor import simulate_candidate
        from pipeline.reactors.equilibrium import run_equilibrium_sweep
        from pipeline.reactors.scorecard import (
            NI_JUDGE_CATALYST, NI_JUDGE_HEADLINE_T_MAX, NI_JUDGE_HEADLINE_T_MIN,
            build_solids_scorecard, log_solids_scorecard)

        eq_result = run_equilibrium_sweep()
        print(f"  Equilibrium check within_tol={eq_result.get('within_tolerance')} "
              f"worst_err={eq_result.get('worst_abs_error')}")

        reactor_results = []
        reactor_sweep_records = []
        n_reactor = min(MAX_REACTOR_CANDIDATES, len(top_catalysts))
        for i, (_, row) in enumerate(top_catalysts.head(n_reactor).iterrows()):
            candidate_key = _reactor_candidate_key(row, i)
            print(f"  Reactor sim {i+1}/{n_reactor}: E_act={row.get('E_act', 1.0):.3f} eV")
            try:
                reactor_results.append(_simulate_reactor_candidate(
                    args, row, candidate_key, simulate_candidate, reactor_sweep_records))
            except Exception as e:
                print(f"    Reactor error: {e}")

        tea_results = _tea_estimates(reactor_results)
        _simulate_judge_reference(args, simulate_candidate, reactor_sweep_records)
        scorecard = build_solids_scorecard(
            reactor_sweep_records, judge_catalyst=NI_JUDGE_CATALYST,
            headline_t_min=NI_JUDGE_HEADLINE_T_MIN,
            headline_t_max=NI_JUDGE_HEADLINE_T_MAX)
        save_json(scorecard, 'phase2_solids_scorecard.json', subdir='reactor')
        log_solids_scorecard(scorecard, _PrintLogger())
        pipeline_state['phase2'] = {
            'catalysts_attempted': n_reactor,
            'catalysts_simulated': len(reactor_results),
            'complete_sweeps': sum(
                r.get('sweep_status') == 'complete' for r in reactor_results),
            'partial_sweeps': sum(
                r.get('sweep_status') == 'partial' for r in reactor_results),
            'failed_sweeps': sum(
                r.get('sweep_status') == 'failed' for r in reactor_results),
            'elapsed_s': time.time() - t2,
            'equilibrium_check': {
                'within_tolerance': eq_result.get('within_tolerance'),
                'worst_abs_error': eq_result.get('worst_abs_error'),
            },
            'solids_scorecard': scorecard,
        }
        from pipeline.validation.viability import evaluate_turquoise
        viability = [evaluate_turquoise(r) for r in reactor_results]
        pipeline_state['phase2']['industrial_viability'] = {
            'pass': sum(v['status'] == 'pass' for v in viability),
            'fail': sum(v['status'] == 'fail' for v in viability),
            'unknown': sum(v['status'] == 'unknown' for v in viability),
        }
        # Same meaning as the orchestrator: the judge's single-pass solids X in
        # the headline band. A raw max across catalysts has no ranking
        # authority before B5 closure.
        pipeline_state['phase2']['best_conversion'] = scorecard.get(
            'headline_solids_conversion')
        pipeline_state['phase2']['best_conversion_scope'] = (
            'solids_single_pass_judge_'
            + str(scorecard.get('judge_catalyst') or 'unavailable'))
        if tea_results:
            best_tea = min(tea_results, key=lambda x: x['h2_cost_usd_kg'])
            pipeline_state['phase2']['best_h2_cost_usd_kg'] = best_tea['h2_cost_usd_kg']
    except (ImportError, RuntimeError) as e:
        # Imports are Cantera-guarded; a missing Cantera surfaces as a
        # RuntimeError from the equilibrium check, not an ImportError.
        print(f"  Phase 2 skipped (Cantera not available): {e}")
        pipeline_state['phase2'] = {'skipped': True, 'reason': str(e)}
    save_json(pipeline_state, "pipeline_state.json")


def _calibrate_ranker_from_dft(
    results_dir: Path, dft_candidates, converged_dft: list[dict]
) -> dict:
    """Calibrate TreeRanker with converged Phase 3 DFT validation residuals."""
    if not converged_dft:
        return {'n_validation_samples': 0}
    try:
        import pandas as pd
        from pipeline.screening.small_data_ranker import fit_tree_ranker
        from pipeline.search.adaptive_validation import record_validation
        from pipeline.search.discovery import candidate_id, discovery_region
        from pipeline.utils import save_json

        ev_path = results_dir / 'screening' / 'branch_ranker_evidence.csv'
        scan_db = results_dir / 'screening' / 'indexed_scan.sqlite'
        if not ev_path.is_file():
            return {'n_validation_samples': 0}
        ranker = fit_tree_ranker(pd.read_csv(ev_path), 'turquoise_hydrogen')
        baseline_bias = dict(ranker.class_bias)
        by_cid = {}
        for _, row in dft_candidates.iterrows():
            try:
                g = (
                    ast.literal_eval(row['genome'])
                    if isinstance(row['genome'], str)
                    else tuple(row['genome'])
                )
                by_cid[str(row.get('candidate_id') or candidate_id(g))] = row
            except Exception:
                pass
        recorded = 0
        observations = []
        for r in converged_dft:
            raw_genome = r.get('genome')
            genome = (
                ast.literal_eval(raw_genome)
                if isinstance(raw_genome, str)
                else tuple(raw_genome)
            )
            cid = str(r.get('candidate_id') or candidate_id(genome))
            src = by_cid.get(cid, {})
            raw_force = r.get('max_force_Ry_bohr')
            force = float(raw_force) if pd.notna(raw_force) else 0.0
            raw_base = src.get('E_act', 1.0) if hasattr(src, 'get') else 1.0
            base_e = float(raw_base) if pd.notna(raw_base) else 1.0
            raw_dft = r.get('E_act_dft')
            dft_e_act = (
                float(raw_dft)
                if pd.notna(raw_dft)
                else max(0.05, base_e + 0.10 * force)
            )
            provenance = {'source_id': f'screening_dft:{cid}'}
            record_validation(
                str(scan_db),
                'turquoise_hydrogen',
                genome,
                base_e,
                dft_e_act,
                'screening_dft',
                dft_e_act <= 1.0,
                provenance,
            )
            observations.append({
                'application': 'turquoise_hydrogen',
                'candidate_id': cid,
                'material_class': str(genome[0]),
                'region': '|'.join(discovery_region(genome)),
                'genome': repr(genome),
                'predicted_eV': base_e,
                'observed_eV': dft_e_act,
                'fidelity': 'screening_dft',
                'productive': bool(dft_e_act <= 1.0),
                'provenance': provenance,
            })
            recorded += 1
        class_bias = ranker.apply_validation_calibration(str(scan_db))
        class_delta = {
            cls: float(val - baseline_bias.get(cls, 0.0))
            for cls, val in class_bias.items()
        }
        calibration = {
            'n_validation_samples': recorded,
            'class_bias_eV': class_bias,
            'class_bias_delta_eV': class_delta,
            'region_bias_eV': dict(ranker.region_bias),
            'observations': observations,
        }
        save_json(calibration, 'active_learning_calibration.json', subdir='dft')
        return calibration
    except Exception as exc:
        return {'error': str(exc)[:200]}


def _run_phase3_dft(args, pipeline_state: dict, dft_candidates, results_dir: Path) -> None:
    """Queue-backed, concurrency-bounded screening DFT on the validation slate."""
    from pipeline.utils import print_banner, save_json

    print_banner("PHASE 3: DFT VALIDATION (Quantum ESPRESSO)")
    t3 = time.time()
    try:
        from pipeline.validation.dft_validator import validate_catalyst
        from pipeline.validation.task_queue import ValidationTaskQueue
        from pipeline.search.discovery import candidate_id

        represented_classes = (
            dft_candidates['material_class'].nunique()
            if 'material_class' in dft_candidates else 0)
        n_dft = min(max(14, represented_classes), len(dft_candidates))
        dft_tasks = []
        protocol_id = (
            f'screening-dft-v2:sssp-1.3.0:physical-realization-v2:'
            f'mpi={args.qe_mpi_ranks}:'
            f'omp={args.qe_omp_threads}')
        task_queue = ValidationTaskQueue(
            results_dir / 'dft' / 'validation_tasks.sqlite')
        # The campaign is the only writer: any 'running' row at start is
        # from a killed run and must be reclaimable immediately.
        task_queue.recover_stale(stale_after_s=0)
        for idx, (_, row) in enumerate(dft_candidates.head(n_dft).iterrows()):
            try:
                genome = ast.literal_eval(row['genome'])
                cid = candidate_id(genome)
                name = f"campaign_{cid}"
                task_queue.enqueue(
                    'turquoise_hydrogen', genome,
                    'screening_dft', protocol_id)
                dft_tasks.append((name, cid, genome))
            except Exception as e:
                print(f"    DFT input failed for cat_{idx}: {e}")
        from concurrent.futures import ThreadPoolExecutor

        from pipeline.validation.qe_workflows import QEExecutionConfig
        qe_exec = QEExecutionConfig(
            mpi_ranks=args.qe_mpi_ranks,
            omp_threads=args.qe_omp_threads,
        )

        def _validate_task(task):
            name, cid, genome = task
            if not task_queue.claim(
                    'turquoise_hydrogen', cid,
                    'screening_dft', protocol_id):
                state = task_queue.task(
                    'turquoise_hydrogen', cid, 'screening_dft', protocol_id)
                if (state and state['status'] == 'converged'
                        and state['result_path']
                        and Path(state['result_path']).is_file()):
                    # Already-converged work from an earlier run still counts.
                    return json.loads(Path(state['result_path']).read_text())
                return {'catalyst_name': name, 'candidate_id': cid,
                        'skipped': True, 'converged': False,
                        'reason': 'task_not_claimable'}
            try:
                result = validate_catalyst(
                    name, genome, run_dft=True, restart_incomplete=True,
                    execution=qe_exec)
                task_queue.finish(
                    'turquoise_hydrogen', cid, 'screening_dft',
                    protocol_id, bool(result.get('converged')),
                    result_path=str(results_dir / 'dft' / f'{name}_dft.json'),
                    error=result.get('error'))
                return result
            except Exception as exc:
                task_queue.finish(
                    'turquoise_hydrogen', cid, 'screening_dft',
                    protocol_id, False, error=str(exc))
                return {'catalyst_name': name, 'error': str(exc)[:200],
                        'converged': False}
        with ThreadPoolExecutor(
                max_workers=min(args.qe_max_concurrent, len(dft_tasks) or 1)
        ) as executor:
            dft_results = list(executor.map(_validate_task, dft_tasks))
        converged_dft = [r for r in dft_results if r.get('converged')]
        al_calibration = _calibrate_ranker_from_dft(
            results_dir, dft_candidates, converged_dft
        )
        pipeline_state['phase3'] = {
            'catalysts_validated': len(dft_results),
            'converged': len(converged_dft),
            'active_learning_calibration': al_calibration,
            'qe_max_concurrent': args.qe_max_concurrent,
            'qe_mpi_ranks': args.qe_mpi_ranks,
            'qe_omp_threads': args.qe_omp_threads,
            'task_queue': task_queue.summary(
                'turquoise_hydrogen', 'screening_dft'),
            'evidence_level': 'screening_dft',
            'elapsed_s': time.time() - t3,
        }
    except (ImportError, Exception) as e:
        print(f"  Phase 3 skipped: {e}")
        pipeline_state['phase3'] = {'skipped': True, 'reason': str(e)[:200]}
    save_json(pipeline_state, "pipeline_state.json")


def _candidate_hamiltonian_path(results_dir: Path, candidate_key: str) -> Path:
    """Conventional location of a candidate's provenance-bound FCIDUMP."""
    return results_dir / 'vqe' / 'hamiltonians' / candidate_key / 'FCIDUMP'


def _run_vqe(name: str, results_dir: Path, hamiltonian: Path) -> dict:
    """Run real CUDA-Q in quantum-env and accept only benchmarked evidence.

    CUDA-Q lives in quantum-env and is intentionally not imported from
    fairchem-env. The module's toy/mock paths are never campaign validation,
    so a candidate-specific FCIDUMP (with its JSON sidecar) is required.
    """
    code = (
        "import json; from pipeline.validation.vqe_transition_state import "
        "validate_transition_state; r=validate_transition_state("
        f"{name!r}, 'CH_split', target='nvidia', "
        f"candidate_hamiltonian={str(hamiltonian)!r}); print(json.dumps(r))"
    )
    proc = subprocess.run(
        ['conda', 'run', '-n', 'quantum-env', 'python', '-c', code],
        cwd=str(REPO_DIR), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"CUDA-Q failed for {name}: {proc.stderr[-1000:]}")
    result_path = results_dir / 'vqe' / f'vqe_{name}_CH_split.json'
    if not result_path.exists():
        raise RuntimeError(f"CUDA-Q produced no result for {name}")
    result = json.loads(result_path.read_text())
    if result.get('mock'):
        raise RuntimeError(f"CUDA-Q returned mock evidence for {name}")
    if not result.get('catalyst_specific_hamiltonian') or not result.get('benchmarked'):
        raise RuntimeError(
            f"CUDA-Q result for {name} is not catalyst-specific benchmarked evidence")
    return result


def _run_phase4_vqe(pipeline_state: dict, top_catalysts, results_dir: Path) -> None:
    """VQE transition-state refinement of the leading champions.

    Only champions with a candidate-specific Hamiltonian are run; the rest
    are listed as missing inputs instead of burning CUDA-Q time on a toy
    model that the evidence gate would reject anyway.
    """
    from pipeline.utils import print_banner, save_json

    print_banner("PHASE 4: VQE TRANSITION STATES (CUDA-Q)")
    t4 = time.time()
    try:
        n_vqe = min(MAX_VQE_CANDIDATES, len(top_catalysts))
        ready, missing = [], []
        for i, (_, row) in enumerate(top_catalysts.head(n_vqe).iterrows()):
            candidate_key = _reactor_candidate_key(row, i)
            hamiltonian = _candidate_hamiltonian_path(results_dir, candidate_key)
            (ready if hamiltonian.is_file() else missing).append(
                (candidate_key, hamiltonian))
        if not ready:
            pipeline_state['phase4'] = {
                'skipped': True,
                'reason': 'no candidate-specific Hamiltonians (FCIDUMP + sidecar)',
                'missing_hamiltonians': [str(path) for _, path in missing],
            }
        else:
            vqe_results = [_run_vqe(key, results_dir, path) for key, path in ready]
            pipeline_state['phase4'] = {
                'catalysts_refined': len(vqe_results),
                'candidates': [key for key, _ in ready],
                'missing_hamiltonians': [str(path) for _, path in missing],
                'elapsed_s': time.time() - t4,
            }
    except (ImportError, Exception) as e:
        print(f"  Phase 4 skipped: {e}")
        pipeline_state['phase4'] = {'skipped': True, 'reason': str(e)[:200]}
    save_json(pipeline_state, "pipeline_state.json")


def _model_pemfc_stack(pipeline_state: dict, top_fc) -> None:
    """Membrane sweep on the ORR slate, then a 400-cell stack for the best cell."""
    from pipeline.utils import print_banner
    from pipeline.fuel_cell.pemfc import sweep_membranes
    from pipeline.fuel_cell.stack import StackConfig, model_stack

    print_banner("PHASE 5B: PEMFC STACK MODELING")
    pemfc_results = []
    for _, row in top_fc.iterrows():
        name = str(row.get('candidate_id', ''))
        if not name or Path(name).name != name:
            raise ValueError('PEMFC candidate_id is missing or unsafe')
        eta = row.get('orr_overpotential_V', 0.4)
        mat_cls = row.get('material_class', None)
        mem = sweep_membranes(name, eta, material_class=mat_cls)
        pemfc_results.extend(mem)
    if not pemfc_results:
        return

    from pipeline.fuel_cell import composite_score

    best = max(pemfc_results, key=composite_score)
    stack = model_stack(StackConfig(n_cells=400,
        cell_voltage_V=best.get('peak_voltage_V', 0.65),
        current_density_A_cm2=best.get('peak_current_A_cm2', 1.5)))
    pipeline_state['phase5_stack'] = {
        'best_power_W_cm2': best.get('peak_power_W_cm2', 0),
        'best_efficiency': best.get('efficiency_at_peak', 0),
        'min_overpotential_V': best.get('orr_overpotential_V', 1.0),
        'best_catalyst': best.get('cathode_catalyst', 'unknown'),
        'best_membrane': best.get('membrane', 'unknown'),
        'stack_net_kW': stack.get('net_power_kW', 0),
        'stack_efficiency': stack.get('system_efficiency', 0),
    }
    from pipeline.validation.viability import evaluate_fuel_cell
    viability_record = dict(best)
    viability_record['system_efficiency'] = stack.get('system_efficiency', 0)
    pipeline_state['phase5_stack']['industrial_viability'] = evaluate_fuel_cell(viability_record)


def _run_phase5_fuel_cell(args, pipeline_state: dict, t_deadline: float,
                         results_dir: Path) -> None:
    """ORR branch-and-bound over the same space, then PEMFC stack modelling."""
    from pipeline.utils import print_banner, save_json

    print_banner("PHASE 5: FUEL CELL ORR — BRANCH-AND-BOUND (21.1B ENCODED SPACE)")
    t5 = time.time()

    remaining_hours = (t_deadline - time.time()) / 3600 if t_deadline != float('inf') else float('inf')
    time_str = f"{remaining_hours:.1f}h" if remaining_hours != float('inf') else "unlimited"
    print(f"  Remaining time: {time_str}")
    print("  FC search: catalyst-guided, coverage-safe branch-and-bound")
    print(f"  Same 21.1B encoded design space, ORR-specific objectives")

    from pipeline.screening.fc_genetic_optimizer import run_fc_branch_discovery, FCBranchDiscoveryConfig
    from pipeline.screening.stage_selection import (
        annotate_evidence, select_for_reactor, select_for_validation)

    fc_config = FCBranchDiscoveryConfig(
        **_branch_config_kwargs(args, t_deadline)
    )

    fc_pareto, fc_screening_db = run_fc_branch_discovery(fc_config)

    from pipeline.search.scope import scope_pemfc_pool
    fc_valid = fc_screening_db[fc_screening_db['valid'] == True].copy()
    fc_evidence = annotate_evidence(
        fc_screening_db, 'orr_overpotential_V')
    fc_pool, fc_scope = scope_pemfc_pool(fc_screening_db)
    top_fc = select_for_reactor(
        fc_pool, PEMFC_CANDIDATES, 'orr_overpotential_V',
        min_per_class=1)
    fc_validation = select_for_validation(
        fc_pool, min(args.validation_batch, len(fc_pool)),
        'orr_overpotential_V', min_per_class=args.min_validation_per_class)
    fc_validation_path = results_dir / 'fuel_cell' / 'validation_slate.csv'
    fc_validation_path.parent.mkdir(parents=True, exist_ok=True)
    fc_validation.to_csv(fc_validation_path, index=False)

    pipeline_state['phase5_branch'] = {
        'pareto_size': len(fc_pareto),
        'total_evaluated': len(fc_screening_db),
        'valid_count': len(fc_valid),
        'pemfc_model_count': len(top_fc),
        'validation_resolution_count': len(fc_validation),
        'validation_slate': str(fc_validation_path),
        'pemfc_scope': fc_scope,
        'candidate_dispositions': fc_evidence[
            'candidate_disposition'].value_counts().to_dict(),
        'elapsed_s': time.time() - t5,
    }
    if len(fc_valid) > 0 and 'orr_overpotential_V' in fc_valid.columns:
        pipeline_state['phase5_branch']['best_overpotential_V'] = float(fc_valid['orr_overpotential_V'].min())

    if time.time() < t_deadline and len(top_fc) > 0:
        _model_pemfc_stack(pipeline_state, top_fc)

    pipeline_state['phase5'] = {
        **pipeline_state.get('phase5_branch', {}),
        **pipeline_state.get('phase5_stack', {}),
        'n_cathodes_screened': len(fc_screening_db),
        'n_valid': len(fc_valid),
    }
    save_json(pipeline_state, "pipeline_state.json")
    print(f"\n  Phase 5: {time.time()-t5:.0f}s | FC catalysts evaluated: {len(fc_screening_db):,}")


def _run_phase6_report(pipeline_state: dict, t_start: float) -> dict:
    """Write the full report and print the dual-campaign summary."""
    from pipeline.utils import print_banner, save_json, load_json
    from pipeline.stages.report import run_report_stage

    print_banner("PHASE 6: REPORT")
    t6 = time.time()
    pipeline_state = load_json("pipeline_state.json") or pipeline_state
    pipeline_state['total_elapsed_s'] = time.time() - t_start
    outcome = run_report_stage(pipeline_state)
    pipeline_state['phase6'] = {**outcome.state, 'elapsed_s': time.time() - t6}
    pipeline_state['total_elapsed_s'] = time.time() - t_start
    save_json(pipeline_state, "pipeline_state.json")

    total = pipeline_state['total_elapsed_s']
    h2_eval = pipeline_state.get('phase1', {}).get('total_evaluated', 0)
    fc_eval = pipeline_state.get('phase5', {}).get(
        'total_evaluated',
        pipeline_state.get('phase5_branch', {}).get('total_evaluated', 0))
    print("\n" + "=" * 80)
    print(f"  DUAL-CAMPAIGN COMPLETE: {total:.0f}s ({total/3600:.2f} hours)")
    print(f"  H₂ Production catalysts evaluated: {h2_eval:,}")
    print(f"  Fuel Cell catalysts evaluated:      {fc_eval:,}")
    print(f"  Total catalysts screened:           {h2_eval + fc_eval:,}")
    print("=" * 80)
    return pipeline_state


def _assess_readiness(args, results_dir: Path) -> None:
    """Coverage, prior-art and six-point readiness; fail closed for final campaigns."""
    from pipeline.utils import save_json
    from pipeline.evidence.readiness import campaign_readiness

    h2_ready = campaign_readiness(
        results_dir / 'screening' / 'turquoise_hydrogen_coverage_certificate.json', args.prior_art_db,
        evidence_manifest=args.evidence_manifest if args.final_campaign else None,
        application='turquoise_hydrogen', pyrolysis_mode=args.mode)
    fc_ready = campaign_readiness(
        results_dir / 'fuel_cell' / 'coverage_certificate.json', args.prior_art_db,
        evidence_manifest=args.evidence_manifest if args.final_campaign else None,
        application='fuel_cell')
    readiness = {'turquoise_hydrogen': h2_ready, 'fuel_cell': fc_ready,
                 'ready': h2_ready['ready'] and fc_ready['ready']}
    from pipeline.evidence.campaign_status import assess_campaign
    readiness['six_point_status'] = assess_campaign(
        results_dir, pyrolysis_mode=args.mode)
    readiness['ready'] = readiness['ready'] and readiness['six_point_status']['ready']
    save_json(readiness, 'campaign_readiness.json')
    if args.final_campaign and not readiness['ready']:
        raise SystemExit(f"Final campaign readiness failed: {readiness}")


def main():
    """Run the module command-line entry point.
    """
    parser = _build_parser()
    args = parser.parse_args()
    results_dir = _resolve_paths(args)
    _validate_resources(parser, args)
    _, n_gpus, gpu_names = _configure_environment(args)
    # Final-campaign mode is fail-closed after both application searches have
    # had an opportunity to update their coverage certificates.
    _load_hf_token()
    _check_design_space(args)

    t_start = time.time()
    t_deadline = t_start + args.hours * 3600 if args.hours > 0 else float('inf')
    pipeline_state = _record_provenance(args, n_gpus, gpu_names)

    top_catalysts, dft_candidates = _run_phase1_discovery(args, pipeline_state, t_deadline)
    if time.time() < t_deadline and len(top_catalysts) > 0:
        _run_phase2_reactors(args, pipeline_state, top_catalysts)
    if time.time() < t_deadline and not args.no_dft:
        _run_phase3_dft(args, pipeline_state, dft_candidates, results_dir)
    elif args.no_dft:
        pipeline_state['phase3'] = {'skipped': True, 'reason': '--no-dft flag'}
    if time.time() < t_deadline and not args.no_vqe:
        _run_phase4_vqe(pipeline_state, top_catalysts, results_dir)
    elif args.no_vqe:
        pipeline_state['phase4'] = {'skipped': True, 'reason': '--no-vqe flag'}
    if time.time() < t_deadline:
        _run_phase5_fuel_cell(args, pipeline_state, t_deadline, results_dir)
    _run_phase6_report(pipeline_state, t_start)
    _assess_readiness(args, results_dir)


if __name__ == '__main__':
    main()
