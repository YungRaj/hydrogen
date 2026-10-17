"""Argument parsing, environment setup, and provenancehelpers for production campaigns."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

from pipeline.data_models.campaigns import (
    CampaignReadinessReport,
    PipelineStateDocument,
    validate_pipeline_state,
)
from pipeline.reactors.modes import DEFAULT_MODE, MODE_CHOICES
from typing import cast

REPO_DIR = Path(__file__).resolve().parents[2]


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


def _build_parser() -> argparse.ArgumentParser:
    """Return the campaign command-line interface."""
    parser = argparse.ArgumentParser(
        description='GPU-Saturated Turquoise H₂ Catalyst Discovery'
    )
    parser.add_argument(
        '--calibration-probes',
        type=int,
        default=500,
        help='Deterministic binary-tree probes used to calibrate the surrogate',
    )
    parser.add_argument(
        '--validation-batch',
        type=int,
        default=500,
        help='Branch/archive champions sent to the atomistic model',
    )
    parser.add_argument(
        '--min-validation-per-class',
        type=int,
        default=2,
        help='Reserved atomistic validation quota per represented class',
    )
    parser.add_argument(
        '--hours',
        type=float,
        default=0,
        help='Max wall-clock hours. 0 = unlimited (default: 0)',
    )
    parser.add_argument(
        '--top-k',
        type=int,
        default=200,
        help='Top-K for reactor simulation (default: 200)',
    )
    parser.add_argument(
        '--no-dft',
        dest='no_dft',
        action='store_true',
        default=True,
        help='Skip Phase 3 DFT validation (default: True)',
    )
    parser.add_argument(
        '--run-dft',
        dest='no_dft',
        action='store_false',
        help='Opt in to Phase 3 Quantum ESPRESSO DFT validation',
    )
    parser.add_argument(
        '--no-vqe',
        dest='no_vqe',
        action='store_true',
        default=True,
        help='Skip Phase 4 CUDA-Q VQE refinement (default: True)',
    )
    parser.add_argument(
        '--run-vqe',
        dest='no_vqe',
        action='store_false',
        help='Opt in to Phase 4 CUDA-Q VQE transition-state refinement',
    )
    parser.add_argument(
        '--mode',
        choices=MODE_CHOICES,
        default=DEFAULT_MODE,
        help=f'Methane-conversion pathway (default: {DEFAULT_MODE})',
    )
    parser.add_argument(
        '--ntec-conditions-json',
        default='',
        help=(
            'JSON with measured operating conditions and paired '
            'NTEC/control effect calibration'
        ),
    )
    parser.add_argument(
        '--electrochemical-conditions-json',
        default='',
        help=(
            'JSON measured electrochemical operating point; '
            'electrolyte_phase is aqueous or molten'
        ),
    )
    parser.add_argument(
        '--workers-per-gpu',
        type=int,
        default=32,
        help=(
            'Candidate threads sharing each GPU in the batched '
            'eSEN engine (relaxation CPU work overlaps inference)'
        ),
    )
    parser.add_argument('--scan-batch-size', type=int, default=65536)
    parser.add_argument(
        '--scan-workers',
        type=int,
        default=8,
        help='Independent deterministic scanner shards per resolved leaf',
    )
    parser.add_argument('--branch-leaf-size', type=int, default=1_000_000)
    parser.add_argument('--branch-probes', type=int, default=9)
    parser.add_argument(
        '--branch-max-leaves',
        type=int,
        default=0,
        help='Process at most this many leaves per run; 0 continues until complete',
    )
    parser.add_argument(
        '--branch-class-floor',
        type=int,
        default=1,
        help='Resolve at least this many leaves per class before pure exploitation',
    )
    parser.add_argument(
        '--branch-exploration-interval',
        type=int,
        default=4,
        help='After class floors, reserve every Nth leaf for the least-covered class',
    )
    parser.add_argument(
        '--branch-priority-refresh',
        type=int,
        default=10000,
        help='Re-score up to this many pending nodes when resuming',
    )
    parser.add_argument(
        '--expected-space-size',
        type=int,
        default=21_092_657_271,
        help='Fail if the indexed population denominator differs',
    )
    parser.add_argument(
        '--results-dir',
        default='results',
        help='Isolated output root for this campaign (default: results)',
    )
    parser.add_argument(
        '--mechanisms-dir',
        default=None,
        help='Generated Cantera mechanisms (default: <results-dir>/mechanisms)',
    )
    parser.add_argument(
        '--prior-art-db',
        default=None,
        help='Prior-art SQLite database (default: <results-dir>/prior_art.sqlite)',
    )
    parser.add_argument(
        '--prior-art-csv',
        action='append',
        default=[],
        help='CSV registry to import; repeat for multiple sources',
    )
    parser.add_argument(
        '--calibration-json',
        default=None,
        help=(
            'Active-learning DFT calibration JSON '
            '(default: <results-dir>/dft/active_learning_calibration.json)'
        ),
    )
    parser.add_argument(
        '--final-campaign',
        action='store_true',
        help='Fail closed unless coverage and prior-art readiness requirements pass',
    )
    parser.add_argument(
        '--sota-campaign',
        action='store_true',
        help=(
            'Configure high-volume SOTA discovery preset (56+ sub-class leaves, '
            '4 leaves/class floor, 1000 calibration probes, 500 validation batch)'
        ),
    )
    parser.add_argument(
        '--evidence-manifest',
        default=None,
        help='Measured/validated evidence counts required by --final-campaign',
    )
    parser.add_argument(
        '--kinetics-validation-dir',
        default=None,
        help=(
            'Directory containing <candidate_id>/pyrolysis_validation.json; '
            'complete matching campaigns replace Cantera template barriers'
        ),
    )
    parser.add_argument(
        '--multiphysics-results-dir',
        default=None,
        help=(
            'Required OpenFOAM/MOOSE artifacts (default: '
            '<results-dir>/multiphysics)'
        ),
    )
    parser.add_argument('--qe-mpi-ranks', type=int, default=4)
    parser.add_argument('--qe-omp-threads', type=int, default=1)
    default_qe_concurrency = max(1, min(4, (os.cpu_count() or 1) // 4))
    parser.add_argument(
        '--qe-max-concurrent',
        type=int,
        default=default_qe_concurrency,
        help='Maximum independent candidate DFT jobs launched together',
    )
    return parser


def _apply_sota_preset(args: argparse.Namespace) -> None:
    """Apply high-volume SOTA branch-and-bound discovery defaults when requested."""
    if not getattr(args, 'sota_campaign', False):
        return
    args.branch_class_floor = max(int(args.branch_class_floor), 4)
    if int(args.branch_max_leaves) == 0:
        args.branch_max_leaves = 56
    args.calibration_probes = max(int(args.calibration_probes), 1000)
    args.validation_batch = max(int(args.validation_batch), 500)
    args.min_validation_per_class = max(int(args.min_validation_per_class), 4)


def _resolve_paths(args: argparse.Namespace) -> Path:
    """Make every output path absolute and default it under the results root."""
    _apply_sota_preset(args)
    results_dir = Path(args.results_dir).expanduser().resolve()
    args.results_dir = str(results_dir)
    mechanisms_dir = (
        Path(args.mechanisms_dir).expanduser().resolve()
        if args.mechanisms_dir
        else results_dir / 'mechanisms'
    )
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
        args.kinetics_validation_dir = str(
            results_dir / 'dft' / 'pyrolysis_kinetics'
        )
    if args.multiphysics_results_dir is None:
        args.multiphysics_results_dir = str(results_dir / 'multiphysics')
    else:
        args.multiphysics_results_dir = str(
            Path(args.multiphysics_results_dir).expanduser().resolve()
        )
    os.environ['HYDROGEN_RESULTS_DIR'] = str(results_dir)
    os.environ['HYDROGEN_MECHANISMS_DIR'] = str(mechanisms_dir)
    return results_dir


def _validate_resources(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Reject budgets the branch search or the QE allocation cannot honour."""
    if (
        args.scan_workers < 1
        or args.workers_per_gpu < 1
        or args.qe_mpi_ranks < 1
        or args.qe_omp_threads < 1
        or args.qe_max_concurrent < 1
    ):
        parser.error('scanner and QE resource dimensions must be positive')
    if args.calibration_probes < 20:
        parser.error('--calibration-probes must be at least 20 for the tree ranker')
    required_validation = 14 * args.min_validation_per_class
    if args.validation_batch < required_validation:
        parser.error(
            f'--validation-batch must be at least {required_validation} to reserve '
            f'{args.min_validation_per_class} candidate(s) across all 14 classes'
        )
    if not args.no_dft:
        requested_qe_cpus = (
            args.qe_mpi_ranks * args.qe_omp_threads * args.qe_max_concurrent
        )
        available_cpus = os.cpu_count() or 1
        if requested_qe_cpus > available_cpus:
            parser.error(
                f'QE allocation requests {requested_qe_cpus} CPUs, '
                f'but only {available_cpus} are visible'
            )


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
        os.environ['ELECTROCHEMICAL_CONDITIONS_JSON'] = (
            args.electrochemical_conditions_json
        )
    import torch

    n_gpus = torch.cuda.device_count()
    gpu_names = [torch.cuda.get_device_name(i) for i in range(n_gpus)]
    gpu_mem = [
        torch.cuda.get_device_properties(i).total_memory / 1e9
        for i in range(n_gpus)
    ]

    print('=' * 80)
    print('  GPU-SATURATED CATALYST DISCOVERY CAMPAIGN v2')
    print('=' * 80)
    for i in range(n_gpus):
        print(f'  GPU[{i}] {gpu_names[i]} — {gpu_mem[i]:.1f} GB')
    print(f'  CPUs: {os.cpu_count()} cores')
    print(f'  Mode: {args.mode.upper()}')
    print('  Search: catalyst-guided, coverage-safe branch-and-bound')
    print(f'  Calibration probes: {args.calibration_probes:,}')
    print(f'  Validation batch: {args.validation_batch:,}')
    print(f'  Prior-art records: {prior_registry.count():,}')
    print('=' * 80)
    return prior_registry, n_gpus, gpu_names


def _load_hf_token() -> None:
    """Load the HuggingFace token needed for the gated eSEN surface model.

    Uses the environment, ``$HF_HOME/token`` or the repo ``.hf_token``; prompts
    only on an interactive terminal; and fails fast otherwise, because
    screening has no non-gated fallback model.
    """
    from pipeline.screening.surface_calculator import _ensure_hf_token

    if not _ensure_hf_token() and sys.stdin.isatty():
        print('\n  HuggingFace token needed for the eSEN OC25 surface model')
        print('    Get one at: https://huggingface.co/settings/tokens')
        token = input('  Enter HF token: ').strip()
        if token:
            os.environ['HF_TOKEN'] = token
    if not os.environ.get('HF_TOKEN'):
        raise SystemExit(
            'No HuggingFace token: set HF_TOKEN or run `huggingface-cli login` '
            '(writes $HF_HOME/token); eSEN screening cannot run without it'
        )
    print('  ✓ HuggingFace token loaded (eSEN surface model enabled)')


def _check_design_space(args: argparse.Namespace) -> None:
    """Fail unless the indexed denominator and the canonical-identity audit pass."""
    from pipeline.evidence.design_space_audit import audit_design_space
    from pipeline.search.design_space import (
        ALL_MATERIAL_CLASSES,
        estimate_design_space_size,
    )
    from pipeline.utils import save_json

    sizes = estimate_design_space_size()
    if args.expected_space_size != sizes['TOTAL']:
        raise SystemExit(
            f'Population denominator mismatch: --expected-space-size={args.expected_space_size:,}, '
            f"but this repository indexes {sizes['TOTAL']:,}. Update the design-space "
            'definition or use the verified denominator; do not label it 25.3B.'
        )
    design_audit = audit_design_space(sample_per_class=2048)
    save_json(design_audit, 'design_space_audit.json', subdir='')
    if not design_audit['valid']:
        raise SystemExit(f"Design-space audit failed: {design_audit['failures']}")
    print(f"\n  Design Space: {sizes['TOTAL']:,} ({sizes['TOTAL']/1e9:.1f}B)")
    print(f"  Canonical identities: {design_audit['canonical_total']:,}")
    for cls in ALL_MATERIAL_CLASSES:
        print(f'    {cls:20s}: {sizes[cls]:>15,}')


def _get_git_sha() -> str:
    try:
        return (
            subprocess.check_output(
                ['git', 'rev-parse', '--short', 'HEAD'],
                cwd=str(REPO_DIR),
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
    except Exception:
        return 'unknown'


def _get_conda_env() -> str:
    try:
        return subprocess.check_output(
            ['conda', 'list', '--export'], stderr=subprocess.DEVNULL
        ).decode()[:2000]
    except Exception:
        return 'unavailable'


def _record_provenance(
    args: argparse.Namespace, n_gpus: int, gpu_names: list[str]
) -> PipelineStateDocument:
    """Start ``pipeline_state.json`` with reproducibility metadata."""
    from pipeline.utils import save_json

    pipeline_state = validate_pipeline_state(
        {
            'provenance': {
                'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                'git_sha': _get_git_sha(),
                'cli_args': vars(args),
                'n_gpus': n_gpus,
                'gpu_names': gpu_names,
                'conda_env_snapshot': _get_conda_env(),
            }
        }
    )
    save_json(pipeline_state, 'pipeline_state.json')
    return pipeline_state


def _assess_readiness(
    args: argparse.Namespace, results_dir: Path
) -> CampaignReadinessReport:
    """Coverage, prior-art and six-point readiness; fail closed for final campaigns."""
    from pipeline.evidence.campaign_status import assess_campaign
    from pipeline.evidence.readiness import campaign_readiness
    from pipeline.utils import save_json

    h2_ready = campaign_readiness(
        results_dir / 'screening' / 'turquoise_hydrogen_coverage_certificate.json',
        args.prior_art_db,
        evidence_manifest=args.evidence_manifest if args.final_campaign else None,
        application='turquoise_hydrogen',
        pyrolysis_mode=args.mode,
    )
    fc_ready = campaign_readiness(
        results_dir / 'fuel_cell' / 'coverage_certificate.json',
        args.prior_art_db,
        evidence_manifest=args.evidence_manifest if args.final_campaign else None,
        application='fuel_cell',
    )
    six_point = assess_campaign(results_dir, pyrolysis_mode=args.mode)
    readiness: CampaignReadinessReport = {
        'turquoise_hydrogen': cast(dict, h2_ready),  # type: ignore[typeddict-item]
        'fuel_cell': cast(dict, fc_ready),  # type: ignore[typeddict-item]
        'six_point_status': cast(dict, six_point),  # type: ignore[typeddict-item]
        'ready': bool(
            h2_ready['ready'] and fc_ready['ready'] and six_point['ready']
        ),
    }
    save_json(readiness, 'campaign_readiness.json')
    if args.final_campaign and not readiness['ready']:
        raise SystemExit(f'Final campaign readiness failed: {readiness}')
    return readiness
