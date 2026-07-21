#!/usr/bin/env python3
"""
GPU-SATURATED DETERMINISTIC DISCOVERY CAMPAIGN

The only production candidate-search strategy is persistent, exhaustive
branch-and-bound over all 14 classes in the 21.1B indexed encoded space.
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


def main():
    parser = argparse.ArgumentParser(
        description='GPU-Saturated Turquoise H₂ Catalyst Discovery v2'
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
    parser.add_argument('--mode', type=str, choices=['ntec', 'thermocatalytic'], default='ntec',
                        help='Pyrolysis screening mode (default: ntec)')
    parser.add_argument('--ntec-conditions-json', default='',
                        help=('JSON with measured operating conditions and paired '
                              'NTEC/control effect calibration'))
    parser.add_argument('--scan-batch-size', type=int, default=65536)
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
    parser.add_argument('--expected-space-size', type=int, default=21_092_645_031,
                        help='Fail if the indexed population denominator differs')
    parser.add_argument('--prior-art-db', default='results/prior_art.sqlite')
    parser.add_argument('--prior-art-csv', action='append', default=[],
                        help='CSV registry to import; repeat for multiple sources')
    parser.add_argument('--final-campaign', action='store_true',
                        help='Fail closed unless coverage and prior-art readiness requirements pass')
    parser.add_argument('--evidence-manifest', default='results/evidence_manifest.json',
                        help='Measured/validated evidence counts required by --final-campaign')
    args = parser.parse_args()

    # The campaign runs under fairchem-env, while Quantum ESPRESSO is installed
    # in qe-env. Resolve it explicitly so DFT does not depend on the caller's PATH.
    if not os.environ.get('PW_X'):
        envs_dir = Path(sys.executable).resolve().parents[2]
        qe_binary = envs_dir / 'qe-env' / 'bin' / 'pw.x'
        if qe_binary.is_file():
            os.environ['PW_X'] = str(qe_binary)

    from pipeline.evidence.prior_art import PriorArtRegistry
    prior_registry = PriorArtRegistry(args.prior_art_db)
    for prior_csv in args.prior_art_csv:
        prior_registry.import_csv(prior_csv)

    # ─── Environment ─────────────────────────────────────────────────────────
    os.environ['PYROLYSIS_MODE'] = args.mode
    if args.ntec_conditions_json:
        os.environ['NTEC_CONDITIONS_JSON'] = args.ntec_conditions_json
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
    print(f"  Search: deterministic branch-and-bound only")
    print(f"  Calibration probes: {args.calibration_probes:,}")
    print(f"  Validation batch: {args.validation_batch:,}")
    print(f"  Prior-art records: {prior_registry.count():,}")
    print("=" * 80)

    # Final-campaign mode is fail-closed after both application searches have
    # had an opportunity to update their coverage certificates.

    # ─── HuggingFace Token (for OC20 surface models) ─────────────────────────
    hf_token = os.environ.get('HF_TOKEN', '')
    token_file = Path(__file__).parent / '.hf_token'
    if not hf_token and token_file.exists():
        hf_token = token_file.read_text().strip()
    if not hf_token and sys.stdin.isatty():
        print("\n  ⚡ HuggingFace token needed for OC20 surface models (eSen/UMA)")
        print("    Get one at: https://huggingface.co/settings/tokens")
        hf_token = input("  Enter HF token (or press Enter to skip): ").strip()
    if hf_token:
        os.environ['HF_TOKEN'] = hf_token
        # Save for future runs
        if not token_file.exists():
            token_file.write_text(hf_token)
            token_file.chmod(0o600)
        print(f"  ✓ HuggingFace token loaded (surface models enabled)")
    else:
        print(f"  ⚠ No HF token — using MACE-MP-0 (bulk model) only")

    # ─── Design space ────────────────────────────────────────────────────────
    from pipeline.common.catalyst_spaces import estimate_design_space_size, ALL_MATERIAL_CLASSES
    from pipeline.common.utils import print_banner, save_json, load_json

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

    t_start = time.time()
    t_deadline = t_start + args.hours * 3600 if args.hours > 0 else float('inf')

    # ─── Campaign Provenance (reproducibility metadata) ──────────────────
    def _get_git_sha():
        try:
            return subprocess.check_output(
                ['git', 'rev-parse', '--short', 'HEAD'],
                cwd=str(Path(__file__).parent), stderr=subprocess.DEVNULL
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

    # ─── Phase 1: BRANCH DISCOVERY + META ESEN-SM ───────────────────────────
    print_banner("PHASE 1: DETERMINISTIC BRANCH-AND-BOUND + META ESEN-SM")
    t1 = time.time()

    from pipeline.screening.genetic_optimizer import run_branch_discovery, BranchDiscoveryConfig

    branch_config = BranchDiscoveryConfig(
        initial_fairchem_samples=args.calibration_probes,
        fairchem_eval_top_k=args.validation_batch,
        exhaustive_batch_size=args.scan_batch_size,
        branch_leaf_size=args.branch_leaf_size,
        branch_probe_count=args.branch_probes,
        branch_max_leaves=None if args.branch_max_leaves == 0 else args.branch_max_leaves,
        expected_space_size=args.expected_space_size,
        # One wall-clock budget covers the entire dual-application campaign.
        max_runtime_s=None if args.hours == 0 else max(0.0, t_deadline - time.time()),
        prior_art_db=args.prior_art_db,
        min_validation_per_class=args.min_validation_per_class,
        min_resolved_leaves_per_class=args.branch_class_floor,
        branch_exploration_interval=args.branch_exploration_interval,
        refresh_pending_priorities=args.branch_priority_refresh,
    )

    pareto_genomes, screening_db = run_branch_discovery(branch_config)

    valid_db = screening_db[screening_db['valid'] == True].copy()
    ranking_db = valid_db
    if 'E_act_censored' in ranking_db.columns:
        uncensored = ranking_db[ranking_db['E_act_censored'] != True]
        if len(uncensored):
            ranking_db = uncensored
    if 'E_act' in valid_db.columns:
        top_catalysts = ranking_db.nsmallest(args.top_k, 'E_act')
    else:
        top_catalysts = valid_db.head(args.top_k)

    pipeline_state['phase1'] = {
        'pareto_size': len(pareto_genomes),
        'total_evaluated': len(screening_db),
        'valid_count': len(valid_db),
        'top_catalysts_count': len(top_catalysts),
        'elapsed_s': time.time() - t1,
        'search_strategy': 'deterministic_branch_and_bound',
    }
    if len(valid_db) > 0 and 'E_act' in valid_db.columns:
        pipeline_state['phase1']['best_E_act'] = float(valid_db['E_act'].min())
    save_json(pipeline_state, "pipeline_state.json")
    print(f"\n  Phase 1: {time.time()-t1:.0f}s | Evaluated: {len(screening_db):,}")

    # ─── Phase 2: Cantera Reactor Simulation ─────────────────────────────────
    if time.time() < t_deadline and len(top_catalysts) > 0:
        print_banner("PHASE 2: CANTERA REACTOR SIMULATION")
        t2 = time.time()
        try:
            from pipeline.process.reactor_mechanisms import write_full_mechanism
            from pipeline.process.reactor_models import run_reactor_sweep

            reactor_temps = [773.15, 900.0, 1100.0, 1300.0]

            reactor_results = []
            n_reactor = min(20, len(top_catalysts))
            for i, (_, row) in enumerate(top_catalysts.head(n_reactor).iterrows()):
                e_act = row.get('E_act', 1.0)
                cat_name = f"catalyst_{i}"
                print(f"  Reactor sim {i+1}/{n_reactor}: E_act={e_act:.3f} eV")
                try:
                    # Generate Cantera YAML mechanism from E_act
                    mech_file = write_full_mechanism(cat_name, e_act)
                    sweep = run_reactor_sweep(cat_name, str(mech_file),
                                              temperatures=reactor_temps,
                                              catalyst_E_act_eV=e_act)
                    if any(result.get('mock') for result in sweep):
                        raise RuntimeError('mock reactor output is forbidden in production')
                    best_condition = max(sweep, key=lambda r: r.get('CH4_conversion', 0)) if sweep else {}
                    best_conv = best_condition.get('CH4_conversion', 0)
                    reactor_results.append({
                        'catalyst': cat_name,
                        'E_act': e_act,
                        'best_conversion': best_conv,
                        'n_conditions': len(sweep),
                        **{k: best_condition.get(k) for k in (
                            'temperature_K', 'H2_selectivity', 'CH4_conversion',
                            'deactivation_fraction_per_h', 'coke_fraction')},
                    })
                except Exception as e:
                    print(f"    Reactor error: {e}")

            # ── Pyrolysis TEA ($/kg H₂) ──────────────────────────────────
            tea_results = []
            for r in reactor_results:
                conv = r.get('best_conversion', 0)
                if conv > 0.01:
                    # Simplified TEA: natural gas + energy + capex
                    ng_cost = 3.50  # $/MMBtu natural gas
                    energy_input = 8.5 / conv  # kWh/kg_H2 (endothermic)
                    electricity_cost = 0.06  # $/kWh
                    capex_amortized = 0.50  # $/kg_H2 (capex over 20 yr)
                    carbon_credit = -0.80   # $/kg_H2 (solid C revenue)
                    h2_cost = (ng_cost * 0.05 / conv +
                               energy_input * electricity_cost +
                               capex_amortized + carbon_credit)
                    tea_results.append({
                        'catalyst': r['catalyst'],
                        'h2_cost_usd_kg': round(max(0.5, h2_cost), 2),
                        'conversion': conv,
                    })

            pipeline_state['phase2'] = {
                'catalysts_simulated': len(reactor_results),
                'elapsed_s': time.time() - t2,
            }
            from pipeline.validation.viability import evaluate_turquoise
            viability = [evaluate_turquoise(r) for r in reactor_results]
            pipeline_state['phase2']['industrial_viability'] = {
                'pass': sum(v['status'] == 'pass' for v in viability),
                'fail': sum(v['status'] == 'fail' for v in viability),
                'unknown': sum(v['status'] == 'unknown' for v in viability),
            }
            if reactor_results:
                pipeline_state['phase2']['best_conversion'] = max(
                    r.get('best_conversion', 0) for r in reactor_results)
            if tea_results:
                best_tea = min(tea_results, key=lambda x: x['h2_cost_usd_kg'])
                pipeline_state['phase2']['best_h2_cost_usd_kg'] = best_tea['h2_cost_usd_kg']
        except ImportError as e:
            print(f"  Phase 2 skipped (Cantera not available): {e}")
            pipeline_state['phase2'] = {'skipped': True, 'reason': str(e)}
        save_json(pipeline_state, "pipeline_state.json")

    # ─── Phase 3: DFT Validation ─────────────────────────────────────────────
    if time.time() < t_deadline and not args.no_dft:
        print_banner("PHASE 3: DFT VALIDATION (Quantum ESPRESSO)")
        t3 = time.time()
        try:
            from pipeline.validation.dft_validator import validate_catalyst

            n_dft = min(10, len(top_catalysts))
            dft_results = []
            for idx, (_, row) in enumerate(top_catalysts.head(n_dft).iterrows()):
                try:
                    genome = ast.literal_eval(row['genome'])
                    result = validate_catalyst(f"campaign_cat_{idx}", genome, run_dft=True)
                    dft_results.append(result)
                except Exception as e:
                    print(f"    DFT failed for cat_{idx}: {e}")
            pipeline_state['phase3'] = {
                'catalysts_validated': len(dft_results),
                'elapsed_s': time.time() - t3,
            }
        except (ImportError, Exception) as e:
            print(f"  Phase 3 skipped: {e}")
            pipeline_state['phase3'] = {'skipped': True, 'reason': str(e)[:200]}
        save_json(pipeline_state, "pipeline_state.json")
    elif args.no_dft:
        pipeline_state['phase3'] = {'skipped': True, 'reason': '--no-dft flag'}

    # ─── Phase 4: VQE Transition States ──────────────────────────────────────
    if time.time() < t_deadline and not args.no_vqe:
        print_banner("PHASE 4: VQE TRANSITION STATES (CUDA-Q)")
        t4 = time.time()
        try:
            n_vqe = min(5, len(top_catalysts))
            vqe_results = []
            for i in range(n_vqe):
                # CUDA-Q lives in quantum-env and is intentionally not imported
                # from fairchem-env. Run the real backend there; never accept the
                # module's mock fallback as campaign validation.
                name = f"champion_{i}"
                code = (
                    "import json; from pipeline.validation.vqe_transition_state import "
                    "validate_transition_state; r=validate_transition_state("
                    f"{name!r}, 'CH_split', target='nvidia'); print(json.dumps(r))"
                )
                proc = subprocess.run(
                    ['conda', 'run', '-n', 'quantum-env', 'python', '-c', code],
                    cwd=str(Path(__file__).parent), capture_output=True, text=True)
                if proc.returncode != 0:
                    raise RuntimeError(f"CUDA-Q failed for {name}: {proc.stderr[-1000:]}")
                result_path = Path('results/vqe') / f'vqe_{name}_CH_split.json'
                if not result_path.exists():
                    raise RuntimeError(f"CUDA-Q produced no result for {name}")
                result = json.loads(result_path.read_text())
                if result.get('mock'):
                    raise RuntimeError(f"CUDA-Q returned mock evidence for {name}")
                if not result.get('catalyst_specific_hamiltonian') or not result.get('benchmarked'):
                    raise RuntimeError(
                        f"CUDA-Q result for {name} is not catalyst-specific benchmarked evidence")
                vqe_results.append(result)
            pipeline_state['phase4'] = {
                'catalysts_refined': len(vqe_results),
                'elapsed_s': time.time() - t4,
            }
        except (ImportError, Exception) as e:
            print(f"  Phase 4 skipped: {e}")
            pipeline_state['phase4'] = {'skipped': True, 'reason': str(e)[:200]}
        save_json(pipeline_state, "pipeline_state.json")
    elif args.no_vqe:
        pipeline_state['phase4'] = {'skipped': True, 'reason': '--no-vqe flag'}

    # ─── Phase 5: Fuel Cell ORR branch search ────────────────────────────────
    if time.time() < t_deadline:
        print_banner("PHASE 5: FUEL CELL ORR — BRANCH-AND-BOUND (21.1B ENCODED SPACE)")
        t5 = time.time()

        remaining_hours = (t_deadline - time.time()) / 3600 if t_deadline != float('inf') else float('inf')
        # Allocate 60% of remaining time to FC screening, 40% to PEMFC/stack
        time_str = f"{remaining_hours:.1f}h" if remaining_hours != float('inf') else "unlimited"
        print(f"  Remaining time: {time_str}")
        print(f"  FC search: deterministic branch-and-bound")
        print(f"  Same 21.1B encoded design space, ORR-specific objectives")

        from pipeline.screening.fc_genetic_optimizer import run_fc_branch_discovery, FCBranchDiscoveryConfig

        fc_config = FCBranchDiscoveryConfig(
            initial_fairchem_samples=args.calibration_probes,
            fairchem_eval_top_k=args.validation_batch,
            exhaustive_batch_size=args.scan_batch_size,
            branch_leaf_size=args.branch_leaf_size,
            branch_probe_count=args.branch_probes,
            branch_max_leaves=None if args.branch_max_leaves == 0 else args.branch_max_leaves,
            expected_space_size=args.expected_space_size,
            max_runtime_s=None if args.hours == 0 else max(0.0, t_deadline - time.time()),
            prior_art_db=args.prior_art_db,
            min_validation_per_class=args.min_validation_per_class,
            min_resolved_leaves_per_class=args.branch_class_floor,
            branch_exploration_interval=args.branch_exploration_interval,
            refresh_pending_priorities=args.branch_priority_refresh,
        )

        fc_pareto, fc_screening_db = run_fc_branch_discovery(fc_config)

        fc_valid = fc_screening_db[fc_screening_db['valid'] == True].copy()
        if 'orr_overpotential_V' in fc_valid.columns:
            top_fc = fc_valid.nsmallest(30, 'orr_overpotential_V')
        else:
            top_fc = fc_valid.head(30)

        pipeline_state['phase5_branch'] = {
            'pareto_size': len(fc_pareto),
            'total_evaluated': len(fc_screening_db),
            'valid_count': len(fc_valid),
            'elapsed_s': time.time() - t5,
        }
        if len(fc_valid) > 0 and 'orr_overpotential_V' in fc_valid.columns:
            pipeline_state['phase5_branch']['best_overpotential_V'] = float(fc_valid['orr_overpotential_V'].min())

        # ─── PEMFC Stack Modeling on top ORR catalysts ────────────────────
        if time.time() < t_deadline and len(top_fc) > 0:
            print_banner("PHASE 5B: PEMFC STACK MODELING")
            from pipeline.process.pemfc_model import sweep_membranes
            from pipeline.process.fuel_cell_stack import StackConfig, model_stack

            pemfc_results = []
            for _, row in top_fc.iterrows():
                name = row.get('name', str(row.get('genome', ''))[:30])
                eta = row.get('orr_overpotential_V', 0.4)
                mat_cls = row.get('material_class', None)
                mem = sweep_membranes(name, eta, material_class=mat_cls)
                pemfc_results.extend(mem)

            if pemfc_results:
                # Optimize for highest efficiency and least overvoltage first and foremost,
                # while maximizing peak power output as much as possible.
                def fc_composite_score(r):
                    eff = r.get('efficiency_at_peak', 0.0)
                    power = r.get('peak_power_W_cm2', 0.0)
                    eta = max(r.get('orr_overpotential_V', 0.4), 0.01)
                    return (eff * power) / eta

                best = max(pemfc_results, key=fc_composite_score)
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

        save_json(pipeline_state, "pipeline_state.json")
        print(f"\n  Phase 5: {time.time()-t5:.0f}s | FC catalysts evaluated: {len(fc_screening_db):,}")

    # ─── Phase 6: Report ─────────────────────────────────────────────────────
    print_banner("PHASE 6: REPORT")
    from pipeline.evidence.report_generator import generate_full_report
    pipeline_state = load_json("pipeline_state.json") or pipeline_state
    generate_full_report(pipeline_state)
    save_json(pipeline_state, "pipeline_state.json")

    total = time.time() - t_start
    h2_eval = pipeline_state.get('phase1', {}).get('total_evaluated', 0)
    fc_eval = pipeline_state.get('phase5_branch', {}).get('total_evaluated', 0)
    print("\n" + "=" * 80)
    print(f"  DUAL-CAMPAIGN COMPLETE: {total:.0f}s ({total/3600:.2f} hours)")
    print(f"  H₂ Production catalysts evaluated: {h2_eval:,}")
    print(f"  Fuel Cell catalysts evaluated:      {fc_eval:,}")
    print(f"  Total catalysts screened:           {h2_eval + fc_eval:,}")
    print("=" * 80)

    from pipeline.evidence.readiness import campaign_readiness
    h2_ready = campaign_readiness(
        'results/screening/turquoise_hydrogen_coverage_certificate.json', args.prior_art_db,
        evidence_manifest=args.evidence_manifest if args.final_campaign else None,
        application='turquoise_hydrogen', pyrolysis_mode=args.mode)
    fc_ready = campaign_readiness(
        'results/fuel_cell/coverage_certificate.json', args.prior_art_db,
        evidence_manifest=args.evidence_manifest if args.final_campaign else None,
        application='fuel_cell')
    readiness = {'turquoise_hydrogen': h2_ready, 'fuel_cell': fc_ready,
                 'ready': h2_ready['ready'] and fc_ready['ready']}
    from pipeline.evidence.campaign_status import assess_campaign
    readiness['six_point_status'] = assess_campaign('results', pyrolysis_mode=args.mode)
    readiness['ready'] = readiness['ready'] and readiness['six_point_status']['ready']
    save_json(readiness, 'campaign_readiness.json')
    if args.final_campaign and not readiness['ready']:
        raise SystemExit(f"Final campaign readiness failed: {readiness}")


if __name__ == '__main__':
    main()
