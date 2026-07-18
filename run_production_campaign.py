#!/usr/bin/env python3
"""
GPU-SATURATED PRODUCTION CAMPAIGN v2

Key changes from v1:
  - Smaller GA population (1000) for faster NSGA-II O(N²) sort
  - MACE validation every 5 generations (not 50) to keep GPUs hot
  - Larger MACE batch (500 per round) for better GPU saturation
  - 3000 generations × 500 MACE/round = 300,000 MACE evaluations
  - All 10 material classes with 25.3B design space
"""

import os
import sys
import time
import json
import argparse
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def main():
    parser = argparse.ArgumentParser(
        description='GPU-Saturated Turquoise H₂ Catalyst Discovery v2'
    )
    parser.add_argument('--pop', type=int, default=1000,
                        help='GA population (default: 1000)')
    parser.add_argument('--gens', type=int, default=3000,
                        help='Total generations (default: 3000)')
    parser.add_argument('--mace-batch', type=int, default=500,
                        help='Initial MACE batch (default: 500)')
    parser.add_argument('--mace-per-round', type=int, default=500,
                        help='MACE evaluations per validation round (default: 500)')
    parser.add_argument('--mace-interval', type=int, default=5,
                        help='Generations between MACE rounds (default: 5)')
    parser.add_argument('--hours', type=float, default=48.0,
                        help='Max wall-clock hours (default: 48)')
    parser.add_argument('--top-k', type=int, default=200,
                        help='Top-K for reactor simulation (default: 200)')
    parser.add_argument('--no-dft', action='store_true')
    parser.add_argument('--no-vqe', action='store_true')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    # ─── Environment ─────────────────────────────────────────────────────────
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
    print(f"  Pop: {args.pop} | Gens: {args.gens}")
    print(f"  MACE: {args.mace_per_round}/round every {args.mace_interval} gens")
    total_mace = args.mace_batch + args.mace_per_round * (args.gens // args.mace_interval)
    print(f"  Estimated total MACE evaluations: {total_mace:,}")
    print(f"  Estimated surrogate evaluations: {args.pop * args.gens:,}")
    print("=" * 80)

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
    from pipeline.catalyst_spaces import estimate_design_space_size, ALL_MATERIAL_CLASSES
    from pipeline.utils import print_banner, save_json, load_json

    sizes = estimate_design_space_size()
    print(f"\n  Design Space: {sizes['TOTAL']:,} ({sizes['TOTAL']/1e9:.1f}B)")
    for cls in ALL_MATERIAL_CLASSES:
        print(f"    {cls:20s}: {sizes[cls]:>15,}")

    t_start = time.time()
    t_deadline = t_start + args.hours * 3600
    pipeline_state = {}

    # ─── Phase 1: GA + META ESEN-SM ──────────────────────────────────────────
    print_banner("PHASE 1: GPU-ACCELERATED META ESEN-SM + NSGA-II")
    t1 = time.time()

    from pipeline.genetic_optimizer import run_genetic_algorithm, GAConfig

    ga_config = GAConfig(
        pop_size=args.pop,
        n_generations=args.gens,
        initial_mace_samples=args.mace_batch,
        mace_eval_interval=args.mace_interval,
        mace_eval_top_k=args.mace_per_round,
        surrogate_retrain_interval=args.mace_interval,
        mutation_rate=0.35,
        crossover_rate=0.7,
        seed=args.seed,
    )

    pareto_genomes, screening_db = run_genetic_algorithm(ga_config)

    valid_db = screening_db[screening_db['valid'] == True].copy()
    if 'E_act' in valid_db.columns:
        top_catalysts = valid_db.nsmallest(args.top_k, 'E_act')
    else:
        top_catalysts = valid_db.head(args.top_k)

    pipeline_state['phase1'] = {
        'pareto_size': len(pareto_genomes),
        'total_evaluated': len(screening_db),
        'valid_count': len(valid_db),
        'top_catalysts_count': len(top_catalysts),
        'elapsed_s': time.time() - t1,
        'total_mace_target': total_mace,
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
            from pipeline.reactor_mechanisms import write_full_mechanism
            from pipeline.reactor_models import run_reactor_sweep

            reactor_results = []
            n_reactor = min(20, len(top_catalysts))
            for i, (_, row) in enumerate(top_catalysts.head(n_reactor).iterrows()):
                e_act = row.get('E_act', 1.0)
                cat_name = f"catalyst_{i}"
                print(f"  Reactor sim {i+1}/{n_reactor}: E_act={e_act:.3f} eV")
                try:
                    # Generate Cantera YAML mechanism from E_act
                    mech_file = write_full_mechanism(cat_name, e_act)
                    sweep = run_reactor_sweep(cat_name, str(mech_file))
                    best_conv = max(r.get('CH4_conversion', 0) for r in sweep) if sweep else 0
                    reactor_results.append({
                        'catalyst': cat_name,
                        'E_act': e_act,
                        'best_conversion': best_conv,
                        'n_conditions': len(sweep),
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
            from pipeline.dft_validator import run_dft_validation

            n_dft = min(10, len(top_catalysts))
            dft_results = run_dft_validation(top_catalysts.head(n_dft))
            pipeline_state['phase3'] = {
                'catalysts_validated': len(dft_results) if dft_results else 0,
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
            from pipeline.vqe_transition_state import run_vqe_refinement

            n_vqe = min(5, len(top_catalysts))
            vqe_results = run_vqe_refinement(top_catalysts.head(n_vqe))
            pipeline_state['phase4'] = {
                'catalysts_refined': len(vqe_results) if vqe_results else 0,
                'elapsed_s': time.time() - t4,
            }
        except (ImportError, Exception) as e:
            print(f"  Phase 4 skipped: {e}")
            pipeline_state['phase4'] = {'skipped': True, 'reason': str(e)[:200]}
        save_json(pipeline_state, "pipeline_state.json")
    elif args.no_vqe:
        pipeline_state['phase4'] = {'skipped': True, 'reason': '--no-vqe flag'}

    # ─── Phase 5: Fuel Cell ORR GA (SAME 25.3B DESIGN SPACE) ─────────────────
    if time.time() < t_deadline:
        print_banner("PHASE 5: FUEL CELL ORR — GA + META ESEN-SM (25.3B DESIGN SPACE)")
        t5 = time.time()

        remaining_hours = (t_deadline - time.time()) / 3600
        # Allocate 60% of remaining time to FC screening, 40% to PEMFC/stack
        fc_gens = max(100, int(args.gens * 0.5))  # Half the gens of methane
        print(f"  Remaining time: {remaining_hours:.1f}h")
        print(f"  FC-GA: {fc_gens} generations, pop={args.pop}")
        print(f"  Same 25.3B design space, ORR-specific objectives")

        from pipeline.fc_genetic_optimizer import run_fc_genetic_algorithm, FCGAConfig

        fc_config = FCGAConfig(
            pop_size=args.pop,
            n_generations=fc_gens,
            initial_mace_samples=args.mace_batch,
            mace_eval_interval=args.mace_interval,
            mace_eval_top_k=args.mace_per_round,
            surrogate_retrain_interval=args.mace_interval,
            mutation_rate=0.35,
            crossover_rate=0.7,
            seed=args.seed + 1000,  # different seed for diversity
        )

        fc_pareto, fc_screening_db = run_fc_genetic_algorithm(fc_config)

        fc_valid = fc_screening_db[fc_screening_db['valid'] == True].copy()
        if 'orr_overpotential_V' in fc_valid.columns:
            top_fc = fc_valid.nsmallest(30, 'orr_overpotential_V')
        else:
            top_fc = fc_valid.head(30)

        pipeline_state['phase5_ga'] = {
            'pareto_size': len(fc_pareto),
            'total_evaluated': len(fc_screening_db),
            'valid_count': len(fc_valid),
            'elapsed_s': time.time() - t5,
        }
        if len(fc_valid) > 0 and 'orr_overpotential_V' in fc_valid.columns:
            pipeline_state['phase5_ga']['best_overpotential_V'] = float(fc_valid['orr_overpotential_V'].min())

        # ─── PEMFC Stack Modeling on top ORR catalysts ────────────────────
        if time.time() < t_deadline and len(top_fc) > 0:
            print_banner("PHASE 5B: PEMFC STACK MODELING")
            from pipeline.pemfc_model import sweep_membranes
            from pipeline.fuel_cell_stack import StackConfig, model_stack

            pemfc_results = []
            for _, row in top_fc.iterrows():
                name = row.get('name', str(row.get('genome', ''))[:30])
                eta = row.get('orr_overpotential_V', 0.4)
                mem = sweep_membranes(name, eta)
                pemfc_results.extend(mem)

            if pemfc_results:
                # Optimize for highest efficiency and least overvoltage first and foremost,
                # while maximizing peak power output as much as possible.
                def fc_composite_score(r):
                    eff = r.get('efficiency_at_rated', 0.0)
                    power = r.get('peak_power_W_cm2', 0.0)
                    eta = max(r.get('orr_overpotential_V', 0.4), 0.01)
                    return (eff * power) / eta

                best = max(pemfc_results, key=fc_composite_score)
                stack = model_stack(StackConfig(n_cells=400,
                    cell_voltage_V=best.get('peak_voltage_V', 0.65),
                    current_density_A_cm2=best.get('peak_current_A_cm2', 1.5)))
                pipeline_state['phase5_stack'] = {
                    'best_power_W_cm2': best.get('peak_power_W_cm2', 0),
                    'best_efficiency': best.get('efficiency_at_rated', 0),
                    'min_overpotential_V': best.get('orr_overpotential_V', 1.0),
                    'best_catalyst': best.get('catalyst', 'unknown'),
                    'best_membrane': best.get('membrane', 'unknown'),
                    'stack_net_kW': stack.get('net_power_kW', 0),
                    'stack_efficiency': stack.get('efficiency', 0),
                }

        save_json(pipeline_state, "pipeline_state.json")
        print(f"\n  Phase 5: {time.time()-t5:.0f}s | FC catalysts evaluated: {len(fc_screening_db):,}")

    # ─── Phase 6: Report ─────────────────────────────────────────────────────
    print_banner("PHASE 6: REPORT")
    from pipeline.report_generator import generate_full_report
    pipeline_state = load_json("pipeline_state.json") or pipeline_state
    generate_full_report(pipeline_state)
    save_json(pipeline_state, "pipeline_state.json")

    total = time.time() - t_start
    h2_eval = pipeline_state.get('phase1', {}).get('total_evaluated', 0)
    fc_eval = pipeline_state.get('phase5_ga', {}).get('total_evaluated', 0)
    print("\n" + "=" * 80)
    print(f"  DUAL-CAMPAIGN COMPLETE: {total:.0f}s ({total/3600:.2f} hours)")
    print(f"  H₂ Production catalysts evaluated: {h2_eval:,}")
    print(f"  Fuel Cell catalysts evaluated:      {fc_eval:,}")
    print(f"  Total catalysts screened:           {h2_eval + fc_eval:,}")
    print("=" * 80)


if __name__ == '__main__':
    main()
