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

    # ─── Phase 1: GA + MACE ──────────────────────────────────────────────────
    print_banner("PHASE 1: GPU-ACCELERATED MACE + NSGA-II")
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

    # ─── Phase 5: Fuel Cell ORR GA (SAME 25.3B DESIGN SPACE) ─────────────────
    if time.time() < t_deadline:
        print_banner("PHASE 5: FUEL CELL ORR — GA + MACE (25.3B DESIGN SPACE)")
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
                best = max(pemfc_results, key=lambda r: r.get('peak_power_W_cm2', 0))
                stack = model_stack(StackConfig(n_cells=400,
                    cell_voltage_V=best.get('peak_voltage_V', 0.65),
                    current_density_A_cm2=best.get('peak_current_A_cm2', 1.5)))
                pipeline_state['phase5_stack'] = {
                    'best_power_W_cm2': best.get('peak_power_W_cm2', 0),
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
