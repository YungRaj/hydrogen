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

    # ─── Phase 5: Fuel Cell ──────────────────────────────────────────────────
    if time.time() < t_deadline:
        print_banner("PHASE 5: FUEL CELL + PEMFC STACK")
        from pipeline.fuel_cell_cathode_screener import run_cathode_screening
        from pipeline.pemfc_model import sweep_membranes
        from pipeline.fuel_cell_stack import StackConfig, model_stack
        t5 = time.time()

        cathode_df = run_cathode_screening()
        valid_c = cathode_df[cathode_df['valid'] == True].copy()
        top_c = valid_c.nsmallest(30, 'orr_overpotential_V') if 'orr_overpotential_V' in valid_c.columns else valid_c.head(30)

        pemfc_results = []
        for _, row in top_c.iterrows():
            mem = sweep_membranes(row['name'], row.get('orr_overpotential_V', 0.4))
            pemfc_results.extend(mem)
        if pemfc_results:
            best = max(pemfc_results, key=lambda r: r.get('peak_power_W_cm2', 0))
            stack = model_stack(StackConfig(n_cells=400,
                cell_voltage_V=best.get('peak_voltage_V', 0.65),
                current_density_A_cm2=best.get('peak_current_A_cm2', 1.5)))
            pipeline_state['phase5'] = {
                'best_power': best.get('peak_power_W_cm2', 0),
                'stack_kW': stack.get('net_power_kW', 0),
                'elapsed_s': time.time() - t5,
            }
        save_json(pipeline_state, "pipeline_state.json")

    # ─── Phase 6: Report ─────────────────────────────────────────────────────
    print_banner("PHASE 6: REPORT")
    from pipeline.report_generator import generate_full_report
    pipeline_state = load_json("pipeline_state.json") or pipeline_state
    generate_full_report(pipeline_state)
    save_json(pipeline_state, "pipeline_state.json")

    total = time.time() - t_start
    print("\n" + "=" * 80)
    print(f"  CAMPAIGN v2 COMPLETE: {total:.0f}s ({total/3600:.2f} hours)")
    print(f"  Catalysts evaluated: {pipeline_state.get('phase1', {}).get('total_evaluated', '?'):,}")
    print("=" * 80)


if __name__ == '__main__':
    main()
