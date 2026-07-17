#!/usr/bin/env python3
"""
Live Dashboard — Real-time monitoring of the catalyst discovery campaign.

Refreshes every 5 seconds showing:
  - GPU utilization per device
  - GA generation progress + best metrics
  - MACE screening throughput
  - Material class distribution
  - Pareto front summary

Usage:
    conda run -n battery-env python live_dashboard.py
    # or just:
    python live_dashboard.py
"""

import os
import sys
import time
import glob
import subprocess
from pathlib import Path
from datetime import datetime, timedelta

BASE = Path(__file__).parent
SCREENING = BASE / "results" / "screening"
GA_LOG = SCREENING / "genetic_optimizer.log"
MACE_LOG = SCREENING / "mace_screening.log"


def clear():
    os.system('clear')


def gpu_status():
    """Get GPU utilization via nvidia-smi."""
    try:
        out = subprocess.check_output([
            'nvidia-smi',
            '--query-gpu=name,utilization.gpu,memory.used,memory.total,temperature.gpu,power.draw',
            '--format=csv,noheader,nounits'
        ], text=True, timeout=5)
        gpus = []
        for line in out.strip().split('\n'):
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 6:
                gpus.append({
                    'name': parts[0][:30],
                    'util': int(parts[1]),
                    'mem_used': int(parts[2]),
                    'mem_total': int(parts[3]),
                    'temp': int(parts[4]),
                    'power': float(parts[5]),
                })
        return gpus
    except Exception:
        return []


def parse_ga_log():
    """Parse the GA log for generation progress."""
    if not GA_LOG.exists():
        return []
    entries = []
    try:
        with open(GA_LOG) as f:
            for line in f:
                if 'Gen ' in line and 'best_E_act' in line:
                    # [2026-07-16 15:07:06] genetic_optimizer | INFO |   Gen 1/3000: best_E_act=...
                    parts = line.split('Gen ')[1]
                    gen_str = parts.split('/')[0].strip()
                    total_str = parts.split('/')[1].split(':')[0].strip()
                    e_act = float(parts.split('best_E_act=')[1].split(' ')[0])
                    pareto = int(parts.split('pareto_size=')[1].split(',')[0])
                    diversity = int(parts.split('pop_diversity=')[1].split(',')[0])
                    elapsed = parts.split('(')[1].split('s)')[0]
                    timestamp = line.split(']')[0].strip('[')
                    entries.append({
                        'gen': int(gen_str),
                        'total': int(total_str),
                        'e_act': e_act,
                        'pareto': pareto,
                        'diversity': diversity,
                        'elapsed_s': float(elapsed),
                        'timestamp': timestamp,
                    })
                elif 'Running MACE' in line:
                    timestamp = line.split(']')[0].strip('[')
                    entries.append({'type': 'mace', 'timestamp': timestamp})
                elif 'Retraining surrogate' in line:
                    samples = int(line.split('on ')[1].split(' ')[0])
                    entries.append({'type': 'retrain', 'samples': samples})
    except Exception:
        pass
    return entries


def parse_mace_log():
    """Parse the MACE log for throughput info."""
    if not MACE_LOG.exists():
        return {}
    info = {}
    try:
        with open(MACE_LOG) as f:
            lines = f.readlines()
        for line in reversed(lines):
            if 'Progress:' in line and 'candidates/sec' in line:
                rate = float(line.split('(')[1].split(' candidates')[0])
                done = line.split('Progress: ')[1].split('/')[0]
                total = line.split('/')[1].split(' ')[0]
                valid = line.split(', ')[1].split(' valid')[0]
                info['rate'] = rate
                info['done'] = int(done)
                info['total'] = int(total)
                info['valid'] = int(valid)
                break
            elif 'Screening complete' in line:
                info['status'] = 'complete'
                break
            elif 'Screening ' in line and 'catalyst' in line:
                n = int(line.split('Screening ')[1].split(' catalyst')[0])
                info['screening'] = n
    except Exception:
        pass
    return info


def count_csvs():
    """Count total MACE evaluations from CSV files."""
    total = 0
    n_files = 0
    class_counts = {}
    try:
        import pandas as pd
        for f in glob.glob(str(SCREENING / "ga_*.csv")):
            df = pd.read_csv(f)
            total += len(df)
            n_files += 1
            if 'material_class' in df.columns:
                for cls, count in df['material_class'].value_counts().items():
                    class_counts[cls] = class_counts.get(cls, 0) + count
    except ImportError:
        # No pandas — just count files
        for f in glob.glob(str(SCREENING / "ga_*.csv")):
            n_files += 1
            with open(f) as fh:
                total += sum(1 for _ in fh) - 1
    return total, n_files, class_counts


def bar(value, max_val, width=30, filled='█', empty='░'):
    """Render a progress bar."""
    if max_val <= 0:
        return empty * width
    ratio = min(value / max_val, 1.0)
    n = int(ratio * width)
    return filled * n + empty * (width - n)


def colorize(text, code):
    """ANSI color wrapper."""
    return f"\033[{code}m{text}\033[0m"


def gpu_bar(util):
    """Colored GPU utilization bar."""
    b = bar(util, 100, width=20)
    if util > 60:
        return colorize(b, '92')  # green
    elif util > 20:
        return colorize(b, '93')  # yellow
    else:
        return colorize(b, '91')  # red


def main():
    print("Starting live dashboard... (Ctrl+C to exit)\n")
    time.sleep(1)

    while True:
        try:
            clear()
            now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

            # ─── Header ──────────────────────────────────────────────────
            print(colorize("═" * 80, '96'))
            print(colorize("  ⚡ TURQUOISE H₂ CATALYST DISCOVERY — LIVE DASHBOARD", '1;96'))
            print(colorize(f"  {now}", '90'))
            print(colorize("═" * 80, '96'))

            # ─── GPU Status ──────────────────────────────────────────────
            gpus = gpu_status()
            print(colorize("\n  ┌─ GPU STATUS ─────────────────────────────────────────────┐", '94'))
            for i, g in enumerate(gpus):
                util_bar = gpu_bar(g['util'])
                mem_pct = 100 * g['mem_used'] / g['mem_total'] if g['mem_total'] > 0 else 0
                status = "🔥 ACTIVE" if g['util'] > 10 else "💤 idle"
                print(f"  │ GPU {i}: {g['name']:25s} {util_bar} {g['util']:>3}%  "
                      f"{g['mem_used']:>5}/{g['mem_total']} MB  {g['temp']}°C  {status} │")
            if not gpus:
                print("  │  No GPUs detected                                       │")
            print(colorize("  └─────────────────────────────────────────────────────────┘", '94'))

            # ─── GA Progress ─────────────────────────────────────────────
            ga_entries = parse_ga_log()
            gen_entries = [e for e in ga_entries if 'gen' in e]
            print(colorize("\n  ┌─ GENETIC ALGORITHM ───────────────────────────────────────┐", '93'))

            if gen_entries:
                latest = gen_entries[-1]
                gen = latest['gen']
                total = latest['total']
                pct = 100 * gen / total if total > 0 else 0
                progress = bar(gen, total, width=30)

                print(f"  │ Generation:  {colorize(f'{gen:>5}', '1')}/{total}  {progress} {pct:5.1f}%")
                e_act_color = '92' if latest['e_act'] < 0 else '91'
                e_act_str = colorize(f"{latest['e_act']:>8.4f} eV", e_act_color)
                print(f"  │ Best E_act:  {e_act_str}")
                print(f"  │ Pareto front: {latest['pareto']:>5} candidates")
                print(f"  │ Diversity:    {latest['diversity']:>5} unique genomes")
                print(f"  │ Gen time:     {latest['elapsed_s']:>7.1f}s")

                # ETA
                if len(gen_entries) >= 2:
                    t0 = gen_entries[0]
                    avg_s = (latest['elapsed_s'] + gen_entries[-2].get('elapsed_s', latest['elapsed_s'])) / 2
                    remaining_gens = total - gen
                    eta_s = remaining_gens * avg_s
                    eta_str = str(timedelta(seconds=int(eta_s)))
                    print(f"  │ ETA:          {colorize(eta_str, '96')}")

                # E_act history (last 5)
                if len(gen_entries) > 1:
                    hist = gen_entries[-min(5, len(gen_entries)):]
                    trend = " → ".join([f"{e['e_act']:.3f}" for e in hist])
                    print(f"  │ E_act trend:  {trend}")
            else:
                print("  │ Waiting for first generation...")

            print(colorize("  └─────────────────────────────────────────────────────────┘", '93'))

            # ─── MACE Screening ──────────────────────────────────────────
            mace = parse_mace_log()
            total_eval, n_files, class_counts = count_csvs()

            print(colorize("\n  ┌─ MACE SCREENING ──────────────────────────────────────────┐", '92'))
            print(f"  │ Total evaluated: {colorize(f'{total_eval:>6,}', '1')} across {n_files} rounds")
            if 'rate' in mace:
                print(f"  │ Current round:   {mace.get('done', '?')}/{mace.get('total', '?')}  "
                      f"({mace['rate']:.1f} cand/sec, {mace.get('valid', '?')} valid)")
            elif mace.get('status') == 'complete':
                print(f"  │ Status:          {colorize('Round complete', '92')}")
            print(colorize("  └─────────────────────────────────────────────────────────┘", '92'))

            # ─── Class Distribution ──────────────────────────────────────
            if class_counts:
                print(colorize("\n  ┌─ MATERIAL CLASS DISTRIBUTION ─────────────────────────────┐", '95'))
                max_count = max(class_counts.values()) if class_counts else 1
                all_classes = ['MoltenMetal', 'SolidCatalyst', 'SAC', 'DAC', 'MOF',
                               'COF', 'Perovskite', 'MetalHydride', 'MAXPhase', 'HEA']
                for cls in all_classes:
                    count = class_counts.get(cls, 0)
                    pct = 100 * count / total_eval if total_eval > 0 else 0
                    cls_bar = bar(count, max_count, width=20)
                    print(f"  │ {cls:18s} {cls_bar} {count:>5} ({pct:4.1f}%)")
                print(colorize("  └─────────────────────────────────────────────────────────┘", '95'))

            # ─── Process Info ────────────────────────────────────────────
            try:
                ps = subprocess.check_output(
                    "ps aux | grep run_production | grep -v grep | awk '{print $3, $4, $10}'",
                    shell=True, text=True, timeout=3
                ).strip()
                if ps:
                    parts = ps.split('\n')[0].split()
                    if len(parts) >= 3:
                        print(f"\n  Process: CPU={parts[0]}%  MEM={parts[1]}%  Time={parts[2]}")
            except Exception:
                pass

            print(colorize(f"\n  Refreshing every 5s... Press Ctrl+C to exit", '90'))
            time.sleep(5)

        except KeyboardInterrupt:
            print("\n\n  Dashboard stopped.")
            break


if __name__ == '__main__':
    main()
