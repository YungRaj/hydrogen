"""B6-7 joint-band search: same row in both Ermakova yield and TOS lifetime.

Exhausts the B6-6 star / 2-D records plus an optional A_γ × θ* × s0 cube
at 923.15 / 973.15 K. A hit is one complete, non-overshooting row whose
filament yield is in 8–10 gC/(gNi·h) AND whose encapsulation lifetime is
in 4–50 h. No new physics is added to force an overlap.

Writes results/sweeps/b67_joint_band.json.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Optional

from pipeline.common.utils import SWEEPS_DIR, repo_relative
from pipeline.process.b66_criteria import (
    NI_FILAMENT_YIELD_BAND,
    NI_TOS_LIFETIME_BAND_H,
    SWEEP_NAMES,
    _finite,
    _row_brief,
    in_band,
    is_scorable,
    load_runs,
    sweep_value,
)

JOINT_SWEEP = 'ni_np_b67_joint'
YIELD_WIDTH = NI_FILAMENT_YIELD_BAND[1] - NI_FILAMENT_YIELD_BAND[0]
LIFE_WIDTH = NI_TOS_LIFETIME_BAND_H[1] - NI_TOS_LIFETIME_BAND_H[0]
ROI_T_MIN = 923.15
ROI_T_MAX = 973.15


def _distance_to_band(value: Optional[float], lo: float, hi: float) -> Optional[float]:
    number = _finite(value)
    if number is None:
        return None
    if lo <= number <= hi:
        return 0.0
    return lo - number if number < lo else number - hi


def _normalized_miss(record: dict) -> Optional[float]:
    """0 if both bands hit; else hypot of yield/life distances over band widths."""
    y = _distance_to_band(
        record.get('filament_yield_gC_per_gMetal_h'), *NI_FILAMENT_YIELD_BAND)
    t = _distance_to_band(
        record.get('encapsulation_lifetime_h'), *NI_TOS_LIFETIME_BAND_H)
    if y is None or t is None:
        return None
    return math.hypot(y / YIELD_WIDTH, t / LIFE_WIDTH)


def search_joint_band(payloads: dict) -> dict:
    records = []
    for name, payload in payloads.items():
        for rec in payload.get('records') or []:
            row = dict(rec)
            row['_sweep'] = name
            records.append(row)
    scorable = [rec for rec in records if is_scorable(rec)]
    yield_hits = [
        rec for rec in scorable
        if in_band(rec.get('filament_yield_gC_per_gMetal_h'),
                   *NI_FILAMENT_YIELD_BAND)]
    life_hits = [
        rec for rec in scorable
        if in_band(rec.get('encapsulation_lifetime_h'), *NI_TOS_LIFETIME_BAND_H)]
    both = [
        rec for rec in scorable
        if in_band(rec.get('filament_yield_gC_per_gMetal_h'),
                   *NI_FILAMENT_YIELD_BAND)
        and in_band(rec.get('encapsulation_lifetime_h'), *NI_TOS_LIFETIME_BAND_H)]
    both_roi = [
        rec for rec in both
        if _finite(rec.get('T_K')) is not None
        and ROI_T_MIN - 0.2 <= float(rec['T_K']) <= ROI_T_MAX + 0.2]

    ranked = []
    for rec in scorable:
        miss = _normalized_miss(rec)
        if miss is None:
            continue
        ranked.append((miss, rec))
    ranked.sort(key=lambda item: item[0])
    closest = [
        {**_row_brief(rec), 'normalized_miss': miss}
        for miss, rec in ranked[:8]
    ]

    declaration = (
        'no_simultaneous_hit_in_filament_ROI' if not both_roi else 'simultaneous_hit'
    )
    return {
        'schema_version': 1,
        'definition': (
            'A joint hit is one scorable row (complete, not X>X_eq) with '
            f'filament_yield in {list(NI_FILAMENT_YIELD_BAND)} gC/(gNi·h) '
            f'AND encapsulation_lifetime_h in {list(NI_TOS_LIFETIME_BAND_H)}. '
            'Same cell, reactor, T, and (A_γ, θ*, s0).'
        ),
        'declaration': declaration,
        'n_records': len(records),
        'n_scorable': len(scorable),
        'n_yield_only': len(yield_hits) - len(both),
        'n_lifetime_only': len(life_hits) - len(both),
        'n_both': len(both),
        'n_both_roi': len(both_roi),
        'both': [_row_brief(rec) for rec in both],
        'both_roi': [_row_brief(rec) for rec in both_roi],
        'closest_misses': closest,
        'read': (
            'In-band yield is a high-Cδ / low-A_γ row (surface already '
            'encapsulating). In-band lifetime is a high-A_γ / low-T / low-s0 '
            'row where Cδ is a trickle. The arrival rate that fills the '
            'Ermakova band also fills θ* inside the pass; the transport '
            'that keeps θ_encap a trickle over-produces filaments. No extra '
            'DOF was added to force an overlap.'
        ),
    }


def load_joint_payloads(sweeps_dir: Optional[Path] = None) -> dict:
    root = Path(sweeps_dir) if sweeps_dir is not None else SWEEPS_DIR
    payloads = load_runs(root)
    joint = root / JOINT_SWEEP / 'run.json'
    if joint.is_file():
        payloads[JOINT_SWEEP] = json.loads(joint.read_text(encoding='utf-8'))
    return payloads


def write_summary(summary: dict, path: Optional[Path] = None) -> Path:
    out = Path(path) if path is not None else SWEEPS_DIR / 'b67_joint_band.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return out


def print_table(summary: dict) -> None:
    print('B6-7 joint band')
    print(f"declaration: {summary['declaration']}")
    print(
        f"scorable {summary['n_scorable']}: "
        f"{summary['n_yield_only']} yield-only, "
        f"{summary['n_lifetime_only']} lifetime-only, "
        f"{summary['n_both']} both "
        f"({summary['n_both_roi']} in 923–973 K ROI)")
    if summary['both']:
        print('hits (all T):')
        for row in summary['both']:
            print(
                f"  {row.get('sweep')} {row.get('cell')} "
                f"{row.get('reactor_type')} {row.get('T_K')} "
                f"Y={row.get('filament_yield_gC_per_gMetal_h')} "
                f"τ={row.get('encapsulation_lifetime_h')}")
    print('closest misses (normalized):')
    for row in summary.get('closest_misses') or []:
        print(
            f"  miss={row['normalized_miss']:.3g}  "
            f"{row.get('sweep')} {row.get('cell')} "
            f"{row.get('reactor_type')} T={row.get('T_K')} "
            f"Y={row.get('filament_yield_gC_per_gMetal_h')} "
            f"τ={row.get('encapsulation_lifetime_h')} "
            f"{row.get('sweep_values')}")


def main() -> dict:
    payloads = load_joint_payloads()
    summary = search_joint_band(payloads)
    out = write_summary(summary)
    print_table(summary)
    print(f'\nWrote {repo_relative(out)}')
    return summary


if __name__ == '__main__':
    main()
