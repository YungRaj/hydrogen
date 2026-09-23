"""B6-6 acceptance criteria over the five Ni star / 2-D sweep run.json files.

Reads results/sweeps/ni_np_b66_*/run.json and writes
results/sweeps/b66_summary.json. Prints a short table.

The five gates are the user-approved set. 1300 K rows with X > X_eq are
flagged and never counted as successes.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

from pipeline.utils import SWEEPS_DIR, repo_relative
from pipeline.reactors.eligibility import is_usable_result

SWEEP_NAMES = (
    'ni_np_b66_eact',
    'ni_np_b66_agamma',
    'ni_np_b66_theta',
    'ni_np_b66_sticking',
    'ni_np_b66_agamma_theta',
)

BASE_KINETICS = {
    'E_act': 1.0,
    'carbon_transfer_prefactor_1_s': 1.0e13,
    'encapsulation_crossover_coverage': 0.5,
    'ch4_sticking_coefficient': 0.01,
}

FLAT_THRESHOLD_PER_EV = 5.0
LINEARITY_THRESHOLD = 0.20
ENCAP_ONSET_THETA = 0.1
ENCAP_ONSET_RATIO = 10.0
NI_TOS_LIFETIME_BAND_H = (4.0, 50.0)
NI_FILAMENT_YIELD_BAND = (8.0, 10.0)

# Neighbours used for the central-difference flat test at E_act = 1.0.
FLAT_E_LO = 0.9
FLAT_E_HI = 1.1
FLAT_E_BASE = 1.0


def _finite(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _close(a: Any, b: Any, rel: float = 1e-6, abs_tol: float = 1e-12) -> bool:
    left, right = _finite(a), _finite(b)
    if left is None or right is None:
        return False
    return abs(left - right) <= max(abs_tol, rel * max(abs(left), abs(right)))


def sweep_value(record: dict, key: str) -> Any:
    """Grid-point value, falling back to the record / job-fixed field.

    Args:
        key: Input controlling key.

    Returns:
        Validated Any output for this operation.
    """
    swept = record.get('sweep') or {}
    if key in swept and swept[key] is not None:
        return swept[key]
    return record.get(key)


def is_complete(record: dict) -> bool:
    """Return whether a sweep record completed successfully.

    Args:
        record: Input controlling record.

    Returns:
        Validated bool output for this operation.
    """
    return record.get('status') == 'complete'


def exceeds_equilibrium(record: dict) -> bool:
    """Return whether a record violates its equilibrium conversion limit.

    Args:
        record: Input controlling record.

    Returns:
        Validated bool output for this operation.
    """
    return bool(record.get('exceeds_equilibrium'))


def is_scorable(record: dict) -> bool:
    """Shared usable baseline: complete, not mock, no X>X_eq, carbon OK.

    Args:
        record: Input controlling record.

    Returns:
        Validated bool output for this operation.
    """
    return is_usable_result(record)


def match_kinetics(record: dict, kinetics: Optional[dict] = None) -> bool:
    """Return whether a record matches every requested kinetic parameter.

    Args:
        record: Input controlling record.
        kinetics: Input controlling kinetics.

    Returns:
        Validated bool output for this operation.
    """
    target = BASE_KINETICS if kinetics is None else kinetics
    for key, expected in target.items():
        got = sweep_value(record, key)
        if got is None:
            continue
        if not _close(got, expected):
            return False
    return True


def select_records(
    records: Iterable[dict],
    *,
    cell: Optional[str] = None,
    reactor: Optional[str] = None,
    T_K: Optional[float] = None,
    regen: Optional[int] = None,
    kinetics: Optional[dict] = None,
    complete_only: bool = True,
    scorable_only: bool = False,
) -> list:
    """Filter sweep records by operating cell, reactor, kinetics, and validity.

    Args:
        records: Input controlling records.
        cell: Input controlling cell.
        reactor: Input controlling reactor.
        T_K: Input controlling T K.
        regen: Input controlling regen.
        kinetics: Input controlling kinetics.
        complete_only: Input controlling complete only.
        scorable_only: Input controlling scorable only.

    Returns:
        Validated list output for this operation.
    """
    out = []
    for rec in records:
        if complete_only and not is_complete(rec):
            continue
        if scorable_only and not is_scorable(rec):
            continue
        if cell is not None and rec.get('cell') != cell:
            continue
        if reactor is not None and rec.get('reactor_type') != reactor:
            continue
        if T_K is not None and not _close(rec.get('T_K'), T_K, rel=0.0, abs_tol=0.2):
            continue
        if regen is not None:
            got = sweep_value(rec, 'max_regen_cycles')
            if got is None or int(got) != int(regen):
                continue
        if kinetics is not None and not match_kinetics(rec, kinetics):
            continue
        out.append(rec)
    return out


def d_lnX_d_Eact(x_lo: float, x_hi: float, e_lo: float = FLAT_E_LO,
                 e_hi: float = FLAT_E_HI) -> float:
    """Central difference |d ln X / d E_act| from two neighbour conversions.

    Args:
        x_lo: Input controlling x lo.
        x_hi: Input controlling x hi.
        e_lo: Input controlling e lo.
        e_hi: Input controlling e hi.

    Returns:
        Validated float output for this operation.
    """
    if x_lo <= 0 or x_hi <= 0:
        raise ValueError('conversions must be positive for ln X')
    denom = float(e_hi) - float(e_lo)
    if denom == 0:
        raise ValueError('E_act neighbours must differ')
    return abs((math.log(x_hi) - math.log(x_lo)) / denom)


def linearity_relative_diff(turnovers_a: float, turnovers_b: float) -> float:
    """|t_a - t_b| / mean. Pure a-scaling gives 0 (equal turnovers).

    Args:
        turnovers_a: Input controlling turnovers a.
        turnovers_b: Input controlling turnovers b.

    Returns:
        Validated float output for this operation.
    """
    mean = 0.5 * (float(turnovers_a) + float(turnovers_b))
    if mean == 0:
        raise ValueError('mean turnovers is zero')
    return abs(float(turnovers_a) - float(turnovers_b)) / mean


def in_band(value: Optional[float], lo: float, hi: float) -> bool:
    """Return whether a finite value lies inside an inclusive interval.

    Args:
        value: Input controlling value.
        lo: Input controlling lo.
        hi: Input controlling hi.

    Returns:
        Validated bool output for this operation.
    """
    number = _finite(value)
    return False if number is None else (lo <= number <= hi)


def first_where(
    records: Iterable[dict],
    predicate: Callable[[dict], bool],
    sort_key: Callable[[dict], float],
) -> Optional[dict]:
    """Return the first sorted record satisfying a predicate.

    Args:
        records: Input controlling records.
        predicate: Input controlling predicate.
        sort_key: Input controlling sort key.

    Returns:
        Validated Optional[dict] output for this operation.
    """
    ranked = sorted(
        (rec for rec in records if predicate(rec)),
        key=sort_key)
    return ranked[0] if ranked else None


def _all_records(payloads: dict) -> list:
    records = []
    for name, payload in payloads.items():
        for rec in payload.get('records') or []:
            row = dict(rec)
            row['_sweep'] = name
            records.append(row)
    return records


def _row_brief(record: Optional[dict]) -> Optional[dict]:
    if record is None:
        return None
    return {
        'sweep': record.get('_sweep'),
        'cell': record.get('cell'),
        'reactor_type': record.get('reactor_type'),
        'T_K': record.get('T_K'),
        'sweep_values': dict(record.get('sweep') or {}),
        'CH4_conversion': record.get('CH4_conversion'),
        'exit_theta_C_encap': record.get('exit_theta_C_encap'),
        'c_gamma_to_c_delta_ratio': record.get('c_gamma_to_c_delta_ratio'),
        'encapsulation_lifetime_h': record.get('encapsulation_lifetime_h'),
        'filament_yield_gC_per_gMetal_h': record.get(
            'filament_yield_gC_per_gMetal_h'),
        'exceeds_equilibrium': exceeds_equilibrium(record),
    }


def evaluate_criteria(payloads: dict) -> dict:
    """Score the five B6-6 gates on already-loaded run payloads.

    Returns:
        Validated dict output for this operation.
    """
    records = _all_records(payloads)
    flagged_overshoot = [
        _row_brief(rec) for rec in records
        if is_complete(rec) and exceeds_equilibrium(rec)]

    # Flat: |d ln X / d E_act| at production / PFR / 923 K / regen 0 / E=1.0.
    eact_rows = select_records(
        records, cell='production', reactor='PFR', T_K=923.15, regen=0)
    x_by_e = {}
    for rec in eact_rows:
        e_act = _finite(sweep_value(rec, 'E_act'))
        x = _finite(rec.get('CH4_conversion'))
        if e_act is None or x is None:
            continue
        if match_kinetics(rec, {k: v for k, v in BASE_KINETICS.items()
                                if k != 'E_act'}):
            x_by_e[round(e_act, 6)] = rec
    lo = x_by_e.get(FLAT_E_LO)
    hi = x_by_e.get(FLAT_E_HI)
    base = x_by_e.get(FLAT_E_BASE)
    slope = None
    if lo is not None and hi is not None:
        slope = d_lnX_d_Eact(
            float(lo['CH4_conversion']), float(hi['CH4_conversion']))
    flat = {
        'name': 'flat',
        'description': (
            '|d ln X / d E_act| at production PFR 923 K regen 0, '
            'central difference from E_act 0.9 and 1.1'),
        'value': slope,
        'threshold': FLAT_THRESHOLD_PER_EV,
        'unit': '1/eV',
        'pass': bool(slope is not None and slope >= FLAT_THRESHOLD_PER_EV),
        'X_0.9': None if lo is None else lo.get('CH4_conversion'),
        'X_1.0': None if base is None else base.get('CH4_conversion'),
        'X_1.1': None if hi is None else hi.get('CH4_conversion'),
    }

    # Linearity: turnovers on the two cells at identical base kinetics.
    prod = select_records(
        records, cell='production', reactor='PFR', T_K=923.15, regen=0,
        kinetics=BASE_KINETICS)
    large = select_records(
        records, cell='large_particle_ni', reactor='PFR', T_K=923.15, regen=0,
        kinetics=BASE_KINETICS)
    t_prod = _finite(prod[0].get('carbon_turnovers_per_site')) if prod else None
    t_large = _finite(large[0].get('carbon_turnovers_per_site')) if large else None
    rel = None
    if t_prod is not None and t_large is not None:
        rel = linearity_relative_diff(t_prod, t_large)
    linearity = {
        'name': 'linearity',
        'description': (
            '|turnovers_production - turnovers_large| / mean at identical '
            'base kinetics, production vs large_particle_ni, PFR 923 K regen 0'),
        'value': rel,
        'threshold': LINEARITY_THRESHOLD,
        'turnovers_production': t_prod,
        'turnovers_large_particle_ni': t_large,
        'pass': bool(rel is not None and rel > LINEARITY_THRESHOLD),
    }

    # Encapsulation onset: first scorable T (base kinetics) and first A_gamma
    # (923 K, other knobs at base) where exit theta_encap >= 0.1; same for
    # C_gamma/C_delta < 10.
    def _onset_pred_theta(rec):
        th = _finite(rec.get('exit_theta_C_encap'))
        return th is not None and th >= ENCAP_ONSET_THETA

    def _onset_pred_ratio(rec):
        ratio = _finite(rec.get('c_gamma_to_c_delta_ratio'))
        return ratio is not None and ratio < ENCAP_ONSET_RATIO

    t_axis = select_records(
        records, cell='production', reactor='PFR', regen=0,
        kinetics=BASE_KINETICS, scorable_only=True)
    a_axis = select_records(
        records, cell='production', reactor='PFR', T_K=923.15, regen=0,
        kinetics={k: v for k, v in BASE_KINETICS.items()
                  if k != 'carbon_transfer_prefactor_1_s'},
        scorable_only=True)
    first_T = first_where(t_axis, _onset_pred_theta, lambda r: float(r['T_K']))
    first_A = first_where(
        a_axis, _onset_pred_theta,
        lambda r: float(sweep_value(r, 'carbon_transfer_prefactor_1_s')))
    first_T_ratio = first_where(
        t_axis, _onset_pred_ratio, lambda r: float(r['T_K']))
    first_A_ratio = first_where(
        a_axis, _onset_pred_ratio,
        lambda r: float(sweep_value(r, 'carbon_transfer_prefactor_1_s')))
    onset = {
        'name': 'encapsulation_onset',
        'description': (
            f'first T or A_gamma where exit theta_encap >= {ENCAP_ONSET_THETA}; '
            f'secondary C_gamma/C_delta < {ENCAP_ONSET_RATIO} '
            '(scorable rows only)'),
        'first_T_theta_encap': _row_brief(first_T),
        'first_A_gamma_theta_encap': _row_brief(first_A),
        'first_T_ratio': _row_brief(first_T_ratio),
        'first_A_gamma_ratio': _row_brief(first_A_ratio),
        'pass': bool(first_T is not None or first_A is not None),
    }

    # Lifetime and yield: any scorable row in the literature band.
    scorable = [rec for rec in records if is_scorable(rec)]
    life_hits = [
        rec for rec in scorable
        if in_band(rec.get('encapsulation_lifetime_h'), *NI_TOS_LIFETIME_BAND_H)]
    yield_hits = [
        rec for rec in scorable
        if in_band(rec.get('filament_yield_gC_per_gMetal_h'),
                   *NI_FILAMENT_YIELD_BAND)]
    lifetime = {
        'name': 'lifetime',
        'description': (
            f'encapsulation_lifetime_h vs NI TOS band '
            f'{NI_TOS_LIFETIME_BAND_H[0]:g}–{NI_TOS_LIFETIME_BAND_H[1]:g} h'),
        'band_h': list(NI_TOS_LIFETIME_BAND_H),
        'n_hits': len(life_hits),
        'example': _row_brief(life_hits[0]) if life_hits else None,
        'pass': bool(life_hits),
    }
    yield_crit = {
        'name': 'yield',
        'description': (
            f'filament yield vs Ermakova '
            f'{NI_FILAMENT_YIELD_BAND[0]:g}–{NI_FILAMENT_YIELD_BAND[1]:g} '
            'gC/(gNi·h)'),
        'band_gC_per_gNi_h': list(NI_FILAMENT_YIELD_BAND),
        'n_hits': len(yield_hits),
        'example': _row_brief(yield_hits[0]) if yield_hits else None,
        'pass': bool(yield_hits),
    }

    criteria = [flat, linearity, onset, lifetime, yield_crit]
    return {
        'schema_version': 1,
        'sweeps': {
            name: {
                'n_records': len((payloads.get(name) or {}).get('records') or []),
                'n_complete': sum(
                    1 for r in (payloads.get(name) or {}).get('records') or []
                    if is_complete(r)),
                'n_overshoot': sum(
                    1 for r in (payloads.get(name) or {}).get('records') or []
                    if is_complete(r) and exceeds_equilibrium(r)),
            }
            for name in SWEEP_NAMES
        },
        'n_flagged_X_gt_Xeq': len(flagged_overshoot),
        'flagged_X_gt_Xeq_note': (
            'X > X_eq rows (typically 1300 K with fast C_gamma) are flagged '
            'and are not scored as successes'),
        'criteria': {c['name']: c for c in criteria},
        'n_pass': sum(1 for c in criteria if c['pass']),
        'n_fail': sum(1 for c in criteria if not c['pass']),
    }


def load_runs(sweeps_dir: Optional[Path] = None) -> dict:
    """Load the B6-6 sweep payloads required by the criteria evaluator.

    Args:
        sweeps_dir: Input controlling sweeps dir.

    Returns:
        Validated dict output for this operation.
    """
    root = Path(sweeps_dir) if sweeps_dir is not None else SWEEPS_DIR
    payloads = {}
    missing = []
    for name in SWEEP_NAMES:
        path = root / name / 'run.json'
        if not path.is_file():
            missing.append(repo_relative(path) if path.is_absolute() else str(path))
            continue
        payloads[name] = json.loads(path.read_text(encoding='utf-8'))
    if missing:
        raise FileNotFoundError(
            'missing B6-6 run.json file(s): ' + ', '.join(missing))
    return payloads


def print_table(summary: dict) -> None:
    """Print a compact human-readable B6-6 criteria summary.

    Args:
        summary: Input controlling summary.
    """
    print('B6-6 criteria')
    print(f"flagged X>X_eq rows: {summary['n_flagged_X_gt_Xeq']} "
          '(not scored as successes)')
    for name, stats in summary['sweeps'].items():
        print(
            f"  {name:<24} {stats['n_complete']}/{stats['n_records']} complete "
            f"({stats['n_overshoot']} overshoot)")
    print()
    print(f"{'criterion':<22} {'value':>12} {'gate':>10} {'result':>8}")
    for name, crit in summary['criteria'].items():
        value = crit.get('value')
        if name == 'encapsulation_onset':
            first_A = crit.get('first_A_gamma_theta_encap') or {}
            first_T = crit.get('first_T_theta_encap') or {}
            a_val = (first_A.get('sweep_values') or {}).get(
                'carbon_transfer_prefactor_1_s')
            t_val = first_T.get('T_K')
            if a_val is not None:
                shown = f'Aγ={a_val:.2g}'
            elif t_val is not None:
                shown = f'T={t_val:g}K'
            else:
                shown = 'none'
            gate = f'θ≥{ENCAP_ONSET_THETA:g}'
        elif name == 'lifetime':
            shown = f"{crit['n_hits']} hits"
            gate = f'{NI_TOS_LIFETIME_BAND_H[0]:g}–{NI_TOS_LIFETIME_BAND_H[1]:g}h'
        elif name == 'yield':
            shown = f"{crit['n_hits']} hits"
            gate = (f'{NI_FILAMENT_YIELD_BAND[0]:g}–'
                    f'{NI_FILAMENT_YIELD_BAND[1]:g}')
        elif name == 'flat':
            shown = '—' if value is None else f'{value:.3g}'
            gate = f'≥{FLAT_THRESHOLD_PER_EV:g}/eV'
        else:
            shown = '—' if value is None else f'{value:.3g}'
            gate = f'>{LINEARITY_THRESHOLD:.0%}'
        mark = 'PASS' if crit['pass'] else 'FAIL'
        print(f'{name:<22} {shown:>12} {gate:>10} {mark:>8}')
    print(f"\n{summary['n_pass']} pass, {summary['n_fail']} fail")


def write_summary(summary: dict, path: Optional[Path] = None) -> Path:
    """Persist a B6-6 criteria summary as JSON.

    Args:
        summary: Input controlling summary.
        path: Input controlling path.

    Returns:
        Validated Path output for this operation.
    """
    out = Path(path) if path is not None else SWEEPS_DIR / 'b66_summary.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return out


def main() -> dict:
    """Run the module command-line workflow.

    Returns:
        Validated dict output for this operation.
    """
    payloads = load_runs()
    summary = evaluate_criteria(payloads)
    out = write_summary(summary)
    print_table(summary)
    print(f'\nWrote {repo_relative(out)}')
    return summary


if __name__ == '__main__':
    main()
