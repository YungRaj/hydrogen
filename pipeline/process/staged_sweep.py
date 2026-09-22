"""General two-stage quantized sweep: gridded, then targeted.

Any lever set can use this. The coarse (gridded) JSON is write-once so a
later targeted pass cannot erase it. Archive layout:

    results/sweeps/<name>/manifest.json
    results/sweeps/<name>/coarse.json
    results/sweeps/<name>/targeted.json
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Callable, Mapping, Optional, Sequence

from pipeline.common.utils import SWEEPS_DIR, repo_relative

RecordFilter = Callable[[dict], bool]


@dataclass
class SweepSpec:
    """Describe one family of quantized levers and how to score cells."""
    name: str
    levers: Mapping[str, Sequence[float]]
    score_key: str = 'CH4_conversion'
    constraint_key: Optional[str] = None
    record_filter: Optional[RecordFilter] = None
    keep_fraction_of_max: float = 0.30
    n_targeted_per_lever: int = 6
    hard_bounds: Mapping[str, tuple] = field(default_factory=dict)
    extend_edge: bool = True


def cartesian_cells(levers: Mapping[str, Sequence]) -> list[dict]:
    keys = list(levers)
    if not keys:
        return [{}]
    cells = []
    for values in product(*(levers[k] for k in keys)):
        cells.append({k: v for k, v in zip(keys, values)})
    return cells


def sweep_dir(name: str) -> Path:
    return SWEEPS_DIR / name


def stage_path(name: str, stage: str) -> Path:
    if stage not in ('coarse', 'targeted'):
        raise ValueError(f'stage must be coarse or targeted, got {stage!r}')
    return sweep_dir(name) / f'{stage}.json'


def load_stage(name: str, stage: str) -> Optional[dict]:
    path = stage_path(name, stage)
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding='utf-8'))


def load_manifest(name: str) -> Optional[dict]:
    path = sweep_dir(name) / 'manifest.json'
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding='utf-8'))


def load_sweep(name: str) -> dict:
    """Both stages plus the manifest. Coarse is never dropped when targeted exists."""
    return {
        'name': name,
        'coarse': load_stage(name, 'coarse'),
        'targeted': load_stage(name, 'targeted'),
        'manifest': load_manifest(name),
    }


def write_stage(name: str, stage: str, payload: dict,
                overwrite: bool = False) -> Path:
    """Persist one stage. Coarse refuses overwrite unless explicitly forced."""
    if stage not in ('coarse', 'targeted'):
        raise ValueError(f'stage must be coarse or targeted, got {stage!r}')
    d = sweep_dir(name)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f'{stage}.json'
    if stage == 'coarse' and path.exists() and not overwrite:
        raise FileExistsError(
            f'{path} already exists; the gridded sweep is write-once. '
            'Pass overwrite=True only if you intend to replace it.')
    body = dict(payload)
    body.setdefault('stage', stage)
    body.setdefault('sweep_name', name)
    path.write_text(json.dumps(body, indent=2), encoding='utf-8')
    _update_manifest(name, stage, path, body)
    return path


def _update_manifest(name: str, stage: str, path: Path, payload: dict) -> None:
    man_path = sweep_dir(name) / 'manifest.json'
    manifest = {}
    if man_path.exists():
        manifest = json.loads(man_path.read_text(encoding='utf-8'))
    stages = dict(manifest.get('stages') or {})
    for entry in stages.values():
        if isinstance(entry, dict) and entry.get('path'):
            entry['path'] = repo_relative(entry['path'])
    stages[stage] = {
        'path': repo_relative(path),
        'written_at': datetime.now(timezone.utc).isoformat(),
        'n_records': len(payload.get('records') or []),
        'n_grid_cells': payload.get('n_grid_cells'),
        'levels': payload.get('levels'),
    }
    manifest.update({
        'name': name,
        'stages': stages,
        'coarse_preserved': (sweep_dir(name) / 'coarse.json').exists(),
    })
    if payload.get('roi'):
        manifest['roi'] = payload['roi']
    man_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')


def _feasible_records(records: Sequence[dict], spec: SweepSpec) -> list[dict]:
    out = []
    for row in records:
        if spec.record_filter is not None and not spec.record_filter(row):
            continue
        if spec.constraint_key is not None and row.get(spec.constraint_key) is False:
            continue
        score = row.get(spec.score_key)
        if score is None:
            continue
        out.append(row)
    return out


def _extend_bounds(coarse_levels: Sequence[float], lo: float, hi: float) -> tuple:
    levels = sorted(float(x) for x in coarse_levels)
    if len(levels) < 2:
        return lo, hi
    if lo <= levels[0] * 1.0000001:
        if levels[0] > 0 and levels[1] > 0:
            lo = levels[0] * (levels[0] / levels[1])
        else:
            lo = levels[0] - (levels[1] - levels[0])
    if hi >= levels[-1] * 0.9999999:
        if levels[-1] > 0 and levels[-2] > 0:
            hi = levels[-1] * (levels[-1] / levels[-2])
        else:
            hi = levels[-1] + (levels[-1] - levels[-2])
    return lo, hi


def _clip_hard(name: str, lo: float, hi: float,
               hard_bounds: Mapping[str, tuple]) -> tuple:
    if name not in hard_bounds:
        return lo, hi
    hlo, hhi = hard_bounds[name]
    return max(lo, float(hlo)), min(hi, float(hhi))


def densify_levels(coarse_levels: Sequence[float], lo: float, hi: float,
                   n: int) -> list[float]:
    kept = [float(x) for x in coarse_levels if lo - 1e-12 <= float(x) <= hi + 1e-12]
    if n < 2:
        n = 2
    step = (hi - lo) / (n - 1) if hi != lo else 0.0
    extra = [lo + i * step for i in range(n)]
    merged = sorted(set(round(x, 10) for x in kept + extra))
    return merged


def propose_roi(coarse_payload: dict, spec: SweepSpec) -> dict:
    """Bound each lever to the high-scoring feasible band, then densify.

    If the best cells sit on a coarse-grid edge, extend one step past that
    edge (so the targeted sweep can find a constraint wall the grid missed).
    """
    records = coarse_payload.get('records') or []
    feasible = _feasible_records(records, spec)
    if not feasible:
        raise ValueError('no feasible coarse records to propose an ROI from')
    scores = [float(r[spec.score_key]) for r in feasible]
    peak = max(scores)
    floor = peak * spec.keep_fraction_of_max
    kept = [r for r in feasible if float(r[spec.score_key]) >= floor]
    levels = {}
    bounds = {}
    notes = []
    for lever, coarse in spec.levers.items():
        values = [float(r[lever]) for r in kept if r.get(lever) is not None]
        if not values:
            raise ValueError(f'kept cells have no values for lever {lever!r}')
        lo, hi = min(values), max(values)
        if spec.extend_edge:
            lo, hi = _extend_bounds(coarse, lo, hi)
        lo, hi = _clip_hard(lever, lo, hi, spec.hard_bounds)
        if hi < lo:
            lo, hi = hi, lo
        bounds[lever] = [lo, hi]
        levels[lever] = densify_levels(coarse, lo, hi, spec.n_targeted_per_lever)
        notes.append(
            f'{lever} in [{lo:g}, {hi:g}] from {len(values)} cells '
            f'scoring >= {spec.keep_fraction_of_max:.0%} of max {peak:g}'
        )
    return {
        'levers': levels,
        'bounds': bounds,
        'n_grid_cells': len(cartesian_cells(levels)),
        'keep_fraction_of_max': spec.keep_fraction_of_max,
        'score_peak': peak,
        'score_floor': floor,
        'n_kept_coarse_cells': len(kept),
        'reason': '; '.join(notes),
    }


def import_existing_stage(name: str, stage: str, source: Path,
                          overwrite: bool = False) -> Path:
    """Copy an already-run JSON into the write-once archive without resimulating."""
    payload = json.loads(Path(source).read_text(encoding='utf-8'))
    return write_stage(name, stage, payload, overwrite=overwrite)
