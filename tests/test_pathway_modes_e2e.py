#!/usr/bin/env python3
"""Every pathway mode runs end to end through the campaign reactor stage.

Uses only tracked literature rows (no screening database, no GPU): the Ni
nanoparticle reference for solids modes, a Bi-Ni melt for MMBCR, and the same
Ni kinetics for NTEC / electrochemical. Requires Cantera and, for the NTEC and
electrochemical modes, the MOOSE executable (opt-in ``e2e`` suite).
"""

from pathlib import Path
import sys
import tempfile

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

REFERENCE_CSV = REPO_ROOT / 'sweeps' / 'ni_np_lit_screening_row.csv'

# mode -> (reactor types, temperatures K, material-class override)
MODES = {
    'thermocatalytic_pfr': (['PFR'], [923.15, 973.15], None),
    'thermocatalytic_fluidized': (['Fluidized'], [923.15, 973.15], None),
    'thermocatalytic': (['PFR', 'Fluidized'], [973.15], None),
    'mmbcr': (['MMBCR'], [1173.15], 'MoltenMetal'),
    'ntec': (['NTEC'], [473.15], None),
    'electrochemical': (['Electrochemical'], [353.15, 823.15], None),
}


def _reference_row(material_class: str | None) -> pd.Series:
    row = pd.read_csv(REFERENCE_CSV).iloc[0].copy()
    if material_class == 'MoltenMetal':
        row['candidate_id'] = 'bini_melt_reference'
        row['material_class'] = 'MoltenMetal'
        row['genome'] = "('MoltenMetal', 'Bi', 'Ni', 10.0, 1000)"
    return row


def _run_mode(mode: str, results: Path) -> list[dict]:
    from pipeline.stages.reactor import simulate_candidate

    reactors, temperatures, material_class = MODES[mode]
    row = _reference_row(material_class)
    result = simulate_candidate(
        row, str(row['candidate_id']), temperatures, reactor_types=reactors,
        forbid_mock=True, pathway_mode=mode,
        multiphysics_results_dir=str(results / 'multiphysics'),
    )
    return result['sweep']


def _assert_complete(mode: str, sweep: list[dict]) -> None:
    reactors, temperatures, _ = MODES[mode]
    assert len(sweep) == len(reactors) * len(temperatures), (mode, len(sweep))
    for record in sweep:
        status = record.get('status')
        assert status == 'complete', (
            f"{mode}/{record.get('reactor_type')} at {record.get('T_K')} K: "
            f"{status} [{record.get('reactor_evidence_tier')}] "
            f"{record.get('screening_error') or record.get('reason') or ''}")
        conversion = float(record['CH4_conversion'])
        assert 0.0 <= conversion <= 1.0, (mode, conversion)
        assert record.get('can_exclude_candidate') is False, mode


def main() -> None:
    """Run every pathway mode and fail on any non-complete record."""
    import os

    with tempfile.TemporaryDirectory(prefix='modes_e2e_') as tmp:
        results = Path(tmp)
        os.environ.setdefault('HYDROGEN_RESULTS_DIR', str(results))
        failures = []
        for mode in MODES:
            try:
                sweep = _run_mode(mode, results)
                _assert_complete(mode, sweep)
                best = max(float(r['CH4_conversion']) for r in sweep)
                tiers = sorted({str(r.get('reactor_evidence_tier')) for r in sweep})
                print(f'PASS {mode:28s} max X={best:.4f} tiers={tiers}')
            except Exception as exc:  # report every mode, then fail
                failures.append(mode)
                print(f'FAIL {mode:28s} {type(exc).__name__}: {exc}')
        print(f'{len(MODES) - len(failures)}/{len(MODES)} pathway modes complete')
        if failures:
            raise SystemExit(1)


if __name__ == '__main__':
    main()
