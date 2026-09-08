"""Generate case skeletons and report batch multiphysics readiness."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pipeline.process.multiphysics_contract import mode_preflight
from pipeline.process.pathway_modes import reactor_types_for_mode
from pipeline.process.physical_case import case_template, load_physical_case


def prepare_manifest(manifest_path: str | Path, *, create: bool = False) -> dict:
    manifest_file = Path(manifest_path).expanduser().resolve()
    manifest = json.loads(manifest_file.read_text())
    entries = manifest.get('cases', [])
    results = []
    for entry in entries:
        candidate = str(entry['candidate_id'])
        mode = str(entry['mode'])
        temperature = float(entry['temperature_K'])
        reactor = str(entry.get('reactor_type') or reactor_types_for_mode(mode)[0])
        directory = Path(entry['case_dir']).expanduser()
        if not directory.is_absolute():
            directory = manifest_file.parent / directory
        case_path = directory / 'hydrogen_case.json'
        if create and not case_path.exists():
            directory.mkdir(parents=True, exist_ok=True)
            case_path.write_text(json.dumps(case_template(
                candidate_id=candidate, mode=mode, reactor_type=reactor,
                temperature_K=temperature,
                electrolyte_phase=entry.get('electrolyte_phase')), indent=2) + '\n')
        failures = []
        preflight = mode_preflight(mode)
        if preflight['missing']:
            failures.append('missing_solvers:' + ','.join(preflight['missing']))
        try:
            load_physical_case(
                case_path, candidate_id=candidate, mode=mode,
                reactor_type=reactor, temperature_K=temperature)
        except (ValueError, OSError) as exc:
            failures.append(str(exc))
        results.append({
            'candidate_id': candidate, 'mode': mode, 'reactor_type': reactor,
            'temperature_K': temperature, 'case_dir': str(directory),
            'ready': not failures, 'failures': failures,
        })
    return {
        'total': len(results),
        'ready': sum(row['ready'] for row in results),
        'not_ready': sum(not row['ready'] for row in results),
        'cases': results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest')
    parser.add_argument('--create', action='store_true',
                        help='create missing non-runnable case templates')
    parser.add_argument('--output')
    args = parser.parse_args()
    report = prepare_manifest(args.manifest, create=args.create)
    rendered = json.dumps(report, indent=2) + '\n'
    if args.output:
        Path(args.output).write_text(rendered)
    print(rendered, end='')
    raise SystemExit(0 if report['not_ready'] == 0 else 2)


if __name__ == '__main__':
    main()
