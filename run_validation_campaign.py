#!/usr/bin/env python3
"""Inspect or safely resume candidate-specific QE validation stages."""

import argparse
import json
from dataclasses import fields
from pathlib import Path

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--pyro-dir')
    parser.add_argument('--pyro-prefix', default='candidate')
    parser.add_argument('--pyro-manifest',
                        help='JSON describing explicit geometries for all methane elementary steps')
    parser.add_argument('--prepare', action='store_true',
                        help='Prepare inputs declared by --pyro-manifest')
    parser.add_argument('--orr-dir')
    parser.add_argument('--orr-name')
    parser.add_argument('--orr-structure',
                        help='ASE-readable surface used to emit a site/coverage validation plan')
    parser.add_argument('--orr-coverages', default='0.25,0.5,1.0')
    parser.add_argument('--orr-plan-output', default='orr_validation_plan.json')
    parser.add_argument('--orr-ensemble-results',
                        help='JSON list of converged site/coverage/adsorbate dG records')
    parser.add_argument('--orr-corrections-json',
                        help='JSON list of provenance-bearing ORR correction models')
    parser.add_argument('--orr-expected-cases', type=int)
    parser.add_argument('--advance', action='store_true',
                        help='Run the next eligible stages; otherwise report only')
    parser.add_argument('--restart-incomplete', action='store_true',
                        help='Overwrite incomplete ORR outputs (never use while a job is active)')
    parser.add_argument('--timeout-s', type=int, default=86400)
    parser.add_argument('--mpi-ranks', type=int, default=4)
    parser.add_argument('--omp-threads', type=int, default=1)
    parser.add_argument('--kpoint-pools', type=int, default=1)
    parser.add_argument('--neb-image-groups', type=int, default=1)
    args = parser.parse_args()
    if not any((args.pyro_dir, args.pyro_manifest, args.orr_dir,
                args.orr_structure, args.orr_ensemble_results)):
        parser.error('provide a pyrolysis or ORR validation input')
    if args.orr_dir and not args.orr_name:
        parser.error('--orr-name is required with --orr-dir')

    # Keep argument discovery usable before the optional scientific stack is
    # installed. Actual validation operations still import their dependencies
    # eagerly here and fail with a useful missing-package error when needed.
    from pipeline.validation.qe_workflows import QEExecutionConfig
    from pipeline.validation.orr_workflows import (
        ORRCorrections, build_orr_validation_plan, evaluate_orr_ensemble)
    from pipeline.validation.production_workflow import (
        advance_methane_neb,
        advance_pyrolysis_campaign,
        methane_neb_status,
        orr_campaign_status,
        prepare_pyrolysis_campaign,
        pyrolysis_campaign_status,
        run_orr_sequence,
    )
    execution = QEExecutionConfig(
        mpi_ranks=args.mpi_ranks,
        omp_threads=args.omp_threads,
        kpoint_pools=args.kpoint_pools,
        image_groups=args.neb_image_groups,
    )

    result = {}
    if args.pyro_manifest:
        if args.prepare:
            result['pyrolysis_campaign_preparation'] = prepare_pyrolysis_campaign(
                args.pyro_manifest)
        result['pyrolysis_campaign'] = (
            advance_pyrolysis_campaign(
                args.pyro_manifest, timeout_s=args.timeout_s,
                execution=execution,
                restart_incomplete=args.restart_incomplete)
            if args.advance else pyrolysis_campaign_status(args.pyro_manifest)
        )
    if args.pyro_dir:
        result['methane_neb'] = (
            advance_methane_neb(
                args.pyro_dir, args.pyro_prefix,
                timeout_s=args.timeout_s, execution=execution,
                restart_incomplete=args.restart_incomplete)
            if args.advance else methane_neb_status(args.pyro_dir)
        )
    if args.orr_dir:
        result['orr'] = (
            run_orr_sequence(args.orr_dir, args.orr_name,
                             timeout_s=args.timeout_s,
                             restart_incomplete=args.restart_incomplete,
                             execution=QEExecutionConfig(
                                 mpi_ranks=args.mpi_ranks,
                                 omp_threads=args.omp_threads,
                                 kpoint_pools=args.kpoint_pools,
                             ))
            if args.advance else orr_campaign_status(args.orr_dir, args.orr_name)
        )
    if args.orr_structure:
        from ase.io import read as ase_read

        coverages = tuple(float(value) for value in args.orr_coverages.split(','))
        plan = build_orr_validation_plan(
            ase_read(args.orr_structure), coverages=coverages)
        target = Path(args.orr_plan_output)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(plan, indent=2, sort_keys=True))
        result['orr_validation_plan'] = {
            'path': str(target), 'tasks': len(plan),
            'sites': len({item['site_id'] for item in plan}),
            'coverages': sorted({item['coverage_ML'] for item in plan}),
        }
    if args.orr_ensemble_results:
        if not args.orr_corrections_json:
            parser.error('--orr-corrections-json is required with --orr-ensemble-results')
        rows = json.loads(Path(args.orr_ensemble_results).read_text())
        raw_corrections = json.loads(Path(args.orr_corrections_json).read_text())
        allowed = {field.name for field in fields(ORRCorrections)}
        corrections = [ORRCorrections(**{
            key: value for key, value in record.items() if key in allowed})
            for record in raw_corrections]
        result['orr_ensemble'] = evaluate_orr_ensemble(
            rows, corrections, expected_cases=args.orr_expected_cases)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
