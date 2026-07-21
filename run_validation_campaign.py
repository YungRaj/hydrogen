#!/usr/bin/env python3
"""Inspect or safely resume candidate-specific QE validation stages."""

import argparse
import json

from pipeline.validation.production_workflow import (
    methane_neb_status,
    orr_campaign_status,
    run_methane_neb,
    run_orr_sequence,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--pyro-dir')
    parser.add_argument('--pyro-prefix', default='candidate')
    parser.add_argument('--orr-dir')
    parser.add_argument('--orr-name')
    parser.add_argument('--advance', action='store_true',
                        help='Run the next eligible stages; otherwise report only')
    parser.add_argument('--restart-incomplete', action='store_true',
                        help='Overwrite incomplete ORR outputs (never use while a job is active)')
    parser.add_argument('--timeout-s', type=int, default=86400)
    args = parser.parse_args()
    if not args.pyro_dir and not args.orr_dir:
        parser.error('provide --pyro-dir and/or --orr-dir')
    if args.orr_dir and not args.orr_name:
        parser.error('--orr-name is required with --orr-dir')

    result = {}
    if args.pyro_dir:
        result['methane_neb'] = (
            run_methane_neb(args.pyro_dir, args.pyro_prefix,
                            timeout_s=args.timeout_s)
            if args.advance else methane_neb_status(args.pyro_dir)
        )
    if args.orr_dir:
        result['orr'] = (
            run_orr_sequence(args.orr_dir, args.orr_name,
                             timeout_s=args.timeout_s,
                             restart_incomplete=args.restart_incomplete)
            if args.advance else orr_campaign_status(args.orr_dir, args.orr_name)
        )
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
