#!/usr/bin/env python3
"""Run a reactor-cell sweep from a YAML spec.

    python runsweep.py sweeps/headline_cat9_1300K.yaml

Schema: docs/sweep-template.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pipeline.process.yaml_sweep import run_sweep


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Run a Cantera reactor-cell sweep from YAML')
    parser.add_argument(
        'input_yaml', type=Path,
        help='Sweep specification (see docs/sweep-template.md)')
    args = parser.parse_args()
    run_sweep(args.input_yaml)


if __name__ == '__main__':
    main()
