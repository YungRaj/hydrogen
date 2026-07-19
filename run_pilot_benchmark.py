#!/usr/bin/env python3
"""Run the candidate-held-out pilot selection benchmark."""

from pipeline.pilot_benchmark import default_specs, run_pilot, write_report


def main() -> None:
    results = [run_pilot(spec) for spec in default_specs()]
    json_path, md_path = write_report(results)
    for result in results:
        print(f"{result['application']}: learned={result['learned']['hit_rate']:.1%}, "
              f"expert={result['expert']['hit_rate']:.1%}, "
              f"random={result['random']['mean_hit_rate']:.1%}, "
              f"enrichment={result['enrichment_vs_random']:.2f}x")
    print(f"Reports: {json_path}, {md_path}")


if __name__ == "__main__":
    main()
