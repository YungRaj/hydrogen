"""Independent evidence-report generation stage."""

from __future__ import annotations

from typing import Callable, Mapping

from pipeline.stages.contracts import StageOutcome


def run_report_stage(pipeline_state: Mapping, *,
                     generator: Callable | None = None) -> StageOutcome:
    """Generate a report through an injectable renderer boundary.

    Args:
        pipeline_state: Mapping supplying pipeline state.
        generator: Injected callable used to perform generator.

    Returns:
        Computed `StageOutcome` result.
    """
    if generator is None:
        from pipeline.evidence.report_generator import generate_full_report
        generator = generate_full_report
    report_path = generator(dict(pipeline_state))
    return StageOutcome(
        state={'report_path': str(report_path)},
        products={'report_path': report_path})
