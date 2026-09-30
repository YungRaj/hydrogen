"""Fast, fail-closed qualification before a bounded campaign is launched."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sqlite3
import subprocess
import tempfile
from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path

from pipeline.evidence.design_space_audit import audit_design_space
from pipeline.reactors.modes import reactor_types_for_mode, validate_mode_reactors
from pipeline.search.discovery import candidate_id
from pipeline.search.exhaustive_search import ScanConfig, run_streaming_scan
from pipeline.search.indexed_space import (
    CLASS_ORDER,
    deterministic_tree_probes,
    is_physically_admissible,
)
from pipeline.validation.task_queue import ValidationTaskQueue
from pipeline.utils import check_element_safety


@dataclass(frozen=True, slots=True)
class QualificationCheck:
    """One independently inspectable campaign-launch condition."""

    name: str
    passed: bool
    blocking: bool
    detail: str


@dataclass(frozen=True, slots=True)
class CampaignQualification:
    """Complete bounded-campaign qualification result."""

    schema_version: int
    mode: str
    ready: bool
    checks: tuple[QualificationCheck, ...]


Preflight = Callable[[str], Mapping[str, object]]


@dataclass(frozen=True, slots=True)
class HostResources:
    """Measured resources required before campaign work is admitted."""

    writable: bool
    free_disk_gb: float
    gpu_count: int


@dataclass(frozen=True, slots=True)
class CandidateTableIdentity:
    """Typed identity fields read from one external candidate-table row."""

    candidate_id: str
    material_class: str


ResourceProbe = Callable[[Path], HostResources]


def probe_host_resources(results_dir: Path) -> HostResources:
    """Measure writable storage and visible NVIDIA devices without workloads.

    Args:
        results_dir: Campaign output directory or its nearest existing parent.

    Returns:
        Typed snapshot of writable storage and visible GPU count.
    """
    target = results_dir.expanduser().resolve()
    existing = target
    while not existing.exists() and existing != existing.parent:
        existing = existing.parent
    free_disk_gb = shutil.disk_usage(existing).free / 1_000_000_000
    writable = False
    try:
        target.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=target, prefix=".qualification-", delete=True
        ):
            writable = True
    except OSError:
        writable = False
    try:
        completed = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        gpu_count = (
            len([line for line in completed.stdout.splitlines() if line.strip()])
            if completed.returncode == 0
            else 0
        )
    except (OSError, subprocess.TimeoutExpired):
        gpu_count = 0
    return HostResources(writable, free_disk_gb, gpu_count)


def _resource_check(
    results_dir: Path,
    minimum_free_disk_gb: float,
    resource_probe: ResourceProbe,
) -> QualificationCheck:
    resources = resource_probe(results_dir)
    passed = (
        resources.writable
        and resources.free_disk_gb >= minimum_free_disk_gb
        and resources.gpu_count > 0
    )
    return QualificationCheck(
        "host_resources",
        passed,
        True,
        (
            f"writable={resources.writable}; free_disk_gb={resources.free_disk_gb:.1f}; "
            f"visible_nvidia_gpus={resources.gpu_count}; "
            f"minimum_free_disk_gb={minimum_free_disk_gb:.1f}"
        ),
    )


def _design_space_check(samples_per_class: int) -> QualificationCheck:
    audit = audit_design_space(sample_per_class=samples_per_class)
    return QualificationCheck(
        name="design_space",
        passed=bool(audit["valid"])
        and audit["classes_represented"] == len(CLASS_ORDER),
        blocking=True,
        detail=(
            f"{audit['classes_represented']}/{len(CLASS_ORDER)} classes represented; "
            f"{len(audit['failures'])} failure(s)"
        ),
    )


def _identity_control_check() -> QualificationCheck:
    probes = deterministic_tree_probes(len(CLASS_ORDER) * 2)
    classes = {genome[0] for genome in probes}
    identities = [candidate_id(genome) for genome in probes]
    unsafe_control, _ = check_element_safety(["Hg"])
    malformed_control = ("Perovskite", "La", "Fe", "Fe", 0.1, "O_vac")
    malformed_admissible, _ = is_physically_admissible(malformed_control)
    passed = (
        classes == set(CLASS_ORDER)
        and len(identities) == len(set(identities))
        and not unsafe_control
        and not malformed_admissible
    )
    return QualificationCheck(
        name="identity_and_controls",
        passed=passed,
        blocking=True,
        detail=(
            f"{len(classes)}/{len(CLASS_ORDER)} classes; "
            f"{len(identities)} unique deterministic identities; "
            f"unsafe control rejected={not unsafe_control}; "
            f"malformed control rejected={not malformed_admissible}"
        ),
    )


def _routing_check(mode: str) -> QualificationCheck:
    try:
        reactors = reactor_types_for_mode(mode)
        validate_mode_reactors(mode, reactors)
    except ValueError as exc:
        return QualificationCheck("reactor_routing", False, True, str(exc))
    return QualificationCheck(
        "reactor_routing", True, True, f"{mode}: {', '.join(reactors)}"
    )


def _solver_check(mode: str, preflight: Preflight) -> QualificationCheck:
    try:
        result = preflight(mode)
    except Exception as exc:
        return QualificationCheck(
            "solver_preflight",
            False,
            True,
            f"preflight failed: {type(exc).__name__}: {exc}",
        )
    missing_value = result.get("missing")
    missing = (
        tuple(str(value) for value in missing_value)
        if isinstance(missing_value, (list, tuple))
        else ("invalid_preflight_result",)
    )
    return QualificationCheck(
        name="solver_preflight",
        passed=not missing,
        blocking=True,
        detail=(
            "all required solvers available"
            if not missing
            else "missing: " + ", ".join(missing)
        ),
    )


def _physical_case_check(mode: str, manifest: Path | None) -> QualificationCheck:
    external_reactors = {"Fluidized", "MMBCR", "NTEC", "Electrochemical"}
    try:
        required = external_reactors.intersection(reactor_types_for_mode(mode))
    except ValueError as exc:
        return QualificationCheck("physical_case_inputs", False, True, str(exc))
    if not required:
        return QualificationCheck(
            "physical_case_inputs", True, True, "not required for PFR mode"
        )
    if manifest is None or not manifest.is_file():
        return QualificationCheck(
            "physical_case_inputs",
            False,
            True,
            "validated case manifest required for: " + ", ".join(sorted(required)),
        )
    try:
        from pipeline.simulation.case_preparation import prepare_manifest

        report = prepare_manifest(manifest, create=False)
        cases = report.get("cases", [])
        covered = {
            str(case.get("reactor_type"))
            for case in cases
            if case.get("mode") == mode and case.get("ready") is True
        }
        passed = report.get("not_ready") == 0 and required.issubset(covered)
        detail = (
            f"ready_cases={report.get('ready', 0)}/{report.get('total', 0)}; "
            f"covered={sorted(covered)}; required={sorted(required)}"
        )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        passed = False
        detail = f"case manifest invalid: {exc}"
    return QualificationCheck("physical_case_inputs", passed, True, detail)


def _resume_check() -> QualificationCheck:
    genome = ("SAC", "Fe", "N4", "N-graphene", "none")
    with tempfile.TemporaryDirectory(prefix="hydrogen-qualification-") as temporary:
        queue = ValidationTaskQueue(Path(temporary) / "tasks.sqlite")
        identity = queue.enqueue("qualification", genome, "dft", "canary-v1")
        first_claim = queue.claim("qualification", identity, "dft", "canary-v1")
        recovered = queue.recover_stale(stale_after_s=0)
        second_claim = queue.claim("qualification", identity, "dft", "canary-v1")
        queue.finish("qualification", identity, "dft", "canary-v1", True, "canary")
        summary = queue.summary("qualification", "dft")
    passed = (
        first_claim and recovered == 1 and second_claim and summary == {"converged": 1}
    )
    return QualificationCheck(
        "resume_and_recovery",
        passed,
        True,
        f"first_claim={first_claim}; recovered={recovered}; final={summary}",
    )


def _search_resume_check() -> QualificationCheck:
    import numpy as np

    def scorer(genomes: list[tuple]) -> np.ndarray:
        return np.asarray(
            [
                [
                    float(int(candidate_id(genome)[:8], 16)),
                    float(int(candidate_id(genome)[8:16], 16)),
                ]
                for genome in genomes
            ]
        )

    with tempfile.TemporaryDirectory(
        prefix="hydrogen-search-qualification-"
    ) as temporary:
        interrupted = Path(temporary) / "interrupted.sqlite"
        reference = Path(temporary) / "reference.sqlite"
        common = {
            "application": "qualification",
            "start": 0,
            "stop": 64,
            "batch_size": 16,
            "global_archive_size": 16,
        }
        run_streaming_scan(
            ScanConfig(database=str(interrupted), max_batches=1, **common), scorer
        )
        resumed = run_streaming_scan(
            ScanConfig(database=str(interrupted), **common), scorer
        )
        direct = run_streaming_scan(
            ScanConfig(database=str(reference), **common), scorer
        )

        def archive(path: Path) -> list[tuple[str, float]]:
            with sqlite3.connect(path) as connection:
                return connection.execute(
                    "SELECT candidate_id, primary_score FROM global_archive "
                    "WHERE application='qualification' ORDER BY candidate_id"
                ).fetchall()

        archives_match = archive(interrupted) == archive(reference)
    passed = (
        resumed["next_index"] == direct["next_index"] == 64
        and resumed["processed_this_run"] + 16 == direct["processed_this_run"]
        and archives_match
    )
    return QualificationCheck(
        "search_resume_equivalence",
        passed,
        True,
        f"cursor={resumed['next_index']}; archive_match={archives_match}",
    )


def _typed_cache_check() -> QualificationCheck:
    from pipeline.screening.result_cache import (
        load_cached_results,
        store_completed_results,
    )

    genome = ("SAC", "Fe", "N4", "N-graphene", "none")
    protocol = "qualification-cache-v1"
    result = {
        "screening_protocol": protocol,
        "valid": True,
        "relax_clean_converged": True,
        "energy_eV": -1.25,
    }
    with tempfile.TemporaryDirectory(
        prefix="hydrogen-cache-qualification-"
    ) as temporary:
        path = Path(temporary) / "cache.sqlite"
        stored = store_completed_results(
            path, "qualification", protocol, "implementation-a", [genome], [result]
        )
        restored = load_cached_results(
            path, "qualification", protocol, "implementation-a", [genome]
        )
        changed_protocol = load_cached_results(
            path,
            "qualification",
            "qualification-cache-v2",
            "implementation-a",
            [genome],
        )
    passed = stored == 1 and restored.get(0) == result and not changed_protocol
    return QualificationCheck(
        "typed_cache_identity",
        passed,
        True,
        f"stored={stored}; exact_round_trip={restored.get(0) == result}; invalidated={not changed_protocol}",
    )


def _scientific_canary_check() -> QualificationCheck:
    from pipeline.reactors.equilibrium import ch4_conversion_from_argon_tracer
    from pipeline.utils import arrhenius_rate, orr_overpotential

    prefactor = 1.0e13
    cold = arrhenius_rate(prefactor, 0.8, 800.0)
    hot = arrhenius_rate(prefactor, 0.8, 900.0)
    higher_barrier = arrhenius_rate(prefactor, 0.9, 900.0)
    zero_barrier = arrhenius_rate(prefactor, 0.0, 900.0)
    overpotential, _ = orr_overpotential(1.23, 2.46, 3.69)
    conversion = ch4_conversion_from_argon_tracer(0.45, 0.10, 0.90, 0.10)
    passed = (
        0.0 < cold < hot < prefactor
        and higher_barrier < hot
        and zero_barrier == prefactor
        and abs(overpotential) < 1.0e-12
        and abs(conversion - 0.5) < 1.0e-12
    )
    return QualificationCheck(
        "scientific_equation_canaries",
        passed,
        True,
        (
            f"arrhenius_temperature_monotonic={hot > cold}; "
            f"arrhenius_barrier_monotonic={higher_barrier < hot}; "
            f"ideal_orr_V={overpotential:.3g}; tracer_conversion={conversion:.3g}"
        ),
    )


def _candidate_table_check(path: Path | None) -> QualificationCheck:
    if path is None or not path.is_file():
        return QualificationCheck(
            "candidate_table", True, False, "not present; discovery will create it"
        )
    try:
        with path.open(newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, None)
            malformed_rows = 0
            if header is None:
                rows: list[CandidateTableIdentity] = []
                missing_columns = {"candidate_id", "material_class"}
            else:
                missing_columns = {"candidate_id", "material_class"}.difference(header)
                rows = []
                if not missing_columns:
                    candidate_index = header.index("candidate_id")
                    class_index = header.index("material_class")
                    for row in reader:
                        if len(row) != len(header):
                            malformed_rows += 1
                            continue
                        rows.append(
                            CandidateTableIdentity(
                                row[candidate_index], row[class_index]
                            )
                        )
    except (csv.Error, OSError, UnicodeError) as exc:
        return QualificationCheck(
            "candidate_table",
            False,
            True,
            f"candidate table unreadable: {type(exc).__name__}: {exc}",
        )
    identities = [row.candidate_id for row in rows]
    classes = {row.material_class for row in rows}
    passed = (
        bool(rows)
        and not missing_columns
        and malformed_rows == 0
        and all(identities)
        and len(identities) == len(set(identities))
        and classes.issubset(set(CLASS_ORDER))
    )
    return QualificationCheck(
        "candidate_table",
        passed,
        True,
        f"rows={len(rows)}; classes={len(classes)}; malformed_rows={malformed_rows}; "
        f"missing_columns={sorted(missing_columns)}",
    )


def qualify_campaign(
    *,
    mode: str = "thermocatalytic_pfr",
    candidate_table: str | Path | None = None,
    case_manifest: str | Path | None = None,
    results_dir: str | Path = "results",
    minimum_free_disk_gb: float = 10.0,
    samples_per_class: int = 64,
    preflight: Preflight | None = None,
    resource_probe: ResourceProbe | None = None,
) -> CampaignQualification:
    """Run fast launch checks without executing scientific solvers.

    Args:
        mode: Methane-conversion pathway to qualify.
        candidate_table: Optional existing screening-table artifact.
        case_manifest: Optional external-physics case manifest for non-PFR modes.
        results_dir: Campaign output directory whose storage is checked.
        minimum_free_disk_gb: Minimum free storage admitted for a bounded campaign.
        samples_per_class: Deterministic design-space audit sample per class.
        preflight: Injectable solver-availability check used by tests and adapters.
        resource_probe: Injectable host-resource probe used by tests and adapters.

    Returns:
        Typed go/no-go result containing every independently visible check.
    """
    if samples_per_class < 3:
        raise ValueError("samples_per_class must be at least 3")
    if minimum_free_disk_gb <= 0:
        raise ValueError("minimum_free_disk_gb must be positive")
    if preflight is None:
        from pipeline.simulation.result_contract import mode_preflight

        preflight = mode_preflight
    if resource_probe is None:
        resource_probe = probe_host_resources
    checks = (
        _resource_check(Path(results_dir), minimum_free_disk_gb, resource_probe),
        _design_space_check(samples_per_class),
        _identity_control_check(),
        _routing_check(mode),
        _solver_check(mode, preflight),
        _physical_case_check(mode, Path(case_manifest) if case_manifest else None),
        _resume_check(),
        _search_resume_check(),
        _typed_cache_check(),
        _scientific_canary_check(),
        _candidate_table_check(Path(candidate_table) if candidate_table else None),
    )
    ready = all(check.passed for check in checks if check.blocking)
    return CampaignQualification(1, mode, ready, checks)


def main() -> None:
    """Run qualification and emit an external JSON report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", default="thermocatalytic_pfr")
    parser.add_argument("--candidate-table")
    parser.add_argument("--case-manifest")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--minimum-free-disk-gb", type=float, default=10.0)
    parser.add_argument("--samples-per-class", type=int, default=64)
    parser.add_argument("--output")
    args = parser.parse_args()
    report = qualify_campaign(
        mode=args.mode,
        candidate_table=args.candidate_table,
        case_manifest=args.case_manifest,
        results_dir=args.results_dir,
        minimum_free_disk_gb=args.minimum_free_disk_gb,
        samples_per_class=args.samples_per_class,
    )
    rendered = json.dumps(asdict(report), indent=2, sort_keys=True) + "\n"
    if args.output:
        target = Path(args.output)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(rendered)
    print(rendered, end="")
    raise SystemExit(0 if report.ready else 2)


if __name__ == "__main__":
    main()
