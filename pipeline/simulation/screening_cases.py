"""Candidate-specific, UNCALIBRATED NTEC/electrochemical screening cases.

NTEC and electrochemical production results require experimentally
calibrated MOOSE + Cantera artifacts. Until those exist, this module lets the
campaign still run the same coupled MOOSE <-> Cantera physics per candidate
under one shared, *assumed* representative operating point
(``moose/operating_points``). The candidate-specific input is the Cantera
interface mechanism written from the candidate's screening kinetics.

Results are labelled ``uncalibrated_multiphysics_screening``, never carry a
calibration block, never become production artifacts, and never authorize
candidate exclusion (analogous to the labelled analytical closures used by
the Fluidized and MMBCR reactors).
"""

from __future__ import annotations

import json
import re
import shutil
import tempfile
from pathlib import Path

from pipeline.simulation.moose_backend import REPOSITORY_ROOT, resolve_moose_executable
from pipeline.simulation.physical_case import CASE_SCHEMA_VERSION
from pipeline.simulation.solver_handoff import require_pristine_case

OPERATING_POINT_DIR = REPOSITORY_ROOT / 'moose' / 'operating_points'
SCREENING_MODES = {'NTEC': 'ntec', 'Electrochemical': 'electrochemical'}
# Above this temperature water-based electrolytes are impractical at 1 atm
# (and the case contract rejects aqueous T >= 647.096 K); use a molten salt.
AQUEOUS_MAX_TEMPERATURE_K = 600.0
KINETICS_SOURCE = 'candidate screening kinetics (write_full_mechanism)'
TEMPERATURE_SOURCE = ('campaign reactor temperature grid '
                      '(pipeline.reactors.modes.reactor_temperatures_for_mode)')
_CASE_SECTIONS = ('geometry', 'operating', 'properties', 'model_parameters',
                  'models', 'feed')


def _operating_point_files(reactor_type: str) -> tuple[Path, ...]:
    if reactor_type == 'NTEC':
        return (OPERATING_POINT_DIR / 'ntec.json',)
    if reactor_type == 'Electrochemical':
        return tuple(OPERATING_POINT_DIR / f'electrochemical_{phase}.json'
                     for phase in ('aqueous', 'molten'))
    return ()


def screening_available(reactor_type: str) -> bool:
    """Report whether the uncalibrated screening tier can run for a reactor.

    Args:
        reactor_type: Physical reactor implementation identifier.

    Returns:
        True only for NTEC/Electrochemical when every operating-point file
        exists and a MOOSE executable resolves.
    """
    files = _operating_point_files(reactor_type)
    if not files or not all(path.is_file() for path in files):
        return False
    return resolve_moose_executable(required=False) is not None


def electrolyte_phase_for(temperature_K: float) -> str:
    """Select the electrolyte phase for an electrochemical screening point.

    A phase declared in ``ELECTROCHEMICAL_CONDITIONS_JSON`` wins; otherwise
    temperatures below 600 K use the aqueous point and higher ones the molten
    point.

    Args:
        temperature_K: Absolute operating temperature in kelvin.

    Returns:
        ``'aqueous'`` or ``'molten'``.
    """
    from pipeline.electrochemistry.model import conditions_from_environment

    declared = conditions_from_environment().electrolyte_phase
    if declared in {'aqueous', 'molten'}:
        return str(declared)
    return 'aqueous' if float(temperature_K) < AQUEOUS_MAX_TEMPERATURE_K else 'molten'


def _surface_phase_name(mechanism_text: str) -> str:
    """Return the single ideal-surface phase declared in a Cantera YAML file."""
    names = re.findall(r'^- name:\s*(\S+)\s*\n\s+thermo:\s*ideal-surface\s*$',
                       mechanism_text, flags=re.MULTILINE)
    if len(names) != 1:
        raise ValueError(
            f'mechanism must declare exactly one ideal-surface phase; found {names}')
    return names[0]


def _load_operating_point(reactor_type: str, phase: str | None) -> dict:
    files = _operating_point_files(reactor_type)
    path = files[0] if reactor_type == 'NTEC' else \
        OPERATING_POINT_DIR / f'electrochemical_{phase}.json'
    if not path.is_file():
        raise FileNotFoundError(f'operating point is missing: {path}')
    value = json.loads(path.read_text())
    if value.get('reactor_type') != reactor_type or 'temperature_K' in value.get(
            'operating', {}):
        raise ValueError(f'operating point {path.name} does not match {reactor_type} '
                         'or declares temperature_K')
    value['_path'] = path
    return value


def build_screening_case(case_dir: str | Path, *, candidate_id: str, mode: str,
                         reactor_type: str, temperature_K: float,
                         mechanism_path: str | Path) -> Path:
    """Create a pristine, uncalibrated screening case for one candidate.

    The case combines the tracked representative operating point with the
    campaign temperature and a copy of the candidate's Cantera mechanism. It
    deliberately has no ``calibration`` block, so only
    ``load_physical_case(..., require_calibration=False)`` accepts it.

    Args:
        case_dir: Directory to create; it must be absent or empty.
        candidate_id: Stable candidate identifier (a safe path component).
        mode: Pathway mode routing ``reactor_type``.
        reactor_type: ``NTEC`` or ``Electrochemical``.
        temperature_K: Absolute operating temperature in kelvin.
        mechanism_path: Candidate Cantera YAML from ``write_full_mechanism``.

    Returns:
        Path of the written ``hydrogen_case.json``.
    """
    if SCREENING_MODES.get(reactor_type) != mode:
        raise ValueError(f'{reactor_type} has no screening case for mode {mode!r}')
    if Path(str(candidate_id)).name != str(candidate_id):
        raise ValueError('candidate_id is not a safe path component')
    mechanism = Path(mechanism_path)
    if not str(mechanism_path) or not mechanism.is_file():
        raise FileNotFoundError(f'candidate mechanism file is missing: {mechanism_path}')
    text = mechanism.read_text(encoding='utf-8')
    interface_phase = _surface_phase_name(text)
    phase = electrolyte_phase_for(temperature_K) \
        if reactor_type == 'Electrochemical' else None
    point = _load_operating_point(reactor_type, phase)

    case = Path(case_dir)
    case.mkdir(parents=True, exist_ok=True)
    if any(case.iterdir()):
        raise ValueError(f'screening case directory is not empty: {case}')
    (case / mechanism.name).write_text(text, encoding='utf-8')
    value: dict = {
        'schema_version': CASE_SCHEMA_VERSION,
        'candidate_id': str(candidate_id),
        'pathway_mode': mode,
        'reactor_type': reactor_type,
    }
    for section in _CASE_SECTIONS:
        value[section] = json.loads(json.dumps(point[section]))
    value['operating']['temperature_K'] = float(temperature_K)
    value['parameter_sources'] = {
        **point['parameter_sources'], 'operating.temperature_K': TEMPERATURE_SOURCE}
    value['kinetics'] = {
        'source': KINETICS_SOURCE, 'mechanism': mechanism.name,
        'gas_phase': 'gas', 'interface_phase': interface_phase,
    }
    value['screening'] = {
        'tier': 'uncalibrated_multiphysics_screening', 'calibrated': False,
        'operating_point': point['_path'].relative_to(REPOSITORY_ROOT).as_posix(),
        'operating_point_status': point.get('status'),
    }
    if reactor_type == 'Electrochemical':
        value['electrolyte_phase'] = phase
        value['electrolyte'] = dict(point['electrolyte'])
    target = case / 'hydrogen_case.json'
    target.write_text(json.dumps(value, indent=1, sort_keys=True) + '\n')
    require_pristine_case(case)
    return target


def _cached_result(target: Path, *, candidate_id: str, mode: str,
                   reactor_type: str, temperature_K: float, digest: str) -> dict | None:
    from pipeline.simulation.external_runner import SCREENING_RESULT, SCREENING_TIER

    path = target / SCREENING_RESULT
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text())
        matches = (
            value.get('valid') is True and value.get('tier') == SCREENING_TIER
            and value.get('candidate_exclusion_authorized') is False
            and value.get('candidate_id') == str(candidate_id)
            and value.get('pathway_mode') == mode
            and value.get('reactor_type') == reactor_type
            and abs(float(value.get('temperature_K')) - float(temperature_K)) < 1e-6
            and value.get('input_sha256') == digest)
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return None
    return value if matches else None


def run_or_reuse_screening_case(results_root: str | Path, *, candidate_id: str,
                                mode: str, reactor_type: str, temperature_K: float,
                                mechanism_path: str | Path,
                                timeout_s: int = 3600, execution=None) -> dict:
    """Build and run (or reuse) the uncalibrated screening case for a condition.

    The case lives in
    ``<results_root>/screening_cases/<candidate_id>/<reactor_type>/<T:.2f>K``.
    A previous ``hydrogen_screening_result.json`` there is reused only when
    its identity and input digest match a freshly staged case; otherwise the
    directory is replaced and the solvers rerun.

    Args:
        results_root: Multiphysics results root for the campaign.
        candidate_id: Stable candidate identifier (a safe path component).
        mode: Pathway mode routing ``reactor_type``.
        reactor_type: ``NTEC`` or ``Electrochemical``.
        temperature_K: Absolute operating temperature in kelvin.
        mechanism_path: Candidate Cantera YAML from ``write_full_mechanism``.
        timeout_s: Wall-time bound for each solver process in seconds.
        execution: Optional ``SolverExecutionServices`` for tests.

    Returns:
        The labelled screening record from ``run_screening_case``, with
        ``reused`` set when a matching cached result was returned.
    """
    from pipeline.simulation.external_runner import (
        run_screening_case, screening_case_digest)

    if Path(str(candidate_id)).name != str(candidate_id):
        raise ValueError('candidate_id is not a safe path component')
    if not str(mechanism_path) or not Path(mechanism_path).is_file():
        raise FileNotFoundError(f'candidate mechanism file is missing: {mechanism_path}')
    target = (Path(results_root) / 'screening_cases' / str(candidate_id) /
              reactor_type / f'{float(temperature_K):.2f}K')
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.staging-', dir=target.parent))
    try:
        build_screening_case(
            staging, candidate_id=candidate_id, mode=mode,
            reactor_type=reactor_type, temperature_K=temperature_K,
            mechanism_path=mechanism_path)
        digest = screening_case_digest(staging, reactor_type)
        cached = _cached_result(
            target, candidate_id=candidate_id, mode=mode,
            reactor_type=reactor_type, temperature_K=temperature_K, digest=digest)
        if cached is not None:
            return {**cached, 'reused': True}
        if target.exists():
            shutil.rmtree(target)
        staging.rename(target)
    finally:
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
    result = run_screening_case(
        mode=mode, reactor_type=reactor_type, candidate_id=str(candidate_id),
        temperature_K=float(temperature_K), case_dir=target,
        timeout_s=timeout_s, execution=execution)
    return {**result, 'reused': False}
