"""Load a reactor-cell sweep from YAML and run it.

Input specs live under sweeps/*.yaml (git-tracked). Run products go to
results/sweeps/<name>/ (gitignored). See docs/sweep-template.md.

This is not a Cantera mechanism file. Root keys are name / catalyst /
conditions / cells. Generated mechanisms live in each job's output directory.
"""

from __future__ import annotations

import itertools
import json
import shutil
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, List, Optional

from pipeline.utils import BASE_DIR, SWEEPS_DIR, repo_relative, setup_logger

logger = setup_logger('yaml_sweep', 'reactor/yaml_sweep.log')

ALLOWED_REACTORS = ('PFR', 'Fluidized', 'MMBCR')
# catalyst.kinetics keys (kinetics-only sweeps). Adsorption energies use the
# screener's conventions (dE_H vs 1/2 H2, dE_CH3 vs CH3 radical, dE_C vs
# CH4 - 2 H2); the writer puts them on Cantera's absolute scale. The three
# B6 keys are only written for gated (Ni/Fe/Co nanoparticle) genomes.
KINETICS_KEYS = (
    'E_act', 'dE_H', 'dE_CH3', 'dE_C',
    'carbon_transfer_eV', 'carbon_encapsulation_eV',
    'encapsulation_crossover_coverage', 'carbon_transfer_prefactor_1_s',
    'carbon_transfer_particle_nm', 'carbon_diffusion_prefactor_m2_s',
    'ch4_sticking_coefficient',
)
KINETICS_FIELD_FOR_KEY = {
    'E_act': 'methane_activation_eV',
    'dE_H': 'h_adsorption_eV',
    'dE_CH3': 'ch3_adsorption_eV',
    'dE_C': 'c_adsorption_eV',
    'carbon_transfer_eV': 'carbon_transfer_eV',
    'carbon_encapsulation_eV': 'carbon_encapsulation_eV',
    'encapsulation_crossover_coverage': 'encapsulation_crossover_coverage',
    'carbon_transfer_prefactor_1_s': 'carbon_transfer_prefactor_1_s',
    'carbon_transfer_particle_nm': 'carbon_transfer_particle_nm',
    'carbon_diffusion_prefactor_m2_s': 'carbon_diffusion_prefactor_m2_s',
    'ch4_sticking_coefficient': 'ch4_sticking_coefficient',
}
DEFAULT_POLICY = {
    'co2_permitted': False,
    'fluidized_mode': 'circulating',
    'max_regen_cycles': 3,
    'regen_mechanism': 'mechanical',
}
# `sweep:` block. Kinetics keys take a list of numbers; policy keys take a
# list of values of the policy field's type. The grid is the cartesian
# product over every listed key. co2_permitted is deliberately not
# sweepable (turquoise policy, not a design variable).
SWEEPABLE_POLICY_KEYS = ('max_regen_cycles', 'fluidized_mode', 'regen_mechanism')


@dataclass
class SweepCell:
    """Define one named reactor-geometry and catalyst-inventory sweep cell.
    """
    name: str
    catalyst_particle_mm: float
    metal_loading: float
    metal_dispersion: float


@dataclass
class SweepJob:
    """Define a validated YAML sweep job and its operating conditions.
    """
    name: str
    description: str
    catalyst_name: str
    temperatures_K: List[float]
    reactor_types: List[str]
    policy: dict
    cells: List[SweepCell]
    source_file: str
    screening_csv: Optional[str] = None
    screening_index: Optional[int] = None
    kinetics: dict = field(default_factory=dict)
    # Explicit class / genome for kinetics-only sweeps (B6 gate, reactor
    # applicability). Screening sweeps take both from the CSV row.
    material_class: Optional[str] = None
    genome: Optional[str] = None
    # `sweep:` block: key -> list of values. Empty when the spec is a single
    # point. Kinetics keys write one mechanism per grid point.
    sweep_kinetics: dict = field(default_factory=dict)
    sweep_policy: dict = field(default_factory=dict)


@dataclass(frozen=True)
class GridPoint:
    """Bind one deterministic kinetic/policy combination to its index.
    """
    index: int
    kinetics: dict
    policy: dict

    @property
    def values(self) -> dict:
        """Merge the kinetic and policy values represented by this grid point.

        Returns:
            Validated dict output for this operation.
        """
        return {**self.kinetics, **self.policy}


def grid_points(job: SweepJob) -> List[GridPoint]:
    """Cartesian product of the `sweep:` lists; one point for a single spec.

    Args:
        job: Input controlling job.

    Returns:
        Validated List[GridPoint] output for this operation.
    """
    keys = list(job.sweep_kinetics) + list(job.sweep_policy)
    if not keys:
        return [GridPoint(0, {}, {})]
    lists = [job.sweep_kinetics.get(k, job.sweep_policy.get(k)) for k in keys]
    points = []
    for i, combo in enumerate(itertools.product(*lists)):
        chosen = dict(zip(keys, combo))
        points.append(GridPoint(
            i,
            {k: chosen[k] for k in job.sweep_kinetics},
            {k: chosen[k] for k in job.sweep_policy}))
    return points


def _require_mapping(value: Any, label: str) -> dict:
    if not isinstance(value, dict):
        raise ValueError(f'{label} must be a mapping')
    return value


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        text = value.strip().lower()
        if text in ('true', '1', 'yes'):
            return True
        if text in ('false', '0', 'no'):
            return False
    raise ValueError(f'expected boolean, got {value!r}')


def _as_float_list(value: Any, label: str) -> List[float]:
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, str):
        parts = [p for p in value.replace(',', ' ').split() if p]
        if not parts:
            raise ValueError(f'{label}: expected at least one number')
        return [float(p) for p in parts]
    if isinstance(value, list):
        if not value:
            raise ValueError(f'{label}: expected at least one number')
        return [float(v) for v in value]
    raise ValueError(f'{label} must be a number or a list of numbers')


def _as_name_list(value: Any, label: str) -> List[str]:
    if isinstance(value, str):
        names = [p for p in value.replace(',', ' ').split() if p]
    elif isinstance(value, list):
        names = [str(v).strip() for v in value if str(v).strip()]
    else:
        raise ValueError(f'{label} must be a string or a list of names')
    if not names:
        raise ValueError(f'{label}: expected at least one name')
    return names


def _resolve(path_text: str) -> Path:
    path = Path(path_text)
    if not path.is_absolute():
        path = BASE_DIR / path
    return path


def _load_yaml(path: Path) -> dict:
    text = path.read_text(encoding='utf-8')
    try:
        from ruamel.yaml import YAML
        data = YAML(typ='safe', pure=True).load(text)
    except ImportError:
        try:
            import yaml
        except ImportError as exc:
            raise ImportError(
                'sweep specs need ruamel.yaml (Cantera) or PyYAML') from exc
        data = yaml.safe_load(text)
    return _require_mapping(data, 'sweep YAML')


def parse_sweep(yaml_path: Path) -> SweepJob:
    """Parse one sweep spec. Does not touch Cantera or write mechanisms.

    Args:
        yaml_path: Input controlling yaml path.

    Returns:
        Validated SweepJob output for this operation.
    """
    yaml_path = Path(yaml_path)
    if yaml_path.suffix.lower() == '.xml':
        raise ValueError(
            'sweep specs are YAML; convert the file and see docs/sweep-template.md')
    if not yaml_path.is_file():
        raise FileNotFoundError(f'no such sweep file: {yaml_path}')
    root = _load_yaml(yaml_path)
    name = str(root.get('name') or yaml_path.stem).strip()
    if not name:
        raise ValueError('sweep name is required')

    catalyst = _require_mapping(root.get('catalyst'), 'catalyst')
    catalyst_name = str(catalyst.get('name') or name).strip()
    screening = catalyst.get('screening')
    kinetics_raw = catalyst.get('kinetics')
    screening_csv = None
    screening_index = None
    kinetics = {}
    if screening is not None:
        screening = _require_mapping(screening, 'catalyst.screening')
        csv_text = screening.get('csv')
        index_val = screening.get('index')
        if not csv_text or index_val is None:
            raise ValueError('catalyst.screening requires csv and index')
        screening_csv = str(csv_text)
        screening_index = int(index_val)
    if kinetics_raw is not None:
        kinetics_raw = _require_mapping(kinetics_raw, 'catalyst.kinetics')
        for key in KINETICS_KEYS:
            raw = kinetics_raw.get(key)
            if raw is not None and raw != '':
                kinetics[key] = float(raw)
        unknown_keys = set(kinetics_raw) - set(KINETICS_KEYS) - {'provenance'}
        if unknown_keys:
            raise ValueError(
                f'unknown catalyst.kinetics key(s) {sorted(unknown_keys)}; '
                f'allowed {list(KINETICS_KEYS)}')
        provenance = kinetics_raw.get('provenance')
        if provenance is not None:
            kinetics['provenance'] = {
                str(k): str(v) for k, v in _require_mapping(
                    provenance, 'catalyst.kinetics.provenance').items()}
    sweep_kinetics: dict = {}
    sweep_policy: dict = {}
    sweep_raw = root.get('sweep')
    if sweep_raw is not None:
        sweep_raw = _require_mapping(sweep_raw, 'sweep')
        unknown_blocks = set(sweep_raw) - {'kinetics', 'policy'}
        if unknown_blocks:
            raise ValueError(
                f'unknown sweep block(s) {sorted(unknown_blocks)}; '
                "allowed ['kinetics', 'policy']")
        kin_block = sweep_raw.get('kinetics')
        if kin_block is not None:
            kin_block = _require_mapping(kin_block, 'sweep.kinetics')
            for key, raw in kin_block.items():
                if key not in KINETICS_KEYS:
                    raise ValueError(
                        f'unknown sweep.kinetics key {key!r}; allowed {list(KINETICS_KEYS)}')
                if key in kinetics:
                    raise ValueError(
                        f'{key} is fixed in catalyst.kinetics and swept in '
                        'sweep.kinetics; choose one')
                sweep_kinetics[key] = _as_float_list(raw, f'sweep.kinetics.{key}')
            if screening is not None and sweep_kinetics:
                raise ValueError(
                    'sweep.kinetics needs a kinetics-only catalyst (no screening row)')
        pol_block = sweep_raw.get('policy')
        if pol_block is not None:
            pol_block = _require_mapping(pol_block, 'sweep.policy')
            for key, raw in pol_block.items():
                if key not in SWEEPABLE_POLICY_KEYS:
                    raise ValueError(
                        f'unknown or non-sweepable sweep.policy key {key!r}; '
                        f'allowed {list(SWEEPABLE_POLICY_KEYS)}')
                if not isinstance(raw, list) or not raw:
                    raise ValueError(f'sweep.policy.{key} must be a non-empty list')
                if key == 'max_regen_cycles':
                    sweep_policy[key] = [int(v) for v in raw]
                else:
                    sweep_policy[key] = [str(v).strip() for v in raw]
    if screening is None and 'E_act' not in kinetics and 'E_act' not in sweep_kinetics:
        raise ValueError(
            'catalyst needs screening: {csv, index} or kinetics: {E_act}')
    material_class = catalyst.get('material_class')
    material_class = str(material_class).strip() if material_class else None
    genome = catalyst.get('genome')
    genome = str(genome).strip() if genome else None
    if screening is None and material_class is None:
        raise ValueError(
            'kinetics-only sweeps need catalyst.material_class (reactor '
            'applicability and the B6 gate depend on it)')

    conditions = _require_mapping(root.get('conditions'), 'conditions')
    if 'temperatures_K' in conditions:
        temperatures = _as_float_list(conditions['temperatures_K'], 'temperatures_K')
    elif 'temperatures' in conditions:
        temperatures = _as_float_list(conditions['temperatures'], 'temperatures')
    else:
        raise ValueError('conditions.temperatures_K is required')
    reactor_types = _as_name_list(conditions.get('reactors'), 'reactors')
    unknown = [r for r in reactor_types if r not in ALLOWED_REACTORS]
    if unknown:
        raise ValueError(f'unknown reactor type(s) {unknown}; allowed {ALLOWED_REACTORS}')

    policy = dict(DEFAULT_POLICY)
    policy_raw = root.get('policy')
    if policy_raw is not None:
        policy_raw = _require_mapping(policy_raw, 'policy')
        if 'co2_permitted' in policy_raw:
            policy['co2_permitted'] = _as_bool(policy_raw['co2_permitted'])
        if 'fluidized_mode' in policy_raw:
            policy['fluidized_mode'] = str(policy_raw['fluidized_mode']).strip()
        if 'max_regen_cycles' in policy_raw:
            policy['max_regen_cycles'] = int(policy_raw['max_regen_cycles'])
        if 'regen_mechanism' in policy_raw:
            policy['regen_mechanism'] = str(policy_raw['regen_mechanism']).strip()

    cells_raw = root.get('cells')
    if not isinstance(cells_raw, list) or not cells_raw:
        raise ValueError('cells must be a non-empty list')
    cells = []
    for i, cell_raw in enumerate(cells_raw):
        cell_raw = _require_mapping(cell_raw, f'cells[{i}]')
        cell_name = str(cell_raw.get('name') or f'cell_{i}').strip()
        try:
            d_p = float(cell_raw['catalyst_particle_mm'])
            loading = float(cell_raw['metal_loading'])
            dispersion = float(cell_raw['metal_dispersion'])
        except KeyError as exc:
            raise ValueError(
                f'{cell_name}: missing {exc.args[0]}') from exc
        if d_p <= 0:
            raise ValueError(f'{cell_name}: catalyst_particle_mm must be positive')
        if not (0.0 < loading <= 1.0 and 0.0 < dispersion <= 1.0):
            raise ValueError(
                f'{cell_name}: metal_loading and metal_dispersion must be in (0, 1]')
        cells.append(SweepCell(cell_name, d_p, loading, dispersion))

    description = root.get('description') or ''
    if not isinstance(description, str):
        raise ValueError('description must be a string')
    return SweepJob(
        name=name,
        description=' '.join(description.split()),
        catalyst_name=catalyst_name,
        temperatures_K=temperatures,
        reactor_types=reactor_types,
        policy=policy,
        cells=cells,
        source_file=str(yaml_path.resolve()),
        screening_csv=screening_csv,
        screening_index=screening_index,
        kinetics=kinetics,
        material_class=material_class,
        genome=genome,
        sweep_kinetics=sweep_kinetics,
        sweep_policy=sweep_policy,
    )


def _load_kinetics_row(job: SweepJob, overrides: Optional[dict] = None,
                       candidate_name: Optional[str] = None):
    """Return (row, kinetics, material_class, candidate_id).

    ``overrides`` are grid-point kinetics values (sweep.kinetics) layered on
    the fixed catalyst.kinetics; ``candidate_name`` names that point's
    mechanism. Screening-row sweeps take neither.
    """
    import pandas as pd
    from pipeline.reactors.mechanisms import (
        CandidateKinetics, parse_catalyst_genome)

    if job.screening_csv is not None:
        path = _resolve(job.screening_csv)
        if not path.is_file():
            raise FileNotFoundError(f'screening CSV not found: {path}')
        frame = pd.read_csv(path)
        if job.screening_index not in frame.index:
            raise KeyError(
                f'screening row {job.screening_index} missing from {path}')
        row = frame.loc[job.screening_index]
        kinetics = CandidateKinetics.from_screening_row(
            row, candidate_id=job.catalyst_name)
        material_class = row.get('material_class')
        material_class = (str(material_class) if isinstance(material_class, str)
                          else None)
        candidate_id = row.get('candidate_id')
        candidate_id = (str(candidate_id) if isinstance(candidate_id, str)
                        else job.catalyst_name)
        return row, kinetics, material_class, candidate_id
    values = {k: v for k, v in job.kinetics.items() if k != 'provenance'}
    overrides = dict(overrides or {})
    values.update(overrides)
    declared = job.kinetics.get('provenance', {})
    sources = {}
    for key, field_name in KINETICS_FIELD_FOR_KEY.items():
        if key in overrides:
            sources[field_name] = 'sweep_yaml_grid: sweep.kinetics'
        elif key in values:
            sources[field_name] = 'sweep_yaml' + (
                f': {declared[key]}' if key in declared else '')
    kwargs = {
        field_name: values.get(key)
        for key, field_name in KINETICS_FIELD_FOR_KEY.items()
        if key != 'E_act'}
    kinetics = CandidateKinetics(
        methane_activation_eV=float(values['E_act']),
        candidate_id=candidate_name or job.catalyst_name,
        screening_protocol='sweep_yaml_literature',
        sources=sources,
        material_class=job.material_class,
        genome=parse_catalyst_genome(job.genome) if job.genome else None,
        **kwargs,
    )
    return values, kinetics, job.material_class, job.catalyst_name


def _reactors_by_mode(reactor_types: List[str]) -> List[tuple]:
    """Group a sweep's reactor list into (pathway_mode, [reactors]) runs.

    Upstream routes reactors through pathway modes and rejects a sweep that
    mixes solids and melt modes; a spec may list both, so we split.
    """
    from pipeline.reactors.models import SINGLE_REACTOR_MODE
    groups: dict = {}
    for reactor in reactor_types:
        groups.setdefault(SINGLE_REACTOR_MODE[reactor], []).append(reactor)
    return list(groups.items())


# Short column headers for swept keys.
_SWEEP_COLUMN = {
    'E_act': 'E_act', 'dE_H': 'dE_H', 'dE_CH3': 'dE_CH3', 'dE_C': 'dE_C',
    'carbon_transfer_eV': 'Eγ', 'carbon_encapsulation_eV': 'Eδ',
    'encapsulation_crossover_coverage': 'θ*',
    'carbon_transfer_prefactor_1_s': 'A_γ',
    'carbon_transfer_particle_nm': 'L_nm',
    'carbon_diffusion_prefactor_m2_s': 'D0',
    'ch4_sticking_coefficient': 's0',
    'max_regen_cycles': 'regen', 'fluidized_mode': 'fl_mode',
    'regen_mechanism': 'regen_by',
}


def _print_table(job: SweepJob, records: list) -> None:
    print(f'\nSweep: {job.name}')
    if job.description:
        print(job.description)
    swept = list(job.sweep_kinetics) + list(job.sweep_policy)
    if swept:
        n_points = len(grid_points(job))
        print(f'grid: {n_points} points over {swept}; {len(records)} records')
    off_site = any(rec.get('off_site_carbon_active') for rec in records)
    extra_hdr = (f" {'X_bound':>8} {'turnov':>7} {'γ/δ':>8} {'gC/gM/h':>8}"
                 f" {'θ_enc':>7} {'life_h':>8}"
                 if off_site else '')
    sweep_hdr = ''.join(f" {_SWEEP_COLUMN.get(k, k)[:9]:>9}" for k in swept)
    print(
        f"{'cell':<22} {'reactor':<10} {'T_K':>7}{sweep_hdr} {'X':>10} "
        f"{'a_1/m':>12} {'WHSV':>8} {'dP_bar':>8}{extra_hdr}  status"
    )

    def fmt(value, spec):
        return format(value, spec) if value is not None else '—'

    def fmt_swept(value):
        if isinstance(value, float):
            return format(value, '.3g')
        return str(value)

    for rec in records:
        status = rec.get('status') or '?'
        if rec.get('reason'):
            status = f"{status} ({rec['reason']})"
        if rec.get('exceeds_equilibrium'):
            status += ' [X > X_eq]'
        extra = ''
        if off_site:
            extra = (
                f" {fmt(rec.get('site_inventory_bound_X'), '.4%'):>8}"
                f" {fmt(rec.get('carbon_turnovers_per_site'), '.2f'):>7}"
                f" {fmt(rec.get('c_gamma_to_c_delta_ratio'), '.3g'):>8}"
                f" {fmt(rec.get('filament_yield_gC_per_gMetal_h'), '.3g'):>8}"
                f" {fmt(rec.get('exit_theta_C_encap'), '.2g'):>7}"
                f" {fmt(rec.get('encapsulation_lifetime_h'), '.3g'):>8}")
        sweep_cols = ''.join(
            f" {fmt_swept(rec.get('sweep', {}).get(k)):>9}" for k in swept)
        print(
            f"{rec['cell']:<22} {rec['reactor_type']:<10} "
            f"{fmt(rec.get('T_K'), '.1f'):>7}{sweep_cols} "
            f"{fmt(rec.get('CH4_conversion'), '.4%'):>10} "
            f"{fmt(rec.get('active_sv_1_m'), '.1f'):>12} "
            f"{fmt(rec.get('WHSV_h-1'), '.1f'):>8} "
            f"{fmt(rec.get('ergun_delta_p_bar'), '.3f'):>8}{extra}  {status}"
        )


def _closure_source(result: dict) -> Optional[str]:
    closure = result.get('reactor_closure_evidence')
    return closure.get('source') if isinstance(closure, dict) else None


def run_sweep(yaml_path: Path) -> dict:
    """Parse YAML, write one mechanism, run every cell × T × reactor, persist JSON.

    Reactors are grouped by pathway mode (solids vs melt) because upstream
    rejects mixed-mode sweeps. Non-applicable cells (e.g. MMBCR on a
    SolidCatalyst) are kept as ``not_applicable`` records, not dropped.

    Args:
        yaml_path: Input controlling yaml path.

    Returns:
        Validated dict output for this operation.
    """
    from pipeline.reactors.mechanisms import write_full_mechanism
    from pipeline.reactors.models import run_reactor_sweep

    job = parse_sweep(yaml_path)
    out_dir = SWEEPS_DIR / job.name
    points = grid_points(job)
    swept = bool(job.sweep_kinetics or job.sweep_policy)

    # One mechanism per distinct kinetics combination (policy points share it).
    mech_for_kinetics: dict = {}
    mechanism_files: List[str] = []
    records = []
    material_class = candidate_id = None
    for point in points:
        kin_key = tuple(sorted(point.kinetics.items()))
        if kin_key not in mech_for_kinetics:
            name = (f'{job.catalyst_name}_p{len(mech_for_kinetics):03d}'
                    if job.sweep_kinetics else job.catalyst_name)
            row, kinetics, material_class, candidate_id = _load_kinetics_row(
                job, overrides=point.kinetics, candidate_name=name)
            try:
                e_act = float(row.get('E_act', kinetics.methane_activation_eV))
                dE_H = float(row.get('dE_H', 0.0) or 0.0)
            except (TypeError, ValueError, AttributeError):
                e_act = float(kinetics.methane_activation_eV)
                dE_H = float(job.kinetics.get('dE_H') or 0.0)
            mech = write_full_mechanism(
                name, kinetics=kinetics, output_dir=out_dir / 'mechanisms')
            mech_for_kinetics[kin_key] = (name, mech, e_act, dE_H)
            mechanism_files.append(repo_relative(mech))
        # The surface phase is '<name>_surface'; the reactor must load by
        # the same name the writer used for this grid point.
        mech_name, mech, e_act, dE_H = mech_for_kinetics[kin_key]
        policy = {**job.policy, **point.policy}
        for cell in job.cells:
            for pathway_mode, reactors in _reactors_by_mode(job.reactor_types):
                results = run_reactor_sweep(
                    mech_name, str(mech),
                    temperatures=list(job.temperatures_K),
                    reactor_types=reactors,
                    catalyst_E_act_eV=e_act,
                    pathway_mode=pathway_mode,
                    material_class=material_class,
                    candidate_id=candidate_id,
                    catalyst_dE_H_eV=dE_H,
                    reactor_config_kwargs={
                        **policy,
                        'catalyst_particle_mm': cell.catalyst_particle_mm,
                        'metal_loading': cell.metal_loading,
                        'metal_dispersion': cell.metal_dispersion,
                    },
                )
                for result in results:
                    records.append(_record(
                        cell, point, policy, pathway_mode, material_class,
                        repo_relative(mech), result))

    job_record = asdict(job)
    job_record['source_file'] = repo_relative(job.source_file)
    payload = {
        'job': job_record,
        'material_class': material_class,
        'candidate_id': candidate_id,
        'mechanism_file': mechanism_files[0] if mechanism_files else None,
        'mechanism_files': mechanism_files,
        'grid_points': [
            {'index': p.index, 'kinetics': p.kinetics, 'policy': p.policy}
            for p in points] if swept else [],
        'written_at': datetime.now(timezone.utc).isoformat(),
        'records': records,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / 'run.json'
    out_json.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    shutil.copy2(job.source_file, out_dir / 'input.yaml')
    logger.info(f'Wrote {repo_relative(out_json)} ({len(records)} records)')
    _print_table(job, records)
    print(f'\nWrote {repo_relative(out_json)}')
    return payload


def _record(cell: SweepCell, point: GridPoint, policy: dict, pathway_mode: str,
            material_class: Optional[str], mechanism_file: str,
            result: dict) -> dict:
    return {
        'cell': cell.name,
        'catalyst_particle_mm': cell.catalyst_particle_mm,
        'metal_loading': cell.metal_loading,
        'metal_dispersion': cell.metal_dispersion,
        'reactor_type': result.get('reactor_type'),
        'pathway_mode': pathway_mode,
        'material_class': material_class,
        'point': point.index,
        'sweep': point.values,
        'max_regen_cycles': policy.get('max_regen_cycles'),
        'fluidized_mode': policy.get('fluidized_mode'),
        'regen_mechanism': policy.get('regen_mechanism'),
        'mechanism_file': mechanism_file,
        'T_K': result.get('T_K'),
        'status': result.get('status'),
        'reason': result.get('reason') or result.get('error'),
        'CH4_conversion': result.get('CH4_conversion'),
        'emulsion_CH4_conversion': result.get('emulsion_CH4_conversion'),
        'active_sv_1_m': result.get('active_sv_1_m'),
        'WHSV_h-1': result.get('WHSV_h-1'),
        'residence_time_s': result.get('residence_time_s'),
        'ergun_delta_p_bar': result.get('ergun_delta_p_bar'),
        'closure_source': _closure_source(result),
        'can_exclude_candidate': result.get('can_exclude_candidate'),
        'surface_loaded': result.get('surface_loaded'),
        'graphite_loaded': result.get('graphite_loaded'),
        # B5/B6 closure accounting (PFR and Fluidized emulsion); None elsewhere.
        'site_inventory_bound_X': result.get('site_inventory_bound_X'),
        'site_inventory_bound_basis': result.get('site_inventory_bound_basis'),
        'X_eq_table': result.get('X_eq_table'),
        'exceeds_equilibrium': result.get('exceeds_equilibrium'),
        'carbon_balance_ok': result.get('carbon_balance_ok'),
        'carbon_balance_residual_mol': result.get('carbon_balance_residual_mol'),
        'mock': result.get('mock'),
        'carbon_turnovers_per_site': result.get('carbon_turnovers_per_site'),
        'turnover_factor_vs_bound': result.get('turnover_factor_vs_bound'),
        'off_site_carbon_active': result.get('off_site_carbon_active'),
        'c_gamma_to_c_delta_ratio': result.get('c_gamma_to_c_delta_ratio'),
        'exit_theta_C': result.get('exit_theta_C'),
        'exit_theta_C_encap': result.get('exit_theta_C_encap'),
        'filament_yield_gC_per_gMetal_h': result.get(
            'filament_yield_gC_per_gMetal_h'),
        'filament_yield_within_band': result.get('filament_yield_within_band'),
        # B6-6 encapsulation onset and lifetime.
        'encapsulation_onset': result.get('encapsulation_onset'),
        'c_delta_competitive': result.get('c_delta_competitive'),
        'encapsulation_lifetime_h': result.get('encapsulation_lifetime_h'),
        'lifetime_within_tos_band': result.get('lifetime_within_tos_band'),
        'outfeed_carbon_mol_per_pass': result.get('outfeed_carbon_mol_per_pass'),
        'regen_cycles_completed': result.get('regen_cycles_completed'),
    }
