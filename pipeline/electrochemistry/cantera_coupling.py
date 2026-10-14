"""Cantera interface step of the runner-owned MOOSE <-> Cantera coupling loop.

One invocation follows each MOOSE run. It reads the MOOSE interface state,
evaluates the candidate's Cantera interface mechanism at that state, and
writes the exchanged parameter MOOSE consumes in the next iteration. It also
packages the outputs, conservation budgets, mesh-refinement series and the
hashed coupling proof that ``external_runner`` validates.

Coupled quantities
------------------
Electrochemical: the anode exchange current density
    i0_a = n F r_CH4(T, x_CH4,surface),
where r_CH4 is Cantera's CH4 destruction rate on the candidate interface at
the MOOSE anode-surface CH4 concentration (via Henry's law). The electrode
potential enters MOOSE's Butler-Volmer law, not the Cantera mechanism.

NTEC: two exchanged quantities.
    sigma = sigma_sat (1 - exp(-gamma_wall / gamma_c))
is the contact-electrification charge driven by MOOSE's wall shear rate, and
    k_s = r_CH4(T, x_CH4,surface) exp(min(beta |E|, cap) / k_B T) / c_surface
is the first-order surface rate constant MOOSE applies at the interface,
where |E| is MOOSE's interfacial field from sigma. sigma_sat, gamma_c,
beta [eV m/V] and cap [eV] are case inputs that require paired-control
calibration. Carbon detachment uses the impact/adhesion estimate in
``ntec_mechanics``.

Both choices are documented modeling assumptions; neither is validated until
the case's held-out calibration passes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.electrochemistry.ntec_mechanics import (  # noqa: E402
    detachment_probability, impact_energy_J)

FARADAY = 96485.33212
K_B_EV = 8.617333262145e-5
M_CH4 = 0.01604246
M_C = 0.0120107
M_H = 0.00100794
MOOSE_CSV = 'hydrogen_moose.csv'
PARAMETER_RECORD = 'hydrogen_moose_parameters.json'
FEEDBACK = 'hydrogen_interface_rates.json'
HISTORY = 'hydrogen_coupling_history.json'
RATE_EXCHANGE = 'hydrogen_rate_exchange.json'
CANTERA_LOG = 'hydrogen_cantera.log'
REFINEMENT_LEVELS = (1, 2)
# Exchanged MOOSE variable -> case field holding its iteration-1 guess.
EXCHANGED = {
    'Electrochemical': {'anode_exchange_current_A_m2': 'anode_exchange_current_guess_A_m2'},
    'NTEC': {'interface_charge_C_m2': 'interface_charge_guess_C_m2',
             'interface_rate_constant_m_s': 'interface_rate_constant_guess_m_s'},
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=1, sort_keys=True) + '\n')


def _last_row(path: Path) -> dict[str, float]:
    if not path.is_file():
        raise RuntimeError(f'MOOSE postprocessor output missing: {path.name}')
    rows = list(csv.DictReader(path.open()))
    if not rows:
        raise RuntimeError(f'MOOSE postprocessor output is empty: {path.name}')
    values = {key: float(value) for key, value in rows[-1].items()}
    if not all(math.isfinite(value) for value in values.values()):
        raise RuntimeError('MOOSE postprocessors contain non-finite values')
    return values


def _parameter(case: dict, name: str, default: float | None = None) -> float:
    for section in ('geometry', 'operating', 'properties', 'model_parameters'):
        value = case.get(section, {}).get(name)
        if value is not None:
            return float(value)
    if default is None:
        raise RuntimeError(f'case lacks required parameter {name!r}')
    return float(default)


def _fraction(name: str, value: float) -> float:
    """Return a [0, 1] fraction, absorbing only floating-point round-off."""
    if not -1e-9 <= value <= 1.0 + 1e-9:
        raise RuntimeError(f'{name}={value} is outside [0, 1]')
    return min(max(value, 0.0), 1.0)


def _used_values(record: dict, case: dict, reactor_type: str) -> dict[str, float]:
    overrides = dict(token.partition('=')[::2] for token in record.get('overrides', []))
    return {name: float(overrides[name]) if name in overrides else _parameter(case, guess)
            for name, guess in EXCHANGED[reactor_type].items()}


def interface_ch4_destruction_mol_m2_s(case_dir: Path, case: dict,
                                       x_ch4: float) -> tuple[float, str]:
    """Evaluate the candidate interface mechanism at the MOOSE surface state.

    Args:
        case_dir: Case directory containing the declared mechanism.
        case: Physical case document.
        x_ch4: CH4 mole fraction in the gas in equilibrium with the surface.

    Returns:
        CH4 destruction rate in mol m^-2 s^-1 and a log line.
    """
    import cantera as ct

    kinetics = case.get('kinetics', {})
    mechanism = case_dir / str(kinetics.get('mechanism', ''))
    if not mechanism.is_file() or not mechanism.resolve().is_relative_to(case_dir.resolve()):
        raise RuntimeError('kinetics.mechanism must name a Cantera YAML file in the case')
    gas = ct.Solution(str(mechanism), str(kinetics['gas_phase']))
    surface = ct.Interface(str(mechanism), str(kinetics['interface_phase']), [gas])
    balance = str(kinetics.get('balance_species', 'Ar'))
    x_ch4 = min(max(x_ch4, 0.0), 1.0)
    temperature = _parameter(case, 'temperature_K')
    pressure = _parameter(case, 'pressure_Pa')
    gas.TPX = temperature, pressure, {'CH4': max(x_ch4, 1e-30), balance: max(1.0 - x_ch4, 1e-30)}
    surface.TP = temperature, pressure
    surface.advance_coverages(float(kinetics.get('coverage_relaxation_s', 1.0e2)))
    index = surface.kinetics_species_index('CH4')
    rate = float(surface.destruction_rates[index]) * 1.0e3
    if not math.isfinite(rate) or rate < 0:
        raise RuntimeError('Cantera returned an invalid CH4 destruction rate')
    line = (f'cantera {ct.__version__} mechanism={mechanism.name} T={temperature} '
            f'P={pressure} x_CH4={x_ch4:.6e} r_CH4={rate:.6e} mol/m2/s')
    return rate, line


def _electrochemical(case: dict, row: dict, rate: float) -> tuple[dict, dict, dict]:
    n = _parameter(case, 'electrons_per_ch4', 4)
    share = _parameter(case, 'ch4_current_fraction', 1.0)
    fe_h2 = _parameter(case, 'cathode_h2_faradaic_efficiency', 1.0)
    area = _parameter(case, 'electrode_area_m2')
    voltage = _parameter(case, 'applied_potential_V')
    feed = _parameter(case, 'methane_mass_flow_kg_s') / M_CH4
    i_anode = row['anode_current_density_A_m2']
    i_cathode = row['cathode_current_density_A_m2']
    i_ionic = row['ionic_current_density_A_m2']
    consumed = row['ch4_consumption_mol_m2_s'] * area
    supplied = -row['ch4_supply_mol_m2_s']
    if i_anode < 0:
        raise RuntimeError('anode current is cathodic at the applied potential')
    conversion = consumed / feed
    if conversion > 1.0:
        raise RuntimeError('CH4 consumption exceeds the declared feed')
    h2 = fe_h2 * i_cathode * area / (2.0 * FARADAY)
    outputs = {
        'current_density_A_cm2': i_anode / 1.0e4,
        'cell_voltage_V': voltage,
        'electrical_power_density_W_cm2': voltage * i_anode / 1.0e4,
        'faradaic_efficiency_H2': fe_h2,
        'CH4_conversion': _fraction('CH4_conversion', conversion),
        'H2_selectivity': _fraction('H2_selectivity', h2 / (2.0 * consumed)) if consumed > 0 else 0.0,
        'anode_overpotential_V': row['anode_overpotential_V'],
        'cathode_overpotential_V': row['cathode_overpotential_V'],
        'ohmic_drop_V': row['ohmic_drop_V'],
        'anode_surface_ch4_mol_m3': row['anode_surface_ch4_mol_m3'],
    }
    h_released = share * i_ionic / FARADAY
    budgets = {
        'charge': {'inlet': i_anode, 'outlet': i_cathode},
        'carbon': {'inlet': supplied, 'outlet': row['ch4_consumption_mol_m2_s']},
        'hydrogen': {'inlet': 4.0 * supplied, 'outlet': h_released * 4.0 / n},
        'mass': {'inlet': supplied * M_CH4,
                 'outlet': row['ch4_consumption_mol_m2_s'] * M_C
                 + h_released * 4.0 / n * M_H},
    }
    return {'anode_exchange_current_A_m2': n * FARADAY * rate}, outputs, budgets


def _ntec(case: dict, row: dict, rate: float) -> tuple[dict, dict, dict]:
    temperature = _parameter(case, 'temperature_K')
    shear = abs(row['wall_shear_rate_s_inv'])
    field = abs(row['interface_field_V_m'])
    beta = _parameter(case, 'field_barrier_coefficient_eV_m_V')
    cap = _parameter(case, 'field_barrier_cap_eV', 0.25)
    enhancement = math.exp(min(beta * field, cap) / (K_B_EV * temperature))
    surface_c = row['interface_ch4_mol_m3']
    if surface_c <= 0:
        raise RuntimeError('interface CH4 concentration is not positive')
    sigma = _parameter(case, 'contact_charge_saturation_C_m2') * (
        1.0 - math.exp(-shear / _parameter(case, 'contact_shear_scale_s_inv')))
    area = _parameter(case, 'interface_area_m2')
    feed = _parameter(case, 'methane_mass_flow_kg_s') / M_CH4
    consumed = row['ch4_consumption_mol_m2_s'] * area
    conversion = consumed / feed
    if conversion > 1.0:
        raise RuntimeError('CH4 consumption exceeds the declared feed')
    detached = detachment_probability(
        _parameter(case, 'adhesion_energy_J_m2'),
        impact_energy_J(_parameter(case, 'effective_mass_kg'), shear,
                        _parameter(case, 'particle_diameter_m')),
        _parameter(case, 'contact_area_m2'))
    h2_kg_s = 2.0 * consumed * 2.0 * M_H
    mechanical_w = row['wall_work_W_m2'] * area
    outputs = {
        'CH4_conversion': _fraction('CH4_conversion', conversion),
        # The candidate interface mechanism has no gas-phase carbon product.
        'H2_selectivity': 1.0,
        'solid_C_selectivity': 1.0,
        'specific_energy_kWh_kg_H2': mechanical_w / h2_kg_s / 3.6e6 if h2_kg_s > 0 else math.inf,
        'carbon_detachment_probability': detached,
        'interface_field_V_m': field,
        'interface_charge_C_m2': row['interface_charge_C_m2'],
        'field_rate_enhancement': enhancement,
        'wall_shear_rate_s_inv': shear,
    }
    supplied = -row['ch4_supply_mol_m2_s']
    budgets = {
        'carbon': {'inlet': supplied, 'outlet': row['ch4_consumption_mol_m2_s']},
        'hydrogen': {'inlet': 4.0 * supplied, 'outlet': 4.0 * row['ch4_consumption_mol_m2_s']},
        'mass': {'inlet': supplied * M_CH4, 'outlet': row['ch4_consumption_mol_m2_s'] * M_CH4},
        'charge': {'inlet': row['interface_charge_C_m2'], 'outlet': row['gauss_charge_C_m2']},
        'energy': {'inlet': row['wall_work_W_m2'], 'outlet': row['viscous_dissipation_W_m2']},
    }
    exchanged = {'interface_charge_C_m2': sigma,
                 'interface_rate_constant_m_s': rate * enhancement / surface_c}
    return exchanged, outputs, budgets


def _mesh_series(case_dir: Path, record: dict, base_row: dict, observable,
                 timeout_s: int) -> list[dict]:
    series = [{'cells': int(base_row['elements']), 'observable': observable(base_row)}]
    command = list(record['command'])
    for level in REFINEMENT_LEVELS:
        stem = f'hydrogen_moose_refine{level}'
        refined = [token for token in command if not token.startswith('Outputs/file_base=')]
        refined += [f'Outputs/file_base={case_dir / stem}', f'Mesh/uniform_refine={level}',
                    'Outputs/exodus=false']
        completed = subprocess.run(refined, cwd=case_dir, capture_output=True,
                                   text=True, timeout=timeout_s)
        (case_dir / f'{stem}.log').write_text(completed.stdout + completed.stderr)
        if completed.returncode:
            raise RuntimeError(f'MOOSE refinement level {level} failed')
        row = _last_row(case_dir / f'{stem}.csv')
        series.append({'cells': int(row['elements']), 'observable': observable(row)})
    return series


def run_coupling_step(case_dir: Path, reactor_type: str, timeout_s: int = 3600) -> dict:
    """Execute one Cantera interface step and write all exchange artifacts.

    Args:
        case_dir: Case directory after the current MOOSE iteration.
        reactor_type: ``Electrochemical`` or ``NTEC``.
        timeout_s: Wall-time bound for each mesh-refinement MOOSE run.

    Returns:
        The coupling state written for the runner.
    """
    if reactor_type not in EXCHANGED:
        raise ValueError(f'{reactor_type} has no Cantera interface coupling')
    case = json.loads((case_dir / 'hydrogen_case.json').read_text())
    request = json.loads((case_dir / 'hydrogen_coupling_request.json').read_text())
    record = json.loads((case_dir / PARAMETER_RECORD).read_text())
    iteration = int(request['iteration'])
    row = _last_row(case_dir / MOOSE_CSV)
    henry = _parameter(case, 'methane_henry_mol_m3_Pa')
    pressure = _parameter(case, 'pressure_Pa')
    surface_c = row['anode_surface_ch4_mol_m3' if reactor_type == 'Electrochemical'
                    else 'interface_ch4_mol_m3']
    rate, log_line = interface_ch4_destruction_mol_m2_s(
        case_dir, case, surface_c / (henry * pressure))
    model = _electrochemical if reactor_type == 'Electrochemical' else _ntec
    computed, outputs, budgets = model(case, row, rate)
    used = _used_values(record, case, reactor_type)
    residual = max(abs(computed[name] - used[name]) / max(abs(computed[name]), 1e-300)
                   for name in computed)
    tolerance = _parameter(case, 'coupling_tolerance_relative', 1e-6)
    relaxation = _parameter(case, 'coupling_relaxation', 1.0)
    if not 0 < relaxation <= 1:
        raise RuntimeError('coupling_relaxation must lie in (0, 1]')
    converged = residual <= tolerance
    next_values = {name: used[name] + relaxation * (computed[name] - used[name])
                   for name in computed}
    feedback = case_dir / FEEDBACK
    _write_json(feedback, {'schema_version': 1, 'iteration': iteration,
                           'computed': computed, 'moose_values_used': used,
                           'moose_parameters': next_values})
    log = case_dir / CANTERA_LOG
    with log.open('a') as handle:
        handle.write(f'iteration {iteration}: {log_line}\n')
    history_path = case_dir / HISTORY
    history = json.loads(history_path.read_text()) if history_path.is_file() else {'iterations': []}
    history['iterations'] = [row_ for row_ in history['iterations']
                             if int(row_['iteration']) < iteration]
    history['iterations'].append({'iteration': iteration, 'residual_relative': residual,
                                  'exchanged': computed})
    _write_json(history_path, history)
    temperature = _parameter(case, 'temperature_K')
    thickness = _parameter(case, 'electrolyte_thickness_m') if reactor_type == 'Electrochemical' \
        else _parameter(case, 'film_thickness_m')
    _write_json(case_dir / RATE_EXCHANGE, {
        'schema_version': 1, 'candidate_id': case['candidate_id'],
        'temperature_K': temperature,
        'reaction_rates_mol_m3_s': {'CH4': -rate / thickness},
        'interface_rate_mol_m2_s': {'CH4': -rate}})
    _write_json(case_dir / 'hydrogen_outputs.json', outputs)
    convergence = {
        'converged': converged,
        'mesh_series': [],
        'mesh_tolerance_relative': _parameter(case, 'mesh_tolerance_relative', 0.01),
        'conservation_budgets': budgets,
    }
    if converged and iteration >= 2:
        index = {'Electrochemical': 'anode_current_density_A_m2',
                 'NTEC': 'ch4_consumption_mol_m2_s'}[reactor_type]
        convergence['mesh_series'] = _mesh_series(
            case_dir, record, row, lambda values: values[index], timeout_s)
    _write_json(case_dir / 'hydrogen_convergence.json', convergence)
    if reactor_type == 'NTEC':
        _write_ntec_handoff(case_dir, case, row, iteration)
    _write_metadata(case_dir, case, reactor_type, iteration, residual, tolerance)
    state = {
        'schema_version': 1, 'candidate_id': case['candidate_id'],
        'reactor_type': reactor_type, 'temperature_K': temperature,
        'iteration': iteration, 'residual_relative': residual,
        'tolerance_relative': tolerance, 'converged': converged,
        'feedback_artifact': {'path': feedback.name, 'sha256': _sha256(feedback)},
    }
    _write_json(case_dir / 'hydrogen_coupling_state.json', state)
    return state


def _write_ntec_handoff(case_dir: Path, case: dict, row: dict, iteration: int) -> None:
    field = case_dir / 'hydrogen_moose.e'
    if not field.is_file():
        raise RuntimeError('NTEC MOOSE run did not write its Exodus field file')
    _write_json(case_dir / 'hydrogen_hydrodynamics.json', {
        'schema_version': 1, 'candidate_id': case['candidate_id'],
        'pathway_mode': case['pathway_mode'], 'reactor_type': 'NTEC',
        'temperature_K': _parameter(case, 'temperature_K'), 'iteration': iteration,
        'fields': {
            'velocity_m_s': row['mean_velocity_m_s'],
            'pressure_Pa': _parameter(case, 'pressure_Pa'),
            'temperature_K': _parameter(case, 'temperature_K'),
            'liquid_volume_fraction': 1.0,
            'shear_rate_s_inv': row['wall_shear_rate_s_inv'],
        },
        'field_artifact': {'path': field.name, 'sha256': _sha256(field),
                           'format': 'EXODUS', 'mesh_id': f'moose-{int(row["elements"])}',
                           'coordinate_system': 'cartesian-m'},
    })


def _write_metadata(case_dir: Path, case: dict, reactor_type: str, iteration: int,
                    residual: float, tolerance: float) -> None:
    kinetics = case.get('kinetics', {})
    mechanism = case_dir / str(kinetics['mechanism'])
    fields = ['species', 'temperature', 'reaction_heat', 'reaction_rates',
              'charge', 'potential']
    if reactor_type == 'NTEC':
        fields.append('momentum')
    proof = {
        'schema_version': 1, 'cantera_used': True,
        'coupling_method': 'iterative_two_way',
        'coupling_iterations': iteration,
        'coupling_residual_relative': residual,
        'coupling_tolerance_relative': tolerance,
        'exchanged_fields': fields,
    }
    for stem, path in (('mechanism', mechanism), ('cantera_log', case_dir / CANTERA_LOG),
                       ('rate_exchange', case_dir / RATE_EXCHANGE),
                       ('coupling_history', case_dir / HISTORY)):
        proof[f'{stem}_path'] = path.name
        proof[f'{stem}_sha256'] = _sha256(path)
    metadata = {'solver_coupling': proof,
                'mechanism': {'complete': True, 'source': kinetics.get('source')}}
    if reactor_type == 'Electrochemical':
        metadata['electrolyte_phase'] = case['electrolyte_phase']
    else:
        calibration = case.get('calibration', {})
        metadata['calibration'] = {
            'paired_control': calibration.get('paired_control') is True,
            'paired_control_source': calibration.get('paired_control_source'),
        }
    _write_json(case_dir / 'hydrogen_metadata.json', metadata)


def main() -> None:
    """Run one coupling step from the command line."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--case-dir', required=True)
    parser.add_argument('--reactor-type', required=True, choices=sorted(EXCHANGED))
    parser.add_argument('--timeout-s', type=int, default=3600)
    args = parser.parse_args()
    state = run_coupling_step(Path(args.case_dir).resolve(), args.reactor_type, args.timeout_s)
    print(json.dumps(state, sort_keys=True))


if __name__ == '__main__':
    main()
