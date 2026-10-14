"""Analytic verification of the MOOSE NTEC and electrochemical models.

Opt-in suite (``python tests/run_tests.py --include-moose``): it runs the real
``moose-opt`` on the repository input templates and compares the solutions
with closed-form results. A missing MOOSE installation is a failure, not a
skip, because this suite exists to verify that installation.
"""

from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from pipeline.electrochemistry.cantera_coupling import _last_row  # noqa: E402
from pipeline.electrochemistry.kinetics import butler_volmer_current_density  # noqa: E402
from pipeline.electrochemistry.transport import (  # noqa: E402
    limiting_current_A, ohmic_resistance_ohm)
from pipeline.simulation.moose_backend import (  # noqa: E402
    build_moose_command, default_moose_input, moose_environment,
    resolve_moose_executable)

FARADAY = 96485.33212
GAS_CONSTANT = 8.314462618


def _solve(reactor_type: str, parameters: dict, refine: int = 0) -> dict[str, float]:
    executable = resolve_moose_executable(required=True)
    with tempfile.TemporaryDirectory() as tmp:
        case = Path(tmp)
        command = build_moose_command(
            executable, default_moose_input(reactor_type), case, parameters)
        if refine:
            command.append(f'Mesh/uniform_refine={refine}')
        completed = subprocess.run(command, cwd=case, capture_output=True,
                                   text=True, timeout=600,
                                   env=moose_environment(executable))
        assert completed.returncode == 0, completed.stdout[-3000:] + completed.stderr[-3000:]
        return _last_row(case / 'hydrogen_moose.csv')


def _close(actual: float, expected: float, rel: float, label: str) -> None:
    assert math.isclose(actual, expected, rel_tol=rel, abs_tol=0.0), (
        f'{label}: MOOSE {actual!r} vs analytic {expected!r}')


ELECTROCHEMICAL = {
    'electrolyte_thickness_m': 1e-3, 'temperature_K': 353.15, 'pressure_Pa': 101325.0,
    'ionic_conductivity_S_m': 10.0, 'methane_diffusivity_m2_s': 1e-9,
    'methane_henry_mol_m3_Pa': 1.4e-5, 'anode_equilibrium_potential_V': 0.2,
    'cathode_equilibrium_potential_V': 0.0, 'anode_transfer_coefficient': 0.5,
    'cathode_transfer_coefficient': 0.5, 'cathode_exchange_current_A_m2': 10.0,
    'ch4_current_fraction': 1.0, 'electrons_per_ch4': 4, 'mesh_elements': 100,
}


def test_electrochemical_cell_reaches_the_analytic_limiting_current():
    parameters = dict(ELECTROCHEMICAL, applied_potential_V=1.4,
                      anode_exchange_current_guess_A_m2=1e-2)
    row = _solve('Electrochemical', parameters)
    c_bulk = parameters['methane_henry_mol_m3_Pa'] * parameters['pressure_Pa']
    limit = limiting_current_A(4, parameters['methane_diffusivity_m2_s'], c_bulk, 1.0,
                               parameters['electrolyte_thickness_m'])
    _close(row['anode_current_density_A_m2'], limit, 1e-4, 'limiting current')
    _close(row['cathode_current_density_A_m2'], row['anode_current_density_A_m2'],
           1e-9, 'charge balance')
    ohmic = row['anode_current_density_A_m2'] * ohmic_resistance_ohm(
        parameters['electrolyte_thickness_m'], parameters['ionic_conductivity_S_m'], 1.0)
    _close(row['ohmic_drop_V'], ohmic, 1e-9, 'ohmic drop')
    _close(-row['ch4_supply_mol_m2_s'], row['ch4_consumption_mol_m2_s'], 1e-9,
           'CH4 supply equals consumption')


def test_electrochemical_kinetic_regime_matches_butler_volmer():
    parameters = dict(ELECTROCHEMICAL, applied_potential_V=0.5,
                      anode_exchange_current_guess_A_m2=1e-6)
    row = _solve('Electrochemical', parameters)
    c_bulk = parameters['methane_henry_mol_m3_Pa'] * parameters['pressure_Pa']
    surface_ratio = row['anode_surface_ch4_mol_m3'] / c_bulk
    assert surface_ratio > 0.999, 'case must be kinetically limited'
    eta = row['anode_overpotential_V']
    frt = FARADAY / (GAS_CONSTANT * parameters['temperature_K'])
    expected = 1e-6 * (surface_ratio * math.exp(0.5 * frt * eta) - math.exp(-0.5 * frt * eta))
    _close(row['anode_current_density_A_m2'], expected, 1e-8, 'Butler-Volmer current')
    reference = butler_volmer_current_density(1e-6, eta, parameters['temperature_K'], 0.5, 1)
    _close(row['anode_current_density_A_m2'], reference * surface_ratio, 1e-3,
           'agreement with pipeline.electrochemistry.kinetics')


def test_electrochemical_solution_is_mesh_independent():
    parameters = dict(ELECTROCHEMICAL, applied_potential_V=0.9,
                      anode_exchange_current_guess_A_m2=1e-4)
    coarse = _solve('Electrochemical', parameters)
    fine = _solve('Electrochemical', parameters, refine=2)
    _close(fine['anode_current_density_A_m2'], coarse['anode_current_density_A_m2'], 1e-6,
           'mesh refinement')


NTEC = {
    'temperature_K': 373.15, 'pressure_Pa': 101325.0, 'shear_rate_s_inv': 150.0,
    'liquid_viscosity_Pa_s': 2e-3, 'permittivity_F_m': 7e-11,
    'liquid_density_kg_m3': 6100.0, 'film_thickness_m': 1e-3,
    'dielectric_thickness_m': 1e-5, 'methane_diffusivity_m2_s': 2e-9,
    'methane_henry_mol_m3_Pa': 1e-5, 'interface_charge_guess_C_m2': 5e-5,
    'interface_rate_constant_guess_m_s': 1e-6, 'mesh_elements_x': 4,
    'mesh_elements_film': 20, 'mesh_elements_dielectric': 5,
}


def test_ntec_layer_matches_couette_gauss_and_robin_solutions():
    row = _solve('NTEC', NTEC)
    h = NTEC['film_thickness_m']
    shear = NTEC['shear_rate_s_inv']
    _close(row['wall_shear_rate_s_inv'], shear, 1e-8, 'Couette wall shear rate')
    _close(row['mean_velocity_m_s'], shear * h / 2.0, 1e-8, 'Couette mean velocity')
    _close(abs(row['interface_field_V_m']),
           NTEC['interface_charge_guess_C_m2'] / NTEC['permittivity_F_m'], 1e-8,
           'Gauss interfacial field')
    _close(row['viscous_dissipation_W_m2'], NTEC['liquid_viscosity_Pa_s'] * shear**2 * h,
           1e-8, 'viscous dissipation')
    _close(row['wall_work_W_m2'], row['viscous_dissipation_W_m2'], 1e-8, 'energy balance')
    c_bulk = NTEC['methane_henry_mol_m3_Pa'] * NTEC['pressure_Pa']
    k = NTEC['interface_rate_constant_guess_m_s']
    expected_surface = c_bulk / (1.0 + k * h / NTEC['methane_diffusivity_m2_s'])
    _close(row['interface_ch4_mol_m3'], expected_surface, 1e-8, 'Robin surface CH4')
    _close(row['ch4_consumption_mol_m2_s'], k * expected_surface, 1e-8, 'surface consumption')
    _close(-row['ch4_supply_mol_m2_s'], row['ch4_consumption_mol_m2_s'], 1e-8,
           'CH4 supply equals consumption')


def test_ntec_solution_is_mesh_independent():
    coarse = _solve('NTEC', NTEC)
    fine = _solve('NTEC', NTEC, refine=1)
    for name in ('ch4_consumption_mol_m2_s', 'wall_shear_rate_s_inv', 'interface_field_V_m'):
        _close(fine[name], coarse[name], 1e-6, f'mesh refinement: {name}')


def main() -> None:
    tests = [value for name, value in sorted(globals().items())
             if name.startswith('test_') and callable(value)]
    for test in tests:
        test()
        print('PASS', test.__name__)
    print(f'{len(tests)}/{len(tests)} MOOSE model verifications passed')


if __name__ == '__main__':
    main()
