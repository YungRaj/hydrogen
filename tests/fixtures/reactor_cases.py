"""Representative reactor fixtures that exercise real pipeline contracts.

External OpenFOAM and FEniCSx processes are replaced by deterministic writers.
Everything around that boundary remains real. These are test fixtures, never
scientific or production evidence.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from pipeline.simulation.solver_execution import SolverExecutionServices

MODES = {"Fluidized": "thermocatalytic_fluidized", "MMBCR": "mmbcr",
         "NTEC": "ntec", "Electrochemical": "electrochemical"}


def _write_json(path: Path, value: dict) -> None:
    """Serialize a fixture object deterministically."""
    path.write_text(json.dumps(value, sort_keys=True) + "\n")


def _sha256(path: Path) -> str:
    """Return the digest recorded by artifact contracts."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_case(case: Path, reactor_type: str, candidate_id: str,
               temperature_K: float) -> None:
    """Write a complete, sourced, non-template physical reactor case."""
    mode = MODES[reactor_type]
    sections = {
        "Fluidized": {"geometry": {"column_diameter_m": .1, "bed_height_m": .8},
            "operating": {"temperature_K": temperature_K, "pressure_Pa": 101325,
                          "methane_mass_flow_kg_s": 1e-4},
            "properties": {"particle_diameter_m": 1e-3,
                           "particle_density_kg_m3": 2500,
                           "gas_viscosity_Pa_s": 2e-5},
            "models": {"drag_model": "Gidaspow",
                       "heat_transfer_model": "Ranz-Marshall"}},
        "MMBCR": {"geometry": {"column_diameter_m": .1, "liquid_height_m": 1.,
                                "sparger_orifice_diameter_m": 1e-3},
            "operating": {"temperature_K": temperature_K, "pressure_Pa": 101325,
                          "methane_mass_flow_kg_s": 1e-4},
            "properties": {"liquid_density_kg_m3": 6000,
                           "liquid_viscosity_Pa_s": 2e-3,
                           "surface_tension_N_m": .7},
            "models": {"bubble_breakup_model": "Lehr",
                       "bubble_coalescence_model": "Prince-Blanch"}},
        "NTEC": {"geometry": {"reactor_volume_m3": 1e-3,
                                "interface_area_m2": .1},
            "operating": {"temperature_K": temperature_K, "pressure_Pa": 101325,
                          "methane_mass_flow_kg_s": 1e-6,
                          "shear_rate_s_inv": 100., "mechanical_power_W_kg": 20.},
            "properties": {"liquid_viscosity_Pa_s": 1e-3,
                           "permittivity_F_m": 7e-10,
                           "ionic_conductivity_S_m": 1.},
            "models": {"contact_electrification_model": "measured-source",
                       "species_transport_model": "Nernst-Planck"}},
        "Electrochemical": {"geometry": {"electrode_area_m2": .01,
                                           "electrolyte_thickness_m": 1e-3},
            "operating": {"temperature_K": temperature_K, "pressure_Pa": 101325,
                          "methane_mass_flow_kg_s": 1e-6,
                          "applied_potential_V": 1.4},
            "properties": {"ionic_conductivity_S_m": 1.,
                           "electronic_conductivity_S_m": 100.,
                           "methane_diffusivity_m2_s": 1e-9},
            "models": {"charge_transfer_model": "Butler-Volmer",
                       "species_transport_model": "Nernst-Planck"}},
    }[reactor_type]
    value = {"schema_version": 1, "candidate_id": candidate_id,
        "pathway_mode": mode, "reactor_type": reactor_type, **sections,
        "feed": {"composition": {"CH4": 1.0}, "source": "fixture:feed-v1"},
        "kinetics": {"source": "fixture:kinetics-v1"},
        "calibration": {"training_ids": ["fixture-train-1", "fixture-train-2"],
            "validation_ids": ["fixture-holdout-1"],
            "source": "fixture:synthetic-validation-v1", "metric": "relative_rmse",
            "acceptance_threshold": .1}}
    value["parameter_sources"] = {
        f"{section}.{name}": "fixture:representative-input-v1"
        for section in ("geometry", "operating", "properties")
        for name in value[section]}
    if reactor_type == "NTEC":
        value["calibration"]["paired_control"] = True
    if reactor_type == "Electrochemical":
        value.update({"electrolyte_phase": "aqueous", "electrolyte": {
            "identity": "fixture 1 M KOH", "source": "fixture:electrolyte-v1"}})
    _write_json(case / "hydrogen_case.json", value)
    _write_json(case / "hydrogen_validation_records.json", {
        "source": "fixture:synthetic-validation-v1", "records": [
            {"id": "fixture-train-1", "predicted": .1, "observed": .1},
            {"id": "fixture-train-2", "predicted": .2, "observed": .2},
            {"id": "fixture-holdout-1", "predicted": .105, "observed": .1}]})


def _convergence() -> dict:
    """Return mesh-stable, conservative data independently rechecked by code."""
    return {"converged": True,
        "mesh_series": [{"cells": 100, "observable": .9},
                        {"cells": 400, "observable": .99},
                        {"cells": 1600, "observable": 1.}],
        "mesh_tolerance_relative": .02,
        "conservation_budgets": {name: {"inlet": 1., "outlet": 1.}
            for name in ("mass", "carbon", "hydrogen", "energy", "charge")}}


def _coupling_proof(case: Path, candidate_id: str,
                    temperature_K: float) -> dict:
    """Write hashed Cantera exchange records for a two-iteration fixture."""
    files = {"mechanism": ("fixture-mechanism.yaml", "# test fixture mechanism\n"),
        "cantera_log": ("fixture-cantera.log", "test fixture Cantera execution\n"),
        "rate_exchange": ("fixture-rates.json", json.dumps({
            "schema_version": 1, "candidate_id": candidate_id,
            "temperature_K": temperature_K,
            "reaction_rates_mol_m3_s": {"CH4": -.1, "H2": .2}})),
        "coupling_history": ("fixture-coupling-history.json", json.dumps({
            "iterations": [{"iteration": 1, "residual_relative": .01},
                           {"iteration": 2, "residual_relative": 1e-5}]}))}
    proof = {"schema_version": 1, "cantera_used": True,
        "coupling_method": "iterative_two_way", "coupling_iterations": 2,
        "coupling_residual_relative": 1e-5, "coupling_tolerance_relative": 1e-4,
        "exchanged_fields": ["species", "temperature", "reaction_heat",
                             "reaction_rates", "momentum", "charge", "potential"]}
    for stem, (name, content) in files.items():
        path = case / name
        path.write_text(content)
        proof[f"{stem}_path"] = name
        proof[f"{stem}_sha256"] = _sha256(path)
    return proof


class FixtureSolver:
    """Deterministic process-boundary substitute for representative cases."""

    def __init__(self, reactor_type: str, candidate_id: str,
                 temperature_K: float):
        self.reactor_type = reactor_type
        self.candidate_id = candidate_id
        self.temperature_K = temperature_K

    def preflight(self, _mode: str) -> dict:
        """Declare fixture executables without probing the host machine."""
        return {"missing": [], "solvers": {
            "openfoam": {"executable": "fixture-openfoam"},
            "fenicsx": {"available": True}, "cantera": {"available": True}}}

    def execute(self, _command: list[str], case: Path, _timeout: int,
                backend: str = "solver") -> None:
        """Write the same files an external solver must return."""
        request = case / "hydrogen_coupling_request.json"
        iteration = int(json.loads(request.read_text())["iteration"]) if request.exists() else 1
        if backend.startswith("openfoam") and self.reactor_type == "NTEC":
            field = case / "fixture-hydrodynamics.xdmf"
            field.write_text("<Xdmf><!-- deterministic test fixture --></Xdmf>\n")
            _write_json(case / "hydrogen_hydrodynamics.json", {
                "schema_version": 1, "candidate_id": self.candidate_id,
                "pathway_mode": "ntec", "reactor_type": "NTEC",
                "temperature_K": self.temperature_K, "iteration": iteration,
                "fields": {"velocity_m_s": .1, "pressure_Pa": 101325.,
                           "temperature_K": self.temperature_K,
                           "liquid_volume_fraction": .9, "shear_rate_s_inv": 100.},
                "field_artifact": {"path": field.name, "sha256": _sha256(field),
                    "format": "XDMF", "mesh_id": "fixture-mesh-v1",
                    "coordinate_system": "cartesian-m"}})
        if self.reactor_type in {"Fluidized", "MMBCR"}:
            outputs = ({"gas_velocity_m_s": .1, "u_mf_m_s": .02,
                        "bubble_fraction": .2} if self.reactor_type == "Fluidized"
                else {"gas_velocity_m_s": .05, "gas_holdup_fraction": .1,
                      "bubble_diameter_mm": 5.})
        else:
            residual = .01 if iteration == 1 else 1e-5
            feedback = case / "hydrogen_feedback.json"
            _write_json(feedback, {"iteration": iteration, "source": "test-fixture"})
            _write_json(case / "hydrogen_coupling_state.json", {
                "schema_version": 1, "candidate_id": self.candidate_id,
                "reactor_type": self.reactor_type, "temperature_K": self.temperature_K,
                "iteration": iteration, "residual_relative": residual,
                "tolerance_relative": 1e-4, "converged": iteration >= 2,
                "feedback_artifact": {"path": feedback.name,
                                      "sha256": _sha256(feedback)}})
            outputs = ({"CH4_conversion": .2, "H2_selectivity": .9,
                        "solid_C_selectivity": .95,
                        "specific_energy_kWh_kg_H2": 15.}
                if self.reactor_type == "NTEC" else
                {"CH4_conversion": .12, "H2_selectivity": .88,
                 "faradaic_efficiency_H2": .91, "current_density_A_cm2": .2,
                 "cell_voltage_V": 1.4, "electrical_power_density_W_cm2": .28})
        _write_json(case / "hydrogen_outputs.json", outputs)
        _write_json(case / "hydrogen_convergence.json", _convergence())
        if self.reactor_type in {"NTEC", "Electrochemical"}:
            metadata = {"solver_coupling": _coupling_proof(
                case, self.candidate_id, self.temperature_K),
                "mechanism": {"complete": True, "source": "fixture:mechanism-v1"}}
            if self.reactor_type == "Electrochemical":
                metadata["electrolyte_phase"] = "aqueous"
            else:
                metadata["calibration"] = {"paired_control": True,
                    "paired_control_source": "fixture:paired-control-v1"}
            _write_json(case / "hydrogen_metadata.json", metadata)

    def services(self) -> SolverExecutionServices:
        """Expose fixture operations through the production service interface."""
        return SolverExecutionServices(preflight=self.preflight, execute=self.execute,
            openfoam_version=lambda _path: "test-fixture-openfoam",
            fenics_command=lambda script: (["fixture-fenicsx", str(script)],
                                           "test-fixture-fenicsx"),
            cantera_version=lambda: "test-fixture-cantera")
