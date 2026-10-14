# Physical Multiphysics Cases

An OpenFOAM or MOOSE process exit code is not sufficient evidence. Every
external reactor calculation must contain `hydrogen_case.json` in its case
directory. The runner validates this file before launching a solver and hashes
the pristine directory so the numerical result is tied to the exact inputs.

## Shared contract

Every case declares:

- `schema_version`, canonical `candidate_id`, `pathway_mode`, and
  `reactor_type`;
- `geometry`, `operating`, and `properties` values in the units encoded in
  their field names;
- a nonempty model choice for each required closure;
- feed composition and its source;
- kinetic/mechanism source;
- a source for every numerical parameter in `parameter_sources`;
- nonempty and disjoint calibration and holdout-validation experiment IDs;
- the calibration dataset source.

Training and holdout IDs may not overlap. Placeholder values and templates are
never runnable. NTEC additionally requires a paired
control declaration. Electrochemical cases additionally select `aqueous` or
`molten`. Missing, non-finite, zero/negative, unsourced, mismatched, or leaking
inputs stop the external calculation before compute is spent.

## Required fields by reactor

| Reactor | Geometry | Operating values | Properties | Closures |
|---|---|---|---|---|
| Fluidized | column diameter, bed height | T, P, methane mass flow | particle diameter/density, gas viscosity | drag, heat transfer |
| MMBCR | column diameter, liquid height, sparger orifice | T, P, methane mass flow | liquid density/viscosity, surface tension | bubble breakup/coalescence |
| NTEC | reactor volume, interface area | T, P, methane flow, shear rate, mechanical power | liquid viscosity, permittivity, ionic conductivity | contact electrification, species transport |
| Electrochemical | electrode area, electrolyte thickness | T, P, methane flow, applied potential | ionic/electronic conductivity, methane diffusivity | charge transfer, species transport |

These are the minimum admissibility fields, not a claim that every problem is
fully identified. A case may and normally should contain additional parameters
such as walls, thermal properties, diffusion coefficients, turbulence models,
porosity, tortuosity, particle distributions, solubilities, and boundary
conditions. Case-specific solver files remain authoritative for discretization.

## MOOSE model inputs (NTEC and electrochemical)

Both pathways run repository MOOSE templates (`moose/inputs/*.i`). Every
finite number in `geometry`, `operating`, `properties` and `model_parameters`
is passed to MOOSE under its own unit-bearing name, so these names are part of
the contract. `kinetics` must also name the candidate Cantera interface
mechanism in the case directory (`mechanism`, `gas_phase`, `interface_phase`).

| Model | Physics | `model_parameters` |
|---|---|---|
| Electrochemical (`electrochemical_cell.i`, 1-D) | Ohm's law in the electrolyte (excess supporting electrolyte), CH4 diffusion, concentration-dependent Butler-Volmer at both electrodes; anode exchange current exchanged with Cantera | `methane_henry_mol_m3_Pa`, `anode_equilibrium_potential_V`, `cathode_equilibrium_potential_V`, `anode_transfer_coefficient`, `cathode_transfer_coefficient`, `cathode_exchange_current_A_m2`, `anode_exchange_current_guess_A_m2`, `ch4_current_fraction`, `cathode_h2_faradaic_efficiency`, `electrons_per_ch4`, `mesh_elements` |
| NTEC (`ntec_shear_layer.i`, 2-D) | Incompressible Navier-Stokes shear layer, Poisson in the dielectric with contact charge sigma, CH4 diffusion with a first-order surface sink; sigma(wall shear) and the field-enhanced rate constant exchanged with Cantera | `liquid_density_kg_m3`, `film_thickness_m`, `dielectric_thickness_m`, `methane_diffusivity_m2_s`, `methane_henry_mol_m3_Pa`, `contact_charge_saturation_C_m2`, `contact_shear_scale_s_inv`, `field_barrier_coefficient_eV_m_V`, `field_barrier_cap_eV`, `interface_charge_guess_C_m2`, `interface_rate_constant_guess_m_s`, `adhesion_energy_J_m2`, `effective_mass_kg`, `particle_diameter_m`, `contact_area_m2`, `mesh_elements_x`, `mesh_elements_film`, `mesh_elements_dielectric` |

Optional for both: `coupling_tolerance_relative` (default 1e-6),
`coupling_relaxation` (default 1.0) and `mesh_tolerance_relative` (default
0.01). The NTEC contact-electrification and field-barrier parameters are
hypothesis closures: they must be fitted to paired-control data, and results
remain non-excluding until the held-out check passes.

`python tests/run_tests.py --include-moose` runs both templates with the real
`moose-opt` and checks them against closed-form solutions (limiting current,
Butler-Volmer, ohmic drop, Couette shear, Gauss's law, viscous dissipation,
Robin surface concentration) and mesh refinement.

## Example skeleton

```json
{
  "schema_version": 1,
  "candidate_id": "canonical-candidate-id",
  "pathway_mode": "thermocatalytic_fluidized",
  "reactor_type": "Fluidized",
  "geometry": {"column_diameter_m": 0.1, "bed_height_m": 0.8},
  "operating": {
    "temperature_K": 900.0,
    "pressure_Pa": 101325.0,
    "methane_mass_flow_kg_s": 0.0001
  },
  "properties": {
    "particle_diameter_m": 0.001,
    "particle_density_kg_m3": 2500.0,
    "gas_viscosity_Pa_s": 0.00002
  },
  "models": {
    "drag_model": "Gidaspow",
    "heat_transfer_model": "Ranz-Marshall"
  },
  "feed": {"composition": {"CH4": 0.95, "Ar": 0.05}, "source": "experiment DOI or dataset record"},
  "kinetics": {"source": "candidate mechanism digest or publication"},
  "parameter_sources": {
    "geometry.column_diameter_m": "apparatus drawing",
    "geometry.bed_height_m": "apparatus drawing",
    "operating.temperature_K": "calibrated thermocouple record",
    "operating.pressure_Pa": "pressure-transducer record",
    "operating.methane_mass_flow_kg_s": "MFC calibration record",
    "properties.particle_diameter_m": "sieve/laser measurement",
    "properties.particle_density_kg_m3": "pycnometry record",
    "properties.gas_viscosity_Pa_s": "property database citation"
  },
  "calibration": {
    "training_ids": ["experiment-001", "experiment-002"],
    "validation_ids": ["experiment-101"],
    "source": "immutable dataset path or DOI",
    "metric": "relative_rmse",
    "acceptance_threshold": 0.10
  }
}
```

## Solver outputs and acceptance

The case must emit `hydrogen_outputs.json`,
`hydrogen_convergence.json`, and, where required,
`hydrogen_metadata.json`. The runner accepts the result only when the artifact
passes identity, required-backend, version, output, numerical-bound,
convergence, mesh-independence, conservation, physical-case, calibration, and
provenance checks. The case must also emit raw prediction/observation records
in `hydrogen_validation_records.json`. The runner—not the external model—then
calculates `model_validation.metric`, `holdout_error`, `acceptance_threshold`,
and `passed`; acceptance requires a finite held-out error no larger than the
declared threshold. It likewise recomputes mesh stability from at least three
successively refined meshes and recomputes mass, carbon, hydrogen, energy,
and/or charge closure from inlet/outlet budgets appropriate to the mode.
Failure removes the would-be accepted artifact. Downstream
screening records `validation_required` and retains the candidate.

## Solver handoffs and coupling proof

An NTEC OpenFOAM stage must write `hydrogen_hydrodynamics.json`. The runner
validates schema version, candidate, mode, reactor, temperature, finite
unit-bearing summary fields (`velocity_m_s`, `pressure_Pa`, `temperature_K`,
`liquid_volume_fraction`, and `shear_rate_s_inv`), mesh identity, coordinate
system, spatial-field format, and the checksum of an in-case XDMF/HDF5/VTU
field artifact. Existence alone is not accepted.

NTEC and electrochemical metadata must replace a boolean Cantera assertion with
a `solver_coupling` proof containing relative in-case paths and SHA-256 hashes
for the exact mechanism, nonempty Cantera log, candidate-specific reaction-rate
exchange, and coupling-iteration history. The rate exchange identifies the
candidate and temperature and records finite `reaction_rates_mol_m3_s`.
Production artifacts require `coupling_method=iterative_two_way`, at least two
sequential iterations, a finite final residual no larger than its declared
tolerance, and exchange of species, temperature, reaction heat, reaction rates,
charge, and potential; NTEC additionally exchanges momentum.

The runner owns this outer loop. Before each pass it writes
`hydrogen_coupling_request.json` with the candidate, reactor, temperature, and
iteration number; it then reruns MOOSE and the repository Cantera interface
step (`pipeline/electrochemistry/cantera_coupling.py`), validates
`hydrogen_coupling_state.json` and its checksummed feedback artifact, and stops
only after convergence or the configured maximum. The final coupling receipt
must match the iteration count, residual, and tolerance observed by the runner.

The MOOSE templates in `moose/inputs` consume the exchanged values as
command-line overrides; the model-specific physics and its analytic
verification are documented in those files and `tests/test_moose_models.py`. The repository
verifies their inputs, outputs, observed history, residual, and conservation
evidence; it cannot infer correct boundary conditions or material laws from a
log file.

Before hashing or execution, the runner rejects old hydrogen outputs/logs and
nonzero OpenFOAM time directories. OpenFOAM's `0` initial-condition directory
remains valid input. The MOOSE input file outside the case directory is hashed as an
explicit external input, and its digest is retained in the accepted artifact.

This contract establishes reproducible computational machinery. Predictive
validation still requires real measurements and candidate-specific kinetics;
the repository intentionally cannot synthesize either as ground truth.

Use `python -m pipeline.simulation.physical_case init ...` to generate a complete
key skeleton and `python -m pipeline.simulation.physical_case validate ...` for a
standalone preflight. Use
`python -m pipeline.simulation.case_preparation MANIFEST.json --create` to
create missing batch skeletons and report why each case is not ready.

A batch manifest has this minimal form:

```json
{
  "cases": [
    {
      "candidate_id": "CANONICAL_ID",
      "mode": "mmbcr",
      "reactor_type": "MMBCR",
      "temperature_K": 900.0,
      "case_dir": "cases/CANONICAL_ID/mmbcr"
    }
  ]
}
```
