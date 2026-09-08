# Physical Multiphysics Cases

An OpenFOAM or FEniCSx process exit code is not sufficient evidence. Every
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

This contract establishes reproducible computational machinery. Predictive
validation still requires real measurements and candidate-specific kinetics;
the repository intentionally cannot synthesize either as ground truth.

Use `python -m pipeline.process.physical_case init ...` to generate a complete
key skeleton and `python -m pipeline.process.physical_case validate ...` for a
standalone preflight. Use
`python -m pipeline.process.multiphysics_prepare MANIFEST.json --create` to
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
