# Representative screening operating points (assumed, uncalibrated)

These files hold the **assumed** geometry, operating conditions, material
properties and model parameters used by the uncalibrated NTEC and
electrochemical MOOSE + Cantera screening tier
(`pipeline/simulation/screening_cases.py` ->
`pipeline/simulation/external_runner.run_screening_case`).

| File | Reactor | Electrolyte |
| --- | --- | --- |
| `ntec.json` | NTEC liquid-metal shear layer | - |
| `electrochemical_aqueous.json` | Electrochemical cell | aqueous (1 M KOH), used below 600 K |
| `electrochemical_molten.json` | Electrochemical cell | molten (NaOH-KOH eutectic), used at/above 600 K |

## What these values are

- They are representative values chosen so that the repository MOOSE models
  (`moose/inputs/*.i`) solve a physically plausible case. They are **not
  measurements** and have **not been calibrated** against any experiment.
- `temperature_K` is deliberately absent: the campaign temperature grid
  supplies it per condition (`pipeline.reactors.modes.reactor_temperatures_for_mode`).
- The only candidate-specific input is the Cantera interface mechanism written
  from the candidate's screening kinetics (`write_full_mechanism`).

Results produced from these files are labelled
`uncalibrated_multiphysics_screening`, carry
`candidate_exclusion_authorized: false`, and are never written as production
multiphysics artifacts. They rank candidates under one shared assumed
condition; they do not predict absolute reactor performance.

## Replacing them with measured conditions

1. Copy the relevant file and replace each value with a measured or
   literature value.
2. Replace the matching `parameter_sources` entry (`section.field`) with the
   citation or measurement identifier. Every `geometry`, `operating` and
   `properties` field must keep a source; `REPLACE_*` placeholders are rejected.
3. Update `feed.source` and, for electrochemistry, `electrolyte.identity` /
   `electrolyte.source`.
4. For a calibrated (production) result, do not edit these files: build a
   full case with a `calibration` block and held-out validation records
   (`python -m pipeline.simulation.physical_case init ...`) and run
   `pipeline/simulation/external_runner.py`. Only that path produces
   artifacts that the reactor stage treats as calibrated evidence.
