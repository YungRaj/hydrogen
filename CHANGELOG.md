# Changelog

## Unreleased

### Added

- Six explicit methane-conversion pathway modes, with thermocatalytic operation
  as the default.
- Portable OpenFOAM, FEniCSx, and Cantera solver discovery and execution.
- Strict external-result artifacts gated by identity, backend versions,
  convergence, mesh independence, conservation, physical bounds, provenance,
  and held-out model error.
- Mode-specific `hydrogen_case.json` validation for fluidized-bed, MMBCR, NTEC,
  and electrochemical calculations, including per-parameter sources and
  disjoint calibration/validation records.
- Reproducible `openfoam-env` and `fenicsx-env` definitions.
- Aqueous/molten electrochemical configuration and paired-control NTEC evidence
  boundaries.

### Changed

- Packed-bed PFR geometry distinguishes total catalyst-bed area from void gas
  residence volume.
- Fluidized screening uses a reacting emulsion plus conservative bubble bypass
  populated by validated OpenFOAM hydrodynamics.
- MMBCR screening uses physical column flow, gas holdup, bubble diameter, and
  interfacial area populated by validated OpenFOAM hydrodynamics.
- Missing specialized physics remains visible and non-excluding instead of
  silently falling back to an unrelated reactor.
- Production and orchestrator CLIs propagate pathway and multiphysics artifact
  configuration end to end.

### Validation

- 74 pipeline tests, 34 scientific contracts, and 24 exclusion-audit checks
  pass locally.
- Real Cantera chemistry loading, OpenFOAM executable startup, and a FEniCSx
  finite-element solve were exercised successfully.

### Known scientific boundaries

- PFR is operational for qualified screening, not experimentally validated
  prediction; incomplete elementary kinetics remain explicitly labeled.
- Fluidized, MMBCR, NTEC, and electrochemical modes require real case inputs,
  candidate mechanisms, calibration data, and held-out validation.
- The legacy ideal-gas carbon tracer remains pending the separate graphite and
  MMBCR work intended for integration from the collaborating fork.
