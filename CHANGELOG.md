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
- Checksum- and geometry-bound candidate Hamiltonian construction from standard
  FCIDUMP integrals using PySCF and OpenFermion.
- Pre-write validation for curated, checksum-bound prior-art manifests with an
  explicit blinded time split.
- A portable CUDA-Q/PySCF/OpenFermion environment definition.
- A checksum-bound experimental dataset contract for reactor, paired NTEC,
  MEA, durability, and hydrogen-impurity calibration/holdout evidence, including
  preregistered sample minima and experimental-unit leakage prevention.
- Strict OpenFOAM-to-FEniCSx field handoffs and Cantera coupling receipts,
  including mechanism/log/rate/history hashes and two-way residual convergence.
- Runner-owned, numbered OpenFOAM/FEniCSx outer iterations whose observed
  feedback states must match the final coupling receipt.
- Pristine-case enforcement and external FEniCSx script hashing for reproducible
  multiphysics input identities.

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

- 74 pipeline tests, 41 scientific contracts, and 24 exclusion-audit checks
  pass locally.
- Real Cantera chemistry loading, OpenFOAM executable startup, and a FEniCSx
  finite-element solve were exercised successfully.
- External modes now include safe case-template generation, batch readiness,
  independently computed raw holdout error, quantitative mesh refinement,
  recomputed elemental/energy/charge balances, phase consistency, physical
  output identities, backend-specific logs, and an explicit NTEC handoff.

### Known scientific boundaries

- PFR is operational for qualified screening, not experimentally validated
  prediction; incomplete elementary kinetics remain explicitly labeled.
- Fluidized, MMBCR, NTEC, and electrochemical modes require real case inputs,
  candidate mechanisms, calibration data, and held-out validation.
- The legacy ideal-gas carbon tracer remains pending the separate graphite and
  MMBCR work intended for integration from the collaborating fork.
