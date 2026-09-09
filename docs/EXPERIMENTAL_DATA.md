# Experimental Calibration and Holdout Data

This contract makes experimental evidence reproducible and prevents calibration
data, replicates, or paired controls from leaking into reported holdout results.
It does not manufacture measurements or decide an apparatus design.

## Experimental design

Assign each physical catalyst batch, MEA, or paired NTEC/control preparation an
`experimental_unit_id`. Assign that unit to `calibration` or `holdout` before a
holdout measurement begins. Every repeat from the same unit remains in that
split. Holdout records set `blinded: true` and preserve the timestamp and source
of the assignment. Analysts should not inspect holdout outcomes while fitting
kinetics, transport closures, surrogate models, thresholds, or preprocessing.

Use independent catalyst preparations—not repeated instrument reads from one
preparation—as the primary independent replicates. Record failed runs and their
failure reasons in the checksum-bound `run_log_path`; do not silently replace
them. A checksum-bound preregistered analysis plan declares the minimum number
of independent calibration and holdout units. The plan must be locked before
holdout measurements are made.
Randomization, blocking, power analysis, instrument calibration, detection
limits, censoring rules, and uncertainty propagation belong in the referenced
`protocol_id` and raw source artifacts.

## Required record types

| Type | Required purpose and outputs |
|---|---|
| `pyrolysis_reactor` | CH4 conversion, H2 selectivity, solid-carbon yield, carbon/hydrogen balance closure, net energy per kg H2, and deactivation rate under recorded T, P, complete feed, duration, apparatus, and candidate identity |
| `ntec_pair` | The same outcomes for exactly one NTEC treatment and one control using distinct experimental units but sharing candidate, apparatus, split, T, P, feed, and duration |
| `mea_performance` | Peak power, system efficiency, and ORR overpotential with temperature, pressure, anode/cathode humidity, duration, and MEA identity |
| `durability` | Voltage degradation and power retention over a stated duration |
| `hydrogen_impurity` | Power and retention under an explicit nonempty impurity mole-fraction composition |

Every observation has an uncertainty entry with the identical key and declares
whether it is a standard deviation, standard error, or expanded uncertainty.
Fractions
use `[0, 1]`; temperature is K, pressure Pa, flow mol/s, duration h, power
W/cm2, degradation microvolt/h, and energy kWh/kg-H2. The field names encode
these units. Conversion and yield must use carbon-inclusive inlet/outlet
bookkeeping; solid-carbon yield cannot exceed methane conversion.

## Raw evidence and identity

Each summarized record references an immutable raw instrument export using
`raw_source_path` and `raw_source_sha256`. The validator resolves relative paths
from the dataset directory and verifies the bytes. `candidate_id`,
`apparatus_id`, `protocol_id`, `measurement_id`, `experimental_unit_id`, and
`replicate_id` must refer to laboratory records, not display labels invented
after analysis.

Validate a dataset with:

```bash
python -m pipeline.evidence.experimental_dataset \
  data/experimental_dataset.json \
  --output results/experimental_dataset_validation.json \
  --evidence-output results/experimental_evidence_records.json
```

The report lists independent calibration and holdout units for every
candidate/type group. A group becomes `split_complete` only when both counts
meet the preregistered minima. Schema validity therefore cannot be mistaken for
adequate sample size or scientific acceptance.
The optional evidence output contains manifest-ready records only for
split-complete groups. It never promotes a calibration-only dataset.

## Minimal structure

See `docs/experimental_dataset.example.json`. Its paths and checksums are
illustrative and must be replaced with real immutable raw exports. For an NTEC
pair, use two distinct experimental units, add the same `control_pair_id` to
both records, and set
`conditions.treatment` to `ntec` and `control`, respectively.
The control declares zero shear/power intervention; the NTEC record declares
positive `shear_rate_s_inv` and `mechanical_power_W`.

## Relationship to model validation

This dataset establishes trustworthy observations. Predictions are generated
separately by the frozen computational model. The existing
`pipeline.process.model_validation.score_holdout` routine compares those
predictions with declared holdout IDs and computes the error itself. Predictions
must never be written into the raw experimental source or used to change the
holdout assignment.

Only a checksum-verified dataset validation report should be attached to the
project evidence manifest. Passing schema validation proves provenance and
split integrity; it does not prove that the experiment was unbiased, the model
was calibrated, or performance targets were achieved.
