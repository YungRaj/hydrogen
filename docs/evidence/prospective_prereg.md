# Prospective catalyst-guidance preregistration

This protocol is locked before any new outcomes are generated. Batches 01–09
are exploratory because the policy changed between them. They cannot confirm
this test. ORR is the sole confirmatory endpoint; pyrolysis is reported as an
exploratory endpoint because its 30-batch CPU replay estimated only 4.2%
probability of passing this design.

- **Frozen policy:** the enforced digest is stored in
  `docs/evidence/prospective_policy_lock.json`. `prepare()` fails before writing
  a manifest if any frozen source differs.
- **Pool and budget:** generate 10 fresh candidates for each of 14 material
  classes. The confirmatory campaign screens only ORR, whose scope excludes
  MetalHydride and MoltenMetal, so each batch evaluates up to 120 candidates.
  The 20-candidate policy selection retains one floor slot for every eligible
  class and reallocates the remaining slots by policy score. The entire eligible
  pool is screened because its outcomes define both the deterministic top-20%
  hit set and the policy-matched random null.
- **Replicates:** exactly 35 completed batches. Failed infrastructure runs may
  be retried under the same locked manifest but are not replaced selectively.
- **Hit definition:** deterministic best `ceil(0.20 * n)` eligible outcomes,
  ordered by metric then candidate ID. Pyrolysis eligibility requires finite
  `E_act`, `valid=True`, `E_act_censored=False`, and
  `pyrolysis_viable=True`.
- **Primary statistic:** summed ORR catalyst-policy hits across all 35 batches,
  compared with 50,000 policy-matched random draws using a one-sided
  permutation p-value.
- **Controls:** summed uncertainty- and validity-policy hits remain reported.
  Confirmatory ORR catalyst guidance must exceed both controls.
- **Error rate:** one confirmatory endpoint, one-sided `alpha = 0.05`.
- **Stopping:** no efficacy stop or policy edit is allowed before all 35
  batches resolve. Safety, corruption, or protocol failures stop the campaign
  without counting as evidence.

Pyrolysis is excluded from confirmatory batches and may be evaluated only in
separately labelled exploratory or smoke work. The batch count is fixed from a
CPU chronological replay that includes the completed geometry smoke and must
provide at least 80% estimated ORR pass probability when oversized historical
batch 04 is removed. At 35 batches, the full-data estimate is 98.9% and the
leave-batch-04-out estimate is 83.5%. These estimates use 2,000 bootstraps and
remain indicative because historical pools were smaller than the planned pools.

The policy is sequentially adaptive: every finalized batch available before a
new manifest is prepared contributes eligible outcomes to that manifest's
training data. This includes historical, smoke, exploratory, and earlier
confirmatory batches. A manifest is immutable after preparation, so outcomes
from its own batch and later batches cannot influence its selections.

The test evaluates computational enrichment under the current eSen screening
protocol. It does not establish DFT accuracy or experimental catalyst
performance.

## Locked result

All 35 batches completed under policy digest
`bd68883dcbb4d45d24462bcefe84e85fbd8a97aeb62c99cc5812c3df8b109920`.
Catalyst guidance found 289 pooled hits, uncertainty ranking found 149,
validity ranking found 153, and policy-matched random selection averaged
120.2459. The one-sided random p-value was 0.0000199996 at alpha 0.05, so the
preregistered ORR endpoint passed. The immutable machine-readable result is
`results/prospective_search/confirmatory_analysis.json`.
