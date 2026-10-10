# Prospective catalyst-guidance preregistration

This protocol is locked before any new outcomes are generated. Batches 01–09
are exploratory because the policy changed between them. They cannot confirm
this test. ORR is the sole confirmatory endpoint; pyrolysis is reported as an
exploratory endpoint because its 30-batch CPU replay estimated only 4.2%
probability of passing this design.

- **Frozen policy:** the enforced digest is stored in
  `docs/evidence/prospective_policy_lock.json`. `prepare()` fails before writing
  a manifest if any frozen source differs.
- **Pool and budget:** generate 20 fresh candidates per each of 14 material
  classes, then apply application scope before screening. Pyrolysis excludes
  COF, MOF, MXene, MetalHydride, and Perovskite under ADR 0001; ORR excludes
  MetalHydride and MoltenMetal. Each application evaluates exactly 20
  candidates, retaining one floor slot for every eligible class and reallocating
  the remaining slots by policy score.
- **Replicates:** exactly 30 completed batches. Failed infrastructure runs may
  be retried under the same locked manifest but are not replaced selectively.
- **Hit definition:** deterministic best `ceil(0.20 * n)` eligible outcomes,
  ordered by metric then candidate ID. Pyrolysis eligibility requires finite
  `E_act`, `valid=True`, `E_act_censored=False`, and
  `pyrolysis_viable=True`.
- **Primary statistic:** summed ORR catalyst-policy hits across all 30 batches,
  compared with 50,000 policy-matched random draws using a one-sided
  permutation p-value.
- **Controls:** summed uncertainty- and validity-policy hits remain reported.
  Confirmatory ORR catalyst guidance must exceed both controls.
- **Error rate:** one confirmatory endpoint, one-sided `alpha = 0.05`.
- **Stopping:** no efficacy stop or policy edit is allowed before all 30
  batches resolve. Safety, corruption, or protocol failures stop the campaign
  without counting as evidence.

Pyrolysis uses the same frozen policy and is always reported, but it cannot
pass or fail the confirmatory campaign. With 30 batches, the CPU chronological
replay estimated 98.4% pass probability for ORR and 4.2% for pyrolysis. The
conservative ORR estimate with oversized batch 04 removed was 83.3%; at only
10 batches it was 57.9%. These estimates use 2,000 bootstraps and remain
indicative because historical pools were smaller than the planned pools.

The test evaluates computational enrichment under the current eSen screening
protocol. It does not establish DFT accuracy or experimental catalyst
performance.
