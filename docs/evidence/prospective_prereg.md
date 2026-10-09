# Prospective catalyst-guidance preregistration

This protocol is locked before any new outcomes are generated. Batches 01–09
are exploratory because the policy changed between them. They cannot confirm
this test.

- **Frozen policy:** `policy_source_sha256 = 855fb9816555bace82a3ccb56878645fd2bd0a8713c6a645f45099fa6b35f269`.
- **Pool and budget:** 20 fresh candidates per each of 14 material classes
  (280 total) for turquoise hydrogen. The ORR pool is its PEMFC-eligible
  subset. Each application receives one class-floor slot plus six additional
  slots; the budget is not increased with pool size.
- **Replicates:** exactly 10 completed batches. Failed infrastructure runs may
  be retried under the same locked manifest but are not replaced selectively.
- **Hit definition:** deterministic best `ceil(0.20 * n)` eligible outcomes,
  ordered by metric then candidate ID. Pyrolysis eligibility requires finite
  `E_act`, `valid=True`, `E_act_censored=False`, and
  `pyrolysis_viable=True`.
- **Primary statistic:** summed catalyst-policy hits across all 10 batches,
  compared with 50,000 policy-matched random draws using a one-sided
  permutation p-value, independently for each application.
- **Controls:** summed uncertainty- and validity-policy hits remain reported.
  Catalyst guidance must exceed both controls in each application.
- **Error rate:** family-wise alpha 0.05, Bonferroni split to one-sided
  `p <= 0.025` for each of the two applications.
- **Stopping:** no efficacy stop or policy edit is allowed before all 10
  batches resolve. Safety, corruption, or protocol failures stop the campaign
  without counting as evidence.

The test evaluates computational enrichment under the current eSen screening
protocol. It does not establish DFT accuracy or experimental catalyst
performance.
