# Project status

- **Discovery policy:** catalyst-guided, coverage-safe branch-and-bound is the
  production standard. Uncertainty and validity are controls and validation
  signals, not production branch priorities.
- **Confirmatory evidence:** the locked 35-batch ORR campaign passed. Catalyst
  guidance found 289 pooled hits versus 149 for uncertainty, 153 for validity,
  and a policy-matched random mean of 120.25. The preregistered one-sided
  random p-value is 0.000020. Catalyst-guided, coverage-safe branch-and-bound
  is therefore the evidence-backed production search policy for ORR.
- **Project focus:** turquoise hydrogen is the primary application. The ORR
  result validates the search machinery; further ORR validation is deferred
  while the pyrolysis descriptor and reactor evidence are hardened.
- **ORR quality diagnostic:** across all 4,200 confirmatory rows, 16.0% were
  invalid and 32.27% of valid rows were capped at 3.0 V. OOH non-convergence is
  the largest invalidity family. This is descriptive evidence for a future
  protocol revision; the locked v4 policy and verdict remain unchanged.
- **Pyrolysis redesign:** the external benchmark is complete and locked at 12
  independent primary sources per label. Correct surface-layer composition
  scaling made all 18 unique genomes and all 24 conditions scoreable. The
  current raw barrier achieved AUC 0.500; the proposed energetic-span score
  achieved AUC 0.333, an improvement of −0.167, and failed every acceptance
  threshold. Neither descriptor is externally validated, and no new pyrolysis
  campaign is authorized. The result and model boundary are
  in `docs/evidence/pyrolysis_descriptor_benchmark_result.json` and
  `docs/evidence/pyrolysis_descriptor_benchmark_plan.md`.
- **Pyrolysis selection:** reactor handoffs and headline barriers use only
  finite, valid, uncensored, phase-viable screening rows. Unresolved or
  censored rows remain eligible for higher-fidelity rescue, never reactor
  ranking.
- **Geometry smoke:** `smoke-relax-v4-20261009` completed with SAC/DAC validity
  of 4/4 and 3/4 for pyrolysis and 4/4 and 4/4 for ORR. There were no
  `atomic_overlap` failures; the sole DAC rejection was an out-of-bounds
  adsorption energy caught by the physical validity guard.
- **Phase 2:** B5 remains closed pending B6. Solids conversion has no authority
  to rank catalysts, choose DFT work, or choose the experimental handoff slate.
- **Claim boundary:** screening outcomes are eSen GNN estimates, not DFT or
  experimental validation.
- **Readiness:** the six-point campaign status remains fail-closed until its
  required calibrated and held-out evidence is present.
- **Preservation:** tag `orr-confirmatory-v1` identifies the confirmed ORR
  result and its checksum-locked DFT slate.
