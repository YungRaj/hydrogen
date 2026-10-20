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
  the largest invalidity family. This motivated the current v4 protocol. The
  completed confirmatory verdict remains the historical v3 campaign and is not
  silently reclassified.
- **Pyrolysis redesign:** the external benchmark is complete and locked at 12
  independent primary sources per label. Correct surface-layer composition
  scaling and five-replica configurational averaging made all 18 unique genomes
  and all 24 conditions scoreable. The current raw barrier achieved AUC 0.549;
  the proposed energetic-span score achieved AUC 0.299, an improvement of
  −0.250, and failed every acceptance
  threshold. Neither descriptor is externally validated, and no new pyrolysis
  campaign is authorized. The result and model boundary are
  in `docs/evidence/pyrolysis_descriptor_benchmark_result.json` and
  `docs/evidence/pyrolysis_descriptor_benchmark_plan.md`.
- **Pyrolysis rate/free-energy follow-up:** rate averaging and the NIST 1-bar
  H2 chemical-potential correction are implemented and locked, but the inherited
  molten BEP model predicts only 1.29–1.63 eV for Cu45Bi55 versus published NEB
  barriers of 2.53–2.60 eV. The independent DFT gate therefore failed before
  any new held-out experimental conditions were scored. No campaign is open.
- **Pyrolysis MLIP-NEB follow-up:** a checksum-locked eSEN climbing-image NEB
  check attempted five Cu45Bi55 configurations against Palmer's 2.53–2.60 eV
  barriers. Zero paths passed the frozen convergence and continuity gates: one
  met the 0.08 eV/A force target but was discontinuous, while four exhausted
  400 steps or remained discontinuous. MLIP-NEB is therefore not validated as
  a barrier replacement, the held-out conditions remain untouched, and no
  pyrolysis campaign is authorized.
- **Molten-alloy NEB endpoint-alignment follow-up:** `methane_endpoints` in
  `pipeline/evidence/molten_neb.py` now relaxes the product (`*CH3 + *H`) first
  and constructs the reactant (`CH4(g) + melt`) from the relaxed product slab
  with `FixAtoms` on the melt and rigid `CH3` lift, guaranteeing
  `max_surface_displacement_A == 0.0` across all amorphous seeds before NEB
  interpolation.
- **Pyrolysis interfacial carbon displacement follow-up:** across all 24
  external primary-literature conditions (`docs/evidence/pyrolysis_benchmark_dataset.json`,
  12 independent sources per label), zero-parameter interfacial carbon
  displacement `|c_height_above_surface_A|` achieved AUC 0.9931 (`95% CI
  [0.9580, 1.0000]`) on eSEN and AUC 0.9583 (`95% CI [0.8667, 1.0000]`) on
  UMA-s-1p1, while the coupled activation + interfacial retention score
  `interfacial_pyrolysis_score_eV = max(E_act_raw, 0) + |c_height_above_surface_A|`
  achieved AUC 0.9722 (`95% CI [0.8929, 1.0000]`, +0.4236 over `E_act_raw`),
  clearing all preregistered benchmark gates (`AUC >= 0.70`, `ci95_low > 0.50`,
  `auc_improvement >= +0.10`) and wired into secondary objective training and
  reactor stage selection.
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
- **Ni(111) NEB workflow check:** the corrected 11-image staged path was tested
  against Bengaard's 1.05 eV reference. Both endpoints and plain NEB converged,
  but climbing exhausted 1,000 steps at 0.557 eV/A with a discontinuous path.
  This is a workflow failure, not a GNN accuracy result; molten re-testing and
  the held-out set remain blocked.
- **Ni(111) interpolation follow-up:** replacing IDPP with linear interpolation
  and adding a constrained reactant endpoint reduced the final path force to
  0.449 eV/A, but plain NEB still exhausted 1,000 steps and the largest image
  jump remained 3.85 A. The setup is therefore still not validated.
- **Ni(111) endpoint-alignment follow-up:** the reactant is now constructed from
  the relaxed product with only C and the breaking H changed, plus an explicit
  atom-matching gate. The product endpoint failed to reconverge under this
  stricter construction; that intermediate blocked result was superseded by the
  corrected replay below.
- **Ni(111) corrected replay:** after attaching the missing product calculator,
  rigidly lifting CH3, and rebuilding the reactant from the relaxed product,
  endpoint matching passed, the path was continuous, and both NEB stages
  converged at 0.0798 eV/A. eSEN predicted 0.386 eV versus Bengaard's 1.05 eV
  (0.664 eV absolute error), so workflow validity passed but model accuracy
  failed. Molten and held-out scoring remain blocked.
- **Claim boundary:** screening outcomes are eSen GNN estimates, not DFT or
  experimental validation.
- **Readiness:** the six-point campaign status remains fail-closed until its
  required calibrated and held-out evidence is present.
- **Preservation:** tag `orr-confirmatory-v1` identifies the confirmed ORR
  result and its checksum-locked DFT slate.
- **Supported nanoparticles:** the opt-in genome records diameter and loading for
  B6 studies and now varies slab capillary/loading strain, dispersion, and
  sub-5 nm step-edge vacancy morphology with `diameter_nm` and `loading_wt_pct`;
  experimental calibration is still required before treating them as validated
  physical controls.
