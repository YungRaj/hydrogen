# Project status

- **Discovery policy:** catalyst-guided, coverage-safe branch-and-bound is the
  production standard. Uncertainty and validity are controls and validation
  signals, not production branch priorities.
- **Prospective evidence:** batches 01–09 are exploratory. After excluding
  censored and non-viable pyrolysis outcomes, pooled hits are pyrolysis
  14 catalyst / 10 uncertainty / 18 validity / 15.90 random mean, and ORR
  40 / 36 / 23 / 24.73. Only ORR exceeds policy-matched random in this
  retrospective analysis.
- **Next campaign:** the frozen 35-batch protocol screens only ORR, with 10
  candidates per class and a 20-candidate policy budget. CPU replay including
  the completed geometry smoke estimates 98.9% ORR pass probability; the
  conservative estimate excluding oversized batch 04 is 83.5%. Do not edit
  the policy while it runs.
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
