# Pyrolysis descriptor benchmark plan

Status: **benchmark executed and failed its descriptor acceptance criteria; no
pyrolysis campaign is authorized**.

## Result (2026-10-10)

The locked run screened all 18 unique genomes represented by the 24 independent
experimental conditions. A molten-alloy construction defect had scaled atomic
percent against all 36 slab atoms while placing every substitution into the
9-atom surface layer. Correcting that in the canonical implementation made the
concentrated Ni33Sn67 condition physically scoreable and cleared the coverage
gate: all 12 active and all 12 poor sources are eligible.

Five fixed-topology replicas showed that a single alloy termination is noisy:
for Bi–Ni 27 at% at 1273 K, the raw barrier spanned −1.336 to 1.595 eV and
`dE_H` spanned −3.255 to 0.663 eV. Temperature had also been part of the
substitution seed, confounding comparisons. The canonical screener now keeps
topology fixed across temperature and averages five deterministic
configurational snapshots, requiring at least four valid relaxations and
retaining their spread. Short MD was not substituted for this ensemble because
it would not sample promoter rearrangements on a screening trajectory.

The frozen analysis then rejected the proposed descriptor. The current raw
barrier achieved AUC 0.549 (grouped-bootstrap 95% CI 0.292–0.800); the
energetic-span score achieved AUC 0.299 (0.081–0.543), an improvement of
−0.250. It therefore missed the required AUC, lower-confidence-bound, and
improvement thresholds.

The apparent `dE_H ≈ −5 eV` boundary is not a clamp: pure Te replicas converge
there with negligible spread. Only the derived operational barrier is clamped.
Positive `dE_C` follows from the declared thermochemical carbon reference
`E(CH4) − 2E(H2)` and is therefore carbon formation relative to methane and
hydrogen, not isolated-atom adsorption. CH3 remained intact in every accepted
ensemble; carbon relaxed below the top surface plane in 10 of 18 candidates,
confirming that the carbon term overlaps the unresolved B6 transport physics.

The immutable summary is
`docs/evidence/pyrolysis_descriptor_benchmark_result.json`; its raw input is
`results/screening/pyrolysis_benchmark.csv`. This is a useful negative result: the
current GNN/BEP descriptors do not reliably recover external activity labels.
The result rules out deploying the tested energetic-span score as the
pyrolysis discovery objective; it does not validate the raw barrier, whose AUC
is indistinguishable from random ranking on this benchmark.

The machine-readable curation file is
`docs/evidence/pyrolysis_benchmark_dataset.json`. The loader in
`pipeline/evidence/pyrolysis_benchmark.py` refuses to score it until every
provenance, independence, exact-encoding, sample-size, and digest gate passes.
Admitted and rejected primary sources, including the exact reason for every
rejection, are recorded in `docs/evidence/pyrolysis_source_audit.md`.

## Purpose

The current screening objective is a class-specific BEP transform of methane
split energy. Within a class it is monotone in `dE_CH3 + dE_H`, so it can reward
strong binding even when hydrogen removal or carbon transport limits turnover.
Before a two-sided energetic-span descriptor can steer discovery, it must
recover externally labelled active and poor/poisoning catalyst families better
than the current objective.

## Dataset contract

- Include at least 12 active conditions and 12 poor, deactivating, or
  poisoning conditions from 12 distinct primary publications per label.
- Treat conditions from one paper and catalyst family as one resampling group;
  compositional variants from the same series are not independent replicates.
- Each row requires catalyst identity, temperature, reactor context, activity
  label, label basis, DOI or stable source URL, and one exact representable
  genome. Approximate encodings must state the mismatch and cannot enter the
  primary analysis.
- Experimental condition IDs must be unique. The same catalyst may appear in
  independently sourced conditions; candidate identity is not mistaken for
  experimental-condition identity.
- Supported-Ni/Fe/Co, molten Ni-Bi/Cu-Bi/Ni-In, carbon catalysts, and explicit
  non-catalytic or deactivating controls should be represented where comparable
  evidence exists. Phase-inadmissible encodings remain excluded under ADR 0001.
- Labels must be fixed from external evidence before either descriptor is
  evaluated. Screening results cannot determine or revise a label.

Tracked starting references include Upham et al. 2017
([DOI 10.1126/science.aao5023](https://doi.org/10.1126/science.aao5023)),
Palmer et al. 2019
([DOI 10.1021/acscatal.9b01833](https://doi.org/10.1021/acscatal.9b01833)),
Xu et al. 2021
([DOI 10.1016/j.apcata.2020.117967](https://doi.org/10.1016/j.apcata.2020.117967)),
and the mechanistic carbon-transport sources in ADR 0001 and backlog B6. These
references identify candidate families; they do not by themselves assign every
needed binary benchmark label.

## Frozen analysis once the dataset is complete

- Primary statistic: grouped-bootstrap ROC AUC for known-active versus
  poor/poisoning conditions, with 10,000 resamples by source-family group.
- Comparator: the current uncensored `E_act`, evaluated on exactly the same
  rows. Lower descriptor values rank as more active.
- Minimum success rule: energetic-span AUC point estimate at least 0.70, its
  grouped-bootstrap 95% lower bound above 0.50, and an AUC improvement of at
  least 0.10 over current `E_act`.
- Missing, invalid, censored, or descriptor-uncertain rows are reported and
  fail the benchmark if they reduce either label group below 12 independent
  source-family conditions.
- Secondary outputs: limiting-step confusion by catalyst family and a
  descriptive re-ranking of historical batches. Neither is confirmatory.

## Locking sequence

1. Commit the complete cited dataset and its SHA-256 manifest.
2. Commit the versioned descriptor parameter table with a citation or explicit
   uncertainty flag for every parameter.
3. Record both digests and the implementation digest in a preregistration lock.
4. Only then compute either descriptor's benchmark AUC.

No pyrolysis campaign should be prepared unless this benchmark passes and a
separate replay estimates at least 80% power under the proposed campaign design.
