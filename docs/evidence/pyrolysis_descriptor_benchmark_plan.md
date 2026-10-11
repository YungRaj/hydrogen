# Pyrolysis descriptor benchmark plan

Status: **dataset curation in progress (3 active, 3 poor exact conditions); not
locked; no descriptor scoring allowed**.

The machine-readable curation file is
`docs/evidence/pyrolysis_benchmark_dataset.json`. The loader in
`pipeline/evidence/pyrolysis_benchmark.py` refuses to score it until every
provenance, independence, exact-encoding, sample-size, and digest gate passes.

## Purpose

The current screening objective is a class-specific BEP transform of methane
split energy. Within a class it is monotone in `dE_CH3 + dE_H`, so it can reward
strong binding even when hydrogen removal or carbon transport limits turnover.
Before a two-sided energetic-span descriptor can steer discovery, it must
recover externally labelled active and poor/poisoning catalyst families better
than the current objective.

## Dataset contract

- Include at least 12 independently sourced active conditions and 12
  independently sourced poor, deactivating, or poisoning conditions.
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
