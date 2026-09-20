# Phase 2 backlog

Issue bodies are self-contained. File them on the working fork (`naryandesai/hydrogen`). `gh` is installed but this machine is not logged in (`gh auth login`).

Reference numbering matches [ADR 0001](../adr/0001-pyrolysis-phase-admissibility.md).

| ID | Title | Type |
|---|---|---|
| [B1](B1-solid-site-inventory.md) | Solid site inventory without inventing Γ | implementation (done) |
| [B2](B2-solids-outfeed-clearing.md) | Solids outfeed / clearing during produce | implementation |
| [B3](B3-mmbcr-bubble-flotation.md) | MMBCR k·a·τ is a sanity check, not kinetic closure | docs + tests |
| [B4](B4-pfr-residence-time-limits.md) | PFR τ limits | implementation |
| [B5](B5-phase2-x-gate-before-dft.md) | **Gate:** do not send Phase 2 X to DFT until PFR/fluidized closure | milestone blocker |
| [B6](B6-off-site-carbon-nucleation.md) | Off-site carbon nucleation (Cα → Cγ / Cδ) | implementation (**blocks B5**; B6-1–6 done, PFR + Fluidized Ni; Cδ ∝ θ_C²; five B6-6 gates pass; judge still `cat_9`; Ni surrogate row open) |

After `gh auth login`:

```powershell
gh milestone create "Phase 2 solids closure" --repo naryandesai/hydrogen --description "B6 (off-site C nucleation) then B5: PFR/fluidized must approach X_eq on a supported-Ni-like catalyst before Phase 2 X is used for DFT. B2 does not close this."
gh issue create --repo naryandesai/hydrogen --title "B1: Solid-catalyst site inventory without inventing Γ" --body-file docs/backlog/B1-solid-site-inventory.md --label "enhancement"
gh issue create --repo naryandesai/hydrogen --title "B2: Solids produce → outfeed/clear → return" --body-file docs/backlog/B2-solids-outfeed-clearing.md --label "enhancement"
gh issue create --repo naryandesai/hydrogen --title "B3: MMBCR k·a·τ is a sanity check, not kinetic closure" --body-file docs/backlog/B3-mmbcr-bubble-flotation.md --label "documentation"
gh issue create --repo naryandesai/hydrogen --title "B4: Practical limits on lengthening the PFR" --body-file docs/backlog/B4-pfr-residence-time-limits.md --label "enhancement"
gh issue create --repo naryandesai/hydrogen --title "B5: Do not send Phase 2 X to DFT until PFR/fluidized closure" --body-file docs/backlog/B5-phase2-x-gate-before-dft.md --label "gate" --milestone "Phase 2 solids closure"
gh issue create --repo naryandesai/hydrogen --title "B6: Off-site carbon nucleation (Cα → Cγ / Cδ)" --body-file docs/backlog/B6-off-site-carbon-nucleation.md --label "enhancement" --milestone "Phase 2 solids closure"
```
