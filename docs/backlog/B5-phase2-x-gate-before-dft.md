# B5 — Do not send Phase 2 conversion to DFT until solids closure passes

**Type:** Gate. Blocks the DFT / Phase 3 milestone. **Blocked by [B6](B6-off-site-carbon-nucleation.md).** B6-7 moved the named solids judge to `ni_np_lit` at the 650–700 °C ROI; that does **not** lift this gate. Ni at 923 K is 14.6 % / 16.9 % (PFR / Fluidized), not X_eq.

**Acceptance**

- Phase 2 `CH4_conversion` is not used to pick DFT candidates until a **supported-Ni-like** catalyst (nanoparticle metal with an off-site Cγ path, modest H desorption) approaches tabulated X_eq on **PFR and fluidized**, where the Langmuir surface is actually integrated.
- That experiment is at **650–700 °C**, where Ni filaments actually win [1, 29]. Not `cat_9` (SAC Rh) at 1300 K: a SAC has no Baker/Helveg path [33], and Alves [1] says Ni at that T encapsulates.
- MMBCR 0.1 eV → X_eq is **not** this gate. That form relaxes to X_eq by construction (B3). The melt–bed gap is **turnovers vs no turnovers**, not continuous-C-removal vs coking.
- Header language: **closure sanity-checked on MMBCR; pending on PFR and fluidized (needs B6).**
- Observed counterexample to ranking: 0.01 eV CH₄ + |dE_H| ≈ 5 eV → ~0.02% X (H parked as H_s) on the old solids path.
- B2 / decoke will **not** pass this gate. Outfeed restores sites between passes; the monolayer identity is within a pass (B6).

**Why**

Until the surface path can turn sites over inside a pass (B6), Phase 2 X on solids is `n_sites / n_CH4` (linear in `a`, independent of E_act once the monolayer fills). That ranks inventory, not catalysts. Do not raise Γ or add a universal `C_s => C(gr) + site` to every genome to force the number.

**Refs (same numbering as ADR 0001)**

1. Alves et al., *Renew. Sustain. Energy Rev.* **2021**, *137*, 110465.
29. Xu et al., *Appl. Catal. A* **2021**, *611*, 117967.
33. Akri et al., *Nat. Commun.* **2019**, *10*, 5181.
