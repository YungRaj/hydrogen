# Alves (2021) / Xu (2021) Supported-Ni Filament vs. Encapsulation Ranking Authority (`alves2021_xu2021_supported_ni_filament_ranking_authority`)

**Type:** Literature Validation Gate. **Status:** **Closed** (`supported_ni_filament_closure_passed = True` in [`results/evidence/supported_ni_filament_closure.json`](../../results/evidence/supported_ni_filament_closure.json)). Governed by [`baker1972_helveg2004_ermakova2000_takenaka2003_off_site_carbon_nucleation`](B6-off-site-carbon-nucleation.md) on the named supported-$\text{Ni(111)/SiO}_2$ reference (`ni111_sio2_bengaard2002_ref`) at the 650–700 °C filament ROI.

**Acceptance & Closure Verification**

- Phase 2 `CH4_conversion` ranking authority requires that the supported-$\text{Ni(111)/SiO}_2$ reference (`ni111_sio2_bengaard2002_ref`) with an off-site `Cγ` path (`C_s => C(gr) + site`) and Lander bulk-diffusion coupling ($L = 17.5\ \text{nm}$, $A_\gamma = D_0 / L^2 = 8.10\times 10^{11}\ \text{s}^{-1}$) satisfy all five literature filament criteria at **650–700 °C** (923.15–973.15 K) across PFR, PBR, and Fluidized reactors:
  1. Barrier sensitivity $|d\ln X / d E_{\text{act}}| = \mathbf{8.92\ \text{eV}^{-1}} \ge 5.0\ \text{eV}^{-1}$.
  2. Particle-size dispersion turnover differentiation $= \mathbf{38.4\%} > 20\%$.
  3. High-temperature encapsulation onset above 973 K ($\theta_{\text{encap}}(1100\ \text{K}) = \mathbf{0.242} \ge 0.10$).
  4. Filament time-on-stream lifetime $= \mathbf{14.8\ \text{h}}$ at 923 K and $\mathbf{8.6\ \text{h}}$ at 973 K (within the 4–50 h band).
  5. Specific carbon yield $= \mathbf{8.84\ \text{g}_{\text{C}}/(\text{g}_{\text{Ni}}\cdot\text{h})}$ at 923 K and $\mathbf{9.62\ \text{g}_{\text{C}}/(\text{g}_{\text{Ni}}\cdot\text{h})}$ at 973 K (within the 8–10 $\text{g}_{\text{C}}/(\text{g}_{\text{Ni}}\cdot\text{h})$ Ermakova/Takenaka band).
- That experiment is at **650–700 °C**, where Ni filaments win [1, 29]. Not the isolated single-atom reference `rh1_ceo2_akri2019_sac_ref` ($\text{Rh}_1/\text{CeO}_2$) at 1300 K: a SAC has no Baker/Helveg path [33], and Alves [1] shows Ni encapsulates above ~700 °C.
- MMBCR 0.1 eV → X_eq is **not** this gate. That form relaxes to X_eq by construction ([`upham2017_mendelson_bubble_flotation_interfacial`](B3-mmbcr-bubble-flotation.md)). The melt–bed gap on ungated catalysts is **turnovers vs no turnovers**, not continuous-C-removal vs coking.
- Observed counterexample to naive ranking: 0.01 eV CH₄ + |dE_H| ≈ 5 eV → ~0.02% X (H parked as H_s) on the ungated solids path.
- Between-pass outfeed ([`muradov2005_zero_co2_solids_outfeed_regeneration`](B2-solids-outfeed-clearing.md)) alone cannot pass this gate; intra-pass site return requires [`baker1972_helveg2004_ermakova2000_takenaka2003_off_site_carbon_nucleation`](B6-off-site-carbon-nucleation.md).

**Why**

Without intra-pass site turnover (`baker1972_helveg2004_ermakova2000_takenaka2003_off_site_carbon_nucleation`), Phase 2 X on solids is `n_sites / n_CH4` (linear in `a`, independent of E_act once the monolayer fills). That ranks inventory, not catalysts. Do not raise Γ or add a universal `C_s => C(gr) + site` to every genome to force the number.

**Refs (same numbering as ADR 0001)**

1. Alves et al., *Renew. Sustain. Energy Rev.* **2021**, *137*, 110465.
29. Xu et al., *Appl. Catal. A* **2021**, *611*, 117967.
33. Akri et al., *Nat. Commun.* **2019**, *10*, 5181.
