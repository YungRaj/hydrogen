# Muradov (2005) / Pinilla (2011) Zero-CO₂ Mechanical Solids Outfeed & Regeneration (`muradov2005_zero_co2_solids_outfeed_regeneration`)

**Acceptance**

- Mechanical / consumable outfeed is the default solids carbon policy. Oxidative regen remains locked unless `co2_permitted=True`.
- Circulating fluidized carbon removal acts **during** integration, not after `net.advance()`. **Shipped:** `simulate_fluidized_bed` substeps and applies `circulating_carbon_removal_rate_1_s` between them. PFR is discrete regen only.
- Reported metrics distinguish single-pass X from cycle-averaged production.
- Default `max_regen_cycles >= 1` is not sufficient by itself: outfeed restores sites **between** passes; inventory **during** a pass is [`deutschmann2000_monolayer_site_inventory`](B1-solid-site-inventory.md). Intra-pass site return is [`baker1972_helveg2004_ermakova2000_takenaka2003_off_site_carbon_nucleation`](B6-off-site-carbon-nucleation.md) (`Cα → Cγ` chemistry), not a coverage reset.

**Why**

`max_regen_cycles = 0` was a side effect of blocking CO₂ burnoff. Fixed beds clog and ΔP rises [1,2]. Fluidized / circulating beds exist so carbon-laden particles can leave while stripped particles return [1,2]. Pinilla et al. [4] is the fixed-vs-fluidized **deactivation** comparison (higher X and faster coke in a fixed bed at the same WHSV), not the existence proof for circulating beds.

Between-pass solids outfeed (`muradov2005_zero_co2_solids_outfeed_regeneration`) does **not** by itself close [`alves2021_xu2021_supported_ni_filament_ranking_authority`](B5-phase2-x-gate-before-dft.md). Decoking restores sites for the next pass. The ungated solids identity `X ≈ n_sites / n_CH4` is **within** a pass: the lattice is consumed unless off-site filament nucleation ([`baker1972_helveg2004_ermakova2000_takenaka2003_off_site_carbon_nucleation`](B6-off-site-carbon-nucleation.md)) is active. Circulating “removal” is a coverage lump, not Baker/Helveg transport.

**Refs**

1. Alves et al., *Renew. Sustain. Energy Rev.* **2021**, *137*, 110465.
2. Sánchez-Bastardo et al., *Ind. Eng. Chem. Res.* **2021**, *60*, 11855–11881.
4. Pinilla et al., *Int. J. Hydrogen Energy* **2012**.
6. Abdollahi et al., *Energies* **2024**, *17*, 290 (melt skim — analogous continuous C takeoff).
