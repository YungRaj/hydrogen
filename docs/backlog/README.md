# Phase 2 Literature Validation Gates

Reference numbering matches [ADR 0001](../adr/0001-pyrolysis-phase-admissibility.md). The six literature-backed validation gates (`LITERATURE_VALIDATION_GATES` in [`pipeline/search/discovery.py`](../../pipeline/search/discovery.py)) and named reference catalysts (`NAMED_REFERENCE_CATALYSTS`) are summarized below:

| Validation Gate ID | Literature Name & Primary References | Status |
|---|---|---|
| [`deutschmann2000_monolayer_site_inventory`](B1-solid-site-inventory.md) | **Deutschmann (2000) Monolayer Active-Site Inventory & Geometric Bed Dispersion** (`Γ = 2.5×10⁻⁹ mol/cm²`, `a = a_geom × loading × dispersion`) | **closed** |
| [`muradov2005_zero_co2_solids_outfeed_regeneration`](B2-solids-outfeed-clearing.md) | **Muradov (2005) / Pinilla (2011) Zero-CO₂ Mechanical Solids Outfeed & Regeneration** | **closed** |
| [`upham2017_mendelson_bubble_flotation_interfacial`](B3-mmbcr-bubble-flotation.md) | **Upham (2017) / Mendelson (1967) Molten-Metal Bubble Hydrodynamics & Carbon Flotation** | **closed** |
| [`ergun1952_packed_bed_pressure_drop_pilot_residence`](B4-pfr-residence-time-limits.md) | **Ergun (1952) / Froment & Bischoff (2011) Packed-Bed Pressure-Drop & Residence Envelope** | **closed** |
| [`alves2021_xu2021_supported_ni_filament_ranking_authority`](B5-phase2-x-gate-before-dft.md) | **Alves (2021) / Xu (2021) Supported-Ni Filament vs. Encapsulation Ranking Authority** | **closed** (`supported_ni_filament_closure_passed = True`) |
| [`baker1972_helveg2004_ermakova2000_takenaka2003_off_site_carbon_nucleation`](B6-off-site-carbon-nucleation.md) | **Baker (1972) / Helveg (2004) / Ermakova (2002) / Takenaka (2003) Off-Site Carbon Filament Nucleation vs. Encapsulation** (`Cα → Cγ / Cδ`) | **closed** (`ni111_sio2_bengaard2002_ref` at 650–700 °C; `results/evidence/supported_ni_filament_closure.json`) |

## Named Literature Reference Catalysts (`NAMED_REFERENCE_CATALYSTS`)

- **`ni111_sio2_bengaard2002_ref`** — **Supported $\text{Ni(111)/SiO}_2$ Nanoparticle Catalyst** (Bengaard et al., *J. Catal.* 2002, 209, 365; Snoeck et al., *J. Catal.* 1997, 169, 250; Ermakova et al., *Catal. Today* 2002, 77, 225; Takenaka et al., *J. Catal.* 2003, 220, 468): primary solids filament reference at 650–700 °C (923.15–973.15 K).
- **`rh1_ceo2_akri2019_sac_ref`** — **$\text{Rh}_1/\text{CeO}_2(111)$ Single-Atom Reference Catalyst** (Akri et al., *Nat. Commun.* 2019, 10, 5181): isolated single-atom negative control for off-site carbon filament nucleation (`off_site_carbon_included = False`).
