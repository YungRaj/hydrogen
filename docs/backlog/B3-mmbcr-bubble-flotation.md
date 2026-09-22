# B3 — MMBCR: bubble area + flotation (not a Langmuir slab)

**Status:** Rate form shipped. Not kinetic closure.

**Acceptance (remaining)**

- Interfacial `k0` and `a = 6/d` are documented as a **calibrated melt-side prefactor**, not a DFT or Upham-fitted mass-transfer coefficient. Validator: `0 < k0 ≤ 1` m/s. Default `0.01` m/s (Upham 2017; Chen 2023; Abdollahi 2024).
- Flotation frequency `mmbcr_carbon_removal_rate_1_s` is wired: `None` = unconstrained (production), `0` = fouled interface, finite = `η = k_float / (k_float + k_if · a)` on the ODE. It is not a Langmuir site rate.
- Tests assert: X ≤ X_eq always; X → X_eq only when `k·a·τ` is large; X is kinetics-limited at Chen-scale E_a (~0.84 eV) on the current column.
- Do **not** treat “E_act = 0.1 eV → X = X_eq” as evidence the kinetics are right. That is a property of `dX/dt ∝ (X_eq − X)`.
- Cross-reactor comparison: MMBCR is pinned toward X_eq while PFR/fluidized have no intra-pass turnovers, so the melt will systematically look better until **B6** gives the bed a Cγ path. That gap is turnovers vs no turnovers, not flotation vs coking. B2 (between-pass outfeed) does not close it.

**Shipped 2026-09-06**

`simulate_mmbcr` uses `k(E_act,T) × a_bubble × (X_eq − X)`; carbon is not a site occupant. Sanity: 0.1 eV, 1300 K, k0 = 0.01 m/s → X = 0.985 = X_eq.

**Refs**

5. Upham et al., *Science* **2017**, *358*, 917–921.
6. Abdollahi et al., *Energies* **2024**, *17*, 290.
7. Chen et al., *Science* **2023**.
8. Meloni et al., *J. Anal. Appl. Pyrolysis* **2024**.
9. Le et al., *AIChE J.* **2024**.
10. Rahimi et al., *Carbon* **2019**, *151*, 181–191.
