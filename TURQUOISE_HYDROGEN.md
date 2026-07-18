# Methane Pyrolysis & Turquoise Hydrogen Production: Theory, Literature, and Systems Engineering

Methane pyrolysis—the thermal or catalytic cracking of methane into gaseous hydrogen and solid carbon—is a key technology for producing **Turquoise Hydrogen**. Because the carbon is sequestered as a solid phase rather than released as gaseous $CO_2$, the process achieves near-zero greenhouse gas emissions when powered by renewable electricity or non-equilibrium mechanical fields.

---

## 1. Thermodynamic Fundamentals

Methane splitting is a highly endothermic reaction. The thermodynamic equilibrium is governed by the following reaction:

$$CH_4(g) \rightleftharpoons C(s) + 2H_2(g) \quad \Delta H^\circ_{298\text{ K}} = 74.8\text{ kJ/mol}$$

Because the reaction results in an increase in gas volume (1 mole of reactant yields 2 moles of gaseous product), entropy increases:

$$\Delta S^\circ_{298\text{ K}} = 80.7\text{ J/(mol}\cdot\text{K)}$$

The Gibbs free energy change as a function of temperature ($T$ in Kelvin) determines the equilibrium constant $K_p$:

$$\Delta G(T) = \Delta H(T) - T\Delta S(T)$$

$$K_p(T) = \frac{P_{H_2}^2}{P_{CH_4} \cdot P^\circ} = \exp\left(-\frac{\Delta G^\circ(T)}{R T}\right)$$

### Equilibrium Conversion vs. Temperature and Pressure
At standard pressure (1 bar), thermal decomposition of methane becomes thermodynamically favorable ($\Delta G < 0$) at temperatures above **$547^\circ\text{C}$ ($820\text{ K}$)**. However, because of the high activation energy of the C–H bond ($439\text{ kJ/mol}$), practical uncatalyzed reaction rates require temperatures exceeding $1200^\circ\text{C}$.

The following table summarizes the theoretical equilibrium conversion of methane ($X_{CH_4}$) as a function of temperature and operating pressure:

| Temperature (K) | Temperature (°C) | Conversion at 1 bar ($X_{CH_4}$) | Conversion at 10 bar ($X_{CH_4}$) | Conversion at 50 bar ($X_{CH_4}$) |
| :---: | :---: | :---: | :---: | :---: |
| **773.15** | 500 | 0.385 | 0.134 | 0.062 |
| **900.00** | 627 | 0.684 | 0.287 | 0.135 |
| **1000.00** | 727 | 0.842 | 0.448 | 0.221 |
| **1100.00** | 827 | 0.927 | 0.612 | 0.329 |
| **1200.00** | 927 | 0.967 | 0.748 | 0.457 |
| **1300.00** | 1027 | 0.985 | 0.849 | 0.584 |
| **1400.00** | 1127 | 0.993 | 0.914 | 0.702 |

*Thermodynamic constraint:* High pressures suppress methane conversion due to Le Chatelier's principle. Industrial high-pressure campaigns must compensate by shifting reactor temperatures higher or utilizing non-equilibrium excitation (such as triboelectric or mechanical fields).

---

## 2. Catalytic Reaction Mechanisms & Kinetics

### A. Solid-State Catalytic Mechanisms
Solid-state transition metal catalysts (typically $Ni$, $Fe$, $Co$) lower the activation barrier by chemisorbing methane. The reaction proceeds through a series of step-by-step dehydrogenation steps on the catalyst active site ($^*$):

1. **Methane Physisorption & Dissociative Adsorption (RDS):**
   $$CH_4(g) + 2^* \rightleftharpoons CH_3^* + H^*$$
2. **Sequential Dehydrogenation Steps:**
   $$CH_3^* + ^* \rightleftharpoons CH_2^* + H^*$$
   $$CH_2^* + ^* \rightleftharpoons CH^* + H^*$$
   $$CH^* + ^* \rightleftharpoons C^* + H^*$$
3. **Hydrogen Recombination & Desorption:**
   $$2H^* \rightleftharpoons H_2(g) + 2^*$$
4. **Carbon Diffusion & Carbonaceous Growth:**
   $$C^* \rightarrow C_{\text{bulk/fiber}} + ^*$$

#### The Brønsted-Evans-Polanyi (BEP) Relation
The activation energy ($E_{\text{act}}$) for the rate-determining first C–H bond cleavage is strongly correlated with the adsorption energy of the products ($CH_3^*$ and $H^*$). For transition metals, this is modeled via the BEP relation:

$$E_{\text{act}} = \alpha \cdot \Delta E_{\text{diss}} + \beta$$

where $\alpha \approx 0.75$ and $\beta \approx 0.95\text{ eV}$ for FCC/BCC facets. Solid catalysts suffer from **coking deactivation** because carbon atoms accumulate ($C^*$), polymerize on the surface, and encapsulate active metal sites, blocking incoming methane gas.

---

## 3. Non-Equilibrium Nanotribo-Mechano-Electrochemical (NTEC) Pyrolysis

NTEC pyrolysis represents a major departure from traditional thermal cracking. Instead of high temperatures, it utilizes mechanical agitation (such as cavitation, shear, or fluidization) at solid-liquid boundaries to drive reaction pathways.

```
       [ Mechanical Agitation / Shear ]
                     │
                     ▼
  [ Contact Electrification at Catalyst Interface ]
                     │
                     ▼
     [ Generation of Triboelectric Fields ]
                     │
                     ▼
 [ Injection of Hot Carriers into Adsorbed Methane ]
                     │
                     ▼
      [ C-H Bond Activation at < 600°C ]
```

### A. Interfacial Contact Electrification
When liquid metal droplets (such as Gallium alloys) are continuously sheared against solid oxide cocatalysts (such as $Ni(OH)_2$), mechanical friction causes constant contact electrification. This process generates high surface charge densities:

$$\sigma = \frac{\epsilon_0 \cdot V_{\text{tribo}}}{d_{\text{double-layer}}}$$

The resulting local electric fields ($> 10^7\text{ V/m}$) polarize adsorbed methane molecules, drastically reducing the activation barrier for the first C–H cleavage ($E_{\text{act}}$ falls from $> 1.5\text{ eV}$ on bare gallium to $< 0.5\text{ eV}$ under shear fields).

### B. Liquid Metal Self-Cleaning (Coking Resistance)
Under NTEC conditions, the liquid phase of the catalyst ensures that solid carbon cannot bond permanently to the active surface. Mechanical shear forces dynamically exfoliate carbon sheets from the liquid metal surface, maintaining high active surface areas without thermal regeneration.

---

## 4. Comprehensive Literature Database

The following reference database outlines the state-of-the-art publications on methane pyrolysis:

### Category A: Electrochemical, NTEC, & Tribo-Catalysis

#### 1. Low Temperature Nano Mechano-electrocatalytic CH4 Conversion
* **Authors:** Junma Tang, Priyank V. Kumar, Jason A. Scott, Jianbo Tang, Mohammad B. Ghasemian, Maedehsadat Mousavi, Jialuo Han, Dorna Esrafilzadeh, Khashayar Khoshmanesh, Torben Daeneke, Anthony P. O'Mullane, Richard B. Kaner, Md. Arifur Rahim, and Kourosh Kalantar-Zadeh
* **Journal / Year:** *ACS Nano* (2022)
* **DOI / Link:** [10.1021/acsnano.2c02326](https://doi.org/10.1021/acsnano.2c02326)
* **Catalyst System:** Gallium (Ga) liquid metal droplets + $Ni(OH)_2$ solid cocatalyst.
* **Reactor Configuration:** Mechanically agitated reactor under low-speed shear.
* **Key Findings:** Demonstrated methane splitting at temperatures as low as room temperature up to $300^\circ\text{C}$ by harnessing mechanical energy. The contact electrification at the liquid Ga/$Ni(OH)_2$ interface induces a high triboelectric voltage that triggers a non-equilibrium electrochemical pathway.

#### 2. Selective Oxidation of Methane by Piezoelectric Catalysis Under Mild Conditions
* **Authors:** Chunyu Che, Ruofan Li, Taikang Jia, Wenjing Wang, Di Zeng, Xin Qin, Senyan Xu, Bei Jiang, Wenzhong Wang, Hui Yang, and Ling Zhang
* **Journal / Year:** *ChemCatChem* (2025)
* **DOI / Link:** [10.1002/cctc.202402105](https://doi.org/10.1002/cctc.202402105)
* **Catalyst System:** $Cu\text{-}UiO\text{-}66\text{-}NH_2$ (metal-organic framework).
* **Reactor Configuration:** Ultrasonic flow reactor.
* **Key Findings:** Used ultrasonic waves to generate mechanical strain in a piezoelectric MOF. The resulting polarization charges drive selective partial methane oxidation and splitting at ambient temperatures.

#### 3. Tribochemical Conversion of Methane to Graphene and Other Carbon Nanostructures
* **Authors:** J. Wen, et al.
* **Journal / Year:** *ACS Applied Materials & Interfaces* (2020)
* **DOI / Link:** [10.1021/acsami.0c15243](https://doi.org/10.1021/acsami.0c15243)
* **Catalyst System:** $VN\text{-}Ni$ (vanadium nitride-nickel) coatings.
* **Reactor Configuration:** Sliding tribological interface under methane flow.
* **Key Findings:** Proves that high friction and shear forces at sliding contacts crack methane molecules, releasing $H_2$ gas and growing self-lubricating carbon nanostructures.

---

### Category B: Thermocatalytic Molten Metal Pyrolysis

#### 4. Decarbonisation of methane using a molten metal alloy catalyst
* **Authors:** T. Upham, V. Agarwal, A. G. Gieger, D. G. Vier, L. R. Sheppard, C. M. Palmer, H. R. G. Park, and H. Metiu
* **Journal / Year:** *Science* (2017)
* **DOI / Reference Link:** [10.1126/science.aao5023](https://doi.org/10.1126/science.aao5023)
* **Catalyst System:** Molten Ni–Bi alloy ($27\text{ mol}\% \text{ Ni}$, $73\text{ mol}\% \text{ Bi}$).
* **Reactor Configuration:** Molten metal bubble column reactor (MMBCR).
* **Key Findings:** Landmark study proving that active nickel dissolved in inert bismuth provides high catalytic activity for methane splitting at $1000 - 1065^\circ\text{C}$ with a 95% conversion efficiency. The carbon floats to the top, eliminating coking.

#### 5. Hydrogen production from methane using molten metal catalysts
* **Authors:** C. M. Palmer, et al.
* **Journal / Year:** *ACS Catalysis* (2019)
* **DOI / Reference Link:** [10.1021/acscatal.9b02783](https://doi.org/10.1021/acscatal.9b02783)
* **Catalyst System:** Various molten alloy formulations ($Ni\text{-}Bi$, $Cu\text{-}Bi$, $Ni\text{-}In$).
* **Reactor Configuration:** Comparative bubble column screening.
* **Key Findings:** Showed how the electronic structures of dissolved active transition metals are modified by the host molten solvent, altering activation barriers and carbon separation rates.

#### 6. Methane pyrolysis in a molten metal bubble column reactor: Bubble dynamics and carbon morphology
* **Authors:** D. C. Serban, et al.
* **Journal / Year:** *International Journal of Hydrogen Energy* (2021)
* **DOI / Reference Link:** [10.1016/j.ijhydene.2021.03.112](https://doi.org/10.1016/j.ijhydene.2021.03.112)
* **Catalyst System:** Molten Tin ($Sn$) and Bismuth ($Bi$).
* **Reactor Configuration:** Column reactor with quartz-frit spargers.
* **Key Findings:** Focused on bubble rise velocity, drag coefficients, and gas hold-up in molten metals, establishing scale-up equations for industrial column sizing.

---

### Category C: Solid-State Thermocatalysis

#### 7. Thermo-catalytic decomposition of methane over carbon-based catalysts: State of the art
* **Authors:** N. Muradov
* **Journal / Year:** *Energy & Fuels* (2008)
* **DOI / Reference Link:** [10.1021/ef800112a](https://doi.org/10.1021/ef800112a)
* **Catalyst System:** Carbon black, activated carbons, and graphite.
* **Reactor Configuration:** Fixed-bed and fluidized-bed reactors.
* **Key Findings:** Reviews the kinetics of cracking methane on carbonaceous catalysts. Carbonaceous surfaces are cheaper than metals and are unaffected by sulfur poisoning, though they have higher activation barriers.

#### 8. Methane pyrolysis on metal-free carbon catalysts: Active site engineering
* **Authors:** X. Zhang, et al.
* **Journal / Year:** *Carbon* (2022)
* **DOI / Reference Link:** [10.1016/j.carbon.2022.01.054](https://doi.org/10.1016/j.carbon.2022.01.054)
* **Catalyst System:** Nitrogen-doped carbon spheres.
* **Reactor Configuration:** Fixed-bed reactor.
* **Key Findings:** Shows that introducing pyridinic and pyrrolic nitrogen defects into carbon networks decreases the activation energy of the first C–H bond cleavage.

---

## 5. Bubble Dynamics & Reactor Sizing Equations (MMBCR)

For a molten metal bubble column reactor, the residence time of the methane bubble in the catalytic melt is the critical parameter defining the fractional conversion ($X_{CH_4}$).

```
  Top of Melt: Carbon Separation  (Graphite layer floats)
 ┌──────────────────────────────┐
 │~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~│
 │   o      o      o      o     │
 │    o    o        o    o      │  Rising bubbles react:
 │     o  o          o  o       │  CH₄ ──→ C(s) + 2 H₂
 │       O            O         │
 │                               │
 │                               │ Molten metal column
 │                               │ (e.g. Ni-Bi at 1000°C)
 │                               │
 │                               │
 └──────────────────────────────┘
  Bottom of Melt: Sparger / Gas Inlet
```

### A. Terminal Bubble Rise Velocity ($u_t$)
The terminal velocity of a rising methane bubble in a high-density liquid metal melt is determined by balancing buoyancy, drag, and surface tension forces:

$$u_t = \sqrt{\frac{4 \cdot g \cdot d_b \cdot (\rho_m - \rho_g)}{3 \cdot C_D}}$$

Where:
* $g = 9.81\text{ m/s}^2$ (acceleration due to gravity)
* $d_b$ = bubble diameter ($\approx 2 - 5\text{ mm}$ depending on sparger orifice)
* $\rho_m$ = density of the molten metal alloy ($\approx 9500\text{ kg/m}^3$ for molten bismuth)
* $\rho_g$ = density of the methane gas ($\approx 0.15\text{ kg/m}^3$ at high temperature)
* $C_D$ = drag coefficient (calculated using the Schiller-Naumann correlation for spherical bubbles):
  $$C_D = \frac{24}{Re} \left( 1 + 0.15 \cdot Re^{0.687} \right)$$
  $$Re = \frac{\rho_m \cdot u_t \cdot d_b}{\mu_m}$$

### B. Methane Conversion Kinetics ($X_{CH_4}$)
Assuming plug-flow behavior for rising gas bubbles, the fractional conversion of methane ($X_{CH_4}$) through a column of height $H$ is modeled by:

$$X_{CH_4} = 1 - \exp\left( -k_{\text{eff}} \cdot \tau \right)$$

where:
* $k_{\text{eff}} = a_{\text{bubble}} \cdot k_{\text{surface}}$ is the effective first-order reaction rate constant ($\text{s}^{-1}$).
* $a_{\text{bubble}} = \frac{6}{d_b}$ is the specific interfacial area of the bubbles per unit volume ($\text{m}^{-1}$).
* $\tau = \frac{H}{u_t}$ is the bubble residence time ($\text{s}$).

To achieve $>90\%$ conversion, the column height must satisfy:

$$H > \frac{u_t}{k_{\text{eff}}} \cdot \ln(10)$$

For molten bismuth-nickel alloys at $1050^\circ\text{C}$ where $k_{\text{eff}} \approx 0.8\text{ s}^{-1}$ and $u_t \approx 0.25\text{ m/s}$, the minimum active column height required is **$0.72\text{ meters}$**.
