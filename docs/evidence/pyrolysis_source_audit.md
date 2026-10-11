# Pyrolysis benchmark source audit

This audit records why a primary experimental source is admitted to, or held
out of, the locked descriptor benchmark. It prevents an inaccessible catalyst
identity, reactor effect, or approximate genome from silently becoming a
training label.

## Admitted conditions

| Primary source | Conditions admitted | Evidence boundary |
|---|---:|---|
| [Upham et al. (2017)](https://doi.org/10.1126/science.aao5023) | 2 | Ni17Bi83 active and pure Bi poor in the same 1000 °C differential comparison |
| [Palmer et al. (2019)](https://doi.org/10.1021/acscatal.9b01833) | 2 | Cu45Bi55 active and pure Bi inactive in the reported molten-alloy comparison |
| [Zeng et al. (2020)](https://doi.org/10.1021/acscatal.0c00805) | 1 | Liquid Te converts methane at 1173 K while the catalyst-free control does not |
| [Msheik et al. (2022)](https://doi.org/10.1016/j.energy.2022.124943) | 1 | Pure Sn has negligible catalytic activity in a matched gas-phase comparison |
| [Plevan et al. (2015)](https://doi.org/10.1016/j.ijhydene.2015.04.062) | 1 | Pure Sn has no catalytic effect relative to the matched blank tube |
| [Zaghloul et al. (2021)](https://doi.org/10.1002/ceat.202100210) | 2 | Exact 5 wt% Ni-Sn (9.62 at% Ni) is the best tested medium; pure Sn is the lowest-performing matched medium |
| [Wi et al. (2024)](https://doi.org/10.1007/s40831-024-00962-8) | 2 | Corrected measurements put pure Ga at the fast end and pure Sn near the slow end after subtracting the hydrogen-noise baseline |
| [Baltrusaitis et al. (2023)](https://doi.org/10.1021/acscatal.3c02955) | 1 | The exact single-phase Ni33Sn67 control has a 355 kJ/mol barrier; the active two-phase material remains excluded |
| [Kim et al. (2025)](https://doi.org/10.1016/j.cej.2025.161304) | 2 | Zn45Bi55 is the top screened alloy (157.0 kJ/mol); matched pure Bi is poor (308.7 kJ/mol) |
| [Song et al. (2026)](https://doi.org/10.1016/j.cej.2026.177345) | 1 | Ga11Bi89 sustains hydrogen production for 158 h at 800 °C with a 67.2 kJ/mol barrier |
| [Rahimi et al. (2019)](https://doi.org/10.1016/j.carbon.2019.05.041) | 1 | An exact single-phase Ni27Bi73 control independently produces hydrogen and recoverable solid carbon |
| [Son et al. (2025)](https://doi.org/10.1016/j.apcatb.2024.125009) | 1 | Exact Ni27Bi73 is an active measured baseline for the Se-promoted alloy series at 1000 °C |
| [Chen et al. (2025)](https://doi.org/10.1021/acscatal.5c02867) | 1 | Exact Ni27Bi73 independently gives 8.2% conversion and 4.02 mmol g-1 h-1 hydrogen at 850 °C |
| [Liu et al. (2026)](https://doi.org/10.1039/d6sc05383a) | 1 | Exact Ni27Bi73 independently gives 33.3% conversion at 1000 °C in the quaternary-alloy baseline |
| [Wang et al. (2008)](https://doi.org/10.1016/j.molcata.2007.12.018) | 1 | Pure molten Mg gives approximately 20% methane conversion at 700 °C |
| [Serban et al. (2003)](https://doi.org/10.1021/ef020271q) | 1 | Pure molten Pb gives only approximately 4–10% methane conversion at 750 °C |
| [Geißler et al. (2016)](https://doi.org/10.1016/j.cej.2016.04.066) | 1 | Pure Sn at 930 °C remains far below equilibrium in a packed bubble column and is used as heat-transfer medium, not an active catalyst |
| [Neuschitzer et al. (2023)](https://doi.org/10.3390/en16207058) | 1 | Ascending-temperature trials over pure Sn show no methane conversion at the exact 900 °C condition |
| [Busillo et al. (2024)](https://doi.org/10.1007/s11705-024-2437-x) | 1 | The exact 950 °C pure-Sn condition belongs to a direct molten-versus-empty-reactor comparison and is poor under its short residence time |

These are twelve active and twelve poor conditions from twelve distinct active
and twelve distinct poor primary publications. Multiple conditions from one
publication do not increase the independent-source count.

## Audited but not admitted

| Primary source | Decision | Reason |
|---|---|---|
| [Wang et al. (2008)](https://doi.org/10.1016/j.molcata.2007.12.018) | Admit after narrow schema review | Mg is now an allowed molten host because the primary experiment supplies an exact in-band condition; no ternary or morphology approximation was introduced. |
| [Leal Pérez et al. (2021)](https://doi.org/10.1016/j.ijhydene.2020.11.079) | Hold | The Ga experiment is at 1119 °C, above the 1300 K turquoise screening band. |
| [Baltrusaitis et al. (2023)](https://doi.org/10.1021/acscatal.3c02955) | Partial | The exact liquid-only Ni33Sn67 control is admitted as poor; the active two-phase Sn–Ni/Ni3Sn2 material remains unrepresentable. |
| [Chen et al. (2023)](https://doi.org/10.1126/science.adh8872) | Hold | The active Ni–Mo–Bi ternary liquid cannot be represented by the binary molten genome. |
| [Xu et al. (2021)](https://doi.org/10.1016/j.apcata.2020.117967) | Hold | Activity depends on measured Ni particle size and support; the current solid genome does not encode particle size or loading. |
| [Li and Smith (2008)](https://doi.org/10.1016/j.apcata.2008.07.011) | Hold | Supported Co and Ni identities lack the unique facet, loading, and particle state required for an exact screening genome. |
| [Chein and Ou (2025)](https://doi.org/10.1016/j.ijhydene.2024.07.134) | Hold | Nickel-foam PPI controls stability, but foam morphology/PPI is absent from the genome. |
| [Ding et al. (2020)](https://doi.org/10.1016/j.crcon.2020.12.002) | Hold | Activated-carbon-supported alloy loading and morphology are not exactly represented by the current solid genome. |
| [Kim et al. (2023)](https://doi.org/10.1016/j.carbon.2023.02.053) | Hold | The 15% to 97% change is caused by gas-sparger geometry at fixed Ni–Sn composition; the catalyst descriptor cannot predict that reactor variable. |
| [Scheiblehner et al. (2023)](https://doi.org/10.1016/j.ijhydene.2022.08.115) | Hold | The binary Cu-alloy experiments operate at 1160 °C (1433 K), outside the 1300 K screening band. |
| [Scheiblehner et al. (2023)](https://doi.org/10.3390/met13071310) | Hold | The Cu-Ni-Sn materials are ternary and were tested at 1160 °C, so both composition and temperature are out of scope. |
| [Wei et al. (2026)](https://doi.org/10.1016/j.cej.2026.177959) | Hold | Nickel metal fog in molten carbonate is not represented by a binary molten-metal genome. |
| [Kerosirov et al. (2023)](https://doi.org/10.1051/e3sconf/202341101009) | Hold | Pure Sn is exact, but the reported product concentrations do not establish an externally matched poor-catalyst label. |
| [Reddy et al. (2024)](https://doi.org/10.1007/s42823-024-00716-2) | Hold | Pure Sn is exact, but the paper characterizes carbon products rather than establishing catalyst activity against a matched comparator. |

## Stopping rule

Do not loosen exact encoding, count multiple catalyst families from one paper
as independent publications, infer unreported facets, or label a catalyst from
absolute conversion across unmatched reactors. Curation resumes when primary
data add exact in-scope contrasts, or after a separately reviewed design-space
extension represents the missing physical variable.
