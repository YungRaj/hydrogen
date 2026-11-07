"""Multi-hydrocarbon and biomass feedstock thermochemistry, stoichiometry, and BEP kinetics registry.

Methane (`CH4`) is the strict default across all reactor models and TEA calculations.
Additional hydrocarbon, waste polyolefin, and biomass/biogenic feedstocks (`ASSOC_GAS`,
`C2H6`, `C3H8`, `C4H10`, `C5H12`, `C6H14`, `C6H12`, `C8H18`, `C12H26` / `diesel`,
`C16H34` / `cetane`, `C10H20_PE` / `waste_plastic`, `BIOMETHANE`, `BIOGAS`,
`C18H34O2_LIPID` / `waste_cooking_oil`, `C3H8O3_GLYCEROL`, `C6H10O2_BIOOIL`,
and `C6H10O4_BIOMASS` / `biomass`) enable comparative simulation of natural-gas,
light-alkane, naphtha, distillate, plastic-wax, and biomass pyrolysis:

    Pure hydrocarbons (n_O == 0): C_n H_m -> n C(s, graphite) + (m / 2) H2(g)
    Oxygenated biomass (n_O > 0): C_n H_m O_k -> (n - k) C(s) + k CO(g) + (m / 2) H2(g)

References:
- NIST Chemistry WebBook (2024): Standard gas-phase enthalpies and entropies of formation.
- Upham et al., Science 358, 917-921 (2017): Catalytic molten metals for the direct
  conversion of methane, ethane, propane, hexane, and cetane/diesel to hydrogen and
  separable carbon.
- Palmer et al., Energy & Fuels 35, 14065-14074 (2021): Pyrolysis of hydrocarbon
  feedstocks in molten metal bubble column reactors.
- Sinfelt, Catal. Rev. 3, 175-205 (1969) & Nørskov et al., J. Catal. 209, 275 (2002):
  Brønsted-Evans-Polanyi (BEP) scaling of C-C and secondary C-H activation barriers
  relative to methane primary C-H activation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import cast

import numpy as np

from pipeline.data_models.reactors import (
    FeedstockCategory,
    FeedstockName,
    ReactorResult,
)
from pipeline.utils import R_J_MOL_K

DEFAULT_FEEDSTOCK: FeedstockName = 'CH4'

MW_C_G_MOL: float = 12.011
MW_H_G_MOL: float = 1.00794
MW_O_G_MOL: float = 15.999
MW_H2_G_MOL: float = 2.01588


@dataclass(frozen=True)
class HydrocarbonFeedstockSpec:
    """Thermochemical, stoichiometric, kinetic, and economic specification for a hydrocarbon or biogenic feedstock.

    Supports both pure hydrocarbons (`n_O == 0`: `C_n H_m -> n C(s) + (m/2) H2`) and
    oxygen-bearing biomass/biogenic feedstocks (`n_O > 0`:
    `C_n H_m O_k -> (n - k) C(s) + k CO(g) + (m/2) H2(g)`).

    Attributes:
        formula: Canonical feedstock key (`CH4` default, `ASSOC_GAS`, `C2H6`, `C3H8`, `C4H10`,
            `C5H12`, `C6H14`, `C6H12`, `C8H18`, `C12H26`, `C16H34`, `C10H20_PE`, `BIOMETHANE`,
            `BIOGAS`, `C18H34O2_LIPID`, `C3H8O3_GLYCEROL`, `C6H10O2_BIOOIL`, `C6H10O4_BIOMASS`).
        name: Human-readable feedstock name.
        n_C: Number of carbon atoms per formula unit (`n`).
        n_H: Number of hydrogen atoms per formula unit (`m`).
        dh_rxn_298_kJ_mol_feed: Standard enthalpy of complete pyrolysis to graphite (+ CO) + H2 at 298.15 K [kJ/mol feed].
        ds_rxn_298_J_mol_K_feed: Standard entropy of complete pyrolysis at 298.15 K [J/(mol·K)].
        dcp_rxn_J_mol_K_feed: Mean heat-capacity change ΔCp of pyrolysis across 300-1300 K [J/(mol·K)].
        bep_e_act_offset_eV: Brønsted-Evans-Polanyi C-C / secondary C-H / C-O activation barrier offset
            relative to CH4 primary C-H scission [eV] (negative means lower barrier than CH4).
        hom_a_1_s: Homogeneous gas-phase thermal cracking prefactor [1/s].
        hom_ea_kJ_mol: Homogeneous gas-phase thermal cracking activation energy [kJ/mol].
        hhv_MJ_kg: Higher heating value of the feedstock [MJ/kg].
        base_price_usd_mmbtu_multiplier: Industrial bulk feedstock price multiplier relative to pipeline natural gas ($/MMBtu basis).
        literature_reference: Primary literature reference for pyrolysis of this feedstock.
        n_O: Number of oxygen atoms per formula unit (`k`, default `0` for pure hydrocarbons).
        category: Feedstock family (`natural_gas`, `light_alkane`, `liquid_hydrocarbon`,
            `cycloalkane_olefin`, `waste_plastic`, `biomass_biogenic`).
    """

    formula: FeedstockName
    name: str
    n_C: int
    n_H: int
    dh_rxn_298_kJ_mol_feed: float
    ds_rxn_298_J_mol_K_feed: float
    dcp_rxn_J_mol_K_feed: float
    bep_e_act_offset_eV: float
    hom_a_1_s: float
    hom_ea_kJ_mol: float
    hhv_MJ_kg: float
    base_price_usd_mmbtu_multiplier: float
    literature_reference: str
    n_O: int = 0
    category: FeedstockCategory = 'light_alkane'

    @property
    def h_to_c_ratio(self) -> float:
        """Atomic H/C ratio of the feedstock (m / n).

        Returns:
            Atomic hydrogen-to-carbon ratio.
        """
        return float(self.n_H) / max(float(self.n_C), 1.0)

    @property
    def nu_H2(self) -> float:
        """Moles of H2 produced per mole of feedstock completely cracked (m / 2).

        Returns:
            Stoichiometric moles of H2 per mole of feedstock.
        """
        return 0.5 * float(self.n_H)

    @property
    def nu_CO(self) -> float:
        """Moles of CO(g) produced per mole of oxygenated feedstock cracked (`0.0` for pure hydrocarbons).

        Returns:
            Stoichiometric moles of CO per mole of feedstock.
        """
        return float(self.n_O)

    @property
    def nu_C(self) -> float:
        """Moles of solid carbon C(s) produced per mole of feedstock completely cracked.

        Returns:
            Stoichiometric moles of solid carbon per mole of feedstock.
        """
        if self.n_O <= 0:
            return float(self.n_C)
        return max(float(self.n_C - self.n_O), 0.25 * float(self.n_C))

    @property
    def gas_expansion_delta_nu(self) -> float:
        """Net increase in gas moles per mole of feedstock cracked (nu_H2 + nu_CO - 1).

        Returns:
            Net change in gas-phase moles per mole of feedstock converted.
        """
        return self.nu_H2 + self.nu_CO - 1.0

    @property
    def molar_mass_g_mol(self) -> float:
        """Molecular weight in g/mol.

        Returns:
            Molar mass of the feedstock formula unit in g/mol.
        """
        return (
            float(self.n_C) * MW_C_G_MOL
            + float(self.n_H) * MW_H_G_MOL
            + float(self.n_O) * MW_O_G_MOL
        )

    @property
    def kg_feed_per_kg_h2(self) -> float:
        """Stoichiometric mass of feedstock required per kg of H2 produced [kg feed / kg H2].

        Returns:
            Mass of feedstock per unit mass of H2 produced in kg/kg.
        """
        return self.molar_mass_g_mol / (self.nu_H2 * MW_H2_G_MOL)

    @property
    def kg_c_per_kg_h2(self) -> float:
        """Stoichiometric mass of solid carbon co-produced per kg of H2 [kg C / kg H2].

        Returns:
            Mass of solid carbon co-product per unit mass of H2 in kg/kg.
        """
        return (self.nu_C * MW_C_G_MOL) / (self.nu_H2 * MW_H2_G_MOL)

    @property
    def dh_rxn_298_kJ_mol_h2(self) -> float:
        """Standard pyrolysis enthalpy at 298.15 K normalized per mole of H2 [kJ/mol H2].

        Returns:
            Standard reaction enthalpy at 298.15 K in kJ per mole of H2.
        """
        return self.dh_rxn_298_kJ_mol_feed / self.nu_H2


HYDROCARBON_FEEDSTOCKS: dict[FeedstockName, HydrocarbonFeedstockSpec] = {
    # ------------------------------------------------------------------
    # 1. Natural Gas & Associated / Flare Gas Family (CH4 is DEFAULT)
    # ------------------------------------------------------------------
    'CH4': HydrocarbonFeedstockSpec(
        formula='CH4',
        name='Methane (Pipeline Natural Gas)',
        n_C=1,
        n_H=4,
        dh_rxn_298_kJ_mol_feed=74.85,
        ds_rxn_298_J_mol_K_feed=80.80,
        dcp_rxn_J_mol_K_feed=15.5,
        bep_e_act_offset_eV=0.00,
        hom_a_1_s=4.5e13,
        hom_ea_kJ_mol=356.0,
        hhv_MJ_kg=55.50,
        base_price_usd_mmbtu_multiplier=1.00,
        literature_reference='Upham et al., Science 2017; Palmer et al., ACS Catal. 2020',
        n_O=0,
        category='natural_gas',
    ),
    'ASSOC_GAS': HydrocarbonFeedstockSpec(
        formula='ASSOC_GAS',
        name='Associated / Flare Gas (Wet C1-C4 Blend)',
        n_C=2,
        n_H=7,
        dh_rxn_298_kJ_mol_feed=81.50,
        ds_rxn_298_J_mol_K_feed=148.00,
        dcp_rxn_J_mol_K_feed=17.5,
        bep_e_act_offset_eV=-0.10,
        hom_a_1_s=8.5e13,
        hom_ea_kJ_mol=340.0,
        hhv_MJ_kg=53.10,
        base_price_usd_mmbtu_multiplier=0.35,
        literature_reference='Palmer et al., Energy Fuels 2021 (Associated Gas & C1-C3 Co-Pyrolysis)',
        n_O=0,
        category='natural_gas',
    ),
    # ------------------------------------------------------------------
    # 2. Light Alkanes, NGLs, and Wholesale LPG / Refinery Raffinates
    # ------------------------------------------------------------------
    'C2H6': HydrocarbonFeedstockSpec(
        formula='C2H6',
        name='Ethane (Shale NGL)',
        n_C=2,
        n_H=6,
        dh_rxn_298_kJ_mol_feed=83.82,
        ds_rxn_298_J_mol_K_feed=174.10,
        dcp_rxn_J_mol_K_feed=18.2,
        bep_e_act_offset_eV=-0.14,
        hom_a_1_s=1.1e14,
        hom_ea_kJ_mol=336.0,
        hhv_MJ_kg=51.90,
        base_price_usd_mmbtu_multiplier=0.95,
        literature_reference='Upham et al., Science 2017; Sinfelt, Catal. Rev. 1969',
        n_O=0,
        category='light_alkane',
    ),
    'C3H8': HydrocarbonFeedstockSpec(
        formula='C3H8',
        name='Propane (Wholesale LPG)',
        n_C=3,
        n_H=8,
        dh_rxn_298_kJ_mol_feed=103.85,
        ds_rxn_298_J_mol_K_feed=269.50,
        dcp_rxn_J_mol_K_feed=21.0,
        bep_e_act_offset_eV=-0.20,
        hom_a_1_s=2.0e14,
        hom_ea_kJ_mol=318.0,
        hhv_MJ_kg=50.35,
        base_price_usd_mmbtu_multiplier=1.25,
        literature_reference='Takenaka et al., J. Catal. 2004; Palmer et al., Energy Fuels 2021',
        n_O=0,
        category='light_alkane',
    ),
    'C4H10': HydrocarbonFeedstockSpec(
        formula='C4H10',
        name='n-Butane (Refinery C4 Raffinate / LPG)',
        n_C=4,
        n_H=10,
        dh_rxn_298_kJ_mol_feed=125.79,
        ds_rxn_298_J_mol_K_feed=365.90,
        dcp_rxn_J_mol_K_feed=24.0,
        bep_e_act_offset_eV=-0.23,
        hom_a_1_s=2.5e14,
        hom_ea_kJ_mol=310.0,
        hhv_MJ_kg=49.50,
        base_price_usd_mmbtu_multiplier=1.20,
        literature_reference='Palmer et al., Energy Fuels 2021',
        n_O=0,
        category='light_alkane',
    ),
    # ------------------------------------------------------------------
    # 3. Practical Liquid Hydrocarbons (Condensates, Naphthas, Cycloalkanes, Distillates)
    # ------------------------------------------------------------------
    'C5H12': HydrocarbonFeedstockSpec(
        formula='C5H12',
        name='n-Pentane (Natural Gasoline / C5+ Condensate)',
        n_C=5,
        n_H=12,
        dh_rxn_298_kJ_mol_feed=146.76,
        ds_rxn_298_J_mol_K_feed=462.10,
        dcp_rxn_J_mol_K_feed=27.2,
        bep_e_act_offset_eV=-0.25,
        hom_a_1_s=3.0e14,
        hom_ea_kJ_mol=304.0,
        hhv_MJ_kg=48.80,
        base_price_usd_mmbtu_multiplier=1.25,
        literature_reference='NIST WebBook 2024; Gas-Plant Condensate Catalytic Cracking',
        n_O=0,
        category='liquid_hydrocarbon',
    ),
    'C6H14': HydrocarbonFeedstockSpec(
        formula='C6H14',
        name='n-Hexane (Straight-Run Light Naphtha)',
        n_C=6,
        n_H=14,
        dh_rxn_298_kJ_mol_feed=167.19,
        ds_rxn_298_J_mol_K_feed=558.40,
        dcp_rxn_J_mol_K_feed=30.5,
        bep_e_act_offset_eV=-0.26,
        hom_a_1_s=3.5e14,
        hom_ea_kJ_mol=298.0,
        hhv_MJ_kg=48.30,
        base_price_usd_mmbtu_multiplier=1.35,
        literature_reference='Upham et al., Science 2017 (Liquid Hydrocarbon Melt Pyrolysis)',
        n_O=0,
        category='liquid_hydrocarbon',
    ),
    'C6H12': HydrocarbonFeedstockSpec(
        formula='C6H12',
        name='Cyclohexane (Naphthenic Ring Cut / LOHC)',
        n_C=6,
        n_H=12,
        dh_rxn_298_kJ_mol_feed=123.40,
        ds_rxn_298_J_mol_K_feed=518.50,
        dcp_rxn_J_mol_K_feed=29.0,
        bep_e_act_offset_eV=-0.28,
        hom_a_1_s=3.2e14,
        hom_ea_kJ_mol=294.0,
        hhv_MJ_kg=46.60,
        base_price_usd_mmbtu_multiplier=1.38,
        literature_reference='Sinfelt, Catal. Rev. 1969; Cycloalkane Dehydrogenation & Pyrolysis',
        n_O=0,
        category='cycloalkane_olefin',
    ),
    'C8H18': HydrocarbonFeedstockSpec(
        formula='C8H18',
        name='n-Octane (Heavy Naphtha / Gasoline Surrogate)',
        n_C=8,
        n_H=18,
        dh_rxn_298_kJ_mol_feed=208.45,
        ds_rxn_298_J_mol_K_feed=751.20,
        dcp_rxn_J_mol_K_feed=36.0,
        bep_e_act_offset_eV=-0.29,
        hom_a_1_s=4.5e14,
        hom_ea_kJ_mol=290.0,
        hhv_MJ_kg=47.90,
        base_price_usd_mmbtu_multiplier=1.45,
        literature_reference='Palmer et al., Energy Fuels 2021; Logistics Fuel Pyrolysis',
        n_O=0,
        category='liquid_hydrocarbon',
    ),
    'C12H26': HydrocarbonFeedstockSpec(
        formula='C12H26',
        name='n-Dodecane (Kerosene / Middle Distillate Diesel)',
        n_C=12,
        n_H=26,
        dh_rxn_298_kJ_mol_feed=290.90,
        ds_rxn_298_J_mol_K_feed=1136.50,
        dcp_rxn_J_mol_K_feed=48.0,
        bep_e_act_offset_eV=-0.32,
        hom_a_1_s=6.0e14,
        hom_ea_kJ_mol=280.0,
        hhv_MJ_kg=47.20,
        base_price_usd_mmbtu_multiplier=1.55,
        literature_reference='Upham et al., Science 2017 (Diesel/Distillate Molten Metal Pyrolysis)',
        n_O=0,
        category='liquid_hydrocarbon',
    ),
    'C16H34': HydrocarbonFeedstockSpec(
        formula='C16H34',
        name='n-Hexadecane / Cetane (Heavy Gas Oil / Slack Wax)',
        n_C=16,
        n_H=34,
        dh_rxn_298_kJ_mol_feed=373.30,
        ds_rxn_298_J_mol_K_feed=1521.80,
        dcp_rxn_J_mol_K_feed=60.0,
        bep_e_act_offset_eV=-0.34,
        hom_a_1_s=7.5e14,
        hom_ea_kJ_mol=274.0,
        hhv_MJ_kg=46.80,
        base_price_usd_mmbtu_multiplier=1.30,
        literature_reference='Upham et al., Science 2017 (Cetane C16H34 Molten Bi-Ni Pyrolysis)',
        n_O=0,
        category='liquid_hydrocarbon',
    ),
    # ------------------------------------------------------------------
    # 4. Waste Polyolefin Plastics (PE/PP Pyrolysis Wax-Oil)
    # ------------------------------------------------------------------
    'C10H20_PE': HydrocarbonFeedstockSpec(
        formula='C10H20_PE',
        name='Waste Polyolefin Pyrolysis Oil (PE/PP -(CH2)-n)',
        n_C=10,
        n_H=20,
        dh_rxn_298_kJ_mol_feed=165.00,
        ds_rxn_298_J_mol_K_feed=890.00,
        dcp_rxn_J_mol_K_feed=41.0,
        bep_e_act_offset_eV=-0.31,
        hom_a_1_s=5.2e14,
        hom_ea_kJ_mol=282.0,
        hhv_MJ_kg=46.20,
        base_price_usd_mmbtu_multiplier=0.45,
        literature_reference='Jie et al., Nature Catal. 2020; Williams et al. (Waste Polyolefin to H2 + CNTs)',
        n_O=0,
        category='waste_plastic',
    ),
    # ------------------------------------------------------------------
    # 5. Biomass & Biogenic Feedstocks (Carbon-Negative / Renewable H2)
    # ------------------------------------------------------------------
    'BIOMETHANE': HydrocarbonFeedstockSpec(
        formula='BIOMETHANE',
        name='Biomethane / RNG (Upgraded Anaerobic Biogas)',
        n_C=1,
        n_H=4,
        dh_rxn_298_kJ_mol_feed=74.85,
        ds_rxn_298_J_mol_K_feed=80.80,
        dcp_rxn_J_mol_K_feed=15.5,
        bep_e_act_offset_eV=0.00,
        hom_a_1_s=4.5e13,
        hom_ea_kJ_mol=356.0,
        hhv_MJ_kg=55.50,
        base_price_usd_mmbtu_multiplier=0.85,
        literature_reference='timmerberg et al., Energy Convers. Manag. 2020 (Carbon-Negative Biomethane Pyrolysis)',
        n_O=0,
        category='biomass_biogenic',
    ),
    'BIOGAS': HydrocarbonFeedstockSpec(
        formula='BIOGAS',
        name='Raw Anaerobic Digester / Landfill Biogas (CH4 + CO2)',
        n_C=2,
        n_H=6,
        dh_rxn_298_kJ_mol_feed=136.00,
        ds_rxn_298_J_mol_K_feed=212.00,
        dcp_rxn_J_mol_K_feed=20.0,
        bep_e_act_offset_eV=-0.06,
        hom_a_1_s=6.0e13,
        hom_ea_kJ_mol=345.0,
        hhv_MJ_kg=33.50,
        base_price_usd_mmbtu_multiplier=0.45,
        literature_reference='Palmer et al., Renew. Sustain. Energy Rev. 2021 (Direct Biogas Pyrolysis)',
        n_O=1,
        category='biomass_biogenic',
    ),
    'C18H34O2_LIPID': HydrocarbonFeedstockSpec(
        formula='C18H34O2_LIPID',
        name='Waste Cooking Oil / Tall Oil Bio-Lipid (C18H34O2)',
        n_C=18,
        n_H=34,
        dh_rxn_298_kJ_mol_feed=348.00,
        ds_rxn_298_J_mol_K_feed=1680.00,
        dcp_rxn_J_mol_K_feed=64.0,
        bep_e_act_offset_eV=-0.33,
        hom_a_1_s=7.0e14,
        hom_ea_kJ_mol=275.0,
        hhv_MJ_kg=39.60,
        base_price_usd_mmbtu_multiplier=0.80,
        literature_reference='Czernik & French, Energy Fuels 2006 (Catalytic Pyrolysis of Waste Lipids to H2)',
        n_O=2,
        category='biomass_biogenic',
    ),
    'C3H8O3_GLYCEROL': HydrocarbonFeedstockSpec(
        formula='C3H8O3_GLYCEROL',
        name='Crude Biodiesel Glycerol Byproduct (C3H8O3)',
        n_C=3,
        n_H=8,
        dh_rxn_298_kJ_mol_feed=138.00,
        ds_rxn_298_J_mol_K_feed=410.00,
        dcp_rxn_J_mol_K_feed=28.0,
        bep_e_act_offset_eV=-0.27,
        hom_a_1_s=3.8e14,
        hom_ea_kJ_mol=288.0,
        hhv_MJ_kg=19.20,
        base_price_usd_mmbtu_multiplier=0.50,
        literature_reference='Valliyappan et al., Int. J. Hydrogen Energy 2008 (Glycerol Pyrolysis to H2 + Carbon)',
        n_O=2,
        category='biomass_biogenic',
    ),
    'C6H10O2_BIOOIL': HydrocarbonFeedstockSpec(
        formula='C6H10O2_BIOOIL',
        name='Fast-Pyrolysis Bio-Oil (Hydro-Deoxygenated C6H10O2)',
        n_C=6,
        n_H=10,
        dh_rxn_298_kJ_mol_feed=122.00,
        ds_rxn_298_J_mol_K_feed=540.00,
        dcp_rxn_J_mol_K_feed=32.0,
        bep_e_act_offset_eV=-0.26,
        hom_a_1_s=3.4e14,
        hom_ea_kJ_mol=292.0,
        hhv_MJ_kg=26.50,
        base_price_usd_mmbtu_multiplier=0.65,
        literature_reference='Wang et al., Ind. Eng. Chem. Res. 1997 (Catalytic Cracking of Bio-Oil to H2 + Biochar)',
        n_O=2,
        category='biomass_biogenic',
    ),
    'C6H10O4_BIOMASS': HydrocarbonFeedstockSpec(
        formula='C6H10O4_BIOMASS',
        name='Torrefied Lignocellulosic Biomass (Woody/Ag C6H10O4)',
        n_C=6,
        n_H=10,
        dh_rxn_298_kJ_mol_feed=118.00,
        ds_rxn_298_J_mol_K_feed=590.00,
        dcp_rxn_J_mol_K_feed=34.0,
        bep_e_act_offset_eV=-0.24,
        hom_a_1_s=2.8e14,
        hom_ea_kJ_mol=296.0,
        hhv_MJ_kg=20.80,
        base_price_usd_mmbtu_multiplier=0.50,
        literature_reference='Chen et al., Prog. Energy Combust. Sci. 2021 (Biomass Catalytic Pyrolysis for H2 + Biochar)',
        n_O=4,
        category='biomass_biogenic',
    ),
}

_FEEDSTOCK_ALIASES: dict[str, FeedstockName] = {
    'ch4': 'CH4',
    'methane': 'CH4',
    'natural_gas': 'CH4',
    'assoc_gas': 'ASSOC_GAS',
    'associated_gas': 'ASSOC_GAS',
    'flare_gas': 'ASSOC_GAS',
    'wet_gas': 'ASSOC_GAS',
    'c2h6': 'C2H6',
    'ethane': 'C2H6',
    'ngl': 'C2H6',
    'c3h8': 'C3H8',
    'propane': 'C3H8',
    'lpg': 'C3H8',
    'c4h10': 'C4H10',
    'butane': 'C4H10',
    'n-butane': 'C4H10',
    'raffinate': 'C4H10',
    'c5h12': 'C5H12',
    'pentane': 'C5H12',
    'n-pentane': 'C5H12',
    'condensate': 'C5H12',
    'natural_gasoline': 'C5H12',
    'c6h14': 'C6H14',
    'hexane': 'C6H14',
    'n-hexane': 'C6H14',
    'naphtha': 'C6H14',
    'light_naphtha': 'C6H14',
    'c6h12': 'C6H12',
    'cyclohexane': 'C6H12',
    'naphthene': 'C6H12',
    'c8h18': 'C8H18',
    'octane': 'C8H18',
    'n-octane': 'C8H18',
    'gasoline': 'C8H18',
    'heavy_naphtha': 'C8H18',
    'c12h26': 'C12H26',
    'dodecane': 'C12H26',
    'n-dodecane': 'C12H26',
    'diesel': 'C12H26',
    'kerosene': 'C12H26',
    'distillate': 'C12H26',
    'c16h34': 'C16H34',
    'hexadecane': 'C16H34',
    'cetane': 'C16H34',
    'heavy_gas_oil': 'C16H34',
    'slack_wax': 'C16H34',
    'c10h20_pe': 'C10H20_PE',
    'waste_plastic': 'C10H20_PE',
    'polyolefin': 'C10H20_PE',
    'pe_wax': 'C10H20_PE',
    'plastic_oil': 'C10H20_PE',
    'biomethane': 'BIOMETHANE',
    'rng': 'BIOMETHANE',
    'renewable_natural_gas': 'BIOMETHANE',
    'biogas': 'BIOGAS',
    'landfill_gas': 'BIOGAS',
    'digester_gas': 'BIOGAS',
    'c18h34o2_lipid': 'C18H34O2_LIPID',
    'waste_cooking_oil': 'C18H34O2_LIPID',
    'yellow_grease': 'C18H34O2_LIPID',
    'tall_oil': 'C18H34O2_LIPID',
    'bio_lipid': 'C18H34O2_LIPID',
    'c3h8o3_glycerol': 'C3H8O3_GLYCEROL',
    'glycerol': 'C3H8O3_GLYCEROL',
    'crude_glycerol': 'C3H8O3_GLYCEROL',
    'c6h10o2_biooil': 'C6H10O2_BIOOIL',
    'bio_oil': 'C6H10O2_BIOOIL',
    'biooil': 'C6H10O2_BIOOIL',
    'pyrolysis_bio_oil': 'C6H10O2_BIOOIL',
    'c6h10o4_biomass': 'C6H10O4_BIOMASS',
    'biomass': 'C6H10O4_BIOMASS',
    'lignocellulose': 'C6H10O4_BIOMASS',
    'wood_biomass': 'C6H10O4_BIOMASS',
    'torrefied_biomass': 'C6H10O4_BIOMASS',
}


def resolve_feedstock(
    feedstock: FeedstockName | str | None = None,
) -> HydrocarbonFeedstockSpec:
    """Resolve a feedstock formula or alias into its `HydrocarbonFeedstockSpec`.

    Defaults strictly to `CH4` when `feedstock` is `None` or empty.

    Args:
        feedstock: Feedstock formula or common name (e.g., `'CH4'`, `'propane'`, `'diesel'`).

    Returns:
        Matching `HydrocarbonFeedstockSpec`.

    Raises:
        ValueError: If `feedstock` is not recognized.
    """
    if feedstock is None or not str(feedstock).strip():
        return HYDROCARBON_FEEDSTOCKS[DEFAULT_FEEDSTOCK]
    key: str = str(feedstock).strip()
    if key in HYDROCARBON_FEEDSTOCKS:
        return HYDROCARBON_FEEDSTOCKS[cast(FeedstockName, key)]
    canonical: FeedstockName | None = _FEEDSTOCK_ALIASES.get(key.lower())
    if canonical is not None:
        return HYDROCARBON_FEEDSTOCKS[canonical]
    raise ValueError(
        f'Unsupported hydrocarbon feedstock {feedstock!r}; supported feedstocks: '
        f'{sorted(HYDROCARBON_FEEDSTOCKS.keys())} (or aliases: {sorted(_FEEDSTOCK_ALIASES.keys())})'
    )


def feedstock_reaction_enthalpy_kJ_mol(
    T_K: float,
    feedstock: FeedstockName | str | HydrocarbonFeedstockSpec = DEFAULT_FEEDSTOCK,
) -> float:
    """Compute endothermic reaction enthalpy ΔH_rxn(T) [kJ/mol feed] at temperature `T_K`.

    Args:
        T_K: Temperature in kelvin.
        feedstock: Feedstock identifier or `HydrocarbonFeedstockSpec`.

    Returns:
        Reaction enthalpy in kJ per mole of feedstock converted.
    """
    spec: HydrocarbonFeedstockSpec = (
        feedstock
        if isinstance(feedstock, HydrocarbonFeedstockSpec)
        else resolve_feedstock(feedstock)
    )
    dt: float = max(float(T_K), 298.15) - 298.15
    return float(spec.dh_rxn_298_kJ_mol_feed + 1e-3 * spec.dcp_rxn_J_mol_K_feed * dt)


def feedstock_equilibrium_conversion(
    T_K: float,
    P_Pa: float = 101325.0,
    feedstock: FeedstockName | str | HydrocarbonFeedstockSpec = DEFAULT_FEEDSTOCK,
    y_feed_inlet: float = 0.95,
) -> float:
    """Compute thermodynamic equilibrium conversion X_eq(T, P) for C_n H_m -> n C(s) + (m/2) H2(g).

    For `CH4`, callers continue to use `tabulated_x_eq(T_K)` or `local_tp_x_eq(T_K, P_Pa)`.
    For higher hydrocarbons (`n >= 2`), solves the exact ideal-gas reaction equilibrium
    with activity of condensed graphite `a_C(s) = 1`.

    Args:
        T_K: Reactor temperature in kelvin.
        P_Pa: Reactor pressure in pascals.
        feedstock: Feedstock identifier or `HydrocarbonFeedstockSpec`.
        y_feed_inlet: Inlet hydrocarbon mole fraction (balance inert Ar).

    Returns:
        Equilibrium conversion fraction in `[0.0, 0.999]`.
    """
    spec: HydrocarbonFeedstockSpec = (
        feedstock
        if isinstance(feedstock, HydrocarbonFeedstockSpec)
        else resolve_feedstock(feedstock)
    )
    t_eff: float = max(float(T_K), 300.0)
    dt: float = t_eff - 298.15
    dh_J: float = spec.dh_rxn_298_kJ_mol_feed * 1000.0 + spec.dcp_rxn_J_mol_K_feed * dt
    ds_J_K: float = spec.ds_rxn_298_J_mol_K_feed + spec.dcp_rxn_J_mol_K_feed * math.log(
        t_eff / 298.15
    )
    dg_J: float = dh_J - t_eff * ds_J_K
    ln_k_eq: float = -dg_J / (R_J_MOL_K * t_eff)
    ln_k_eq = min(max(ln_k_eq, -60.0), 60.0)
    k_eq: float = math.exp(ln_k_eq)

    p_bar: float = max(float(P_Pa) / 1e5, 1e-3)
    y0: float = min(max(float(y_feed_inlet), 1e-4), 1.0)
    nu_gas: float = spec.nu_H2 + spec.nu_CO

    # Bisection on X in (0, 0.999) for:
    # K_eq = ((nu_gas * y0 * X)^nu_gas / ((y0 * (1 - X)) * (1 + (nu_gas - 1) * y0 * X)^(nu_gas - 1))) * p_bar^(nu_gas - 1)
    lo: float = 1e-9
    hi: float = 0.999
    for _ in range(48):
        mid: float = 0.5 * (lo + hi)
        n_tot: float = 1.0 + (nu_gas - 1.0) * y0 * mid
        lhs: float = (
            ((nu_gas * y0 * mid) ** nu_gas)
            * (p_bar ** (nu_gas - 1.0))
            / max((y0 * (1.0 - mid)) * (n_tot ** (nu_gas - 1.0)), 1e-30)
        )
        if lhs < k_eq:
            lo = mid
        else:
            hi = mid
    return float(np.clip(0.5 * (lo + hi), 0.0, 0.999))


def effective_catalyst_barrier_for_feedstock(
    catalyst_E_act_CH4_eV: float,
    feedstock: FeedstockName | str | HydrocarbonFeedstockSpec = DEFAULT_FEEDSTOCK,
) -> float:
    """Apply the BEP C-C / secondary C-H / C-O barrier offset relative to methane activation.

    For `CH4` (`bep_e_act_offset_eV == 0.0`), returns `catalyst_E_act_CH4_eV` unchanged.

    Args:
        catalyst_E_act_CH4_eV: Methane C-H activation barrier on the catalyst [eV].
        feedstock: Feedstock identifier or `HydrocarbonFeedstockSpec`.

    Returns:
        Effective rate-limiting activation barrier for the feedstock [eV].
    """
    spec: HydrocarbonFeedstockSpec = (
        feedstock
        if isinstance(feedstock, HydrocarbonFeedstockSpec)
        else resolve_feedstock(feedstock)
    )
    if spec.formula == DEFAULT_FEEDSTOCK or spec.bep_e_act_offset_eV == 0.0:
        return float(catalyst_E_act_CH4_eV)
    return float(max(0.18, float(catalyst_E_act_CH4_eV) + spec.bep_e_act_offset_eV))


def adapt_reactor_result_to_feedstock(
    base_ch4_result: ReactorResult,
    feedstock: FeedstockName | str | HydrocarbonFeedstockSpec = DEFAULT_FEEDSTOCK,
    *,
    y_feed_inlet: float = 0.95,
    y_ar_inlet: float = 0.05,
) -> ReactorResult:
    """Adapt a baseline CH4 reactor simulation result to the selected hydrocarbon or biomass feedstock.

    Architecture note:
        Cantera gas/surface mechanisms in the screening sweep integrate the
        baseline methane (``CH4``) reaction network first (`base_ch4_result`).
        When `feedstock == 'CH4'` (the strict default), this function returns a
        shallow copy of `base_ch4_result` enriched with feedstock-general
        metadata and explicit ``feedstock_*`` keys so all existing methane
        fields remain 100% identical. When a non-methane alkane, cycloalkane,
        polyolefin, or biomass feedstock is selected, this function maps the
        baseline CH4 reference metrics (`*_ch4_ref`) to the active feedstock
        (`*_feed`) via exact stoichiometry, Kirchhoff ΔG°(T, P) equilibrium,
        and BEP C-C / secondary C-H / C-O barrier scaling.

    Args:
        base_ch4_result: Baseline methane reactor simulation dictionary produced by `simulate_*`.
        feedstock: Feedstock formula or alias (`'CH4'`, `'ASSOC_GAS'`, `'C3H8'`, `'diesel'`, `'biomass'`, etc.).
        y_feed_inlet: Inlet hydrocarbon/feedstock mole fraction.
        y_ar_inlet: Inlet argon tracer mole fraction.

    Returns:
        Reactor result dictionary adapted to `feedstock`, carrying both explicit
        `feedstock_*` fields and synchronized legacy `CH4_*` keys.
    """
    spec: HydrocarbonFeedstockSpec = (
        feedstock
        if isinstance(feedstock, HydrocarbonFeedstockSpec)
        else resolve_feedstock(feedstock)
    )
    out: ReactorResult = cast(ReactorResult, dict(base_ch4_result))
    out['feedstock'] = spec.formula
    out['feedstock_name'] = spec.name
    out['feedstock_category'] = spec.category
    out['feedstock_n_C'] = spec.n_C
    out['feedstock_n_H'] = spec.n_H
    out['feedstock_n_O'] = spec.n_O
    out['feedstock_h_to_c_ratio'] = round(spec.h_to_c_ratio, 4)
    out['feedstock_nu_H2'] = spec.nu_H2
    out['feedstock_nu_C'] = round(spec.nu_C, 4)
    out['feedstock_nu_CO'] = round(spec.nu_CO, 4)
    out['feedstock_kg_per_kg_h2'] = round(spec.kg_feed_per_kg_h2, 6)
    out['feedstock_kg_c_per_kg_h2'] = round(spec.kg_c_per_kg_h2, 6)

    if (
        spec.formula == DEFAULT_FEEDSTOCK
        or (spec.n_C == 1 and spec.n_H == 4 and spec.n_O == 0)
        or out.get('status') not in (None, 'complete')
    ):
        has_ch4_conv: bool = (
            'CH4_conversion' in out or 'single_pass_CH4_conversion' in out
        )
        default_conv: float = float(
            (
                out.get('CH4_conversion')
                if out.get('CH4_conversion') is not None
                else out.get('single_pass_CH4_conversion')
            )
            if has_ch4_conv
            else (
                out.get('feedstock_conversion')
                or out.get('single_pass_feedstock_conversion')
                or 0.0
            )
        )
        default_single_pass_conv: float = float(
            (
                out.get('single_pass_CH4_conversion')
                if out.get('single_pass_CH4_conversion') is not None
                else default_conv
            )
            if has_ch4_conv
            else (out.get('single_pass_feedstock_conversion') or default_conv)
        )
        default_cycle_conv: float = float(
            (
                out.get('cycle_averaged_CH4_conversion')
                if out.get('cycle_averaged_CH4_conversion') is not None
                else default_conv
            )
            if has_ch4_conv
            else (out.get('cycle_averaged_feedstock_conversion') or default_conv)
        )
        out['feedstock_conversion'] = default_conv
        out['single_pass_feedstock_conversion'] = default_single_pass_conv
        out['cycle_averaged_feedstock_conversion'] = default_cycle_conv
        out['CH4_conversion'] = default_conv
        out['single_pass_CH4_conversion'] = default_single_pass_conv
        if 'exit_x_CH4' in out and out['exit_x_CH4'] is not None:
            out['exit_x_feedstock'] = float(out['exit_x_CH4'])
        if (
            'ch4_conversion_rate_mol_m3_s' in out
            and out['ch4_conversion_rate_mol_m3_s'] is not None
        ):
            out['feedstock_conversion_rate_mol_m3_s'] = float(
                out['ch4_conversion_rate_mol_m3_s']
            )
        if (
            'endothermic_heat_duty_kJ_per_mol_CH4_fed' in out
            and out['endothermic_heat_duty_kJ_per_mol_CH4_fed'] is not None
        ):
            out['endothermic_heat_duty_kJ_per_mol_feed_fed'] = float(
                out['endothermic_heat_duty_kJ_per_mol_CH4_fed']
            )
        if (
            'endothermic_heat_duty_T_kJ_per_mol_CH4_fed' in out
            and out['endothermic_heat_duty_T_kJ_per_mol_CH4_fed'] is not None
        ):
            out['endothermic_heat_duty_T_kJ_per_mol_feed_fed'] = float(
                out['endothermic_heat_duty_T_kJ_per_mol_CH4_fed']
            )
        out.setdefault('exit_x_CO', 0.0)
        return out

    t_k: float = float(out.get('T_K') or 1073.15)
    p_pa: float = float(out.get('P_Pa') or 101325.0)
    rtype: str = str(out.get('reactor_type') or 'PFR')
    e_act_ch4_ref: float = float(out.get('catalyst_E_act_eV') or 0.80)
    e_act_feed: float = effective_catalyst_barrier_for_feedstock(e_act_ch4_ref, spec)
    tau_s: float = float(out.get('residence_time_s') or 1.25)

    # Baseline methane reference conversion and equilibrium (`*_ch4_ref`) from the Cantera CH4 pass:
    x_ch4_ref: float = float(
        out.get('single_pass_CH4_conversion')
        or out.get('CH4_conversion')
        or 0.0
    )
    x_eq_ch4_ref: float = float(
        out.get('X_eq_local_tp') or out.get('X_eq_table') or 0.985
    )
    x_eq_feed: float = feedstock_equilibrium_conversion(
        t_k, p_pa, spec, y_feed_inlet=y_feed_inlet
    )

    # Extract baseline effective Damköhler number Da_CH4_ref from X_CH4_ref / X_eq_CH4_ref
    frac_ch4_ref: float = min(max(x_ch4_ref / max(x_eq_ch4_ref, 1e-4), 0.0), 0.995)
    da_ch4_ref: float = -math.log(max(1.0 - frac_ch4_ref, 1e-4))

    # BEP rate acceleration from weaker C-C / secondary C-H / C-O bonds:
    # moderated on solid beds by solid-carbon stoichiometry (nu_C^-0.35 surface site competition)
    # whereas liquid-metal MMBCR floats carbon continuously without site blocking.
    ddg_j_mol: float = (e_act_ch4_ref - e_act_feed) * 96485.33212
    bep_rate_mult: float = math.exp(min(ddg_j_mol / (R_J_MOL_K * max(t_k, 300.0)), 6.0))
    coking_stoich_factor: float = (
        1.0
        if rtype == 'MMBCR'
        else (max(float(spec.nu_C), 1.0) ** -0.35)
    )
    if str(out.get('direct_cantera_feedstock') or '').upper() == spec.formula:
        x_feed = float(np.clip(x_ch4_ref, 0.0, x_eq_feed))
        frac_direct = min(max(x_feed / max(x_eq_feed, 1e-4), 0.0), 0.995)
        da_feed = -math.log(max(1.0 - frac_direct, 1e-4))
    else:
        da_feed = da_ch4_ref * bep_rate_mult * coking_stoich_factor
        x_feed = float(np.clip(x_eq_feed * (1.0 - math.exp(-da_feed)), 0.0, x_eq_feed))

    # Competitive homogeneous thermal cracking + secondary CH4/C2H4 intermediate slip:
    eps_val: float = float(
        out.get('gas_holdup_fraction')
        or out.get('bed_porosity')
        or (0.22 if rtype == 'MMBCR' else 0.40)
    )
    k_hom_s: float = spec.hom_a_1_s * math.exp(
        -(spec.hom_ea_kJ_mol * 1000.0) / (R_J_MOL_K * max(t_k, 300.0))
    )
    k_het_eff_s: float = max(da_feed / max(tau_s, 1e-3), 1e-4)
    eff_c_chain: float = max(float(spec.nu_C), 1.0)
    slip_ceiling: float = min(0.08 + 0.018 * math.log(eff_c_chain), 0.14)
    s_light_hc_carbon: float = (
        float(
            np.clip(
                slip_ceiling
                * (eps_val * k_hom_s)
                / (k_het_eff_s + eps_val * k_hom_s),
                0.0,
                slip_ceiling,
            )
        )
        if x_feed > 1e-6 and t_k >= 750.0
        else 0.0
    )

    solid_c_sel: float = float(np.clip(1.0 - s_light_hc_carbon, 0.0, 1.0))
    h_per_c_feed: float = float(spec.n_H) / max(float(spec.nu_C), 1.0)
    h_retained_frac: float = min(2.56 / max(h_per_c_feed, 2.0), 1.0) * s_light_hc_carbon
    h2_sel: float = float(np.clip(1.0 - h_retained_frac, 0.0, 1.0))
    h2_yield: float = float(x_feed * h2_sel)

    # Molar gas expansion and outlet mole fractions for active feedstock:
    # 1 mol feed -> nu_H2 * h2_sel mol H2 + nu_CO mol CO + light hydrocarbon byproducts
    dn_feed_converted: float = y_feed_inlet * x_feed
    n_h2_out: float = spec.nu_H2 * dn_feed_converted * h2_sel
    n_co_out: float = spec.nu_CO * dn_feed_converted
    n_c2_tot: float = 0.5 * float(spec.nu_C) * dn_feed_converted * s_light_hc_carbon
    n_feed_out: float = y_feed_inlet * (1.0 - x_feed)
    n_tot_out: float = max(
        n_feed_out + n_h2_out + n_co_out + n_c2_tot + y_ar_inlet, 1e-12
    )

    # Volumetric feedstock conversion & H2 production rates
    c_total_mol_m3: float = p_pa / (R_J_MOL_K * max(t_k, 1.0))
    eps_eff: float = min(max(eps_val, 0.01), 0.99)
    r_feed_mol_m3_s: float = float(
        eps_eff * c_total_mol_m3 * y_feed_inlet * x_feed / max(tau_s, 1e-6)
    )
    r_h2_mol_m3_s: float = float(spec.nu_H2 * r_feed_mol_m3_s * h2_sel)
    r_h2_kg_m3_h: float = float(r_h2_mol_m3_s * 0.00201588 * 3600.0)

    raw_sv: float | None = out.get('active_sv_1_m') or out.get('interfacial_sv_ratio_1_m')
    sv_val: float = (
        float(raw_sv) * (eps_eff if rtype == 'MMBCR' and 'active_sv_1_m' not in out else 1.0)
        if isinstance(raw_sv, (int, float)) and float(raw_sv) > 0.0
        else (900.0 * eps_eff if rtype == 'MMBCR' else 3600.0)
    )
    site_density_mol_cm2: float = float(out.get('site_density_mol_cm2') or 2.5e-9)
    sites_mol_m3: float = site_density_mol_cm2 * 1e4 * sv_val
    tof_h2_1_s: float = float(r_h2_mol_m3_s / max(sites_mol_m3, 1e-12))

    dh_t_kJ_mol_feed: float = feedstock_reaction_enthalpy_kJ_mol(t_k, spec)
    heat_duty_298_kJ_mol_feed: float = round(spec.dh_rxn_298_kJ_mol_feed * x_feed, 4)
    heat_duty_t_kJ_mol_feed: float = round(dh_t_kJ_mol_feed * x_feed, 4)
    exit_x_feed: float = round(n_feed_out / n_tot_out, 6)

    # Populate both feedstock-general keys (`feedstock_*`) and synchronized
    # legacy `CH4_*` keys so downstream scorecard/TEA/reporting consumers work
    # seamlessly whether they inspect `feedstock_*` or `CH4_*`.
    out.update(
        {
            'feedstock_conversion': round(x_feed, 6),
            'single_pass_feedstock_conversion': round(x_feed, 6),
            'cycle_averaged_feedstock_conversion': round(x_feed, 6),
            'CH4_conversion': round(x_feed, 6),
            'single_pass_CH4_conversion': round(x_feed, 6),
            'cycle_averaged_CH4_conversion': round(x_feed, 6),
            'X_eq_local_tp': round(x_eq_feed, 6),
            'catalyst_E_act_eV': round(e_act_feed, 6),
            'catalyst_E_act_CH4_reference_eV': round(e_act_ch4_ref, 6),
            'H2_selectivity': round(h2_sel, 6),
            'H2_atom_balance': round(h2_sel if x_feed > 0.01 else 0.0, 6),
            'H2_yield': round(h2_yield, 6),
            'solid_C_selectivity': round(solid_c_sel, 6),
            'C2_carbon_selectivity': round(s_light_hc_carbon, 6),
            'exit_x_H2': round(n_h2_out / n_tot_out, 6),
            'exit_x_CO': round(n_co_out / n_tot_out, 6),
            'exit_x_feedstock': exit_x_feed,
            'exit_x_CH4': exit_x_feed,
            'exit_x_C2H6': round(0.20 * n_c2_tot / n_tot_out, 8),
            'exit_x_C2H4': round(0.65 * n_c2_tot / n_tot_out, 8),
            'exit_x_C2H2': round(0.15 * n_c2_tot / n_tot_out, 8),
            'feedstock_conversion_rate_mol_m3_s': round(r_feed_mol_m3_s, 6),
            'ch4_conversion_rate_mol_m3_s': round(r_feed_mol_m3_s, 6),
            'h2_production_rate_mol_m3_s': round(r_h2_mol_m3_s, 6),
            'h2_production_rate_kg_m3_h': round(r_h2_kg_m3_h, 6),
            'h2_turnover_frequency_1_s': round(tof_h2_1_s, 6),
            'endothermic_heat_duty_kJ_per_mol_feed_fed': heat_duty_298_kJ_mol_feed,
            'endothermic_heat_duty_T_kJ_per_mol_feed_fed': heat_duty_t_kJ_mol_feed,
            'endothermic_heat_duty_kJ_per_mol_CH4_fed': heat_duty_298_kJ_mol_feed,
            'endothermic_heat_duty_T_kJ_per_mol_CH4_fed': heat_duty_t_kJ_mol_feed,
            'c2_tracked': True,
        }
    )
    return out
