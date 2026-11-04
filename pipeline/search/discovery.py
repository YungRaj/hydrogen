"""Deterministic, coverage-aware search candidate acquisition.

The evolutionary optimizers are useful local exploiters, but a GA alone cannot
show which parts of a combinatorial space were never visited.  This module
provides canonical candidate identities, chemistry-region coverage accounting,
and a deterministic quality/novelty acquisition rule shared by both campaigns.
"""

from __future__ import annotations

import hashlib
import json
import ast
from typing import Iterable, List, Optional, Sequence, Tuple

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]


def fast_non_dominated_sort(objectives: np.ndarray) -> List[List[int]]:
    """Partition minimization objectives into Pareto fronts.

    Args:
        objectives: Candidate rows with one or more minimized objectives.

    Returns:
        Candidate indices grouped from the best to worst Pareto front.
    """
    remaining = np.arange(len(objectives))
    fronts: List[List[int]] = []
    while len(remaining):
        subset = objectives[remaining]
        efficient = np.ones(len(subset), dtype=bool)
        for index in range(len(subset)):
            if efficient[index]:
                dominated = np.all(subset[index] <= subset, axis=1) & np.any(
                    subset[index] < subset, axis=1
                )
                efficient[dominated] = False
        front_indices = np.flatnonzero(efficient)
        fronts.append(remaining[front_indices].tolist())
        remaining = np.delete(remaining, front_indices)
    return fronts


def canonicalize_genome(genome: tuple) -> tuple:
    """Return a stable representation without changing site semantics.

    Args:
        genome: Encoded catalyst composition and structural configuration.

    Returns:
        Ordered tuple of computed values.
    """
    g = list(genome)
    if not g:
        raise ValueError("empty catalyst genome")
    if g[0] == "SolidCatalyst":
        # Dopant order is not represented by the generated slab; permutations
        # therefore describe the same candidate.  Preserve multiplicity.
        g[5] = tuple(sorted(g[5]))
        g[4] = round(float(g[4]), 6)
    elif g[0] == "MoltenMetal":
        # Zero loading and self-promotion both describe the pure host.
        if g[2] == "None" or float(g[3]) == 0.0 or g[2] == g[1]:
            g[2], g[3] = "None", 0.0
    elif g[0] == "Perovskite":
        # Dopant identity has no physical meaning at zero fraction.
        if float(g[4]) == 0.0:
            g[3] = "None"
    elif g[0] == "HEA":
        g[1] = tuple(sorted(g[1]))
    elif g[0] == "MetalHydride" and g[3] == g[1]:
        # A repeated secondary metal is the single-metal composition.
        g[3] = "None"
    elif g[0] == "MAXPhase" and g[5] == g[1]:
        g[5] = "None"
    elif g[0] == "MXene" and g[5] == g[1]:
        g[5] = "None"
    return tuple(g)


def candidate_id(genome: tuple) -> str:
    """Content-addressed ID suitable for resumable, sharded campaigns.

    Args:
        genome: Encoded catalyst composition and structural configuration.

    Returns:
        Computed `str` result.
    """
    payload = json.dumps(canonicalize_genome(genome), separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


def discovery_region(genome: tuple) -> Tuple[str, ...]:
    """A coarse chemistry cell used to prevent global-rank mode collapse.

    Args:
        genome: Encoded catalyst composition and structural configuration.

    Returns:
        Ordered tuple of computed values.
    """
    g = canonicalize_genome(genome)
    cls = g[0]
    if cls == "SolidCatalyst":
        return cls, str(g[1]), str(g[2]), str(g[3])
    if cls == "MoltenMetal":
        return cls, str(g[1]), str(g[2])
    if cls in ("SAC", "DAC"):
        return tuple([cls] + [str(x) for x in g[1:4]])
    if cls in ("MOF", "COF"):
        return cls, str(g[1]), str(g[2]), str(g[3])
    if cls == "HEA":
        return cls, "-".join(g[1]), str(g[2])
    return tuple([cls] + [str(x) for x in g[1:3]])



LITERATURE_VALIDATION_GATES: dict[str, dict[str, object]] = {
    'deutschmann2000_monolayer_site_inventory': {
        'literature_gate_name': (
            'Deutschmann2000_MonolayerSiteDensity_PelletInventory_Gate'
        ),
        'short_title': 'Monolayer Site Density & Pellet Inventory Scaling (Deutschmann 2000; Alves 2021)',
        'governing_equation': (
            'Gamma = 2.5e-9 mol/cm^2; a_active = a_geom * metal_loading * metal_dispersion'
        ),
        'literature_citations': [
            'Deutschmann et al., Catal. Today 59 (2000) 141-150',
            'Alves et al., Renew. Sustain. Energy Rev. 137 (2021) 110465',
            'Sanchez-Bastardo et al., Ind. Eng. Chem. Res. 60 (2021) 11855-11881',
        ],
    },
    'muradov2005_zero_co2_solids_outfeed_regeneration': {
        'literature_gate_name': (
            'TurquoisePyrolysis_ZeroCO2_CarbonRegeneration_Policy'
        ),
        'short_title': 'Zero-CO2 Mechanical & Circulating Carbon Outfeed Policy (Muradov 2005; Parkinson 2019)',
        'governing_equation': (
            'co2_permitted = False; mechanical/circulating removal clears C_s without clearing encapsulating C_encap_s'
        ),
        'literature_citations': [
            'Muradov & Veziroglu, Int. J. Hydrogen Energy 30 (2005) 225-237',
            'Parkinson et al., Energy Environ. Sci. 12 (2019) 19-40',
        ],
    },
    'upham2017_mendelson_bubble_flotation_interfacial': {
        'literature_gate_name': (
            'Upham2017_MendelsonBubbleRise_MoltenInterfacialMassTransfer_Closure'
        ),
        'short_title': 'Molten-Metal Mendelson Bubble Rise & Interfacial Flotation Closure (Upham 2017; Chen 2023)',
        'governing_equation': (
            'u_b = sqrt(2*sigma/(rho*d_b) + g*d_b/2); eps_g = u_sup/u_b <= 0.30; k0 in (0, 1.0] m/s'
        ),
        'literature_citations': [
            'Upham et al., Science 358 (2017) 917-921',
            'Palmer et al., ACS Catal. 10 (2020) 9945-9952',
            'Chen et al., Ind. Eng. Chem. Res. 62 (2023) 1832-1844',
        ],
    },
    'ergun1952_packed_bed_pressure_drop_pilot_residence': {
        'literature_gate_name': (
            'Ergun1952_PackedBedPressureDrop_ThieleEffectiveness_PilotResidence_Gate'
        ),
        'short_title': 'Packed-Bed Ergun Pressure Drop, Thiele Effectiveness & Pilot Residence Window (Ergun 1952)',
        'governing_equation': (
            'DeltaP_Ergun <= 1.01325 bar; tau_pilot in [2.0, 12.0] s; Omega = eta_internal / (1 + Da_ext * eta_internal)'
        ),
        'literature_citations': [
            'Ergun, Chem. Eng. Prog. 48 (1952) 89-94',
            'Froment, De Wilde & Bischoff, Chemical Reactor Analysis and Design (2011)',
        ],
    },
    'alves2021_xu2021_supported_ni_filament_ranking_authority': {
        'literature_gate_name': (
            'Alves2021_Xu2021_SupportedNi_FilamentousCarbon_ConversionClosure_Gate'
        ),
        'short_title': (
            'Supported-Ni(111)/SiO2 Filamentous-Carbon ROI Conversion & Ranking Authority Gate '
            '(Alves 2021; Xu 2021; Bengaard 2002)'
        ),
        'governing_equation': (
            'X_single_pass(650-700 C) >= 0.50 on PFR, PBR, and Fluidized with X <= X_eq and '
            'Baker-Helveg-Ermakova-Takenaka off-site carbon criteria passed'
        ),
        'literature_citations': [
            'Alves et al., Renew. Sustain. Energy Rev. 137 (2021) 110465',
            'Xu et al., Appl. Catal. A 611 (2021) 117967',
            'Bengaard et al., J. Catal. 209 (2002) 365-384',
        ],
    },
    'baker1972_helveg2004_ermakova2000_takenaka2003_off_site_carbon_nucleation': {
        'literature_gate_name': (
            'Baker1972_Helveg2004_AbildPedersen2006_Ermakova2000_Takenaka2003_'
            'OffSiteFilamentousCarbon_vs_Encapsulation_Gate'
        ),
        'short_title': (
            'Baker-Helveg-Abild-Pedersen Off-Site Carbon Diffusion (C_gamma) vs. Graphitic '
            'Encapsulation (C_delta) Calibration (Ermakova 2000; Takenaka 2003)'
        ),
        'governing_equation': (
            'r_gamma = A_gamma * exp(-E_gamma/RT) * theta_C; '
            'r_delta = (A_gamma/theta*) * exp(-E_delta/RT) * theta_C^2; '
            'Y_C in [8, 10] g_C/(g_Ni*h); tau_encap in [4, 50] h'
        ),
        'literature_citations': [
            'Baker et al., J. Catal. 26 (1972) 51-62',
            'Helveg et al., Nature 427 (2004) 426-429',
            'Abild-Pedersen et al., Phys. Rev. B 73 (2006) 115419',
            'Ermakova et al., Appl. Catal. A 201 (2000) 61-74',
            'Takenaka et al., J. Catal. 217 (2003) 79-87',
            'Amin et al., Ind. Eng. Chem. Res. 50 (2011) 12460-12470',
        ],
    },
}

NAMED_REFERENCE_CATALYSTS: dict[str, dict[str, str]] = {
    'ni111_sio2_bengaard2002_ref': {
        'catalyst_id': 'ni111_sio2_bengaard2002_ref',
        'chemical_formula': 'Ni(111)/SiO2',
        'literature_catalyst_name': (
            'Ni(111)/SiO2 Supported Nanoparticle Reference '
            '(Bengaard 2002 / Ermakova 2000 / Takenaka 2003 Filamentous-Carbon Judge)'
        ),
        'material_class': 'SolidCatalyst',
        'physical_mechanism': (
            'Baker-Helveg-Abild-Pedersen subsurface/step carbon diffusion to CNT/CNF '
            '(C_alpha -> C_gamma) competing with second-order graphene encapsulation (C_delta)'
        ),
    },
    'rh1_ceo2_akri2019_sac_ref': {
        'catalyst_id': 'rh1_ceo2_akri2019_sac_ref',
        'chemical_formula': 'Rh1/CeO2(111)',
        'literature_catalyst_name': (
            'Rh1/CeO2(111) Single-Atom Catalyst Control '
            '(Akri 2019 Non-Filamentous Stoichiometric Monolayer Baseline)'
        ),
        'material_class': 'SAC',
        'physical_mechanism': (
            'Isolated single-atom active site lacking extended metal bulk/step edges '
            'for Baker-Helveg carbon filament nucleation (Akri et al., Nat. Commun. 2019)'
        ),
    },
}

CANDIDATE_CLASS_LITERATURE_ARCHETYPES: dict[str, dict[str, str]] = {
    'MoltenMetal': {
        'archetype_name': 'Liquid Metal Alloy Bubble-Column Catalyst',
        'primary_citations': 'Upham et al., Science 2017; Palmer et al., ACS Catal. 2020; Chen et al., Nat. Catal. 2023',
    },
    'SolidCatalyst': {
        'archetype_name': 'Supported Extended-Surface Transition Metal Catalyst',
        'primary_citations': 'Bengaard et al., J. Catal. 2002; Alves et al., Appl. Catal. A 2021',
    },
    'SupportedNanoparticle': {
        'archetype_name': 'Supported Bimetallic Step-Edge Nanoparticle Catalyst',
        'primary_citations': 'Snoeck et al., J. Catal. 1997; Reshetenko et al., Carbon 2003; Ashok et al., Int. J. Hydrogen Energy 2007',
    },
    'SAC': {
        'archetype_name': 'Atomically Dispersed Single-Atom Catalyst',
        'primary_citations': 'Zitolo et al., Nat. Mater. 2015; Akri et al., Nat. Commun. 2019',
    },
    'DAC': {
        'archetype_name': 'Synergistic Bimetallic Dual-Atom Catalyst',
        'primary_citations': 'Wang et al., J. Am. Chem. Soc. 2018; Li et al., Nat. Catal. 2020',
    },
    'SAA': {
        'archetype_name': 'Dilute Single-Atom Alloy Catalyst',
        'primary_citations': 'Kyriakou et al., Science 2012; Darby et al., Nat. Mater. 2018',
    },
    'HEA': {
        'archetype_name': 'Multi-Principal-Element High-Entropy Alloy Catalyst',
        'primary_citations': 'Miracle & Senkov, Acta Mater. 2017; Yao et al., Science 2018',
    },
    'Intermetallic': {
        'archetype_name': 'Ordered Intermetallic Compound Catalyst',
        'primary_citations': 'Stamenkovic et al., Science 2007; Wang et al., Nat. Mater. 2013',
    },
    'Perovskite': {
        'archetype_name': 'ABO3 Perovskite Oxide Catalyst',
        'primary_citations': 'Suntivich et al., Science 2011',
    },
    'Spinel': {
        'archetype_name': 'AB2O4 Spinel Oxide Catalyst',
        'primary_citations': 'Liang et al., Nat. Mater. 2011',
    },
    'MAXPhase': {
        'archetype_name': 'Mn+1AXn Nanolaminate MAX Phase Catalyst',
        'primary_citations': 'Barsoum, Prog. Solid State Chem. 2000',
    },
    'MXene': {
        'archetype_name': '2D Transition-Metal Carbide/Nitride MXene Catalyst',
        'primary_citations': 'Naguib et al., Adv. Mater. 2011; Seh et al., ACS Energy Lett. 2016',
    },
    'TMDC': {
        'archetype_name': 'Transition-Metal Dichalcogenide Edge/Basal Catalyst',
        'primary_citations': 'Jaramillo et al., Science 2007',
    },
    'TransitionMetalCompound': {
        'archetype_name': 'Interstitial Transition-Metal Carbide/Nitride/Phosphide Catalyst',
        'primary_citations': 'Levy & Boudart, Science 1973',
    },
    'MetalHydride': {
        'archetype_name': 'Subsurface Transition-Metal Hydride Catalyst',
        'primary_citations': 'Fukuzumi, Eur. J. Inorg. Chem. 2008; Wang et al., Nat. Chem. 2016',
    },
    'MOF': {
        'archetype_name': 'Coordinatively Unsaturated Metal-Organic Framework Catalyst',
        'primary_citations': 'Furukawa et al., Science 2013',
    },
    'COF': {
        'archetype_name': 'Metallated Covalent Organic Framework Catalyst',
        'primary_citations': 'Furukawa et al., Science 2013; Lin et al., Science 2015',
    },
    'MetalFreeCarbon': {
        'archetype_name': 'Heteroatom-Doped Defect-Rich Metal-Free Carbon Catalyst',
        'primary_citations': 'Muradov et al., Catal. Today 2005; Gong et al., Science 2009',
    },
}


def _format_genome_tuple_label(genome: tuple) -> str:
    """Format a canonical catalyst genome tuple into a literature-backed chemical label."""
    if not genome:
        return 'Unknown Catalyst'
    cls = str(genome[0])
    if cls == 'SolidCatalyst':
        # Distinguish pyrolysis SolidCatalyst: ('SolidCatalyst', metal, support, facet, strain, dopants, n_dop, n_vac)
        # from fuel-cell SolidCatalyst: ('SolidCatalyst', m1, m2, ratio_pct, support, facet, strain, n_vac)
        if len(genome) >= 6 and isinstance(genome[5], (tuple, list)):
            metal = str(genome[1])
            support = str(genome[2])
            facet = str(genome[3])
            strain = float(genome[4]) if isinstance(genome[4], (int, float)) else 0.0
            dopants = [str(d) for d in genome[5]] if genome[5] else []
            n_dop = int(genome[6]) if len(genome) > 6 and isinstance(genome[6], (int, float)) else 0
            dop_prefix = f"{','.join(dopants)}-doped " if (dopants and n_dop > 0) else ""
            strain_str = f", strain={strain:+.2f}%" if abs(strain) > 1e-4 else ""
            return (
                f"{dop_prefix}{metal}({facet})/{support}{strain_str} Supported Solid Catalyst "
                f"(Bengaard 2002 / Alves 2021)"
            )
        if len(genome) >= 6:
            m1, m2, ratio, support, facet = (
                str(genome[1]),
                str(genome[2]),
                genome[3],
                str(genome[4]),
                str(genome[5]),
            )
            alloy = f"{m1}{ratio}{m2}{100 - int(ratio)}" if isinstance(ratio, (int, float)) and m2 != m1 else m1
            return f"{alloy}({facet})/{support} Supported Alloy Catalyst (Stamenkovic 2007 / Greeley 2009)"
    if cls == 'SupportedNanoparticle' and len(genome) >= 6:
        m1, m2 = str(genome[1]), str(genome[2])
        counts = genome[3] if isinstance(genome[3], (tuple, list)) and len(genome[3]) == 2 else (1, 0)
        support, facet = str(genome[4]), str(genome[5])
        comp = f"{m1}{counts[0]}{m2}{counts[1]}" if (m2 not in ('None', 'none', m1) and int(counts[1]) > 0) else f"{m1}{counts[0]}"
        return (
            f"{comp}/{support}({facet}) Supported Bimetallic Nanoparticle "
            f"(Reshetenko 2003 / Ashok 2007 / Snoeck 1997 Step-Edge Alloy)"
        )
    if cls == 'MoltenMetal' and len(genome) >= 4:
        m1, m2, frac = str(genome[1]), str(genome[2]), genome[3]
        if m2 in ('None', 'none', m1) or not isinstance(frac, (int, float)) or float(frac) == 0.0:
            comp = f"Liquid {m1}"
        else:
            val = float(frac)
            if val <= 1.0:
                comp = f"{m2}{val:.2f}{m1}{1.0 - val:.2f} Liquid Metal Alloy"
            else:
                comp = f"{m1}{100 - int(val)}{m2}{int(val)} Liquid Metal Alloy"
        return f"{comp} (Upham et al., Science 2017 / Palmer et al., ACS Catal. 2020 Bubble-Column Melt)"
    if cls == 'HEA' and len(genome) >= 4:
        metals = "".join(str(m) for m in genome[1]) if isinstance(genome[1], (tuple, list)) else str(genome[1])
        struct = str(genome[2])
        facet = str(genome[3])
        return (
            f"{metals}-{struct}({facet}) High-Entropy Alloy "
            f"(Miracle & Senkov 2017 / Yao et al., Science 2018 Multi-Principal-Element Surface)"
        )
    if cls == 'SAA' and len(genome) >= 4:
        m_dop, m_host, facet = str(genome[1]), str(genome[2]), str(genome[3])
        return (
            f"{m_dop}1/{m_host}({facet}) Single-Atom Alloy "
            f"(Kyriakou et al., Science 2012 / Darby et al., Nat. Mater. 2018)"
        )
    if cls == 'SAC' and len(genome) >= 4:
        metal, coord, axial = str(genome[1]), str(genome[2]), str(genome[3])
        support = str(genome[4]) if len(genome) > 4 else 'Carbon'
        ax = f"-{axial}" if axial not in ('None', 'none') else ""
        return (
            f"{metal}-{coord}{ax}/{support} Single-Atom Catalyst "
            f"(Zitolo et al., Nat. Mater. 2015 / Akri et al., Nat. Commun. 2019)"
        )
    if cls == 'DAC' and len(genome) >= 4:
        m1, m2, coord = str(genome[1]), str(genome[2]), str(genome[3])
        support = str(genome[5]) if len(genome) > 5 else 'Graphene'
        return (
            f"{m1}{m2}-{coord}/{support} Dual-Atom Catalyst "
            f"(Wang et al., J. Am. Chem. Soc. 2018 Synergistic Bimetallic Site)"
        )
    if cls == 'Intermetallic' and len(genome) >= 4:
        m1, m2, phase, facet = str(genome[1]), str(genome[2]), str(genome[3]), str(genome[4])
        return (
            f"{phase}-{m1}{m2}({facet}) Ordered Intermetallic "
            f"(Stamenkovic et al., Science 2007 / Wang et al., Nat. Mater. 2013)"
        )
    if cls == 'Perovskite' and len(genome) >= 4:
        a_site, b_site, dop = str(genome[1]), str(genome[2]), str(genome[3])
        frac = float(genome[4]) if len(genome) > 4 and isinstance(genome[4], (int, float)) else 0.0
        term = str(genome[5]) if len(genome) > 5 else '001'
        formula = f"{a_site}{1.0 - frac:.2g}{dop}{frac:.2g}{b_site}O3" if (dop not in ('None', 'none') and frac > 0) else f"{a_site}{b_site}O3"
        return f"{formula}({term}) Perovskite Oxide (Suntivich et al., Science 2011 eg-Filling Oxide)"
    if cls == 'Spinel' and len(genome) >= 5:
        a_site, b_site, dop, morph = str(genome[1]), str(genome[2]), str(genome[3]), str(genome[4])
        sup = f"/{genome[5]}" if len(genome) > 5 and str(genome[5]) not in ('None', 'none') else ""
        dop_str = f"{dop}-doped " if dop not in ('None', 'none') else ""
        return f"{dop_str}{a_site}{b_site}2O4{sup} ({morph}) Spinel Oxide (Liang et al., Nat. Mater. 2011)"
    if cls == 'MAXPhase' and len(genome) >= 6:
        m_el, a_el, x_el, n_val, dop = str(genome[1]), str(genome[2]), str(genome[3]), int(genome[4]), str(genome[5])
        facet = str(genome[6]) if len(genome) > 6 else '0001'
        m_str = f"({m_el},{dop})" if dop not in ('None', 'none', m_el) else m_el
        return f"{m_str}{n_val + 1}{a_el}{x_el}{n_val}({facet}) MAX Phase Nanolaminate (Barsoum, Prog. Solid State Chem. 2000)"
    if cls == 'MXene' and len(genome) >= 5:
        m_el, x_el, n_val, term = str(genome[1]), str(genome[2]), int(genome[3]), str(genome[4])
        dop = str(genome[5]) if len(genome) > 5 else 'None'
        m_str = f"({m_el},{dop})" if dop not in ('None', 'none', m_el) else m_el
        return f"{m_str}{n_val + 1}{x_el}{n_val}{term}x 2D MXene (Naguib et al., Adv. Mater. 2011 / Seh 2016)"
    if cls == 'TMDC' and len(genome) >= 5:
        m_el, x_el, phase, dop = str(genome[1]), str(genome[2]), str(genome[3]), str(genome[4])
        defect = str(genome[5]) if len(genome) > 5 else 'edge'
        dop_str = f"{dop}-doped " if dop not in ('None', 'none') else ""
        return f"{dop_str}{phase}-{m_el}{x_el}2 ({defect}) Dichalcogenide (Jaramillo et al., Science 2007)"
    if cls == 'TransitionMetalCompound' and len(genome) >= 5:
        m_el, x_el, stoich, dop = str(genome[1]), str(genome[2]), str(genome[3]), str(genome[4])
        facet = str(genome[5]) if len(genome) > 5 else '001'
        dop_str = f"{dop}-doped " if dop not in ('None', 'none') else ""
        return f"{dop_str}{m_el}{x_el} ({stoich}, {facet}) Interstitial Compound (Levy & Boudart, Science 1973)"
    if cls == 'MetalHydride' and len(genome) >= 4:
        m_el, h_st, dop = str(genome[1]), str(genome[2]), str(genome[3])
        dop_str = f"{dop}-doped " if dop not in ('None', 'none', m_el) else ""
        return f"{dop_str}{m_el}H{h_st} Subsurface Metal Hydride (Fukuzumi 2008 / Wang 2016)"
    if cls in ('MOF', 'COF') and len(genome) >= 4:
        node, linker, func = str(genome[1]), str(genome[2]), str(genome[3])
        topo = str(genome[4]) if len(genome) > 4 else 'net'
        return f"{node}-{linker}-{func} ({topo}) {cls} Network (Furukawa et al., Science 2013)"
    if cls == 'MetalFreeCarbon' and len(genome) >= 4:
        d1, d2, defect = str(genome[1]), str(genome[2]), str(genome[3])
        scaffold = str(genome[4]) if len(genome) > 4 else 'nanocarbon'
        dops = ",".join(d for d in (d1, d2) if d not in ('None', 'none'))
        dop_str = f"{dops}-doped " if dops else "Undoped "
        return f"{dop_str}{defect} {scaffold} Metal-Free Carbon (Muradov 2005 / Gong et al., Science 2009)"
    return f"{cls}({'-'.join(str(x) for x in genome[1:])})"


_CANDIDATE_GENOME_CACHE: dict[str, tuple] = {}


def _populate_candidate_genome_cache() -> None:
    if _CANDIDATE_GENOME_CACHE:
        return
    from pathlib import Path
    from pipeline.utils import BASE_DIR

    results_root = BASE_DIR / 'results'
    if not results_root.is_dir():
        return
    for kin_file in results_root.glob('**/mechanism_*.kinetics.json'):
        try:
            data = json.loads(kin_file.read_text(encoding='utf-8'))
            inputs = data.get('inputs') or {}
            cid = str(inputs.get('candidate_id') or data.get('catalyst_name') or '')
            raw_g = inputs.get('genome')
            if cid and raw_g:
                parsed = (
                    ast.literal_eval(raw_g)
                    if isinstance(raw_g, str)
                    else tuple(
                        tuple(x) if isinstance(x, list) else x for x in raw_g
                    )
                )
                if isinstance(parsed, tuple):
                    _CANDIDATE_GENOME_CACHE[cid] = parsed
        except Exception:
            continue
    for csv_file in results_root.glob('**/*.csv'):
        try:
            text = csv_file.read_text(encoding='utf-8', errors='ignore')
            if 'genome' not in text:
                continue
            import csv as _csv
            reader = _csv.DictReader(text.splitlines())
            for row in reader:
                raw_g = row.get('genome')
                if not raw_g:
                    continue
                try:
                    parsed = ast.literal_eval(raw_g)
                    if isinstance(parsed, tuple):
                        cid = row.get('candidate_id') or candidate_id(parsed)
                        _CANDIDATE_GENOME_CACHE[str(cid)] = parsed
                except Exception:
                    continue
        except Exception:
            continue


def format_literature_catalyst_label(
    genome_or_id: object,
    *,
    material_class: object = None,
    genome: object = None,
) -> str:
    """Resolve a catalyst genome or candidate ID into a self-explanatory literature-backed name.

    Args:
        genome_or_id: Catalyst genome tuple/list/string or candidate identifier
            (`ni111_sio2_bengaard2002_ref`, `rh1_ceo2_akri2019_sac_ref`, or SHA-1 ID).
        material_class: Optional material class string when `genome` is not provided.
        genome: Optional explicit genome tuple/list/string to format when
            `genome_or_id` is a candidate ID.

    Returns:
        Human-readable literature-backed chemical formula, structure, and citation label.
    """
    if genome is not None and str(genome).strip():
        base_label = format_literature_catalyst_label(genome)
        if isinstance(genome_or_id, str) and genome_or_id.strip():
            cid = genome_or_id.strip()
            if cid not in base_label:
                return f"{base_label} [id={cid}]"
        return base_label
    if genome_or_id is None:
        return 'Unspecified Catalyst'
    if isinstance(genome_or_id, str):
        text = genome_or_id.strip()
        if text in NAMED_REFERENCE_CATALYSTS:
            return NAMED_REFERENCE_CATALYSTS[text]['literature_catalyst_name']
        if text.startswith('(') or text.startswith('['):
            try:
                parsed = ast.literal_eval(text)
                if isinstance(parsed, (tuple, list)):
                    return _format_genome_tuple_label(
                        tuple(tuple(x) if isinstance(x, list) else x for x in parsed)
                    )
            except (ValueError, SyntaxError):
                pass
        _populate_candidate_genome_cache()
        if text in _CANDIDATE_GENOME_CACHE:
            return (
                f"{_format_genome_tuple_label(_CANDIDATE_GENOME_CACHE[text])} "
                f"[id={text}]"
            )
        if material_class and str(material_class) in CANDIDATE_CLASS_LITERATURE_ARCHETYPES:
            arch = CANDIDATE_CLASS_LITERATURE_ARCHETYPES[str(material_class)]
            return f"{arch['archetype_name']} ({arch['primary_citations']}) [id={text}]"
        return text
    if isinstance(genome_or_id, (tuple, list)):
        return _format_genome_tuple_label(
            tuple(tuple(x) if isinstance(x, list) else x for x in genome_or_id)
        )
    return str(genome_or_id)


def add_discovery_metadata(frame):
    """Add stable identity, coverage-cell, and literature-backed name columns to a screening DataFrame.

    Args:
        frame: Tabular candidate or result records.

    Returns:
        Computed result described above.
    """
    if frame is None or "genome" not in frame.columns:
        return frame
    ids, regions, labels = [], [], []
    for raw in frame["genome"]:
        try:
            genome = ast.literal_eval(raw) if isinstance(raw, str) else tuple(raw)
            ids.append(candidate_id(genome))
            regions.append("|".join(discovery_region(genome)))
            labels.append(_format_genome_tuple_label(genome))
        except (ValueError, SyntaxError, TypeError):
            ids.append(None)
            regions.append(None)
            labels.append(None)
    frame = frame.copy()
    frame["candidate_id"] = ids
    frame["discovery_region"] = regions
    frame["literature_catalyst_name"] = labels
    return frame



def _quality_score(
    objectives: np.ndarray,
    active_indices: Optional[Sequence[int]] = None,
) -> np.ndarray:
    """Robustly scale minimization objectives and combine without units."""
    obj = np.asarray(objectives, dtype=float)
    if obj.ndim != 2 or len(obj) == 0:
        return np.zeros(len(obj), dtype=float)
    if active_indices is not None:
        cols = [int(i) for i in active_indices if 0 <= int(i) < obj.shape[1]]
        if cols:
            obj = obj[:, cols]
    finite = np.where(np.isfinite(obj), obj, np.nan)
    lo = np.nanpercentile(finite, 5, axis=0)
    hi = np.nanpercentile(finite, 95, axis=0)
    scaled = np.clip((finite - lo) / np.maximum(hi - lo, 1e-12), 0.0, 1.0)
    scaled = np.nan_to_num(scaled, nan=1.0, posinf=1.0, neginf=0.0)
    return 1.0 - scaled.mean(axis=1)


def select_discovery_batch(
    candidates: Sequence[tuple],
    objectives: np.ndarray,
    n_select: int,
    evaluated: Iterable[tuple] = (),
    uncertainties: Optional[Sequence[float]] = None,
    confidence: Optional[Sequence[float]] = None,
) -> List[int]:
    """Select a deterministic batch balancing viability, novelty, and OOD value.

        First take the best candidate from as many unseen chemistry cells as the
        budget permits.  Remaining slots use a combined score.  Low confidence is
        treated as a reason to *validate* a candidate, not as proof it is poor.

    Args:
        candidates: Candidate records to process.
        objectives: Objectives used by this operation.
        n_select: Number of select to use.
        evaluated: Evaluated used by this operation.
        uncertainties: Ordered values supplying uncertainties.
        confidence: Ordered values supplying confidence.

    Returns:
        List of computed or validated records.
    """
    if n_select <= 0 or not candidates:
        return []
    n_select = min(n_select, len(candidates))
    quality = _quality_score(objectives)
    seen_ids = {candidate_id(g) for g in evaluated}
    seen_regions = {discovery_region(g) for g in evaluated}
    uncertainty = (
        np.zeros(len(candidates))
        if uncertainties is None
        else np.asarray(uncertainties, float)
    )
    if len(uncertainty) != len(candidates):
        raise ValueError("uncertainties must match candidates")
    uncertainty_span = np.ptp(uncertainty)
    if uncertainty_span > 0:
        uncertainty = (uncertainty - uncertainty.min()) / uncertainty_span
    conf = (
        np.ones(len(candidates))
        if confidence is None
        else np.asarray(confidence, float)
    )
    if len(conf) != len(candidates):
        raise ValueError("confidence must match candidates")

    ids = [candidate_id(g) for g in candidates]
    regions = [discovery_region(g) for g in candidates]
    eligible = [i for i, cid in enumerate(ids) if cid not in seen_ids]
    by_region = {}
    for i in eligible:
        by_region.setdefault(regions[i], []).append(i)

    # Region champions: quality remains important, with deterministic ID ties
    # resolved by the same total order as the ranking sort below.
    champions = [
        min(indices, key=lambda i: (-quality[i], -uncertainty[i], ids[i]))
        for region, indices in sorted(by_region.items())
        if region not in seen_regions
    ]
    champions.sort(key=lambda i: (-quality[i], -uncertainty[i], ids[i]))
    selected = champions[:n_select]
    selected_set = set(selected)

    # OOD value is capped: unfamiliar chemistry receives validation opportunity
    # but cannot swamp a batch solely because the model knows little about it.
    score = 0.60 * quality + 0.20 * uncertainty + 0.20 * np.clip(1.0 - conf, 0.0, 1.0)
    remainder = sorted(
        (i for i in eligible if i not in selected_set),
        key=lambda i: (-score[i], ids[i]),
    )
    selected.extend(remainder[: n_select - len(selected)])
    return selected
