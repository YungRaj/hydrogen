"""Immutable scientific protocol definitions for atomistic screening."""

from dataclasses import dataclass


@dataclass(frozen=True)
class RelaxationBudget:
    fmax_eV_A: float
    steps: int


@dataclass(frozen=True)
class ScreeningProtocol:
    protocol_id: str
    reference: RelaxationBudget
    clean: RelaxationBudget
    adsorbate: RelaxationBudget


PYROLYSIS_PROTOCOL = ScreeningProtocol(
    protocol_id='esen-sm-conserving-all-oc25:relax-v3:pyrolysis-v2',
    reference=RelaxationBudget(0.05, 200),
    clean=RelaxationBudget(0.08, 150),
    adsorbate=RelaxationBudget(0.08, 100),
)

ORR_PROTOCOL = ScreeningProtocol(
    protocol_id='esen-sm-conserving-all-oc25:relax-v3:orr-che-v2',
    reference=RelaxationBudget(0.05, 200),
    clean=RelaxationBudget(0.08, 150),
    adsorbate=RelaxationBudget(0.08, 100),
)
