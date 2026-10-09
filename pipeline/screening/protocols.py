"""Immutable scientific protocol definitions for atomistic screening."""

from dataclasses import dataclass


@dataclass(frozen=True)
class RelaxationBudget:
    """Define optimizer limits for one screening relaxation tier.

    Attributes:
        fmax_eV_A: Configured fmax eV A value.
        steps: Configured steps value.
    """

    fmax_eV_A: float
    steps: int


@dataclass(frozen=True)
class ScreeningProtocol:
    """Define a named, auditable hierarchy of relaxation budgets.

    Attributes:
        protocol_id: Configured protocol id value.
        reference: Configured reference value.
        clean: Configured clean value.
        adsorbate: Configured adsorbate value.
    """

    protocol_id: str
    reference: RelaxationBudget
    clean: RelaxationBudget
    adsorbate: RelaxationBudget


PYROLYSIS_PROTOCOL = ScreeningProtocol(
    protocol_id='esen-sm-conserving-all-oc25:relax-v4:pyrolysis-v3',
    reference=RelaxationBudget(0.05, 200),
    clean=RelaxationBudget(0.08, 150),
    adsorbate=RelaxationBudget(0.08, 100),
)

ORR_PROTOCOL = ScreeningProtocol(
    protocol_id='esen-sm-conserving-all-oc25:relax-v4:orr-che-v3',
    reference=RelaxationBudget(0.05, 200),
    clean=RelaxationBudget(0.08, 150),
    adsorbate=RelaxationBudget(0.08, 100),
)
