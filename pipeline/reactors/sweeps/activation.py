"""Phase 2 continuous carbon-detachment ablation response map by reactor family."""

from __future__ import annotations

# PFR and PBR clear C_s only by discrete regen. They do not read a continuous
# detachment rate; a flat PFR/PBR pair is "not wired", not "detachment = 0".
ABLATION_RESPONSE: dict[str, bool] = {
    'MMBCR': True,
    'Fluidized': True,
    'PFR': False,
    'PBR': False,
}
