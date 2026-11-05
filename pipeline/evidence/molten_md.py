"""Deterministic MLIP molecular-dynamics protocol contract for liquid alloys."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LiquidMDProtocol:
    """Controls equilibration and decorrelated snapshot collection."""

    temperature_K: float
    timestep_fs: float = 1.0
    equilibration_steps: int = 2_000
    stride_steps: int = 200
    snapshots: int = 10
    friction_per_fs: float = 0.01
    seed: int = 0

    def validate(self) -> None:
        """Reject physically or numerically invalid MD settings."""
        if self.temperature_K <= 0.0:
            raise ValueError('temperature_K must be positive')
        if self.timestep_fs <= 0.0 or self.friction_per_fs <= 0.0:
            raise ValueError('timestep and friction must be positive')
        if (
            self.equilibration_steps < 0
            or self.stride_steps <= 0
            or self.snapshots <= 0
        ):
            raise ValueError(
                'MD step counts must be non-negative and snapshots positive'
            )


def run_liquid_alloy_md(
    genome: tuple[object, ...],
    calc: object,
    protocol: LiquidMDProtocol | None = None,
) -> dict[str, object]:
    """Run finite-temperature Langevin NVT molecular dynamics on a molten-metal slab.

    Args:
        genome: MoltenMetal catalyst genome tuple `('MoltenMetal', host, promoter, at_pct, T_K)`.
        calc: Attached ASE-compatible interatomic potential calculator.
        protocol: Optional `LiquidMDProtocol` controlling equilibration and snapshot collection.

    Returns:
        Dictionary with decorrelated snapshot energies, mean promoter vertical offset,
        and finite-temperature ensemble statistics.
    """
    import numpy as np
    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from pipeline.screening.structure_builders import generate_structure

    temp_default = (
        float(str(genome[4]))
        if len(genome) > 4 and isinstance(genome[4], (int, float, str))
        else 1273.15
    )
    cfg = (
        protocol
        if protocol is not None
        else LiquidMDProtocol(
            temperature_K=temp_default,
            timestep_fs=1.0,
            equilibration_steps=20,
            stride_steps=5,
            snapshots=5,
            friction_per_fs=0.01,
            seed=42,
        )
    )
    cfg.validate()
    if calc is None:
        raise ValueError('run_liquid_alloy_md requires an attached calculator')

    slab, top_idx, _ = generate_structure(genome, structure_seed=cfg.seed)
    slab.pbc = True
    slab.calc = calc
    rng = np.random.default_rng(cfg.seed)
    MaxwellBoltzmannDistribution(
        slab, temperature_K=cfg.temperature_K, rng=rng
    )
    dyn = Langevin(
        slab,
        timestep=cfg.timestep_fs * units.fs,
        temperature_K=cfg.temperature_K,
        friction=cfg.friction_per_fs / units.fs,
        rng=rng,
    )
    if cfg.equilibration_steps > 0:
        dyn.run(cfg.equilibration_steps)

    host_sym = str(genome[1]) if len(genome) > 1 else 'Ga'
    promoter_sym = str(genome[2]) if len(genome) > 2 else 'None'
    snapshots_data: list[dict[str, float | int]] = []
    for snap_idx in range(cfg.snapshots):
        dyn.run(cfg.stride_steps)
        e_pot = float(slab.get_potential_energy())
        e_kin = float(slab.get_kinetic_energy())
        n_free = max(1, len(slab) - len(slab.constraints[0].get_indices()) if slab.constraints else len(slab))
        inst_temp = e_kin / (1.5 * units.kB * float(n_free))
        prom_z = [
            float(slab.positions[i, 2])
            for i in top_idx
            if 0 <= int(i) < len(slab) and slab[int(i)].symbol == promoter_sym
        ]
        host_z = [
            float(slab.positions[i, 2])
            for i in top_idx
            if 0 <= int(i) < len(slab) and slab[int(i)].symbol == host_sym
        ]
        dz_prom = (
            float(np.mean(prom_z) - np.mean(host_z))
            if (prom_z and host_z)
            else 0.0
        )
        snapshots_data.append(
            {
                'snapshot_index': int(snap_idx),
                'potential_energy_eV': e_pot,
                'kinetic_energy_eV': e_kin,
                'instantaneous_temperature_K': float(inst_temp),
                'promoter_vertical_offset_A': dz_prom,
            }
        )

    pot_energies = [float(s['potential_energy_eV']) for s in snapshots_data]
    dz_values = [float(s['promoter_vertical_offset_A']) for s in snapshots_data]
    return {
        'genome': str(genome),
        'temperature_K': float(cfg.temperature_K),
        'equilibration_steps': int(cfg.equilibration_steps),
        'stride_steps': int(cfg.stride_steps),
        'n_snapshots': len(snapshots_data),
        'mean_potential_energy_eV': float(np.mean(pot_energies)),
        'std_potential_energy_eV': float(np.std(pot_energies)),
        'mean_promoter_vertical_offset_A': float(np.mean(dz_values)),
        'snapshots': snapshots_data,
        'complete': len(snapshots_data) == cfg.snapshots,
    }
