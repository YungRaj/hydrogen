"""Independent fuel-cell screening, PEMFC, and stack stage."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from typing import Any, Protocol, cast

from pipeline.data_models.fuel_cells import EmptyStackResult, PEMFCResult, StackResult
from pipeline.data_models.stages import FuelCellProducts, FuelCellState, TableLike
from pipeline.fuel_cell import (
    composite_score as fuel_cell_composite_score,  # pyright: ignore[reportUnknownVariableType]
)
from pipeline.stages.contracts import StageOutcome


class CathodeTable(Protocol):
    """Small dataframe surface required by the fuel-cell coordinator."""

    columns: Any

    def __len__(self) -> int: ...
    def __getitem__(self, key: Any) -> Any: ...
    def copy(self) -> "CathodeTable":
        """Return an isolated table preserving the cathode schema.

        Returns:
            An independent cathode table with the same columns.
        """
        ...

    def nsmallest(self, n: int, column: str) -> "CathodeTable":
        """Return the lowest rows ranked by one numeric column.

        Args:
            n: Maximum number of rows to retain.
            column: Numeric column used for ascending selection.

        Returns:
            A cathode table containing the selected rows.
        """
        ...

    def head(self, n: int) -> "CathodeTable":
        """Return the first cathode records.

        Args:
            n: Maximum number of rows to retain.

        Returns:
            A cathode table containing the leading rows.
        """
        ...

    def iterrows(self) -> Iterator[tuple[Any, Any]]:
        """Yield index and row pairs for selected cathodes.

        Returns:
            An iterator of ``(index, row)`` pairs.
        """
        ...


@dataclass(frozen=True)
class FuelCellServices:
    """Bundle replaceable fuel-cell-stage operations.

    Attributes:
        screen_cathodes: Configured screen cathodes value.
        sweep_membranes: Configured sweep membranes value.
        build_stack_config: Configured build stack config value.
        model_stack: Configured model stack value.
    """

    screen_cathodes: Callable[[], Any]
    sweep_membranes: Callable[..., list[PEMFCResult]]
    build_stack_config: Callable[..., Any]
    model_stack: Callable[..., StackResult]


def default_fuel_cell_services() -> FuelCellServices:
    """Construct production fuel-cell-stage dependencies.

    Returns:
        A `FuelCellServices` containing the default fuel cell services result.
    """
    from pipeline.fuel_cell.pemfc import (
        sweep_membranes,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.fuel_cell.stack import (
        StackConfig,  # pyright: ignore[reportUnknownVariableType]
        model_stack,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.screening.fc_cathode_screener import (
        run_cathode_screening,  # pyright: ignore[reportUnknownVariableType]
    )

    return FuelCellServices(
        cast(Callable[[], Any], run_cathode_screening),
        cast(Callable[..., list[PEMFCResult]], sweep_membranes),
        cast(Callable[..., Any], StackConfig),
        cast(Callable[..., StackResult], model_stack),
    )


def _has_mock_rows(cathodes: Any) -> bool:
    """Return whether any screening row is flagged as mock evidence."""
    columns = getattr(cathodes, 'columns', None)
    if columns is not None:
        if 'mock' not in columns:
            return False
        return bool(cathodes['mock'].fillna(False).astype(bool).any())
    rows: list[Any] = list(cathodes) if cathodes is not None else []
    return any(bool(cast(dict[str, Any], row).get('mock')) for row in rows)


def run_fuel_cell_stage(
    *,
    top_k_pemfc: int,
    stack_cells: int,
    services: FuelCellServices | None = None,
    allow_mock_inputs: bool = False,
) -> StageOutcome[FuelCellState, FuelCellProducts]:
    """Run the complete fuel-cell phase through independently replaceable tools.

    Args:
        top_k_pemfc: Bound controlling top k pemfc.
        stack_cells: Stack cells used by this operation.
        services: Services used by this operation.
        allow_mock_inputs: Explicit test-only permission to consume cathode
            screening rows flagged ``mock`` (e.g. the random fallback used when
            the atomistic calculator is unavailable).

    Returns:
        Computed `StageOutcome` result.

    Raises:
        RuntimeError: If the screening table contains mock rows and
            ``allow_mock_inputs`` is False.
    """
    if top_k_pemfc < 0 or stack_cells <= 0:
        raise ValueError('fuel-cell limits must be physically positive')
    services = services or default_fuel_cell_services()
    cathodes: CathodeTable = services.screen_cathodes()
    if not allow_mock_inputs and _has_mock_rows(cathodes):
        raise RuntimeError(
            'cathode screening returned mock rows; refusing to model PEMFC/stack '
            'performance from them (mock fallback is disabled)'
        )
    from pipeline.search.scope import (
        scope_pemfc_pool,  # pyright: ignore[reportUnknownVariableType]
    )

    cathode_pool, _ = cast(tuple[CathodeTable, Any], scope_pemfc_pool(cathodes))
    valid: CathodeTable = cathode_pool[cathode_pool['valid'] == True].copy()
    top: CathodeTable = (
        valid.nsmallest(top_k_pemfc, 'orr_overpotential_V')
        if 'orr_overpotential_V' in valid.columns
        else valid.head(top_k_pemfc)
    )
    pemfc: list[PEMFCResult] = []
    for _, row in top.iterrows():
        pemfc.extend(
            services.sweep_membranes(
                row['name'],
                # 0.40 V default ORR overpotential (Pt/C benchmark baseline)
                row.get('orr_overpotential_V', 0.4),
                material_class=row.get('material_class', None),
            )
        )
    stack: StackResult | EmptyStackResult = {}
    score_fn = cast(Callable[[PEMFCResult], float], fuel_cell_composite_score)
    best: PEMFCResult | None = None
    if pemfc:
        best = max(pemfc, key=score_fn)
        # Nominal automotive PEMFC rated point: 0.65 V/cell at 1.5 A/cm^2 (DOE 80 kW stack baseline)
        stack = services.model_stack(
            services.build_stack_config(
                n_cells=stack_cells,
                cell_voltage_V=best.get('peak_voltage_V', 0.65),
                current_density_A_cm2=best.get('peak_current_A_cm2', 1.5),
            )
        )
    state: FuelCellState = {
        'n_cathodes_screened': len(cathodes),
        'n_valid': len(valid),
        'n_pemfc_simulations': len(pemfc),
    }
    if 'pourbaix_acid_stable' in valid.columns:
        state['pourbaix_acid_stable_count'] = int(
            valid['pourbaix_acid_stable'].fillna(False).astype(bool).sum()
        )
    if 'orr_overpotential_gc_V' in valid.columns and len(valid) > 0:
        gc_series = valid['orr_overpotential_gc_V'].dropna()
        if len(gc_series) > 0:
            state['best_overpotential_gc_V'] = float(gc_series.min())
    if pemfc and best is not None:
        min_eta = float(min(row.get('orr_overpotential_V', 1.0) for row in pemfc))
        state.update(
            {
                'best_power_W_cm2': float(
                    max(row.get('peak_power_W_cm2', 0.0) for row in pemfc)
                ),
                'best_efficiency': float(
                    max(row.get('efficiency_at_peak', 0.0) for row in pemfc)
                ),
                'min_overpotential_V': min_eta,
                'best_overpotential_V': float(
                    best.get('orr_overpotential_V', min_eta)
                ),
            }
        )
        if best.get('cathode_catalyst'):
            state['best_catalyst'] = str(best['cathode_catalyst'])
        if best.get('membrane'):
            state['best_membrane'] = str(best['membrane'])
        stack_map = cast(dict[str, Any], stack)
        if stack_map.get('net_power_kW') is not None:
            state['stack_net_kW'] = float(stack_map['net_power_kW'])
        if stack_map.get('system_efficiency') is not None:
            state['stack_efficiency'] = float(stack_map['system_efficiency'])
    products: FuelCellProducts = {
        'cathode_database': cast(TableLike, cathodes),
        'valid_cathodes': cast(TableLike, valid),
        'pemfc_results': pemfc,
        'stack_result': stack,
    }
    return StageOutcome(state=state, products=products)
