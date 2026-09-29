"""Independent fuel-cell screening, PEMFC, and stack stage."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Protocol

from pipeline.data_models.stages import FuelCellProducts, FuelCellState
from pipeline.data_models.fuel_cells import (
    EmptyStackResult, PEMFCResult, StackResult)
from pipeline.stages.contracts import StageOutcome


class CathodeTable(Protocol):
    """Small dataframe surface required by the fuel-cell coordinator."""

    columns: object

    def __len__(self) -> int: ...
    def __getitem__(self, key): ...
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

    def iterrows(self):
        """Yield index and row pairs for selected cathodes."""
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
    screen_cathodes: Callable
    sweep_membranes: Callable
    build_stack_config: Callable
    model_stack: Callable


def default_fuel_cell_services() -> FuelCellServices:
    """Construct production fuel-cell-stage dependencies.

    Returns:
        A `FuelCellServices` containing the default fuel cell services result.
    """
    from pipeline.screening.fc_cathode_screener import run_cathode_screening
    from pipeline.fuel_cell.pemfc import sweep_membranes
    from pipeline.fuel_cell.stack import StackConfig, model_stack
    return FuelCellServices(
        run_cathode_screening, sweep_membranes, StackConfig, model_stack)


def fuel_cell_composite_score(result: PEMFCResult) -> float:
    """Preserve the established efficiency/power/overpotential priority.

    Args:
        result: Mapping supplying result.

    Returns:
        Computed `float` value in the units documented above.
    """
    efficiency = result.get('efficiency_at_peak', 0.0)
    power = result.get('peak_power_W_cm2', 0.0)
    overpotential = max(result.get('orr_overpotential_V', 0.4), 0.01)
    return (efficiency * power) / overpotential


def run_fuel_cell_stage(*, top_k_pemfc: int, stack_cells: int,
                        services: FuelCellServices | None = None
                        ) -> StageOutcome[FuelCellState, FuelCellProducts]:
    """Run the complete fuel-cell phase through independently replaceable tools.

    Args:
        top_k_pemfc: Bound controlling top k pemfc.
        stack_cells: Stack cells used by this operation.
        services: Services used by this operation.

    Returns:
        Computed `StageOutcome` result.
    """
    if top_k_pemfc < 0 or stack_cells <= 0:
        raise ValueError('fuel-cell limits must be physically positive')
    services = services or default_fuel_cell_services()
    cathodes = services.screen_cathodes()
    from pipeline.search.scope import scope_pemfc_pool
    cathode_pool, _ = scope_pemfc_pool(cathodes)
    valid = cathode_pool[cathode_pool['valid'] == True].copy()
    top = valid.nsmallest(top_k_pemfc, 'orr_overpotential_V') \
        if 'orr_overpotential_V' in valid.columns else valid.head(top_k_pemfc)
    pemfc = []
    for _, row in top.iterrows():
        pemfc.extend(services.sweep_membranes(
            row['name'], row.get('orr_overpotential_V', 0.4),
            material_class=row.get('material_class', None)))
    stack: StackResult | EmptyStackResult = {}
    if pemfc:
        best = max(pemfc, key=fuel_cell_composite_score)
        stack = services.model_stack(services.build_stack_config(
            n_cells=stack_cells,
            cell_voltage_V=best.get('peak_voltage_V', 0.65),
            current_density_A_cm2=best.get('peak_current_A_cm2', 1.5)))
    state: FuelCellState = {
        'n_cathodes_screened': len(cathodes), 'n_valid': len(valid),
        'n_pemfc_simulations': len(pemfc)}
    if pemfc:
        state.update({
            'best_power_W_cm2': max(
                row.get('peak_power_W_cm2', 0) for row in pemfc),
            'best_efficiency': max(
                row.get('efficiency_at_peak', 0) for row in pemfc),
            'min_overpotential_V': min(
                row.get('orr_overpotential_V', 1.0) for row in pemfc),
        })
    products: FuelCellProducts = {
        'cathode_database': cathodes, 'valid_cathodes': valid,
        'pemfc_results': pemfc, 'stack_result': stack}
    return StageOutcome(state=state, products=products)
