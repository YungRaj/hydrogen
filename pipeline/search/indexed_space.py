"""O(1), deterministic search addressing of the encoded catalyst space.

No population is materialized.  A global integer maps to exactly one encoded
configuration, enabling disjoint workers, restarts, and auditable coverage.
The size intentionally matches ``estimate_design_space_size``; physically
invalid Cartesian entries are emitted and must be recorded as rejected rather
than silently disappearing from the denominator.
"""

from __future__ import annotations

from math import comb
from typing import Iterator, Sequence, Tuple
import heapq

import numpy as np

from pipeline.common import catalyst_spaces as cs


def _decode(index: int, dimensions: Sequence[Sequence]) -> list:
    values = [None] * len(dimensions)
    for pos in range(len(dimensions) - 1, -1, -1):
        dim = dimensions[pos]
        index, digit = divmod(index, len(dim))
        values[pos] = dim[digit]
    if index:
        raise IndexError("mixed-radix index outside space")
    return values


def _unrank_combination(n: int, k: int, rank: int) -> Tuple[int, ...]:
    """Lexicographically unrank one k-combination without enumeration."""
    if rank < 0 or rank >= comb(n, k):
        raise IndexError("combination rank outside space")
    result, start = [], 0
    for remaining in range(k, 0, -1):
        for value in range(start, n):
            block = comb(n - value - 1, remaining - 1)
            if rank < block:
                result.append(value)
                start = value + 1
                break
            rank -= block
    return tuple(result)


CLASS_SIZES = cs.estimate_design_space_size().copy()
CLASS_SIZES.pop("TOTAL")
CLASS_ORDER = tuple(cs.ALL_MATERIAL_CLASSES)
CLASS_OFFSETS = {}
_offset = 0
for _cls in CLASS_ORDER:
    CLASS_OFFSETS[_cls] = _offset
    _offset += CLASS_SIZES[_cls]
TOTAL_SIZE = _offset
SOLID_STRAINS = tuple(float(x) for x in np.linspace(
    cs.SOLID_STRAIN_RANGE[0], cs.SOLID_STRAIN_RANGE[1], 20
))


def candidate_at_class(material_class: str, index: int) -> tuple:
    """Return the candidate at a class-local index in O(number of genes)."""
    size = CLASS_SIZES[material_class]
    if index < 0 or index >= size:
        raise IndexError(f"{material_class} index {index} outside [0, {size})")

    if material_class == "MoltenMetal":
        return (material_class, *_decode(index, [cs.MOLTEN_HOSTS, cs.MOLTEN_PROMOTERS,
            cs.MOLTEN_PROMOTER_AT_PCT, cs.MOLTEN_TEMPERATURES_K]))
    if material_class == "SolidCatalyst":
        metal, support, facet, strain, d1, d2, nsub, vac = _decode(index, [
            cs.SOLID_ACTIVE_METALS, cs.SOLID_SUPPORTS, cs.SOLID_FACETS, SOLID_STRAINS,
            cs.SOLID_DOPANTS, cs.SOLID_DOPANTS, range(1, 5), range(3)])
        return material_class, metal, support, facet, strain, (d1, d2), nsub, vac
    if material_class == "SAC":
        return (material_class, *_decode(index, [cs.SAC_METALS, cs.SAC_COORDINATIONS,
            cs.SAC_SUBSTRATES, cs.SAC_AXIAL_LIGANDS]))
    if material_class == "DAC":
        return (material_class, *_decode(index, [cs.DAC_METALS_1, cs.DAC_METALS_2,
            cs.DAC_COORDINATIONS, cs.SAC_SUBSTRATES]))
    if material_class == "MOF":
        return (material_class, *_decode(index, [cs.MOF_METAL_NODES, cs.MOF_LINKERS,
            cs.MOF_CAVITIES, cs.MOF_PORE_SIZES]))
    if material_class == "COF":
        return (material_class, *_decode(index, [cs.MOF_METAL_NODES + ['None'], cs.COF_LINKAGES,
            cs.MOF_CAVITIES, cs.MOF_PORE_SIZES]))
    if material_class == "Perovskite":
        return (material_class, *_decode(index, [cs.PEROVSKITE_A_SITE, cs.PEROVSKITE_B_SITE,
            cs.PEROVSKITE_B_SITE + ['None'], cs.PEROVSKITE_DOPANT_FRAC, cs.PEROVSKITE_DEFECTS]))
    if material_class == "MetalHydride":
        temps = [300, 350, 400, 450, 500, 550, 600, 700, 800]
        return (material_class, *_decode(index, [cs.HYDRIDE_METALS, cs.HYDRIDE_TYPES,
            cs.HYDRIDE_SECOND_METAL, cs.HYDRIDE_ADDITIVES, temps]))
    if material_class == "MAXPhase":
        return (material_class, *_decode(index, [cs.MAX_M_ELEMENTS, cs.MAX_A_ELEMENTS,
            cs.MAX_X_ELEMENTS, cs.MAX_N_VALUES, cs.MAX_M_ELEMENTS + ['None'],
            ['basal_0001', 'edge_1010', 'edge_1120']]))
    if material_class == "HEA":
        tail_size = len(cs.HEA_STRUCTURES) * 4 * 6
        combo_rank, tail = divmod(index, tail_size)
        for k in (4, 5, 6):
            count = comb(len(cs.HEA_ELEMENTS), k)
            if combo_rank < count:
                components = tuple(cs.HEA_ELEMENTS[i] for i in _unrank_combination(len(cs.HEA_ELEMENTS), k, combo_rank))
                structure, facet, temp = _decode(tail, [cs.HEA_STRUCTURES, ['111', '100', '110', '211'],
                                                       [800, 900, 1000, 1100, 1200, 1300]])
                return material_class, components, structure, facet, temp
            combo_rank -= count
    if material_class == "Spinel":
        return (material_class, *_decode(index, [cs.SPINEL_A_METALS, cs.SPINEL_B_METALS,
            cs.SPINEL_DOPANTS, cs.SPINEL_MORPHOLOGIES, cs.SPINEL_SUPPORT_CARBONS]))
    if material_class == "MXene":
        return (material_class, *_decode(index, [cs.MXENE_M_ELEMENTS, cs.MXENE_X_ELEMENTS,
            cs.MXENE_N_VALUES, cs.MXENE_TERMINATIONS, cs.MXENE_SAC_METALS]))
    if material_class == "SAA":
        return (material_class, *_decode(index, [cs.SAA_TRACE_METALS, cs.SAA_HOST_METALS,
            cs.SAA_FACETS, cs.SAA_LOADINGS_PPM]))
    if material_class == "MetalFreeCarbon":
        return (material_class, *_decode(index, [cs.MFC_N_TYPES, cs.MFC_N_FRACTIONS,
            cs.MFC_DEFECT_TYPES, cs.MFC_SUBSTRATES, cs.MFC_DOPANTS]))
    raise KeyError(material_class)


def candidate_at(global_index: int) -> tuple:
    """Map a global index in ``[0, TOTAL_SIZE)`` to a candidate."""
    if global_index < 0 or global_index >= TOTAL_SIZE:
        raise IndexError(f"global index {global_index} outside [0, {TOTAL_SIZE})")
    # Fourteen entries: a linear lookup is faster and clearer than allocating a
    # 21-billion-entry map.
    for cls in reversed(CLASS_ORDER):
        if global_index >= CLASS_OFFSETS[cls]:
            return candidate_at_class(cls, global_index - CLASS_OFFSETS[cls])
    raise AssertionError("unreachable")


def iter_shard(start: int, stop: int, worker_id: int = 0,
               num_workers: int = 1) -> Iterator[Tuple[int, tuple]]:
    """Yield one deterministic, disjoint strided worker shard."""
    if not 0 <= start <= stop <= TOTAL_SIZE:
        raise ValueError("invalid shard bounds")
    if num_workers <= 0 or not 0 <= worker_id < num_workers:
        raise ValueError("invalid worker assignment")
    for index in range(start + worker_id, stop, num_workers):
        yield index, candidate_at(index)


def is_physically_admissible(genome: tuple) -> Tuple[bool, str]:
    """Conservative, cheap rejection rules; uncertainty is never rejected."""
    cls = genome[0]
    if cls == "MoltenMetal" and not cs.validate_molten_metal(genome):
        return False, "molten_temperature_or_composition"
    if cls == "Perovskite":
        _, _, base, dopant, frac, _ = genome
        if dopant == 'None' and frac != 0.0:
            return False, "dopant_fraction_without_dopant"
        if dopant == base and frac > 0:
            return False, "dopant_equals_base_site"
    if cls == "MAXPhase" and genome[5] == genome[1]:
        return False, "dopant_equals_host"
    if cls == "SAA" and genome[1] == genome[2]:
        return False, "trace_equals_host"
    return True, "accepted"


def deterministic_tree_probes(count: int) -> list:
    """Return calibration points by recursively bisecting all class ranges.

    These points calibrate the surrogate; they are not used to claim population
    coverage. Population discovery remains the exhaustive branch scanner.
    """
    if count <= 0:
        return []
    # Visit every class root before recursively refining the largest intervals.
    queues = {cls: [] for cls in CLASS_ORDER}
    probes, seen = [], set()

    def visit(cls, start, stop):
        mid = start + (stop - start) // 2
        # Search deterministically outward when the midpoint is hard-invalid.
        chosen = None
        for delta in range(min(stop - start, 1024)):
            for index in (mid + delta, mid - delta):
                if start <= index < stop and index not in seen:
                    genome = candidate_at(index)
                    if is_physically_admissible(genome)[0]:
                        chosen = (index, genome)
                        break
            if chosen:
                break
        if chosen:
            seen.add(chosen[0]); probes.append(chosen[1])
        if mid > start:
            heapq.heappush(queues[cls], (-(mid - start), start, mid))
        if stop > mid + 1:
            heapq.heappush(queues[cls], (-(stop - mid - 1), mid + 1, stop))
    for cls in CLASS_ORDER:
        if len(probes) >= count:
            break
        start = CLASS_OFFSETS[cls]
        visit(cls, start, start + CLASS_SIZES[cls])
    # Round-robin refinement prevents the 20.9B SolidCatalyst class from
    # monopolizing calibration evidence solely because it is largest.
    while len(probes) < count and any(queues.values()):
        for cls in CLASS_ORDER:
            if len(probes) >= count:
                break
            if queues[cls]:
                _, start, stop = heapq.heappop(queues[cls])
                visit(cls, start, stop)
    return probes
