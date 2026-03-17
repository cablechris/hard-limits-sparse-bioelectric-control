from __future__ import annotations

import random

from .geometry import edge_key
from .model import DisorderKind, SimulationConfig


def apply_disorder(
    config: SimulationConfig,
    voltages: list[float],
    conductances: dict[tuple[int, int], float],
    rng: random.Random,
) -> tuple[list[float], dict[tuple[int, int], float]]:
    site_offsets = [0.0 for _ in voltages]
    edge_scales = {edge: 1.0 for edge in conductances}
    spec = config.disorder
    if spec.kind == DisorderKind.NONE or spec.strength == 0.0:
        return site_offsets, edge_scales
    for site in range(len(site_offsets)):
        if spec.kind in (DisorderKind.SITE, DisorderKind.RANDOM_FIELD):
            site_offsets[site] = rng.gauss(0.0, spec.strength)
    for edge in list(edge_scales):
        if spec.kind in (DisorderKind.BOND, DisorderKind.RANDOM_BOND):
            edge_scales[edge] = max(0.0, 1.0 + rng.gauss(0.0, spec.strength))
        if spec.kind == DisorderKind.DILUTION and rng.random() < spec.strength:
            edge_scales[edge] = 0.0
    return site_offsets, edge_scales


def apply_quenched_edge_removals(
    conductances: dict[tuple[int, int], float],
    neighbors: dict[int, tuple[int, ...]],
    edge_scales: dict[tuple[int, int], float],
) -> tuple[dict[tuple[int, int], float], dict[int, tuple[int, ...]]]:
    removed_edges = {edge for edge, scale in edge_scales.items() if scale <= 0.0}
    if not removed_edges:
        return conductances, neighbors

    updated_conductances = dict(conductances)
    for edge in removed_edges:
        updated_conductances.pop(edge, None)

    updated_neighbors: dict[int, tuple[int, ...]] = {}
    for site, adjacent in neighbors.items():
        filtered = tuple(other for other in adjacent if edge_key(site, other) not in removed_edges)
        updated_neighbors[site] = filtered
    return updated_conductances, updated_neighbors
