from __future__ import annotations

from .geometry import build_square_neighbors, edge_key
from .model import InitialConditionKind, LesionKind, SimulationConfig


def build_initial_voltages(config: SimulationConfig, active_nodes: set[int]) -> list[float]:
    node_count = config.lattice_size * config.lattice_size
    voltages = [0.0 for _ in range(node_count)]
    if config.initial_condition != InitialConditionKind.ANTERIOR_POSTERIOR:
        return voltages

    for node in range(node_count):
        if node not in active_nodes:
            continue
        row, col = divmod(node, config.lattice_size)
        coordinate = col if config.polarity_axis == "x" else row
        if coordinate < config.lattice_size / 2:
            voltages[node] = config.posterior_voltage
        else:
            voltages[node] = config.anterior_voltage
    return voltages


def build_polarity_target(config: SimulationConfig, active_nodes: set[int]) -> list[float]:
    node_count = config.lattice_size * config.lattice_size
    target = [0.0 for _ in range(node_count)]
    if not config.polarity_field.enabled:
        return target

    for node in range(node_count):
        if node not in active_nodes:
            continue
        row, col = divmod(node, config.lattice_size)
        coordinate = col if config.polarity_axis == "x" else row
        if coordinate < config.lattice_size / 2:
            target[node] = -config.polarity_field.amplitude
        else:
            target[node] = config.polarity_field.amplitude
    return target


def initial_state(
    config: SimulationConfig,
) -> tuple[
    list[float],
    list[float],
    list[float],
    dict[tuple[int, int], float],
    set[int],
    set[int],
    dict[int, tuple[int, ...]],
]:
    node_count = config.lattice_size * config.lattice_size
    neighbors = build_square_neighbors(config.lattice_size)
    active_nodes = set(range(node_count))
    sink_nodes: set[int] = set()

    if config.lesion.kind == LesionKind.SINK:
        sink_nodes = set(config.lesion.nodes)
    elif config.lesion.kind == LesionKind.TOPOLOGICAL_CUT:
        active_nodes.difference_update(config.lesion.nodes)

    severed = {edge_key(i, j) for i, j in config.lesion.edges}
    conductances: dict[tuple[int, int], float] = {}
    updated_neighbors: dict[int, tuple[int, ...]] = {}

    for site, raw_neighbors in neighbors.items():
        filtered: list[int] = []
        for other in raw_neighbors:
            edge = edge_key(site, other)
            removed_by_cut = (
                config.lesion.kind == LesionKind.TOPOLOGICAL_CUT
                and (site not in active_nodes or other not in active_nodes)
            )
            severed_edge = config.lesion.kind == LesionKind.EDGE_SEVER and edge in severed
            if removed_by_cut or severed_edge:
                continue
            filtered.append(other)
            if edge not in conductances:
                conductances[edge] = config.base_conductance
        updated_neighbors[site] = tuple(filtered)

    voltages = build_initial_voltages(config, active_nodes)
    polarity_target = build_polarity_target(config, active_nodes)
    polarity_field = list(polarity_target)
    return voltages, polarity_field, polarity_target, conductances, active_nodes, sink_nodes, updated_neighbors
