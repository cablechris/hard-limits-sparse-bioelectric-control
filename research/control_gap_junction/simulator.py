from __future__ import annotations

from dataclasses import dataclass
import random

from .disorder import apply_disorder, apply_quenched_edge_removals
from .dynamics import step
from .initialization import initial_state
from .model import SimulationConfig


@dataclass(slots=True)
class SimulationState:
    voltages: list[float]
    polarity_field: list[float]
    conductances: dict[tuple[int, int], float]
    active_nodes: set[int]
    sink_nodes: set[int]
    neighbors: dict[int, tuple[int, ...]]


def run_simulation(config: SimulationConfig) -> list[SimulationState]:
    rng = random.Random(config.disorder.seed)
    voltages, polarity_field, polarity_target, conductances, active_nodes, sink_nodes, neighbors = initial_state(config)
    site_offsets, edge_scales = apply_disorder(config, voltages, conductances, rng)
    conductances, neighbors = apply_quenched_edge_removals(conductances, neighbors, edge_scales)
    history: list[SimulationState] = []
    for _ in range(config.steps):
        voltages, polarity_field, conductances = step(
            config=config,
            voltages=voltages,
            polarity_field=polarity_field,
            polarity_target=polarity_target,
            conductances=conductances,
            active_nodes=active_nodes,
            sink_nodes=sink_nodes,
            neighbors=neighbors,
            site_offsets=site_offsets,
            edge_scales=edge_scales,
            rng=rng,
        )
        history.append(
            SimulationState(
                voltages=list(voltages),
                polarity_field=list(polarity_field),
                conductances=dict(conductances),
                active_nodes=set(active_nodes),
                sink_nodes=set(sink_nodes),
                neighbors=dict(neighbors),
            )
        )
    return history
