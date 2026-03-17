from __future__ import annotations

from dataclasses import dataclass
import random

from .analysis import (
    control_operator_matrix,
    linearized_control_operator,
    linearized_voltage_jacobian,
    project_control_to_sites,
)
from .disorder import apply_disorder, apply_quenched_edge_removals
from .dynamics import step
from .initialization import initial_state
from .model import LinearizationPoint, LinearizedPair, SimulationConfig


@dataclass(slots=True)
class SimulationState:
    voltages: list[float]
    polarity_field: list[float]
    conductances: dict[tuple[int, int], float]
    active_nodes: set[int]
    sink_nodes: set[int]
    neighbors: dict[int, tuple[int, ...]]
    control_channels: list[float]
    site_control: list[float]


def linearization_point_from_state(state: SimulationState) -> LinearizationPoint:
    """Freeze a simulation snapshot into the default Phase 2 operating point."""
    return LinearizationPoint(
        voltages=tuple(state.voltages),
        conductances=tuple(sorted(state.conductances.items())),
        polarity_field=tuple(state.polarity_field),
    )


def linearized_pair_from_state(
    config: SimulationConfig,
    state: SimulationState,
    *,
    site_offsets: list[float] | None = None,
    edge_scales: dict[tuple[int, int], float] | None = None,
    site_control: list[float] | None = None,
    epsilon: float = 1e-6,
) -> LinearizedPair:
    """Package the default Phase 1 local linearized pair from one frozen state."""
    point = linearization_point_from_state(state)
    return LinearizedPair(
        point=point,
        a_matrix=linearized_voltage_jacobian(
            config,
            point,
            state.neighbors,
            site_offsets=site_offsets,
            edge_scales=edge_scales,
            site_control=site_control,
            epsilon=epsilon,
        ),
        b_matrix=linearized_control_operator(config),
    )


def _control_channels_at_step(config: SimulationConfig, step_index: int, channel_count: int) -> list[float]:
    if channel_count == 0 or step_index >= len(config.control.control_trace):
        return [0.0 for _ in range(channel_count)]

    control = list(config.control.control_trace[step_index])
    if len(control) != channel_count:
        raise ValueError("control_trace entry width must match the actuator channel count")
    return control


def run_simulation(config: SimulationConfig) -> list[SimulationState]:
    rng = random.Random(config.disorder.seed)
    voltages, polarity_field, polarity_target, conductances, active_nodes, sink_nodes, neighbors = initial_state(config)
    site_offsets, edge_scales = apply_disorder(config, voltages, conductances, rng)
    conductances, neighbors = apply_quenched_edge_removals(conductances, neighbors, edge_scales)
    operator = control_operator_matrix(config)
    channel_count = len(operator[0]) if operator else 0
    history: list[SimulationState] = []
    for step_index in range(config.steps):
        control_channels = _control_channels_at_step(config, step_index, channel_count)
        site_control = (
            project_control_to_sites(operator, control_channels)
            if operator
            else [0.0 for _ in voltages]
        )
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
            site_control=site_control,
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
                control_channels=list(control_channels),
                site_control=list(site_control),
            )
        )
    return history
