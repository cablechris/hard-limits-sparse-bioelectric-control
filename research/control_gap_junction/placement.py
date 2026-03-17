from __future__ import annotations

from collections.abc import Iterable
from itertools import combinations

from .analysis import (
    ProjectedSteeringReport,
    linear_minimum_energy_symmetric_projected,
    linear_minimum_energy_symmetric_reachable,
    linearized_control_operator,
)
from .model import ActuatorFamily, ControlMode, ControlSpec, SimulationConfig


def side_patch_displacement(
    lattice_size: int,
    *,
    width: int = 2,
    amplitude: float = 1.0,
    axis: str = "x",
    side: str = "low",
) -> list[float]:
    """Return a uniform rectangular depolarization patch on one side of the lattice."""
    if lattice_size <= 0:
        raise ValueError("lattice_size must be positive")
    if width <= 0 or width > lattice_size:
        raise ValueError("width must lie in [1, lattice_size]")
    if side not in {"low", "high"}:
        raise ValueError("side must be 'low' or 'high'")

    values = [0.0 for _ in range(lattice_size * lattice_size)]
    for node in range(len(values)):
        row, col = divmod(node, lattice_size)
        coordinate = col if axis == "x" else row
        selected = coordinate < width if side == "low" else coordinate >= lattice_size - width
        if selected:
            values[node] = amplitude
    return values


def corner_patch_displacement(
    lattice_size: int,
    *,
    width: int = 2,
    amplitude: float = 1.0,
    corner: str = "low_low",
) -> list[float]:
    """Return a square depolarization patch anchored at a chosen corner."""
    if corner not in {"low_low", "low_high", "high_low", "high_high"}:
        raise ValueError("corner must be one of low_low, low_high, high_low, high_high")
    if width <= 0 or width > lattice_size:
        raise ValueError("width must lie in [1, lattice_size]")

    values = [0.0 for _ in range(lattice_size * lattice_size)]
    for node in range(len(values)):
        row, col = divmod(node, lattice_size)
        row_selected = row < width if corner.startswith("low") else row >= lattice_size - width
        col_selected = col < width if corner.endswith("low") else col >= lattice_size - width
        if row_selected and col_selected:
            values[node] = amplitude
    return values


def central_patch_displacement(
    lattice_size: int,
    *,
    width: int = 2,
    amplitude: float = 1.0,
) -> list[float]:
    """Return a centered square depolarization patch."""
    if width <= 0 or width > lattice_size:
        raise ValueError("width must lie in [1, lattice_size]")

    start = max(0, (lattice_size - width) // 2)
    stop = start + width
    values = [0.0 for _ in range(lattice_size * lattice_size)]
    for node in range(len(values)):
        row, col = divmod(node, lattice_size)
        if start <= row < stop and start <= col < stop:
            values[node] = amplitude
    return values


def two_site_sparse_b_matrix(lattice_size: int, tau_v: float, nodes: tuple[int, int]) -> list[list[float]]:
    """Return the linearized sparse two-site actuator matrix for the chosen nodes."""
    return sparse_b_matrix(lattice_size, tau_v, nodes)


def sparse_b_matrix(lattice_size: int, tau_v: float, nodes: tuple[int, ...]) -> list[list[float]]:
    """Return the linearized sparse actuator matrix for an arbitrary node set."""
    config = SimulationConfig(
        lattice_size=lattice_size,
        tau_v=tau_v,
        control=ControlSpec(
            mode=ControlMode.ADDITIVE_VOLTAGE,
            actuator_family=ActuatorFamily.SPARSE_SITE,
            actuator_nodes=nodes,
        ),
    )
    return linearized_control_operator(config)


def exhaustive_two_site_placement_energies(
    a_matrix: list[list[float]],
    *,
    lattice_size: int,
    tau_v: float,
    horizon: float,
    delta_v0: list[float],
    delta_v_target: list[float] | None = None,
    candidate_nodes: Iterable[int] | None = None,
    tolerance: float = 1e-12,
) -> list[dict[str, float | tuple[int, int]]]:
    """Evaluate minimum energy for every two-site sparse actuator placement."""
    nodes = tuple(candidate_nodes) if candidate_nodes is not None else tuple(range(lattice_size * lattice_size))
    target = delta_v_target or [0.0 for _ in delta_v0]
    results: list[dict[str, float | tuple[int, int]]] = []
    for pair in combinations(nodes, 2):
        b_matrix = two_site_sparse_b_matrix(lattice_size, tau_v, pair)
        energy = linear_minimum_energy_symmetric_reachable(
            a_matrix,
            b_matrix,
            delta_v0=delta_v0,
            delta_v_target=target,
            horizon=horizon,
            tolerance=tolerance,
        )
        results.append({"nodes": pair, "energy": energy})
    return results


def placement_projected_report(
    a_matrix: list[list[float]],
    *,
    lattice_size: int,
    tau_v: float,
    horizon: float,
    delta_v0: list[float],
    nodes: tuple[int, ...],
    delta_v_target: list[float] | None = None,
    tolerance: float = 1e-12,
) -> ProjectedSteeringReport:
    """Return projected residual-energy metrics for a sparse actuator set."""
    target = delta_v_target or [0.0 for _ in delta_v0]
    b_matrix = sparse_b_matrix(lattice_size, tau_v, nodes)
    return linear_minimum_energy_symmetric_projected(
        a_matrix,
        b_matrix,
        delta_v0=delta_v0,
        delta_v_target=target,
        horizon=horizon,
        tolerance=tolerance,
    )


def greedy_sparse_placement_sequence(
    a_matrix: list[list[float]],
    *,
    lattice_size: int,
    tau_v: float,
    horizon: float,
    delta_v0: list[float],
    max_count: int,
    initial_nodes: tuple[int, ...] = (),
    delta_v_target: list[float] | None = None,
    tolerance: float = 1e-12,
) -> list[dict[str, object]]:
    """Greedily add actuator sites minimizing projected residual then energy."""
    if max_count <= 0:
        raise ValueError("max_count must be positive")

    chosen = tuple(sorted(initial_nodes))
    node_count = lattice_size * lattice_size
    if len(set(chosen)) != len(chosen):
        raise ValueError("initial_nodes must be unique")

    sequence: list[dict[str, object]] = []
    if chosen:
        report = placement_projected_report(
            a_matrix,
            lattice_size=lattice_size,
            tau_v=tau_v,
            horizon=horizon,
            delta_v0=delta_v0,
            nodes=chosen,
            delta_v_target=delta_v_target,
            tolerance=tolerance,
        )
        sequence.append(
            {
                "nodes": chosen,
                "projected_energy": report.reachable_energy,
                "residual_norm": report.unreachable_residual_norm,
                "reachable_rank": report.reachable_rank,
            }
        )

    while len(chosen) < max_count:
        best_candidate: dict[str, object] | None = None
        for node in range(node_count):
            if node in chosen:
                continue
            candidate_nodes = tuple(sorted(chosen + (node,)))
            report = placement_projected_report(
                a_matrix,
                lattice_size=lattice_size,
                tau_v=tau_v,
                horizon=horizon,
                delta_v0=delta_v0,
                nodes=candidate_nodes,
                delta_v_target=delta_v_target,
                tolerance=tolerance,
            )
            record = {
                "nodes": candidate_nodes,
                "added_node": node,
                "projected_energy": report.reachable_energy,
                "residual_norm": report.unreachable_residual_norm,
                "reachable_rank": report.reachable_rank,
            }
            if best_candidate is None:
                best_candidate = record
                continue
            if (
                record["residual_norm"],
                record["projected_energy"],
                -record["reachable_rank"],
                record["added_node"],
            ) < (
                best_candidate["residual_norm"],
                best_candidate["projected_energy"],
                -best_candidate["reachable_rank"],
                best_candidate["added_node"],
            ):
                best_candidate = record

        assert best_candidate is not None
        chosen = tuple(best_candidate["nodes"])
        sequence.append(best_candidate)
    return sequence


def cellwise_best_energy_map(
    lattice_size: int,
    placement_results: list[dict[str, float | tuple[int, int]]],
) -> list[list[float]]:
    """Return the best achievable energy when each cell is forced to be one actuator site."""
    node_count = lattice_size * lattice_size
    best = [float("inf") for _ in range(node_count)]
    for result in placement_results:
        pair = result["nodes"]
        energy = float(result["energy"])
        assert isinstance(pair, tuple)
        for node in pair:
            best[node] = min(best[node], energy)
    return [
        best[row * lattice_size:(row + 1) * lattice_size]
        for row in range(lattice_size)
    ]
