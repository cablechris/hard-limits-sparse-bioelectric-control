from research.control_gap_junction.analysis import (
    control_energy,
    control_operator_matrix,
    default_sparse_actuator_nodes,
)
from research.control_gap_junction.model import ActuatorFamily, ControlMode, ControlSpec, SimulationConfig


def test_default_sparse_actuator_nodes_returns_two_sites_for_standard_lattice() -> None:
    config = SimulationConfig(lattice_size=16)
    nodes = default_sparse_actuator_nodes(config)
    assert len(nodes) == 2
    assert nodes[0] != nodes[1]


def test_control_operator_matrix_full_site_is_identity() -> None:
    config = SimulationConfig(
        lattice_size=2,
        control=ControlSpec(mode=ControlMode.ADDITIVE_VOLTAGE, actuator_family=ActuatorFamily.FULL_SITE),
    )
    operator = control_operator_matrix(config)
    assert operator == [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]


def test_control_operator_matrix_sparse_site_uses_selected_nodes() -> None:
    config = SimulationConfig(
        lattice_size=2,
        control=ControlSpec(
            mode=ControlMode.ADDITIVE_VOLTAGE,
            actuator_family=ActuatorFamily.SPARSE_SITE,
            actuator_nodes=(0, 3),
        ),
    )
    operator = control_operator_matrix(config)
    assert operator == [
        [1.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 1.0],
    ]


def test_control_energy_matches_quadratic_time_sum() -> None:
    trace = [
        [1.0, 0.0],
        [0.0, 2.0],
    ]
    assert control_energy(trace, dt=0.5) == 2.5
