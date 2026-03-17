from research.control_gap_junction.analysis import (
    control_energy,
    control_operator_matrix,
    default_sparse_actuator_nodes,
    linearized_control_operator,
    linearized_voltage_jacobian,
    project_control_to_sites,
)
from research.control_gap_junction.model import ActuatorFamily, ControlMode, ControlSpec, SimulationConfig
from research.control_gap_junction.simulator import linearization_point_from_state, linearized_pair_from_state, run_simulation


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


def test_linearized_control_operator_rescales_by_tau_v() -> None:
    config = SimulationConfig(
        lattice_size=2,
        tau_v=2.0,
        control=ControlSpec(
            mode=ControlMode.ADDITIVE_VOLTAGE,
            actuator_family=ActuatorFamily.SPARSE_SITE,
            actuator_nodes=(0, 3),
        ),
    )
    operator = linearized_control_operator(config)
    assert operator == [
        [0.5, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.5],
    ]


def test_linearized_control_operator_matches_full_site_identity_over_tau_v() -> None:
    config = SimulationConfig(
        lattice_size=1,
        tau_v=4.0,
        control=ControlSpec(
            mode=ControlMode.ADDITIVE_VOLTAGE,
            actuator_family=ActuatorFamily.FULL_SITE,
        ),
    )
    assert linearized_control_operator(config) == [[0.25]]


def test_control_energy_matches_quadratic_time_sum() -> None:
    trace = [
        [1.0, 0.0],
        [0.0, 2.0],
    ]
    assert control_energy(trace, dt=0.5) == 2.5


def test_project_control_to_sites_maps_sparse_channels_into_site_forcing() -> None:
    operator = [
        [1.0, 0.0],
        [0.0, 0.0],
        [0.0, 1.0],
    ]
    assert project_control_to_sites(operator, [2.0, -1.5]) == [2.0, 0.0, -1.5]


def test_run_simulation_applies_open_loop_control_trace() -> None:
    config = SimulationConfig(
        lattice_size=2,
        dt=1.0,
        steps=1,
        noise_amplitude=0.0,
        base_conductance=0.0,
        well_linear=0.0,
        well_cubic=1.0,
        stimulation_amplitude=0.0,
        injury_amplitude=0.0,
        control=ControlSpec(
            mode=ControlMode.ADDITIVE_VOLTAGE,
            actuator_family=ActuatorFamily.SPARSE_SITE,
            actuator_nodes=(0,),
            control_trace=((2.0,),),
        ),
    )
    config.polarity_field.enabled = False

    history = run_simulation(config)

    assert history[0].site_control == [2.0, 0.0, 0.0, 0.0]
    assert history[0].control_channels == [2.0]
    assert history[0].voltages[0] == 2.0


def test_linearization_point_from_state_freezes_snapshot_data() -> None:
    config = SimulationConfig(lattice_size=2, steps=1, noise_amplitude=0.0)
    history = run_simulation(config)

    point = linearization_point_from_state(history[0])

    assert tuple(history[0].voltages) == point.voltages
    assert tuple(history[0].polarity_field) == point.polarity_field
    assert tuple(sorted(history[0].conductances.items())) == point.conductances


def test_linearized_voltage_jacobian_matches_decoupled_on_site_derivative() -> None:
    config = SimulationConfig(
        lattice_size=2,
        steps=1,
        noise_amplitude=0.0,
        base_conductance=0.0,
        well_linear=1.2,
        well_cubic=0.5,
    )
    config.polarity_field.enabled = False
    config.plasticity.adaptation_strength = 0.0
    history = run_simulation(config)
    point = linearization_point_from_state(history[0])

    jacobian = linearized_voltage_jacobian(config, point, history[0].neighbors)

    expected_diagonal = [
        (config.well_linear - 3.0 * config.well_cubic * voltage * voltage) / config.tau_v
        for voltage in point.voltages
    ]
    for index, expected in enumerate(expected_diagonal):
        assert abs(jacobian[index][index] - expected) < 1e-5
        off_diagonal_sum = sum(abs(jacobian[index][j]) for j in range(len(point.voltages)) if j != index)
        assert off_diagonal_sum < 1e-6


def test_linearized_pair_from_state_packages_point_a_and_b() -> None:
    config = SimulationConfig(
        lattice_size=2,
        tau_v=2.0,
        steps=1,
        noise_amplitude=0.0,
    )
    history = run_simulation(config)

    pair = linearized_pair_from_state(config, history[0])

    assert pair.point.voltages == tuple(history[0].voltages)
    assert len(pair.a_matrix) == 4
    assert all(len(row) == 4 for row in pair.a_matrix)
    assert len(pair.b_matrix) == 4
    assert all(len(row) == 2 for row in pair.b_matrix)
