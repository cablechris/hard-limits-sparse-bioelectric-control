from research.control_gap_junction.analysis import (
    finite_horizon_controllability_report,
    finite_horizon_controllability_gramian,
    finite_horizon_controllability_gramian_symmetric,
    linear_minimum_energy,
    linear_minimum_energy_symmetric_projected,
    linear_minimum_energy_symmetric,
    linear_minimum_energy_symmetric_reachable,
    linear_minimum_energy_control_samples,
    mode_overlap_matrix,
    symmetric_modal_control_weights,
    symmetric_eigensystem,
)


def test_finite_horizon_gramian_matches_identity_case() -> None:
    gramian = finite_horizon_controllability_gramian(
        [[0.0]],
        [[1.0]],
        horizon=2.5,
        quadrature_steps=64,
    )
    assert abs(gramian[0][0] - 2.5) < 1e-10


def test_finite_horizon_symmetric_gramian_matches_identity_case() -> None:
    gramian = finite_horizon_controllability_gramian_symmetric(
        [[0.0]],
        [[1.0]],
        horizon=2.5,
    )
    assert abs(gramian[0][0] - 2.5) < 1e-12


def test_linear_minimum_energy_matches_identity_case() -> None:
    energy = linear_minimum_energy(
        [[0.0]],
        [[1.0]],
        delta_v0=[3.0],
        delta_v_target=[0.0],
        horizon=2.0,
        quadrature_steps=64,
    )
    assert abs(energy - 4.5) < 1e-10


def test_linear_minimum_energy_symmetric_matches_identity_case() -> None:
    energy = linear_minimum_energy_symmetric(
        [[0.0]],
        [[1.0]],
        delta_v0=[3.0],
        delta_v_target=[0.0],
        horizon=2.0,
    )
    assert abs(energy - 4.5) < 1e-12


def test_linear_minimum_energy_symmetric_reachable_detects_unreachable_target() -> None:
    energy = linear_minimum_energy_symmetric_reachable(
        [[0.0, 0.0], [0.0, 0.0]],
        [[1.0], [0.0]],
        delta_v0=[0.0, 0.0],
        delta_v_target=[0.0, 1.0],
        horizon=1.0,
        tolerance=1e-12,
    )
    assert energy == float("inf")


def test_linear_minimum_energy_symmetric_projected_reports_unreachable_residual() -> None:
    report = linear_minimum_energy_symmetric_projected(
        [[0.0, 0.0], [0.0, 0.0]],
        [[1.0], [0.0]],
        delta_v0=[0.0, 0.0],
        delta_v_target=[0.0, 1.0],
        horizon=1.0,
        tolerance=1e-12,
    )
    assert report.reachable_rank == 1
    assert report.reachable_energy == 0.0
    assert abs(report.unreachable_residual_norm - 1.0) < 1e-12


def test_linear_minimum_energy_control_samples_are_constant_for_identity_case() -> None:
    controls = linear_minimum_energy_control_samples(
        [[0.0]],
        [[1.0]],
        delta_v0=[3.0],
        delta_v_target=[0.0],
        horizon=2.0,
        sample_times=[0.0, 1.0, 2.0],
        quadrature_steps=64,
    )
    for control in controls:
        assert abs(control[0] + 1.5) < 1e-10


def test_symmetric_eigensystem_and_mode_overlap_match_diagonal_case() -> None:
    eigenvalues, eigenvectors = symmetric_eigensystem([[2.0, 0.0], [0.0, 1.0]], descending=True)
    assert eigenvalues == [2.0, 1.0]
    overlaps = mode_overlap_matrix(eigenvectors, eigenvectors)
    assert abs(overlaps[0][0] - 1.0) < 1e-12
    assert abs(overlaps[1][1] - 1.0) < 1e-12


def test_symmetric_modal_control_weights_match_identity_case() -> None:
    eigenvalues, weights = symmetric_modal_control_weights([[0.0]], [[1.0]], horizon=2.5)
    assert eigenvalues == [0.0]
    assert abs(weights[0] - 2.5) < 1e-12


def test_finite_horizon_controllability_report_matches_identity_case() -> None:
    report = finite_horizon_controllability_report(
        [[0.0]],
        [[1.0]],
        horizon=2.5,
        quadrature_steps=64,
        tolerance=1e-12,
    )
    assert report.state_dimension == 1
    assert report.input_dimension == 1
    assert report.numerical_rank == 1
    assert report.controllable_at_tolerance is True
    assert abs(report.smallest_positive_eigenvalue - 2.5) < 1e-10
    assert abs(report.reachable_condition_number - 1.0) < 1e-12
