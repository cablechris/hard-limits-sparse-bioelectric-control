from research.control_gap_junction.placement import (
    cellwise_best_energy_map,
    central_patch_displacement,
    corner_patch_displacement,
    exhaustive_two_site_placement_energies,
    greedy_sparse_placement_sequence,
    placement_projected_report,
    side_patch_displacement,
)


def test_side_patch_displacement_marks_requested_side_columns() -> None:
    values = side_patch_displacement(4, width=1, amplitude=2.0, axis="x", side="low")
    assert values[0] == 2.0
    assert values[4] == 2.0
    assert values[1] == 0.0


def test_corner_patch_displacement_marks_requested_corner() -> None:
    values = corner_patch_displacement(4, width=1, amplitude=2.0, corner="low_high")
    assert values[3] == 2.0
    assert values[0] == 0.0


def test_central_patch_displacement_marks_center_square() -> None:
    values = central_patch_displacement(4, width=2, amplitude=1.0)
    assert values[5] == 1.0
    assert values[6] == 1.0
    assert values[0] == 0.0


def test_exhaustive_two_site_placement_energies_counts_all_pairs_on_small_lattice() -> None:
    results = exhaustive_two_site_placement_energies(
        [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ],
        lattice_size=2,
        tau_v=1.0,
        horizon=1.0,
        delta_v0=[0.0, 0.0, 0.0, 0.0],
        candidate_nodes=[0, 1, 2, 3],
    )
    assert len(results) == 6


def test_cellwise_best_energy_map_uses_best_pair_for_each_cell() -> None:
    energy_map = cellwise_best_energy_map(
        2,
        [
            {"nodes": (0, 1), "energy": 3.0},
            {"nodes": (0, 2), "energy": 2.0},
            {"nodes": (1, 3), "energy": 1.0},
        ],
    )
    assert energy_map == [
        [2.0, 1.0],
        [2.0, 1.0],
    ]


def test_placement_projected_report_matches_identity_case() -> None:
    report = placement_projected_report(
        [[0.0]],
        lattice_size=1,
        tau_v=1.0,
        horizon=2.0,
        delta_v0=[3.0],
        nodes=(0,),
    )
    assert abs(report.reachable_energy - 4.5) < 1e-12
    assert abs(report.unreachable_residual_norm) < 1e-12
    assert report.reachable_rank == 1


def test_greedy_sparse_placement_sequence_adds_requested_number_of_sites() -> None:
    sequence = greedy_sparse_placement_sequence(
        [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ],
        lattice_size=2,
        tau_v=1.0,
        horizon=1.0,
        delta_v0=[1.0, 0.0, 0.0, 0.0],
        max_count=2,
    )
    assert len(sequence) == 2
    assert len(sequence[-1]["nodes"]) == 2
