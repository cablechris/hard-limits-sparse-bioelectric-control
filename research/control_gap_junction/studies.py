from __future__ import annotations

from dataclasses import dataclass

from .placement import (
    central_patch_displacement,
    corner_patch_displacement,
    greedy_sparse_placement_sequence,
    side_patch_displacement,
)


@dataclass(slots=True)
class DamagePatternSpec:
    name: str
    displacement: list[float]


def default_damage_patterns_8x8() -> list[DamagePatternSpec]:
    return [
        DamagePatternSpec(
            name="side_patch",
            displacement=side_patch_displacement(8, width=2, amplitude=1.0, axis="x", side="low"),
        ),
        DamagePatternSpec(
            name="corner_wound",
            displacement=corner_patch_displacement(8, width=2, amplitude=1.0, corner="low_low"),
        ),
        DamagePatternSpec(
            name="central_lesion",
            displacement=central_patch_displacement(8, width=2, amplitude=1.0),
        ),
    ]


def scaled_damage_patterns_16x16() -> list[DamagePatternSpec]:
    return [
        DamagePatternSpec(
            name="side_patch",
            displacement=side_patch_displacement(16, width=4, amplitude=1.0, axis="x", side="low"),
        ),
        DamagePatternSpec(
            name="corner_wound",
            displacement=corner_patch_displacement(16, width=4, amplitude=1.0, corner="low_low"),
        ),
        DamagePatternSpec(
            name="central_lesion",
            displacement=central_patch_displacement(16, width=4, amplitude=1.0),
        ),
    ]


def residual_threshold_hit(
    sequence: list[dict[str, object]],
    threshold: float,
) -> int | None:
    for record in sequence:
        if float(record["residual_norm"]) <= threshold:
            return len(record["nodes"])
    return None


def greedy_rank_trajectory(
    a_matrix: list[list[float]],
    *,
    lattice_size: int,
    tau_v: float,
    horizon: float,
    displacement: list[float],
    max_count: int = 16,
    tolerance: float = 1e-12,
) -> list[dict[str, object]]:
    return greedy_sparse_placement_sequence(
        a_matrix,
        lattice_size=lattice_size,
        tau_v=tau_v,
        horizon=horizon,
        delta_v0=displacement,
        max_count=max_count,
        tolerance=tolerance,
    )
