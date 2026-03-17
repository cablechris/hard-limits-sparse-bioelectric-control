from __future__ import annotations


def final_abs_mean_voltage(voltages: list[float]) -> float:
    return abs(sum(voltages) / len(voltages))


def ap_voltage_variance(
    voltages: list[float],
    lattice_size: int,
    axis: str = "x",
) -> tuple[float, float, float, float]:
    low_side: list[float] = []
    high_side: list[float] = []
    for node, value in enumerate(voltages):
        row, col = divmod(node, lattice_size)
        coordinate = col if axis == "x" else row
        if coordinate < lattice_size / 2:
            low_side.append(value)
        else:
            high_side.append(value)
    if not low_side or not high_side:
        return 0.0, 0.0, 0.0, 0.0
    low_mean = sum(low_side) / len(low_side)
    high_mean = sum(high_side) / len(high_side)
    combined = low_side + high_side
    mean_value = sum(combined) / len(combined)
    variance = sum((value - mean_value) ** 2 for value in combined) / len(combined)
    return high_mean, low_mean, high_mean - low_mean, variance


def spatial_variance(voltages: list[float]) -> float:
    if not voltages:
        return 0.0
    mean_value = sum(voltages) / len(voltages)
    return sum((value - mean_value) ** 2 for value in voltages) / len(voltages)


def consensus_fraction(voltages: list[float]) -> float:
    positive = sum(1 for value in voltages if value >= 0.0)
    negative = len(voltages) - positive
    return max(positive, negative) / len(voltages)


def polarity_contrast(voltages: list[float], lattice_size: int, axis: str = "x") -> float:
    low_side: list[float] = []
    high_side: list[float] = []
    for node, value in enumerate(voltages):
        row, col = divmod(node, lattice_size)
        coordinate = col if axis == "x" else row
        if coordinate < lattice_size / 2:
            low_side.append(value)
        else:
            high_side.append(value)
    if not low_side or not high_side:
        return 0.0
    return (sum(high_side) / len(high_side)) - (sum(low_side) / len(low_side))
