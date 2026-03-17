from __future__ import annotations

from dataclasses import dataclass

from .model import MacrostateKind, SimulationConfig


@dataclass(slots=True)
class MorphologyFeatures:
    anterior_mean: float
    posterior_mean: float
    polarity_contrast: float
    dh_propensity: float


@dataclass(slots=True)
class MorphologyAssessment:
    macrostate: MacrostateKind
    features: MorphologyFeatures
    note: str


def side_means(voltages: list[float], lattice_size: int, axis: str = "x") -> tuple[float, float]:
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
        return 0.0, 0.0
    return sum(low_side) / len(low_side), sum(high_side) / len(high_side)


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def extract_morphology_features(voltages: list[float], config: SimulationConfig) -> MorphologyFeatures:
    posterior_mean, anterior_mean = side_means(
        voltages,
        lattice_size=config.lattice_size,
        axis=config.polarity_axis,
    )
    polarity_contrast = anterior_mean - posterior_mean
    reference_span = max(config.anterior_voltage - config.posterior_voltage, 1e-9)
    dh_propensity = clamp01((posterior_mean - config.posterior_voltage) / reference_span)
    return MorphologyFeatures(
        anterior_mean=anterior_mean,
        posterior_mean=posterior_mean,
        polarity_contrast=polarity_contrast,
        dh_propensity=dh_propensity,
    )


def classify_macrostate(
    voltages: list[float],
    config: SimulationConfig,
    wt_cutoff: float = 0.125,
    dh_cutoff: float = 0.875,
) -> MorphologyAssessment:
    features = extract_morphology_features(voltages, config)
    if features.dh_propensity <= wt_cutoff:
        macrostate = MacrostateKind.WT
    elif features.dh_propensity >= dh_cutoff:
        macrostate = MacrostateKind.DH
    else:
        macrostate = MacrostateKind.CRYPTIC
    return MorphologyAssessment(
        macrostate=macrostate,
        features=features,
        note=(
            "Three-state readout is currently a coarse-grained surrogate based on posterior "
            "depolarization propensity, not a full three-well dynamical realization."
        ),
    )
