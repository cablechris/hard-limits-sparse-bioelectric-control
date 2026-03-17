from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class PotentialKind(StrEnum):
    # WT / Cryptic / DH support will require an explicit three-well or hybrid
    # slow-variable extension here. The current double-well model is only a
    # voltage-surrogate foundation.
    DOUBLE_WELL = "double_well"


class InitialConditionKind(StrEnum):
    ANTERIOR_POSTERIOR = "anterior_posterior"


class ControlMode(StrEnum):
    NONE = "none"
    ADDITIVE_VOLTAGE = "additive_voltage"


class ActuatorFamily(StrEnum):
    FULL_SITE = "full_site"
    SPARSE_SITE = "sparse_site"
    MODE_RESTRICTED = "mode_restricted"


class MacrostateKind(StrEnum):
    WT = "wt"
    CRYPTIC = "cryptic"
    DH = "dh"


class DisorderKind(StrEnum):
    NONE = "none"
    SITE = "site"
    BOND = "bond"
    RANDOM_FIELD = "random_field"
    RANDOM_BOND = "random_bond"
    DILUTION = "dilution"


class LesionKind(StrEnum):
    NONE = "none"
    SINK = "sink"
    EDGE_SEVER = "edge_sever"
    TOPOLOGICAL_CUT = "topological_cut"


@dataclass(slots=True)
class DisorderSpec:
    kind: DisorderKind = DisorderKind.NONE
    distribution: str = "none"
    strength: float = 0.0
    seed: int = 0


@dataclass(slots=True)
class LesionSpec:
    kind: LesionKind = LesionKind.NONE
    nodes: tuple[int, ...] = ()
    edges: tuple[tuple[int, int], ...] = ()
    clamp_voltage: float = 0.0


@dataclass(slots=True)
class PlasticitySpec:
    adaptation_strength: float = 0.1
    mismatch_target: float = 0.0
    mismatch_slope: float = 1.0
    decay_rate: float = 0.05
    field_scale: float = 0.0


@dataclass(slots=True)
class PolarityFieldSpec:
    enabled: bool = True
    amplitude: float = 0.15
    tau_h: float = 2_000.0


@dataclass(slots=True)
class ControlSpec:
    mode: ControlMode = ControlMode.ADDITIVE_VOLTAGE
    actuator_family: ActuatorFamily = ActuatorFamily.SPARSE_SITE
    actuator_nodes: tuple[int, ...] = ()
    mode_vectors: tuple[tuple[float, ...], ...] = ()
    control_trace: tuple[tuple[float, ...], ...] = ()
    energy_weight: float = 1.0
    amplitude_limit: float | None = None


@dataclass(slots=True)
class LinearizationPoint:
    voltages: tuple[float, ...]
    conductances: tuple[tuple[tuple[int, int], float], ...]
    polarity_field: tuple[float, ...]
    note: str = "Frozen damaged pre-intervention snapshot for local voltage-sector linearization."


@dataclass(slots=True)
class LinearizedPair:
    point: LinearizationPoint
    a_matrix: list[list[float]]
    b_matrix: list[list[float]]


@dataclass(slots=True)
class SimulationConfig:
    lattice_size: int = 16
    tau_v: float = 1.0
    tau_g: float = 20.0
    potential_kind: PotentialKind = PotentialKind.DOUBLE_WELL
    initial_condition: InitialConditionKind = InitialConditionKind.ANTERIOR_POSTERIOR
    polarity_axis: str = "x"
    anterior_voltage: float = 1.0
    posterior_voltage: float = -1.0
    well_linear: float = 1.0
    well_cubic: float = 1.0
    base_conductance: float = 0.35
    dt: float = 0.02
    steps: int = 200
    noise_amplitude: float = 0.02
    stimulation_amplitude: float = 0.0
    injury_amplitude: float = 0.0
    disorder: DisorderSpec = field(default_factory=DisorderSpec)
    lesion: LesionSpec = field(default_factory=LesionSpec)
    plasticity: PlasticitySpec = field(default_factory=PlasticitySpec)
    polarity_field: PolarityFieldSpec = field(default_factory=PolarityFieldSpec)
    control: ControlSpec = field(default_factory=ControlSpec)

    @property
    def epsilon(self) -> float:
        return self.tau_v / self.tau_g
