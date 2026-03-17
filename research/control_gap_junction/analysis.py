from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

from .dynamics import deterministic_voltage_drift
from .initialization import build_polarity_target
from .model import ActuatorFamily, ControlMode, MacrostateKind, SimulationConfig
from .morphology import extract_morphology_features


def _copy_matrix(matrix: list[list[float]]) -> list[list[float]]:
    return [list(row) for row in matrix]


def _solve_linear_system(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Solve a dense linear system by Gaussian elimination with pivoting."""
    size = len(matrix)
    if size == 0:
        return []
    augmented = [list(matrix[row]) + [rhs[row]] for row in range(size)]

    for pivot in range(size):
        pivot_row = max(range(pivot, size), key=lambda idx: abs(augmented[idx][pivot]))
        if abs(augmented[pivot_row][pivot]) <= 1e-12:
            raise ValueError("linear system is singular or ill-conditioned at the requested tolerance")
        if pivot_row != pivot:
            augmented[pivot], augmented[pivot_row] = augmented[pivot_row], augmented[pivot]

        pivot_value = augmented[pivot][pivot]
        for column in range(pivot, size + 1):
            augmented[pivot][column] /= pivot_value

        for row in range(size):
            if row == pivot:
                continue
            factor = augmented[row][pivot]
            if factor == 0.0:
                continue
            for column in range(pivot, size + 1):
                augmented[row][column] -= factor * augmented[pivot][column]

    return [augmented[row][size] for row in range(size)]


def _max_symmetric_eigenvalue(matrix: list[list[float]], iterations: int = 64) -> float:
    """Estimate the largest eigenvalue of a symmetric matrix by power iteration."""
    size = len(matrix)
    if size == 0:
        return 0.0
    vector = [1.0 for _ in range(size)]
    norm = math.sqrt(sum(value * value for value in vector))
    vector = [value / norm for value in vector]

    eigenvalue = 0.0
    for _ in range(iterations):
        image = [
            sum(matrix[row][col] * vector[col] for col in range(size))
            for row in range(size)
        ]
        norm = math.sqrt(sum(value * value for value in image))
        if norm <= 1e-15:
            return 0.0
        vector = [value / norm for value in image]
        eigenvalue = sum(
            vector[row] * sum(matrix[row][col] * vector[col] for col in range(size))
            for row in range(size)
        )
    return eigenvalue


def _indices_from_mask(mask: list[bool] | tuple[bool, ...] | list[int] | tuple[int, ...], size: int) -> list[int]:
    if len(mask) == 0:
        return []
    if all(isinstance(entry, bool) for entry in mask):
        bool_mask = [bool(entry) for entry in mask]
        if len(bool_mask) != size:
            raise ValueError("boolean mask length must match the state dimension")
        return [index for index, selected in enumerate(bool_mask) if selected]
    return [int(entry) for entry in mask]


def _dot(vec_a: list[float], vec_b: list[float]) -> float:
    return sum(a * b for a, b in zip(vec_a, vec_b, strict=True))


def _orthonormalize(vectors: list[list[float]], tolerance: float = 1e-12) -> list[list[float]]:
    basis: list[list[float]] = []
    for vector in vectors:
        working = list(vector)
        for basis_vector in basis:
            overlap = _dot(working, basis_vector)
            working = [
                value - overlap * basis_value
                for value, basis_value in zip(working, basis_vector, strict=True)
            ]
        norm = math.sqrt(_dot(working, working))
        if norm > tolerance:
            basis.append([value / norm for value in working])
    return basis


def default_sparse_actuator_nodes(config: SimulationConfig) -> tuple[int, ...]:
    """Return a compact default sparse-site actuator set near the lattice midline."""
    size = config.lattice_size
    if size < 2:
        return (0,)
    center_row = size // 2
    left_col = max(0, size // 2 - 1)
    right_col = min(size - 1, left_col + 1)
    return (
        center_row * size + left_col,
        center_row * size + right_col,
    )


def first_ap_mode_laplacian_eigenvalue(lattice_size: int) -> float:
    """Return the first non-constant Neumann-grid mode along the AP axis.

    For the square lattice with free boundaries, the lowest x-varying mode has
    eigenvalue lambda_1 = 2 * (1 - cos(pi / L)).
    """
    if lattice_size < 2:
        return 0.0
    return 2.0 * (1.0 - math.cos(math.pi / lattice_size))


def effective_smoothing_rate(
    average_conductance: float,
    lattice_size: int,
    bond_occupation: float,
) -> float:
    """Estimate AP-mode Laplacian smoothing under bond dilution.

    This mean-field estimate replaces the diluted graph by an occupied-edge
    fraction p multiplying the intact-grid conductance scale.
    """
    return bond_occupation * average_conductance * first_ap_mode_laplacian_eigenvalue(lattice_size)


def ap_mode_drive_threshold_linearized(
    *,
    average_conductance: float,
    lattice_size: int,
    bond_occupation: float,
    target_amplitude: float,
    well_linear: float,
) -> float:
    """Projected drive amplitude needed to hold the first AP mode.

    Linearizing around small AP-mode amplitude A gives

        tau_V dA/dt ~= (a - p G lambda_1) A + h_1,

    so the drive needed to hold A is

        h_1^* = (p G lambda_1 - a) A.

    Negative values mean the on-site instability already supports the mode in
    the linearized approximation, so the minimum sustaining external drive is
    clipped to zero.
    """
    smoothing = effective_smoothing_rate(
        average_conductance=average_conductance,
        lattice_size=lattice_size,
        bond_occupation=bond_occupation,
    )
    return max(0.0, (smoothing - well_linear) * target_amplitude)


def homogenization_timescale(
    *,
    average_conductance: float,
    lattice_size: int,
    bond_occupation: float,
    tau_v: float,
) -> float:
    """Estimate AP-mode decay time after drive removal.

    Neglecting on-site destabilization and conductance adaptation, the first AP
    mode decays on the Laplacian timescale

        tau_relax = tau_v / (p G lambda_1).

    Infinite timescale means the projected smoothing rate vanishes.
    """
    smoothing = effective_smoothing_rate(
        average_conductance=average_conductance,
        lattice_size=lattice_size,
        bond_occupation=bond_occupation,
    )
    if smoothing <= 0.0:
        return math.inf
    return tau_v / smoothing


def minimum_existence_condition(a: float, average_conductance: float, lambda_max_laplacian: float) -> bool:
    """Return the linear-response invertibility condition for ``-a I + G L``.

    In the present Case 1 convention this checks

        a > G * lambda_max(L).

    This is the condition under which the linear operator used in the driven
    fixed-point estimate is negative definite and therefore invertible. It is a
    local linear-response condition, not a standalone proof of global uniqueness
    for the full quartic double-well functional.
    """
    return a > average_conductance * lambda_max_laplacian


def driven_fixed_point_linear(
    h: list[float],
    a: float,
    average_conductance: float,
    laplacian_matrix: list[list[float]],
) -> list[float]:
    """Solve the linearized driven fixed-point equation.

        (-a I + G_bar L) V_* = h

    The helper emits a warning when ``a <= G_bar * lambda_max(L)`` because that
    falls outside the negative-definite regime of the linear-response operator
    used in the current Case 1 derivation.
    """
    size = len(h)
    if size != len(laplacian_matrix):
        raise ValueError("drive vector and Laplacian matrix must have the same dimension")
    if any(len(row) != size for row in laplacian_matrix):
        raise ValueError("laplacian_matrix must be square")

    lambda_max = _max_symmetric_eigenvalue(laplacian_matrix)
    if not minimum_existence_condition(a, average_conductance, lambda_max):
        warnings.warn(
            "linear-response condition a > G_bar * lambda_max(L) is violated; "
            "the fixed-point estimate may not represent a stable driven minimum",
            stacklevel=2,
        )

    operator = _copy_matrix(laplacian_matrix)
    for row in range(size):
        for column in range(size):
            operator[row][column] *= average_conductance
        operator[row][row] -= a
    return _solve_linear_system(operator, h)


def ap_differential_at_minimum(
    h: list[float],
    a: float,
    average_conductance: float,
    laplacian_matrix: list[list[float]],
    anterior_mask: list[bool] | tuple[bool, ...] | list[int] | tuple[int, ...],
    posterior_mask: list[bool] | tuple[bool, ...] | list[int] | tuple[int, ...],
) -> float:
    """Return the AP differential of the linearized driven fixed point."""
    fixed_point = driven_fixed_point_linear(h, a, average_conductance, laplacian_matrix)
    anterior_indices = _indices_from_mask(anterior_mask, len(fixed_point))
    posterior_indices = _indices_from_mask(posterior_mask, len(fixed_point))
    if not anterior_indices or not posterior_indices:
        raise ValueError("anterior and posterior selections must both be non-empty")

    anterior_mean = sum(fixed_point[index] for index in anterior_indices) / len(anterior_indices)
    posterior_mean = sum(fixed_point[index] for index in posterior_indices) / len(posterior_indices)
    return anterior_mean - posterior_mean


@dataclass(slots=True)
class GenerativeModelSpec:
    macrostate_priors: dict[MacrostateKind, float]
    macrostate_targets: dict[MacrostateKind, tuple[float, float, float, float]]
    feature_sigmas: tuple[float, float, float, float]


def feature_vector_from_voltages(voltages: list[float], config: SimulationConfig) -> tuple[float, float, float, float]:
    """Return the current coarse feature map phi(V)."""
    features = extract_morphology_features(voltages, config)
    return (
        features.anterior_mean,
        features.posterior_mean,
        features.polarity_contrast,
        features.dh_propensity,
    )


def default_generative_model(config: SimulationConfig) -> GenerativeModelSpec:
    """Return the smallest explicit Phase 3 generative family.

    This is intentionally simple:

    - uniform prior over WT/Cryptic/DH
    - diagonal Gaussian likelihood on the current coarse feature map
    - macrostate targets inherited from the present surrogate classifier
    """
    anterior = config.anterior_voltage
    posterior = config.posterior_voltage
    midpoint = 0.5 * (anterior + posterior)
    reference_span = anterior - posterior
    return GenerativeModelSpec(
        macrostate_priors={
            MacrostateKind.WT: 1.0 / 3.0,
            MacrostateKind.CRYPTIC: 1.0 / 3.0,
            MacrostateKind.DH: 1.0 / 3.0,
        },
        macrostate_targets={
            MacrostateKind.WT: (anterior, posterior, reference_span, 0.0),
            MacrostateKind.CRYPTIC: (anterior, midpoint, 0.5 * reference_span, 0.5),
            MacrostateKind.DH: (anterior, anterior, 0.0, 1.0),
        },
        feature_sigmas=(0.25, 0.25, 0.5, 0.2),
    )


def negative_log_joint_from_features(
    feature_vector: tuple[float, float, float, float],
    macrostate: MacrostateKind,
    model: GenerativeModelSpec,
) -> float:
    """Return -log P(m) - log P(phi(V) | m) up to an additive constant."""
    prior = model.macrostate_priors[macrostate]
    if prior <= 0.0:
        raise ValueError("macrostate prior must be positive")
    target = model.macrostate_targets[macrostate]
    quadratic = 0.0
    for observed, center, sigma in zip(feature_vector, target, model.feature_sigmas, strict=True):
        if sigma <= 0.0:
            raise ValueError("feature sigma must be positive")
        quadratic += 0.5 * ((observed - center) / sigma) ** 2
    return -math.log(prior) + quadratic


def negative_log_joint_from_voltages(
    voltages: list[float],
    config: SimulationConfig,
    macrostate: MacrostateKind,
    model: GenerativeModelSpec | None = None,
) -> float:
    """Evaluate the explicit coarse-grained generative objective on voltages."""
    generative_model = model or default_generative_model(config)
    feature_vector = feature_vector_from_voltages(voltages, config)
    return negative_log_joint_from_features(feature_vector, macrostate, generative_model)


def classify_macrostate_from_objective(
    voltages: list[float],
    config: SimulationConfig,
    model: GenerativeModelSpec | None = None,
) -> tuple[MacrostateKind, dict[MacrostateKind, float]]:
    """Classify by the smallest coarse-grained negative log joint objective."""
    generative_model = model or default_generative_model(config)
    scores = {
        macrostate: negative_log_joint_from_voltages(voltages, config, macrostate, generative_model)
        for macrostate in (MacrostateKind.WT, MacrostateKind.CRYPTIC, MacrostateKind.DH)
    }
    best = min(scores, key=scores.get)
    return best, scores


def correspondence_status_base_case() -> str:
    """Return the current Phase 3 classification for the base system."""
    return "no_exact_correspondence_base_case"


def feature_jacobian(voltages: list[float], config: SimulationConfig) -> list[list[float]]:
    """Return the Jacobian of the current coarse feature map phi(V).

    Rows are ordered as:

    1. anterior_mean
    2. posterior_mean
    3. AP_differential
    4. dh_propensity

    The `dh_propensity` derivative is piecewise because the current surrogate
    clamps it into `[0, 1]`. We use zero derivative in the saturated regime and
    the linear posterior-mean derivative in the interior.
    """
    node_count = len(voltages)
    anterior_indices: list[int] = []
    posterior_indices: list[int] = []
    for node in range(node_count):
        row, col = divmod(node, config.lattice_size)
        coordinate = col if config.polarity_axis == "x" else row
        if coordinate < config.lattice_size / 2:
            posterior_indices.append(node)
        else:
            anterior_indices.append(node)

    if not anterior_indices or not posterior_indices:
        raise ValueError("feature_jacobian requires both anterior and posterior partitions")

    anterior_weight = 1.0 / len(anterior_indices)
    posterior_weight = 1.0 / len(posterior_indices)
    reference_span = max(config.anterior_voltage - config.posterior_voltage, 1e-9)
    features = extract_morphology_features(voltages, config)
    dh_interior = 0.0 < features.dh_propensity < 1.0

    jacobian = [[0.0 for _ in range(node_count)] for _ in range(4)]
    for index in anterior_indices:
        jacobian[0][index] = anterior_weight
        jacobian[2][index] = anterior_weight
    for index in posterior_indices:
        jacobian[1][index] = posterior_weight
        jacobian[2][index] = -posterior_weight
        if dh_interior:
            jacobian[3][index] = posterior_weight / reference_span
    return jacobian


def coarse_objective_gradient(
    voltages: list[float],
    config: SimulationConfig,
    macrostate: MacrostateKind,
    model: GenerativeModelSpec | None = None,
) -> list[float]:
    """Pull back the coarse generative objective gradient into voltage space."""
    generative_model = model or default_generative_model(config)
    feature_vector = feature_vector_from_voltages(voltages, config)
    target = generative_model.macrostate_targets[macrostate]
    residual = [
        (observed - center) / (sigma * sigma)
        for observed, center, sigma in zip(feature_vector, target, generative_model.feature_sigmas, strict=True)
    ]
    jacobian = feature_jacobian(voltages, config)
    node_count = len(voltages)
    return [
        sum(jacobian[row][node] * residual[row] for row in range(len(jacobian)))
        for node in range(node_count)
    ]


def voltage_drift_vector(
    config: SimulationConfig,
    voltages: list[float],
    polarity_field: list[float],
    conductances: dict[tuple[int, int], float],
    neighbors: dict[int, tuple[int, ...]],
    site_offsets: list[float],
    edge_scales: dict[tuple[int, int], float],
) -> list[float]:
    """Return the deterministic voltage drift vector for the current state."""
    return [
        deterministic_voltage_drift(
            config,
            site=site,
            voltages=voltages,
            polarity_field=polarity_field,
            conductances=conductances,
            neighbors=neighbors,
            site_offsets=site_offsets,
            edge_scales=edge_scales,
        )
        for site in range(len(voltages))
    ]


def cosine_similarity(vec_a: list[float], vec_b: list[float]) -> float:
    """Return cosine similarity between two vectors, or 0 for zero-norm input."""
    norm_a = math.sqrt(sum(value * value for value in vec_a))
    norm_b = math.sqrt(sum(value * value for value in vec_b))
    if norm_a <= 1e-15 or norm_b <= 1e-15:
        return 0.0
    return sum(a * b for a, b in zip(vec_a, vec_b, strict=True)) / (norm_a * norm_b)


def drift_variance_in_feature_subspace(drift_v: list[float], jacobian: list[list[float]]) -> float:
    """Return the fraction of drift variance captured by the feature subspace.

    The relevant subspace is the column space of ``J_phi^T``, equivalently the
    row space of ``J_phi``. We compute the orthogonal projection directly using
    an orthonormal basis of the row space so the result remains well-defined
    even when feature rows are linearly dependent.
    """
    drift_norm_sq = _dot(drift_v, drift_v)
    if drift_norm_sq <= 1e-15:
        return 0.0
    basis = _orthonormalize(jacobian)
    if not basis:
        return 0.0
    projected = [0.0 for _ in drift_v]
    for basis_vector in basis:
        weight = _dot(drift_v, basis_vector)
        projected = [
            value + weight * basis_value
            for value, basis_value in zip(projected, basis_vector, strict=True)
        ]
    projected_norm_sq = _dot(projected, projected)
    return projected_norm_sq / drift_norm_sq


def weighted_laplacian_matrix(
    node_count: int,
    conductances: dict[tuple[int, int], float],
) -> list[list[float]]:
    """Build the weighted graph Laplacian from edge conductances."""
    matrix = [[0.0 for _ in range(node_count)] for _ in range(node_count)]
    for (i, j), weight in conductances.items():
        matrix[i][i] += weight
        matrix[j][j] += weight
        matrix[i][j] -= weight
        matrix[j][i] -= weight
    return matrix


def control_operator_matrix(config: SimulationConfig) -> list[list[float]]:
    """Build the additive actuator matrix B for the current control specification."""
    node_count = config.lattice_size * config.lattice_size
    if config.control.mode == ControlMode.NONE:
        return []
    if config.control.mode != ControlMode.ADDITIVE_VOLTAGE:
        raise ValueError(f"Unsupported control mode: {config.control.mode}")

    if config.control.actuator_family == ActuatorFamily.FULL_SITE:
        return [
            [1.0 if row == col else 0.0 for col in range(node_count)]
            for row in range(node_count)
        ]

    if config.control.actuator_family == ActuatorFamily.SPARSE_SITE:
        nodes = config.control.actuator_nodes or default_sparse_actuator_nodes(config)
        columns: list[list[float]] = []
        for node in nodes:
            column = [0.0 for _ in range(node_count)]
            column[node] = 1.0
            columns.append(column)
        return [
            [column[row] for column in columns]
            for row in range(node_count)
        ]

    if config.control.actuator_family == ActuatorFamily.MODE_RESTRICTED:
        if not config.control.mode_vectors:
            raise ValueError("mode_restricted control requires explicit mode_vectors")
        if any(len(vector) != node_count for vector in config.control.mode_vectors):
            raise ValueError("every mode vector must have length equal to the state dimension")
        return [
            [vector[row] for vector in config.control.mode_vectors]
            for row in range(node_count)
        ]

    raise ValueError(f"Unsupported actuator family: {config.control.actuator_family}")


def control_energy(control_trace: list[list[float]], dt: float, energy_weight: float = 1.0) -> float:
    """Quadratic intervention energy integral approximated by a time sum."""
    total = 0.0
    for control in control_trace:
        total += sum(value * value for value in control)
    return energy_weight * dt * total


def static_ap_drive_vector(config: SimulationConfig, amplitude: float | None = None) -> list[float]:
    """Return the static AP drive vector used in the Case 1 comparison.

    By default this uses the configured polarity-field amplitude and target
    geometry but treats the drive as static rather than dynamically relaxing.
    """
    drive_config = SimulationConfig(
        lattice_size=config.lattice_size,
        polarity_axis=config.polarity_axis,
        anterior_voltage=config.anterior_voltage,
        posterior_voltage=config.posterior_voltage,
        polarity_field=config.polarity_field,
    )
    target = build_polarity_target(drive_config, set(range(config.lattice_size * config.lattice_size)))
    if amplitude is None or abs(config.polarity_field.amplitude) <= 1e-15:
        return target
    scale = amplitude / config.polarity_field.amplitude
    return [scale * value for value in target]
