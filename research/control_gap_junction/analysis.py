from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

from .dynamics import deterministic_voltage_drift
from .initialization import build_polarity_target
from .model import ActuatorFamily, ControlMode, LinearizationPoint, MacrostateKind, SimulationConfig
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


@dataclass(slots=True)
class ControllabilityReport:
    state_dimension: int
    input_dimension: int
    horizon: float
    numerical_rank: int
    tolerance: float
    min_eigenvalue: float
    max_eigenvalue: float
    smallest_positive_eigenvalue: float
    reachable_condition_number: float
    controllable_at_tolerance: bool


@dataclass(slots=True)
class ProjectedSteeringReport:
    reachable_energy: float
    unreachable_residual_norm: float
    reachable_rank: int


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
    site_control: list[float] | None = None,
) -> list[float]:
    """Return the deterministic voltage drift vector for the current state."""
    local_site_control = site_control or [0.0 for _ in voltages]
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
            site_control=local_site_control,
        )
        for site in range(len(voltages))
    ]


def frozen_conductances_from_point(point: LinearizationPoint) -> dict[tuple[int, int], float]:
    """Return the frozen conductance field stored in a linearization point."""
    return dict(point.conductances)


def linearized_voltage_jacobian(
    config: SimulationConfig,
    point: LinearizationPoint,
    neighbors: dict[int, tuple[int, ...]],
    *,
    site_offsets: list[float] | None = None,
    edge_scales: dict[tuple[int, int], float] | None = None,
    site_control: list[float] | None = None,
    epsilon: float = 1e-6,
) -> list[list[float]]:
    """Return a finite-difference Jacobian for the frozen voltage-sector drift.

    This computes the local matrix

        A ~= d/dV [(1 / tau_V) F(V, G_bar)] |_(V_bar)

    around the supplied frozen operating point. Conductances are held fixed at
    ``G_bar`` from ``point``; no control-law, Gramian, or finite-horizon solver
    is introduced here.
    """
    if epsilon <= 0.0:
        raise ValueError("epsilon must be positive")

    voltages = list(point.voltages)
    polarity_field = list(point.polarity_field)
    conductances = frozen_conductances_from_point(point)
    node_count = len(voltages)

    if len(polarity_field) != node_count:
        raise ValueError("linearization point polarity field must match the voltage dimension")

    local_site_offsets = site_offsets or [0.0 for _ in range(node_count)]
    local_edge_scales = edge_scales or {edge: 1.0 for edge in conductances}
    local_site_control = site_control or [0.0 for _ in range(node_count)]

    if len(local_site_offsets) != node_count:
        raise ValueError("site_offsets must match the voltage dimension")
    if len(local_site_control) != node_count:
        raise ValueError("site_control must match the voltage dimension")

    jacobian = [[0.0 for _ in range(node_count)] for _ in range(node_count)]
    for column in range(node_count):
        plus_voltages = list(voltages)
        minus_voltages = list(voltages)
        plus_voltages[column] += epsilon
        minus_voltages[column] -= epsilon

        drift_plus = voltage_drift_vector(
            config,
            plus_voltages,
            polarity_field,
            conductances,
            neighbors,
            local_site_offsets,
            local_edge_scales,
            local_site_control,
        )
        drift_minus = voltage_drift_vector(
            config,
            minus_voltages,
            polarity_field,
            conductances,
            neighbors,
            local_site_offsets,
            local_edge_scales,
            local_site_control,
        )
        for row in range(node_count):
            jacobian[row][column] = (drift_plus[row] - drift_minus[row]) / (2.0 * epsilon)
    return jacobian


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


def linearized_control_operator(config: SimulationConfig) -> list[list[float]]:
    """Return the linearized additive control operator B_tilde = B / tau_v."""
    if abs(config.tau_v) <= 1e-15:
        raise ValueError("tau_v must be nonzero to define the linearized control operator")
    operator = control_operator_matrix(config)
    return [
        [value / config.tau_v for value in row]
        for row in operator
    ]


def symmetric_eigensystem(matrix: list[list[float]], *, descending: bool = False) -> tuple[list[float], list[list[float]]]:
    """Return eigenvalues and eigenvectors for a symmetric matrix."""
    import numpy as np

    array = _numpy_real_array(matrix)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("matrix must be square")
    if not np.allclose(array, array.T, atol=1e-9, rtol=1e-9):
        raise ValueError("matrix must be symmetric")

    eigenvalues, eigenvectors = np.linalg.eigh(array)
    order = np.argsort(eigenvalues)
    if descending:
        order = order[::-1]
    sorted_values = eigenvalues[order]
    sorted_vectors = eigenvectors[:, order]
    return sorted_values.tolist(), [sorted_vectors[:, index].tolist() for index in range(sorted_vectors.shape[1])]


def mode_overlap_matrix(left_vectors: list[list[float]], right_vectors: list[list[float]]) -> list[list[float]]:
    """Return absolute overlap magnitudes |<v_i, w_j>| between two orthonormal bases."""
    import numpy as np

    left = _numpy_real_array(left_vectors)
    right = _numpy_real_array(right_vectors)
    if left.ndim != 2 or right.ndim != 2:
        raise ValueError("left_vectors and right_vectors must be 2D")
    if left.shape[1] != right.shape[1]:
        raise ValueError("vector dimensions must match")
    overlaps = np.abs(left @ right.T)
    return overlaps.tolist()


def symmetric_modal_control_weights(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    horizon: float,
    descending: bool = True,
) -> tuple[list[float], list[float]]:
    """Return exact diagonal Gramian weights in the symmetric A-eigenbasis."""
    import numpy as np

    if horizon < 0.0:
        raise ValueError("horizon must be non-negative")
    eigenvalues, eigenvectors = symmetric_eigensystem(a_matrix, descending=descending)
    q = _numpy_real_array(eigenvectors)
    b = _numpy_real_array(b_matrix)
    if b.ndim != 2 or b.shape[0] != q.shape[1]:
        raise ValueError("b_matrix row count must match the state dimension")

    modal_b = q @ b
    weights: list[float] = []
    for index, eigenvalue in enumerate(eigenvalues):
        projection_sq = float(np.sum(modal_b[index, :] ** 2))
        if abs(eigenvalue) <= 1e-12:
            factor = horizon
        else:
            factor = (math.exp(2.0 * eigenvalue * horizon) - 1.0) / (2.0 * eigenvalue)
        weights.append(projection_sq * factor)
    return eigenvalues, weights


def _numpy_real_array(values: list[list[float]] | list[float]) -> "numpy.ndarray":
    import numpy as np

    return np.asarray(values, dtype=float)


def _spectral_exponential_action(
    a_matrix: list[list[float]],
    rhs: list[list[float]] | list[float],
    time: float,
) -> "numpy.ndarray":
    import numpy as np

    a = _numpy_real_array(a_matrix)
    b = _numpy_real_array(rhs)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a_matrix must be square")

    if np.allclose(a, a.T, atol=1e-9, rtol=1e-9):
        eigenvalues, eigenvectors = np.linalg.eigh(a)
        transformed = eigenvectors.T @ b
        weighted = np.exp(eigenvalues * time)
        if transformed.ndim == 1:
            image = eigenvectors @ (weighted * transformed)
        else:
            image = eigenvectors @ (weighted[:, None] * transformed)
        return image

    eigenvalues, eigenvectors = np.linalg.eig(a)
    inverse = np.linalg.inv(eigenvectors)
    transformed = inverse @ b
    weighted = np.exp(eigenvalues * time)
    if transformed.ndim == 1:
        image = eigenvectors @ (weighted * transformed)
    else:
        image = eigenvectors @ (weighted[:, None] * transformed)
    if np.max(np.abs(np.imag(image))) > 1e-8:
        raise ValueError("matrix exponential produced a materially complex result")
    return np.real(image)


def state_transition_matrix(a_matrix: list[list[float]], time: float) -> list[list[float]]:
    """Return exp(A t) using spectral decomposition."""
    import numpy as np

    size = len(a_matrix)
    if size == 0:
        return []
    identity = np.eye(size, dtype=float)
    return _spectral_exponential_action(a_matrix, identity, time).tolist()


def finite_horizon_controllability_gramian(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    horizon: float,
    quadrature_steps: int = 256,
) -> list[list[float]]:
    """Return the finite-horizon controllability Gramian W(T)."""
    import numpy as np

    if horizon < 0.0:
        raise ValueError("horizon must be non-negative")
    if quadrature_steps <= 0:
        raise ValueError("quadrature_steps must be positive")

    a = _numpy_real_array(a_matrix)
    b = _numpy_real_array(b_matrix)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a_matrix must be square")
    if b.ndim != 2 or b.shape[0] != a.shape[0]:
        raise ValueError("b_matrix row count must match the state dimension")

    size = a.shape[0]
    if horizon == 0.0:
        return np.zeros((size, size), dtype=float).tolist()

    dt = horizon / quadrature_steps
    gramian = np.zeros((size, size), dtype=float)
    for index in range(quadrature_steps + 1):
        tau = index * dt
        propagated = _spectral_exponential_action(a_matrix, b_matrix, tau)
        weight = 0.5 if index in (0, quadrature_steps) else 1.0
        gramian += weight * (propagated @ propagated.T)
    gramian *= dt
    gramian = 0.5 * (gramian + gramian.T)
    return gramian.tolist()


def finite_horizon_controllability_gramian_symmetric(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    horizon: float,
) -> list[list[float]]:
    """Return the exact finite-horizon Gramian when A is symmetric."""
    import numpy as np

    if horizon < 0.0:
        raise ValueError("horizon must be non-negative")

    a = _numpy_real_array(a_matrix)
    b = _numpy_real_array(b_matrix)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("a_matrix must be square")
    if not np.allclose(a, a.T, atol=1e-9, rtol=1e-9):
        raise ValueError("a_matrix must be symmetric")
    if b.ndim != 2 or b.shape[0] != a.shape[0]:
        raise ValueError("b_matrix row count must match the state dimension")

    eigenvalues, eigenvectors = np.linalg.eigh(a)
    modal_b = eigenvectors.T @ b
    pair_sums = eigenvalues[:, None] + eigenvalues[None, :]
    kernel = np.empty_like(pair_sums)
    near_zero = np.abs(pair_sums) <= 1e-12
    kernel[near_zero] = horizon
    kernel[~near_zero] = (np.exp(pair_sums[~near_zero] * horizon) - 1.0) / pair_sums[~near_zero]
    modal_gramian = (modal_b @ modal_b.T) * kernel
    gramian = eigenvectors @ modal_gramian @ eigenvectors.T
    gramian = 0.5 * (gramian + gramian.T)
    return gramian.tolist()


def linear_minimum_energy_symmetric(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    delta_v0: list[float],
    delta_v_target: list[float],
    horizon: float,
) -> float:
    """Return the exact minimum energy when A is symmetric."""
    import numpy as np

    gramian = _numpy_real_array(
        finite_horizon_controllability_gramian_symmetric(
            a_matrix,
            b_matrix,
            horizon=horizon,
        )
    )
    mismatch = _numpy_real_array(
        endpoint_mismatch(
            a_matrix,
            delta_v0=delta_v0,
            delta_v_target=delta_v_target,
            horizon=horizon,
        )
    )
    if float(np.linalg.norm(mismatch)) <= 1e-15:
        return 0.0
    solution = np.linalg.solve(gramian, mismatch)
    return float(mismatch @ solution)


def linear_minimum_energy_symmetric_reachable(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    delta_v0: list[float],
    delta_v_target: list[float],
    horizon: float,
    tolerance: float = 1e-12,
) -> float:
    """Return the exact minimum energy if reachable, else ``math.inf``."""
    import numpy as np

    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")

    gramian = _numpy_real_array(
        finite_horizon_controllability_gramian_symmetric(
            a_matrix,
            b_matrix,
            horizon=horizon,
        )
    )
    mismatch = _numpy_real_array(
        endpoint_mismatch(
            a_matrix,
            delta_v0=delta_v0,
            delta_v_target=delta_v_target,
            horizon=horizon,
        )
    )
    if float(np.linalg.norm(mismatch)) <= 1e-15:
        return 0.0

    eigenvalues, eigenvectors = np.linalg.eigh(gramian)
    coordinates = eigenvectors.T @ mismatch
    positive = eigenvalues > tolerance
    unreachable = np.abs(coordinates[~positive])
    if unreachable.size and float(np.max(unreachable)) > tolerance:
        return math.inf
    if not np.any(positive):
        return math.inf
    return float(np.sum((coordinates[positive] ** 2) / eigenvalues[positive]))


def linear_minimum_energy_symmetric_projected(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    delta_v0: list[float],
    delta_v_target: list[float],
    horizon: float,
    tolerance: float = 1e-12,
) -> ProjectedSteeringReport:
    """Return reachable-projection energy and unreachable residual norm."""
    import numpy as np

    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")

    gramian = _numpy_real_array(
        finite_horizon_controllability_gramian_symmetric(
            a_matrix,
            b_matrix,
            horizon=horizon,
        )
    )
    mismatch = _numpy_real_array(
        endpoint_mismatch(
            a_matrix,
            delta_v0=delta_v0,
            delta_v_target=delta_v_target,
            horizon=horizon,
        )
    )
    eigenvalues, eigenvectors = np.linalg.eigh(gramian)
    coordinates = eigenvectors.T @ mismatch
    positive = eigenvalues > tolerance
    reachable_rank = int(np.sum(positive))
    unreachable_norm = float(np.linalg.norm(coordinates[~positive]))
    if reachable_rank == 0:
        energy = 0.0 if float(np.linalg.norm(mismatch)) <= 1e-15 else math.inf
    else:
        energy = float(np.sum((coordinates[positive] ** 2) / eigenvalues[positive]))
    return ProjectedSteeringReport(
        reachable_energy=energy,
        unreachable_residual_norm=unreachable_norm,
        reachable_rank=reachable_rank,
    )


def endpoint_mismatch(
    a_matrix: list[list[float]],
    *,
    delta_v0: list[float],
    delta_v_target: list[float],
    horizon: float,
) -> list[float]:
    """Return delta_v_target - exp(A T) delta_v0."""
    import numpy as np

    propagated = _spectral_exponential_action(a_matrix, delta_v0, horizon)
    mismatch = _numpy_real_array(delta_v_target) - propagated
    return np.asarray(mismatch, dtype=float).tolist()


def linear_minimum_energy(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    delta_v0: list[float],
    delta_v_target: list[float],
    horizon: float,
    quadrature_steps: int = 256,
) -> float:
    """Return the exact linear minimum energy for the chosen endpoint."""
    import numpy as np

    gramian = _numpy_real_array(
        finite_horizon_controllability_gramian(
            a_matrix,
            b_matrix,
            horizon=horizon,
            quadrature_steps=quadrature_steps,
        )
    )
    mismatch = _numpy_real_array(
        endpoint_mismatch(
            a_matrix,
            delta_v0=delta_v0,
            delta_v_target=delta_v_target,
            horizon=horizon,
        )
    )
    solution = np.linalg.solve(gramian, mismatch)
    return float(mismatch @ solution)


def linear_minimum_energy_control_samples(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    delta_v0: list[float],
    delta_v_target: list[float],
    horizon: float,
    sample_times: list[float],
    quadrature_steps: int = 256,
) -> list[list[float]]:
    """Return samples of the minimum-energy control input u*(t)."""
    import numpy as np

    gramian = _numpy_real_array(
        finite_horizon_controllability_gramian(
            a_matrix,
            b_matrix,
            horizon=horizon,
            quadrature_steps=quadrature_steps,
        )
    )
    mismatch = _numpy_real_array(
        endpoint_mismatch(
            a_matrix,
            delta_v0=delta_v0,
            delta_v_target=delta_v_target,
            horizon=horizon,
        )
    )
    multiplier = np.linalg.solve(gramian, mismatch)
    b = _numpy_real_array(b_matrix)

    controls: list[list[float]] = []
    for time in sample_times:
        if time < 0.0 or time > horizon:
            raise ValueError("sample_times must lie inside [0, horizon]")
        propagated = _spectral_exponential_action(
            _numpy_real_array(a_matrix).T.tolist(),
            multiplier.tolist(),
            horizon - time,
        )
        control = b.T @ propagated
        controls.append(np.asarray(control, dtype=float).tolist())
    return controls


def finite_horizon_controllability_report(
    a_matrix: list[list[float]],
    b_matrix: list[list[float]],
    *,
    horizon: float,
    quadrature_steps: int = 256,
    tolerance: float = 1e-10,
) -> ControllabilityReport:
    """Return a numerical controllability summary from the Gramian spectrum."""
    import numpy as np

    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")

    gramian = np.asarray(
        finite_horizon_controllability_gramian(
            a_matrix,
            b_matrix,
            horizon=horizon,
            quadrature_steps=quadrature_steps,
        ),
        dtype=float,
    )
    eigenvalues = np.linalg.eigvalsh(gramian)
    numerical_rank = int(np.sum(eigenvalues > tolerance))
    positive = eigenvalues[eigenvalues > tolerance]
    state_dimension = gramian.shape[0]
    input_dimension = len(b_matrix[0]) if b_matrix else 0
    return ControllabilityReport(
        state_dimension=state_dimension,
        input_dimension=input_dimension,
        horizon=horizon,
        numerical_rank=numerical_rank,
        tolerance=tolerance,
        min_eigenvalue=float(eigenvalues[0]),
        max_eigenvalue=float(eigenvalues[-1]),
        smallest_positive_eigenvalue=float(positive[0]) if positive.size else 0.0,
        reachable_condition_number=(
            float(positive[-1] / positive[0])
            if positive.size
            else math.inf
        ),
        controllable_at_tolerance=(numerical_rank == state_dimension),
    )


def control_energy(control_trace: list[list[float]], dt: float, energy_weight: float = 1.0) -> float:
    """Quadratic intervention energy integral approximated by a time sum."""
    total = 0.0
    for control in control_trace:
        total += sum(value * value for value in control)
    return energy_weight * dt * total


def project_control_to_sites(operator: list[list[float]], control: list[float]) -> list[float]:
    """Map control-channel amplitudes into site-level additive forcing."""
    if not operator:
        return []
    if any(len(row) != len(control) for row in operator):
        raise ValueError("control dimension must match actuator operator column count")
    return [
        sum(weight * value for weight, value in zip(row, control, strict=True))
        for row in operator
    ]


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
