"""Auditable finite-dimensional information-map toy model.

This module is a *semiclassical code-subspace demonstration only*.  It
constructs a finite isometry and checks its quantum-information bookkeeping.
It does not factorize the physical Hilbert space of gravity, supply a
covariant black-hole-to-cosmology transition, or derive horizon thermodynamics
or H-SAT.  In particular, the ``R`` and ``C`` tensor factors below are an
explicit approximation whose physical dressing, constraints, and error budget
would have to be supplied by a gravitational completion.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atan2, cos, isfinite, log, sin, sqrt
from numbers import Real
from typing import Sequence

Matrix = list[list[complex]]
_DEFAULT_TOLERANCE = 1.0e-10


def _finite_real(name: str, value: object) -> float:
    """Return a finite real value without silently treating booleans as data."""

    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a real number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _require_finite_matrix(matrix: Sequence[Sequence[complex]]) -> None:
    """Reject matrix entries whose real or imaginary part is non-finite."""

    shape(matrix)
    for row in matrix:
        for entry in row:
            try:
                value = complex(entry)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(
                    "matrix entries must be finite complex numbers"
                ) from exc
            if not isfinite(value.real) or not isfinite(value.imag):
                raise ValueError("matrix entries must be finite")


@dataclass(frozen=True, slots=True)
class TransitionCodeSpec:
    """Finite code and coherent location-erasure probability.

    ``d`` is the dimension of the input code.  Each output subsystem has
    dimension ``d + 1``; its final basis state is the orthogonal blank state.
    ``p`` is a toy-map parameter, not a derived branching probability.
    """

    d: int
    p: float

    def __post_init__(self) -> None:
        if isinstance(self.d, bool) or not isinstance(self.d, int) or self.d < 1:
            raise ValueError("d must be an integer greater than or equal to 1")
        probability = _finite_real("p", self.p)
        if not 0.0 <= probability <= 1.0:
            raise ValueError("p must be finite and lie in [0, 1]")
        object.__setattr__(self, "p", probability)

    @property
    def local_dimension(self) -> int:
        """Dimension of either output factor, including its blank state."""

        return self.d + 1

    @property
    def output_dimension(self) -> int:
        """Dimension of the effective R tensor C output space."""

        return self.local_dimension**2

    @property
    def blank_index(self) -> int:
        """The basis index of the blank state in either output factor."""

        return self.d


def zeros(rows: int, columns: int) -> Matrix:
    """Return a complex zero matrix, rejecting non-positive dimensions."""

    if rows < 1 or columns < 1:
        raise ValueError("matrix dimensions must be positive")
    return [[0.0j for _ in range(columns)] for _ in range(rows)]


def identity(dimension: int) -> Matrix:
    """Return an identity matrix."""

    result = zeros(dimension, dimension)
    for index in range(dimension):
        result[index][index] = 1.0
    return result


def shape(matrix: Sequence[Sequence[complex]]) -> tuple[int, int]:
    """Return a non-empty rectangular matrix shape."""

    if not matrix or not matrix[0]:
        raise ValueError("matrix must be non-empty")
    rows, columns = len(matrix), len(matrix[0])
    if any(len(row) != columns for row in matrix):
        raise ValueError("matrix must be rectangular")
    return rows, columns


def dagger(matrix: Sequence[Sequence[complex]]) -> Matrix:
    """Return the conjugate transpose."""

    rows, columns = shape(matrix)
    return [
        [complex(matrix[row][column]).conjugate() for row in range(rows)]
        for column in range(columns)
    ]


def matmul(
    left: Sequence[Sequence[complex]], right: Sequence[Sequence[complex]]
) -> Matrix:
    """Return a dense matrix product using only the Python standard library."""

    left_rows, left_columns = shape(left)
    right_rows, right_columns = shape(right)
    if left_columns != right_rows:
        raise ValueError("matrix dimensions do not compose")
    result = zeros(left_rows, right_columns)
    for row in range(left_rows):
        for pivot in range(left_columns):
            coefficient = complex(left[row][pivot])
            if coefficient == 0.0j:
                continue
            for column in range(right_columns):
                result[row][column] += coefficient * complex(right[pivot][column])
    return result


def trace(matrix: Sequence[Sequence[complex]]) -> complex:
    """Return a square matrix trace."""

    rows, columns = shape(matrix)
    if rows != columns:
        raise ValueError("trace requires a square matrix")
    return sum((complex(matrix[index][index]) for index in range(rows)), 0.0j)


def max_abs_difference(
    left: Sequence[Sequence[complex]], right: Sequence[Sequence[complex]]
) -> float:
    """Return the entrywise max norm of a matrix difference."""

    if shape(left) != shape(right):
        raise ValueError("matrix shapes differ")
    return max(
        abs(complex(left[row][column]) - complex(right[row][column]))
        for row in range(len(left))
        for column in range(len(left[0]))
    )


def is_hermitian(
    matrix: Sequence[Sequence[complex]], tolerance: float = _DEFAULT_TOLERANCE
) -> bool:
    """Check Hermiticity to a supplied absolute tolerance."""

    return max_abs_difference(matrix, dagger(matrix)) <= tolerance


def state_density(
    amplitudes: Sequence[complex], tolerance: float = _DEFAULT_TOLERANCE
) -> Matrix:
    """Construct |psi><psi| after checking that |psi> is normalized."""

    if not amplitudes:
        raise ValueError("a state must contain at least one amplitude")
    try:
        vector = [complex(value) for value in amplitudes]
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("state amplitudes must be finite complex numbers") from exc
    if any(not isfinite(value.real) or not isfinite(value.imag) for value in vector):
        raise ValueError("state amplitudes must be finite complex numbers")
    norm = sum((abs(value) ** 2 for value in vector), 0.0)
    if not isfinite(norm) or abs(norm - 1.0) > tolerance:
        raise ValueError("state amplitudes must be normalized")
    return [
        [
            vector[row] * vector[column].conjugate()
            for column in range(len(vector))
        ]
        for row in range(len(vector))
    ]


def coherent_erasure_isometry(spec: TransitionCodeSpec) -> Matrix:
    """Build V_p|i> = sqrt(1-p)|i,blank> + sqrt(p)|blank,i>.

    Rows use the product index ``r * local_dimension + c``.  The blank state
    makes the two terms orthogonal, so this matrix is an isometry for all
    permitted ``p``.
    """

    matrix = zeros(spec.output_dimension, spec.d)
    to_r = sqrt(1.0 - spec.p)
    to_c = sqrt(spec.p)
    local, blank = spec.local_dimension, spec.blank_index
    for code_index in range(spec.d):
        matrix[code_index * local + blank][code_index] = to_r
        matrix[blank * local + code_index][code_index] = to_c
    return matrix


def apply_isometry(spec: TransitionCodeSpec, rho: Sequence[Sequence[complex]]) -> Matrix:
    """Apply the channel Phi(rho)=V rho V-dagger to a d-dimensional density matrix."""

    if shape(rho) != (spec.d, spec.d):
        raise ValueError("rho must have shape (d, d)")
    _require_finite_matrix(rho)
    # This validates Hermiticity, unit trace, and positive semidefiniteness in
    # addition to finiteness.  Its return value is not needed here.
    von_neumann_entropy(rho)
    isometry = coherent_erasure_isometry(spec)
    return matmul(matmul(isometry, rho), dagger(isometry))


def recover_code_state(
    spec: TransitionCodeSpec, rho_rc: Sequence[Sequence[complex]]
) -> Matrix:
    """Compress a joint state with V-dagger and V.

    This exactly recovers inputs of ``apply_isometry`` because they lie in the
    range of ``V``.  It is not a trace-preserving recovery channel for an
    arbitrary joint state outside that code range, and it is not recovery from
    either marginal alone.
    """

    if shape(rho_rc) != (spec.output_dimension, spec.output_dimension):
        raise ValueError("rho_rc has the wrong output shape")
    _require_finite_matrix(rho_rc)
    isometry = coherent_erasure_isometry(spec)
    return matmul(matmul(dagger(isometry), rho_rc), isometry)


def partial_trace_child(
    spec: TransitionCodeSpec, rho_rc: Sequence[Sequence[complex]]
) -> Matrix:
    """Trace C from an R tensor C matrix, returning the R marginal."""

    if shape(rho_rc) != (spec.output_dimension, spec.output_dimension):
        raise ValueError("rho_rc has the wrong output shape")
    local = spec.local_dimension
    result = zeros(local, local)
    for r_row in range(local):
        for r_column in range(local):
            result[r_row][r_column] = sum(
                (
                    complex(
                        rho_rc[r_row * local + child][r_column * local + child]
                    )
                    for child in range(local)
                ),
                0.0j,
            )
    return result


def partial_trace_radiation(
    spec: TransitionCodeSpec, rho_rc: Sequence[Sequence[complex]]
) -> Matrix:
    """Trace R from an R tensor C matrix, returning the C marginal."""

    if shape(rho_rc) != (spec.output_dimension, spec.output_dimension):
        raise ValueError("rho_rc has the wrong output shape")
    local = spec.local_dimension
    result = zeros(local, local)
    for c_row in range(local):
        for c_column in range(local):
            result[c_row][c_column] = sum(
                (
                    complex(
                        rho_rc[radiation * local + c_row][
                            radiation * local + c_column
                        ]
                    )
                    for radiation in range(local)
                ),
                0.0j,
            )
    return result


def expected_radiation_marginal(
    spec: TransitionCodeSpec, rho: Sequence[Sequence[complex]]
) -> Matrix:
    """Return (1-p)rho plus p times the R blank-state projector."""

    if shape(rho) != (spec.d, spec.d):
        raise ValueError("rho must have shape (d, d)")
    result = zeros(spec.local_dimension, spec.local_dimension)
    for row in range(spec.d):
        for column in range(spec.d):
            result[row][column] = (1.0 - spec.p) * complex(rho[row][column])
    result[spec.blank_index][spec.blank_index] = spec.p
    return result


def expected_child_marginal(
    spec: TransitionCodeSpec, rho: Sequence[Sequence[complex]]
) -> Matrix:
    """Return p rho plus (1-p) times the C blank-state projector."""

    if shape(rho) != (spec.d, spec.d):
        raise ValueError("rho must have shape (d, d)")
    result = zeros(spec.local_dimension, spec.local_dimension)
    for row in range(spec.d):
        for column in range(spec.d):
            result[row][column] = spec.p * complex(rho[row][column])
    result[spec.blank_index][spec.blank_index] = 1.0 - spec.p
    return result


def binary_entropy(probability: float) -> float:
    """Return h_2(p) in nats, with h_2(0)=h_2(1)=0."""

    probability = _finite_real("probability", probability)
    if not 0.0 <= probability <= 1.0:
        raise ValueError("probability must be finite and lie in [0, 1]")
    if probability in (0.0, 1.0):
        return 0.0
    return -probability * log(probability) - (1.0 - probability) * log(1.0 - probability)


def _real_symmetric_eigenvalues(
    matrix: list[list[float]], tolerance: float = 1.0e-13
) -> list[float]:
    """Diagonalize a real symmetric matrix by Jacobi rotations.

    This small dependency-free routine is sufficient for the low-dimensional
    matrices in this demonstrator.  It is not intended as a large-scale
    numerical-linear-algebra replacement.
    """

    dimension = len(matrix)
    working = [row[:] for row in matrix]
    for _ in range(max(16, 100 * dimension * dimension)):
        largest, row, column = 0.0, 0, 0
        for first in range(dimension):
            for second in range(first + 1, dimension):
                candidate = abs(working[first][second])
                if candidate > largest:
                    largest, row, column = candidate, first, second
        if largest <= tolerance:
            return [working[index][index] for index in range(dimension)]
        angle = 0.5 * atan2(
            2.0 * working[row][column],
            working[column][column] - working[row][row],
        )
        cosine, sine = cos(angle), sin(angle)
        diagonal_row, diagonal_column = working[row][row], working[column][column]
        off_diagonal = working[row][column]
        working[row][row] = (
            cosine * cosine * diagonal_row
            - 2.0 * sine * cosine * off_diagonal
            + sine * sine * diagonal_column
        )
        working[column][column] = (
            sine * sine * diagonal_row
            + 2.0 * sine * cosine * off_diagonal
            + cosine * cosine * diagonal_column
        )
        working[row][column] = working[column][row] = 0.0
        for index in range(dimension):
            if index in (row, column):
                continue
            row_value, column_value = working[index][row], working[index][column]
            working[index][row] = working[row][index] = cosine * row_value - sine * column_value
            working[index][column] = working[column][index] = (
                sine * row_value + cosine * column_value
            )
    raise ValueError("Jacobi eigensolver did not converge")


def von_neumann_entropy(
    rho: Sequence[Sequence[complex]], tolerance: float = _DEFAULT_TOLERANCE
) -> float:
    """Return S_vN(rho) in nats for a low-dimensional Hermitian density matrix."""

    rows, columns = shape(rho)
    _require_finite_matrix(rho)
    if rows != columns or not is_hermitian(rho, tolerance):
        raise ValueError("rho must be a Hermitian square matrix")
    if abs(trace(rho).imag) > tolerance or abs(trace(rho).real - 1.0) > tolerance:
        raise ValueError("rho must have unit trace")
    # A complex Hermitian n-by-n matrix becomes a real symmetric 2n-by-2n
    # matrix with every eigenvalue duplicated.  Half its entropy is S(rho).
    real = [
        [float(complex(rho[row][column]).real) for column in range(rows)]
        for row in range(rows)
    ]
    imaginary = [
        [float(complex(rho[row][column]).imag) for column in range(rows)]
        for row in range(rows)
    ]
    doubled = [
        real[row] + [-imaginary[row][column] for column in range(rows)]
        for row in range(rows)
    ] + [
        imaginary[row] + real[row]
        for row in range(rows)
    ]
    values = _real_symmetric_eigenvalues(doubled)
    entropy = 0.0
    for value in values:
        if value < -tolerance:
            raise ValueError("rho is not positive semidefinite")
        if value > tolerance:
            entropy -= 0.5 * value * log(value)
    return entropy


def analytic_radiation_entropy(
    spec: TransitionCodeSpec, input_entropy_nats: float
) -> float:
    """Return h_2(p)+(1-p)S(rho) for the R erasure marginal."""

    input_entropy_nats = _finite_real("input_entropy_nats", input_entropy_nats)
    if input_entropy_nats < -_DEFAULT_TOLERANCE:
        raise ValueError("entropy cannot be negative")
    return binary_entropy(spec.p) + (1.0 - spec.p) * max(0.0, input_entropy_nats)


def analytic_child_entropy(spec: TransitionCodeSpec, input_entropy_nats: float) -> float:
    """Return h_2(p)+p S(rho) for the C erasure marginal."""

    input_entropy_nats = _finite_real("input_entropy_nats", input_entropy_nats)
    if input_entropy_nats < -_DEFAULT_TOLERANCE:
        raise ValueError("entropy cannot be negative")
    return binary_entropy(spec.p) + spec.p * max(0.0, input_entropy_nats)


def no_cloning_overlap_contrast() -> dict[str, float]:
    """Return the |0>, |+> overlap contradiction for a putative perfect cloner.

    A linear isometry preserves the input overlap 1/sqrt(2); perfect cloning
    would square it to 1/2.  This is a general quantum-information contrast,
    not a resolution of a gravitational no-cloning problem.
    """

    return {
        "input_overlap": 1.0 / sqrt(2.0),
        "putative_cloned_output_overlap": 0.5,
        "absolute_difference": 1.0 / sqrt(2.0) - 0.5,
    }


def validate_transition(
    spec: TransitionCodeSpec,
    rho: Sequence[Sequence[complex]],
    tolerance: float = _DEFAULT_TOLERANCE,
) -> dict[str, float | bool]:
    """Audit the toy map's mathematical invariants for a supplied input state."""

    isometry = coherent_erasure_isometry(spec)
    output = apply_isometry(spec, rho)
    recovered = recover_code_state(spec, output)
    isometry_error = max_abs_difference(
        matmul(dagger(isometry), isometry), identity(spec.d)
    )
    trace_error = abs(trace(output) - trace(rho))
    hermiticity_error = max_abs_difference(output, dagger(output))
    recovery_error = max_abs_difference(recovered, rho)
    # A finite entropy value also checks positive semidefiniteness of the
    # generated output with the same tolerance used elsewhere in this module.
    von_neumann_entropy(output, tolerance)
    return {
        "isometry_error": isometry_error,
        "trace_error": trace_error,
        "hermiticity_error": hermiticity_error,
        "recovery_error": recovery_error,
        "is_isometry": isometry_error <= tolerance,
        "is_trace_preserving": trace_error <= tolerance,
        "output_is_hermitian": hermiticity_error <= tolerance,
        "output_is_positive_semidefinite": True,
        "is_recoverable_on_code_range": recovery_error <= tolerance,
    }
