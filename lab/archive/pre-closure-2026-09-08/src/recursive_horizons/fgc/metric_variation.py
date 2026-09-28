"""Exact algebraic controls for the FGC-1 covariant metric variation.

The action-level convention is the one frozen by :mod:`.action`.  In four
dimensions, with ``R^rho_{ sigma mu nu}`` using the ``+partial_mu Gamma_nu``
sign and with variation performed with respect to ``g^{ab}``, the
Gauss--Bonnet stress of the action term ``+f(phi) G`` is

``T_GB_ab = -8 P_acbd nabla^c nabla^d f``.

This module checks the compact expression against a separately implemented
seven-term expansion on exact algebraic-curvature fixtures.  It also checks
the nonminimal ``F(phi) R / 2`` derivative source, the trace and Noether
curvature-contraction identities, and an independent spatially flat FLRW lapse
variation evaluated by exact automatic differentiation of the reduced action.
It does not perform a spherical reduction, principal-symbol calculation, or
collapse evolution.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from numbers import Integral
from typing import Any, Mapping, Sequence

from .action import ActionParameters, ModelID


DIMENSION = 4
MINKOWSKI_DIAGONAL = (
    Fraction(-1),
    Fraction(1),
    Fraction(1),
    Fraction(1),
)
LCG64_MULTIPLIER = 6364136223846793005
LCG64_INCREMENT = 1442695040888963407
LCG64_MASK = (1 << 64) - 1

Vector = tuple[Fraction, ...]
Matrix = tuple[tuple[Fraction, ...], ...]
Tensor4 = tuple[
    tuple[tuple[tuple[Fraction, ...], ...], ...], ...
]

EXPANDED_GB_COEFFICIENTS: dict[str, Fraction] = {
    "riemann_hessian": Fraction(-8),
    "metric_ricci_hessian": Fraction(8),
    "ricci_b_hessian_a": Fraction(-8),
    "ricci_a_hessian_b": Fraction(-8),
    "ricci_ab_box": Fraction(8),
    "scalar_metric_box": Fraction(-4),
    "scalar_hessian": Fraction(4),
}


def _fraction(name: str, value: Fraction | Integral) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an exact rational number")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value))
    raise ValueError(f"{name} must be an exact rational number")


def _parameter_fraction(name: str, value: float) -> Fraction:
    """Convert a validated finite ACT1 decimal fixture to an exact rational."""

    try:
        return Fraction(str(value))
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} cannot be represented exactly") from exc


def _validate_metric_diagonal(metric_diagonal: Sequence[Fraction | Integral]) -> Vector:
    if len(metric_diagonal) != DIMENSION:
        raise ValueError("metric_diagonal must contain exactly four entries")
    metric = tuple(
        _fraction(f"metric_diagonal[{index}]", value)
        for index, value in enumerate(metric_diagonal)
    )
    if any(value == 0 for value in metric):
        raise ValueError("metric_diagonal must be nonsingular")
    return metric


def _validate_vector(name: str, vector: Sequence[Fraction | Integral]) -> Vector:
    if len(vector) != DIMENSION:
        raise ValueError(f"{name} must contain exactly four entries")
    return tuple(_fraction(f"{name}[{index}]", value) for index, value in enumerate(vector))


def _validate_matrix(name: str, matrix: Sequence[Sequence[Fraction | Integral]]) -> Matrix:
    if len(matrix) != DIMENSION or any(len(row) != DIMENSION for row in matrix):
        raise ValueError(f"{name} must be a four-by-four matrix")
    return tuple(
        tuple(_fraction(f"{name}[{i}][{j}]", matrix[i][j]) for j in range(DIMENSION))
        for i in range(DIMENSION)
    )


def _validate_tensor4(
    name: str,
    tensor: Sequence[Sequence[Sequence[Sequence[Fraction | Integral]]]],
) -> Tensor4:
    if (
        len(tensor) != DIMENSION
        or any(len(tensor[a]) != DIMENSION for a in range(DIMENSION))
        or any(
            len(tensor[a][b]) != DIMENSION
            for a in range(DIMENSION)
            for b in range(DIMENSION)
        )
        or any(
            len(tensor[a][b][c]) != DIMENSION
            for a in range(DIMENSION)
            for b in range(DIMENSION)
            for c in range(DIMENSION)
        )
    ):
        raise ValueError(f"{name} must be a four-by-four-by-four-by-four tensor")
    return tuple(
        tuple(
            tuple(
                tuple(
                    _fraction(f"{name}[{a}][{b}][{c}][{d}]", tensor[a][b][c][d])
                    for d in range(DIMENSION)
                )
                for c in range(DIMENSION)
            )
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )


def _metric_component(metric_diagonal: Vector, a: int, b: int) -> Fraction:
    return metric_diagonal[a] if a == b else Fraction(0)


def _inverse_metric_diagonal(metric_diagonal: Vector) -> Vector:
    return tuple(Fraction(1, value) for value in metric_diagonal)


def _zero_matrix() -> list[list[Fraction]]:
    return [[Fraction(0) for _ in range(DIMENSION)] for _ in range(DIMENSION)]


def _zero_tensor4() -> list[list[list[list[Fraction]]]]:
    return [
        [
            [[Fraction(0) for _ in range(DIMENSION)] for _ in range(DIMENSION)]
            for _ in range(DIMENSION)
        ]
        for _ in range(DIMENSION)
    ]


def _freeze_matrix(matrix: Sequence[Sequence[Fraction]]) -> Matrix:
    return tuple(tuple(value for value in row) for row in matrix)


def _freeze_tensor4(
    tensor: Sequence[Sequence[Sequence[Sequence[Fraction]]]],
) -> Tensor4:
    return tuple(
        tuple(
            tuple(tuple(value for value in line) for line in plane)
            for plane in block
        )
        for block in tensor
    )


def _matrix_is_symmetric(matrix: Matrix) -> bool:
    return all(
        matrix[a][b] == matrix[b][a]
        for a in range(DIMENSION)
        for b in range(DIMENSION)
    )


def _matrix_is_zero(matrix: Matrix) -> bool:
    return all(value == 0 for row in matrix for value in row)


def _tensor4_is_zero(tensor: Tensor4) -> bool:
    return all(
        tensor[a][b][c][d] == 0
        for a in range(DIMENSION)
        for b in range(DIMENSION)
        for c in range(DIMENSION)
        for d in range(DIMENSION)
    )


def _max_abs_matrix(matrix: Matrix) -> Fraction:
    return max((abs(value) for row in matrix for value in row), default=Fraction(0))


def _subtract_matrices(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[a][b] - right[a][b] for b in range(DIMENSION))
        for a in range(DIMENSION)
    )


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _canonical_exact(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if isinstance(value, tuple):
        return [_canonical_exact(item) for item in value]
    if isinstance(value, list):
        return [_canonical_exact(item) for item in value]
    if isinstance(value, dict):
        return {key: _canonical_exact(item) for key, item in value.items()}
    return value


@dataclass(frozen=True, slots=True)
class AlgebraicCurvatureFixture:
    """One exact pointwise curvature/scalar-jet fixture."""

    fixture_id: str
    metric_diagonal: Vector
    riemann_lower: Tensor4
    scalar_field: Fraction
    gradient_scalar_lower: Vector
    hessian_scalar_lower: Matrix

    def __post_init__(self) -> None:
        if not isinstance(self.fixture_id, str) or not self.fixture_id:
            raise ValueError("fixture_id must be a non-empty string")
        object.__setattr__(
            self,
            "metric_diagonal",
            _validate_metric_diagonal(self.metric_diagonal),
        )
        object.__setattr__(
            self,
            "riemann_lower",
            _validate_tensor4("riemann_lower", self.riemann_lower),
        )
        object.__setattr__(
            self,
            "scalar_field",
            _fraction("scalar_field", self.scalar_field),
        )
        object.__setattr__(
            self,
            "gradient_scalar_lower",
            _validate_vector("gradient_scalar_lower", self.gradient_scalar_lower),
        )
        object.__setattr__(
            self,
            "hessian_scalar_lower",
            _validate_matrix("hessian_scalar_lower", self.hessian_scalar_lower),
        )

    def source_sha256(self) -> str:
        payload = {
            "fixture_id": self.fixture_id,
            "metric_diagonal": self.metric_diagonal,
            "riemann_lower": self.riemann_lower,
            "scalar_field": self.scalar_field,
            "gradient_scalar_lower": self.gradient_scalar_lower,
            "hessian_scalar_lower": self.hessian_scalar_lower,
        }
        encoded = json.dumps(
            _canonical_exact(payload), sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
        return sha256(encoded).hexdigest()


@dataclass(frozen=True, slots=True)
class _ExactDual:
    """One exact value/tangent pair for dependency-free differentiation."""

    value: Fraction
    tangent: Fraction

    def __mul__(self, other: _ExactDual | Fraction | Integral) -> _ExactDual:
        if not isinstance(other, _ExactDual):
            other = _ExactDual(_fraction("dual multiplier", other), Fraction(0))
        return _ExactDual(
            self.value * other.value,
            self.tangent * other.value + self.value * other.tangent,
        )

    def __rmul__(self, other: Fraction | Integral) -> _ExactDual:
        return self * other

    def __pow__(self, power: int) -> _ExactDual:
        if isinstance(power, bool) or not isinstance(power, int):
            raise ValueError("dual power must be an integer")
        if self.value == 0 and power < 0:
            raise ValueError("cannot raise a zero dual value to a negative power")
        if power == 0:
            return _ExactDual(Fraction(1), Fraction(0))
        return _ExactDual(
            self.value**power,
            power * self.value ** (power - 1) * self.tangent,
        )


class _LCG64:
    def __init__(self, seed: int) -> None:
        if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        self.state = seed & LCG64_MASK

    def next(self) -> int:
        self.state = (
            LCG64_MULTIPLIER * self.state + LCG64_INCREMENT
        ) & LCG64_MASK
        return self.state

    def bounded_integer(self, bound: int) -> int:
        if isinstance(bound, bool) or not isinstance(bound, int) or bound <= 0:
            raise ValueError("integer_bound must be a positive integer")
        return int((self.next() >> 16) % (2 * bound + 1)) - bound


def _generated_symmetric_matrix(generator: _LCG64, bound: int) -> Matrix:
    matrix = _zero_matrix()
    for a in range(DIMENSION):
        for b in range(a, DIMENSION):
            value = Fraction(generator.bounded_integer(bound))
            matrix[a][b] = value
            matrix[b][a] = value
    return _freeze_matrix(matrix)


def kulkarni_nomizu_product(left: Matrix, right: Matrix) -> Tensor4:
    """Return the Kulkarni--Nomizu product of two symmetric tensors."""

    left = _validate_matrix("left", left)
    right = _validate_matrix("right", right)
    if not _matrix_is_symmetric(left) or not _matrix_is_symmetric(right):
        raise ValueError("Kulkarni--Nomizu inputs must be symmetric")
    tensor = _zero_tensor4()
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            for c in range(DIMENSION):
                for d in range(DIMENSION):
                    tensor[a][b][c][d] = (
                        left[a][c] * right[b][d]
                        + right[a][c] * left[b][d]
                        - left[a][d] * right[b][c]
                        - right[a][d] * left[b][c]
                    )
    return _freeze_tensor4(tensor)


def _add_tensor4(left: Tensor4, right: Tensor4) -> Tensor4:
    return tuple(
        tuple(
            tuple(
                tuple(
                    left[a][b][c][d] + right[a][b][c][d]
                    for d in range(DIMENSION)
                )
                for c in range(DIMENSION)
            )
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )


def generate_algebraic_curvature_fixtures(
    seed: int,
    fixture_count: int,
    integer_bound: int,
) -> tuple[AlgebraicCurvatureFixture, ...]:
    """Generate exact generic Riemann tensors with a frozen LCG contract."""

    if isinstance(fixture_count, bool) or not isinstance(fixture_count, int) or fixture_count <= 0:
        raise ValueError("fixture_count must be a positive integer")
    if isinstance(integer_bound, bool) or not isinstance(integer_bound, int) or integer_bound <= 0:
        raise ValueError("integer_bound must be a positive integer")
    generator = _LCG64(seed)
    accepted: list[AlgebraicCurvatureFixture] = []
    attempts = 0
    while len(accepted) < fixture_count and attempts < fixture_count * 100:
        attempts += 1
        symmetric = [
            _generated_symmetric_matrix(generator, integer_bound) for _ in range(4)
        ]
        riemann = _add_tensor4(
            kulkarni_nomizu_product(symmetric[0], symmetric[1]),
            kulkarni_nomizu_product(symmetric[2], symmetric[3]),
        )
        scalar_field = Fraction(generator.bounded_integer(integer_bound))
        gradient = tuple(
            Fraction(generator.bounded_integer(integer_bound))
            for _ in range(DIMENSION)
        )
        hessian = _generated_symmetric_matrix(generator, integer_bound)
        fixture = AlgebraicCurvatureFixture(
            fixture_id=f"lcg64-{len(accepted):02d}-draw-{attempts:02d}",
            metric_diagonal=MINKOWSKI_DIAGONAL,
            riemann_lower=riemann,
            scalar_field=scalar_field,
            gradient_scalar_lower=gradient,
            hessian_scalar_lower=hessian,
        )
        validation = validate_algebraic_curvature_fixture(fixture)
        if all(
            validation[key]
            for key in (
                "riemann_first_pair_antisymmetric",
                "riemann_second_pair_antisymmetric",
                "riemann_pair_exchange_symmetric",
                "riemann_first_bianchi",
                "ricci_symmetric",
                "nonzero_curvature",
                "nonzero_gauss_bonnet_invariant",
                "nonzero_scalar_field",
                "nonzero_scalar_gradient",
                "nonzero_scalar_hessian",
                "generic_gb_source_nonzero",
            )
        ):
            accepted.append(fixture)
    if len(accepted) != fixture_count:
        raise ValueError("fixture generator could not satisfy exact coverage gates")
    return tuple(accepted)


def ricci_tensor(metric_diagonal: Vector, riemann_lower: Tensor4) -> Matrix:
    metric = _validate_metric_diagonal(metric_diagonal)
    riemann = _validate_tensor4("riemann_lower", riemann_lower)
    inverse = _inverse_metric_diagonal(metric)
    return tuple(
        tuple(
            sum(
                inverse[index] * riemann[index][b][index][d]
                for index in range(DIMENSION)
            )
            for d in range(DIMENSION)
        )
        for b in range(DIMENSION)
    )


def scalar_curvature(metric_diagonal: Vector, ricci_lower: Matrix) -> Fraction:
    metric = _validate_metric_diagonal(metric_diagonal)
    ricci = _validate_matrix("ricci_lower", ricci_lower)
    inverse = _inverse_metric_diagonal(metric)
    return sum(inverse[index] * ricci[index][index] for index in range(DIMENSION))


def gauss_bonnet_invariant(metric_diagonal: Vector, riemann_lower: Tensor4) -> Fraction:
    metric = _validate_metric_diagonal(metric_diagonal)
    riemann = _validate_tensor4("riemann_lower", riemann_lower)
    inverse = _inverse_metric_diagonal(metric)
    ricci = ricci_tensor(metric, riemann)
    scalar = scalar_curvature(metric, ricci)
    riemann_squared = sum(
        inverse[a]
        * inverse[b]
        * inverse[c]
        * inverse[d]
        * riemann[a][b][c][d]
        * riemann[a][b][c][d]
        for a in range(DIMENSION)
        for b in range(DIMENSION)
        for c in range(DIMENSION)
        for d in range(DIMENSION)
    )
    ricci_squared = sum(
        inverse[a] * inverse[b] * ricci[a][b] * ricci[a][b]
        for a in range(DIMENSION)
        for b in range(DIMENSION)
    )
    return riemann_squared - 4 * ricci_squared + scalar * scalar


def double_dual_p_tensor(metric_diagonal: Vector, riemann_lower: Tensor4) -> Tensor4:
    """Return the explicit divergence-free double-dual curvature tensor.

    The component authority is
    ``P_abcd = R_abcd - 2 g_a[c R_d]b + 2 g_b[c R_d]a
                  + R g_a[c g_d]b``
    with antisymmetrization weight one half.
    """

    metric = _validate_metric_diagonal(metric_diagonal)
    riemann = _validate_tensor4("riemann_lower", riemann_lower)
    ricci = ricci_tensor(metric, riemann)
    scalar = scalar_curvature(metric, ricci)
    tensor = _zero_tensor4()
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            for c in range(DIMENSION):
                for d in range(DIMENSION):
                    tensor[a][b][c][d] = (
                        riemann[a][b][c][d]
                        - _metric_component(metric, a, c) * ricci[d][b]
                        + _metric_component(metric, a, d) * ricci[c][b]
                        + _metric_component(metric, b, c) * ricci[d][a]
                        - _metric_component(metric, b, d) * ricci[c][a]
                        + scalar
                        * Fraction(1, 2)
                        * (
                            _metric_component(metric, a, c)
                            * _metric_component(metric, d, b)
                            - _metric_component(metric, a, d)
                            * _metric_component(metric, c, b)
                        )
                    )
    return _freeze_tensor4(tensor)


def _raise_hessian(metric_diagonal: Vector, hessian_lower: Matrix) -> Matrix:
    inverse = _inverse_metric_diagonal(metric_diagonal)
    return tuple(
        tuple(inverse[a] * inverse[b] * hessian_lower[a][b] for b in range(DIMENSION))
        for a in range(DIMENSION)
    )


def _box(metric_diagonal: Vector, hessian_lower: Matrix) -> Fraction:
    inverse = _inverse_metric_diagonal(metric_diagonal)
    return sum(inverse[index] * hessian_lower[index][index] for index in range(DIMENSION))


def compact_gauss_bonnet_metric_source(
    metric_diagonal: Vector,
    riemann_lower: Tensor4,
    hessian_f_lower: Matrix,
) -> Matrix:
    """Evaluate ``-8 P_acbd nabla^c nabla^d f`` exactly."""

    metric = _validate_metric_diagonal(metric_diagonal)
    riemann = _validate_tensor4("riemann_lower", riemann_lower)
    hessian = _validate_matrix("hessian_f_lower", hessian_f_lower)
    if not _matrix_is_symmetric(hessian):
        raise ValueError("hessian_f_lower must be symmetric")
    p_tensor = double_dual_p_tensor(metric, riemann)
    hessian_up = _raise_hessian(metric, hessian)
    return tuple(
        tuple(
            -8
            * sum(
                p_tensor[a][c][b][d] * hessian_up[c][d]
                for c in range(DIMENSION)
                for d in range(DIMENSION)
            )
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )


def expanded_gauss_bonnet_basis(
    metric_diagonal: Vector,
    riemann_lower: Tensor4,
    hessian_f_lower: Matrix,
) -> dict[str, Matrix]:
    """Return seven independently assembled curvature--Hessian basis tensors."""

    metric = _validate_metric_diagonal(metric_diagonal)
    riemann = _validate_tensor4("riemann_lower", riemann_lower)
    hessian = _validate_matrix("hessian_f_lower", hessian_f_lower)
    if not _matrix_is_symmetric(hessian):
        raise ValueError("hessian_f_lower must be symmetric")
    inverse = _inverse_metric_diagonal(metric)
    hessian_up = _raise_hessian(metric, hessian)
    ricci = ricci_tensor(metric, riemann)
    scalar = scalar_curvature(metric, ricci)
    box_f = _box(metric, hessian)
    ricci_hessian = sum(
        ricci[c][d] * hessian_up[c][d]
        for c in range(DIMENSION)
        for d in range(DIMENSION)
    )

    basis = {key: _zero_matrix() for key in EXPANDED_GB_COEFFICIENTS}
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            basis["riemann_hessian"][a][b] = sum(
                riemann[a][c][b][d] * hessian_up[c][d]
                for c in range(DIMENSION)
                for d in range(DIMENSION)
            )
            basis["metric_ricci_hessian"][a][b] = (
                _metric_component(metric, a, b) * ricci_hessian
            )
            basis["ricci_b_hessian_a"][a][b] = sum(
                ricci[b][c] * inverse[c] * hessian[a][c]
                for c in range(DIMENSION)
            )
            basis["ricci_a_hessian_b"][a][b] = sum(
                ricci[a][c] * inverse[c] * hessian[b][c]
                for c in range(DIMENSION)
            )
            basis["ricci_ab_box"][a][b] = ricci[a][b] * box_f
            basis["scalar_metric_box"][a][b] = (
                scalar * _metric_component(metric, a, b) * box_f
            )
            basis["scalar_hessian"][a][b] = scalar * hessian[a][b]
    return {key: _freeze_matrix(value) for key, value in basis.items()}


def contract_expanded_gauss_bonnet_basis(
    basis: Mapping[str, Matrix],
    coefficients: Mapping[str, Fraction | Integral] | None = None,
) -> Matrix:
    """Contract the seven-term basis with an exact coefficient map."""

    selected = EXPANDED_GB_COEFFICIENTS if coefficients is None else coefficients
    if set(selected) != set(EXPANDED_GB_COEFFICIENTS):
        raise ValueError("expanded Gauss-Bonnet coefficient keys differ")
    if set(basis) != set(EXPANDED_GB_COEFFICIENTS):
        raise ValueError("expanded Gauss-Bonnet basis keys differ")
    exact_coefficients = {
        key: _fraction(f"coefficients.{key}", selected[key]) for key in selected
    }
    validated_basis = {
        key: _validate_matrix(f"basis.{key}", basis[key]) for key in basis
    }
    return tuple(
        tuple(
            sum(
                exact_coefficients[key] * validated_basis[key][a][b]
                for key in EXPANDED_GB_COEFFICIENTS
            )
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )


def expanded_gauss_bonnet_metric_source(
    metric_diagonal: Vector,
    riemann_lower: Tensor4,
    hessian_f_lower: Matrix,
) -> Matrix:
    """Evaluate the seven-term expanded source without constructing ``P``."""

    return contract_expanded_gauss_bonnet_basis(
        expanded_gauss_bonnet_basis(
            metric_diagonal, riemann_lower, hessian_f_lower
        )
    )


def composite_hessian(
    scalar_field: Fraction | Integral,
    gradient_scalar_lower: Vector,
    hessian_scalar_lower: Matrix,
    first_derivative: Fraction | Integral,
    second_derivative: Fraction | Integral,
) -> Matrix:
    """Return ``nabla_a nabla_b q(phi)`` from ``q'`` and ``q''``."""

    _fraction("scalar_field", scalar_field)
    gradient = _validate_vector("gradient_scalar_lower", gradient_scalar_lower)
    hessian = _validate_matrix("hessian_scalar_lower", hessian_scalar_lower)
    if not _matrix_is_symmetric(hessian):
        raise ValueError("hessian_scalar_lower must be symmetric")
    first = _fraction("first_derivative", first_derivative)
    second = _fraction("second_derivative", second_derivative)
    return tuple(
        tuple(
            first * hessian[a][b] + second * gradient[a] * gradient[b]
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )


def compact_nonminimal_metric_source(
    metric_diagonal: Vector,
    hessian_f_lower: Matrix,
) -> Matrix:
    """Evaluate ``(g_ab Box - nabla_a nabla_b) F`` exactly."""

    metric = _validate_metric_diagonal(metric_diagonal)
    hessian = _validate_matrix("hessian_F_lower", hessian_f_lower)
    if not _matrix_is_symmetric(hessian):
        raise ValueError("hessian_F_lower must be symmetric")
    box_f = _box(metric, hessian)
    return tuple(
        tuple(
            _metric_component(metric, a, b) * box_f - hessian[a][b]
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )


def expanded_quadratic_nonminimal_metric_source(
    metric_diagonal: Vector,
    scalar_field: Fraction | Integral,
    gradient_scalar_lower: Vector,
    hessian_scalar_lower: Matrix,
    beta: Fraction | Integral,
) -> Matrix:
    """Evaluate the expanded source for ``F=M_Pl^2+beta*phi^2``."""

    metric = _validate_metric_diagonal(metric_diagonal)
    field = _fraction("scalar_field", scalar_field)
    gradient = _validate_vector("gradient_scalar_lower", gradient_scalar_lower)
    hessian = _validate_matrix("hessian_scalar_lower", hessian_scalar_lower)
    coupling = _fraction("beta", beta)
    if not _matrix_is_symmetric(hessian):
        raise ValueError("hessian_scalar_lower must be symmetric")
    inverse = _inverse_metric_diagonal(metric)
    box_phi = _box(metric, hessian)
    gradient_squared = sum(
        inverse[index] * gradient[index] * gradient[index]
        for index in range(DIMENSION)
    )
    return tuple(
        tuple(
            2
            * coupling
            * (
                _metric_component(metric, a, b)
                * (field * box_phi + gradient_squared)
                - field * hessian[a][b]
                - gradient[a] * gradient[b]
            )
            for b in range(DIMENSION)
        )
        for a in range(DIMENSION)
    )


def _branch_derivatives(
    parameters: ActionParameters, scalar_field: Fraction
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    beta = _parameter_fraction("ricci_coupling", parameters.ricci_coupling)
    f_first = Fraction(0)
    f_second = Fraction(0)
    if parameters.model_id is ModelID.SGB_L:
        f_first = _parameter_fraction(
            "linear_gb_coupling", parameters.linear_gb_coupling
        )
    elif parameters.model_id is ModelID.FGC_QR:
        eta = _parameter_fraction(
            "quadratic_gb_coupling", parameters.quadratic_gb_coupling
        )
        f_first = eta * scalar_field / 4
        f_second = eta / 4
    f_derivatives = (f_first, f_second)
    F_derivatives = (2 * beta * scalar_field, 2 * beta)
    return (*f_derivatives, *F_derivatives)


def branch_metric_source_identity(
    parameters: ActionParameters,
    fixture: AlgebraicCurvatureFixture,
) -> dict[str, Any]:
    """Evaluate both metric-source identities for one ACT1 branch and fixture."""

    f_first, f_second, F_first, F_second = _branch_derivatives(
        parameters, fixture.scalar_field
    )
    hessian_f = composite_hessian(
        fixture.scalar_field,
        fixture.gradient_scalar_lower,
        fixture.hessian_scalar_lower,
        f_first,
        f_second,
    )
    hessian_F = composite_hessian(
        fixture.scalar_field,
        fixture.gradient_scalar_lower,
        fixture.hessian_scalar_lower,
        F_first,
        F_second,
    )
    compact_gb = compact_gauss_bonnet_metric_source(
        fixture.metric_diagonal, fixture.riemann_lower, hessian_f
    )
    expanded_gb = expanded_gauss_bonnet_metric_source(
        fixture.metric_diagonal, fixture.riemann_lower, hessian_f
    )
    compact_F = compact_nonminimal_metric_source(
        fixture.metric_diagonal, hessian_F
    )
    beta = _parameter_fraction("ricci_coupling", parameters.ricci_coupling)
    expanded_F = expanded_quadratic_nonminimal_metric_source(
        fixture.metric_diagonal,
        fixture.scalar_field,
        fixture.gradient_scalar_lower,
        fixture.hessian_scalar_lower,
        beta,
    )
    return {
        "model_id": parameters.model_id.value,
        "f_prime": _fraction_text(f_first),
        "f_second": _fraction_text(f_second),
        "F_prime": _fraction_text(F_first),
        "F_second": _fraction_text(F_second),
        "gb_compact_expanded_max_exact_residual": _fraction_text(
            _max_abs_matrix(_subtract_matrices(compact_gb, expanded_gb))
        ),
        "nonminimal_compact_expanded_max_exact_residual": _fraction_text(
            _max_abs_matrix(_subtract_matrices(compact_F, expanded_F))
        ),
        "gb_source_nonzero": not _matrix_is_zero(compact_gb),
        "nonminimal_source_nonzero": not _matrix_is_zero(compact_F),
    }


def _riemann_symmetry_checks(riemann: Tensor4) -> dict[str, bool]:
    return {
        "riemann_first_pair_antisymmetric": all(
            riemann[a][b][c][d] == -riemann[b][a][c][d]
            for a in range(DIMENSION)
            for b in range(DIMENSION)
            for c in range(DIMENSION)
            for d in range(DIMENSION)
        ),
        "riemann_second_pair_antisymmetric": all(
            riemann[a][b][c][d] == -riemann[a][b][d][c]
            for a in range(DIMENSION)
            for b in range(DIMENSION)
            for c in range(DIMENSION)
            for d in range(DIMENSION)
        ),
        "riemann_pair_exchange_symmetric": all(
            riemann[a][b][c][d] == riemann[c][d][a][b]
            for a in range(DIMENSION)
            for b in range(DIMENSION)
            for c in range(DIMENSION)
            for d in range(DIMENSION)
        ),
        "riemann_first_bianchi": all(
            riemann[a][b][c][d]
            + riemann[a][c][d][b]
            + riemann[a][d][b][c]
            == 0
            for a in range(DIMENSION)
            for b in range(DIMENSION)
            for c in range(DIMENSION)
            for d in range(DIMENSION)
        ),
    }


def _p_trace_identity_residual(
    metric_diagonal: Vector, riemann_lower: Tensor4
) -> Matrix:
    inverse = _inverse_metric_diagonal(metric_diagonal)
    p_tensor = double_dual_p_tensor(metric_diagonal, riemann_lower)
    ricci = ricci_tensor(metric_diagonal, riemann_lower)
    scalar = scalar_curvature(metric_diagonal, ricci)
    return tuple(
        tuple(
            sum(
                inverse[index] * p_tensor[index][c][index][d]
                for index in range(DIMENSION)
            )
            + ricci[c][d]
            - Fraction(1, 2) * _metric_component(metric_diagonal, c, d) * scalar
            for d in range(DIMENSION)
        )
        for c in range(DIMENSION)
    )


def _trace_identity_residual(
    metric_diagonal: Vector, riemann_lower: Tensor4, hessian_f_lower: Matrix
) -> Fraction:
    inverse = _inverse_metric_diagonal(metric_diagonal)
    stress = compact_gauss_bonnet_metric_source(
        metric_diagonal, riemann_lower, hessian_f_lower
    )
    trace = sum(inverse[index] * stress[index][index] for index in range(DIMENSION))
    ricci = ricci_tensor(metric_diagonal, riemann_lower)
    scalar = scalar_curvature(metric_diagonal, ricci)
    hessian_up = _raise_hessian(metric_diagonal, hessian_f_lower)
    ricci_hessian = sum(
        ricci[c][d] * hessian_up[c][d]
        for c in range(DIMENSION)
        for d in range(DIMENSION)
    )
    expected = 8 * ricci_hessian - 4 * scalar * _box(
        metric_diagonal, hessian_f_lower
    )
    return trace - expected


def _noether_curvature_contraction_residual(
    metric_diagonal: Vector, riemann_lower: Tensor4
) -> Matrix:
    """Check ``P_acbd R^d_e{}^{ac} = -G g_be/4`` exactly."""

    inverse = _inverse_metric_diagonal(metric_diagonal)
    p_tensor = double_dual_p_tensor(metric_diagonal, riemann_lower)
    invariant = gauss_bonnet_invariant(metric_diagonal, riemann_lower)
    return tuple(
        tuple(
            sum(
                p_tensor[a][c][b][d]
                * inverse[a]
                * inverse[c]
                * inverse[d]
                * riemann_lower[d][e][a][c]
                for a in range(DIMENSION)
                for c in range(DIMENSION)
                for d in range(DIMENSION)
            )
            + Fraction(1, 4)
            * invariant
            * _metric_component(metric_diagonal, b, e)
            for e in range(DIMENSION)
        )
        for b in range(DIMENSION)
    )


def validate_algebraic_curvature_fixture(
    fixture: AlgebraicCurvatureFixture,
) -> dict[str, bool | str]:
    """Return exact structural and source-coverage checks for one fixture."""

    riemann = fixture.riemann_lower
    ricci = ricci_tensor(fixture.metric_diagonal, riemann)
    generic_source = compact_gauss_bonnet_metric_source(
        fixture.metric_diagonal, riemann, fixture.hessian_scalar_lower
    )
    checks: dict[str, bool | str] = _riemann_symmetry_checks(riemann)
    checks.update(
        {
            "metric_is_minkowski_diagonal": fixture.metric_diagonal
            == MINKOWSKI_DIAGONAL,
            "ricci_symmetric": _matrix_is_symmetric(ricci),
            "nonzero_curvature": not _tensor4_is_zero(riemann),
            "nonzero_gauss_bonnet_invariant": gauss_bonnet_invariant(
                fixture.metric_diagonal, riemann
            )
            != 0,
            "nonzero_scalar_field": fixture.scalar_field != 0,
            "nonzero_scalar_gradient": any(
                value != 0 for value in fixture.gradient_scalar_lower
            ),
            "nonzero_scalar_hessian": not _matrix_is_zero(
                fixture.hessian_scalar_lower
            ),
            "generic_gb_source_nonzero": not _matrix_is_zero(generic_source),
            "generic_gb_source_symmetric": _matrix_is_symmetric(generic_source),
            "p_trace_max_exact_residual": _fraction_text(
                _max_abs_matrix(
                    _p_trace_identity_residual(fixture.metric_diagonal, riemann)
                )
            ),
            "noether_curvature_contraction_max_exact_residual": _fraction_text(
                _max_abs_matrix(
                    _noether_curvature_contraction_residual(
                        fixture.metric_diagonal, riemann
                    )
                )
            ),
        }
    )
    return checks


def _flrw_riemann(hubble: Fraction, hubble_derivative: Fraction) -> Tensor4:
    tensor = _zero_tensor4()
    acceleration = hubble_derivative + hubble * hubble
    for i in range(1, DIMENSION):
        tensor[0][i][0][i] = -acceleration
        tensor[i][0][0][i] = acceleration
        tensor[0][i][i][0] = acceleration
        tensor[i][0][i][0] = -acceleration
    for i in range(1, DIMENSION):
        for j in range(1, DIMENSION):
            for k in range(1, DIMENSION):
                for l in range(1, DIMENSION):
                    delta_ik = Fraction(1) if i == k else Fraction(0)
                    delta_jl = Fraction(1) if j == l else Fraction(0)
                    delta_il = Fraction(1) if i == l else Fraction(0)
                    delta_jk = Fraction(1) if j == k else Fraction(0)
                    tensor[i][j][k][l] = hubble * hubble * (
                        delta_ik * delta_jl - delta_il * delta_jk
                    )
    return _freeze_tensor4(tensor)


def flrw_compact_gb_energy_density(
    hubble: Fraction | Integral,
    hubble_derivative: Fraction | Integral,
    f_time_derivative: Fraction | Integral,
    f_second_time_derivative: Fraction | Integral,
) -> Fraction:
    """Evaluate the compact Gauss--Bonnet source's FLRW ``00`` component."""

    h = _fraction("hubble", hubble)
    h_dot = _fraction("hubble_derivative", hubble_derivative)
    f_dot = _fraction("f_time_derivative", f_time_derivative)
    f_ddot = _fraction("f_second_time_derivative", f_second_time_derivative)
    hessian = _zero_matrix()
    hessian[0][0] = f_ddot
    for i in range(1, DIMENSION):
        hessian[i][i] = -h * f_dot
    source = compact_gauss_bonnet_metric_source(
        MINKOWSKI_DIAGONAL,
        _flrw_riemann(h, h_dot),
        _freeze_matrix(hessian),
    )
    return source[0][0]


def _flrw_reduced_gb_lagrangian_expression(
    lapse: Fraction | _ExactDual,
    scale_factor_velocity: Fraction,
    f_time_derivative: Fraction,
) -> Fraction | _ExactDual:
    """Evaluate ``-8 f_dot a_dot^3 N^-3`` on rationals or exact duals."""

    coefficient = -8 * f_time_derivative * scale_factor_velocity**3
    return (lapse**-3) * coefficient


def flrw_reduced_gb_lagrangian(
    lapse: Fraction | Integral,
    scale_factor_velocity: Fraction | Integral,
    f_time_derivative: Fraction | Integral,
) -> Fraction:
    """Return the boundary-reduced homogeneous Gauss--Bonnet Lagrangian."""

    n = _fraction("lapse", lapse)
    a_dot = _fraction("scale_factor_velocity", scale_factor_velocity)
    f_dot = _fraction("f_time_derivative", f_time_derivative)
    if n <= 0:
        raise ValueError("lapse must be positive")
    value = _flrw_reduced_gb_lagrangian_expression(n, a_dot, f_dot)
    if not isinstance(value, Fraction):
        raise AssertionError("rational lapse expression returned a dual value")
    return value


def flrw_reduced_gb_lapse_derivative(
    lapse: Fraction | Integral,
    scale_factor_velocity: Fraction | Integral,
    f_time_derivative: Fraction | Integral,
) -> Fraction:
    """Differentiate the reduced FLRW action exactly with respect to ``N``."""

    n = _fraction("lapse", lapse)
    a_dot = _fraction("scale_factor_velocity", scale_factor_velocity)
    f_dot = _fraction("f_time_derivative", f_time_derivative)
    if n <= 0:
        raise ValueError("lapse must be positive")
    differentiated = _flrw_reduced_gb_lagrangian_expression(
        _ExactDual(n, Fraction(1)), a_dot, f_dot
    )
    if not isinstance(differentiated, _ExactDual):
        raise AssertionError("dual lapse expression returned a rational value")
    return differentiated.tangent


def flrw_lapse_variation_gb_t00(
    lapse: Fraction | Integral,
    scale_factor: Fraction | Integral,
    scale_factor_velocity: Fraction | Integral,
    f_time_derivative: Fraction | Integral,
) -> Fraction:
    """Convert the exact reduced-action lapse derivative to ``T_GB_00``.

    Per unit comoving volume and up to a boundary term,
    ``sqrt(-g) f G = -8 f_dot a_dot^3/N^3``.  Varying ``N`` and then setting
    ``N=1`` gives ``T_00=-24 H^3 f_dot``. The derivative is obtained from the
    same reduced-action expression by exact dual-number differentiation.
    """

    n = _fraction("lapse", lapse)
    a = _fraction("scale_factor", scale_factor)
    a_dot = _fraction("scale_factor_velocity", scale_factor_velocity)
    f_dot = _fraction("f_time_derivative", f_time_derivative)
    if n <= 0:
        raise ValueError("lapse must be positive")
    if a <= 0:
        raise ValueError("scale_factor must be positive")
    lapse_derivative = flrw_reduced_gb_lapse_derivative(n, a_dot, f_dot)
    return -(n * n / (a**3)) * lapse_derivative


def required_var1_nonclaims() -> dict[str, bool]:
    """Return every physical conclusion withheld by the VAR1 certificate."""

    return {
        "tensor_cas_or_second_independent_functional_variation_performed": False,
        "differential_bianchi_or_full_divergence_fixture_evaluated": False,
        "boundary_or_junction_variation_completed": False,
        "coupled_metric_chi_background_verified": False,
        "spherical_field_equations_derived": False,
        "spherical_principal_symbol_derived": False,
        "nonlinear_strong_hyperbolicity_proven": False,
        "full_coupled_quadratic_action_diagonalized": False,
        "nonlinear_ghost_freedom_proven": False,
        "constraint_propagation_proven": False,
        "metric_null_defocusing_derived": False,
        "finite_collapse_solution_constructed": False,
        "singularity_resolution_proven": False,
        "child_spacetime_constructed": False,
        "shear_robustness_proven": False,
        "global_flux_entropy_closure_derived": False,
        "dark_matter_derived": False,
        "dark_energy_derived": False,
        "variable_speed_of_light_derived": False,
        "particle_spectrum_derived": False,
        "theory_of_everything_derived": False,
        "observational_confirmation_obtained": False,
    }


def metric_variation_certificate(
    models: Sequence[ActionParameters],
    fixtures: Sequence[AlgebraicCurvatureFixture],
) -> dict[str, Any]:
    """Build the exact VAR1 algebraic certificate body."""

    if [model.model_id for model in models] != list(ModelID):
        raise ValueError("models must retain the GR-0, SGB-L, FGC-QR order")
    if not fixtures:
        raise ValueError("at least one algebraic-curvature fixture is required")

    fixture_records: list[dict[str, Any]] = []
    residuals_by_identity: dict[str, list[Fraction]] = {
        "compact_expanded_gb_source": [],
        "quadratic_nonminimal_source": [],
        "p_trace": [],
        "trace": [],
        "noether_curvature_contraction": [],
        "flrw_lapse": [],
    }
    all_structural_checks = True
    all_gb_stresses_symmetric = True
    branch_coverage = {
        "GR-0": {"gb_zero": True, "nonminimal_zero": True},
        "SGB-L": {"gb_nonzero": False, "nonminimal_zero": True},
        "FGC-QR": {"gb_nonzero": False, "nonminimal_nonzero": False},
    }

    for fixture in fixtures:
        validation = validate_algebraic_curvature_fixture(fixture)
        structural_boolean_keys = (
            "metric_is_minkowski_diagonal",
            "riemann_first_pair_antisymmetric",
            "riemann_second_pair_antisymmetric",
            "riemann_pair_exchange_symmetric",
            "riemann_first_bianchi",
            "ricci_symmetric",
            "nonzero_curvature",
            "nonzero_gauss_bonnet_invariant",
            "nonzero_scalar_field",
            "nonzero_scalar_gradient",
            "nonzero_scalar_hessian",
            "generic_gb_source_nonzero",
            "generic_gb_source_symmetric",
        )
        all_structural_checks = all_structural_checks and all(
            validation[key] is True for key in structural_boolean_keys
        )
        residuals_by_identity["p_trace"].append(
            Fraction(validation["p_trace_max_exact_residual"])
        )
        residuals_by_identity["noether_curvature_contraction"].append(
            Fraction(
                validation[
                    "noether_curvature_contraction_max_exact_residual"
                ]
            )
        )

        generic_compact = compact_gauss_bonnet_metric_source(
            fixture.metric_diagonal,
            fixture.riemann_lower,
            fixture.hessian_scalar_lower,
        )
        generic_expanded = expanded_gauss_bonnet_metric_source(
            fixture.metric_diagonal,
            fixture.riemann_lower,
            fixture.hessian_scalar_lower,
        )
        generic_residual = _max_abs_matrix(
            _subtract_matrices(generic_compact, generic_expanded)
        )
        all_gb_stresses_symmetric = all_gb_stresses_symmetric and (
            _matrix_is_symmetric(generic_compact)
            and _matrix_is_symmetric(generic_expanded)
        )
        trace_residual = abs(
            _trace_identity_residual(
                fixture.metric_diagonal,
                fixture.riemann_lower,
                fixture.hessian_scalar_lower,
            )
        )
        residuals_by_identity["compact_expanded_gb_source"].append(
            generic_residual
        )
        residuals_by_identity["trace"].append(trace_residual)

        branches = [branch_metric_source_identity(model, fixture) for model in models]
        for branch in branches:
            residuals_by_identity["compact_expanded_gb_source"].append(
                Fraction(branch["gb_compact_expanded_max_exact_residual"])
            )
            residuals_by_identity["quadratic_nonminimal_source"].append(
                Fraction(
                    branch["nonminimal_compact_expanded_max_exact_residual"]
                )
            )
        gr, linear, quadratic = branches
        branch_coverage["GR-0"]["gb_zero"] &= not gr["gb_source_nonzero"]
        branch_coverage["GR-0"]["nonminimal_zero"] &= not gr[
            "nonminimal_source_nonzero"
        ]
        branch_coverage["SGB-L"]["gb_nonzero"] |= linear["gb_source_nonzero"]
        branch_coverage["SGB-L"]["nonminimal_zero"] &= not linear[
            "nonminimal_source_nonzero"
        ]
        branch_coverage["FGC-QR"]["gb_nonzero"] |= quadratic["gb_source_nonzero"]
        branch_coverage["FGC-QR"]["nonminimal_nonzero"] |= quadratic[
            "nonminimal_source_nonzero"
        ]

        fixture_records.append(
            {
                "fixture_id": fixture.fixture_id,
                "fixture_source_sha256": fixture.source_sha256(),
                "gauss_bonnet_invariant": _fraction_text(
                    gauss_bonnet_invariant(
                        fixture.metric_diagonal, fixture.riemann_lower
                    )
                ),
                "validation": validation,
                "generic_compact_expanded_max_exact_residual": _fraction_text(
                    generic_residual
                ),
                "trace_identity_exact_residual": _fraction_text(trace_residual),
                "branches": branches,
            }
        )

    flrw_values = {
        "lapse": Fraction(1),
        "scale_factor": Fraction(7, 5),
        "hubble": Fraction(3, 2),
        "hubble_derivative": Fraction(-2, 3),
        "f_time_derivative": Fraction(5, 7),
        "f_second_time_derivative": Fraction(-11, 13),
    }
    flrw_compact = flrw_compact_gb_energy_density(
        flrw_values["hubble"],
        flrw_values["hubble_derivative"],
        flrw_values["f_time_derivative"],
        flrw_values["f_second_time_derivative"],
    )
    scale_factor_velocity = (
        flrw_values["lapse"]
        * flrw_values["scale_factor"]
        * flrw_values["hubble"]
    )
    reduced_lagrangian = flrw_reduced_gb_lagrangian(
        flrw_values["lapse"],
        scale_factor_velocity,
        flrw_values["f_time_derivative"],
    )
    exact_lapse_derivative = flrw_reduced_gb_lapse_derivative(
        flrw_values["lapse"],
        scale_factor_velocity,
        flrw_values["f_time_derivative"],
    )
    flrw_lapse = flrw_lapse_variation_gb_t00(
        flrw_values["lapse"],
        flrw_values["scale_factor"],
        scale_factor_velocity,
        flrw_values["f_time_derivative"],
    )
    flrw_residual = abs(flrw_compact - flrw_lapse)
    residuals_by_identity["flrw_lapse"].append(flrw_residual)
    maximum_residuals_by_identity = {
        key: max(values, default=Fraction(0))
        for key, values in residuals_by_identity.items()
    }
    maximum_residual = max(
        maximum_residuals_by_identity.values(), default=Fraction(0)
    )
    all_coverage = all(
        value for branch in branch_coverage.values() for value in branch.values()
    )
    nonclaims = required_var1_nonclaims()

    return {
        "convention_locked_equations": {
            "riemann": "R^rho_sigma_mu_nu=partial_mu_Gamma^rho_nu_sigma-partial_nu_Gamma^rho_mu_sigma+Gamma^rho_mu_lambda*Gamma^lambda_nu_sigma-Gamma^rho_nu_lambda*Gamma^lambda_mu_sigma",
            "gauss_bonnet": "G=R_abcd*R^abcd-4*R_ab*R^ab+R^2",
            "antisymmetrization_weight": "one_half",
            "double_dual_p": "P_abcd=R_abcd-2*g_a[c*R_d]b+2*g_b[c*R_d]a+R*g_a[c*g_d]b",
            "stress_definition": "T_GB_ab=-2/sqrt(-g)*delta(int_sqrt(-g)*f*G)/delta_g^ab",
            "compact_gb_source": "T_GB_ab=-8*P_acbd*nabla^c*nabla^d*f",
            "expanded_gb_coefficients": {
                key: _fraction_text(value)
                for key, value in EXPANDED_GB_COEFFICIENTS.items()
            },
            "nonminimal_source": "(g_ab*Box-nabla_a*nabla_b)F",
            "full_metric_equation": "F*G_ab+(g_ab*Box-nabla_a*nabla_b)F=T_phi_ab+T_chi_ab-8*P_acbd*nabla^c*nabla^d*f",
        },
        "fixture_generator": {
            "algorithm": "lcg64_v1",
            "multiplier": LCG64_MULTIPLIER,
            "increment": LCG64_INCREMENT,
            "modulus": "2^64",
            "curvature_construction": "sum_of_two_kulkarni_nomizu_products",
        },
        "fixture_records": fixture_records,
        "branch_coverage": branch_coverage,
        "flrw_lapse_crosscheck": {
            "inputs": {
                key: _fraction_text(value) for key, value in flrw_values.items()
            },
            "scale_factor_velocity": _fraction_text(scale_factor_velocity),
            "boundary_reduced_lagrangian": _fraction_text(reduced_lagrangian),
            "exact_automatic_lapse_derivative": _fraction_text(
                exact_lapse_derivative
            ),
            "stress_conversion_factor": "-N^2/a^3",
            "compact_T_GB_00": _fraction_text(flrw_compact),
            "independent_lapse_variation_T_GB_00": _fraction_text(flrw_lapse),
            "expected_formula": "T_GB_00=-24*H^3*f_dot",
            "max_exact_residual": _fraction_text(flrw_residual),
        },
        "verified_exact_checks": {
            "all_algebraic_curvature_fixtures_valid": all_structural_checks,
            "gb_metric_source_symmetric": all_gb_stresses_symmetric,
            "compact_expanded_gb_source_identity": maximum_residuals_by_identity[
                "compact_expanded_gb_source"
            ]
            == 0,
            "quadratic_nonminimal_source_identity": maximum_residuals_by_identity[
                "quadratic_nonminimal_source"
            ]
            == 0,
            "p_trace_identity": maximum_residuals_by_identity["p_trace"] == 0,
            "trace_identity": maximum_residuals_by_identity["trace"] == 0,
            "noether_curvature_contraction_identity": maximum_residuals_by_identity[
                "noether_curvature_contraction"
            ]
            == 0,
            "flrw_lapse_sign_and_factor": maximum_residuals_by_identity[
                "flrw_lapse"
            ]
            == 0,
            "all_branch_source_coverage": all_coverage,
            "constant_f_topological_bulk_source_zero": _matrix_is_zero(
                compact_gauss_bonnet_metric_source(
                    MINKOWSKI_DIAGONAL,
                    fixtures[0].riemann_lower,
                    _freeze_matrix(_zero_matrix()),
                )
            ),
            "all_var1_nonclaims_false": all(
                value is False for value in nonclaims.values()
            ),
        },
        "aggregate": {
            "fixture_count": len(fixtures),
            "maximum_exact_residual": _fraction_text(maximum_residual),
            "maximum_exact_residual_by_identity": {
                key: _fraction_text(value)
                for key, value in maximum_residuals_by_identity.items()
            },
            "all_required_coverage": all_coverage,
        },
        "nonclaims": nonclaims,
    }
