"""Exact parametric Krawczyk/Banach inclusion certificates.

This bounded utility consumes *already valid* interval enclosures for a
residual at an exact variable centre and for its variable Jacobian on a box.
It does not evaluate a residual, establish differentiability, or infer the
meaning of caller-supplied parameters.  Conditional on those inputs enclosing
a ``C1`` map ``F(x, p)``, strict inclusion and a strict infinity-norm
contraction prove one root in the declared variable box for each admitted
parameter point.
"""

from __future__ import annotations

from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping, Sequence, TypeAlias

from .exact_interval import (
    Interval,
    coerce_interval,
    interval_matrix_multiply,
    interval_matrix_subtract,
    interval_matrix_vector,
    point_matrix_interval_multiply,
)


Rational: TypeAlias = Fraction | Integral
IntervalMatrix: TypeAlias = tuple[tuple[Interval, ...], ...]
IntervalVector: TypeAlias = tuple[Interval, ...]
PointMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]
MAXIMUM_DIMENSION = 16


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Fraction(int(value))


def _point_matrix(name: str, value: Sequence[Sequence[Rational]]) -> PointMatrix:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{name} must be a nonempty matrix")
    rows: list[tuple[Fraction, ...]] = []
    width: int | None = None
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"{name}[{row_index}] must be a nonempty row")
        converted = tuple(_fraction(f"{name}[{row_index}] entry", item) for item in row)
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ValueError(f"{name} must be rectangular")
        rows.append(converted)
    if len(rows) != width:
        raise ValueError(f"{name} must be square")
    if not 0 < len(rows) <= MAXIMUM_DIMENSION:
        raise ValueError(f"{name} dimension must lie in 1..{MAXIMUM_DIMENSION}")
    return tuple(rows)


def _interval_vector(name: str, value: Sequence[object], *, dimension: int) -> IntervalVector:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a vector")
    converted = tuple(coerce_interval(item) for item in value)
    if len(converted) != dimension:
        raise ValueError(f"{name} has the wrong dimension")
    return converted


def _interval_matrix(name: str, value: Sequence[Sequence[object]], *, dimension: int) -> IntervalMatrix:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != dimension:
        raise ValueError(f"{name} must be a {dimension}-by-{dimension} matrix")
    rows: list[tuple[Interval, ...]] = []
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            raise ValueError(f"{name}[{row_index}] must be a matrix row")
        converted = tuple(coerce_interval(item) for item in row)
        if len(converted) != dimension:
            raise ValueError(f"{name} must be a {dimension}-by-{dimension} matrix")
        rows.append(converted)
    return tuple(rows)


def _identity(dimension: int) -> IntervalMatrix:
    return tuple(
        tuple(Interval.singleton(int(row == column)) for column in range(dimension))
        for row in range(dimension)
    )


def _row_sum_bound(matrix: IntervalMatrix) -> Fraction:
    return max(sum(entry.abs_upper() for entry in row) for row in matrix)


def _point_infinity_norm(matrix: PointMatrix) -> Fraction:
    return max(sum(abs(entry) for entry in row) for row in matrix)


def parametric_krawczyk_inclusion(
    *,
    center_inverse: Sequence[Sequence[Rational]],
    residual_at_center: Sequence[object],
    jacobian_box: Sequence[Sequence[object]],
    displacement_box: Sequence[object],
) -> Mapping[str, Any]:
    """Certify a componentwise parametric Krawczyk inclusion.

    ``center_inverse`` is an exact point preconditioner ``C``.  The supplied
    residual enclosure is ``F(x0, p)`` and ``jacobian_box`` encloses
    ``D_x F(x, p)`` over every ``x=x0+d`` with ``d`` in ``displacement_box``
    and every admitted parameter ``p``.  Each displacement must put zero
    strictly inside its interval, making ``x0`` an interior point of the
    declared variable box.

    The returned image is in displacement coordinates:
    ``-C F(x0,p) + (I-C D_xF(X,p)) displacement_box``.  If the reported
    contraction is strictly below one and this image is strictly inside the
    displacement box, the Banach/Krawczyk conclusion applies *provided the
    caller's residual/Jacobian enclosures are valid for a C1 map*.
    """

    c = _point_matrix("center inverse", center_inverse)
    dimension = len(c)
    residual = _interval_vector(
        "residual at center", residual_at_center, dimension=dimension
    )
    jacobian = _interval_matrix("Jacobian box", jacobian_box, dimension=dimension)
    displacement = _interval_vector(
        "displacement box", displacement_box, dimension=dimension
    )
    if any(not (entry.lower < 0 < entry.upper) for entry in displacement):
        raise ValueError("every displacement interval must contain zero strictly in its interior")

    preconditioned_jacobian = point_matrix_interval_multiply(c, jacobian)
    contraction_matrix = interval_matrix_subtract(
        _identity(dimension), preconditioned_jacobian
    )
    rho = _row_sum_bound(contraction_matrix)
    if rho >= 1:
        raise ValueError("Krawczyk contraction infinity-norm bound must be below one")
    center_correction = tuple(
        -entry for entry in interval_matrix_vector(c, residual)
    )
    contraction_image = interval_matrix_vector(contraction_matrix, displacement)
    image = tuple(
        center_correction[index] + contraction_image[index]
        for index in range(dimension)
    )
    margins = tuple(
        min(
            image[index].lower - displacement[index].lower,
            displacement[index].upper - image[index].upper,
        )
        for index in range(dimension)
    )
    strict_inclusion = all(image[index].strictly_inside(displacement[index]) for index in range(dimension))
    if not strict_inclusion or min(margins) <= 0:
        raise ValueError("Krawczyk image is not strictly inside the displacement box")
    inverse_norm = _point_infinity_norm(c)
    return {
        "dimension": dimension,
        "norm": "unscaled_componentwise_infinity_norm",
        "center_inverse": c,
        "preconditioned_jacobian_box": preconditioned_jacobian,
        "identity_minus_preconditioned_jacobian_box": contraction_matrix,
        "rho_infinity_upper_bound": rho,
        "rho_strictly_below_one": True,
        "center_correction_box": center_correction,
        "displacement_box": displacement,
        "linear_contraction_image_box": contraction_image,
        "krawczyk_displacement_image_box": image,
        "strict_componentwise_inclusion_margins": margins,
        "minimum_strict_componentwise_inclusion_margin": min(margins),
        "krawczyk_image_strictly_inside_displacement_box": True,
        "uniform_inverse_infinity_norm_upper_bound": inverse_norm / (1 - rho),
        "conditional_theorem": (
            "For every admitted parameter, if the supplied residual and Jacobian "
            "enclosures are valid for a C1 map on the declared variable box, "
            "there is a unique variable root in that box."
        ),
        "input_enclosure_proven_here": False,
    }


# Compact alias for artifact adapters.
krawczyk_inclusion = parametric_krawczyk_inclusion
