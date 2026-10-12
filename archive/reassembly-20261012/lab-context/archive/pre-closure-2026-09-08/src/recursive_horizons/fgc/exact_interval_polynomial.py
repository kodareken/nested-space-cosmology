"""Exact interval-polynomial helpers for compact characteristic certificates.

Coefficients are closed rational intervals in ascending degree order.  The
module never treats a zero-containing interval as zero and never divides by a
leading coefficient that may vanish.  Polynomial matrix determinants use a
bounded permutation expansion, avoiding interval Gaussian pivots.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from typing import Any, Sequence, TypeAlias

from .exact_interval import Interval, coerce_interval, interval


IntervalPolynomial: TypeAlias = tuple[Interval, ...]
MAXIMUM_DETERMINANT_SIZE = 6


def _polynomial(value: Sequence[object]) -> IntervalPolynomial:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError("interval polynomial must be a nonempty coefficient sequence")
    coefficients = [coerce_interval(entry) for entry in value]
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients)


def polynomial_add(left: Sequence[object], right: Sequence[object]) -> IntervalPolynomial:
    a, b = _polynomial(left), _polynomial(right)
    size = max(len(a), len(b))
    return _polynomial(
        tuple(
            (a[index] if index < len(a) else interval(0))
            + (b[index] if index < len(b) else interval(0))
            for index in range(size)
        )
    )


def polynomial_negate(value: Sequence[object]) -> IntervalPolynomial:
    return _polynomial(tuple(-entry for entry in _polynomial(value)))


def polynomial_subtract(left: Sequence[object], right: Sequence[object]) -> IntervalPolynomial:
    return polynomial_add(left, polynomial_negate(right))


def polynomial_scale(value: Sequence[object], scalar: object) -> IntervalPolynomial:
    factor = coerce_interval(scalar)
    return _polynomial(tuple(entry * factor for entry in _polynomial(value)))


def polynomial_multiply(left: Sequence[object], right: Sequence[object]) -> IntervalPolynomial:
    a, b = _polynomial(left), _polynomial(right)
    output = [interval(0) for _ in range(len(a) + len(b) - 1)]
    for left_degree, left_coefficient in enumerate(a):
        for right_degree, right_coefficient in enumerate(b):
            output[left_degree + right_degree] = (
                output[left_degree + right_degree]
                + left_coefficient * right_coefficient
            )
    return _polynomial(output)


def polynomial_evaluate(value: Sequence[object], argument: object) -> Interval:
    coefficients = _polynomial(value)
    point = coerce_interval(argument)
    result = interval(0)
    for coefficient in reversed(coefficients):
        result = result * point + coefficient
    return result


def polynomial_derivative(value: Sequence[object]) -> IntervalPolynomial:
    coefficients = _polynomial(value)
    if len(coefficients) == 1:
        return (interval(0),)
    return _polynomial(
        tuple(degree * coefficients[degree] for degree in range(1, len(coefficients)))
    )


def monic_enclosure(value: Sequence[object]) -> IntervalPolynomial:
    """Enclose normalization by the possibly varying nonzero leading term."""

    coefficients = _polynomial(value)
    leading = coefficients[-1]
    if leading.contains_zero():
        raise ValueError("monic interval normalization requires nonzero leading coefficient")
    return _polynomial(tuple(coefficient / leading for coefficient in coefficients))


def divide_by_monic_enclosure(
    dividend: Sequence[object], divisor: Sequence[object]
) -> dict[str, IntervalPolynomial]:
    """Apply long division when the divisor leading coefficient is exactly one.

    The quotient encloses the quotient of every exactly divisible point
    polynomial represented by the inputs.  A non-singleton remainder is only
    an interval dependency remainder; callers need a separate algebraic
    factorization identity before discarding it.
    """

    numerator = list(_polynomial(dividend))
    denominator = _polynomial(divisor)
    if denominator[-1] != 1:
        raise ValueError("interval polynomial divisor must be exactly monic")
    if len(numerator) < len(denominator):
        return {"quotient": (interval(0),), "remainder": tuple(numerator)}
    quotient = [interval(0) for _ in range(len(numerator) - len(denominator) + 1)]
    for degree in range(len(numerator) - 1, len(denominator) - 2, -1):
        factor = numerator[degree]
        quotient_degree = degree - (len(denominator) - 1)
        quotient[quotient_degree] = factor
        for divisor_degree, coefficient in enumerate(denominator):
            numerator[quotient_degree + divisor_degree] = (
                numerator[quotient_degree + divisor_degree] - factor * coefficient
            )
    remainder = numerator[: len(denominator) - 1] or [interval(0)]
    return {"quotient": _polynomial(quotient), "remainder": _polynomial(remainder)}


def polynomial_matrix_determinant(
    matrix: Sequence[Sequence[Sequence[object]]],
) -> IntervalPolynomial:
    """Return a bounded exact interval enclosure of a polynomial determinant."""

    if not isinstance(matrix, Sequence) or isinstance(matrix, (str, bytes)) or not matrix:
        raise ValueError("polynomial matrix must be nonempty")
    size = len(matrix)
    if size > MAXIMUM_DETERMINANT_SIZE or any(len(row) != size for row in matrix):
        raise ValueError(
            f"polynomial determinant requires a square matrix up to {MAXIMUM_DETERMINANT_SIZE}x{MAXIMUM_DETERMINANT_SIZE}"
        )
    converted = tuple(tuple(_polynomial(entry) for entry in row) for row in matrix)
    result: IntervalPolynomial = (interval(0),)
    for ordering in permutations(range(size)):
        inversions = sum(
            ordering[left] > ordering[right]
            for left in range(size)
            for right in range(left + 1, size)
        )
        term: IntervalPolynomial = (interval(-1 if inversions % 2 else 1),)
        for row, column in enumerate(ordering):
            term = polynomial_multiply(term, converted[row][column])
        result = polynomial_add(result, term)
    return result


def quadratic_discriminant(value: Sequence[object]) -> Interval:
    coefficients = _polynomial(value)
    if len(coefficients) != 3 or coefficients[2].contains_zero():
        raise ValueError("interval characteristic factor must remain genuinely quadratic")
    constant, linear, quadratic = coefficients
    return linear**2 - 4 * quadratic * constant


def quadratic_root_bracket(
    value: Sequence[object], root_interval: Interval
) -> dict[str, Any]:
    """Prove one root exists uniquely in a rational interval for every factor.

    Uniform opposite signs at the two endpoints give existence by the
    intermediate value theorem.  A derivative interval of strict one-sided
    sign gives uniqueness for every point polynomial in the coefficient box.
    """

    coefficients = _polynomial(value)
    if len(coefficients) != 3:
        raise ValueError("root bracketing currently requires a quadratic")
    if not isinstance(root_interval, Interval) or root_interval.is_singleton():
        raise ValueError("root bracket must be a nondegenerate exact interval")
    left = polynomial_evaluate(coefficients, interval(root_interval.lower))
    right = polynomial_evaluate(coefficients, interval(root_interval.upper))
    opposite = (
        left.strictly_negative() and right.strictly_positive()
    ) or (
        left.strictly_positive() and right.strictly_negative()
    )
    if not opposite:
        raise ValueError("quadratic endpoint signs do not uniformly bracket a root")
    derivative = polynomial_evaluate(
        polynomial_derivative(coefficients), root_interval
    )
    monotone = derivative.strictly_positive() or derivative.strictly_negative()
    if not monotone:
        raise ValueError("quadratic derivative is not uniformly one-sided on bracket")
    discriminant = quadratic_discriminant(coefficients)
    if not discriminant.strictly_positive():
        raise ValueError("quadratic discriminant is not uniformly positive")
    return {
        "root_interval": root_interval,
        "left_endpoint_value": left,
        "right_endpoint_value": right,
        "derivative_interval": derivative,
        "discriminant_interval": discriminant,
        "uniform_opposite_endpoint_signs": True,
        "uniform_strict_monotonicity": True,
        "exactly_one_real_root_for_every_enclosed_quadratic": True,
        "method": "uniform_endpoint_signs_plus_monotonicity",
    }
