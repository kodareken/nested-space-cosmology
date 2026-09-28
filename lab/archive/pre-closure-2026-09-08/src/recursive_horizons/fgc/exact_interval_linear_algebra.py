"""Bounded exact interval linear algebra for verified Neumann inverses.

This module deliberately keeps the large-matrix operation needed by UHYP1 out
of :mod:`exact_interval`: it certifies an inverse from an exact point
preconditioner and the interval Neumann criterion, rather than expanding a
large interval determinant.  Every endpoint is a ``Fraction`` and every
failure is explicit; no numerical pivot, tolerance, or floating-point
fallback exists.
"""

from __future__ import annotations

from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping, Sequence, TypeAlias

from .exact_interval import Interval, coerce_interval


Rational: TypeAlias = Fraction | Integral
IntervalMatrix: TypeAlias = tuple[tuple[Interval, ...], ...]
PointMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]
MAXIMUM_MATRIX_SIZE = 12


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Fraction(int(value))


def _interval_matrix(
    name: str, value: Sequence[Sequence[object]], *, square: bool = False
) -> IntervalMatrix:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{name} must be a nonempty matrix")
    rows: list[tuple[Interval, ...]] = []
    width: int | None = None
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"{name}[{row_index}] must be a nonempty matrix row")
        converted = tuple(coerce_interval(entry) for entry in row)
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ValueError(f"{name} must be rectangular")
        rows.append(converted)
    if square and len(rows) != width:
        raise ValueError(f"{name} must be square")
    if len(rows) > MAXIMUM_MATRIX_SIZE or width is not None and width > MAXIMUM_MATRIX_SIZE:
        raise ValueError(f"{name} supports at most {MAXIMUM_MATRIX_SIZE}x{MAXIMUM_MATRIX_SIZE}")
    return tuple(rows)


def _point_matrix(
    name: str, value: Sequence[Sequence[Rational]], *, square: bool = False
) -> PointMatrix:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{name} must be a nonempty matrix")
    rows: list[tuple[Fraction, ...]] = []
    width: int | None = None
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"{name}[{row_index}] must be a nonempty matrix row")
        converted = tuple(_fraction(f"{name}[{row_index}] entry", entry) for entry in row)
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ValueError(f"{name} must be rectangular")
        rows.append(converted)
    if square and len(rows) != width:
        raise ValueError(f"{name} must be square")
    if len(rows) > MAXIMUM_MATRIX_SIZE or width is not None and width > MAXIMUM_MATRIX_SIZE:
        raise ValueError(f"{name} supports at most {MAXIMUM_MATRIX_SIZE}x{MAXIMUM_MATRIX_SIZE}")
    return tuple(rows)


def identity(size: int) -> IntervalMatrix:
    """Return the exact interval identity of bounded positive dimension."""

    if isinstance(size, bool) or not isinstance(size, Integral) or not 0 < size <= MAXIMUM_MATRIX_SIZE:
        raise ValueError(f"identity size must lie in 1..{MAXIMUM_MATRIX_SIZE}")
    return tuple(
        tuple(Interval.singleton(int(row == column)) for column in range(int(size)))
        for row in range(int(size))
    )


def interval_matrix_multiply(
    left: Sequence[Sequence[object]], right: Sequence[Sequence[object]]
) -> IntervalMatrix:
    """Return the exact interval enclosure of a bounded matrix product."""

    a, b = _interval_matrix("left matrix", left), _interval_matrix("right matrix", right)
    if len(a[0]) != len(b):
        raise ValueError("interval matrix multiplication dimensions do not agree")
    return tuple(
        tuple(
            sum(
                (a[row][inner] * b[inner][column] for inner in range(len(b))),
                Interval.singleton(0),
            )
            for column in range(len(b[0]))
        )
        for row in range(len(a))
    )


def point_matrix_interval_multiply(
    left: Sequence[Sequence[Rational]], right: Sequence[Sequence[object]]
) -> IntervalMatrix:
    """Return ``left @ right`` for exact point ``left`` and interval ``right``."""

    point = _point_matrix("left point matrix", left)
    return interval_matrix_multiply(
        tuple(tuple(Interval.singleton(value) for value in row) for row in point), right
    )


def interval_matrix_subtract(
    left: Sequence[Sequence[object]], right: Sequence[Sequence[object]]
) -> IntervalMatrix:
    """Subtract equal-shaped bounded interval matrices."""

    a, b = _interval_matrix("left matrix", left), _interval_matrix("right matrix", right)
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        raise ValueError("interval matrix subtraction dimensions do not agree")
    return tuple(
        tuple(a[row][column] - b[row][column] for column in range(len(a[0])))
        for row in range(len(a))
    )


def interval_matrix_infinity_row_sum_bound(matrix: Sequence[Sequence[object]]) -> Fraction:
    """Return ``max_i sum_j sup(abs(matrix[i,j]))`` exactly."""

    converted = _interval_matrix("matrix", matrix)
    return max(sum(entry.abs_upper() for entry in row) for row in converted)


def center_preconditioned_neumann_inverse(
    matrix: Sequence[Sequence[object]],
    center_inverse: Sequence[Sequence[Rational]],
) -> Mapping[str, Any]:
    """Certify a bounded inverse enclosure by the left Neumann criterion.

    With exact point preconditioner ``C`` and interval matrix ``A``, define
    ``E = I - C A``.  When ``||E||_infinity = rho < 1``, every matrix in the
    interval enclosure is invertible and
    ``A^-1 = (I-E)^-1 C``.  The returned entrywise enclosure is
    ``C_ij +/- rho/(1-rho)*||C||_infinity``.  It is intentionally conservative
    but avoids determinant expansion and ambiguous interval pivot choices.
    """

    a = _interval_matrix("matrix", matrix, square=True)
    c = _point_matrix("center inverse", center_inverse, square=True)
    if len(a) != len(c):
        raise ValueError("matrix and center inverse dimensions do not agree")
    product = point_matrix_interval_multiply(c, a)
    residual = interval_matrix_subtract(identity(len(a)), product)
    rho = interval_matrix_infinity_row_sum_bound(residual)
    if rho >= 1:
        raise ValueError("Neumann inverse certificate requires rho < 1")
    center_norm = max(sum(abs(entry) for entry in row) for row in c)
    # The denominator is exact and positive only after the strict rho check.
    tail_bound = rho * center_norm / (1 - rho)
    enclosure: IntervalMatrix = tuple(
        tuple(Interval(entry - tail_bound, entry + tail_bound) for entry in row)
        for row in c
    )
    return {
        "dimension": len(a),
        "preconditioner": c,
        "preconditioned_residual": residual,
        "rho_infinity": rho,
        "center_inverse_infinity_norm": center_norm,
        "inverse_tail_entrywise_bound": tail_bound,
        "inverse_enclosure": enclosure,
        "rho_strictly_below_one": rho < 1,
        "every_enclosed_matrix_invertible": True,
        "method": "exact_left_center_preconditioned_neumann_series",
    }


# Compact names for the UHYP1 artifact adapter.
matrix_multiply = interval_matrix_multiply
point_times_interval = point_matrix_interval_multiply
matrix_subtract = interval_matrix_subtract
matrix_infinity_row_sum_bound = interval_matrix_infinity_row_sum_bound
neumann_inverse_certificate = center_preconditioned_neumann_inverse
