"""Fail-closed exact rational closed intervals for bounded certificates.

The type in this module represents a closed real interval with exact
``Fraction`` endpoints.  It deliberately has no implicit truth value or
ordering: a sign or branch may be used only after one of the explicit
whole-domain predicates has established it over the *entire* interval.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from itertools import permutations
from typing import Any, Sequence, TypeAlias


Rational: TypeAlias = Fraction | Integral
IntervalMatrix: TypeAlias = tuple[tuple["Interval", ...], ...]
IntervalVector: TypeAlias = tuple["Interval", ...]


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Fraction(int(value))


@dataclass(frozen=True, slots=True, eq=False)
class Interval:
    """A finite closed interval with exact rational endpoints."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        lower = _fraction("interval lower endpoint", self.lower)
        upper = _fraction("interval upper endpoint", self.upper)
        if lower > upper:
            raise ValueError("interval lower endpoint must not exceed upper endpoint")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @classmethod
    def singleton(cls, value: Rational) -> "Interval":
        """Return the exact singleton interval containing ``value``."""

        rational = _fraction("interval singleton", value)
        return cls(rational, rational)

    @staticmethod
    def coerce(value: object) -> "Interval" | Any:
        if isinstance(value, Interval):
            return value
        if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
            return NotImplemented
        return Interval.singleton(value)

    def __bool__(self) -> bool:
        raise TypeError("interval truth value is ambiguous; use a whole-domain predicate")

    def _comparison(self, _other: object) -> Any:
        raise TypeError("interval ordering is ambiguous; use a whole-domain predicate")

    __lt__ = _comparison
    __le__ = _comparison
    __gt__ = _comparison
    __ge__ = _comparison

    def __eq__(self, other: object) -> bool:
        value = self.coerce(other)
        return False if value is NotImplemented else (self.lower, self.upper) == (value.lower, value.upper)

    def __hash__(self) -> int:
        return hash((self.lower, self.upper))

    def __neg__(self) -> "Interval":
        return Interval(-self.upper, -self.lower)

    def __add__(self, other: object) -> "Interval" | Any:
        value = self.coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return Interval(self.lower + value.lower, self.upper + value.upper)

    __radd__ = __add__

    def __sub__(self, other: object) -> "Interval" | Any:
        value = self.coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return Interval(self.lower - value.upper, self.upper - value.lower)

    def __rsub__(self, other: object) -> "Interval" | Any:
        value = self.coerce(other)
        return NotImplemented if value is NotImplemented else value - self

    def __mul__(self, other: object) -> "Interval" | Any:
        value = self.coerce(other)
        if value is NotImplemented:
            return NotImplemented
        corners = (
            self.lower * value.lower,
            self.lower * value.upper,
            self.upper * value.lower,
            self.upper * value.upper,
        )
        return Interval(min(corners), max(corners))

    __rmul__ = __mul__

    def reciprocal(self) -> "Interval":
        """Return the reciprocal interval, rejecting a zero-containing input."""

        if self.contains_zero():
            raise ZeroDivisionError("interval reciprocal denominator contains zero")
        return Interval(Fraction(1, 1) / self.upper, Fraction(1, 1) / self.lower)

    def __truediv__(self, other: object) -> "Interval" | Any:
        value = self.coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self * value.reciprocal()

    def __rtruediv__(self, other: object) -> "Interval" | Any:
        value = self.coerce(other)
        return NotImplemented if value is NotImplemented else value / self

    def __pow__(self, exponent: object) -> "Interval" | Any:
        if isinstance(exponent, bool) or not isinstance(exponent, Integral):
            return NotImplemented
        power = int(exponent)
        if power == 0:
            return Interval.singleton(1)
        if power < 0:
            return self.reciprocal() ** (-power)
        if power % 2 or self.lower >= 0:
            return Interval(self.lower**power, self.upper**power)
        if self.upper <= 0:
            return Interval(self.upper**power, self.lower**power)
        return Interval(Fraction(0), max((-self.lower) ** power, self.upper**power))

    def contains_zero(self) -> bool:
        return self.lower <= 0 <= self.upper

    def strictly_positive(self) -> bool:
        return self.lower > 0

    def nonnegative(self) -> bool:
        return self.lower >= 0

    def strictly_negative(self) -> bool:
        return self.upper < 0

    def nonpositive(self) -> bool:
        return self.upper <= 0

    def is_singleton(self) -> bool:
        return self.lower == self.upper

    def abs_upper(self) -> Fraction:
        return max(abs(self.lower), abs(self.upper))

    def width(self) -> Fraction:
        return self.upper - self.lower

    def midpoint(self) -> Fraction:
        return (self.lower + self.upper) / 2

    def radius(self) -> Fraction:
        return self.width() / 2

    def subset_of(self, other: object) -> bool:
        value = self.coerce(other)
        if value is NotImplemented:
            raise TypeError("subset target must be an interval or exact rational")
        return value.lower <= self.lower and self.upper <= value.upper

    def strictly_inside(self, other: object) -> bool:
        value = self.coerce(other)
        if value is NotImplemented:
            raise TypeError("interior target must be an interval or exact rational")
        return value.lower < self.lower and self.upper < value.upper


def interval(lower: Rational, upper: Rational | None = None) -> Interval:
    """Create an interval, treating one argument as an exact singleton."""

    return Interval.singleton(lower) if upper is None else Interval(lower, upper)


def coerce_interval(value: object) -> Interval:
    """Coerce an exact rational or interval, rejecting all other values."""

    result = Interval.coerce(value)
    if result is NotImplemented:
        raise TypeError("value must be an interval or exact rational")
    return result


def interval_matrix_multiply(
    left: Sequence[Sequence[object]], right: Sequence[Sequence[object]]
) -> IntervalMatrix:
    """Multiply nonempty rectangular interval matrices with exact enclosures."""

    a = _matrix("left matrix", left)
    b = _matrix("right matrix", right)
    if len(a[0]) != len(b):
        raise ValueError("interval matrix multiplication dimensions do not agree")
    return tuple(
        tuple(
            sum((a[row][inner] * b[inner][column] for inner in range(len(b))), Interval.singleton(0))
            for column in range(len(b[0]))
        )
        for row in range(len(a))
    )


def interval_matrix_add(
    left: Sequence[Sequence[object]], right: Sequence[Sequence[object]]
) -> IntervalMatrix:
    """Add same-shaped interval matrices entrywise."""

    a = _matrix("left matrix", left)
    b = _matrix("right matrix", right)
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        raise ValueError("interval matrix addition dimensions do not agree")
    return tuple(tuple(a[row][column] + b[row][column] for column in range(len(a[0]))) for row in range(len(a)))


def interval_matrix_subtract(
    left: Sequence[Sequence[object]], right: Sequence[Sequence[object]]
) -> IntervalMatrix:
    """Subtract same-shaped interval matrices entrywise."""

    a = _matrix("left matrix", left)
    b = _matrix("right matrix", right)
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        raise ValueError("interval matrix subtraction dimensions do not agree")
    return tuple(tuple(a[row][column] - b[row][column] for column in range(len(a[0]))) for row in range(len(a)))


def point_matrix_interval_multiply(
    left: Sequence[Sequence[Rational]], right: Sequence[Sequence[object]]
) -> IntervalMatrix:
    """Multiply an exact point matrix by an interval matrix."""

    return interval_matrix_multiply(_point_matrix("left point matrix", left), right)


def interval_matrix_point_multiply(
    left: Sequence[Sequence[object]], right: Sequence[Sequence[Rational]]
) -> IntervalMatrix:
    """Multiply an interval matrix by an exact point matrix."""

    return interval_matrix_multiply(left, _point_matrix("right point matrix", right))


def interval_matrix_determinant(value: Sequence[Sequence[object]]) -> Interval:
    """Enclose a small square determinant by the exact permutation formula.

    Avoiding Gaussian elimination here is intentional: a pivot interval that
    contains zero is not safely invertible.  The permutation expansion has no
    pivot choice and is exact interval arithmetic, but factorial work makes it
    a bounded helper rather than a large-matrix algorithm.
    """

    matrix = _matrix("matrix", value)
    size = len(matrix)
    if size != len(matrix[0]):
        raise ValueError("interval matrix determinant requires a square matrix")
    if size > 6:
        raise ValueError("interval matrix determinant supports at most 6x6 matrices")
    result = Interval.singleton(0)
    for ordering in permutations(range(size)):
        inversions = sum(
            ordering[left] > ordering[right]
            for left in range(size)
            for right in range(left + 1, size)
        )
        term = Interval.singleton(-1 if inversions % 2 else 1)
        for row, column in enumerate(ordering):
            term = term * matrix[row][column]
        result = result + term
    return result


def interval_matrix_vector(
    matrix: Sequence[Sequence[object]], vector: Sequence[object]
) -> IntervalVector:
    """Multiply a nonempty interval matrix by a matching nonempty vector."""

    a = _matrix("matrix", matrix)
    if not isinstance(vector, Sequence) or isinstance(vector, (str, bytes)) or not vector:
        raise ValueError("interval vector must be nonempty")
    b = tuple(coerce_interval(value) for value in vector)
    if len(a[0]) != len(b):
        raise ValueError("interval matrix-vector dimensions do not agree")
    return tuple(sum((entry * b[column] for column, entry in enumerate(row)), Interval.singleton(0)) for row in a)


def interval_matrix_infinity_row_sum_bound(matrix: Sequence[Sequence[object]]) -> Fraction:
    """Return the exact upper bound ``max_i sum_j sup(abs(A_ij))``."""

    value = _matrix("matrix", matrix)
    return max(sum(entry.abs_upper() for entry in row) for row in value)


def _matrix(name: str, value: Sequence[Sequence[object]]) -> IntervalMatrix:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{name} must be a nonempty matrix")
    rows: list[tuple[Interval, ...]] = []
    width: int | None = None
    for index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"{name}[{index}] must be a nonempty matrix row")
        converted = tuple(coerce_interval(entry) for entry in row)
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ValueError(f"{name} must be rectangular")
        rows.append(converted)
    return tuple(rows)


def _point_matrix(name: str, value: Sequence[Sequence[Rational]]) -> tuple[tuple[Fraction, ...], ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{name} must be a nonempty matrix")
    rows: list[tuple[Fraction, ...]] = []
    width: int | None = None
    for index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"{name}[{index}] must be a nonempty matrix row")
        converted = tuple(_fraction(f"{name}[{index}] entry", entry) for entry in row)
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ValueError(f"{name} must be rectangular")
        rows.append(converted)
    return tuple(rows)


# Compact aliases for config/reproducer consumers.
matrix_multiply = interval_matrix_multiply
matrix_vector = interval_matrix_vector
matrix_infinity_row_sum_bound = interval_matrix_infinity_row_sum_bound
matrix_add = interval_matrix_add
matrix_subtract = interval_matrix_subtract
matrix_determinant = interval_matrix_determinant
