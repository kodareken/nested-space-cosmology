"""Small fail-closed exact rational polynomial and matrix primitives.

SYM1 uses these helpers for finite, convention-frozen algebra only.  Values are
always :class:`fractions.Fraction`; no floating-point coercion, numerical pivot
tolerance, pseudo-inverse, or approximate polynomial division is available.
Polynomials are immutable coefficient tuples in low-degree-first order.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from numbers import Integral
from typing import Iterable, Sequence, TypeAlias


Rational: TypeAlias = Fraction | Integral
Polynomial: TypeAlias = tuple[Fraction, ...]
Matrix: TypeAlias = tuple[tuple[Fraction, ...], ...]
PolynomialMatrix: TypeAlias = tuple[tuple[Polynomial, ...], ...]

ZERO: Polynomial = (Fraction(0),)
ONE: Polynomial = (Fraction(1),)


def _fraction(name: str, value: Rational) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Fraction(int(value))


def polynomial(coefficients: Iterable[Rational]) -> Polynomial:
    """Return the canonical low-degree-first exact polynomial."""

    result = tuple(_fraction("polynomial coefficient", value) for value in coefficients)
    if not result:
        raise ValueError("polynomial must contain at least one coefficient")
    index = len(result) - 1
    while index > 0 and result[index] == 0:
        index -= 1
    return result[: index + 1]


def degree(value: Sequence[Rational]) -> int:
    """Return the degree of a canonicalized polynomial."""

    return len(polynomial(value)) - 1


def add(left: Sequence[Rational], right: Sequence[Rational]) -> Polynomial:
    """Return ``left + right`` exactly."""

    a, b = polynomial(left), polynomial(right)
    length = max(len(a), len(b))
    return polynomial(
        (a[index] if index < len(a) else Fraction(0))
        + (b[index] if index < len(b) else Fraction(0))
        for index in range(length)
    )


def subtract(left: Sequence[Rational], right: Sequence[Rational]) -> Polynomial:
    """Return ``left - right`` exactly."""

    return add(left, scale(right, Fraction(-1)))


def scale(value: Sequence[Rational], factor: Rational) -> Polynomial:
    """Return ``factor * value`` exactly."""

    multiplier = _fraction("polynomial scale", factor)
    return polynomial(coefficient * multiplier for coefficient in polynomial(value))


def multiply(left: Sequence[Rational], right: Sequence[Rational]) -> Polynomial:
    """Return the exact polynomial product."""

    a, b = polynomial(left), polynomial(right)
    output = [Fraction(0)] * (len(a) + len(b) - 1)
    for left_index, left_value in enumerate(a):
        for right_index, right_value in enumerate(b):
            output[left_index + right_index] += left_value * right_value
    return polynomial(output)


def evaluate(value: Sequence[Rational], argument: Rational) -> Fraction:
    """Evaluate a polynomial by Horner's rule with exact arithmetic."""

    x = _fraction("polynomial argument", argument)
    accumulator = Fraction(0)
    for coefficient in reversed(polynomial(value)):
        accumulator = accumulator * x + coefficient
    return accumulator


def exact_division(dividend: Sequence[Rational], divisor: Sequence[Rational]) -> Polynomial:
    """Return an exact quotient or raise when the polynomial remainder is nonzero."""

    numerator, denominator = polynomial(dividend), polynomial(divisor)
    if denominator == ZERO:
        raise ZeroDivisionError("polynomial divisor must be nonzero")
    if len(numerator) < len(denominator):
        if numerator == ZERO:
            return ZERO
        raise ValueError("polynomial division has a nonzero remainder")
    remainder = list(numerator)
    quotient = [Fraction(0)] * (len(numerator) - len(denominator) + 1)
    while len(remainder) >= len(denominator):
        leading = remainder[-1] / denominator[-1]
        offset = len(remainder) - len(denominator)
        quotient[offset] = leading
        for index, coefficient in enumerate(denominator):
            remainder[index + offset] -= leading * coefficient
        while remainder and remainder[-1] == 0:
            remainder.pop()
    if remainder:
        raise ValueError("polynomial division has a nonzero remainder")
    return polynomial(quotient)


def _matrix(name: str, value: Sequence[Sequence[Rational]], *, square: bool = False) -> Matrix:
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
    return tuple(rows)


def transpose(value: Sequence[Sequence[Rational]]) -> Matrix:
    """Transpose a nonempty rectangular exact matrix."""

    matrix = _matrix("matrix", value)
    return tuple(tuple(matrix[row][column] for row in range(len(matrix))) for column in range(len(matrix[0])))


def matrix_multiply(left: Sequence[Sequence[Rational]], right: Sequence[Sequence[Rational]]) -> Matrix:
    """Multiply exact rectangular matrices with strict shape validation."""

    a, b = _matrix("left matrix", left), _matrix("right matrix", right)
    if len(a[0]) != len(b):
        raise ValueError("matrix multiplication dimensions do not agree")
    return tuple(
        tuple(
            sum(a[row][inner] * b[inner][column] for inner in range(len(b)))
            for column in range(len(b[0]))
        )
        for row in range(len(a))
    )


def determinant(value: Sequence[Sequence[Rational]]) -> Fraction:
    """Return the exact determinant by fraction-preserving elimination."""

    matrix = [list(row) for row in _matrix("matrix", value, square=True)]
    size = len(matrix)
    sign = Fraction(1)
    determinant_value = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size) if matrix[row][column] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            sign = -sign
        pivot_value = matrix[column][column]
        determinant_value *= pivot_value
        for row in range(column + 1, size):
            factor = matrix[row][column] / pivot_value
            for inner in range(column + 1, size):
                matrix[row][inner] -= factor * matrix[column][inner]
            matrix[row][column] = Fraction(0)
    return sign * determinant_value


def rank(value: Sequence[Sequence[Rational]]) -> int:
    """Return exact row rank without a numerical tolerance."""

    matrix = [list(row) for row in _matrix("matrix", value)]
    rows, columns = len(matrix), len(matrix[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows) if matrix[row][column] != 0), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / pivot_value for entry in matrix[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = matrix[row][column]
            if factor:
                matrix[row] = [
                    entry - factor * pivot_entry
                    for entry, pivot_entry in zip(matrix[row], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def inverse(value: Sequence[Sequence[Rational]]) -> Matrix:
    """Return the exact inverse or reject a singular matrix."""

    original = _matrix("matrix", value, square=True)
    size = len(original)
    augmented = [
        list(row) + [Fraction(int(row_index == column)) for column in range(size)]
        for row_index, row in enumerate(original)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column] != 0), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        augmented[column] = [entry / pivot_value for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    entry - factor * pivot_entry
                    for entry, pivot_entry in zip(augmented[row], augmented[column])
                ]
    return tuple(tuple(row[size:]) for row in augmented)


def adjugate_2x2(value: Sequence[Sequence[Rational]]) -> Matrix:
    """Return the exact adjugate of a 2x2 matrix without inversion."""

    matrix = _matrix("2x2 matrix", value, square=True)
    if len(matrix) != 2:
        raise ValueError("adjugate_2x2 requires a 2x2 matrix")
    return ((matrix[1][1], -matrix[0][1]), (-matrix[1][0], matrix[0][0]))


def schur_numerator_2x2(value: Sequence[Sequence[Rational]]) -> Fraction:
    """Return ``a*d-b*c``, the numerator of the scalar 2x2 Schur complement."""

    matrix = _matrix("2x2 matrix", value, square=True)
    if len(matrix) != 2:
        raise ValueError("schur_numerator_2x2 requires a 2x2 matrix")
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def _polynomial_matrix(value: Sequence[Sequence[Sequence[Rational]]]) -> PolynomialMatrix:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError("polynomial matrix must be nonempty")
    rows: list[tuple[Polynomial, ...]] = []
    width: int | None = None
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"polynomial matrix row {row_index} must be nonempty")
        converted = tuple(polynomial(entry) for entry in row)
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ValueError("polynomial matrix must be rectangular")
        rows.append(converted)
    if len(rows) != width:
        raise ValueError("polynomial matrix must be square")
    if len(rows) > 4:
        raise ValueError("polynomial matrix determinant supports at most 4x4")
    return tuple(rows)


def _permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def polynomial_matrix_determinant(
    value: Sequence[Sequence[Sequence[Rational]]],
) -> Polynomial:
    """Return an exact polynomial determinant for a square matrix up to 4x4.

    The Leibniz expansion intentionally avoids polynomial pivot division, so a
    vanishing leading entry never causes an approximate or singular fallback.
    """

    matrix = _polynomial_matrix(value)
    output = ZERO
    for permutation in permutations(range(len(matrix))):
        term = ONE
        for row, column in enumerate(permutation):
            term = multiply(term, matrix[row][column])
        output = add(output, scale(term, _permutation_sign(permutation)))
    return output


def polynomial_matrix_determinant_bareiss(
    value: Sequence[Sequence[Sequence[Rational]]],
) -> Polynomial:
    """Return an exact determinant by fraction-free Bareiss elimination.

    Unlike :func:`polynomial_matrix_determinant`, this route is suitable for
    the six-by-six modified-harmonic symbol.  Every intermediate division must
    be exact in ``Q[c]``; a violated Bareiss invariant raises instead of
    introducing a rational-function or approximate fallback.
    """

    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError("polynomial matrix must be nonempty")
    rows: list[list[Polynomial]] = []
    width: int | None = None
    for row_index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"polynomial matrix row {row_index} must be nonempty")
        converted = [polynomial(entry) for entry in row]
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ValueError("polynomial matrix must be rectangular")
        rows.append(converted)
    if len(rows) != width:
        raise ValueError("polynomial matrix must be square")
    if len(rows) > 16:
        raise ValueError("Bareiss polynomial determinant supports at most 16x16")
    if len(rows) == 1:
        return rows[0][0]

    sign = Fraction(1)
    previous_pivot = ONE
    size = len(rows)
    for pivot_index in range(size - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_index, size)
                if rows[row][pivot_index] != ZERO
            ),
            None,
        )
        if pivot_row is None:
            return ZERO
        if pivot_row != pivot_index:
            rows[pivot_index], rows[pivot_row] = rows[pivot_row], rows[pivot_index]
            sign = -sign
        pivot = rows[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = subtract(
                    multiply(pivot, rows[row][column]),
                    multiply(rows[row][pivot_index], rows[pivot_index][column]),
                )
                rows[row][column] = (
                    numerator
                    if pivot_index == 0
                    else exact_division(numerator, previous_pivot)
                )
            rows[row][pivot_index] = ZERO
        previous_pivot = pivot
    return scale(rows[-1][-1], sign)


# Compact SYM1-facing names.  The descriptive functions above remain the
# implementation authority and are intentionally aliases rather than wrappers.
poly = polynomial
poly_add = add
poly_sub = subtract
poly_mul = multiply
poly_scale = scale
poly_eval = evaluate
poly_div_exact = exact_division
matrix_mul = matrix_multiply
matrix_det = determinant
matrix_rank = rank
matrix_inverse = inverse
polynomial_matrix_det = polynomial_matrix_determinant
polynomial_matrix_det_bareiss = polynomial_matrix_determinant_bareiss
