"""Exact quotient-ring charts for the physical MHG characteristic sectors.

MODE1 already supplies polynomial bases for the two repeated auxiliary
quadratics.  This module supplies the complementary one-dimensional physical
metric and regulator kernels at the solved COMP1 point.  A unit five-by-five
minor over ``Q[c]/(q)`` proves rank five at both real roots of ``q`` without
requiring either root to be rational.

Everything here is pointwise and exact.  It is a prerequisite for, not a
claim of, a uniform eigenframe or strong hyperbolicity.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from typing import Any, Sequence

from .exact_linear_algebra import Polynomial, poly_div_exact, poly_mul
from .modified_harmonic import (
    _physical_metric,
    inverse_metric_null_polynomial,
    modified_harmonic_symbol,
    monic,
    quadratic_companion_certificate,
    standard_first_order_principal_system,
)
from .modified_harmonic_modes import (
    Residue,
    monic_quadratic_factor,
    quadratic_quotient_parameters,
    quotient_add,
    quotient_inverse,
    quotient_multiply,
    quotient_negate,
    quotient_reduce,
    quotient_unit_norm,
)
from .spherical_reduction import BASE_FIELD_ORDER, SphericalState
from .spherical_symbol import metric_schur_data, quotient_symbol, radial_symbol


Q = Fraction
FIELD_COUNT = len(BASE_FIELD_ORDER)
PHYSICAL_REQUIRED_CHART = (5, 5)


def _sum_residues(values) -> Residue:
    total: Residue = (Q(0), Q(0))
    for value in values:
        total = quotient_add(total, value)
    return total


def _permutation_sign(permutation: Sequence[int]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def quotient_determinant(
    matrix: Sequence[Sequence[Residue]], *, b: Fraction, d: Fraction
) -> Residue:
    """Return a bounded exact determinant over a quadratic quotient ring."""

    size = len(matrix)
    if size == 0 or size > 5 or any(len(row) != size for row in matrix):
        raise ValueError("quotient determinant requires a square matrix of size 1..5")
    total: Residue = (Q(0), Q(0))
    for ordering in permutations(range(size)):
        term: Residue = (Q(_permutation_sign(ordering)), Q(0))
        for row, column in enumerate(ordering):
            term = quotient_multiply(term, matrix[row][column], b=b, d=d)
        total = quotient_add(total, term)
    return total


def quotient_solve(
    matrix: Sequence[Sequence[Residue]],
    rhs: Sequence[Residue],
    *,
    b: Fraction,
    d: Fraction,
) -> tuple[tuple[Residue, ...], Residue]:
    """Solve a bounded square system by exact quotient-ring Cramer rules."""

    size = len(matrix)
    if len(rhs) != size:
        raise ValueError("quotient solve right side has the wrong dimension")
    determinant = quotient_determinant(matrix, b=b, d=d)
    determinant_inverse = quotient_inverse(determinant, b=b, d=d)
    answer: list[Residue] = []
    for column in range(size):
        replacement = [list(row) for row in matrix]
        for row in range(size):
            replacement[row][column] = rhs[row]
        numerator = quotient_determinant(replacement, b=b, d=d)
        answer.append(
            quotient_multiply(numerator, determinant_inverse, b=b, d=d)
        )
    return tuple(answer), determinant


def quadratic_mode_chart(
    symbol: Sequence[Sequence[Sequence[Fraction]]],
    factor: Sequence[Fraction],
    *,
    required_chart: tuple[int, int] | None = None,
) -> dict[str, Any]:
    """Construct one normalized full-symbol kernel modulo a real quadratic.

    A chart is ``(omitted equation row, free field column)``.  All remaining
    rows and columns form a five-by-five pivot.  Its quotient-unit determinant
    is nonzero at both roots, and the normalized full residual proves that the
    kernel has dimension exactly one at each root.
    """

    if len(symbol) != FIELD_COUNT or any(len(row) != FIELD_COUNT for row in symbol):
        raise ValueError("physical mode chart requires a six-by-six symbol")
    monic_factor = monic_quadratic_factor(factor)
    companion = quadratic_companion_certificate(monic_factor)
    if companion["discriminant"] <= 0:
        raise ValueError("physical mode factor must have two distinct real roots")
    b, d = quadratic_quotient_parameters(monic_factor)
    quotient_symbol_matrix = tuple(
        tuple(quotient_reduce(entry, b=b, d=d) for entry in row)
        for row in symbol
    )
    if required_chart is not None:
        if (
            len(required_chart) != 2
            or any(isinstance(value, bool) or not isinstance(value, int) for value in required_chart)
            or any(value not in range(FIELD_COUNT) for value in required_chart)
        ):
            raise ValueError("required chart indices must lie in 0..5")
        choices = (required_chart,)
    else:
        choices = tuple(
            (omitted_row, free_column)
            for omitted_row in range(FIELD_COUNT)
            for free_column in range(FIELD_COUNT)
        )

    for omitted_row, free_column in choices:
        pivot_rows = tuple(row for row in range(FIELD_COUNT) if row != omitted_row)
        pivot_columns = tuple(
            column for column in range(FIELD_COUNT) if column != free_column
        )
        pivot = tuple(
            tuple(quotient_symbol_matrix[row][column] for column in pivot_columns)
            for row in pivot_rows
        )
        rhs = tuple(
            quotient_negate(quotient_symbol_matrix[row][free_column])
            for row in pivot_rows
        )
        try:
            solved, determinant = quotient_solve(pivot, rhs, b=b, d=d)
        except ValueError:
            continue
        unit_norm = quotient_unit_norm(determinant, b=b, d=d)
        if unit_norm == 0:
            continue
        vector: list[Residue] = [(Q(0), Q(0)) for _ in range(FIELD_COUNT)]
        for column, value in zip(pivot_columns, solved, strict=True):
            vector[column] = value
        vector[free_column] = (Q(1), Q(0))
        residual = tuple(
            _sum_residues(
                quotient_multiply(
                    quotient_symbol_matrix[row][column], vector[column], b=b, d=d
                )
                for column in range(FIELD_COUNT)
            )
            for row in range(FIELD_COUNT)
        )
        if any(value != (Q(0), Q(0)) for value in residual):
            continue
        return {
            "factor": monic_factor,
            "quadratic_discriminant": companion["discriminant"],
            "root_isolations": companion["roots"],
            "b": b,
            "d": d,
            "omitted_row": omitted_row,
            "free_column": free_column,
            "pivot_rows": pivot_rows,
            "pivot_columns": pivot_columns,
            "pivot_determinant": determinant,
            "pivot_unit_norm": unit_norm,
            "pivot_is_unit_at_both_roots": True,
            "second_order_vector": tuple(vector),
            "second_order_residue": residual,
            "second_order_residue_zero": True,
            "free_coordinate_normalized": vector[free_column] == (Q(1), Q(0)),
            "kernel_dimension_one_at_both_roots_by_unit_minor": True,
        }
    label = "required" if required_chart is not None else "deterministic"
    raise ValueError(f"no exact unit {label} pivot chart exists")


def _attach_first_order_lift(
    chart: dict[str, Any], first_order: dict[str, Any]
) -> None:
    b, d = chart["b"], chart["d"]
    c: Residue = (Q(0), Q(1))
    second = chart["second_order_vector"]
    first_vector = tuple(
        quotient_negate(quotient_multiply(c, value, b=b, d=d))
        for value in second
    ) + tuple(second)
    residual = tuple(
        _sum_residues(
            quotient_multiply(
                quotient_add(
                    (first_order["A_radial"][row][column], Q(0)),
                    quotient_negate(
                        quotient_multiply(
                            c,
                            (first_order["A_time"][row][column], Q(0)),
                            b=b,
                            d=d,
                        )
                    ),
                ),
                first_vector[column],
                b=b,
                d=d,
            )
            for column in range(2 * FIELD_COUNT)
        )
        for row in range(2 * FIELD_COUNT)
    )
    if any(value != (Q(0), Q(0)) for value in residual):
        raise ValueError("quadratic physical mode does not lift to the first-order pencil")
    chart["first_order_vector"] = first_vector
    chart["first_order_residue"] = residual
    chart["first_order_residue_zero"] = True


def physical_mode_certificate(
    state: SphericalState,
    *,
    tilde_normal_factor: Fraction | int = Q(4),
    hat_normal_factor: Fraction | int = Q(9),
) -> dict[str, Any]:
    """Build the exact physical metric/regulator quotient atlas at one point."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    unredefined = radial_symbol(state)
    quotient = quotient_symbol(unredefined, row_patch="rr")
    schur = metric_schur_data(quotient)
    _, physical_inverse = _physical_metric(state)
    physical_raw = inverse_metric_null_polynomial(physical_inverse)
    physical_factor = monic_quadratic_factor(physical_raw)
    regulator_raw = poly_div_exact(schur["scalar_determinant"], physical_raw)
    regulator_factor = monic_quadratic_factor(regulator_raw)
    scalar_determinant_monic: Polynomial = monic(schur["scalar_determinant"])
    factorization_exact = scalar_determinant_monic == poly_mul(
        physical_factor, regulator_factor
    )
    if not factorization_exact:
        raise ValueError("physical scalar determinant does not have the declared factors")

    full_symbol = modified_harmonic_symbol(
        state,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    charts = {
        "physical_metric": quadratic_mode_chart(
            full_symbol, physical_factor, required_chart=PHYSICAL_REQUIRED_CHART
        ),
        "regulator": quadratic_mode_chart(full_symbol, regulator_factor),
    }
    first_order = standard_first_order_principal_system(full_symbol)
    for chart in charts.values():
        _attach_first_order_lift(chart, first_order)
    if not all(
        chart["second_order_residue_zero"]
        and chart["first_order_residue_zero"]
        and chart["free_coordinate_normalized"]
        and chart["pivot_is_unit_at_both_roots"]
        and chart["kernel_dimension_one_at_both_roots_by_unit_minor"]
        for chart in charts.values()
    ):
        raise ValueError("physical quotient-mode certificate failed")
    return {
        "classification": "exact_pointwise_physical_quadratic_mode_atlas_not_a_uniform_domain_certificate",
        "quotient_patch": "rr",
        "rr_scalar_determinant": schur["scalar_determinant"],
        "rr_scalar_determinant_monic": scalar_determinant_monic,
        "physical_metric_factor": physical_factor,
        "regulator_factor": regulator_factor,
        "scalar_factorization_exact": factorization_exact,
        "charts": charts,
        "all_exact_checks_pass": True,
        "uniform_domain_eigenframe_proven": False,
        "strong_hyperbolicity_proven": False,
    }
