"""Exact MODE1 principal matrices and auxiliary-mode bases.

This module is intentionally pointwise.  It constructs the principal matrix
directly from the complete REF1 evaluator and compares it to MHG1; it does not
claim an open-domain eigenframe or strong hyperbolicity.
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from typing import Any, Sequence

from .exact_linear_algebra import Matrix, matrix_det, matrix_inverse, matrix_mul, matrix_rank, poly_eval
from .modified_harmonic import (
    _physical_metric,
    auxiliary_inverse_metric,
    inverse_metric_null_polynomial,
    modified_harmonic_symbol,
    quadratic_companion_certificate,
    standard_first_order_principal_system,
)
from .modified_harmonic_first_order import exact_full_residual_first_order_argument_jacobian
from .modified_harmonic_implicit import exact_full_residual_acceleration_jacobian
from .modified_harmonic_reference import ReferenceConnection, modified_harmonic_full_residuals
from .spherical_reduction import BASE_FIELD_ORDER, SphericalState
from .spherical_symbol import right_gauge_generators


Q = Fraction
FIELD_COUNT = len(BASE_FIELD_ORDER)
MODE1_DERIVATIVE_STATE_ORDER = tuple(f"dt_{name}" for name in BASE_FIELD_ORDER) + tuple(
    f"dr_{name}" for name in BASE_FIELD_ORDER
)
Residue = tuple[Fraction, Fraction]


def _fraction(value: object, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    return value


def _residue(value: Sequence[Fraction], *, name: str) -> Residue:
    if len(value) != 2:
        raise ValueError(f"{name} must be an (a, e) quadratic-quotient residue")
    return (_fraction(value[0], name=f"{name}[0]"), _fraction(value[1], name=f"{name}[1]"))


def _quotient_parameters(b: Fraction, d: Fraction) -> tuple[Fraction, Fraction]:
    return _fraction(b, name="b"), _fraction(d, name="d")


def quadratic_quotient_parameters(factor: Sequence[Fraction]) -> tuple[Fraction, Fraction]:
    """Validate monic ``c^2+b*c+d`` and return ``(b,d)``."""
    if len(factor) != 3:
        raise ValueError("quadratic quotient factor must have degree exactly two")
    if any(isinstance(value, bool) or not isinstance(value, Fraction) for value in factor):
        raise TypeError("quadratic quotient factor must contain Fractions only")
    if factor[2] != Q(1):
        raise ValueError("quadratic quotient factor must be monic c^2+b*c+d")
    return factor[1], factor[0]


def monic_quadratic_factor(factor: Sequence[Fraction]) -> tuple[Fraction, Fraction, Fraction]:
    """Normalize a nonzero exact quadratic without changing its root ideal."""
    if len(factor) != 3:
        raise ValueError("characteristic factor must have degree exactly two")
    if any(isinstance(value, bool) or not isinstance(value, Fraction) for value in factor):
        raise TypeError("characteristic factor must contain Fractions only")
    leading = factor[2]
    if leading == 0:
        raise ValueError("characteristic factor has no quadratic leading coefficient")
    normalized = tuple(value / leading for value in factor)
    quadratic_quotient_parameters(normalized)
    return normalized  # type: ignore[return-value]


def quotient_reduce(value: Sequence[Fraction], *, b: Fraction, d: Fraction) -> Residue:
    """Reduce a polynomial to ``a+e*c`` in Q[c]/(c^2+b*c+d)."""
    b, d = _quotient_parameters(b, d)
    if any(isinstance(item, bool) or not isinstance(item, Fraction) for item in value):
        raise TypeError("quotient polynomial must contain Fractions only")
    work = list(value) or [Q(0)]
    while len(work) > 2:
        leading = work.pop()
        degree = len(work) - 1
        work[degree - 1] -= leading * d
        work[degree] -= leading * b
    return (work[0] if work else Q(0), work[1] if len(work) > 1 else Q(0))


def quotient_add(left: Residue, right: Residue) -> Residue:
    left, right = _residue(left, name="left"), _residue(right, name="right")
    return (left[0] + right[0], left[1] + right[1])


def quotient_negate(value: Residue) -> Residue:
    value = _residue(value, name="value")
    return (-value[0], -value[1])


def quotient_multiply(left: Residue, right: Residue, *, b: Fraction, d: Fraction) -> Residue:
    left, right = _residue(left, name="left"), _residue(right, name="right")
    b, d = _quotient_parameters(b, d)
    a, e = left
    x, y = right
    return (a * x - e * y * d, a * y + e * x - e * y * b)


def quotient_unit_norm(value: Residue, *, b: Fraction, d: Fraction) -> Fraction:
    value = _residue(value, name="value")
    b, d = _quotient_parameters(b, d)
    a, e = value
    return a * a - a * e * b + e * e * d


def quotient_inverse(value: Residue, *, b: Fraction, d: Fraction) -> Residue:
    value = _residue(value, name="value")
    b, d = _quotient_parameters(b, d)
    a, e = value
    norm = quotient_unit_norm(value, b=b, d=d)
    if norm == 0:
        raise ValueError("quadratic quotient element is not a unit")
    return ((a - e * b) / norm, -e / norm)


def _qdet(matrix: Sequence[Sequence[Residue]], *, b: Fraction, d: Fraction) -> Residue:
    if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
        raise ValueError("quadratic quotient determinant currently requires 3x3")
    result: Residue = (Q(0), Q(0))
    for permutation, sign in (((0, 1, 2), 1), ((0, 2, 1), -1), ((1, 0, 2), -1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((2, 1, 0), -1)):
        term: Residue = (Q(1), Q(0))
        for row, column in enumerate(permutation):
            term = quotient_multiply(term, matrix[row][column], b=b, d=d)
        result = quotient_add(result, term if sign == 1 else quotient_negate(term))
    return result


def quotient_solve_3x3(matrix: Sequence[Sequence[Residue]], rhs: Sequence[Residue], *, b: Fraction, d: Fraction) -> tuple[Residue, Residue, Residue]:
    """Cramer's-rule solve over the validated quadratic quotient ring."""
    if len(rhs) != 3:
        raise ValueError("quadratic quotient right side must have three entries")
    determinant = _qdet(matrix, b=b, d=d)
    inverse_determinant = quotient_inverse(determinant, b=b, d=d)
    answer = []
    for column in range(3):
        replacement = [list(row) for row in matrix]
        for row in range(3):
            replacement[row][column] = rhs[row]
        answer.append(quotient_multiply(_qdet(replacement, b=b, d=d), inverse_determinant, b=b, d=d))
    return tuple(answer)  # type: ignore[return-value]


def _residue_vector_multiply(matrix: Sequence[Sequence[Residue]], vector: Sequence[Residue], *, b: Fraction, d: Fraction) -> tuple[Residue, ...]:
    return tuple(
        _sum_residues((quotient_multiply(entry, vector[column], b=b, d=d) for column, entry in enumerate(row)))
        for row in matrix
    )


def _sum_residues(values) -> Residue:
    total: Residue = (Q(0), Q(0))
    for value in values:
        total = quotient_add(total, value)
    return total


def _residue_symbol(symbol, *, b: Fraction, d: Fraction) -> tuple[tuple[Residue, ...], ...]:
    return tuple(tuple(quotient_reduce(entry, b=b, d=d) for entry in row) for row in symbol)


def _first_order_residual(first: dict[str, Any], vector: Sequence[Residue], *, b: Fraction, d: Fraction) -> tuple[Residue, ...]:
    # Pencil A_r-c A_t; multiplication by c is (0,1) in the quotient.
    c: Residue = (Q(0), Q(1))
    return tuple(
        _sum_residues(
            quotient_multiply(
                quotient_add((first["A_radial"][row][column], Q(0)), quotient_negate(quotient_multiply(c, (first["A_time"][row][column], Q(0)), b=b, d=d))),
                vector[column], b=b, d=d,
            )
            for column in range(2 * FIELD_COUNT)
        )
        for row in range(2 * FIELD_COUNT)
    )


def _identity(size: int) -> Matrix:
    return tuple(tuple(Q(int(row == column)) for column in range(size)) for row in range(size))


def _hstack(left: Matrix, right: Matrix) -> Matrix:
    return tuple(tuple(left[row] + right[row]) for row in range(len(left)))


def _rref_nullspace(matrix: Matrix) -> tuple[tuple[Fraction, ...], ...]:
    """Return an exact RREF null basis, with free variables in increasing order."""
    rows = [list(row) for row in matrix]
    columns = len(rows[0])
    pivots: list[int] = []
    pivot_row = 0
    for column in range(columns):
        candidate = next((row for row in range(pivot_row, len(rows)) if rows[row][column] != 0), None)
        if candidate is None:
            continue
        rows[pivot_row], rows[candidate] = rows[candidate], rows[pivot_row]
        divisor = rows[pivot_row][column]
        rows[pivot_row] = [value / divisor for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row != pivot_row and rows[row][column] != 0:
                factor = rows[row][column]
                rows[row] = [value - factor * pivot for value, pivot in zip(rows[row], rows[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    free = tuple(column for column in range(columns) if column not in pivots)
    basis = []
    for free_column in free:
        vector = [Q(0) for _ in range(columns)]
        vector[free_column] = Q(1)
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -rows[row][free_column]
        basis.append(tuple(vector))
    return tuple(basis)


def _symbol_at(symbol: Sequence[Sequence[Sequence[Fraction]]], speed: Fraction) -> Matrix:
    return tuple(tuple(poly_eval(entry, speed) for entry in row) for row in symbol)


def _pencil_at(first_order: dict[str, Any], speed: Fraction) -> Matrix:
    return tuple(
        tuple(first_order["A_radial"][row][column] - speed * first_order["A_time"][row][column] for column in range(2 * FIELD_COUNT))
        for row in range(2 * FIELD_COUNT)
    )


def ref1_branch_principal_matrix(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Build the exact REF1 branch principal matrix in MHG1 pencil convention.

    The top block is ``J_a^-1 [J_p_r J_q_r]``.  This is the coefficient in
    ``A_time d_t + A_radial d_r``; the evolution right-hand side has its usual
    opposite sign.
    """
    data = exact_full_residual_first_order_argument_jacobian(
        state, reference=reference, coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor, hat_normal_factor=hat_normal_factor,
    )
    blocks = data["blocks"]
    kinetic = blocks["p_t"]
    if matrix_det(kinetic) == 0:
        raise ValueError("REF1 coordinate-time kinetic block is singular")
    radial = _hstack(blocks["p_r"], blocks["q_r"])
    top = matrix_mul(matrix_inverse(kinetic), radial)
    zero = tuple(tuple(Q(0) for _ in range(FIELD_COUNT)) for _ in range(FIELD_COUNT))
    lower = _hstack(tuple(tuple(-entry for entry in row) for row in _identity(FIELD_COUNT)), zero)
    principal = tuple(top) + tuple(lower)
    return {
        "state_order": MODE1_DERIVATIVE_STATE_ORDER,
        "kinetic": kinetic,
        "kinetic_determinant": matrix_det(kinetic),
        "radial_blocks": {"p_r": blocks["p_r"], "q_r": blocks["q_r"]},
        "normal_principal_matrix": principal,
    }


def exact_comp1_acceleration_root(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
    qift_acceleration_half_width: Fraction = Q(1, 128),
) -> dict[str, Any]:
    """Derive and check COMP1's exact REF1 acceleration root at a fixed jet.

    This is a point computation, not a nonlinear acceleration map.  It uses
    the complete exact REF1 residual and its exact acceleration tangent block,
    then explicitly substitutes the Newton candidate back into all six rows.
    """
    qift_acceleration_half_width = _fraction(
        qift_acceleration_half_width, name="qift_acceleration_half_width"
    )
    if qift_acceleration_half_width <= 0:
        raise ValueError("qift_acceleration_half_width must be positive")
    if any(getattr(state, field).dtt != 0 for field in BASE_FIELD_ORDER):
        raise ValueError("COMP1 root derivation requires the declared zero-acceleration base state")
    tangent = exact_full_residual_acceleration_jacobian(
        state,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    acceleration_jacobian = tangent["jacobian"]
    if matrix_det(acceleration_jacobian) == 0:
        raise ValueError("COMP1 acceleration Jacobian is singular")
    residual_at_zero = tangent["base_full_residual_vector"]
    acceleration_inverse = matrix_inverse(acceleration_jacobian)
    root_acceleration = tuple(
        -sum(
            acceleration_inverse[row][column] * residual_at_zero[column]
            for column in range(FIELD_COUNT)
        )
        for row in range(FIELD_COUNT)
    )
    solved_state = replace(
        state,
        **{
            field: replace(getattr(state, field), dtt=root_acceleration[index])
            for index, field in enumerate(BASE_FIELD_ORDER)
        },
    )
    solved = modified_harmonic_full_residuals(
        solved_state,
        reference=reference,
        coordinate_radius=Q(coordinate_radius),
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    solved_residual = tuple(solved["full_residual_vector"])
    if any(value != 0 for value in solved_residual):
        raise ValueError("linearized COMP1 acceleration candidate is not an exact full REF1 root")
    strictly_inside = all(abs(value) < qift_acceleration_half_width for value in root_acceleration)
    if not strictly_inside:
        raise ValueError("COMP1 exact root lies outside the declared QIFT acceleration box")
    return {
        "classification": "exact_fixed_COMP1_REF1_acceleration_root_not_a_closed_form_nonlinear_acceleration_map",
        "acceleration_order": BASE_FIELD_ORDER,
        "residual_at_zero_acceleration": residual_at_zero,
        "acceleration_jacobian": acceleration_jacobian,
        "acceleration_jacobian_determinant": matrix_det(acceleration_jacobian),
        "root_acceleration": root_acceleration,
        "solved_full_residual_vector": solved_residual,
        "solved_full_residual_zero": True,
        "qift_acceleration_half_width": qift_acceleration_half_width,
        "root_strictly_inside_qift_acceleration_box": strictly_inside,
        "solved_state": solved_state,
    }


def auxiliary_mode_bases(
    state: SphericalState,
    *,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Return exact RREF bases for the two repeated auxiliary eigenspaces."""
    _, inverse = _physical_metric(state)
    symbol = modified_harmonic_symbol(state, tilde_normal_factor=tilde_normal_factor, hat_normal_factor=hat_normal_factor)
    first = standard_first_order_principal_system(symbol)
    output: dict[str, Any] = {}
    for name, factor in (
        ("tilde", inverse_metric_null_polynomial(auxiliary_inverse_metric(inverse, tilde_normal_factor))),
        ("hat", inverse_metric_null_polynomial(auxiliary_inverse_metric(inverse, hat_normal_factor))),
    ):
        roots = quadratic_companion_certificate(factor)["roots"]
        if any("exact" not in root for root in roots):
            raise ValueError(f"{name} auxiliary roots must be exact rationals")
        records = []
        for root in (roots[0]["exact"], roots[1]["exact"]):
            second = _symbol_at(symbol, root)
            pencil = _pencil_at(first, root)
            second_basis = _rref_nullspace(second)
            first_basis = _rref_nullspace(pencil)
            records.append({
                "speed": root,
                "second_order_basis": second_basis,
                "first_order_basis": first_basis,
                "second_order_kernel_dimension": len(second_basis),
                "first_order_kernel_dimension": len(first_basis),
                "second_order_residual_zero": all(all(sum(second[row][column] * vector[column] for column in range(FIELD_COUNT)) == 0 for row in range(FIELD_COUNT)) for vector in second_basis),
                "first_order_residual_zero": all(all(sum(pencil[row][column] * vector[column] for column in range(2 * FIELD_COUNT)) == 0 for row in range(2 * FIELD_COUNT)) for vector in first_basis),
            })
        output[name] = {"factor": factor, "roots": tuple(records), "pointwise_semisimple_at_declared_roots": all(item["second_order_kernel_dimension"] == 2 and item["first_order_kernel_dimension"] == 2 and item["second_order_residual_zero"] and item["first_order_residual_zero"] for item in records)}
    return output


def auxiliary_polynomial_mode_bases(
    state: SphericalState,
    *,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Build exact tilde/hat mode atlases over their quadratic quotient rings.

    The hat chart is deliberately nonflat: its selected quotient pivot must be
    a unit.  It is a MODE1 point chart, not a statement of uniformity.
    """
    _, inverse = _physical_metric(state)
    symbol = modified_harmonic_symbol(state, tilde_normal_factor=tilde_normal_factor, hat_normal_factor=hat_normal_factor)
    first = standard_first_order_principal_system(symbol)

    def first_vector(second_vector: Sequence[Residue], *, b: Fraction, d: Fraction) -> tuple[Residue, ...]:
        c: Residue = (Q(0), Q(1))
        return tuple(quotient_negate(quotient_multiply(c, value, b=b, d=d)) for value in second_vector) + tuple(second_vector)

    tilde_characteristic_factor = inverse_metric_null_polynomial(auxiliary_inverse_metric(inverse, tilde_normal_factor))
    tilde_factor = monic_quadratic_factor(tilde_characteristic_factor)
    tilde_b, tilde_d = quadratic_quotient_parameters(tilde_factor)
    tilde_symbol = _residue_symbol(symbol, b=tilde_b, d=tilde_d)
    gauge = right_gauge_generators()
    tilde_vectors = tuple(tuple(quotient_reduce(gauge[row][column], b=tilde_b, d=tilde_d) for row in range(FIELD_COUNT)) for column in range(2))
    tilde_second = tuple(_residue_vector_multiply(tilde_symbol, vector, b=tilde_b, d=tilde_d) for vector in tilde_vectors)
    tilde_first_vectors = tuple(first_vector(vector, b=tilde_b, d=tilde_d) for vector in tilde_vectors)
    tilde_first = tuple(_first_order_residual(first, vector, b=tilde_b, d=tilde_d) for vector in tilde_first_vectors)
    if any(any(value != (Q(0), Q(0)) for value in residual) for residual in tilde_second):
        raise ValueError("SYM1 right gauge generators do not annihilate the full MHG symbol modulo tilde q")
    if any(any(value != (Q(0), Q(0)) for value in residual) for residual in tilde_first):
        raise ValueError("tilde polynomial first-order modes do not annihilate the MHG pencil modulo tilde q")

    hat_characteristic_factor = inverse_metric_null_polynomial(auxiliary_inverse_metric(inverse, hat_normal_factor))
    hat_factor = monic_quadratic_factor(hat_characteristic_factor)
    hat_b, hat_d = quadratic_quotient_parameters(hat_factor)
    hat_symbol = _residue_symbol(symbol, b=hat_b, d=hat_d)
    pivot_rows, pivot_columns, free_columns = (1, 3, 4), (0, 1, 2), (3, 4)
    pivot = tuple(tuple(hat_symbol[row][column] for column in pivot_columns) for row in pivot_rows)
    pivot_determinant = _qdet(pivot, b=hat_b, d=hat_d)
    # The inverse call is also the exact unit test for this nonflat chart.
    pivot_inverse = quotient_inverse(pivot_determinant, b=hat_b, d=hat_d)
    hat_vectors = []
    for free_column in free_columns:
        rhs = tuple(quotient_negate(hat_symbol[row][free_column]) for row in pivot_rows)
        solved = quotient_solve_3x3(pivot, rhs, b=hat_b, d=hat_d)
        vector = [(Q(0), Q(0)) for _ in range(FIELD_COUNT)]
        for column, value in zip(pivot_columns, solved):
            vector[column] = value
        vector[free_column] = (Q(1), Q(0))
        hat_vectors.append(tuple(vector))
    hat_vectors = tuple(hat_vectors)
    hat_second = tuple(_residue_vector_multiply(hat_symbol, vector, b=hat_b, d=hat_d) for vector in hat_vectors)
    if any(any(value != (Q(0), Q(0)) for value in residual) for residual in hat_second):
        raise ValueError("hat polynomial pivot chart does not annihilate the full MHG symbol modulo hat q")
    hat_first_vectors = tuple(first_vector(vector, b=hat_b, d=hat_d) for vector in hat_vectors)
    hat_first = tuple(_first_order_residual(first, vector, b=hat_b, d=hat_d) for vector in hat_first_vectors)
    if any(any(value != (Q(0), Q(0)) for value in residual) for residual in hat_first):
        raise ValueError("hat polynomial first-order modes do not annihilate the MHG pencil modulo hat q")
    free_coordinate_identity = tuple(tuple(vector[column] for vector in hat_vectors) for column in free_columns) == (((Q(1), Q(0)), (Q(0), Q(0))), ((Q(0), Q(0)), (Q(1), Q(0))))
    if not free_coordinate_identity:
        raise ValueError("hat free-coordinate identity failed; the two quotient modes are not independent")
    return {
        "classification": "exact_pointwise_polynomial_auxiliary_mode_atlas_not_a_uniform_domain_certificate",
        "tilde": {
            "characteristic_factor": tilde_characteristic_factor,
            "factor": tilde_factor, "b": tilde_b, "d": tilde_d,
            "basis_source": "SYM1_right_gauge_generators",
            "second_order_vectors": tilde_vectors,
            "second_order_residues": tilde_second,
            "first_order_vectors": tilde_first_vectors,
            "first_order_residues": tilde_first,
            "all_residues_zero": all(all(value == (Q(0), Q(0)) for value in residual) for residual in tilde_second + tilde_first),
        },
        "hat": {
            "characteristic_factor": hat_characteristic_factor,
            "factor": hat_factor, "b": hat_b, "d": hat_d,
            "pivot_rows": pivot_rows, "pivot_columns": pivot_columns, "free_columns": free_columns,
            "pivot_determinant": pivot_determinant,
            "pivot_unit_norm": quotient_unit_norm(pivot_determinant, b=hat_b, d=hat_d),
            "pivot_inverse": pivot_inverse,
            "second_order_vectors": hat_vectors,
            "second_order_residues": hat_second,
            "first_order_vectors": hat_first_vectors,
            "first_order_residues": hat_first,
            "free_coordinate_identity_exact": free_coordinate_identity,
            "rank_two_at_each_hat_root_by_free_coordinate_identity": free_coordinate_identity,
            "all_residues_zero": True,
        },
    }


def mode1_point_certificate(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Exact point certificate; it deliberately makes no open-domain claim."""
    branch = ref1_branch_principal_matrix(state, reference=reference, coordinate_radius=coordinate_radius, tilde_normal_factor=tilde_normal_factor, hat_normal_factor=hat_normal_factor)
    mhg = standard_first_order_principal_system(modified_harmonic_symbol(state, tilde_normal_factor=tilde_normal_factor, hat_normal_factor=hat_normal_factor))
    modes = auxiliary_mode_bases(state, tilde_normal_factor=tilde_normal_factor, hat_normal_factor=hat_normal_factor)
    return {"branch": branch, "mhg1": mhg, "ref1_branch_equals_mhg1": branch["normal_principal_matrix"] == mhg["normal_principal_matrix"], "auxiliary_modes": modes, "all_auxiliary_modes_pointwise_semisimple": all(modes[name]["pointwise_semisimple_at_declared_roots"] for name in ("tilde", "hat")), "classification": "exact_pointwise_MODE1_principal_and_auxiliary_mode_certificate_not_a_uniform_domain_or_evolution_theorem"}
