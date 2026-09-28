"""SYM1 exact radial principal-symbol and ADM formulation preflight.

The module starts from RED1's unredefined six-by-eighteen second-jet matrix.
It constructs the radial covariant polynomial symbol, verifies its exact
diffeomorphism degeneracies, takes one declared quotient chart, and performs
an exact metric Schur elimination.  Separately it maps coordinate second jets
to normal/radial first-order variables in the generalized ADM chart.

Everything here is pointwise and exact.  It is not a nonlinear evolution
system, a constraint-propagation proof, or an open-domain hyperbolicity proof.
"""

from __future__ import annotations

from fractions import Fraction
from math import isqrt
from typing import Any, Mapping, Sequence

from .exact_linear_algebra import (
    ONE,
    ZERO,
    Matrix,
    Polynomial,
    matrix_det,
    matrix_inverse,
    matrix_mul,
    matrix_rank,
    poly,
    poly_add,
    poly_div_exact,
    poly_eval,
    poly_mul,
    poly_scale,
    poly_sub,
    polynomial_matrix_det,
)
from .spherical_reduction import (
    ADM_FIELD_ORDER,
    BASE_FIELD_ORDER,
    INDEPENDENT_EQUATION_ORDER,
    SECOND_DERIVATIVE_ORDER,
    SphericalState,
    adm_pg_principal_matrix,
    direct_4d_curvature,
    principal_matrix,
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction

RADIAL_COVECTOR = "xi_A=(-c,1)"
QUOTIENT_VARIABLE_ORDER = (
    "h_tt_representative",
    "areal_radius",
    "phi",
    "chi",
)
QUOTIENT_RR_EQUATION_ORDER = (
    "metric_rr",
    "metric_theta_theta",
    "scalar_phi",
    "scalar_chi",
)
QUOTIENT_TT_EQUATION_ORDER = (
    "metric_tt",
    "metric_theta_theta",
    "scalar_phi",
    "scalar_chi",
)
NORMAL_VARIABLE_ORDER = (
    "A_perp",
    "A_r",
    "L_r",
    "B_perp",
    "B_r",
    "S_r",
    "P_phi_perp",
    "P_phi_r",
    "Q_phi_r",
    "P_chi_perp",
    "P_chi_r",
    "Q_chi_r",
)
ADM_PROJECTION_ORDER = (
    "E_nn",
    "E_ns",
    "E_ss",
    "E_theta_theta_orthonormal",
    "E_phi",
    "E_chi",
)
KINETIC_VARIABLE_ORDER = (
    "A_perp",
    "B_perp",
    "P_phi_perp",
    "P_chi_perp",
)
GAUGE_SOURCE_SECOND_JET_ORDER = (
    "alpha.dtt",
    "alpha.dtr",
    "alpha.drr",
    "shift.dtt",
    "shift.dtr",
    "shift.drr",
)
FIRST_ORDER_STATE_ORDER = (
    "A",
    "B",
    "P_phi",
    "P_chi",
    "L",
    "S",
    "Q_phi",
    "Q_chi",
)


def _fraction(value: Any, name: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an exact rational value")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str) and value:
        try:
            result = Fraction(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{name} must be an exact rational value") from exc
        canonical = (
            str(result.numerator)
            if result.denominator == 1
            else f"{result.numerator}/{result.denominator}"
        )
        if value != canonical:
            raise ValueError(f"{name} must use canonical rational encoding")
        return result
    raise ValueError(f"{name} must be an exact rational value")


def _fixture_value(fixture: Mapping[str, Any], name: str) -> Fraction:
    state = fixture.get("state", fixture)
    if not isinstance(state, Mapping) or name not in state:
        raise ValueError(f"fixture state is missing {name}")
    value = state[name]
    if isinstance(value, Mapping):
        if "value" not in value:
            raise ValueError(f"fixture state {name} jet is missing value")
        value = value["value"]
    return _fraction(value, f"fixture state {name}")


def _zero_polynomial_matrix(rows: int, columns: int) -> list[list[Polynomial]]:
    return [[ZERO for _ in range(columns)] for _ in range(rows)]


def _polynomial_matrix_multiply(
    left: Sequence[Sequence[Sequence[Fraction]]],
    right: Sequence[Sequence[Sequence[Fraction]]],
) -> tuple[tuple[Polynomial, ...], ...]:
    if not left or not right or not left[0] or not right[0]:
        raise ValueError("polynomial matrices must be nonempty")
    if any(len(row) != len(left[0]) for row in left):
        raise ValueError("left polynomial matrix must be rectangular")
    if any(len(row) != len(right[0]) for row in right):
        raise ValueError("right polynomial matrix must be rectangular")
    if len(left[0]) != len(right):
        raise ValueError("polynomial matrix dimensions do not agree")
    return tuple(
        tuple(
            _poly_sum(
                poly_mul(left[row][inner], right[inner][column])
                for inner in range(len(right))
            )
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def _poly_sum(values) -> Polynomial:
    result = ZERO
    for value in values:
        result = poly_add(result, value)
    return result


def _polynomial_matrix_subtract(left, right):
    if len(left) != len(right) or any(
        len(a) != len(b) for a, b in zip(left, right)
    ):
        raise ValueError("polynomial matrix shapes do not agree")
    return tuple(
        tuple(poly_sub(a, b) for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def _polynomial_matrix_scale(value, factor: Sequence[Fraction]):
    return tuple(
        tuple(poly_mul(entry, factor) for entry in row) for row in value
    )


def _polynomial_matrix_zero(value) -> bool:
    return all(entry == ZERO for row in value for entry in row)


def radial_symbol_from_principal_matrix(
    data: Mapping[str, Any],
) -> tuple[tuple[Polynomial, ...], ...]:
    """Construct ``P(-c,1)`` from one RED1 base-field principal matrix."""

    if tuple(data.get("field_order", ())) != BASE_FIELD_ORDER:
        raise ValueError("principal field order differs from RED1")
    if tuple(data.get("second_derivative_order", ())) != SECOND_DERIVATIVE_ORDER:
        raise ValueError("principal derivative order differs from RED1")
    if tuple(data.get("equation_order", ())) != INDEPENDENT_EQUATION_ORDER:
        raise ValueError("principal equation order differs from RED1")
    columns = tuple(data.get("column_order", ()))
    expected_columns = tuple(
        f"{field}.{derivative}"
        for field in BASE_FIELD_ORDER
        for derivative in SECOND_DERIVATIVE_ORDER
    )
    if columns != expected_columns:
        raise ValueError("principal column order differs from RED1")
    matrix = data.get("matrix")
    if not isinstance(matrix, (tuple, list)) or len(matrix) != 6:
        raise ValueError("principal matrix must have six rows")
    if any(not isinstance(row, (tuple, list)) or len(row) != 18 for row in matrix):
        raise ValueError("principal matrix must be six by eighteen")

    output = _zero_polynomial_matrix(6, 6)
    for equation in range(6):
        for field_index, field in enumerate(BASE_FIELD_ORDER):
            dtt = matrix[equation][columns.index(f"{field}.dtt")]
            dtr = matrix[equation][columns.index(f"{field}.dtr")]
            drr = matrix[equation][columns.index(f"{field}.drr")]
            output[equation][field_index] = poly((drr, -dtr, dtt))
    return tuple(tuple(row) for row in output)


def radial_symbol(state: SphericalState) -> tuple[tuple[Polynomial, ...], ...]:
    return radial_symbol_from_principal_matrix(principal_matrix(state))


def right_gauge_generators() -> tuple[tuple[Polynomial, ...], ...]:
    """Return field-by-generator pure radial diffeomorphism polarizations."""

    return (
        (poly((0, -2)), ZERO),
        (ONE, poly((0, -1))),
        (ZERO, poly((2,))),
        (ZERO, ZERO),
        (ZERO, ZERO),
        (ZERO, ZERO),
    )


def left_bianchi_generators(
    state: SphericalState,
) -> tuple[tuple[Polynomial, ...], ...]:
    """Return the two principal rows implementing ``xi^A E_AB=0``."""

    _riemann, inverse, _connection = direct_4d_curvature(state)
    xi_t_up = poly((inverse[0][1], -inverse[0][0]))
    xi_r_up = poly((inverse[1][1], -inverse[0][1]))
    return (
        (xi_t_up, xi_r_up, ZERO, ZERO, ZERO, ZERO),
        (ZERO, xi_t_up, xi_r_up, ZERO, ZERO, ZERO),
    )


def quotient_symbol(
    symbol, *, row_patch: str = "rr"
) -> tuple[tuple[Polynomial, ...], ...]:
    """Return one declared four-by-four Bianchi-row quotient chart.

    ``rr`` is valid where ``xi^t != 0`` (equivalently ``c+v != 0`` in the
    ADM chart). ``tt`` supplies the overlapping chart needed at ``c=-v``.
    Both use the same globally nonsingular right-gauge representative columns.
    """

    if len(symbol) != 6 or any(len(row) != 6 for row in symbol):
        raise ValueError("radial symbol must be six by six")
    if row_patch == "rr":
        equation_indices = (2, 3, 4, 5)
    elif row_patch == "tt":
        equation_indices = (0, 3, 4, 5)
    else:
        raise ValueError("row_patch must be rr or tt")
    variable_indices = (0, 3, 4, 5)
    return tuple(
        tuple(poly(symbol[row][column]) for column in variable_indices)
        for row in equation_indices
    )


def metric_schur_data(quotient) -> dict[str, Any]:
    """Return exact 2+2 block determinants and Schur numerator data."""

    if len(quotient) != 4 or any(len(row) != 4 for row in quotient):
        raise ValueError("quotient symbol must be four by four")
    pgg = tuple(tuple(quotient[i][j] for j in range(2)) for i in range(2))
    pgs = tuple(tuple(quotient[i][j] for j in range(2, 4)) for i in range(2))
    psg = tuple(tuple(quotient[i][j] for j in range(2)) for i in range(2, 4))
    pss = tuple(tuple(quotient[i][j] for j in range(2, 4)) for i in range(2, 4))
    det_gg = polynomial_matrix_det(pgg)
    adj_gg = (
        (pgg[1][1], poly_scale(pgg[0][1], -1)),
        (poly_scale(pgg[1][0], -1), pgg[0][0]),
    )
    correction = _polynomial_matrix_multiply(
        _polynomial_matrix_multiply(psg, adj_gg), pgs
    )
    schur_numerator = _polynomial_matrix_subtract(
        _polynomial_matrix_scale(pss, det_gg), correction
    )
    det_quotient = polynomial_matrix_det(quotient)
    scalar_determinant = poly_div_exact(det_quotient, det_gg)
    schur_identity = polynomial_matrix_det(schur_numerator) == poly_mul(
        det_gg, det_quotient
    )
    return {
        "metric_block": pgg,
        "metric_to_scalar_block": pgs,
        "scalar_to_metric_block": psg,
        "bare_scalar_block": pss,
        "metric_block_determinant": det_gg,
        "quotient_determinant": det_quotient,
        "schur_numerator": schur_numerator,
        "schur_denominator": det_gg,
        "scalar_determinant": scalar_determinant,
        "schur_determinant_identity": schur_identity,
        "scalar_offdiagonal_numerators_zero": (
            schur_numerator[0][1] == ZERO and schur_numerator[1][0] == ZERO
        ),
    }


def quadratic_discriminant(value: Sequence[Fraction]) -> Fraction:
    coefficients = poly(value)
    if len(coefficients) != 3 or coefficients[2] == 0:
        raise ValueError("characteristic factor must be a genuine quadratic")
    d, b, a = coefficients
    return b * b - 4 * a * d


def _fraction_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator == value.numerator and denominator * denominator == value.denominator:
        return Fraction(numerator, denominator)
    return None


def _same_sign(left: Fraction, right: Fraction) -> bool:
    return (left > 0 and right > 0) or (left < 0 and right < 0)


def _isolate_quadratic_roots(
    value: Sequence[Fraction], *, bisections: int = 24
) -> tuple[dict[str, Fraction], dict[str, Fraction]]:
    coefficients = poly(value)
    discriminant = quadratic_discriminant(coefficients)
    if discriminant <= 0:
        raise ValueError("quadratic roots are not real and distinct")
    d, b, a = coefficients
    root_discriminant = _fraction_square_root(discriminant)
    if root_discriminant is not None:
        roots = sorted(((-b - root_discriminant) / (2 * a), (-b + root_discriminant) / (2 * a)))
        return ({"exact": roots[0]}, {"exact": roots[1]})

    bound_seed = max(abs(d / a), abs(b / a), Fraction(1)) + 1
    bound = (bound_seed.numerator + bound_seed.denominator - 1) // bound_seed.denominator
    while (
        poly_eval(coefficients, -bound) == 0
        or poly_eval(coefficients, bound) == 0
        or not _same_sign(poly_eval(coefficients, -bound), a)
        or not _same_sign(poly_eval(coefficients, bound), a)
    ):
        bound += 1

    # A dyadic scan keeps the public isolating intervals compact.  Because a
    # non-square positive discriminant gives two irrational simple roots, no
    # grid point can equal either root and each is eventually isolated by a
    # separate sign-changing cell.
    for precision in range(4, bisections + 1):
        denominator = 1 << precision
        first = -bound * denominator
        last = bound * denominator
        intervals: list[dict[str, Fraction]] = []
        lower = Fraction(first, denominator)
        lower_value = poly_eval(coefficients, lower)
        for numerator in range(first + 1, last + 1):
            upper = Fraction(numerator, denominator)
            upper_value = poly_eval(coefficients, upper)
            if not _same_sign(lower_value, upper_value):
                intervals.append({"lower": lower, "upper": upper})
            lower, lower_value = upper, upper_value
        if len(intervals) == 2:
            refined: list[dict[str, Fraction]] = []
            for interval in intervals:
                lower = interval["lower"]
                upper = interval["upper"]
                lower_value = poly_eval(coefficients, lower)
                for _ in range(bisections - precision):
                    midpoint = (lower + upper) / 2
                    midpoint_value = poly_eval(coefficients, midpoint)
                    if _same_sign(lower_value, midpoint_value):
                        lower, lower_value = midpoint, midpoint_value
                    else:
                        upper = midpoint
                refined.append({"lower": lower, "upper": upper})
            return refined[0], refined[1]
    raise ValueError("failed to isolate both quadratic roots within dyadic budget")


def quadratic_companion_certificate(value: Sequence[Fraction]) -> dict[str, Any]:
    coefficients = poly(value)
    discriminant = quadratic_discriminant(coefficients)
    if discriminant <= 0:
        raise ValueError("quadratic must have positive discriminant")
    d, b, a = coefficients
    companion = ((Q(0), -d / a), (Q(1), -b / a))
    symmetrizer = (
        (Q(1), -b / (2 * a)),
        (-b / (2 * a), -d / a + b * b / (2 * a * a)),
    )
    product = matrix_mul(symmetrizer, companion)
    symmetrizer_determinant = matrix_det(symmetrizer)
    return {
        "polynomial": coefficients,
        "discriminant": discriminant,
        "roots": _isolate_quadratic_roots(coefficients),
        "companion": companion,
        "symmetrizer": symmetrizer,
        "symmetrizer_times_companion": product,
        "product_is_symmetric": product[0][1] == product[1][0],
        "symmetrizer_first_minor_positive": symmetrizer[0][0] > 0,
        "symmetrizer_determinant": symmetrizer_determinant,
        "symmetrizer_positive_definite": (
            symmetrizer[0][0] > 0 and symmetrizer_determinant > 0
        ),
        "expected_symmetrizer_determinant": discriminant / (4 * a * a),
    }


def _quadratic_resultant(left: Sequence[Fraction], right: Sequence[Fraction]) -> Fraction:
    a0, a1, a2 = poly(left)
    b0, b1, b2 = poly(right)
    return matrix_det(
        (
            (a2, a1, a0, Q(0)),
            (Q(0), a2, a1, a0),
            (b2, b1, b0, Q(0)),
            (Q(0), b2, b1, b0),
        )
    )


def _polynomials_proportional(left: Sequence[Fraction], right: Sequence[Fraction]) -> bool:
    a, b = poly(left), poly(right)
    if len(a) != len(b):
        return False
    pivot = next((index for index, value in enumerate(b) if value), None)
    if pivot is None:
        return a == ZERO
    ratio = a[pivot] / b[pivot]
    return all(x == ratio * y for x, y in zip(a, b))


def adm_normal_principal_matrix(fixture: Mapping[str, Any]) -> dict[str, Any]:
    """Map the RED1 ADM second jets to normal/radial first-order derivatives."""

    raw = adm_pg_principal_matrix(fixture)
    matrix = raw["matrix"]
    columns = tuple(raw["column_order"])
    alpha = _fixture_value(fixture, "alpha")
    shift = _fixture_value(fixture, "shift")
    radial_metric = _fixture_value(fixture, "lambda")
    areal_radius = _fixture_value(fixture, "areal_radius")
    if alpha <= 0 or radial_metric <= 0 or areal_radius <= 0:
        raise ValueError("ADM alpha, lambda, and areal radius must be positive")

    projection = [[Q(0) for _ in range(6)] for _ in range(6)]
    projection[0] = [
        1 / (alpha * alpha),
        -2 * shift / (alpha * alpha),
        shift * shift / (alpha * alpha),
        Q(0),
        Q(0),
        Q(0),
    ]
    projection[1] = [
        Q(0),
        1 / (alpha * radial_metric),
        -shift / (alpha * radial_metric),
        Q(0),
        Q(0),
        Q(0),
    ]
    projection[2][2] = 1 / (radial_metric * radial_metric)
    projection[3][3] = 1 / (areal_radius * areal_radius)
    projection[4][4] = Q(1)
    projection[5][5] = Q(1)

    transform = [[Q(0) for _ in NORMAL_VARIABLE_ORDER] for _ in columns]
    gauge_transform = [
        [Q(0) for _ in GAUGE_SOURCE_SECOND_JET_ORDER] for _ in columns
    ]

    def put(field: str, derivative: str, variable: str, coefficient: Fraction) -> None:
        transform[columns.index(f"{field}.{derivative}")][
            NORMAL_VARIABLE_ORDER.index(variable)
        ] = coefficient

    put("lambda", "dtt", "A_perp", -alpha * radial_metric)
    put("lambda", "dtt", "A_r", -2 * alpha * radial_metric * shift)
    put("lambda", "dtt", "L_r", shift * shift * radial_metric)
    put("lambda", "dtr", "A_r", -alpha * radial_metric)
    put("lambda", "dtr", "L_r", shift * radial_metric)
    put("lambda", "drr", "L_r", radial_metric)

    put("areal_radius", "dtt", "B_perp", -alpha * areal_radius)
    put("areal_radius", "dtt", "B_r", -2 * alpha * shift * areal_radius)
    put(
        "areal_radius",
        "dtt",
        "S_r",
        shift * shift * radial_metric * areal_radius,
    )
    put("areal_radius", "dtr", "B_r", -alpha * areal_radius)
    put(
        "areal_radius",
        "dtr",
        "S_r",
        shift * radial_metric * areal_radius,
    )
    put("areal_radius", "drr", "S_r", radial_metric * areal_radius)

    for field, momentum, gradient in (
        ("phi", "P_phi", "Q_phi"),
        ("chi", "P_chi", "Q_chi"),
    ):
        put(field, "dtt", momentum + "_perp", alpha)
        put(field, "dtt", momentum + "_r", 2 * alpha * shift)
        put(field, "dtt", gradient + "_r", shift * shift * radial_metric)
        put(field, "dtr", momentum + "_r", alpha)
        put(field, "dtr", gradient + "_r", shift * radial_metric)
        put(field, "drr", gradient + "_r", radial_metric)

    # Complete the principal Jacobian for prescribed lapse/shift sources while
    # holding A fixed.  Since
    # A=-(D Lambda/Lambda-v_r)/alpha, varying v_tr or v_rr also changes the
    # coordinate Lambda second jets at fixed A_perp and A_r.
    for source_index, source in enumerate(GAUGE_SOURCE_SECOND_JET_ORDER):
        gauge_transform[columns.index(source)][source_index] = Q(1)
    shift_tr_index = GAUGE_SOURCE_SECOND_JET_ORDER.index("shift.dtr")
    shift_rr_index = GAUGE_SOURCE_SECOND_JET_ORDER.index("shift.drr")
    gauge_transform[columns.index("lambda.dtt")][shift_tr_index] = radial_metric
    gauge_transform[columns.index("lambda.dtt")][shift_rr_index] = (
        shift * radial_metric
    )
    gauge_transform[columns.index("lambda.dtr")][shift_rr_index] = radial_metric

    projected = matrix_mul(projection, matrix)
    normal = matrix_mul(projected, transform)
    gauge_sources = matrix_mul(projected, gauge_transform)
    acceleration_indices = tuple(
        NORMAL_VARIABLE_ORDER.index(name) for name in KINETIC_VARIABLE_ORDER
    )
    constraints_free = all(
        normal[row][column] == 0
        for row in (0, 1)
        for column in acceleration_indices
    )
    constraints_free_of_gauge_second_jets = all(
        gauge_sources[row][column] == 0
        for row in (0, 1)
        for column in range(len(GAUGE_SOURCE_SECOND_JET_ORDER))
    )
    kinetic = tuple(
        tuple(normal[row][column] for column in acceleration_indices)
        for row in (2, 3, 4, 5)
    )
    return {
        "projection_order": ADM_PROJECTION_ORDER,
        "normal_variable_order": NORMAL_VARIABLE_ORDER,
        "matrix": normal,
        "candidate_constraint_rows": normal[:2],
        "candidate_constraints_free_of_normal_accelerations": constraints_free,
        "gauge_source_second_jet_order": GAUGE_SOURCE_SECOND_JET_ORDER,
        "gauge_source_matrix_at_fixed_normal_variables": gauge_sources,
        "candidate_constraints_free_of_gauge_source_second_jets": (
            constraints_free_of_gauge_second_jets
        ),
        "kinetic_variable_order": KINETIC_VARIABLE_ORDER,
        "kinetic_matrix": kinetic,
        "kinetic_determinant": matrix_det(kinetic),
    }


def naive_fixed_gauge_comparator(fixture: Mapping[str, Any]) -> dict[str, Any]:
    """Return the exact defect test for the tempting unconstrained ADM system."""

    normal_data = adm_normal_principal_matrix(fixture)
    matrix = normal_data["matrix"]
    alpha = _fixture_value(fixture, "alpha")
    shift = _fixture_value(fixture, "shift")
    radial_metric = _fixture_value(fixture, "lambda")
    evolution = matrix[2:]
    velocity_perp = (0, 3, 6, 9)
    velocity_r = (1, 4, 7, 10)
    gradient_r = (2, 5, 8, 11)
    kinetic = tuple(
        tuple(row[column] for column in velocity_perp) for row in evolution
    )
    velocity = tuple(
        tuple(row[column] for column in velocity_r) for row in evolution
    )
    gradient = tuple(
        tuple(row[column] for column in gradient_r) for row in evolution
    )
    if matrix_det(kinetic) == 0:
        raise ValueError("naive comparator kinetic matrix is singular")

    at = [[Q(0) for _ in range(8)] for _ in range(8)]
    ar = [[Q(0) for _ in range(8)] for _ in range(8)]
    for row in range(4):
        for column in range(4):
            at[row][column] = kinetic[row][column]
            # The evolution rows are assembled in D=partial_t-v partial_r
            # variables.  Convert K D X+V X_r+G Y_r=0 to coordinate time:
            # K X_t+(V-v K)X_r+G Y_r=0.
            ar[row][column] = (
                velocity[row][column] - shift * kinetic[row][column]
            )
            ar[row][column + 4] = gradient[row][column]
        at[row + 4][row + 4] = Q(1)
        # The kinematic compatibility equations are likewise D Y plus the
        # displayed X_r coupling, so their own radial advection is -v Y_r.
        ar[row + 4][row + 4] = -shift

    compatibility = (alpha, alpha / radial_metric, -alpha / radial_metric, -alpha / radial_metric)
    for index, coefficient in enumerate(compatibility):
        ar[index + 4][index] = coefficient
    principal = matrix_mul(matrix_inverse(at), ar)

    def power(value: Matrix, exponent: int) -> Matrix:
        identity = tuple(
            tuple(Q(int(row == column)) for column in range(len(value)))
            for row in range(len(value))
        )
        result = identity
        for _ in range(exponent):
            result = matrix_mul(result, value)
        return result

    geometric = 8 - matrix_rank(principal)
    generalized = 8 - matrix_rank(power(principal, 8))
    return {
        "state_order": FIRST_ORDER_STATE_ORDER,
        "normal_principal_matrix": principal,
        "zero_speed_geometric_multiplicity": geometric,
        "zero_speed_generalized_multiplicity": generalized,
        "defective_zero_speed_sector": generalized > geometric,
        "strong_hyperbolicity_rejected": generalized > geometric,
    }


def symbol_fixture_certificate(fixture: Mapping[str, Any]) -> dict[str, Any]:
    """Build one exact pointwise SYM1 certificate from a RED1 fixture."""

    state = state_from_generalized_adm_pg_fixture(fixture)
    symbol = radial_symbol(state)
    gauge_residual = _polynomial_matrix_multiply(symbol, right_gauge_generators())
    bianchi_residual = _polynomial_matrix_multiply(
        left_bianchi_generators(state), symbol
    )
    quotient_rr = quotient_symbol(symbol, row_patch="rr")
    quotient_tt = quotient_symbol(symbol, row_patch="tt")
    schur_rr = metric_schur_data(quotient_rr)
    schur_tt = metric_schur_data(quotient_tt)
    shift = _fixture_value(fixture, "shift")
    alpha = _fixture_value(fixture, "alpha")
    radial_metric = _fixture_value(fixture, "lambda")
    patch_factor = poly((shift * shift, 2 * shift, 1))
    patch_quotient = poly_div_exact(
        schur_rr["metric_block_determinant"], patch_factor
    )
    if len(patch_quotient) != 1 or patch_quotient[0] == 0:
        raise ValueError("metric quotient block lacks the nonzero gauge-patch factor")

    metric_null = poly(
        (shift * shift - alpha * alpha / (radial_metric * radial_metric), 2 * shift, 1)
    )
    if schur_rr["scalar_determinant"] != schur_tt["scalar_determinant"]:
        raise ValueError("quotient atlas patches disagree on the scalar determinant")
    regulator = poly_div_exact(schur_rr["scalar_determinant"], metric_null)
    metric_null_certificate = quadratic_companion_certificate(metric_null)
    regulator_certificate = quadratic_companion_certificate(regulator)
    cones_proportional = _polynomials_proportional(regulator, metric_null)
    cone_resultant = _quadratic_resultant(metric_null, regulator)
    metric_patch_point = -shift
    metric_null_at_patch = poly_eval(metric_null, metric_patch_point)
    regulator_at_patch = poly_eval(regulator, metric_patch_point)
    alternate_metric_block_at_patch = poly_eval(
        schur_tt["metric_block_determinant"], metric_patch_point
    )
    adm = adm_normal_principal_matrix(fixture)
    return {
        "fixture_id": str(fixture.get("fixture_id", "unspecified")),
        "branch": state.branch,
        "radial_covector": RADIAL_COVECTOR,
        "raw_symbol": symbol,
        "right_gauge_generators": right_gauge_generators(),
        "right_gauge_residual": gauge_residual,
        "right_gauge_null_identity_exact": _polynomial_matrix_zero(gauge_residual),
        "left_bianchi_generators": left_bianchi_generators(state),
        "left_bianchi_residual": bianchi_residual,
        "left_bianchi_null_identity_exact": _polynomial_matrix_zero(bianchi_residual),
        "quotient_variable_order": QUOTIENT_VARIABLE_ORDER,
        "quotient_atlas": {
            "rr_patch": {
                "equation_order": QUOTIENT_RR_EQUATION_ORDER,
                "validity": "c+shift_nonzero",
                "quotient_symbol": quotient_rr,
                "metric_schur": schur_rr,
            },
            "tt_patch": {
                "equation_order": QUOTIENT_TT_EQUATION_ORDER,
                "validity": "covers_c_equals_minus_shift",
                "quotient_symbol": quotient_tt,
                "metric_schur": schur_tt,
            },
        },
        "quotient_atlas_scalar_determinants_agree_exactly": schur_rr[
            "scalar_determinant"
        ]
        == schur_tt["scalar_determinant"],
        "metric_patch_factor": patch_factor,
        "metric_patch_nonzero_constant": patch_quotient[0],
        "metric_patch_speed": metric_patch_point,
        "metric_null_factor": metric_null,
        "regulator_factor": regulator,
        "scalar_determinant_factorization_exact": schur_rr["scalar_determinant"]
        == poly_mul(metric_null, regulator),
        "metric_null_certificate": metric_null_certificate,
        "regulator_certificate": regulator_certificate,
        "regulator_cone_proportional_to_metric": cones_proportional,
        "cone_resultant": cone_resultant,
        "cones_share_no_root": cone_resultant != 0,
        "metric_null_nonzero_at_patch": metric_null_at_patch,
        "regulator_nonzero_at_patch": regulator_at_patch,
        "alternate_metric_block_nonzero_at_patch": alternate_metric_block_at_patch,
        "quotient_atlas_covers_patch_boundary": alternate_metric_block_at_patch != 0,
        "metric_block_invertible_at_physical_roots": (
            metric_null_at_patch != 0
            and regulator_at_patch != 0
            and alternate_metric_block_at_patch != 0
        ),
        "adm_normal_preflight": adm,
    }


def required_sym1_nonclaims() -> dict[str, bool]:
    return {
        name: False
        for name in (
            "nonlinear_first_order_evolution_system_derived",
            "lower_order_radial_constraints_derived",
            "gauge_preservation_proven",
            "constraint_propagation_proven",
            "regular_center_boundary_system_derived",
            "open_retained_eft_domain_proven",
            "uniform_full_system_symmetrizer_proven",
            "full_strong_hyperbolicity_proven",
            "evolution_authorized",
            "collapse_solution_derived",
            "metric_null_affine_defocusing_derived",
        )
    }


def spherical_symbol_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    """Build the aggregate exact SYM1 certificate from RED1 fixtures."""

    if not isinstance(configuration, Mapping):
        raise ValueError("configuration must be a mapping")
    if set(configuration) != {"fixtures"}:
        raise ValueError("SYM1 configuration adapter accepts only fixtures")
    fixtures = configuration["fixtures"]
    if not isinstance(fixtures, (list, tuple)) or not fixtures:
        raise ValueError("SYM1 requires a nonempty fixture list")
    records = [symbol_fixture_certificate(fixture) for fixture in fixtures]
    gr_records = [record for record in records if record["fixture_id"] == "GR0_regular"]
    if len(gr_records) != 1:
        raise ValueError("SYM1 requires exactly one GR0_regular fixture")
    gr_fixture = next(
        fixture for fixture in fixtures if fixture.get("fixture_id") == "GR0_regular"
    )
    naive = naive_fixed_gauge_comparator(gr_fixture)
    verified = {
        "right_gauge_null_identities": all(
            item["right_gauge_null_identity_exact"] for item in records
        ),
        "left_bianchi_null_identities": all(
            item["left_bianchi_null_identity_exact"] for item in records
        ),
        "metric_schur_identities": all(
            item["quotient_atlas"][patch]["metric_schur"][
                "schur_determinant_identity"
            ]
            for item in records
            for patch in ("rr_patch", "tt_patch")
        ),
        "scalar_schur_offdiagonal_zero": all(
            item["quotient_atlas"][patch]["metric_schur"][
                "scalar_offdiagonal_numerators_zero"
            ]
            for item in records
            for patch in ("rr_patch", "tt_patch")
        ),
        "quotient_atlas_scalar_determinants_agree": all(
            item["quotient_atlas_scalar_determinants_agree_exactly"]
            and item["quotient_atlas_covers_patch_boundary"]
            for item in records
        ),
        "scalar_determinant_factorization": all(
            item["scalar_determinant_factorization_exact"] for item in records
        ),
        "metric_blocks_invertible_at_physical_roots": all(
            item["metric_block_invertible_at_physical_roots"] for item in records
        ),
        "positive_metric_null_discriminants": all(
            item["metric_null_certificate"]["discriminant"] > 0 for item in records
        ),
        "positive_regulator_discriminants": all(
            item["regulator_certificate"]["discriminant"] > 0 for item in records
        ),
        "exact_positive_scalar_companion_symmetrizers": all(
            item["metric_null_certificate"]["product_is_symmetric"]
            and item["metric_null_certificate"]["symmetrizer_positive_definite"]
            and item["regulator_certificate"]["product_is_symmetric"]
            and item["regulator_certificate"]["symmetrizer_positive_definite"]
            for item in records
        ),
        "candidate_constraints_free_of_normal_accelerations": all(
            item["adm_normal_preflight"][
                "candidate_constraints_free_of_normal_accelerations"
            ]
            for item in records
        ),
        "candidate_constraints_free_of_gauge_source_second_jets": all(
            item["adm_normal_preflight"][
                "candidate_constraints_free_of_gauge_source_second_jets"
            ]
            for item in records
        ),
        "pointwise_metric_phi_kinetic_blocks_nonsingular": all(
            item["adm_normal_preflight"]["kinetic_determinant"] != 0
            for item in records
        ),
        "naive_fixed_gauge_gr_comparator_rejected": naive[
            "strong_hyperbolicity_rejected"
        ],
    }
    return {
        "artifact_id": "FGC-1-HYP1-SYM1",
        "classification": "exact_pointwise_covariant_scalar_quotient_and_adm_formulation_preflight_not_full_hyperbolicity_or_evolution",
        "fixture_records": records,
        "naive_fixed_gauge_gr_comparator": naive,
        "verified_exact_checks": verified,
        "all_declared_exact_checks_pass": all(verified.values()),
        "nonclaims": required_sym1_nonclaims(),
    }
