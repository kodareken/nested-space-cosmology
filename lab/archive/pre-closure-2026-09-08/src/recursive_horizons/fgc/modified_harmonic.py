"""Exact spherical principal gate for the ACT1 modified-harmonic formulation.

The physical equations are the unredefined VAR1 equations.  Following
Kovacs--Reall, the metric equation is extended off the gauge surface by a
modified-harmonic gauge-fixing principal term built from two auxiliary inverse
metrics.  This module proves pointwise, exact statements about that complete
spherical second-order symbol and its standard first-order principal
reduction.  It does not yet provide lower-order source terms, an open-domain
symmetrizer, boundary data, or an evolution.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from typing import Any, Mapping, Sequence

from .exact_interval import Interval, interval_matrix_determinant
from .exact_tangent import FirstTangent
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
    polynomial_matrix_det_bareiss,
)
from .interval_tangent import IntervalFirstTangent
from .spherical_reduction import BASE_FIELD_ORDER, SphericalState
from .spherical_symbol import (
    _polynomial_matrix_zero,
    _quadratic_resultant,
    quadratic_companion_certificate,
    radial_symbol,
    symbol_fixture_certificate,
)


Q = Fraction
N = 4

MHG_FIELD_ORDER = BASE_FIELD_ORDER
MHG_EQUATION_ORDER = (
    "metric_tt_mhg",
    "metric_tr_mhg",
    "metric_rr_mhg",
    "metric_theta_theta_mhg",
    "scalar_phi",
    "scalar_chi",
)
MHG_FIRST_ORDER_STATE_ORDER = tuple(
    f"dt_{name}" for name in MHG_FIELD_ORDER
) + tuple(f"dr_{name}" for name in MHG_FIELD_ORDER)
MHG_GAUGE_CONSTRAINT_ORDER = ("H^t", "H^r")


def _generic_four_by_four_determinant(value: Sequence[Sequence[Any]]) -> Any:
    """Return a scalar-generic 4x4 determinant for adapter validation.

    The exact public path continues to use ``matrix_det``.  A numerical
    ``Jet2`` adapter can carry a floating first-tangent whose primal is a
    binary64 scalar; this small Leibniz check avoids routing that scalar through
    the exact-linear-algebra validator while evaluating the same auxiliary
    metric definition.
    """

    if len(value) != N or any(len(row) != N for row in value):
        raise ValueError("generic determinant requires a four by four matrix")
    total: Any = 0
    for permutation in permutations(range(N)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(N)
            for right in range(left + 1, N)
        )
        term: Any = 1
        for row, column in enumerate(permutation):
            term = term * value[row][column]
        total = total - term if inversions % 2 else total + term
    return total


def _physical_metric(state: SphericalState) -> tuple[Matrix, Matrix]:
    """Return the equatorial physical metric and inverse exactly."""

    radius_squared = state.areal_radius.value**2
    metric = (
        (state.h_tt.value, state.h_tr.value, Q(0), Q(0)),
        (state.h_tr.value, state.h_rr.value, Q(0), Q(0)),
        (Q(0), Q(0), radius_squared, Q(0)),
        (Q(0), Q(0), Q(0), radius_squared),
    )
    determinant = (
        state.h_tt.value * state.h_rr.value - state.h_tr.value**2
    )
    inverse = (
        (state.h_rr.value / determinant, -state.h_tr.value / determinant, Q(0), Q(0)),
        (-state.h_tr.value / determinant, state.h_tt.value / determinant, Q(0), Q(0)),
        (Q(0), Q(0), 1 / radius_squared, Q(0)),
        (Q(0), Q(0), Q(0), 1 / radius_squared),
    )
    if matrix_mul(metric, inverse) != tuple(
        tuple(Q(int(row == column)) for column in range(N)) for row in range(N)
    ):
        raise ValueError("physical metric inverse failed exact reconstruction")
    return metric, inverse


def auxiliary_inverse_metric(
    physical_inverse: Sequence[Sequence[Fraction]],
    normal_factor: Fraction | int,
) -> Matrix:
    """Return ``g^ab-(q-1)n^a n^b`` without a square root.

    Here ``n`` is the future unit normal to the coordinate-time slices and
    ``q>1``.  Since ``n^a n^b=-g^{a0}g^{b0}/g^{00}``, the construction is
    exact rational whenever the physical inverse metric is exact rational.
    Its cotangent null cone lies strictly outside the physical one.
    """

    factor = Q(normal_factor)
    if factor <= 1:
        raise ValueError("auxiliary normal factor must be greater than one")
    # IMP1 originally needed tangents only in second-jet slots.  FO1-RC1 also
    # differentiates the complete residual with respect to field values and
    # first derivatives, so the physical inverse can carry an exact
    # ``FirstTangent``.  Preserve that algebraic tangent while retaining the
    # strict rational coercion for the public ordinary-rational path.
    inverse = tuple(
        tuple(
            entry
            if isinstance(entry, (FirstTangent, Interval, IntervalFirstTangent))
            else Q(entry)
            for entry in row
        )
        for row in physical_inverse
    )
    if len(inverse) != N or any(len(row) != N for row in inverse):
        raise ValueError("physical inverse metric must be four by four")
    inverse_tt = inverse[0][0]
    inverse_tt_primal = (
        inverse_tt.primal
        if isinstance(inverse_tt, (FirstTangent, IntervalFirstTangent))
        else inverse_tt
    )
    if isinstance(inverse_tt_primal, Interval):
        if not inverse_tt_primal.strictly_negative():
            raise ValueError("coordinate-time covector must be physical timelike")
    elif inverse_tt_primal >= 0:
        raise ValueError("coordinate-time covector must be physical timelike")
    normal_outer = tuple(
        tuple(-inverse[mu][0] * inverse[nu][0] / inverse[0][0] for nu in range(N))
        for mu in range(N)
    )
    result = tuple(
        tuple(
            inverse[mu][nu] - (factor - 1) * normal_outer[mu][nu]
            for nu in range(N)
        )
        for mu in range(N)
    )
    primal_result = tuple(
        tuple(
            entry.primal
            if isinstance(entry, (FirstTangent, IntervalFirstTangent))
            else entry
            for entry in row
        )
        for row in result
    )
    if any(isinstance(entry, Interval) for row in primal_result for entry in row):
        determinant = interval_matrix_determinant(primal_result)
        if not determinant.strictly_negative():
            raise ValueError("auxiliary inverse metric must remain Lorentzian")
    else:
        determinant = (
            matrix_det(primal_result)
            if all(
                isinstance(entry, Fraction)
                for row in primal_result
                for entry in row
            )
            else _generic_four_by_four_determinant(primal_result)
        )
        if determinant >= 0:
            raise ValueError("auxiliary inverse metric must remain Lorentzian")
    return result


def _trace_reversal_projector(
    inverse: Matrix,
    alpha: int,
    beta: int,
    mu: int,
    nu: int,
) -> Fraction:
    """Return ``P_alpha^(beta mu nu)`` for one inverse metric."""

    return (
        Q(int(alpha == mu)) * inverse[nu][beta]
        + Q(int(alpha == nu)) * inverse[mu][beta]
        - Q(int(alpha == beta)) * inverse[mu][nu]
    ) / 2


def _radial_covector() -> tuple[Polynomial, ...]:
    return (poly((0, -1)), ONE, ZERO, ZERO)


def _metric_perturbation_basis(state: SphericalState) -> tuple[Matrix, ...]:
    radius = state.areal_radius.value
    zero = [[Q(0) for _ in range(N)] for _ in range(N)]
    basis: list[Matrix] = []
    for entries in (
        ((0, 0, Q(1)),),
        ((0, 1, Q(1)), (1, 0, Q(1))),
        ((1, 1, Q(1)),),
        ((2, 2, 2 * radius), (3, 3, 2 * radius)),
    ):
        value = [row[:] for row in zero]
        for row, column, coefficient in entries:
            value[row][column] = coefficient
        basis.append(tuple(tuple(row) for row in value))
    return tuple(basis)


def modified_harmonic_gauge_symbol(
    state: SphericalState,
    *,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> tuple[tuple[Polynomial, ...], ...]:
    """Return the lower-index spherical MHG gauge-fixing principal matrix.

    This is the exact spherical restriction of Kovacs--Reall Eq. (15), scaled
    by the positive ACT1 Einstein coefficient ``F``.  Rows are the four lower
    metric equations followed by two zero scalar rows; columns use RED1's six
    base fields.
    """

    metric, physical_inverse = _physical_metric(state)
    tilde = auxiliary_inverse_metric(physical_inverse, tilde_normal_factor)
    hat = auxiliary_inverse_metric(physical_inverse, hat_normal_factor)
    xi = _radial_covector()
    effective_planck = state.planck_mass**2 + state.beta * state.phi.value**2
    if effective_planck <= 0:
        raise ValueError("modified-harmonic gauge scaling requires F>0")

    # Contract each projector with xi once.  The remaining physical inverse
    # contraction is the g^(alpha beta) in Kovacs--Reall Eq. (15).
    hat_xi = [
        [
            [
                _poly_sum(
                    poly_scale(
                        xi[gamma],
                        _trace_reversal_projector(hat, alpha, gamma, mu, nu),
                    )
                    for gamma in range(N)
                )
                for nu in range(N)
            ]
            for mu in range(N)
        ]
        for alpha in range(N)
    ]
    tilde_xi = [
        [
            [
                _poly_sum(
                    poly_scale(
                        xi[delta],
                        _trace_reversal_projector(
                            tilde, beta, delta, rho, sigma
                        ),
                    )
                    for delta in range(N)
                )
                for sigma in range(N)
            ]
            for rho in range(N)
        ]
        for beta in range(N)
    ]

    contravariant: list[list[list[Polynomial]]] = []
    for perturbation in _metric_perturbation_basis(state):
        response = [[ZERO for _ in range(N)] for _ in range(N)]
        for mu in range(N):
            for nu in range(N):
                terms: list[Polynomial] = []
                for alpha in range(N):
                    for beta in range(N):
                        if physical_inverse[alpha][beta] == 0:
                            continue
                        for rho in range(N):
                            for sigma in range(N):
                                if perturbation[rho][sigma] == 0:
                                    continue
                                terms.append(
                                    poly_scale(
                                        poly_mul(
                                            hat_xi[alpha][mu][nu],
                                            tilde_xi[beta][rho][sigma],
                                        ),
                                        -effective_planck
                                        * physical_inverse[alpha][beta]
                                        * perturbation[rho][sigma],
                                    )
                                )
                response[mu][nu] = _poly_sum(terms)
        contravariant.append(response)

    output = [[ZERO for _ in MHG_FIELD_ORDER] for _ in MHG_EQUATION_ORDER]
    lower_pairs = ((0, 0), (0, 1), (1, 1), (2, 2))
    for column, response in enumerate(contravariant):
        for row, (a, b) in enumerate(lower_pairs):
            output[row][column] = _poly_sum(
                poly_scale(
                    response[mu][nu], metric[a][mu] * metric[b][nu]
                )
                for mu in range(N)
                for nu in range(N)
                if metric[a][mu] and metric[b][nu]
            )
    return tuple(tuple(row) for row in output)


def _poly_sum(values) -> Polynomial:
    result = ZERO
    for value in values:
        result = poly_add(result, value)
    return result


def modified_harmonic_symbol(
    state: SphericalState,
    *,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> tuple[tuple[Polynomial, ...], ...]:
    """Return the complete unredefined-equations-plus-MHG radial symbol."""

    covariant = radial_symbol(state)
    gauge = modified_harmonic_gauge_symbol(
        state,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    return tuple(
        tuple(poly_add(covariant[row][column], gauge[row][column]) for column in range(6))
        for row in range(6)
    )


def inverse_metric_null_polynomial(inverse: Matrix) -> Polynomial:
    """Return ``g^AB xi_A xi_B`` for ``xi_A=(-c,1)``."""

    return poly((inverse[1][1], -2 * inverse[0][1], inverse[0][0]))


def gauge_constraint_propagation_principal_symbol(
    state: SphericalState,
    *,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Derive the principal divergence of the MHG gauge extension exactly.

    At principal order, on the scalar equations, the ACT1 Noether identity
    removes the divergence of the unredefined metric equation.  What remains
    is the double-covector contraction

    ``F xi_mu hat_P_alpha^(beta mu nu) xi_beta H^alpha``.

    The trace-reversal projector makes the two longitudinal terms cancel,
    leaving ``(F/2) hat_g^(rho sigma) xi_rho xi_sigma delta_alpha^nu``.
    This routine performs the full four-index contraction rather than
    inserting that expected result.
    """

    _, physical_inverse = _physical_metric(state)
    hat = auxiliary_inverse_metric(physical_inverse, hat_normal_factor)
    xi = _radial_covector()
    effective_planck = state.planck_mass**2 + state.beta * state.phi.value**2
    if effective_planck <= 0:
        raise ValueError("gauge-constraint propagation requires F>0")

    derived = tuple(
        tuple(
            _poly_sum(
                poly_scale(
                    poly_mul(xi[mu], xi[beta]),
                    effective_planck
                    * _trace_reversal_projector(
                        hat, alpha, beta, mu, nu
                    ),
                )
                for mu in range(N)
                for beta in range(N)
            )
            for alpha in range(N)
        )
        for nu in range(N)
    )
    expected_wave_factor = poly_scale(
        inverse_metric_null_polynomial(hat), effective_planck / 2
    )
    expected = tuple(
        tuple(
            expected_wave_factor if nu == alpha else ZERO
            for alpha in range(N)
        )
        for nu in range(N)
    )
    spherical_indices = (0, 1)
    spherical = tuple(
        tuple(derived[nu][alpha] for alpha in spherical_indices)
        for nu in spherical_indices
    )
    return {
        "input_constraint_components": ("H^t", "H^r", "H^theta", "H^phi"),
        "output_divergence_components": ("nu=t", "nu=r", "nu=theta", "nu=phi"),
        "full_projector_contraction": derived,
        "expected_full_hat_wave_operator": expected,
        "expected_hat_wave_factor": expected_wave_factor,
        "expected_hat_wave_factor_monic": monic(expected_wave_factor),
        "spherical_projector_contraction": spherical,
        "full_projector_contraction_exact": derived == expected,
    }


def monic(value: Sequence[Fraction]) -> Polynomial:
    polynomial_value = poly(value)
    if polynomial_value == ZERO:
        raise ValueError("zero polynomial has no monic normalization")
    return poly_scale(polynomial_value, 1 / polynomial_value[-1])


def _polynomial_coefficient_matrix(
    symbol: Sequence[Sequence[Sequence[Fraction]]], degree: int
) -> Matrix:
    return tuple(
        tuple(entry[degree] if degree < len(entry) else Q(0) for entry in row)
        for row in symbol
    )


def standard_first_order_principal_system(
    symbol: Sequence[Sequence[Sequence[Fraction]]],
) -> dict[str, Any]:
    """Return the standard ``(u_t,u_r)`` first-order principal reduction."""

    if len(symbol) != 6 or any(len(row) != 6 for row in symbol):
        raise ValueError("MHG second-order symbol must be six by six")
    a_tt = _polynomial_coefficient_matrix(symbol, 2)
    a_tr = tuple(
        tuple(-entry for entry in row)
        for row in _polynomial_coefficient_matrix(symbol, 1)
    )
    a_rr = _polynomial_coefficient_matrix(symbol, 0)
    if matrix_det(a_tt) == 0:
        raise ValueError("MHG coordinate-time kinetic matrix is singular")
    zero6 = tuple(tuple(Q(0) for _ in range(6)) for _ in range(6))
    identity6 = tuple(
        tuple(Q(int(row == column)) for column in range(6)) for row in range(6)
    )
    minus_identity6 = tuple(
        tuple(-entry for entry in row) for row in identity6
    )
    at = tuple(
        tuple(a_tt[row] + zero6[row]) for row in range(6)
    ) + tuple(tuple(zero6[row] + identity6[row]) for row in range(6))
    ar = tuple(
        tuple(a_tr[row] + a_rr[row]) for row in range(6)
    ) + tuple(tuple(minus_identity6[row] + zero6[row]) for row in range(6))
    principal = matrix_mul(matrix_inverse(at), ar)
    characteristic_pencil = tuple(
        tuple(poly((ar[row][column], -at[row][column])) for column in range(12))
        for row in range(12)
    )
    characteristic_determinant = polynomial_matrix_det_bareiss(
        characteristic_pencil
    )
    return {
        "state_order": MHG_FIRST_ORDER_STATE_ORDER,
        "A_tt": a_tt,
        "A_tr": a_tr,
        "A_rr": a_rr,
        "A_time": at,
        "A_radial": ar,
        "characteristic_pencil_determinant": characteristic_determinant,
        "coordinate_time_kinetic_determinant": matrix_det(a_tt),
        "A_time_determinant": matrix_det(at),
        "normal_principal_matrix": principal,
    }


def _matrix_at(
    symbol: Sequence[Sequence[Sequence[Fraction]]], value: Fraction
) -> Matrix:
    return tuple(tuple(poly_eval(entry, value) for entry in row) for row in symbol)


def _first_order_pencil_at(
    first_order: Mapping[str, Any], value: Fraction
) -> Matrix:
    a_time = first_order["A_time"]
    a_radial = first_order["A_radial"]
    return tuple(
        tuple(
            a_radial[row][column] - value * a_time[row][column]
            for column in range(12)
        )
        for row in range(12)
    )


def _rational_quadratic_roots(value: Polynomial) -> tuple[Fraction, Fraction]:
    certificate = quadratic_companion_certificate(value)
    roots = certificate["roots"]
    if any("exact" not in root for root in roots):
        raise ValueError("quadratic roots are not rational")
    return roots[0]["exact"], roots[1]["exact"]


def modified_harmonic_fixture_certificate(
    fixture: Mapping[str, Any],
    *,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Build one exact MHG radial-principal fixture certificate."""

    from .spherical_reduction import state_from_generalized_adm_pg_fixture

    state = state_from_generalized_adm_pg_fixture(fixture)
    metric, physical_inverse = _physical_metric(state)
    tilde = auxiliary_inverse_metric(physical_inverse, tilde_normal_factor)
    hat = auxiliary_inverse_metric(physical_inverse, hat_normal_factor)
    gauge = modified_harmonic_gauge_symbol(
        state,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    full = modified_harmonic_symbol(
        state,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    determinant = polynomial_matrix_det_bareiss(full)
    tilde_null_raw = inverse_metric_null_polynomial(tilde)
    hat_null_raw = inverse_metric_null_polynomial(hat)
    physical_null_raw = inverse_metric_null_polynomial(physical_inverse)
    tilde_null = monic(tilde_null_raw)
    hat_null = monic(hat_null_raw)
    physical_null = monic(physical_null_raw)
    sym1 = symbol_fixture_certificate(fixture)
    scalar_determinant = monic(sym1["quotient_atlas"]["rr_patch"]["metric_schur"]["scalar_determinant"])
    regulator = monic(sym1["regulator_factor"])
    expected_factors = poly_mul(
        poly_mul(poly_mul(tilde_null, tilde_null), poly_mul(hat_null, hat_null)),
        scalar_determinant,
    )
    quotient = poly_div_exact(determinant, expected_factors)
    if len(quotient) != 1 or quotient[0] == 0:
        raise ValueError("MHG determinant lacks the declared nonzero factorization")

    first_order = standard_first_order_principal_system(full)
    rational_modes: dict[str, Any] = {}
    for name, factor, expected_kernel in (
        ("tilde_gauge", tilde_null, 2),
        ("hat_constraint", hat_null, 2),
        (
            "physical_metric",
            physical_null,
            2 if sym1["regulator_cone_proportional_to_metric"] else 1,
        ),
    ):
        roots = _rational_quadratic_roots(factor)
        second_order_kernels = tuple(
            6 - matrix_rank(_matrix_at(full, root)) for root in roots
        )
        first_order_kernels = tuple(
            12 - matrix_rank(_first_order_pencil_at(first_order, root))
            for root in roots
        )
        rational_modes[name] = {
            "factor": factor,
            "roots": roots,
            "kernel_dimensions": first_order_kernels,
            "second_order_kernel_dimensions": second_order_kernels,
            "first_order_kernel_dimensions": first_order_kernels,
            "second_to_first_order_kernel_map_exact": (
                second_order_kernels == first_order_kernels
            ),
            "expected_kernel_dimension": expected_kernel,
            "semisimple_at_each_rational_root": all(
                value == expected_kernel for value in first_order_kernels
            )
            and second_order_kernels == first_order_kernels,
        }

    propagation = gauge_constraint_propagation_principal_symbol(
        state,
        hat_normal_factor=hat_normal_factor,
    )
    physical_roots = rational_modes["physical_metric"]["roots"]
    cone_nesting = {
        "physical_roots_are_tilde_timelike": all(
            poly_eval(tilde_null_raw, root) < 0 for root in physical_roots
        ),
        "physical_roots_are_hat_timelike": all(
            poly_eval(hat_null_raw, root) < 0 for root in physical_roots
        ),
        "auxiliary_cones_do_not_intersect": _quadratic_resultant(
            tilde_null, hat_null
        )
        != 0,
    }
    x_zero_regular = (
        state.phi.dt == 0
        and state.phi.dr == 0
        and all(
            all(all(isinstance(coefficient, Fraction) for coefficient in entry) for entry in row)
            for row in full
        )
    )
    cone_resultants = {
        "tilde_vs_hat": _quadratic_resultant(tilde_null, hat_null),
        "tilde_vs_physical": _quadratic_resultant(tilde_null, physical_null),
        "hat_vs_physical": _quadratic_resultant(hat_null, physical_null),
        "tilde_vs_regulator": _quadratic_resultant(tilde_null, regulator),
        "hat_vs_regulator": _quadratic_resultant(hat_null, regulator),
        "physical_vs_regulator": _quadratic_resultant(
            physical_null, regulator
        ),
    }
    auxiliary_separated_from_all_physical_modes = all(
        cone_resultants[name] != 0
        for name in (
            "tilde_vs_hat",
            "tilde_vs_physical",
            "hat_vs_physical",
            "tilde_vs_regulator",
            "hat_vs_regulator",
        )
    )
    rational_repeated_modes_semisimple = all(
        rational_modes[name]["semisimple_at_each_rational_root"]
        for name in ("tilde_gauge", "hat_constraint", "physical_metric")
    )
    scalar_characteristics_real = (
        sym1["metric_null_certificate"]["discriminant"] > 0
        and sym1["regulator_certificate"]["discriminant"] > 0
    )
    physical_mode_multiplicity_closed = (
        cone_resultants["physical_vs_regulator"] == 0
        and sym1["regulator_cone_proportional_to_metric"]
        and rational_modes["physical_metric"][
            "semisimple_at_each_rational_root"
        ]
    ) or (
        cone_resultants["physical_vs_regulator"] != 0
        and not sym1["regulator_cone_proportional_to_metric"]
    )
    return {
        "fixture_id": fixture.get("fixture_id"),
        "model_id": fixture.get("model_id"),
        "effective_planck_coefficient": state.planck_mass**2
        + state.beta * state.phi.value**2,
        "physical_metric": metric,
        "physical_inverse_metric": physical_inverse,
        "tilde_inverse_metric": tilde,
        "hat_inverse_metric": hat,
        "gauge_fixing_symbol": gauge,
        "gauge_fixing_nonzero": not _polynomial_matrix_zero(gauge),
        "complete_mhg_symbol": full,
        "complete_mhg_determinant": determinant,
        "determinant_factorization": {
            "tilde_null_monic": tilde_null,
            "tilde_multiplicity": 2,
            "hat_null_monic": hat_null,
            "hat_multiplicity": 2,
            "physical_scalar_determinant_monic": scalar_determinant,
            "nonzero_constant_quotient": quotient[0],
            "exact": determinant == poly_scale(expected_factors, quotient[0]),
        },
        "rational_characteristic_modes": rational_modes,
        "auxiliary_cone_checks": cone_nesting,
        "cone_resultants": cone_resultants,
        "auxiliary_cones_separated_from_all_physical_modes": (
            auxiliary_separated_from_all_physical_modes
        ),
        "gauge_constraint_propagation_principal_symbol": propagation,
        "gauge_constraint_propagation_hat_null_exact": (
            propagation["full_projector_contraction_exact"]
            and propagation["expected_hat_wave_factor_monic"] == hat_null
        ),
        "standard_first_order_principal_system": first_order,
        "first_order_characteristic_determinant_equals_second_order": (
            first_order["characteristic_pencil_determinant"] == determinant
        ),
        "coordinate_time_kinetic_matches_characteristic_leading_coefficient": (
            first_order["coordinate_time_kinetic_determinant"]
            == determinant[-1]
        ),
        "pointwise_complete_real_first_order_eigenbasis": (
            scalar_characteristics_real
            and auxiliary_separated_from_all_physical_modes
            and rational_repeated_modes_semisimple
            and physical_mode_multiplicity_closed
        ),
        "x_zero_direct_var1_coefficients_regular": x_zero_regular,
        "sym1_physical_regression": {
            "physical_metric_null_factor": sym1["metric_null_factor"],
            "regulator_factor": sym1["regulator_factor"],
            "scalar_determinant_factorization_exact": sym1[
                "scalar_determinant_factorization_exact"
            ],
        },
    }


def required_mhg1_nonclaims() -> dict[str, bool]:
    return {
        name: False
        for name in (
            "complete_lower_order_first_order_sources_derived",
            "gauge_constraint_lower_order_coefficients_serialized",
            "reduction_constraint_propagation_proven",
            "uniform_open_domain_symmetrizer_proven",
            "open_retained_eft_domain_proven",
            "boundary_or_regular_center_system_derived",
            "evolution_authorized",
            "collapse_solution_derived",
            "metric_null_affine_defocusing_derived",
            "singularity_resolution_derived",
        )
    }


def modified_harmonic_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    """Build the aggregate exact MHG1 pointwise formulation certificate."""

    if not isinstance(configuration, Mapping):
        raise ValueError("MHG1 configuration must be a mapping")
    if set(configuration) != {
        "fixtures",
        "tilde_normal_factor",
        "hat_normal_factor",
    }:
        raise ValueError("MHG1 configuration adapter has unexpected keys")
    fixtures = configuration["fixtures"]
    if not isinstance(fixtures, (list, tuple)) or not fixtures:
        raise ValueError("MHG1 requires a nonempty fixture list")
    tilde_factor = Q(configuration["tilde_normal_factor"])
    hat_factor = Q(configuration["hat_normal_factor"])
    if not 1 < tilde_factor < hat_factor:
        raise ValueError("MHG1 requires 1 < tilde factor < hat factor")
    records = [
        modified_harmonic_fixture_certificate(
            fixture,
            tilde_normal_factor=tilde_factor,
            hat_normal_factor=hat_factor,
        )
        for fixture in fixtures
    ]
    x_zero_records = [
        record for record in records if record["x_zero_direct_var1_coefficients_regular"]
    ]
    verified = {
        "positive_effective_planck_coefficients": all(
            record["effective_planck_coefficient"] > 0 for record in records
        ),
        "gauge_fixing_symbols_nonzero": all(
            record["gauge_fixing_nonzero"] for record in records
        ),
        "complete_determinant_factorizations": all(
            record["determinant_factorization"]["exact"] for record in records
        ),
        "auxiliary_cones_nested_and_separated": all(
            all(record["auxiliary_cone_checks"].values())
            and record["auxiliary_cones_separated_from_all_physical_modes"]
            for record in records
        ),
        "gauge_constraint_hat_cone_principal_propagation": all(
            record["gauge_constraint_propagation_hat_null_exact"]
            for record in records
        ),
        "coordinate_time_kinetic_matrices_invertible": all(
            record["standard_first_order_principal_system"][
                "coordinate_time_kinetic_determinant"
            ]
            != 0
            for record in records
        ),
        "first_order_characteristic_reductions_exact": all(
            record[
                "first_order_characteristic_determinant_equals_second_order"
            ]
            and record[
                "coordinate_time_kinetic_matches_characteristic_leading_coefficient"
            ]
            for record in records
        ),
        "pointwise_complete_real_first_order_eigenbases": all(
            record["pointwise_complete_real_first_order_eigenbasis"]
            for record in records
        ),
        "direct_x_zero_var1_regularity_controls": bool(x_zero_records),
        "sym1_physical_factor_regression": all(
            record["sym1_physical_regression"][
                "scalar_determinant_factorization_exact"
            ]
            for record in records
        ),
    }
    return {
        "artifact_id": "FGC-1-HYP1-MHG1",
        "classification": "exact_pointwise_modified_harmonic_spherical_principal_and_gauge_constraint_preflight_not_open_domain_or_evolution",
        "formulation": {
            "physical_equations": "unredefined_ACT1_VAR1",
            "gauge_condition": "H^mu=-tilde_g^rho_sigma*Gamma^mu_rho_sigma=0_up_to_declared_reference_connection",
            "metric_extension": "E_mhg^mu_nu=E^mu_nu+F*hat_P_alpha^beta_mu_nu*partial_beta_H^alpha",
            "tilde_normal_factor": tilde_factor,
            "hat_normal_factor": hat_factor,
            "gauge_surface_equivalence": "H_and_first_derivative_zero_implies_E_mhg_equals_E",
            "gauge_constraint_propagation": "homogeneous_hat_metric_wave_principal_on_scalar_equations_by_Noether_identity",
            "first_order_reduction": "standard_dt_dr_reduction_of_complete_six_field_spherical_second_order_system",
        },
        "fixture_records": records,
        "verified_exact_checks": verified,
        "all_declared_exact_checks_pass": all(verified.values()),
        "nonclaims": required_mhg1_nonclaims(),
    }
