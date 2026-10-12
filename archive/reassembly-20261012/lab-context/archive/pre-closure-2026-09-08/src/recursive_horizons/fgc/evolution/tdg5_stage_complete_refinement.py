"""Exact no-history controls for the TDG5 temporal replacement design.

TDG4 proved that finitely many point values do not determine an unrestricted
smooth history between those points.  TDG5 therefore changes the numerical
object being admitted.  Each accepted time interval is declared to be a
finite-dimensional cubic Hermite continuous extension, determined by endpoint
values and endpoint right-hand sides.  A same-grid one-step/two-half-step pair
then measures temporal refinement without changing the spatial operator.

This module proves only the algebra needed to freeze that design.  It performs
no I/O, reads no campaign history, advances no state, defines no runtime
threshold, and grants no GR-0 or candidate authorization.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from numbers import Rational
from typing import Iterable, Sequence


Q = Fraction

TDG5_METHOD_ORDERS = (("RK4", 4), ("SSPRK3", 3))
TDG5_HERMITE_DATA_ORDER = (
    "left_value",
    "left_scaled_slope",
    "right_value",
    "right_scaled_slope",
)
TDG5_OWNED_ROW_POLICY = (
    "dimensionless_u_p_q_on_noncentre_rows_excluding_four_projector_owned_outer_rows"
)
TDG5_SELECTED_SUFFICIENT_ROUTE = (
    "finite_dimensional_stage_complete_continuous_extension_with_stable_injective_data_map"
)
TDG5_REQUIRED_RUNTIME_LAYERS = (
    "one_full_step_and_two_half_steps_share_one_initial_state_and_spatial_operator",
    "continuous_extension_uses_complete_state_endpoint_values_and_fresh_endpoint_RHS_slopes",
    "every_internal_stage_and_candidate_endpoint_retains_unchanged_source_health_CFL_boundary_and_transaction_guards",
    "coarse_and_fine_proposals_each_pass_unchanged_source_health_CFL_boundary_and_transaction_guards",
    "only_the_fine_two_half_step_path_may_commit_atomically",
    "continuous_coarse_fine_difference_is_maximized_on_each_half_interval",
    "method_order_owned_Richardson_debit_is_recorded_without_becoming_a_rigorous_PDE_bound",
    "both_methods_and_method_owned_spatial_ladders_remain_independent_requirements",
    "projector_owned_rows_are_covered_by_existing_exact_projector_and_constraint_gates",
)


def _exact(value: Rational, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(f"{name} must be an exact rational")
    return Fraction(value)


def _fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _matrix(
    value: Iterable[Iterable[Rational]], *, name: str
) -> tuple[tuple[Fraction, ...], ...]:
    rows = tuple(tuple(_exact(item, name=name) for item in row) for row in value)
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError(f"{name} must be a nonempty rectangular matrix")
    return rows


def _matmul(
    left: Sequence[Sequence[Fraction]],
    right: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("matrix product shapes differ")
    return tuple(
        tuple(
            sum(
                (
                    left[row][inner] * right[inner][column]
                    for inner in range(len(right))
                ),
                Q(0),
            )
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def _determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    size = len(matrix)
    if size < 1 or any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a nonempty square matrix")
    answer = Q(0)
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(size)
            for right in range(left + 1, size)
        )
        product = Q(-1 if inversions % 2 else 1)
        for row, column in enumerate(permutation):
            product *= matrix[row][column]
        answer += product
    return answer


def _infinity_norm(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    return max(sum((abs(item) for item in row), Q(0)) for row in matrix)


def _identity(size: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Q(int(row == column)) for column in range(size))
        for row in range(size)
    )


def hermite_data_map_certificate() -> dict[str, object]:
    """Prove exact injectivity and a finite inverse bound for cubic data.

    For ``p(theta)=a0+a1 theta+a2 theta^2+a3 theta^3``, the data vector is
    ``(p(0), p'(0), p(1), p'(1))``.  Endpoint physical slopes are multiplied
    by the interval width before entering this normalized map.
    """

    sampling = _matrix(
        (
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (1, 1, 1, 1),
            (0, 1, 2, 3),
        ),
        name="Hermite sampling matrix",
    )
    inverse = _matrix(
        (
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (-3, -2, 3, -1),
            (2, 1, -2, 1),
        ),
        name="Hermite inverse matrix",
    )
    identity = _identity(4)
    determinant = _determinant(sampling)
    sampling_norm = _infinity_norm(sampling)
    inverse_norm = _infinity_norm(inverse)
    return {
        "polynomial_class": "P3_on_each_accepted_normalized_time_interval",
        "data_order": list(TDG5_HERMITE_DATA_ORDER),
        "sampling_matrix": [
            [_fraction_text(item) for item in row] for row in sampling
        ],
        "inverse_matrix": [
            [_fraction_text(item) for item in row] for row in inverse
        ],
        "determinant": _fraction_text(determinant),
        "left_inverse_exact": _matmul(inverse, sampling) == identity,
        "right_inverse_exact": _matmul(sampling, inverse) == identity,
        "sampling_infinity_norm": _fraction_text(sampling_norm),
        "inverse_infinity_norm": _fraction_text(inverse_norm),
        "condition_number_infinity_upper_bound": _fraction_text(
            sampling_norm * inverse_norm
        ),
        "data_map_is_injective": determinant != 0,
        "no_between_sample_kernel_inside_declared_class": determinant != 0,
        "physical_endpoint_slopes_are_scaled_by_interval_width": True,
    }


def hermite_coefficients(
    data: Sequence[Rational],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Return exact monomial coefficients from four normalized Hermite data."""

    values = tuple(_exact(item, name="Hermite datum") for item in data)
    if len(values) != 4:
        raise ValueError("Hermite data must contain exactly four entries")
    left, left_slope, right, right_slope = values
    return (
        left,
        left_slope,
        -3 * left - 2 * left_slope + 3 * right - right_slope,
        2 * left + left_slope - 2 * right + right_slope,
    )


def evaluate_cubic(
    coefficients: Sequence[Rational], theta: Rational
) -> Fraction:
    values = tuple(_exact(item, name="cubic coefficient") for item in coefficients)
    if len(values) != 4:
        raise ValueError("cubic coefficients must contain exactly four entries")
    point = _exact(theta, name="theta")
    if point < 0 or point > 1:
        raise ValueError("theta must lie in the closed unit interval")
    return values[0] + point * (
        values[1] + point * (values[2] + point * values[3])
    )


def endpoint_only_alias_control() -> dict[str, object]:
    """Exhibit a cubic missed by endpoint values but found by endpoint slopes."""

    zero = hermite_coefficients((0, 0, 0, 0))
    bump = hermite_coefficients((0, 4, 0, -4))
    midpoint = evaluate_cubic(bump, Q(1, 2))
    return {
        "zero_polynomial_coefficients": [_fraction_text(item) for item in zero],
        "adversarial_polynomial_coefficients": [
            _fraction_text(item) for item in bump
        ],
        "endpoint_values_equal": (
            evaluate_cubic(zero, Q(0)) == evaluate_cubic(bump, Q(0))
            and evaluate_cubic(zero, Q(1)) == evaluate_cubic(bump, Q(1))
        ),
        "endpoint_slopes_differ": True,
        "adversarial_midpoint_value": _fraction_text(midpoint),
        "endpoint_only_sampling_misses_nonzero_interior": midpoint != 0,
        "stage_complete_Hermite_data_detects_it": bump != zero,
    }


def _dot(left: Sequence[Fraction], right: Sequence[Fraction]) -> Fraction:
    if len(left) != len(right):
        raise ValueError("dot-product shapes differ")
    return sum((a * b for a, b in zip(left, right, strict=True)), Q(0))


def _matvec(
    matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    return tuple(_dot(row, vector) for row in matrix)


def _butcher_tableau(
    method: str,
) -> tuple[
    tuple[tuple[Fraction, ...], ...],
    tuple[Fraction, ...],
    tuple[Fraction, ...],
    int,
]:
    if method == "RK4":
        a = _matrix(
            (
                (0, 0, 0, 0),
                (Q(1, 2), 0, 0, 0),
                (0, Q(1, 2), 0, 0),
                (0, 0, 1, 0),
            ),
            name="RK4 A",
        )
        b = (Q(1, 6), Q(1, 3), Q(1, 3), Q(1, 6))
        c = (Q(0), Q(1, 2), Q(1, 2), Q(1))
        formal_order = 4
    elif method == "SSPRK3":
        a = _matrix(
            (
                (0, 0, 0),
                (1, 0, 0),
                (Q(1, 4), Q(1, 4), 0),
            ),
            name="SSPRK3 A",
        )
        b = (Q(1, 6), Q(1, 6), Q(2, 3))
        c = (Q(0), Q(1), Q(1, 2))
        formal_order = 3
    else:
        raise ValueError("unknown TDG5 method")
    return a, b, c, formal_order


def butcher_order_certificate(method: str) -> dict[str, object]:
    """Verify the classical exact Runge--Kutta order conditions in use."""

    a, b, c, formal_order = _butcher_tableau(method)

    ones = tuple(Q(1) for _ in b)
    c2 = tuple(item**2 for item in c)
    conditions = {
        "order_1_b_e": (_dot(b, ones), Q(1)),
        "order_2_b_c": (_dot(b, c), Q(1, 2)),
        "order_3_b_c2": (_dot(b, c2), Q(1, 3)),
        "order_3_b_A_c": (_dot(b, _matvec(a, c)), Q(1, 6)),
    }
    if formal_order == 4:
        ac = _matvec(a, c)
        conditions.update(
            {
                "order_4_b_c3": (
                    _dot(b, tuple(item**3 for item in c)),
                    Q(1, 4),
                ),
                "order_4_b_C_A_c": (
                    _dot(
                        b,
                        tuple(
                            c_i * ac_i
                            for c_i, ac_i in zip(c, ac, strict=True)
                        ),
                    ),
                    Q(1, 8),
                ),
                "order_4_b_A_c2": (_dot(b, _matvec(a, c2)), Q(1, 12)),
                "order_4_b_A_A_c": (
                    _dot(b, _matvec(a, ac)),
                    Q(1, 24),
                ),
            }
        )
    serialized = {
        name: {
            "observed": _fraction_text(observed),
            "required": _fraction_text(required),
            "passed": observed == required,
        }
        for name, (observed, required) in conditions.items()
    }
    return {
        "method": method,
        "formal_order": formal_order,
        "stage_count": len(b),
        "conditions": serialized,
        "all_required_order_conditions_pass": all(
            item["passed"] for item in serialized.values()
        ),
    }


def _stability_polynomial(method: str, value: Fraction) -> Fraction:
    """Evaluate the actual frozen Butcher tableau on ``y'=y`` exactly."""

    a, b, _, _ = _butcher_tableau(method)
    stages: list[Fraction] = []
    for row in range(len(b)):
        stage = Q(1) + value * sum(
            (a[row][column] * stages[column] for column in range(row)),
            Q(0),
        )
        stages.append(stage)
    return Q(1) + value * _dot(b, stages)


def step_doubling_certificate(
    method: str, *, probe_step: Rational = Q(1, 8)
) -> dict[str, object]:
    """Return exact same-grid step-doubling algebra and a scalar control."""

    certificate = butcher_order_certificate(method)
    order = int(certificate["formal_order"])
    step = _exact(probe_step, name="probe_step")
    if step <= 0:
        raise ValueError("probe_step must be positive")

    def difference(width: Fraction) -> Fraction:
        coarse = _stability_polynomial(method, width)
        half = _stability_polynomial(method, width / 2)
        return abs(half * half - coarse)

    coarse_difference = difference(step)
    half_width_difference = difference(step / 2)
    if coarse_difference == 0:
        raise RuntimeError("TDG5 scalar step-doubling control is degenerate")
    contraction = half_width_difference / coarse_difference
    denominator = 2**order - 1
    return {
        "method": method,
        "formal_order": order,
        "Richardson_denominator": denominator,
        "fine_path_error_factor": f"1/{denominator}",
        "probe_equation": "y_prime_equals_y_with_y0_equals_1",
        "probe_step": _fraction_text(step),
        "coarse_vs_two_half_difference": _fraction_text(coarse_difference),
        "half_width_difference": _fraction_text(half_width_difference),
        "exact_contraction_ratio": _fraction_text(contraction),
        "contraction_is_stricter_than_two_to_minus_formal_order": (
            contraction < Q(1, 2**order)
        ),
        "same_initial_state_and_same_spatial_operator_required": True,
        "fine_two_half_step_path_is_the_only_committable_path": True,
        "Richardson_factor_is_not_a_rigorous_PDE_error_bound": True,
    }


def stage_complete_temporal_refinement_preflight() -> dict[str, object]:
    """Assemble the exact prospective TDG5 design controls."""

    hermite = hermite_data_map_certificate()
    endpoint_attack = endpoint_only_alias_control()
    methods = {
        method: {
            "order": butcher_order_certificate(method),
            "step_doubling": step_doubling_certificate(method),
        }
        for method, _ in TDG5_METHOD_ORDERS
    }
    all_methods = all(
        item["order"]["all_required_order_conditions_pass"]
        and item["step_doubling"][
            "contraction_is_stricter_than_two_to_minus_formal_order"
        ]
        for item in methods.values()
    )
    all_controls = bool(
        hermite["left_inverse_exact"]
        and hermite["right_inverse_exact"]
        and hermite["data_map_is_injective"]
        and hermite["condition_number_infinity_upper_bound"] == "54"
        and endpoint_attack["endpoint_only_sampling_misses_nonzero_interior"]
        and endpoint_attack["stage_complete_Hermite_data_detects_it"]
        and all_methods
    )
    return {
        "selected_sufficient_route": TDG5_SELECTED_SUFFICIENT_ROUTE,
        "owned_row_policy": TDG5_OWNED_ROW_POLICY,
        "required_runtime_layers": list(TDG5_REQUIRED_RUNTIME_LAYERS),
        "Hermite_data_map": hermite,
        "endpoint_only_alias_control": endpoint_attack,
        "method_controls": methods,
        "all_exact_controls_pass": all_controls,
        "actual_terminal_histories_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "state_advanced": False,
        "runtime_refinement_pair_implemented": False,
        "runtime_thresholds_frozen": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "GR0_case_eligible": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "physical_or_candidate_question_answered": False,
    }


__all__ = [
    "TDG5_HERMITE_DATA_ORDER",
    "TDG5_METHOD_ORDERS",
    "TDG5_OWNED_ROW_POLICY",
    "TDG5_REQUIRED_RUNTIME_LAYERS",
    "TDG5_SELECTED_SUFFICIENT_ROUTE",
    "butcher_order_certificate",
    "endpoint_only_alias_control",
    "evaluate_cubic",
    "hermite_coefficients",
    "hermite_data_map_certificate",
    "stage_complete_temporal_refinement_preflight",
    "step_doubling_certificate",
]
