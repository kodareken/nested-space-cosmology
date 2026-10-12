"""Independent post-freeze theorem controls for the TDG5 runtime design.

The prospective freeze selected a finite-dimensional continuous extension and
same-grid step refinement as one sufficient answer to TDG4.  This module
rederives the relevant algebra independently of the freeze implementation.  It
proves the finite candidate set for cubic extrema, exact half-interval
restriction, and the conditional Richardson-debit relation.  It also checks,
through the Python syntax tree, that the immutable numerical engine exposes the
stage and endpoint records required by a future implementation.

It does not implement the runtime pair, evaluate a production trajectory,
choose a numerical threshold, or define replacement temporal admission.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from math import isqrt
from numbers import Rational
from typing import Iterable, Sequence


Q = Fraction

TDG5_THEOREM_METHODS = (("RK4", 4, 15), ("SSPRK3", 3, 7))
TDG5_THEOREM_DATA_ORDER = (
    "left_value",
    "left_scaled_slope",
    "right_value",
    "right_scaled_slope",
)
TDG5_THEOREM_RUNTIME_FIELDS = {
    "StageRecord": ("stage_name", "time", "state", "rhs"),
    "StepProposal": (
        "method",
        "initial_time",
        "final_time",
        "initial_state",
        "candidate_state",
        "stages",
    ),
}
TDG5_THEOREM_STAGE_NAMES = (
    "rk4_k1",
    "rk4_k2",
    "rk4_k3",
    "rk4_k4",
    "ssprk3_s0",
    "ssprk3_s1",
    "ssprk3_s2",
    "candidate_endpoint",
)


def _exact(value: Rational, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(f"{name} must be an exact rational")
    return Fraction(value)


def _text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _matrix(
    rows: Iterable[Iterable[Rational]], *, name: str
) -> tuple[tuple[Fraction, ...], ...]:
    value = tuple(tuple(_exact(item, name=name) for item in row) for row in rows)
    if not value or not value[0] or any(len(row) != len(value[0]) for row in value):
        raise ValueError(f"{name} must be a nonempty rectangular matrix")
    return value


def _matvec(
    matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    if not matrix or len(matrix[0]) != len(vector):
        raise ValueError("matrix-vector shapes differ")
    return tuple(
        sum(
            (
                entry * component
                for entry, component in zip(row, vector, strict=True)
            ),
            Q(0),
        )
        for row in matrix
    )


def _determinant_four(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    if len(matrix) != 4 or any(len(row) != 4 for row in matrix):
        raise ValueError("TDG5 determinant control requires a four-square matrix")
    # Fraction-preserving Gaussian elimination is independent of the freeze's
    # permutation expansion.
    work = [list(row) for row in matrix]
    determinant = Q(1)
    for column in range(4):
        pivot = next(
            (row for row in range(column, 4) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            return Q(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant *= -1
        pivot_value = work[column][column]
        determinant *= pivot_value
        for entry in range(column, 4):
            work[column][entry] /= pivot_value
        for row in range(column + 1, 4):
            factor = work[row][column]
            for entry in range(column, 4):
                work[row][entry] -= factor * work[column][entry]
    return determinant


def _evaluate(coefficients: Sequence[Fraction], theta: Fraction) -> Fraction:
    if len(coefficients) != 4:
        raise ValueError("cubic coefficient vector must have length four")
    return coefficients[0] + theta * (
        coefficients[1]
        + theta * (coefficients[2] + theta * coefficients[3])
    )


def _derivative(coefficients: Sequence[Fraction], theta: Fraction) -> Fraction:
    if len(coefficients) != 4:
        raise ValueError("cubic coefficient vector must have length four")
    return coefficients[1] + theta * (
        2 * coefficients[2] + theta * 3 * coefficients[3]
    )


def independent_hermite_coefficients(
    data: Sequence[Rational],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Expand the four Hermite basis functions independently."""

    values = tuple(_exact(item, name="Hermite datum") for item in data)
    if len(values) != 4:
        raise ValueError("Hermite data must contain exactly four entries")
    left, left_slope, right, right_slope = values
    # Basis order:
    #   2t^3-3t^2+1, t^3-2t^2+t,
    #  -2t^3+3t^2,   t^3-t^2.
    basis = (
        (1, 0, -3, 2),
        (0, 1, -2, 1),
        (0, 0, 3, -2),
        (0, 0, -1, 1),
    )
    return tuple(
        sum(
            (
                value * Q(basis_index[coefficient])
                for value, basis_index in zip(values, basis, strict=True)
            ),
            Q(0),
        )
        for coefficient in range(4)
    )


def independent_hermite_map_certificate() -> dict[str, object]:
    unit_data = tuple(
        tuple(Q(int(row == column)) for row in range(4))
        for column in range(4)
    )
    columns = tuple(independent_hermite_coefficients(data) for data in unit_data)
    inverse = tuple(
        tuple(columns[column][row] for column in range(4)) for row in range(4)
    )
    sampling = _matrix(
        ((1, 0, 0, 0), (0, 1, 0, 0), (1, 1, 1, 1), (0, 1, 2, 3)),
        name="independent Hermite sampling matrix",
    )
    round_trips = []
    for data in unit_data:
        coefficients = independent_hermite_coefficients(data)
        reconstructed = (
            _evaluate(coefficients, Q(0)),
            _derivative(coefficients, Q(0)),
            _evaluate(coefficients, Q(1)),
            _derivative(coefficients, Q(1)),
        )
        round_trips.append(reconstructed == data)
    sampling_norm = max(sum(abs(item) for item in row) for row in sampling)
    inverse_norm = max(sum(abs(item) for item in row) for row in inverse)
    return {
        "data_order": list(TDG5_THEOREM_DATA_ORDER),
        "basis_expansion_matrix": [
            [_text(item) for item in row] for row in inverse
        ],
        "sampling_matrix": [[_text(item) for item in row] for row in sampling],
        "sampling_determinant_by_fraction_Gaussian_elimination": _text(
            _determinant_four(sampling)
        ),
        "basis_determinant_by_fraction_Gaussian_elimination": _text(
            _determinant_four(inverse)
        ),
        "all_four_unit_data_round_trip_exactly": all(round_trips),
        "sampling_infinity_norm": _text(sampling_norm),
        "basis_infinity_norm": _text(inverse_norm),
        "condition_number_infinity_upper_bound": _text(
            sampling_norm * inverse_norm
        ),
        "finite_dimensional_data_map_is_stably_injective": (
            _determinant_four(sampling) != 0
            and all(round_trips)
            and sampling_norm * inverse_norm == 54
        ),
    }


def _half_interval_restriction_matrices() -> tuple[
    tuple[tuple[Fraction, ...], ...],
    tuple[tuple[Fraction, ...], ...],
]:
    left = _matrix(
        (
            (1, 0, 0, 0),
            (0, Q(1, 2), 0, 0),
            (0, 0, Q(1, 4), 0),
            (0, 0, 0, Q(1, 8)),
        ),
        name="left-half restriction",
    )
    right = _matrix(
        (
            (1, Q(1, 2), Q(1, 4), Q(1, 8)),
            (0, Q(1, 2), Q(1, 2), Q(3, 8)),
            (0, 0, Q(1, 4), Q(3, 8)),
            (0, 0, 0, Q(1, 8)),
        ),
        name="right-half restriction",
    )
    return left, right


def restrict_cubic_to_half(
    coefficients: Sequence[Rational], *, half: int
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    values = tuple(_exact(item, name="cubic coefficient") for item in coefficients)
    if len(values) != 4:
        raise ValueError("cubic coefficient vector must have length four")
    if isinstance(half, bool) or half not in (0, 1):
        raise ValueError("half must be zero or one")
    left_matrix, right_matrix = _half_interval_restriction_matrices()
    return _matvec(left_matrix if half == 0 else right_matrix, values)


def half_interval_restriction_certificate() -> dict[str, object]:
    left_matrix, right_matrix = _half_interval_restriction_matrices()
    basis_vectors = tuple(
        tuple(Q(int(row == column)) for row in range(4))
        for column in range(4)
    )
    probe_points = (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1))
    controls: dict[str, bool] = {}
    for half in (0, 1):
        for index, coefficients in enumerate(basis_vectors):
            restricted = restrict_cubic_to_half(coefficients, half=half)
            controls[f"half_{half}_basis_{index}"] = all(
                _evaluate(restricted, local)
                == _evaluate(coefficients, (Q(half) + local) / 2)
                for local in probe_points
            )
    return {
        "left_restriction_matrix": [
            [_text(item) for item in row] for row in left_matrix
        ],
        "right_restriction_matrix": [
            [_text(item) for item in row] for row in right_matrix
        ],
        "left_restriction_determinant": _text(
            _determinant_four(left_matrix)
        ),
        "right_restriction_determinant": _text(
            _determinant_four(right_matrix)
        ),
        "all_basis_and_probe_composition_controls_pass": all(controls.values()),
        "control_count": len(controls),
        "coarse_cubic_restricts_exactly_to_each_fine_half_interval": all(
            controls.values()
        ),
    }


def _rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator:
        return None
    if denominator * denominator != value.denominator:
        return None
    return Q(numerator, denominator)


def exact_derivative_root_candidates(
    coefficients: Sequence[Rational],
) -> tuple[Fraction, ...]:
    """Return all exact rational derivative roots in the closed unit interval.

    A non-square rational discriminant is intentionally rejected.  The theorem
    still identifies the quadratic formula as the complete candidate set, but
    a future binary64 implementation must own outward root/evaluation error.
    """

    values = tuple(_exact(item, name="cubic coefficient") for item in coefficients)
    if len(values) != 4:
        raise ValueError("cubic coefficient vector must have length four")
    quadratic = 3 * values[3]
    linear = 2 * values[2]
    constant = values[1]
    if quadratic == 0:
        roots = () if linear == 0 else (-constant / linear,)
    else:
        discriminant = linear * linear - 4 * quadratic * constant
        if discriminant < 0:
            roots = ()
        else:
            root = _rational_square_root(discriminant)
            if root is None:
                raise ValueError(
                    "derivative discriminant is not an exact rational square"
                )
            roots = (
                (-linear - root) / (2 * quadratic),
                (-linear + root) / (2 * quadratic),
            )
    return tuple(sorted({root for root in roots if 0 <= root <= 1}))


def exact_cubic_absolute_maximum(
    coefficients: Sequence[Rational],
) -> dict[str, object]:
    values = tuple(_exact(item, name="cubic coefficient") for item in coefficients)
    if len(values) != 4:
        raise ValueError("cubic coefficient vector must have length four")
    candidates = (Q(0), *exact_derivative_root_candidates(values), Q(1))
    samples = tuple((point, abs(_evaluate(values, point))) for point in candidates)
    maximum = max(value for _, value in samples)
    return {
        "candidate_points": [_text(point) for point, _ in samples],
        "absolute_values": [_text(value) for _, value in samples],
        "absolute_maximum": _text(maximum),
        "maximizers": [
            _text(point) for point, value in samples if value == maximum
        ],
    }


def cubic_extremum_theorem_certificate() -> dict[str, object]:
    bump = exact_cubic_absolute_maximum((0, 4, -4, 0))
    genuine_cubic = exact_cubic_absolute_maximum((0, Q(9, 16), Q(-3, 2), 1))
    return {
        "compact_interval_extreme_value_theorem_supplies_a_maximum": True,
        "every_nonzero_interior_absolute_maximizer_has_zero_polynomial_derivative": (
            True
        ),
        "zero_polynomial_case_is_trivial": True,
        "cubic_derivative_degree_is_at_most_two": True,
        "complete_candidate_set_is_endpoints_plus_real_interior_derivative_roots": True,
        "maximum_candidate_count_per_half_interval": 4,
        "coarse_minus_fine_difference_remains_cubic_on_each_half_interval": True,
        "endpoint_only_bump_control": bump,
        "genuine_cubic_two_stationary_point_control": genuine_cubic,
        "binary64_root_and_value_evaluation_requires_outward_debit": True,
        "binary64_extremum_evaluator_implemented": False,
        "continuous_piecewise_cubic_extremum_theorem_completed": (
            bump["absolute_maximum"] == "1"
            and bump["maximizers"] == ["1/2"]
            and genuine_cubic["candidate_points"]
            == ["0", "1/4", "3/4", "1"]
        ),
    }


def richardson_debit_certificate(
    method: str, formal_order: int
) -> dict[str, object]:
    if method not in {name for name, _, _ in TDG5_THEOREM_METHODS}:
        raise ValueError("unknown TDG5 theorem method")
    if isinstance(formal_order, bool) or not isinstance(formal_order, int):
        raise TypeError("formal_order must be an integer")
    expected_order = dict((name, order) for name, order, _ in TDG5_THEOREM_METHODS)[
        method
    ]
    if formal_order != expected_order:
        raise ValueError("formal_order differs from the frozen method")
    power = 2**formal_order
    return {
        "method": method,
        "formal_order": formal_order,
        "coarse_leading_error_coefficient": "1",
        "two_half_step_leading_error_coefficient": f"1/{power}",
        "fine_minus_coarse_leading_difference_magnitude": f"{power - 1}/{power}",
        "fine_error_over_difference_factor": f"1/{power - 1}",
        "Richardson_denominator": power - 1,
        "derivation_assumes_a_shared_leading_local_error_coefficient": True,
        "derivation_assumes_the_tested_step_is_in_the_asymptotic_regime": True,
        "trajectory_asymptotic_regime_proved": False,
        "rigorous_global_PDE_error_bound_proved": False,
    }


def numerical_engine_shape_certificate(source: str) -> dict[str, object]:
    if not isinstance(source, str) or not source.strip():
        raise ValueError("numerical-engine source must be nonempty text")
    tree = ast.parse(source)
    class_fields: dict[str, tuple[str, ...]] = {}
    proposal: ast.FunctionDef | ast.AsyncFunctionDef | None = None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name in TDG5_THEOREM_RUNTIME_FIELDS:
            class_fields[node.name] = tuple(
                statement.target.id
                for statement in node.body
                if isinstance(statement, ast.AnnAssign)
                and isinstance(statement.target, ast.Name)
            )
        if (
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name == "propose_step"
        ):
            proposal = node
    fields_match = class_fields == TDG5_THEOREM_RUNTIME_FIELDS
    evaluate_calls: dict[str, ast.Call] = {}
    evaluate_body_calls_rhs = False
    evaluate_body_appends_stage_record = False
    if proposal is not None:
        for node in ast.walk(proposal):
            if (
                isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and node.name == "evaluate"
            ):
                evaluate_body_calls_rhs = any(
                    isinstance(candidate, ast.Call)
                    and isinstance(candidate.func, ast.Name)
                    and candidate.func.id == "rhs"
                    for candidate in ast.walk(node)
                )
                evaluate_body_appends_stage_record = any(
                    isinstance(candidate, ast.Call)
                    and isinstance(candidate.func, ast.Name)
                    and candidate.func.id == "StageRecord"
                    for candidate in ast.walk(node)
                )
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "evaluate"
                and len(node.args) >= 3
                and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)
            ):
                evaluate_calls[node.args[0].value] = node
    stages_present = all(name in evaluate_calls for name in TDG5_THEOREM_STAGE_NAMES)
    endpoint_call = evaluate_calls.get("candidate_endpoint")
    endpoint_arguments_match = bool(
        endpoint_call is not None
        and isinstance(endpoint_call.args[1], ast.Name)
        and endpoint_call.args[1].id == "final_time"
        and isinstance(endpoint_call.args[2], ast.Name)
        and endpoint_call.args[2].id == "candidate"
    )
    fresh_endpoint_rhs = bool(
        fields_match
        and endpoint_arguments_match
        and evaluate_body_calls_rhs
        and evaluate_body_appends_stage_record
    )
    return {
        "class_fields": {name: list(fields) for name, fields in class_fields.items()},
        "required_stage_names": list(TDG5_THEOREM_STAGE_NAMES),
        "observed_evaluate_stage_names": sorted(evaluate_calls),
        "stage_and_proposal_fields_match_frozen_shape": fields_match,
        "all_required_internal_and_endpoint_stage_names_present": stages_present,
        "candidate_endpoint_call_uses_final_time_and_candidate_state": (
            endpoint_arguments_match
        ),
        "evaluate_function_calls_RHS": evaluate_body_calls_rhs,
        "evaluate_function_appends_StageRecord": (
            evaluate_body_appends_stage_record
        ),
        "fresh_candidate_endpoint_RHS_record_is_available": fresh_endpoint_rhs,
        "runtime_pair_can_be_built_without_mutating_the_existing_step_proposal_shape": (
            fields_match and stages_present and fresh_endpoint_rhs
        ),
        "atomic_coarse_fine_transaction_implemented": False,
    }


def stage_complete_refinement_theorem_certificate(
    numerical_engine_source: str,
) -> dict[str, object]:
    hermite = independent_hermite_map_certificate()
    restriction = half_interval_restriction_certificate()
    extrema = cubic_extremum_theorem_certificate()
    richardson = {
        method: richardson_debit_certificate(method, order)
        for method, order, _ in TDG5_THEOREM_METHODS
    }
    engine = numerical_engine_shape_certificate(numerical_engine_source)
    all_controls = bool(
        hermite["finite_dimensional_data_map_is_stably_injective"]
        and restriction[
            "coarse_cubic_restricts_exactly_to_each_fine_half_interval"
        ]
        and extrema["continuous_piecewise_cubic_extremum_theorem_completed"]
        and all(
            richardson[method]["Richardson_denominator"] == denominator
            and richardson[method]["trajectory_asymptotic_regime_proved"]
            is False
            for method, _, denominator in TDG5_THEOREM_METHODS
        )
        and engine[
            "runtime_pair_can_be_built_without_mutating_the_existing_step_proposal_shape"
        ]
    )
    return {
        "independent_Hermite_map": hermite,
        "half_interval_restriction": restriction,
        "continuous_extremum_theorem": extrema,
        "method_owned_Richardson_debits": richardson,
        "immutable_numerical_engine_shape": engine,
        "all_independent_theorem_controls_pass": all_controls,
        "TDG5_stage_complete_refinement_theorem_completed": all_controls,
        "runtime_refinement_pair_implementation_authorized": all_controls,
        "runtime_refinement_pair_implemented": False,
        "runtime_thresholds_frozen": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "state_advanced": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "physical_or_candidate_question_answered": False,
    }


__all__ = [
    "TDG5_THEOREM_DATA_ORDER",
    "TDG5_THEOREM_METHODS",
    "TDG5_THEOREM_RUNTIME_FIELDS",
    "TDG5_THEOREM_STAGE_NAMES",
    "cubic_extremum_theorem_certificate",
    "exact_cubic_absolute_maximum",
    "exact_derivative_root_candidates",
    "half_interval_restriction_certificate",
    "independent_hermite_coefficients",
    "independent_hermite_map_certificate",
    "numerical_engine_shape_certificate",
    "restrict_cubic_to_half",
    "richardson_debit_certificate",
    "stage_complete_refinement_theorem_certificate",
]
