"""Focused controls for the two independent TDG9 exact replay evaluators."""

from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg9_exact_temporal_arithmetic as primary,
)
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg9_exact_temporal_arithmetic_independent as independent,
)


Q = Fraction
Cubic = tuple[Fraction, Fraction, Fraction, Fraction]


def _evaluate(coefficients: Cubic, point: Fraction) -> Fraction:
    a0, a1, a2, a3 = coefficients
    return a0 + point * (a1 + point * (a2 + point * a3))


def _derivative(coefficients: Cubic, point: Fraction) -> Fraction:
    return coefficients[1] + point * (2 * coefficients[2] + point * 3 * coefficients[3])


def _segment(coefficients: Cubic, width: Fraction) -> tuple[float, ...]:
    values = (
        _evaluate(coefficients, Q(0)),
        _derivative(coefficients, Q(0)) / width,
        _evaluate(coefficients, Q(1)),
        _derivative(coefficients, Q(1)) / width,
        width,
    )
    answer = tuple(float(item) for item in values)
    if any(Q(*item.as_integer_ratio()) != exact for item, exact in zip(answer, values)):
        raise AssertionError("test fixture is not exactly binary64 representable")
    return answer


def _restrict(coefficients: Cubic, half: int) -> Cubic:
    a0, a1, a2, a3 = coefficients
    if half == 0:
        return a0, a1 / 2, a2 / 4, a3 / 8
    return (
        a0 + a1 / 2 + a2 / 4 + a3 / 8,
        a1 / 2 + a2 / 2 + 3 * a3 / 8,
        a2 / 4 + 3 * a3 / 8,
        a3 / 8,
    )


def _subtract(left: Cubic, right: Cubic) -> Cubic:
    return tuple(a - b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _scaled(coefficients: Cubic, factor: Fraction) -> Cubic:
    return tuple(factor * item for item in coefficients)  # type: ignore[return-value]


def _row(d01: Cubic, d12: Cubic) -> tuple[object, object, object]:
    zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
    medium_polynomials = tuple(_scaled(d01, Q(-1)) for _ in range(2))
    fine_polynomials = tuple(
        _subtract(
            _restrict(medium_polynomials[medium_index], half),
            d12,
        )
        for medium_index in range(2)
        for half in (0, 1)
    )
    return (
        _segment(zero, Q(1)),
        tuple(_segment(item, Q(1, 2)) for item in medium_polynomials),
        tuple(_segment(item, Q(1, 4)) for item in fine_polynomials),
    )


def _both(
    row: object, *, depth: int = 12, nodes: int = 10_000
) -> tuple[object, object]:
    rows = (row,)
    return (
        primary.assess_exact_temporal_refinement(
            rows, max_depth=depth, max_nodes=nodes
        ),
        independent.assess_exact_temporal_refinement_independently(
            rows, max_depth=depth, max_nodes=nodes
        ),
    )


class TDG9ExactTemporalArithmeticTests(unittest.TestCase):
    def test_independent_module_imports_no_primary_or_historical_private_arithmetic(
        self,
    ) -> None:
        path = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/tdg9_exact_temporal_arithmetic_independent.py"
        )
        tree = ast.parse(path.read_text())
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertFalse(
            imported
            & {
                "tdg9_exact_temporal_arithmetic",
                "tdg5_stage_complete_refinement_runtime",
                "tdg6_temporal_admission_runtime",
            }
        )
        self.assertFalse(
            any(
                isinstance(node, ast.ImportFrom)
                and node.module
                and (
                    "tdg5" in node.module
                    or "tdg6" in node.module
                    or node.module.endswith("tdg9_exact_temporal_arithmetic")
                )
                for node in ast.walk(tree)
            )
        )

    def test_exact_binary64_hermite_reconstruction_and_half_restriction(self) -> None:
        coefficients: Cubic = (Q(1), Q(2), Q(3), Q(4))
        segment = _segment(coefficients, Q(1, 2))
        self.assertEqual(primary.exact_hermite_coefficients(segment), coefficients)
        self.assertEqual(
            independent.independent_exact_hermite_coefficients(segment),
            coefficients,
        )
        for half in (0, 1):
            expected = _restrict(coefficients, half)
            first = primary.restrict_exact_cubic_to_half(coefficients, half=half)
            second = independent.independent_restrict_exact_cubic_to_half(
                coefficients, half=half
            )
            self.assertEqual(first, expected)
            self.assertEqual(second, expected)
            for point in (Q(0), Q(1, 4), Q(1, 2), Q(1)):
                self.assertEqual(
                    _evaluate(first, point),
                    _evaluate(coefficients, (Q(half) + point) / 2),
                )

    def test_exact_zero_stream_and_coefficient_hashes_agree(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        first, second = _both(_row(zero, zero))
        self.assertEqual(first.classification, "exact_zero")
        self.assertEqual(second.classification, "exact_zero")
        self.assertFalse(first.sufficient_pass_certified)
        self.assertFalse(second.sufficient_pass_certified)
        self.assertEqual(first.outer_difference.interval.lower, 0)
        self.assertEqual(first.finest_difference.interval.upper, 0)
        self.assertEqual(
            first.outer_difference.coefficient_stream_sha256,
            second.outer_difference.coefficient_stream_sha256,
        )
        self.assertEqual(
            first.finest_difference.coefficient_stream_sha256,
            second.finest_difference.coefficient_stream_sha256,
        )
        self.assertEqual(
            first.combined_coefficient_stream_sha256,
            second.combined_coefficient_stream_sha256,
        )

    def test_coefficient_stream_hash_has_frozen_oracle_and_detects_mutations(
        self,
    ) -> None:
        # These three digests were calculated from the literal ASCII stream
        # documented by the v1 hash domain, independently of either evaluator.
        finest: Cubic = (Q(0), Q(5, 8), Q(-7, 16), Q(9, 32))
        outer = _scaled(finest, Q(16))
        row = _row(outer, finest)
        first, second = _both(row, depth=14, nodes=100_000)
        expected_outer = (
            "e2af6bb6567a4c8338a23061a357e3c990331b6a967b5a29936e0aa1e0e2bec6"
        )
        expected_finest = (
            "2e880d300617262ade9fc189a32c089c80eac510d540ddbb2e3263bbd873fdd8"
        )
        expected_combined = (
            "b4e2e63ab6c895b5be5ee40eae9aea3596067e6aecc32a804b0601f05ccf2db7"
        )
        for evidence in (first, second):
            self.assertEqual(
                evidence.outer_difference.coefficient_stream_sha256,
                expected_outer,
            )
            self.assertEqual(
                evidence.finest_difference.coefficient_stream_sha256,
                expected_finest,
            )
            self.assertEqual(
                evidence.combined_coefficient_stream_sha256,
                expected_combined,
            )

        mutated_finest: Cubic = (Q(1, 32), Q(5, 8), Q(-7, 16), Q(9, 32))
        mutated_first, mutated_second = _both(
            _row(outer, mutated_finest), depth=14, nodes=100_000
        )
        for evidence in (mutated_first, mutated_second):
            self.assertNotEqual(
                evidence.finest_difference.coefficient_stream_sha256,
                expected_finest,
            )

        linear: Cubic = (Q(0), Q(1), Q(0), Q(0))
        other_row = _row(_scaled(linear, Q(16)), linear)
        for evaluator in (
            primary.assess_exact_temporal_refinement,
            independent.assess_exact_temporal_refinement_independently,
        ):
            forward = evaluator((row, other_row), max_depth=14, max_nodes=200_000)
            reverse = evaluator((other_row, row), max_depth=14, max_nodes=200_000)
            self.assertNotEqual(
                forward.outer_difference.coefficient_stream_sha256,
                reverse.outer_difference.coefficient_stream_sha256,
            )
            self.assertNotEqual(
                forward.finest_difference.coefficient_stream_sha256,
                reverse.finest_difference.coefficient_stream_sha256,
            )

    def test_p4_and_p1_replay_controls_resolve_in_both_evaluators(self) -> None:
        linear: Cubic = (Q(0), Q(1), Q(0), Q(0))
        p4_first, p4_second = _both(_row(_scaled(linear, Q(16)), linear))
        for evidence in (p4_first, p4_second):
            self.assertEqual(evidence.classification, "sufficient_contraction_pass")
            self.assertEqual(evidence.outer_difference.interval.lower, 16)
            self.assertEqual(evidence.outer_difference.interval.upper, 16)
            self.assertEqual(evidence.finest_difference.interval.lower, 1)
            self.assertEqual(evidence.finest_difference.interval.upper, 1)
            self.assertEqual(evidence.sufficient_pass_left, 8)
            self.assertEqual(evidence.sufficient_pass_right, 256)

        p1_first, p1_second = _both(_row(_scaled(linear, Q(2)), linear))
        for evidence in (p1_first, p1_second):
            self.assertEqual(evidence.classification, "sufficient_contraction_failure")
            self.assertFalse(evidence.sufficient_pass_certified)
            self.assertEqual(evidence.outer_difference.interval.upper, 2)
            self.assertEqual(evidence.finest_difference.interval.lower, 1)
            self.assertGreater(
                8 * evidence.finest_difference.interval.lower**2,
                evidence.outer_difference.interval.upper**2,
            )

    def test_exact_three_halves_squared_boundary_and_one_quantum_below(self) -> None:
        for classifier in (
            primary.three_halves_order_passes_exact,
            independent.independent_three_halves_order_passes_exact,
        ):
            self.assertTrue(
                classifier(outer_lower_squared=Q(8), finest_upper_squared=Q(1))
            )
            self.assertFalse(
                classifier(
                    outer_lower_squared=Q(8) - Q(1, 2**52),
                    finest_upper_squared=Q(1),
                )
            )

    def test_algebraic_nonzero_boundary_is_fail_closed_without_symbolic_proof(
        self,
    ) -> None:
        # max |4t^3-6t| on [0,1] is 2*sqrt(2), while max |1| is 1.
        # Thus the true squared comparison is exactly 8 == 8.  Rational
        # isolating intervals cannot certify equality at finite depth unless a
        # separate symbolic argument is supplied, so both bounded evaluators
        # must report resource exhaustion rather than round it into a pass.
        boundary: Cubic = (Q(0), Q(-6), Q(0), Q(4))
        constant: Cubic = (Q(1), Q(0), Q(0), Q(0))
        row = _row(boundary, constant)
        with self.assertRaises(
            primary.ExactTemporalArithmeticResourceExhausted
        ) as first:
            primary.assess_exact_temporal_refinement(
                (row,), max_depth=8, max_nodes=10_000
            )
        with self.assertRaises(
            independent.IndependentExactTemporalArithmeticResourceExhausted
        ) as second:
            independent.assess_exact_temporal_refinement_independently(
                (row,), max_depth=8, max_nodes=10_000
            )
        self.assertEqual(
            first.exception.evidence.reason,
            "max_depth_reached_before_contraction_resolved",
        )
        self.assertEqual(
            second.exception.evidence.reason,
            "max_depth_reached_before_contraction_resolved",
        )

    def test_floor_sized_monotone_and_nonmonotone_controls_are_not_zeroed(self) -> None:
        epsilon = Q(1, 2**900)
        monotone: Cubic = (Q(0), Q(1), Q(0), Q(0))
        nonmonotone: Cubic = (Q(0), Q(4), Q(-4), Q(0))
        for shape in (monotone, nonmonotone):
            first, second = _both(
                _row(_scaled(shape, 16 * epsilon), _scaled(shape, epsilon))
            )
            for evidence in (first, second):
                self.assertEqual(evidence.classification, "sufficient_contraction_pass")
                self.assertGreater(evidence.outer_difference.interval.lower, 0)
                self.assertGreater(evidence.finest_difference.interval.lower, 0)
                self.assertEqual(
                    evidence.outer_difference.interval.lower,
                    evidence.outer_difference.interval.upper,
                )
                self.assertEqual(
                    evidence.finest_difference.interval.lower,
                    evidence.finest_difference.interval.upper,
                )

    def test_above_floor_nonconvergent_control_fails_sufficiently(self) -> None:
        scale = Q(1, 2**20)
        shape: Cubic = (Q(0), Q(4), Q(-4), Q(0))
        first, second = _both(_row(_scaled(shape, 2 * scale), _scaled(shape, scale)))
        self.assertEqual(first.classification, "sufficient_contraction_failure")
        self.assertEqual(second.classification, "sufficient_contraction_failure")

    def test_two_irrational_derivative_roots_are_independently_isolated(self) -> None:
        # p'(t)=1-6t+6t^2 has two irrational roots in (0,1).  This exercises
        # both branches of the independent Sturm/bisection route while the
        # primary route sees only Bernstein control polygons.
        shape: Cubic = (Q(0), Q(1), Q(-3), Q(2))
        first, second = _both(
            _row(_scaled(shape, Q(16)), shape), depth=14, nodes=100_000
        )
        self.assertEqual(first.classification, "sufficient_contraction_pass")
        self.assertEqual(second.classification, "sufficient_contraction_pass")
        for first_interval, second_interval in (
            (first.outer_difference.interval, second.outer_difference.interval),
            (first.finest_difference.interval, second.finest_difference.interval),
        ):
            self.assertLessEqual(
                max(first_interval.lower, second_interval.lower),
                min(first_interval.upper, second_interval.upper),
            )
        self.assertEqual(
            first.combined_coefficient_stream_sha256,
            second.combined_coefficient_stream_sha256,
        )

    def test_resource_exhaustion_is_typed_and_never_a_pass(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        row = _row(zero, zero)
        with self.assertRaises(
            primary.ExactTemporalArithmeticResourceExhausted
        ) as first:
            primary.assess_exact_temporal_refinement((row,), max_depth=2, max_nodes=1)
        with self.assertRaises(
            independent.IndependentExactTemporalArithmeticResourceExhausted
        ) as second:
            independent.assess_exact_temporal_refinement_independently(
                (row,), max_depth=2, max_nodes=1
            )
        self.assertEqual(
            first.exception.evidence.reason,
            "max_nodes_exhausted_before_stream_completed",
        )
        self.assertEqual(
            second.exception.evidence.reason,
            "max_nodes_exhausted_before_stream_completed",
        )

    def test_user_runtime_error_cannot_impersonate_private_budget_exhaustion(
        self,
    ) -> None:
        def poisoned_rows() -> object:
            raise RuntimeError("node_budget_exhausted")
            yield  # pragma: no cover - makes this a generator without yielding

        for evaluator in (
            primary.assess_exact_temporal_refinement,
            independent.assess_exact_temporal_refinement_independently,
        ):
            with self.assertRaisesRegex(RuntimeError, "^node_budget_exhausted$"):
                evaluator(poisoned_rows(), max_depth=2, max_nodes=100)

    def test_completed_evidence_relations_and_digests_fail_closed(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        first, second = _both(_row(zero, zero))
        for evidence in (first, second):
            with self.assertRaisesRegex(TypeError, "interval must"):
                replace(evidence.outer_difference, interval=object())
            with self.assertRaisesRegex(ValueError, "lowercase SHA-256"):
                replace(
                    evidence.outer_difference,
                    coefficient_stream_sha256="A" * 64,
                )
            wrong_count = replace(
                evidence.outer_difference,
                polynomial_count=evidence.outer_difference.polynomial_count + 1,
            )
            with self.assertRaisesRegex(ValueError, "polynomial count"):
                replace(evidence, outer_difference=wrong_count)
            with self.assertRaisesRegex(ValueError, "twice row_count"):
                replace(evidence, row_count=evidence.row_count + 1)
            with self.assertRaisesRegex(ValueError, "total nodes"):
                replace(evidence, nodes_visited=evidence.nodes_visited + 1)
            with self.assertRaisesRegex(ValueError, "budget"):
                replace(evidence, max_nodes=evidence.nodes_visited - 1)
            excessive_depth = replace(
                evidence.outer_difference,
                maximum_depth_reached=evidence.max_depth + 1,
            )
            with self.assertRaisesRegex(ValueError, "depth exceeds"):
                replace(evidence, outer_difference=excessive_depth)
            with self.assertRaisesRegex(ValueError, "component hashes"):
                replace(evidence, combined_coefficient_stream_sha256="0" * 64)

    def test_resource_evidence_relations_and_digests_fail_closed(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        row = _row(zero, zero)
        failures: list[object] = []
        for evaluator, exception_type in (
            (
                primary.assess_exact_temporal_refinement,
                primary.ExactTemporalArithmeticResourceExhausted,
            ),
            (
                independent.assess_exact_temporal_refinement_independently,
                independent.IndependentExactTemporalArithmeticResourceExhausted,
            ),
        ):
            with self.assertRaises(exception_type) as caught:
                evaluator((row,), max_depth=2, max_nodes=1)
            failures.append(caught.exception.evidence)

        for evidence in failures:
            with self.assertRaisesRegex(ValueError, "unknown"):
                replace(evidence, reason="not_a_resource_reason")
            with self.assertRaisesRegex(ValueError, "consume the node budget"):
                replace(evidence, nodes_visited=0)
            with self.assertRaisesRegex(ValueError, "completed rows require"):
                replace(evidence, completed_rows=1)
            with self.assertRaisesRegex(ValueError, "lowercase SHA-256"):
                replace(evidence, outer_coefficient_prefix_sha256="F" * 64)

    def test_invalid_inputs_fail_before_any_scientific_classification(self) -> None:
        zero: Cubic = (Q(0), Q(0), Q(0), Q(0))
        row = _row(zero, zero)
        broken = (row[0], row[1], tuple(row[2][:-1]))
        for evaluator in (
            primary.assess_exact_temporal_refinement,
            independent.assess_exact_temporal_refinement_independently,
        ):
            with self.assertRaisesRegex(ValueError, "exactly 4"):
                evaluator((broken,), max_depth=2, max_nodes=100)
            with self.assertRaisesRegex(ValueError, "at least one"):
                evaluator((), max_depth=2, max_nodes=100)
            with self.assertRaisesRegex(ValueError, "positive"):
                evaluator((row,), max_depth=2, max_nodes=0)


if __name__ == "__main__":
    unittest.main()
