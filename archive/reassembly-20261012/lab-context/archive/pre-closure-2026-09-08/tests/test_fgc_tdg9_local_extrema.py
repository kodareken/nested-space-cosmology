"""Pure controls for the two independent LOC1 exact localizers."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg9_local_extrema as primary  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_local_extrema_independent as independent  # noqa: E402


Q = Fraction


def _primary(coefficients: tuple[int, int, int, int], label: str = "x") -> primary.LocalCubic:
    return primary.LocalCubic(coefficients, {"label": label})


def _independent(coefficients: tuple[int, int, int, int], label: str = "x") -> independent.IndependentLocalCubic:
    return independent.IndependentLocalCubic(coefficients, {"label": label})


class TDG9LocalExtremaTests(unittest.TestCase):
    def _both(self, coefficients: tuple[int, int, int, int]) -> tuple[object, object]:
        first = primary.localize_absolute_maximum(
            [_primary(coefficients)], maximum_candidates=4, refinement_depth=80
        )
        second = independent.localize_absolute_maximum_independently(
            [_independent(coefficients)], maximum_candidates=4, refinement_bits=96
        )
        return first, second

    def test_unique_endpoint_agrees(self) -> None:
        first, second = self._both((0, 1, 0, 0))
        self.assertEqual(first.classification, "unique_maximum")
        self.assertEqual(second.classification, "unique_maximum")
        self.assertEqual(first.candidates[0].location, "right_endpoint")
        self.assertEqual(second.candidates[0].location, "right_endpoint")
        self.assertEqual(first.global_absolute_lower, Q(1))
        self.assertEqual(second.global_absolute_upper, Q(1))

    def test_unique_interior_stationary_agrees(self) -> None:
        first, second = self._both((0, 4, -4, 0))
        self.assertEqual(first.classification, "unique_maximum")
        self.assertEqual(second.classification, "unique_maximum")
        self.assertEqual(first.candidates[0].location, "interior_stationary")
        self.assertEqual(second.candidates[0].location, "interior_stationary")
        self.assertEqual(
            (first.candidates[0].parameter_lower, first.candidates[0].parameter_upper),
            (Q(1, 2), Q(1, 2)),
        )
        self.assertEqual(
            (second.candidates[0].parameter_lower, second.candidates[0].parameter_upper),
            (Q(1, 2), Q(1, 2)),
        )

    def test_irrational_stationary_root_isolated_by_both_routes(self) -> None:
        first, second = self._both((0, -1, 0, 1))
        self.assertEqual(first.classification, "unique_maximum")
        self.assertEqual(second.classification, "unique_maximum")
        self.assertEqual(first.candidates[0].location, "interior_stationary")
        self.assertEqual(second.candidates[0].location, "interior_stationary")
        self.assertLess(first.candidates[0].parameter_lower, first.candidates[0].parameter_upper)
        self.assertLess(second.candidates[0].parameter_lower, second.candidates[0].parameter_upper)
        self.assertLessEqual(
            max(first.candidates[0].parameter_lower, second.candidates[0].parameter_lower),
            min(first.candidates[0].parameter_upper, second.candidates[0].parameter_upper),
        )

    def test_nonunique_endpoint_result_is_preserved(self) -> None:
        first, second = self._both((1, 0, 0, 0))
        self.assertEqual(first.classification, "nonunique_or_interval_inconclusive")
        self.assertEqual(second.classification, "nonunique_or_interval_inconclusive")
        self.assertEqual({item.location for item in first.candidates}, {"left_endpoint", "right_endpoint"})
        self.assertEqual({item.location for item in second.candidates}, {"left_endpoint", "right_endpoint"})

    def test_two_interior_absolute_co_maxima_preserve_multiplicity(self) -> None:
        coefficients = (Q(1, 24), Q(5, 12), Q(-3, 2), Q(1))
        first = primary.localize_absolute_maximum(
            [primary.LocalCubic(coefficients, {"label": "two-roots"})],
            maximum_candidates=4,
            refinement_depth=96,
        )
        second = independent.localize_absolute_maximum_independently(
            [independent.IndependentLocalCubic(coefficients, {"label": "two-roots"})],
            maximum_candidates=4,
            refinement_bits=128,
        )
        self.assertEqual(first.classification, "nonunique_or_interval_inconclusive")
        self.assertEqual(second.classification, "nonunique_or_interval_inconclusive")
        self.assertEqual(len(first.candidates), 2)
        self.assertEqual(len(second.candidates), 2)
        self.assertEqual([item.location_ordinal for item in first.candidates], [0, 1])
        self.assertEqual([item.location_ordinal for item in second.candidates], [0, 1])
        for left, right in zip(first.candidates, second.candidates, strict=True):
            self.assertEqual(left.location, "interior_stationary")
            self.assertEqual(right.location, "interior_stationary")
            self.assertLessEqual(
                max(left.parameter_lower, right.parameter_lower),
                min(left.parameter_upper, right.parameter_upper),
            )

    def test_candidate_ceiling_fails_closed(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "candidate_ceiling_exhausted"):
            primary.localize_absolute_maximum([_primary((1, 0, 0, 0))], maximum_candidates=1)
        with self.assertRaisesRegex(RuntimeError, "candidate_ceiling_exhausted"):
            independent.localize_absolute_maximum_independently([_independent((1, 0, 0, 0))], maximum_candidates=1)

    def test_boolean_and_nonrational_inputs_are_rejected(self) -> None:
        with self.assertRaises(TypeError):
            primary.LocalCubic((True, 0, 0, 0), {})
        with self.assertRaises(TypeError):
            independent.IndependentLocalCubic((0.0, 0, 0, 0), {})

    def test_ulp_is_not_an_input_or_tolerance(self) -> None:
        self.assertNotIn("ulp", primary.localize_absolute_maximum.__code__.co_varnames)
        self.assertNotIn("ulp", independent.localize_absolute_maximum_independently.__code__.co_varnames)
        first, second = self._both((0, 1, 0, 0))
        self.assertFalse(first.tolerance_used)
        self.assertFalse(second.tolerance_used)


if __name__ == "__main__":
    unittest.main()
