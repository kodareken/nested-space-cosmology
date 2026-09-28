from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.exact_interval import (
    Interval,
    coerce_interval,
    interval,
    interval_matrix_infinity_row_sum_bound,
    interval_matrix_add,
    interval_matrix_determinant,
    interval_matrix_multiply,
    interval_matrix_point_multiply,
    interval_matrix_subtract,
    interval_matrix_vector,
    point_matrix_interval_multiply,
)
from recursive_horizons.fgc.interval_tangent import IntervalFirstTangent, primal_and_tangent


class ExactIntervalTests(unittest.TestCase):
    def test_exact_arithmetic_and_analytic_enclosure(self) -> None:
        x = interval(Q(-1, 3), Q(1, 2))
        self.assertEqual(x + 2, interval(Q(5, 3), Q(5, 2)))
        self.assertEqual(x - 2, interval(Q(-7, 3), Q(-3, 2)))
        self.assertEqual(x * interval(-2, 3), interval(-1, Q(3, 2)))
        self.assertEqual(x**2, interval(0, Q(1, 4)))
        self.assertEqual(interval(2, 4).reciprocal(), interval(Q(1, 4), Q(1, 2)))
        self.assertEqual(x.width(), Q(5, 6))
        self.assertEqual(x.midpoint(), Q(1, 12))
        self.assertEqual(x.radius(), Q(5, 12))
        self.assertEqual(x.abs_upper(), Q(1, 2))
        polynomial = x**2 + 2 * x - 1
        for value in (Q(-1, 3), Q(0), Q(1, 2)):
            self.assertTrue(polynomial.lower <= value**2 + 2 * value - 1 <= polynomial.upper)

    def test_explicit_sign_and_set_predicates(self) -> None:
        positive = interval(Q(1, 7), Q(3, 7))
        crossing = interval(-1, 1)
        self.assertTrue(positive.strictly_positive())
        self.assertTrue(positive.nonnegative())
        self.assertFalse(positive.contains_zero())
        self.assertTrue(crossing.contains_zero())
        self.assertFalse(crossing.strictly_positive())
        self.assertTrue(interval(1, 2).subset_of(interval(0, 2)))
        self.assertTrue(interval(1, 2).strictly_inside(interval(0, 3)))
        self.assertFalse(interval(1, 2).strictly_inside(interval(1, 3)))
        with self.assertRaisesRegex(TypeError, "truth value"):
            bool(positive)
        with self.assertRaisesRegex(TypeError, "ordering"):
            positive < interval(1, 2)

    def test_matrix_operations_and_infinity_bound(self) -> None:
        left = ((interval(1, 2), interval(-1, 1)),)
        right = ((interval(2, 3),), (interval(4, 5),))
        self.assertEqual(interval_matrix_multiply(left, right), ((interval(-3, 11),),))
        self.assertEqual(interval_matrix_vector(((interval(-1, 2), interval(3, 4)),), (interval(1, 2), 2)), (interval(4, 12),))
        self.assertEqual(
            interval_matrix_infinity_row_sum_bound(
                ((interval(-1, 2), interval(-3, 1)), (interval(0, 1), interval(-1, 1)))
            ),
            Q(5),
        )
        self.assertEqual(
            interval_matrix_add(((interval(1, 2),),), ((interval(-1, 1),),)),
            ((interval(0, 3),),),
        )
        self.assertEqual(
            interval_matrix_subtract(((interval(1, 2),),), ((interval(-1, 1),),)),
            ((interval(0, 3),),),
        )
        self.assertEqual(
            point_matrix_interval_multiply(((2, -1),), ((interval(1, 2),), (interval(3, 4),))),
            ((interval(-2, 1),),),
        )
        self.assertEqual(
            interval_matrix_point_multiply(((interval(1, 2), interval(3, 4)),), ((2,), (-1,))),
            ((interval(-2, 1),),),
        )
        self.assertEqual(
            interval_matrix_determinant(((interval(1, 2), interval(3, 4)), (interval(-1, 1), interval(2, 3)))),
            interval(-2, 10),
        )

    def test_interval_tangent_product_quotient_power_and_enclosure(self) -> None:
        x = IntervalFirstTangent.seed(interval(1, 2), 1)
        expression = (x**3 + 2 * x) / (x + 1)
        # f'(x) = (2x^3 + 3x^2 + 2) / (x + 1)^2.
        for value in (Q(1), Q(3, 2), Q(2)):
            derivative = (2 * value**3 + 3 * value**2 + 2) / (value + 1) ** 2
            self.assertTrue(expression.tangent.lower <= derivative <= expression.tangent.upper)
        self.assertEqual(primal_and_tangent(2), (interval(2), interval(0)))
        self.assertEqual(primal_and_tangent(x), (interval(1, 2), interval(1)))

    def test_fail_closed_invalid_endpoints_division_shapes_and_tangents(self) -> None:
        with self.assertRaisesRegex(TypeError, "exact rational"):
            Interval(0.0, 1)
        with self.assertRaisesRegex(TypeError, "exact rational"):
            interval(True)
        with self.assertRaisesRegex(ValueError, "must not exceed"):
            interval(2, 1)
        with self.assertRaisesRegex(ZeroDivisionError, "contains zero"):
            interval(1, 2) / interval(-1, 1)
        with self.assertRaisesRegex(ValueError, "dimensions"):
            interval_matrix_multiply(((interval(1), interval(2)),), ((interval(1), interval(2)),))
        with self.assertRaisesRegex(ValueError, "rectangular"):
            interval_matrix_vector(((interval(1),), (interval(1), interval(2))), (interval(1),))
        with self.assertRaisesRegex(ValueError, "square"):
            interval_matrix_determinant(((interval(1), interval(2)),))
        with self.assertRaisesRegex(TypeError, "interval or exact rational"):
            coerce_interval("1/2")
        with self.assertRaisesRegex(TypeError, "truth value"):
            bool(IntervalFirstTangent.constant(1))
        with self.assertRaisesRegex(TypeError, "ordering"):
            IntervalFirstTangent.constant(1) < IntervalFirstTangent.constant(2)


if __name__ == "__main__":
    unittest.main()
