from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.exact_interval import interval  # noqa: E402
from recursive_horizons.fgc.exact_interval_polynomial import (  # noqa: E402
    divide_by_monic_enclosure,
    monic_enclosure,
    polynomial_add,
    polynomial_derivative,
    polynomial_evaluate,
    polynomial_matrix_determinant,
    polynomial_multiply,
    quadratic_discriminant,
    quadratic_root_bracket,
)


class ExactIntervalPolynomialTests(unittest.TestCase):
    def test_point_polynomial_arithmetic_and_determinant_collapse_exactly(self) -> None:
        left = (interval(1), interval(2))
        right = (interval(-1), interval(1))
        self.assertEqual(
            polynomial_multiply(left, right),
            (interval(-1), interval(-1), interval(2)),
        )
        self.assertEqual(
            polynomial_add(left, right), (interval(0), interval(3))
        )
        self.assertEqual(polynomial_derivative((1, 2, 3)), (interval(2), interval(6)))
        self.assertEqual(polynomial_evaluate((1, 2, 3), 2), interval(17))
        determinant = polynomial_matrix_determinant(
            (((1, 1), (2,)), ((3,), (4, -1)))
        )
        self.assertEqual(determinant, (interval(-2), interval(3), interval(-1)))

    def test_monic_division_encloses_point_quotient(self) -> None:
        normalized = monic_enclosure((interval(2, 3), interval(4, 5)))
        self.assertTrue(normalized[-1].subset_of(interval(Q(4, 5), Q(5, 4))))
        divided = divide_by_monic_enclosure(
            (interval(-1), interval(0), interval(0), interval(0), interval(1)),
            (interval(-1), interval(0), interval(1)),
        )
        self.assertEqual(
            divided["quotient"], (interval(1), interval(0), interval(1))
        )
        self.assertEqual(divided["remainder"], (interval(0),))

    def test_uniform_quadratic_root_brackets(self) -> None:
        factor = (interval(Q(-101, 100), Q(-99, 100)), interval(0), interval(1))
        self.assertTrue(quadratic_discriminant(factor).strictly_positive())
        negative = quadratic_root_bracket(factor, interval(Q(-11, 10), Q(-9, 10)))
        positive = quadratic_root_bracket(factor, interval(Q(9, 10), Q(11, 10)))
        self.assertTrue(negative["exactly_one_real_root_for_every_enclosed_quadratic"])
        self.assertTrue(positive["exactly_one_real_root_for_every_enclosed_quadratic"])

    def test_fail_closed_on_ambiguous_leading_term_root_or_shape(self) -> None:
        with self.assertRaisesRegex(ValueError, "nonzero leading"):
            monic_enclosure((interval(1), interval(-1, 1)))
        with self.assertRaisesRegex(ValueError, "bracket"):
            quadratic_root_bracket((interval(-1), interval(0), interval(1)), interval(2, 3))
        with self.assertRaisesRegex(ValueError, "square"):
            polynomial_matrix_determinant((((1,),), ((1,), (2,))))


if __name__ == "__main__":
    unittest.main()
