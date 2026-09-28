from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.exact_linear_algebra import (
    ZERO,
    add,
    adjugate_2x2,
    degree,
    determinant,
    evaluate,
    exact_division,
    inverse,
    matrix_multiply,
    matrix_det,
    matrix_inverse,
    matrix_mul,
    matrix_rank,
    multiply,
    polynomial,
    polynomial_matrix_determinant,
    polynomial_matrix_determinant_bareiss,
    polynomial_matrix_det,
    polynomial_matrix_det_bareiss,
    poly,
    poly_add,
    poly_div_exact,
    poly_eval,
    poly_mul,
    poly_scale,
    poly_sub,
    rank,
    scale,
    schur_numerator_2x2,
    subtract,
    transpose,
)


class ExactLinearAlgebraTests(unittest.TestCase):
    def test_polynomial_canonical_arithmetic_and_exact_division(self) -> None:
        self.assertEqual(polynomial([Q(1), Q(2), Q(0), Q(0)]), (Q(1), Q(2)))
        self.assertEqual(degree([Q(0)]), 0)
        self.assertEqual(add((Q(1),), (Q(-1),)), ZERO)
        self.assertEqual(subtract((Q(1), Q(3)), (Q(2),)), (Q(-1), Q(3)))
        self.assertEqual(multiply((Q(1), Q(1)), (Q(-1), Q(1))), (Q(-1), Q(0), Q(1)))
        self.assertEqual(evaluate((Q(-1), Q(0), Q(1)), Q(3)), Q(8))
        self.assertEqual(exact_division((Q(-1), Q(0), Q(1)), (Q(-1), Q(1))), (Q(1), Q(1)))
        with self.assertRaisesRegex(ValueError, "nonzero remainder"):
            exact_division((Q(1),), (Q(1), Q(1)))
        with self.assertRaisesRegex(ZeroDivisionError, "nonzero"):
            exact_division((Q(1),), ZERO)

    def test_exact_matrix_operations_and_singularity_rejection(self) -> None:
        matrix = ((Q(2), Q(1)), (Q(1), Q(1)))
        inverse_matrix = inverse(matrix)
        self.assertEqual(determinant(matrix), Q(1))
        self.assertEqual(rank(matrix), 2)
        self.assertEqual(matrix_multiply(matrix, inverse_matrix), ((Q(1), Q(0)), (Q(0), Q(1))))
        self.assertEqual(transpose(((Q(1), Q(2), Q(3)),)), ((Q(1),), (Q(2),), (Q(3),)))
        singular = ((Q(1), Q(2)), (Q(2), Q(4)))
        self.assertEqual(determinant(singular), Q(0))
        self.assertEqual(rank(singular), 1)
        with self.assertRaisesRegex(ValueError, "singular"):
            inverse(singular)
        with self.assertRaisesRegex(ValueError, "dimensions"):
            matrix_multiply(((Q(1), Q(2)),), ((Q(1), Q(2)),))

    def test_two_by_two_helpers_and_polynomial_matrix_determinant(self) -> None:
        matrix = ((Q(2), Q(3)), (Q(5), Q(7)))
        self.assertEqual(adjugate_2x2(matrix), ((Q(7), Q(-3)), (Q(-5), Q(2))))
        self.assertEqual(schur_numerator_2x2(matrix), Q(-1))
        pencil = (
            ((Q(1), Q(-1)), (Q(2),)),
            ((Q(3),), (Q(4), Q(1))),
        )
        self.assertEqual(polynomial_matrix_determinant(pencil), (Q(-2), Q(-3), Q(-1)))
        self.assertEqual(
            polynomial_matrix_determinant_bareiss(pencil),
            polynomial_matrix_determinant(pencil),
        )
        self.assertIs(poly, polynomial)
        self.assertIs(poly_add, add)
        self.assertIs(poly_sub, subtract)
        self.assertIs(poly_mul, multiply)
        self.assertIs(poly_scale, scale)
        self.assertIs(poly_eval, evaluate)
        self.assertIs(poly_div_exact, exact_division)
        self.assertIs(matrix_mul, matrix_multiply)
        self.assertIs(matrix_det, determinant)
        self.assertIs(matrix_rank, rank)
        self.assertIs(matrix_inverse, inverse)
        self.assertIs(polynomial_matrix_det, polynomial_matrix_determinant)
        self.assertIs(
            polynomial_matrix_det_bareiss,
            polynomial_matrix_determinant_bareiss,
        )
        with self.assertRaisesRegex(ValueError, "at most 4x4"):
            polynomial_matrix_determinant(tuple(tuple((Q(1),) for _ in range(5)) for _ in range(5)))

    def test_bareiss_polynomial_determinant_scales_beyond_four_by_four(self) -> None:
        diagonal = ((Q(1), Q(1)), (Q(2), Q(-1)))
        matrix = tuple(
            tuple(
                diagonal[row % 2] if row == column else ZERO
                for column in range(6)
            )
            for row in range(6)
        )
        expected = (Q(8), Q(12), Q(-6), Q(-11), Q(3), Q(3), Q(-1))
        self.assertEqual(polynomial_matrix_determinant_bareiss(matrix), expected)

        swapped = (
            (ZERO, (Q(1),), ZERO, ZERO, ZERO),
            ((Q(1), Q(1)), ZERO, ZERO, ZERO, ZERO),
            (ZERO, ZERO, (Q(2),), ZERO, ZERO),
            (ZERO, ZERO, ZERO, (Q(3),), ZERO),
            (ZERO, ZERO, ZERO, ZERO, (Q(5),)),
        )
        self.assertEqual(
            polynomial_matrix_determinant_bareiss(swapped),
            (Q(-30), Q(-30)),
        )

    def test_fail_closed_on_nonexact_values_and_bad_shapes(self) -> None:
        with self.assertRaisesRegex(TypeError, "exact rational"):
            polynomial([0.5])
        with self.assertRaisesRegex(TypeError, "exact rational"):
            polynomial([True])
        with self.assertRaisesRegex(ValueError, "rectangular"):
            determinant(((Q(1),), (Q(1), Q(2))))
        with self.assertRaisesRegex(ValueError, "square"):
            determinant(((Q(1), Q(2)),))
        with self.assertRaisesRegex(ValueError, "2x2"):
            adjugate_2x2(((Q(1),),))


if __name__ == "__main__":
    unittest.main()
