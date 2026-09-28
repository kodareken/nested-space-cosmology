from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.exact_interval import interval
from recursive_horizons.fgc.exact_interval_linear_algebra import (
    MAXIMUM_MATRIX_SIZE,
    center_preconditioned_neumann_inverse,
    identity,
    interval_matrix_infinity_row_sum_bound,
    interval_matrix_multiply,
    interval_matrix_subtract,
    point_matrix_interval_multiply,
)


class ExactIntervalLinearAlgebraTests(unittest.TestCase):
    def test_matrix_utilities_point_collapse_and_enclosure(self) -> None:
        left = ((interval(1, 2), interval(-1, 1)),)
        right = ((interval(2, 3),), (interval(4, 5),))
        self.assertEqual(interval_matrix_multiply(left, right), ((interval(-3, 11),),))
        self.assertEqual(identity(2), ((interval(1), interval(0)), (interval(0), interval(1))))
        self.assertEqual(
            point_matrix_interval_multiply(((2, -1),), ((interval(1, 2),), (interval(3, 4),))),
            ((interval(-2, 1),),),
        )
        self.assertEqual(
            interval_matrix_subtract(((interval(1, 2),),), ((interval(-1, 1),),)),
            ((interval(0, 3),),),
        )
        self.assertEqual(
            interval_matrix_infinity_row_sum_bound(
                ((interval(-1, 2), interval(-3, 1)), (interval(0, 1), interval(-1, 1)))
            ),
            Q(5),
        )

    def test_neumann_certificate_exact_identity_and_nontrivial_enclosure(self) -> None:
        exact = ((interval(2), interval(0)), (interval(0), interval(3)))
        exact_certificate = center_preconditioned_neumann_inverse(
            exact, ((Q(1, 2), 0), (0, Q(1, 3)))
        )
        self.assertEqual(exact_certificate["rho_infinity"], 0)
        self.assertEqual(exact_certificate["inverse_tail_entrywise_bound"], 0)
        self.assertEqual(
            exact_certificate["inverse_enclosure"],
            ((interval(Q(1, 2)), interval(0)), (interval(0), interval(Q(1, 3)))),
        )

        matrix = ((interval(Q(19, 10), Q(21, 10)), interval(0)), (interval(0), interval(Q(29, 10), Q(31, 10))))
        certificate = center_preconditioned_neumann_inverse(
            matrix, ((Q(1, 2), 0), (0, Q(1, 3)))
        )
        self.assertEqual(certificate["rho_infinity"], Q(1, 20))
        self.assertTrue(certificate["rho_strictly_below_one"])
        self.assertTrue(certificate["every_enclosed_matrix_invertible"])
        inverse = certificate["inverse_enclosure"]
        for value, reciprocal in ((Q(19, 10), Q(10, 19)), (Q(21, 10), Q(10, 21)), (Q(29, 10), Q(10, 29)), (Q(31, 10), Q(10, 31))):
            entry = inverse[0][0] if value.denominator == 10 and value < Q(5, 2) else inverse[1][1]
            self.assertLessEqual(entry.lower, reciprocal)
            self.assertGreaterEqual(entry.upper, reciprocal)

    def test_failure_closed_shapes_floats_and_noncontractive_mutation(self) -> None:
        with self.assertRaisesRegex(ValueError, "1.."):
            identity(0)
        with self.assertRaisesRegex(ValueError, "1.."):
            identity(MAXIMUM_MATRIX_SIZE + 1)
        with self.assertRaisesRegex(TypeError, "exact rational"):
            center_preconditioned_neumann_inverse(((interval(1.0),),), ((1,),))
        with self.assertRaisesRegex(ValueError, "rectangular"):
            interval_matrix_multiply(((interval(1),), (interval(1), interval(2))), ((interval(1),),))
        with self.assertRaisesRegex(ValueError, "dimensions"):
            center_preconditioned_neumann_inverse(((interval(1),),), ((1, 0), (0, 1)))
        with self.assertRaisesRegex(ValueError, "rho < 1"):
            center_preconditioned_neumann_inverse(((interval(3),),), ((1,),))


if __name__ == "__main__":
    unittest.main()
