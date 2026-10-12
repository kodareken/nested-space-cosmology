from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.exact_interval import interval
from recursive_horizons.fgc.exact_interval_krawczyk import (
    MAXIMUM_DIMENSION,
    parametric_krawczyk_inclusion,
)


class ExactIntervalKrawczykTests(unittest.TestCase):
    def test_exact_point_linear_root(self) -> None:
        certificate = parametric_krawczyk_inclusion(
            center_inverse=((Q(1, 2),),),
            residual_at_center=(interval(0),),
            jacobian_box=((interval(2),),),
            displacement_box=(interval(-1, 1),),
        )
        self.assertEqual(certificate["rho_infinity_upper_bound"], 0)
        self.assertEqual(certificate["center_correction_box"], (interval(0),))
        self.assertEqual(certificate["krawczyk_displacement_image_box"], (interval(0),))
        self.assertEqual(certificate["uniform_inverse_infinity_norm_upper_bound"], Q(1, 2))
        self.assertTrue(certificate["krawczyk_image_strictly_inside_displacement_box"])

    def test_nontrivial_parametric_scalar_enclosure(self) -> None:
        # F(x,p)=2x-p with p in [-1/8,1/8] and x0=0.  The declared
        # displacement box [-1/4,1/4] contains the unique root p/2.
        certificate = parametric_krawczyk_inclusion(
            center_inverse=((Q(1, 2),),),
            residual_at_center=(interval(Q(-1, 8), Q(1, 8)),),
            jacobian_box=((interval(2),),),
            displacement_box=(interval(Q(-1, 4), Q(1, 4)),),
        )
        self.assertEqual(certificate["rho_infinity_upper_bound"], 0)
        self.assertEqual(certificate["center_correction_box"], (interval(Q(-1, 16), Q(1, 16)),))
        self.assertEqual(certificate["minimum_strict_componentwise_inclusion_margin"], Q(3, 16))
        self.assertIn("if the supplied residual", certificate["conditional_theorem"])
        self.assertFalse(certificate["input_enclosure_proven_here"])

    def test_nontrivial_two_by_two_contraction(self) -> None:
        certificate = parametric_krawczyk_inclusion(
            center_inverse=((Q(1, 2), 0), (0, Q(1, 3))),
            residual_at_center=(interval(Q(-1, 20), Q(1, 20)), interval(Q(-1, 30), Q(1, 30))),
            jacobian_box=((interval(Q(19, 10), Q(21, 10)), interval(Q(-1, 20), Q(1, 20))), (interval(Q(-1, 20), Q(1, 20)), interval(Q(29, 10), Q(31, 10)))),
            displacement_box=(interval(Q(-1, 4), Q(1, 4)), interval(Q(-1, 4), Q(1, 4))),
        )
        self.assertLess(certificate["rho_infinity_upper_bound"], 1)
        self.assertEqual(certificate["dimension"], 2)
        self.assertGreater(certificate["minimum_strict_componentwise_inclusion_margin"], 0)

    def test_failure_closed_types_shapes_contraction_and_inclusion(self) -> None:
        valid = {
            "center_inverse": ((Q(1),),),
            "residual_at_center": (interval(0),),
            "jacobian_box": ((interval(1),),),
            "displacement_box": (interval(-1, 1),),
        }
        with self.assertRaisesRegex(TypeError, "exact rational"):
            parametric_krawczyk_inclusion(**(valid | {"center_inverse": ((1.0,),)}))
        with self.assertRaisesRegex(TypeError, "exact rational"):
            parametric_krawczyk_inclusion(**(valid | {"residual_at_center": (interval(True),)}))
        with self.assertRaisesRegex(ValueError, "wrong dimension"):
            parametric_krawczyk_inclusion(**(valid | {"residual_at_center": (interval(0), interval(0))}))
        with self.assertRaisesRegex(ValueError, "strictly"):
            parametric_krawczyk_inclusion(**(valid | {"displacement_box": (interval(0, 1),)}))
        with self.assertRaisesRegex(ValueError, "below one"):
            parametric_krawczyk_inclusion(**(valid | {"jacobian_box": ((interval(3),),)}))
        with self.assertRaisesRegex(ValueError, "strictly inside"):
            parametric_krawczyk_inclusion(**(valid | {"residual_at_center": (interval(1, 2),)}))
        too_large = tuple(tuple(Q(int(row == column)) for column in range(MAXIMUM_DIMENSION + 1)) for row in range(MAXIMUM_DIMENSION + 1))
        with self.assertRaisesRegex(ValueError, "1.."):
            parametric_krawczyk_inclusion(
                center_inverse=too_large,
                residual_at_center=tuple(interval(0) for _ in too_large),
                jacobian_box=tuple(tuple(interval(item) for item in row) for row in too_large),
                displacement_box=tuple(interval(-1, 1) for _ in too_large),
            )


if __name__ == "__main__":
    unittest.main()
