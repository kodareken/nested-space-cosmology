from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.reference_connection import (
    REFERENCE_COORDINATE_ORDER,
    flat_spherical_annulus_connection,
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.reference_connection_second import (
    flat_spherical_annulus_connection_second_derivative,
)


class ReferenceConnectionSecondDerivativeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.reference = flat_spherical_annulus_reference(radial_domain_minimum=Q(1, 2))
        self.t, self.r, self.theta, self.phi = range(4)

    def test_complete_nonzero_second_derivative_table_at_two_radii(self) -> None:
        for radius in (Q(2), Q(5, 2)):
            data = flat_spherical_annulus_connection(
                self.reference, coordinate_radius=radius
            )
            second = flat_spherical_annulus_connection_second_derivative(
                self.reference, coordinate_radius=radius
            )
            expected = {
                (self.r, self.r, self.theta, self.r, self.theta): 2 / radius**3,
                (self.r, self.r, self.theta, self.theta, self.r): 2 / radius**3,
                (self.r, self.r, self.phi, self.r, self.phi): 2 / radius**3,
                (self.r, self.r, self.phi, self.phi, self.r): 2 / radius**3,
                (self.theta, self.theta, self.r, self.phi, self.phi): 2 * radius,
            }
            observed = {
                (first, second_direction, upper, lower_one, lower_two): second[first][second_direction][upper][lower_one][lower_two]
                for first in range(4)
                for second_direction in range(4)
                for upper in range(4)
                for lower_one in range(4)
                for lower_two in range(4)
                if second[first][second_direction][upper][lower_one][lower_two] != 0
            }
            self.assertEqual(observed, expected)
            self.assertEqual(data.coordinate_radius, radius)
            self.assertEqual(data.reference_id, "flat_spherical_annulus")

    def test_index_symmetries_and_mutation_sensitive_values(self) -> None:
        second = flat_spherical_annulus_connection_second_derivative(
            self.reference, coordinate_radius=Q(2)
        )
        self.assertEqual(REFERENCE_COORDINATE_ORDER, ("t", "r", "theta", "phi"))
        for first in range(4):
            for second_direction in range(4):
                for upper in range(4):
                    for lower_one in range(4):
                        for lower_two in range(4):
                            self.assertEqual(
                                second[first][second_direction][upper][lower_one][lower_two],
                                second[second_direction][first][upper][lower_one][lower_two],
                            )
                            self.assertEqual(
                                second[first][second_direction][upper][lower_one][lower_two],
                                second[first][second_direction][upper][lower_two][lower_one],
                            )
        self.assertEqual(second[self.r][self.r][self.theta][self.r][self.theta], Q(1, 4))
        self.assertEqual(second[self.r][self.r][self.phi][self.r][self.phi], Q(1, 4))
        self.assertEqual(second[self.theta][self.theta][self.r][self.phi][self.phi], Q(4))
        # These guard the common sign/power mistakes (+/- 2/r^3 and +/-2r).
        self.assertNotEqual(second[self.r][self.r][self.theta][self.r][self.theta], -Q(1, 4))
        self.assertNotEqual(second[self.theta][self.theta][self.r][self.phi][self.phi], -Q(4))
        self.assertNotEqual(second[self.r][self.r][self.theta][self.r][self.theta], Q(1, 2))

    def test_annulus_boundary_remains_failure_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "strictly inside"):
            flat_spherical_annulus_connection_second_derivative(
                self.reference, coordinate_radius=Q(1, 2)
            )


if __name__ == "__main__":
    unittest.main()
