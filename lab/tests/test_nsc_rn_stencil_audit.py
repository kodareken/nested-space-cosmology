"""Bounded checks for the nodal stencil audit.

The audit exists to answer one question: is the regularity verdict's
throat steepness physics, or an artifact of the stencils the chart
itself reports as failing? These tests pin that the recorded weights are
exact, and that the moment check misfires because it is ill-conditioned
rather than because the weights are wrong.
"""
import unittest

import numpy as np

from recursive_horizons import nsc_rn_stencil_audit as audit


class ExactWeightTests(unittest.TestCase):
    def _apply(self, weights, values):
        return float(np.dot(weights, np.asarray(values, dtype=float)[: weights.size]))

    def test_first_derivative_recovers_a_linear_ramp(self):
        nodes = np.arange(5, dtype=float)
        for at in (0.0, 1.0, 2.0, 3.0, 4.0):
            with self.subTest(at=at):
                weights = audit.fornberg_exact(5, 1, at)
                self.assertAlmostEqual(self._apply(weights, nodes), 1.0, places=12)

    def test_second_derivative_weights_match_the_known_five_point_rule(self):
        expected = np.array([-1.0 / 12.0, 4.0 / 3.0, -2.5, 4.0 / 3.0, -1.0 / 12.0])
        self.assertTrue(np.allclose(audit.fornberg_exact(5, 2, 2.0), expected, atol=1e-15))

    def test_first_derivative_weights_sum_to_zero(self):
        for width in (5, 7, 9):
            with self.subTest(width=width):
                self.assertAlmostEqual(float(audit.fornberg_exact(width, 1, 0.0).sum()), 0.0, places=12)

    def test_derivative_of_a_quadratic_is_exact(self):
        nodes = np.arange(9, dtype=float)
        squares = nodes ** 2
        for at in (0.0, 4.0, 8.0):
            with self.subTest(at=at):
                weights = audit.fornberg_exact(9, 1, at)
                self.assertAlmostEqual(self._apply(weights, squares), 2.0 * at, places=10)

    def test_second_derivative_of_a_cubic_is_exact(self):
        nodes = np.arange(9, dtype=float)
        cubes = nodes ** 3
        for at in (0.0, 4.0, 8.0):
            with self.subTest(at=at):
                weights = audit.fornberg_exact(9, 2, at)
                self.assertAlmostEqual(self._apply(weights, cubes), 6.0 * at, places=9)


class RecordedVersusExactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.audit = audit.audit_stencils()

    def test_recorded_weights_are_exact(self):
        """The whole point: the weights are not the defect."""
        worst = max(row["max_relative_weight_difference"] for row in self.audit["rows"])
        self.assertLess(worst, 1e-12, f"recorded weights differ from exact by {worst:.3e}")

    def test_configured_one_sided_widths_are_covered(self):
        widths = {row["width"] for row in self.audit["configured_one_sided"]}
        self.assertEqual(widths, {9, 10})

    def test_moment_check_misfires_even_on_exact_weights(self):
        """Ill-conditioned check, not wrong weights."""
        failing_exact = [row for row in self.audit["rows"] if not row["extended_passes"]]
        self.assertGreater(len(failing_exact), 0)
        for row in failing_exact[:5]:
            self.assertLess(row["max_relative_weight_difference"], 1e-12)

    def test_extended_precision_is_unavailable_here(self):
        """Why exact rationals rather than a wider accumulator."""
        self.assertEqual(np.longdouble(1.0).dtype, np.dtype(np.float64))


class SlopeValidityTests(unittest.TestCase):
    def _ramp(self, count=869):
        radius = np.linspace(1.04, 32.0, count)
        return radius, radius ** -3.0

    def test_throat_slope_survives_exact_weights(self):
        """The two weight sets must agree; the absolute slope is set by
        the short baseline between the end nodes, which both share."""
        radius, values = self._ramp()
        result = audit.slope_under_both_stencils(radius, values)
        self.assertLess(abs(result["slope_difference"]), 1e-9)
        self.assertLess(result["slope_plain"], 0.0)
        self.assertLess(abs(result["derivative_relative_difference"]), 1e-12)

    def test_short_grid_is_refused(self):
        radius, values = self._ramp(count=8)
        with self.assertRaises(ValueError):
            audit.slope_under_both_stencils(radius, values)

    def test_a_steep_profile_is_reported_as_steep_under_both(self):
        radius = np.linspace(1.04, 32.0, 869)
        values = radius ** -60.0
        result = audit.slope_under_both_stencils(radius, values)
        self.assertLess(result["slope_exact"], -1.0)
        self.assertLess(abs(result["slope_difference"]), 1e-6)


if __name__ == "__main__":
    unittest.main()