"""Bounded checks for the inward-march decision helper."""

import unittest

import numpy as np

from recursive_horizons import nsc_rn_inward as inward


class SlopeDecisionTests(unittest.TestCase):
    def test_regular_power_law_meets_the_bound(self):
        radius = np.linspace(0.8, 32.0, 401)
        pressure = 2.0 * radius ** (-1.0)
        decision = inward.inner_region_decision(radius, pressure, r_minus=0.75)
        self.assertTrue(decision["meets_hayward_on_evaluated_domain"])
        self.assertTrue(decision["decidable_on_this_domain"])
        self.assertAlmostEqual(decision["gap_in_cells"], (0.8 - 0.75) / (radius[1] - radius[0]))

    def test_steep_power_law_fails_the_bound(self):
        radius = np.linspace(0.8, 32.0, 401)
        pressure = 2.0 * radius ** (-7.0)
        decision = inward.inner_region_decision(radius, pressure, r_minus=0.75)
        self.assertFalse(decision["meets_hayward_on_evaluated_domain"])
        self.assertLess(decision["inner_region_slope_min"], inward.REGULAR_SLOPE_LIMIT)

    def test_slope_is_measured_not_extrapolated(self):
        radius = np.linspace(0.8, 2.0, 65)
        pressure = radius ** (-3.0)
        decision = inward.inner_region_decision(radius, pressure, r_minus=0.79)
        self.assertEqual(decision["extrapolation_below_inner_edge"], "none performed")
        self.assertGreater(decision["gap_to_horizon"], 0.0)

    def test_shape_mismatch_is_rejected(self):
        with self.assertRaises(ValueError):
            inward.inner_region_decision(np.linspace(1.0, 2.0, 6), np.ones(5), r_minus=0.9)


class LogSlopeTests(unittest.TestCase):
    def test_power_law_slope_is_exact(self):
        radius = np.linspace(1.0, 32.0, 301)
        for exponent in (-1.0, -3.0, -7.0):
            slope = inward._log_slope(radius, radius ** exponent)
            self.assertTrue(np.allclose(slope[1:-1], exponent, rtol=1e-8, atol=1e-8))


if __name__ == "__main__":
    unittest.main()