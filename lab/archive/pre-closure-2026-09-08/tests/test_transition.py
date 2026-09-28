from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.transition import (  # noqa: E402
    apparent_horizon_chi,
    classify_surface,
    closed_curvature_discriminator,
    curvature_invariants,
    friedmann_residual,
    future_null_expansions,
    hubble_derivative,
    hubble_parameter,
    scale_factor,
)


class ClosedDeSitterBackgroundTests(unittest.TestCase):
    def test_bounce_has_finite_minimum_and_positive_hubble_derivative(self) -> None:
        length = 2.5
        self.assertEqual(scale_factor(0.0, length), length)
        self.assertEqual(hubble_parameter(0.0, length), 0.0)
        self.assertAlmostEqual(hubble_derivative(0.0, length), 1.0 / length**2)
        self.assertGreater(scale_factor(-1.0, length), length)
        self.assertGreater(scale_factor(1.0, length), length)

    def test_constraint_and_invariants_are_regular(self) -> None:
        length = 3.0
        for time in (-20.0, -1.0, 0.0, 1.0, 20.0):
            with self.subTest(time=time):
                self.assertAlmostEqual(friedmann_residual(time, length), 0.0, places=14)
        invariants = curvature_invariants(length)
        self.assertAlmostEqual(invariants["ricci_scalar"], 12.0 / length**2)
        self.assertAlmostEqual(invariants["ricci_tensor_squared"], 36.0 / length**4)
        self.assertAlmostEqual(invariants["riemann_tensor_squared"], 24.0 / length**4)
        self.assertTrue(all(math.isfinite(value) and value > 0.0 for value in invariants.values()))


class NullExpansionTests(unittest.TestCase):
    def test_trapped_normal_and_antitrapped_examples(self) -> None:
        length = 1.0
        equator = math.pi / 2.0
        self.assertEqual(classify_surface(-2.0, equator, length), "trapped")
        self.assertEqual(classify_surface(0.0, math.pi / 4.0, length), "normal")
        self.assertEqual(classify_surface(2.0, equator, length), "anti_trapped")

    def test_apparent_horizons_have_a_zero_expansion(self) -> None:
        time = 1.0
        length = 1.0
        horizons = apparent_horizon_chi(time, length)
        plus, minus = future_null_expansions(time, horizons["theta_plus_zero"], length)
        self.assertAlmostEqual(plus, 0.0, places=13)
        self.assertGreater(minus, 0.0)
        plus, minus = future_null_expansions(time, horizons["theta_minus_zero"], length)
        self.assertGreater(plus, 0.0)
        self.assertAlmostEqual(minus, 0.0, places=13)
        bounce = apparent_horizon_chi(0.0, length)
        self.assertAlmostEqual(bounce["theta_plus_zero"], math.pi / 2.0)
        self.assertAlmostEqual(bounce["theta_minus_zero"], math.pi / 2.0)


class DiscriminatorAndValidationTests(unittest.TestCase):
    def test_closed_curvature_parameter_has_negative_sign(self) -> None:
        omega_k = closed_curvature_discriminator(10.0, 0.2)
        self.assertAlmostEqual(omega_k, -0.25)
        self.assertLess(omega_k, 0.0)
        self.assertAlmostEqual(closed_curvature_discriminator(10.0, -0.2), -0.25)

    def test_invalid_inputs_are_rejected(self) -> None:
        for bad in (0.0, -1.0, math.inf, math.nan):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    scale_factor(0.0, bad)
        for bad in (0.0, math.inf, math.nan):
            with self.subTest(hubble=bad):
                with self.assertRaises(ValueError):
                    closed_curvature_discriminator(1.0, bad)
        for chi in (0.0, math.pi, -0.1, math.inf):
            with self.subTest(chi=chi):
                with self.assertRaises(ValueError):
                    future_null_expansions(0.0, chi, 1.0)


if __name__ == "__main__":
    unittest.main()
