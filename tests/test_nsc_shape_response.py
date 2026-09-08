from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np

from recursive_horizons.nsc_shape_response import (
    StaticAxialMetric, radial_reduction_identities, smooth_metric,
    spin_shape_difference, staggered_operator, stress_summary,
)
from recursive_horizons.nsc_covariant_measure import cylinder_stress

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from check_nsc_shape_response import compare_record  # noqa: E402


class ShapeResponseTests(unittest.TestCase):
    def test_radial_spin_connection_and_hilbert_measure_reduce_together(self):
        exact = radial_reduction_identities()
        self.assertEqual(exact["symmetric_derivative_residual"], "0")
        self.assertEqual(exact["measure_residual"], "0")

    def test_constant_staggered_spectrum_matches_independent_fourier_symbol(self):
        size, length, kappa = 32, 4., 2
        metric = StaticAxialMetric(length, np.ones(size), np.ones(size), np.ones(size))
        for eta in (0., .5):
            matrix, _ = staggered_operator(metric, kappa, eta)
            phase = 2 * np.pi * (np.arange(size) + eta) / size
            expected = np.sqrt(4 / metric.spacing**2 * np.sin(phase / 2)**2
                               + kappa**2 * np.cos(phase / 2)**2)
            np.testing.assert_allclose(np.sort(np.linalg.svd(matrix, compute_uv=False)),
                                       np.sort(expected), atol=1e-12)

    def test_local_conformal_change_preserves_the_same_operator(self):
        metric = smooth_metric(2., 64)
        factor = np.exp(.13 * np.cos(np.pi * metric.x / 2) + .02 * np.sin(np.pi * metric.x))
        changed = replace(metric, lapse=metric.lapse * factor,
                          radial_scale=metric.radial_scale * factor,
                          sphere_radius=metric.sphere_radius * factor)
        for eta in (0., .5):
            np.testing.assert_allclose(staggered_operator(metric, 3, eta)[0],
                                       staggered_operator(changed, 3, eta)[0], atol=1e-13)
        result = spin_shape_difference(changed, 4)
        self.assertLess(np.max(np.abs(result["local_weyl_identity"])), 1e-12)

    def test_all_metric_HF_gradients_match_local_energy_interventions(self):
        base = smooth_metric(2., 64)
        base = replace(base, lapse=1 + .08 * np.cos(np.pi * base.x / 2),
                       radial_scale=1 + .06 * np.sin(np.pi * base.x / 2))
        result = spin_shape_difference(base, 4)
        direction = np.eye(base.points)[base.points // 2]
        for field, gradient in (("lapse", "gradient_lapse"), ("radial_scale", "gradient_radial_scale"),
                                ("sphere_radius", "gradient_sphere_radius")):
            step = .001
            plus = replace(base, **{field: getattr(base, field) + step * direction})
            minus = replace(base, **{field: getattr(base, field) - step * direction})
            finite = (spin_shape_difference(plus, 4, False)["energy_difference"]
                      - spin_shape_difference(minus, 4, False)["energy_difference"]) / (2 * step)
            self.assertAlmostEqual(finite, np.dot(result[gradient], direction), delta=2e-7)

    def test_constant_stress_converges_to_full_angular_Bessel_response(self):
        length = 4.
        p = cylinder_stress(length, eta=0.)
        a = cylinder_stress(length, eta=.5)
        reference = p["radial_null_stress"] - a["radial_null_stress"]
        errors = []
        for size in (64, 128, 256):
            metric = StaticAxialMetric(length, np.ones(size), np.ones(size), np.ones(size))
            value = spin_shape_difference(metric, 8)
            errors.append(abs(np.mean(value["axial_null_difference"]) - reference))
        self.assertTrue(3.5 < errors[0] / errors[1] < 4.5)
        self.assertTrue(3.5 < errors[1] / errors[2] < 4.5)
        self.assertLess(errors[-1], 6e-7)

    def test_varying_radius_trace_and_spatial_conservation(self):
        residuals = []
        for size in (64, 128, 256):
            metric = smooth_metric(4., size)
            response = spin_shape_difference(metric, 6)
            summary = stress_summary(metric, response)
            self.assertLess(summary["maximum_absolute_local_trace"], 1e-12)
            self.assertLess(abs(summary["global_lapse_homogeneity_residual"]), 2e-7)
            self.assertGreater(summary["throat_axial_null_difference"], 0.)
            residuals.append(summary["rms_conservation_residual"])
        self.assertTrue(all(3.4 < a / b < 4.6 for a, b in zip(residuals, residuals[1:])))

    def test_angular_convergence_is_distinct_from_grid_convergence(self):
        metric = smooth_metric(2., 128)
        values = [stress_summary(metric, spin_shape_difference(metric, k))["throat_axial_null_difference"]
                  for k in (2, 6, 8)]
        self.assertGreater(abs(values[0] - values[-1]), 1e-4)
        self.assertLess(abs(values[1] - values[-1]), 2e-7)

    def test_units_and_lapse_homogeneity_of_stress(self):
        metric = smooth_metric(2., 64)
        original = spin_shape_difference(metric, 4)
        scale = 2.3
        changed = replace(metric, length=metric.length * scale,
                          sphere_radius=metric.sphere_radius * scale)
        result = spin_shape_difference(changed, 4)
        self.assertAlmostEqual(result["energy_difference"] * scale,
                               original["energy_difference"], delta=2e-8)
        np.testing.assert_allclose(result["axial_null_difference"] * scale**4,
                                   original["axial_null_difference"], atol=1e-8)

    def test_bad_metric_and_unresolved_angular_stencil_are_rejected(self):
        with self.assertRaises(ValueError):
            StaticAxialMetric(1., np.ones(8), np.zeros(8), np.ones(8))
        with self.assertRaises(ValueError):
            staggered_operator(smooth_metric(4., 16), 4, 0.)
        with self.assertRaises(ValueError):
            staggered_operator(smooth_metric(2., 64), 2, .2)

    def test_reproduction_checks_response_scope_and_every_key(self):
        with self.assertRaises(RuntimeError):
            compare_record({"absolute_AP_stress_derived": False}, {"absolute_AP_stress_derived": True})
        with self.assertRaises(RuntimeError):
            compare_record({"null": .0018}, {"null": -.0018})
        with self.assertRaises(RuntimeError):
            compare_record({"grid": 512}, {"grid": 513})
        with self.assertRaises(RuntimeError):
            compare_record({"gradient": [1., 2.]}, {"gradient": [1.]})


if __name__ == "__main__":
    unittest.main()
