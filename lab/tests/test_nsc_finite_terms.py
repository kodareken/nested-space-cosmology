from dataclasses import replace
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from recursive_horizons.nsc_finite_terms import (
    BASIS, FIELDS, StaticAxialMetric, direct_curvature_identities, exact_identities,
    exact_throat_sensitivities, finite_response, general_control_metric,
    independent_energy, periodic_derivative, response_matrix, smooth_metric, ward_controls,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from check_nsc_finite_terms import compare_record  # noqa: E402


class FiniteTermTests(unittest.TestCase):
    def test_independent_full_metric_curvature_and_exact_boundary_identities(self):
        self.assertEqual(set(direct_curvature_identities().values()), {"0"})
        identities = exact_identities()
        for name in ("weyl_contraction_residual", "euler_contraction_residual", "euler_primitive_residual", "pointwise_C2_local_weyl_residual"):
            self.assertEqual(identities[name], "0")
        self.assertEqual(set(identities["rigid_weyl_density_residuals"]), {"0"})

    def test_constant_cylinder_has_known_energy_and_accidental_rank_two(self):
        metric = StaticAxialMetric(4., np.ones(64), np.ones(64), np.ones(64))
        result = finite_response(metric)
        self.assertAlmostEqual(result["energy"][2], 64*np.pi/3, delta=1e-12)
        np.testing.assert_allclose(result["null"], 0, atol=1e-12)
        matrix = response_matrix(metric, result)
        self.assertEqual(matrix["rank"], 2)
        values = np.asarray(matrix["matrix"])
        np.testing.assert_allclose(values@np.array([0, 0, 1, -1/3, 0, 0]), 0, atol=1e-12)
        np.testing.assert_allclose(values@np.array([1, 1, 0, .25, 0, 0]), 0, atol=1e-12)

    def test_actual_smooth_profiles_resolve_four_bulk_coefficients_under_refinement(self):
        for radius in (2., 4.):
            matrices = []
            for points in (64, 128):
                metric = smooth_metric(radius, points)
                result = response_matrix(metric)
                self.assertEqual(result["rank"], 4)
                self.assertGreater(result["column_normalized_singular_values"][3], .5)
                matrices.append(result["matrix"])
            np.testing.assert_allclose(matrices[0], matrices[1], atol=2e-7)
            np.testing.assert_allclose(np.asarray(matrices[-1])[:, 4:], 0, atol=1e-10)

    def test_all_metric_variations_match_independent_scalar_energy_interventions(self):
        metric = general_control_metric(128)
        response = finite_response(metric)
        f = .4*np.cos(2*np.pi*metric.x/metric.length)+.2*np.sin(4*np.pi*metric.x/metric.length)
        for j, field in enumerate(FIELDS):
            direction = getattr(metric, field)*f
            derivative = metric.spacing*response["gradients"][:, j]@direction
            errors = []
            for step in (.001, .0005):
                plus = replace(metric, **{field: getattr(metric, field)+step*direction})
                minus = replace(metric, **{field: getattr(metric, field)-step*direction})
                observed = (independent_energy(plus)-independent_energy(minus))/(2*step)
                errors.append(np.linalg.norm(observed[:4]-derivative[:4]))
            self.assertLess(errors[-1], 1e-4)
            self.assertTrue(3.7 < errors[0]/errors[1] < 4.3)

    def test_local_weyl_and_radial_coordinate_identities_before_gauge_fixing(self):
        metric = general_control_metric(128)
        result = ward_controls(metric)
        self.assertLess(max(result["maximum_local_weyl_residual_by_basis"]), 2e-6)
        self.assertLess(max(result["maximum_local_radial_diffeomorphism_residual_by_basis"]), 1e-6)
        np.testing.assert_allclose(result["integrated_radial_gauge_response_by_basis"], 0., atol=1e-10)
        np.testing.assert_allclose(result["global_lapse_homogeneity_residual_by_basis"], 0., atol=1e-10)

    def test_each_local_metric_node_gradient_matches_finite_energy_change(self):
        metric = general_control_metric(64)
        response = finite_response(metric)
        direction = np.eye(metric.points)[metric.points//2]
        step = 1e-6
        for j, field in enumerate(FIELDS):
            plus = replace(metric, **{field: getattr(metric, field)+step*direction})
            minus = replace(metric, **{field: getattr(metric, field)-step*direction})
            observed = (independent_energy(plus)-independent_energy(minus))/(2*step)
            analytic = metric.spacing*response["gradients"][:, j]@direction
            np.testing.assert_allclose(observed, analytic, atol=5e-6, rtol=0)

    def test_exact_finite_coordinate_reparametrization_preserves_all_six_integrals(self):
        points, radius, epsilon = 128, 2., .09
        original = general_control_metric(points)
        x, length = original.x, original.length
        angle = 2*np.pi*x/length
        y = x+epsilon*np.sin(angle)
        jacobian = 1+epsilon*2*np.pi/length*np.cos(angle)
        angle_y = 2*np.pi*y/length
        wave = np.pi/(2*radius)
        changed = replace(original, lapse=np.exp(.08*np.cos(angle_y)+.03*np.sin(2*angle_y)),
                          radial_scale=np.exp(.06*np.sin(angle_y)-.02*np.cos(2*angle_y))*jacobian,
                          sphere_radius=np.sqrt(1+(np.sin(wave*y)/wave)**2))
        np.testing.assert_allclose(independent_energy(changed), independent_energy(original), atol=1e-10)

    def test_finite_local_weyl_keeps_C2_but_does_not_force_its_shape_response_to_zero(self):
        metric = general_control_metric(128)
        omega = np.exp(.11*np.cos(2*np.pi*metric.x/metric.length))
        changed = replace(metric, **{field: getattr(metric, field)*omega for field in FIELDS})
        original, conformal = finite_response(metric), finite_response(changed)
        np.testing.assert_allclose(original["densities"][2], conformal["densities"][2], atol=1e-10)
        self.assertGreater(np.linalg.norm(np.asarray(response_matrix(metric)["matrix"])[:, 2]), 1.)

    def test_throat_null_sensitivities_are_exact_and_not_the_cylinder_values(self):
        exact = exact_throat_sensitivities()
        self.assertEqual(exact["C2"]["null"], "32*(a**2*k**2 + 1)/(3*a**4)")
        self.assertEqual(exact["R2"]["null"], "-16*(4*a**2*k**2 + 1)/a**4")
        for radius in (2., 4.):
            result = finite_response(smooth_metric(radius, 128))
            k = np.pi/(2*radius)
            expected = [0, 4, 32*(1+k*k)/3, -16*(1+4*k*k), 0, 0]
            np.testing.assert_allclose(result["null"][:, 64], expected, atol=1e-7)

    def test_normalization_mass_and_spatial_units_are_kept_distinct(self):
        metric = general_control_metric(128)
        original = finite_response(metric)["energy"]
        factor = 1.7
        changed = replace(metric, length=metric.length*factor, sphere_radius=metric.sphere_radius*factor)
        np.testing.assert_allclose(finite_response(changed)["energy"], original*factor**np.array([3, 1, -1, -1, -1, -1]), atol=1e-10)
        np.testing.assert_allclose(finite_response(changed, 1/factor)["energy"], original/factor, atol=1e-10)
        with self.assertRaises(ValueError):
            finite_response(metric, 0)

    def test_spectral_derivative_adjoint_does_not_silently_assume_a_product_rule(self):
        rng = np.random.default_rng(421)
        a, b = rng.normal(size=(2, 64))
        da, db = (periodic_derivative(f, 4.) for f in (a, b))
        self.assertAlmostEqual(np.dot(a, db), -np.dot(da, b), delta=1e-11)
        np.testing.assert_allclose(periodic_derivative(da, 4.), periodic_derivative(a, 4., 2), atol=1e-11)

    def test_reproduction_checks_every_key_scope_null_vector_and_hash(self):
        for left, right in (({"scope": False}, {"scope": True}),
                            ({"rank": 4}, {"rank": 3}),
                            ({"hash": "abc"}, {"hash": "abd"}),
                            ({"matrix": [[1., 0.]]}, {"matrix": [[1., 1.]]}),
                            ({"a": 1}, {"a": 1, "b": 2})):
            with self.assertRaises(RuntimeError):
                compare_record(left, right)

    def test_cli_rejects_combined_modes_and_existing_output_without_overwriting(self):
        script = ROOT/"scripts/check_nsc_finite_terms.py"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"existing.json"
            path.write_text("preserved bytes\n")
            for args in (("--check", "--output", str(path)), ("--output", str(path))):
                result = subprocess.run([sys.executable, str(script), *args], capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(path.read_text(), "preserved bytes\n")


if __name__ == "__main__":
    unittest.main()
