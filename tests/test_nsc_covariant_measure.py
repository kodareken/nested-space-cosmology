import copy
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np
import mpmath as mp

from recursive_horizons.nsc_covariant_measure import (
    antiperiodic_fermi_integral, cylinder_energy, cylinder_stress,
    finite_cutoff_cylinder_stress, finite_weyl_average, product_geometry_identities,
    spatial_heat_weyl_ratio, spin_difference, spin_difference_heat_integral,
    winding_theta_sums,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from check_nsc_scale_closure import compare  # noqa: E402


class CovariantMeasureTests(unittest.TestCase):
    def test_nonzero_shape_response_has_zero_local_weyl_variation(self):
        identities = product_geometry_identities()
        self.assertEqual(identities["weyl_squared"], "4/(3*a**4)")
        self.assertEqual(identities["euler_density"], "0")
        self.assertEqual(identities["local_weyl_weight_residual"], "0")
        self.assertEqual(identities["rigid_four_dimensional_weyl_residual"], "0")
        self.assertEqual(identities["spacing_derivative"], "16*pi*alpha/(3*a**2)")
        self.assertEqual(identities["fixed_volume_spacing_derivative"], "128*pi**2*L*alpha/(3*V)")

    def test_flat_scale_average_retains_positive_measure_at_large_operator(self):
        value = finite_weyl_average((5., 10.))
        first, second = value["integrals"]
        slope = (second["negative_orbit_integral"] - first["negative_orbit_integral"]) / 5
        self.assertAlmostEqual(slope, np.exp(-np.euler_gamma), delta=1e-12)
        self.assertAlmostEqual(second["left_endpoint_integrand"], slope, delta=1e-12)
        with self.assertRaises(ValueError):
            finite_weyl_average((np.inf,))

    def test_spacing_ambiguous_local_term_cancels_in_radial_null_stress(self):
        identities = product_geometry_identities()
        self.assertEqual(identities["local_C_squared_energy_density"], "4*alpha/(3*a**4)")
        self.assertEqual(identities["local_C_squared_radial_pressure"], "-4*alpha/(3*a**4)")
        self.assertEqual(identities["local_C_squared_radial_null_projection"], "0")
        self.assertEqual(identities["local_C_squared_stress_trace"], "0")

    def test_full_angular_count_recovers_four_component_weyl_law(self):
        ratios = [spatial_heat_weyl_ratio(t) for t in (.1, .01, .001)]
        self.assertTrue(ratios[0] < ratios[1] < ratios[2] < 1.)
        self.assertAlmostEqual(ratios[-1], 1., delta=2e-4)
        # The O(t) coefficient also checks the spinor-sphere sum, beyond rank.
        self.assertAlmostEqual((1 - ratios[-1]) / .001, 1 / 6, delta=2e-5)

    def test_spin_energy_matches_independent_gaussian_heat_integral(self):
        for length, radius in ((1., 1.), (2.1, .8), (3.6, 1.4)):
            energy, derivative = spin_difference(length, radius)
            independent, error = spin_difference_heat_integral(length, radius)
            self.assertAlmostEqual(energy, independent, delta=max(2e-11, 3 * error))
            periodic = cylinder_energy(length, radius, eta=0.)
            antiperiodic = cylinder_energy(length, radius, eta=.5)
            np.testing.assert_allclose(np.subtract(periodic, antiperiodic),
                                       [energy, derivative], atol=1e-12)
            self.assertGreater(energy, 0.)
            self.assertLess(derivative, 0.)

    def test_pressure_variation_includes_length_prefactor_of_heat_integral(self):
        length, radius = 2.1, .8
        actual = spin_difference(length, radius)[1]
        step = .0005
        independent = (spin_difference_heat_integral(length + step, radius)[0]
                       - spin_difference_heat_integral(length - step, radius)[0]) / (2 * step)
        self.assertAlmostEqual(actual, independent, delta=3e-7)

    def test_fermi_integral_independently_gives_negative_AP_radial_null_stress(self):
        for length in (1., 2., 4.):
            stress = cylinder_stress(length, eta=.5)
            independent = antiperiodic_fermi_integral(length)
            for key in ("energy", "length_derivative", "radius_derivative"):
                self.assertAlmostEqual(stress[key], independent[key], delta=3e-11)
            self.assertLess(independent["energy"], 0.)
            self.assertGreater(independent["length_derivative"], 0.)
            self.assertLess(stress["radial_null_stress"], 0.)
            self.assertGreater(cylinder_stress(length, eta=0.)["radial_null_stress"], 0.)

    def test_transverse_pressure_comes_from_independent_radius_variation(self):
        length, radius, step = 2., 1., .0005
        stress = cylinder_stress(length, radius, eta=.5)
        radius_fd = (antiperiodic_fermi_integral(length, radius + step)["energy"]
                     - antiperiodic_fermi_integral(length, radius - step)["energy"]) / (2 * step)
        pressure = -radius_fd / (8 * np.pi * radius * length)
        self.assertAlmostEqual(stress["relative_transverse_pressure"], pressure, delta=3e-10)

    def test_compactification_stress_trace_and_null_dilation(self):
        for eta in (0., .5):
            original = cylinder_stress(2., 1.3, eta)
            scaled = cylinder_stress(5.4, 3.51, eta)
            self.assertAlmostEqual(original["relative_stress_trace"], 0., delta=1e-13)
            self.assertAlmostEqual(original["homogeneity_residual"], 0., delta=1e-12)
            self.assertAlmostEqual(scaled["radial_null_stress"] * 2.7**4,
                                   original["radial_null_stress"], delta=1e-12)

    def test_poisson_theta_derivative_resolves_small_negative_AP_tail(self):
        with mp.workdps(100):
            for b in (".02", ".1", "1", "10"):
                q = mp.exp(-mp.mpf(b))
                independent = q * mp.diff(lambda z: mp.jtheta(4, 0, z), q) / 2
                actual = winding_theta_sums(float(b), .5)[1]
                self.assertLess(actual, 0.)
                # Relative comparison is necessary: the first value is ~1e-49.
                self.assertLess(abs(mp.mpf(actual) / independent - 1), 1e-12)

    def test_finite_proper_time_cutoff_keeps_null_sign_and_recovers_continuum(self):
        for eta in (0., .5):
            continuum = cylinder_stress(1., eta=eta)["radial_null_stress"]
            magnitudes = []
            for cutoff in (1., 2., 4., 16.):
                finite = finite_cutoff_cylinder_stress(1., cutoff=cutoff, eta=eta)
                null = finite["axial_null_stress_integral"]
                self.assertGreater(null, 0.) if eta == 0 else self.assertLess(null, 0.)
                self.assertAlmostEqual(null, finite["axial_null_stress_from_energy_variation"], delta=1e-12)
                magnitudes.append(abs(null))
            self.assertTrue(all(a < b for a, b in zip(magnitudes, magnitudes[1:])))
            self.assertAlmostEqual(null, continuum, delta=1e-11)
        # The continuum E_L>0 argument does not apply to finite cutoff: the
        # exact theta/null identity must be used in this regime instead.
        self.assertLess(finite_cutoff_cylinder_stress(1., cutoff=1., eta=.5)["length_derivative_fixed_cutoff"], 0.)

    def test_null_cutoff_integral_matches_independent_energy_difference(self):
        for length, cutoff in ((1., 1.), (1., 4.), (2., 2.)):
            finite = finite_cutoff_cylinder_stress(length, cutoff=cutoff, eta=.5)
            step = .0001
            minus = finite_cutoff_cylinder_stress(length - step, cutoff=cutoff, eta=.5)["energy"]
            plus = finite_cutoff_cylinder_stress(length + step, cutoff=cutoff, eta=.5)["energy"]
            null = finite["relative_energy_density"] - (plus - minus) / (2 * step * 4 * np.pi)
            self.assertAlmostEqual(null, finite["axial_null_stress_integral"], delta=3e-8)

    def test_finite_cutoff_scale_identity_retains_cutoff_response(self):
        length, radius, cutoff = 2., .7, 2.
        finite = finite_cutoff_cylinder_stress(length, radius, cutoff, eta=.5)
        self.assertAlmostEqual(finite["cutoff_homogeneity_residual"], 0., delta=1e-12)
        self.assertAlmostEqual(finite["relative_stress_trace"],
                               -cutoff * finite["cutoff_derivative"] / (4 * np.pi * radius**2 * length),
                               delta=1e-12)
        self.assertGreater(abs(finite["relative_stress_trace"]), 1e-5)
        scale = 2.3
        changed = finite_cutoff_cylinder_stress(length * scale, radius * scale, cutoff / scale, eta=.5)
        self.assertAlmostEqual(changed["axial_null_stress_integral"] * scale**4,
                               finite["axial_null_stress_integral"], delta=1e-12)
        with self.assertRaises(ValueError):
            finite_cutoff_cylinder_stress(1., cutoff=0.)

    def test_angular_and_winding_truncations_converge_independently(self):
        reference = np.array(spin_difference(1.))
        for dimension in ("angular_max", "winding_max"):
            coarse = np.array(spin_difference(1., **{dimension: 7}))
            fine = np.array(spin_difference(1., **{dimension: 31}))
            self.assertGreater(np.linalg.norm(coarse - reference), 1e-6)
            self.assertLess(np.linalg.norm(fine - reference), 2e-9)

    def test_physical_dilation_has_energy_and_pressure_dimensions(self):
        energy, derivative = spin_difference(2., radius=1.3)
        scale = 2.7
        scaled = spin_difference(2. * scale, radius=1.3 * scale)
        np.testing.assert_allclose(scaled, [energy / scale, derivative / scale**2], atol=1e-12)
        pressure = -derivative / (4 * np.pi * 1.3**2)
        scaled_pressure = -scaled[1] / (4 * np.pi * (1.3 * scale)**2)
        self.assertAlmostEqual(scaled_pressure * scale**4, pressure, delta=1e-12)

    def test_domain_rejects_unstated_twist_or_invalid_geometry(self):
        for kwargs in ({"length": 0.}, {"length": 1., "radius": -1.},
                       {"length": 1., "eta": .2}, {"length": 1., "angular_max": 1.5},
                       {"length": 1., "winding_max": True}):
            with self.assertRaises(ValueError):
                cylinder_energy(**kwargs)

    def test_reproduction_compares_nested_numbers_and_scope_booleans(self):
        sample = {"controls": [{"energy": spin_difference(2.)[0]}],
                  "stationary_geometry_derived": False}
        changed = copy.deepcopy(sample)
        changed["controls"][0]["energy"] += .001
        with self.assertRaises(RuntimeError):
            compare(sample, changed)
        changed = copy.deepcopy(sample)
        changed["stationary_geometry_derived"] = True
        with self.assertRaises(RuntimeError):
            compare(sample, changed)

    def test_check_and_output_are_mutually_exclusive(self):
        result = subprocess.run([sys.executable, str(ROOT / "scripts/check_nsc_covariant_measure.py"),
                                 "--check", "--output", "unused.json"],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("not allowed with argument", result.stderr)


if __name__ == "__main__":
    unittest.main()
