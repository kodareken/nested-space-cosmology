"""Focused tests for the seam-free radius, Dirac gap and Einstein residual.

All-field reproduction is `python3 scripts/check_nsc_smooth_geometry.py --check`.
These tests do not rebuild the published record and do not edit older NSC modules.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

from recursive_horizons.nsc_smooth_geometry import (
    areal_radius,
    areal_radius_derivative,
    areal_radius_second,
    areal_radius_third,
    bloch_gap,
    conservation_residual_grid,
    continuum_discriminant,
    continuum_first_band_edge,
    cosmological_constant_null_projection,
    derive_einstein_identities,
    eight_pi_G_stress,
    finite_difference_curvature_check,
    period_null_integrals,
    smoothness_invariants,
    smooth_motif,
    wave_number,
    zero_energy_integral,
)


ROOT = Path(__file__).resolve().parents[1]


class SmoothnessTests(unittest.TestCase):
    def test_throat_and_periodic_endpoint_invariants(self):
        for radius in (2.0, 4.0, 8.0):
            data = smoothness_invariants(1.0, radius)
            self.assertAlmostEqual(data["k_times_R"], np.pi / 2, places=15)
            self.assertAlmostEqual(data["r_min"], 1.0, places=15)
            self.assertAlmostEqual(data["r_second_at_origin"], 1.0, places=12)
            self.assertAlmostEqual(data["r_prime_plus_R"], 0.0, places=12)
            self.assertAlmostEqual(data["r_prime_minus_R"], 0.0, places=12)
            self.assertAlmostEqual(data["r_third_plus_R"], 0.0, places=10)
            self.assertAlmostEqual(data["r_plus_R"], data["r_minus_R"], places=15)
            self.assertAlmostEqual(data["smooth_r_prime_jump"], 0.0, places=15)
            self.assertGreater(data["seamed_sqrt_one_plus_rho_squared_r_prime_jump"], 1.0)

    def test_odd_even_parity_of_the_periodic_extension(self):
        radius = 4.0
        x = np.array([0.25, 1.1, 2.7])
        np.testing.assert_allclose(
            areal_radius(x, 1.0, radius), areal_radius(-x, 1.0, radius), atol=1e-15
        )
        np.testing.assert_allclose(
            areal_radius_derivative(x, 1.0, radius),
            -areal_radius_derivative(-x, 1.0, radius),
            atol=1e-15,
        )
        np.testing.assert_allclose(
            areal_radius_second(x, 1.0, radius),
            areal_radius_second(-x, 1.0, radius),
            atol=1e-14,
        )
        np.testing.assert_allclose(
            areal_radius_third(x, 1.0, radius),
            -areal_radius_third(-x, 1.0, radius),
            atol=1e-13,
        )
        self.assertAlmostEqual(
            float(areal_radius(radius, 1.0, radius)),
            float(areal_radius(-radius, 1.0, radius)),
            places=15,
        )


class CurvatureTests(unittest.TestCase):
    def test_symbolic_einstein_conservation_and_lapse_residuals_vanish(self):
        identities = derive_einstein_identities()
        for key in (
            "G_tt_minus_eight_pi_G_rho",
            "G_xx_minus_eight_pi_G_p_r",
            "G_theta_theta_minus_r_squared_times_eight_pi_G_p_t",
            "conservation_p_r_prime_plus_two_r_prime_over_r_times_p_r_minus_p_t",
            "local_null_identity_rho_plus_p_r_plus_two_r_second_over_r",
            "lapse_normal_projection_minus_eight_pi_G_rho",
            "lapse_euler_plus_r_squared_spatial_scalar",
        ):
            self.assertEqual(identities[key], "0")

    def test_metric_finite_difference_recovers_the_einstein_stress(self):
        check = finite_difference_curvature_check()
        self.assertLess(abs(check["rho_residual"]), 1e-8)
        self.assertLess(abs(check["p_r_residual"]), 1e-8)
        self.assertLess(abs(check["p_t_residual"]), 1e-8)

    def test_period_null_identity_and_throat_nec_violation(self):
        data = period_null_integrals(1.0, 2.0)
        self.assertLess(abs(data["identity_residual"]), 1e-12)
        self.assertLess(data["eight_pi_G_integral_rho_plus_p_r"], 0.0)
        self.assertAlmostEqual(data["throat_eight_pi_G_rho_plus_p_r"], -2.0, places=12)
        conservation = conservation_residual_grid(1.0, 2.0)
        self.assertLess(conservation["max_abs_exact_residual"], 1e-12)

    def test_cosmological_constant_has_vanishing_null_projection(self):
        data = cosmological_constant_null_projection()
        self.assertEqual(data["eight_pi_G_rho_plus_p"], 0.0)
        self.assertFalse(data["can_supply_required_radial_null_stress"])
        throat = eight_pi_G_stress(0.0, 1.0, 4.0)
        self.assertLess(float(np.min(np.asarray(throat["eight_pi_G_rho_plus_p_r"]))), 0.0)


class GapTests(unittest.TestCase):
    def test_zero_energy_monodromy_matches_elliptic_integral(self):
        radius = 2.0
        trace, det = continuum_discriminant(0.0, 1.0, radius)
        integral = zero_energy_integral(1.0, radius)
        self.assertAlmostEqual(trace, 2 * np.cosh(integral), delta=1e-9)
        self.assertAlmostEqual(det, 1.0, delta=1e-9)
        self.assertGreater(integral, 0.0)

    def test_continuum_and_bloch_gaps_agree_at_second_order(self):
        radius = 2.0
        continuum = continuum_first_band_edge(1.0, radius)
        coarse = bloch_gap(1.0, radius, 32)
        fine = bloch_gap(1.0, radius, 64)
        error_coarse = abs(coarse["minimum_sampled_bloch_gap"] - continuum["band_edge"])
        error_fine = abs(fine["minimum_sampled_bloch_gap"] - continuum["band_edge"])
        self.assertGreater(continuum["band_edge"], 0.5)
        self.assertAlmostEqual(coarse["minimum_phase"], 0.0, places=12)
        self.assertGreater(np.log2(error_coarse / error_fine), 1.8)
        self.assertLess(error_fine, 1e-5)

    def test_removing_angular_potential_closes_the_zero_phase_gap(self):
        free = smooth_motif(1.0, 2.0, 32, kappa=0.0)
        confined = smooth_motif(1.0, 2.0, 32, kappa=1.0)
        def gap(motif):
            return np.min(np.abs(np.linalg.eigvalsh(motif["H"] + motif["B"] + motif["B"].T)))
        self.assertLess(gap(free), 1e-12)
        self.assertGreater(gap(confined), 0.6)
        self.assertAlmostEqual(wave_number(2.0), np.pi / 4, places=15)


class PreservationTests(unittest.TestCase):
    def test_published_nsc3_source_bytes_are_unchanged(self):
        geometric = json.loads((ROOT / "results/nsc-3-geometric-chain.json").read_text())
        for path, digest in geometric["source_hashes"].items():
            observed = hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            self.assertEqual(observed, digest, msg=path)
        regulated = json.loads((ROOT / "results/nsc-3-regulated-recursion.json").read_text())
        for path, digest in regulated["source_hashes"].items():
            observed = hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            self.assertEqual(observed, digest, msg=path)
        radial = json.loads((ROOT / "results/nsc-3-radial-spectrum.json").read_text())
        observed = hashlib.sha256((ROOT / "scripts/check_nsc_radial_spectrum.py").read_bytes()).hexdigest()
        self.assertEqual(observed, radial["source_sha256"])


if __name__ == "__main__":
    unittest.main()
