"""Focused probes of the normalization profiles; the runner reproduces the record."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import mpmath as mp
import numpy as np

from recursive_horizons.nsc_measure_normalization import (
    EULER, SMOOTH_KINDS, analytic_cylinder_values, cylinder_metric, dx_profiles,
    exact_hard_cylinder, fiber_operator, g_hard_log, g_heat_rank, g_proper_time,
    g_weighted_log, gauss_hard_cylinder, h_remainder, h_series, hard_cylinder_mode,
    heat_trace_expm, identity_pieces, identity_residual, independent_e1, integrate_profiles,
    independent_h_moments, independent_hard_cylinder_quadrature, lambda_profiles,
    log_ratio, smooth_metric,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from check_nsc_measure_normalization import compare_record  # noqa: E402


class ProfileIdentityTests(unittest.TestCase):
    def test_exact_identity_holds_and_hard_projector_differs(self):
        ratio = log_ratio(2., 1.)
        x = np.array([0.03, 0.2, 0.7, 1., 1.4, 4.])
        residual = identity_residual(x, ratio)
        self.assertLess(float(np.max(np.abs(residual))), 1e-14)
        pieces = identity_pieces(x, ratio)
        np.testing.assert_allclose(
            pieces["weighted_log"],
            g_proper_time(x) + g_heat_rank(x, ratio) + h_remainder(x),
            atol=1e-14,
        )
        self.assertGreater(float(np.max(np.abs(g_hard_log(x, ratio) - g_weighted_log(x, ratio)))), 0.05)

    def test_mpmath_identity_and_low_x_series(self):
        mp.mp.dps = 40
        ratio = mp.log(mp.mpf("1") / mp.mpf("2"))
        for raw in ("1e-8", "0.02", "0.5", "1", "3.25"):
            x = mp.mpf(raw)
            g_pt = mp.mpf("0.5") * mp.e1(x)
            heat = (ratio + mp.euler / 2) * mp.exp(-x)
            remainder = -mp.mpf("0.5") * (mp.exp(-x) * (mp.log(x) + mp.euler) + mp.e1(x))
            g_log = -mp.mpf("0.5") * mp.exp(-x) * (mp.log(x) - 2 * ratio)
            self.assertLess(abs(g_log - g_pt - heat - remainder), mp.mpf("1e-28"))
        xs = np.array([1e-6, 1e-4, 8e-3])
        np.testing.assert_allclose(h_remainder(xs), h_series(xs, 18), atol=1e-15)
        self.assertLess(abs(float(h_remainder(np.array([1e-12]))[0])), 2e-11)

    def test_independent_quadrature_of_e1_and_h_moments(self):
        value, error = independent_e1(1.2)
        self.assertLess(abs(value - float(__import__("scipy.special").special.exp1(1.2))), max(3 * error, 1e-12))
        moments = independent_h_moments()
        self.assertAlmostEqual(moments["integral_h"], -0.5, delta=1e-9)
        self.assertAlmostEqual(moments["integral_x_h"], -0.75, delta=1e-8)
        self.assertEqual(moments["a4_seeley_weight_h0"], 0.)

    def test_independent_profile_x_derivatives(self):
        ratio = log_ratio(3., 1.1)
        xs = np.array([0.18, 0.9, 2.2])
        analytic = dx_profiles(xs, ratio)
        step = 1e-6
        pieces_plus = identity_pieces(xs + step, ratio)
        pieces_minus = identity_pieces(xs - step, ratio)
        for kind in SMOOTH_KINDS:
            observed = (pieces_plus[kind] - pieces_minus[kind]) / (2 * step)
            np.testing.assert_allclose(observed, analytic[kind], atol=4e-10)


class CovariantMatrixTests(unittest.TestCase):
    def test_all_smooth_profiles_obey_uniform_clock_scaling(self):
        metric = smooth_metric(12, general=True)
        base = integrate_profiles(metric, 2., 1., 3, 24, gradients=False)
        shifted = integrate_profiles(replace(metric, lapse=1.3 * metric.lapse),
                                     2., 1., 3, 24, gradients=False)
        for kind in SMOOTH_KINDS:
            self.assertAlmostEqual(shifted["energy"][kind],
                                   1.3 * base["energy"][kind], delta=1e-12)

    def test_cylinder_fiber_matches_analytic_spectrum_and_expm_heat(self):
        metric = cylinder_metric(16)
        fiber = fiber_operator(metric, 0.7, 2, 2., 1., True)
        np.testing.assert_allclose(np.sort(fiber["values"]), np.sort(analytic_cylinder_values(metric, 0.7, 2)), atol=1e-13)
        self.assertAlmostEqual(
            heat_trace_expm(fiber["matrix"], 2.),
            float(np.sum(np.exp(-(fiber["values"] / 2.) ** 2))),
            delta=1e-12,
        )
        profiles, _ = lambda_profiles(fiber["values"], 2., 1.)
        self.assertAlmostEqual(
            float(np.sum(profiles["weighted_log"] - profiles["proper_time"] - profiles["heat_rank"] - profiles["remainder"])),
            0.,
            delta=1e-13,
        )
        self.assertGreater(abs(fiber["energies"]["remainder"]), 0.1)
        self.assertGreater(abs(fiber["energies"]["hard_log"] - fiber["energies"]["weighted_log"]), 0.1)

    def test_independent_differentiation_on_general_covariant_matrix(self):
        metric = smooth_metric(16, general=True)
        omega, kappa, cutoff, mass = 0.7, 2, 2., 1.
        fiber = fiber_operator(metric, omega, kappa, cutoff, mass, True)
        direction = metric.lapse * (0.3 * np.cos(2 * np.pi * metric.x / metric.length)
                                    + 0.2 * np.sin(4 * np.pi * metric.x / metric.length))
        step = 1e-4
        plus = replace(metric, lapse=metric.lapse + step * direction)
        minus = replace(metric, lapse=metric.lapse - step * direction)
        for kind in SMOOTH_KINDS:
            analytic = float(np.dot(fiber["gradients"][kind][0], direction))
            ep = fiber_operator(plus, omega, kappa, cutoff, mass, False)["energies"][kind]
            em = fiber_operator(minus, omega, kappa, cutoff, mass, False)["energies"][kind]
            self.assertAlmostEqual((ep - em) / (2 * step), analytic, delta=2e-9)

    def test_normalization_derivative_is_heat_for_weighted_log_not_rank(self):
        values = np.array([-2.2, -0.8, 0.5, 1.7])
        cutoff, mass = 2., 1.
        profiles, _ = lambda_profiles(values, cutoff, mass)
        step = 1e-5
        plus, _ = lambda_profiles(values, cutoff, mass * np.exp(step))
        minus, _ = lambda_profiles(values, cutoff, mass * np.exp(-step))
        heat = float(np.sum(np.exp(-(values / cutoff) ** 2)))
        for kind, predicted in (("proper_time", 0.), ("remainder", 0.),
                                ("heat_rank", heat), ("weighted_log", heat)):
            observed = float(np.sum(plus[kind] - minus[kind]) / (2 * step))
            self.assertAlmostEqual(observed, predicted, delta=2e-8)
        hard_plus = float(np.sum(plus["hard_log"] - minus["hard_log"]) / (2 * step))
        projected = float(np.sum((values / cutoff) ** 2 <= 1))
        self.assertAlmostEqual(hard_plus, projected, delta=2e-6)
        self.assertGreater(abs(hard_plus - heat), 0.2)

    def test_exact_hard_cylinder_matches_adaptive_split_and_wall(self):
        metric = cylinder_metric(16)
        exact = exact_hard_cylinder(metric, 2., 1., 8)
        quadrature = independent_hard_cylinder_quadrature(metric, 2., 1., 8)
        self.assertAlmostEqual(exact["energy"], -3.2978447242754125, delta=1e-15)
        self.assertAlmostEqual(exact["energy"], exact["interior"] + exact["wall"], delta=1e-15)
        self.assertAlmostEqual(quadrature["energy"], exact["energy"], delta=1e-12)
        self.assertAlmostEqual(quadrature["interior"], exact["interior"], delta=1e-12)
        self.assertAlmostEqual(quadrature["wall"], exact["wall"], delta=1e-12)
        self.assertEqual(len(exact["retained_modes"]), 2)
        reconstructed = 0.
        derivative = 0.
        for mode in exact["retained_modes"]:
            reconstructed += mode["interior"] + mode["wall"]
            derivative += 8 * mode["kappa"] / np.pi * mode["W"]
            self.assertAlmostEqual(mode["energy"], mode["interior"] + mode["wall"], delta=1e-15)
        self.assertAlmostEqual(reconstructed, exact["energy"], delta=1e-15)
        self.assertAlmostEqual(derivative, exact["normalization_log_derivative"], delta=1e-15)

    def test_uniform_lapse_hard_cylinder_is_interior_plus_wall(self):
        mode = hard_cylinder_mode(np.sqrt((np.pi / 4) ** 2 + 1.), 2., 1., 1)
        self.assertGreater(mode["W"], 0.)
        self.assertAlmostEqual(mode["energy"], mode["interior"] + mode["wall"], delta=1e-15)
        self.assertAlmostEqual(
            mode["normalization_log_derivative"],
            8 / np.pi * mode["W"],
            delta=1e-15,
        )
        closed = exact_hard_cylinder(cylinder_metric(16), 2., 1., 8)
        self.assertTrue(closed["uniform_lapse_derivative_equals_interior_plus_wall"])
        self.assertAlmostEqual(closed["interior_plus_wall"], closed["energy"], delta=1e-15)
        coarse = gauss_hard_cylinder(cylinder_metric(16), 2., 1., 8, 32)
        self.assertGreater(abs(coarse["energy"] - closed["energy"]), 0.4)
        self.assertAlmostEqual(coarse["energy"], -2.841515236464148, delta=1e-12)

    def test_invalid_scales_and_zero_modes_are_rejected(self):
        with self.assertRaises(ValueError):
            g_proper_time(0.)
        with self.assertRaises(ValueError):
            h_remainder(np.array([-0.1]))
        with self.assertRaises(ValueError):
            log_ratio(0., 1.)
        with self.assertRaises(ValueError):
            lambda_profiles(np.array([1., 0.]), 2., 1.)
        with self.assertRaises(ValueError):
            hard_cylinder_mode(0.5, 2., 1., 0)


class RecordContractTests(unittest.TestCase):
    def test_reproduction_rejects_changed_prescription_and_hashes(self):
        for left, right in (
            ({"identity_forced_into_Gamma_one": False}, {"identity_forced_into_Gamma_one": True}),
            ({"throat_root_fitted": False}, {"throat_root_fitted": True}),
            ({"hash": "abc"}, {"hash": "abd"}),
            ({"remainder": -14.7}, {"remainder": 0.0}),
            ({"a": 1}, {"a": 1, "b": 2}),
        ):
            with self.assertRaises(RuntimeError):
                compare_record(left, right)

    def test_cli_rejects_combined_modes_and_existing_output_without_overwriting(self):
        script = ROOT / "scripts/check_nsc_measure_normalization.py"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "existing.json"
            path.write_text("preserved bytes\n")
            for args in (("--check", "--output", str(path)), ("--output", str(path))):
                result = subprocess.run([sys.executable, str(script), *args], capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(path.read_text(), "preserved bytes\n")


if __name__ == "__main__":
    unittest.main()
