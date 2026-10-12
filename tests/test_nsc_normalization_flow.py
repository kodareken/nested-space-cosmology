"""Focused probes of the free-Dirac normalization flow; the runner reproduces the record."""
from __future__ import annotations

from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import sympy as sp

from recursive_horizons.nsc_covariant_operator import cylinder_metric, smooth_metric
from recursive_horizons.nsc_finite_terms import BASIS, FIELDS
from recursive_horizons.nsc_normalization_flow import (
    ABS_TOL, ARBITRARY_CONTROL_COEFFICIENTS, FINITE_DIFFERENCE_ATOL,
    SYMBOLIC_COEFFICIENTS, a4_basis_weights, beta_mref_exact, build_record,
    c2_only_coefficient_shift, centered_normalization_derivative, compare_record,
    discrete_euler_refinement, exact_flow_identities, finite_coefficient_shift,
    finite_increment, float_tolerances, matched_cancellation,
    remainder_log_derivatives, shift_coefficients,
)

ROOT = Path(__file__).resolve().parents[1]


class ExactFlowTests(unittest.TestCase):
    def test_exact_identities_and_a4_reduction(self):
        identities = exact_flow_identities()
        for key in (
            "dE_sub_dlnM_residual", "a4_R2_derivative", "a4_weyl_euler_identity_residual",
            "compensated_energy_residual", "arbitrary_initial_coefficients_cancel",
            "euler_lagrange_compensation_residual", "M_basis_classical_M4_residual",
            "M_basis_classical_M2R_residual", "continuum_bulk_energy_residual_on_closed_cell",
            "continuum_bulk_EL_residual_with_vanishing_Euler_box",
        ):
            self.assertEqual(identities[key], "0")
        self.assertEqual(set(identities["composition_residuals"]), {"0"})
        self.assertEqual(set(identities["inverse_residuals"]), {"0"})
        self.assertEqual(identities["dE_sub_dlnM"], "A4/(16*pi**2)")
        self.assertEqual(identities["a4_density"], "(-18*C2+11*E4-12*BoxR)/360")
        self.assertFalse(identities["a4_contains_R2"])
        self.assertTrue(identities["M_ratio_is_not_recursive_Omega"])
        self.assertTrue(identities["continuum_closed_cell_bulk_running_is_C2"])
        self.assertTrue(identities["continuum_Euler_and_box_bulk_EL_vanish"])
        self.assertTrue(identities["discrete_Euler_column_is_product_rule_error"])
        self.assertEqual(identities["unspecified_initial_coefficients"], "None")
        self.assertEqual(identities["symbolic_initial_coefficients"], list(SYMBOLIC_COEFFICIENTS))
        self.assertEqual(identities["marginal_stationary_coefficient"], "c_R2")
        weights = a4_basis_weights()
        self.assertEqual(str(weights[2]), "-1/20")
        self.assertEqual(str(weights[3]), "0")
        self.assertEqual(str(weights[4]), "11/360")
        self.assertEqual(str(weights[5]), "-1/30")

    def test_mref_betas_and_symbolic_initial_coefficients(self):
        betas = beta_mref_exact()
        self.assertEqual(str(betas[0]), "0")
        self.assertEqual(str(betas[1]), "0")
        self.assertEqual(str(sp.simplify(betas[2] - 1 / (320 * sp.pi**2))), "0")
        self.assertEqual(str(betas[3]), "0")
        self.assertEqual(str(sp.simplify(betas[4] + 11 / (5760 * sp.pi**2))), "0")
        self.assertEqual(str(sp.simplify(betas[5] - 1 / (480 * sp.pi**2))), "0")
        scale = sp.symbols("s", positive=True)
        functionals = sp.symbols("F0:6")
        coefficients = sp.symbols("c0:6")
        jump = sum(
            (coeff + beta * sp.log(scale)) * term - coeff * term
            for coeff, beta, term in zip(coefficients, betas, functionals)
        )
        a4 = sum(weight * term for weight, term in zip(a4_basis_weights(), functionals))
        self.assertEqual(sp.simplify(jump + a4 * sp.log(scale) / (16 * sp.pi**2)), 0)

    def test_m_basis_classical_betas_keep_arbitrary_coefficients(self):
        mass, reference, c_m4, c_m2r = sp.symbols("M M_ref c_M4 c_M2R", positive=True)
        hat_m4 = c_m4 * (reference / mass)**4
        hat_m2r = c_m2r * (reference / mass)**2
        self.assertEqual(sp.simplify(mass * sp.diff(hat_m4, mass) + 4 * hat_m4), 0)
        self.assertEqual(sp.simplify(mass * sp.diff(hat_m2r, mass) + 2 * hat_m2r), 0)

    def test_no_newton_cosmological_or_recursive_fixed_point_inference(self):
        identities = exact_flow_identities()
        self.assertTrue(identities["no_physical_Newton_running_inferred"])
        self.assertTrue(identities["no_physical_cosmological_running_inferred"])
        self.assertTrue(identities["no_recursive_fixed_point_inferred"])
        self.assertTrue(identities["M_ratio_is_not_recursive_Omega"])


class CancellationTests(unittest.TestCase):
    def test_cylinder_analytic_derivative(self):
        metric = cylinder_metric(32)
        result = remainder_log_derivatives(metric)
        self.assertAlmostEqual(result["energy"], -1 / (15 * np.pi), delta=1e-14)
        self.assertAlmostEqual(result["A4"], -16 * np.pi / 15, delta=1e-14)

    def test_energy_and_every_metric_gradient_cancel(self):
        metrics = (cylinder_metric(32), smooth_metric(32), smooth_metric(32, general=True))
        scales = (0.5, 2., float(np.e), 10., 1 / np.pi)
        initials = (None, ARBITRARY_CONTROL_COEFFICIENTS)
        for metric in metrics:
            for scale in scales:
                for initial in initials:
                    row = matched_cancellation(metric, 2., 1., scale, 1., initial)
                    self.assertLess(abs(row["energy_residual"]), ABS_TOL)
                    self.assertLess(row["maximum_metric_gradient_residual"], ABS_TOL)
                    for field in FIELDS:
                        self.assertLess(
                            row["maximum_metric_gradient_residual_by_field"][field], ABS_TOL
                        )
                    self.assertGreater(abs(row["unmatched_energy"]), 1e-6)
                    self.assertGreater(row["unmatched_maximum_metric_gradient"], 1e-8)
                    if initial is None:
                        self.assertIsNone(row["initial_coefficients"])
                        self.assertTrue(row["coefficients_unspecified"])
                    else:
                        np.testing.assert_allclose(row["initial_coefficients"], initial)
                        self.assertFalse(row["coefficients_unspecified"])

    def test_unspecified_coefficients_are_not_defaulted_to_zero(self):
        metric = cylinder_metric(32)
        row = finite_increment(metric, 2.)
        self.assertIsNone(row["initial_coefficients"])
        self.assertIsNone(row["shifted_coefficients"])
        self.assertTrue(row["coefficients_unspecified"])
        self.assertFalse(np.allclose(row["delta_coefficients"], 0))
        with self.assertRaises(ValueError):
            shift_coefficients(None, 2.)

    def test_centered_derivative_is_independent_of_the_beta_formula(self):
        metric = smooth_metric(32, general=True)
        analytic = remainder_log_derivatives(metric)
        observed = centered_normalization_derivative(metric, 2., 1., 1e-4)
        self.assertAlmostEqual(observed["energy"], analytic["energy"], delta=FINITE_DIFFERENCE_ATOL)
        np.testing.assert_allclose(
            observed["metric_gradients"], analytic["metric_gradients"],
            atol=FINITE_DIFFERENCE_ATOL,
        )

    def test_composition_and_inverse_shifts(self):
        np.testing.assert_allclose(
            finite_coefficient_shift(2.) + finite_coefficient_shift(3.),
            finite_coefficient_shift(6.),
            atol=1e-15,
        )
        np.testing.assert_allclose(
            finite_coefficient_shift(2.) + finite_coefficient_shift(0.5),
            np.zeros(len(BASIS)),
            atol=1e-15,
        )
        np.testing.assert_allclose(finite_coefficient_shift(1.), np.zeros(len(BASIS)), atol=1e-15)
        initial = np.array([1., 2., 3., 4., 5., 6.])
        np.testing.assert_allclose(
            shift_coefficients(shift_coefficients(initial, 2.), 0.5), initial, atol=1e-15
        )
        c2 = c2_only_coefficient_shift(np.e)
        full = finite_coefficient_shift(np.e)
        self.assertAlmostEqual(c2[2], full[2], delta=1e-16)
        np.testing.assert_allclose(np.delete(c2, 2), 0, atol=1e-16)

    def test_discrete_euler_residual_vanishes_under_refinement(self):
        rows = discrete_euler_refinement()
        self.assertEqual([row["points"] for row in rows], [32, 64, 128])
        self.assertGreater(rows[0]["c2_only_gradient_residual"], 1e-12)
        self.assertLess(rows[1]["c2_only_gradient_residual"], 1e-12)
        self.assertLess(rows[2]["c2_only_gradient_residual"], 1e-12)
        for row in rows:
            self.assertLess(abs(row["full_energy_residual"]), ABS_TOL)
            self.assertLess(row["full_gradient_residual"], ABS_TOL)
            self.assertLess(abs(row["E4_energy"]), 1e-12)

    def test_reference_mass_is_not_the_normalization(self):
        metric = smooth_metric(32, general=True)
        one = remainder_log_derivatives(metric, 1.)
        two = remainder_log_derivatives(metric, 2.)
        self.assertAlmostEqual(one["A4"], two["A4"], delta=1e-14)
        self.assertAlmostEqual(one["energy"], two["energy"], delta=1e-14)
        self.assertNotEqual(one["reference_mass"], two["reference_mass"])

    def test_r2_shift_vanishes_and_c2_runs(self):
        delta = finite_coefficient_shift(np.e)
        self.assertEqual(list(BASIS), ["M4", "M2R", "C2", "R2", "E4", "boxR"])
        self.assertAlmostEqual(delta[0], 0., delta=1e-16)
        self.assertAlmostEqual(delta[1], 0., delta=1e-16)
        self.assertAlmostEqual(delta[3], 0., delta=1e-16)
        self.assertGreater(abs(delta[2]), 1e-4)
        self.assertGreater(abs(delta[4]), 1e-5)
        self.assertGreater(abs(delta[5]), 1e-5)

    def test_invalid_scales_and_coefficient_shapes_are_rejected(self):
        with self.assertRaises(ValueError):
            finite_coefficient_shift(0.)
        with self.assertRaises(ValueError):
            finite_coefficient_shift(-2.)
        with self.assertRaises(ValueError):
            shift_coefficients(np.ones(5), 2.)
        with self.assertRaises(ValueError):
            remainder_log_derivatives(cylinder_metric(32), 0.)


class RecordContractTests(unittest.TestCase):
    def test_comparator_keeps_primitive_fields_tight(self):
        compare_record({"energy_residual": 0.0}, {"energy_residual": 1e-16})
        with self.assertRaises(RuntimeError):
            compare_record({"energy_residual": 0.0}, {"energy_residual": 1e-11})
        with self.assertRaises(RuntimeError):
            compare_record({"A4": 1.0}, {"A4": 1.0 + 1e-11})
        self.assertEqual(float_tolerances("$/families/0/energy_residual"), (ABS_TOL, 2e-12))

    def test_comparator_uses_cancellation_aware_finite_difference_tolerance(self):
        compare_record({"energy_minus_analytic": 0.0}, {"energy_minus_analytic": 1e-11})
        compare_record(
            {"centered_differences": [{"maximum_gradient_minus_analytic": 0.0}]},
            {"centered_differences": [{"maximum_gradient_minus_analytic": 1.5e-11}]},
        )
        with self.assertRaises(RuntimeError):
            compare_record({"energy_minus_analytic": 0.0}, {"energy_minus_analytic": 3e-10})
        self.assertEqual(
            float_tolerances("$/families/0/derivative_controls/centered_differences/2/energy_minus_analytic"),
            (FINITE_DIFFERENCE_ATOL, 2e-12),
        )

    def test_reproduction_rejects_changed_scope_and_hashes(self):
        for left, right in (
            ({"complete_ultraviolet_functional_fixed": False},
             {"complete_ultraviolet_functional_fixed": True}),
            ({"M_ratio_identified_with_Omega": False},
             {"M_ratio_identified_with_Omega": True}),
            ({"physical_Newton_running_inferred": False},
             {"physical_Newton_running_inferred": True}),
            ({"physical_cosmological_running_inferred": False},
             {"physical_cosmological_running_inferred": True}),
            ({"zero_action_selected_as_baseline": False},
             {"zero_action_selected_as_baseline": True}),
            ({"continuum_Euler_bulk_stress_from_this_logarithm": False},
             {"continuum_Euler_bulk_stress_from_this_logarithm": True}),
            ({"finite_coefficients_selected": False},
             {"finite_coefficients_selected": True}),
            ({"hash": "abc"}, {"hash": "abd"}),
            ({"energy_residual": -1.2e-16}, {"energy_residual": 0.5}),
            ({"a": 1}, {"a": 1, "b": 2}),
            ({"initial_coefficients": None}, {"initial_coefficients": [0.0] * 6}),
        ):
            with self.assertRaises(RuntimeError):
                compare_record(left, right)

    def test_cli_rejects_combined_modes_and_existing_output_without_overwriting(self):
        script = ROOT / "scripts/check_nsc_normalization_flow.py"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "existing.json"
            path.write_text("preserved bytes\n")
            for args in (("--check", "--output", str(path)), ("--output", str(path))):
                result = subprocess.run(
                    [sys.executable, str(script), *args], capture_output=True, text=True
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(path.read_text(), "preserved bytes\n")

    def test_frozen_record_reproduces_all_fields(self):
        path = ROOT / "results/nsc-12-normalization-flow.json"
        if not path.exists():
            self.skipTest("exclusive result has not been written")
        compare_record(json.loads(path.read_text()), build_record())

    def test_build_record_flags_remain_unfitted_and_not_omega(self):
        identities = exact_flow_identities()
        self.assertTrue(identities["M_ratio_is_not_recursive_Omega"])
        self.assertEqual(identities["M_ref_basis_M4_M2R_betas"], ["0", "0"])
        self.assertFalse(identities["a4_contains_R2"])
        self.assertEqual(identities["unspecified_initial_coefficients"], "None")


if __name__ == "__main__":
    unittest.main()
