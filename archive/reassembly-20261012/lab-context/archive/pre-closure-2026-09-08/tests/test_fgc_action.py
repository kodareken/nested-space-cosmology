from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc import (  # noqa: E402
    ActionParameters,
    ModelID,
    action_certificate,
    coupling_functions,
    dimensionless_scan_coordinates,
    effective_planck_mass_squared,
    linear_effective_mass_squared,
    required_nonclaims,
    scalar_equation_algebraic_residual,
    schwarzschild_activation_radius,
    schwarzschild_gauss_bonnet,
)


class FGCActionTests(unittest.TestCase):
    def setUp(self) -> None:
        shared = dict(
            planck_mass=2.0,
            scalar_mass=3.0,
            quartic_coupling=0.5,
            pulse_width=4.0,
        )
        self.gr = ActionParameters(ModelID.GR_0, **shared)
        self.linear = ActionParameters(
            ModelID.SGB_L, **shared, linear_gb_coupling=-0.25
        )
        self.qr = ActionParameters(
            ModelID.FGC_QR,
            **shared,
            scalar_field=1.0,
            ricci_coupling=-0.25,
            quadratic_gb_coupling=0.5,
        )

    def test_exact_branch_couplings_and_scalar_equation(self) -> None:
        gr = coupling_functions(self.gr, 0.0)
        linear = coupling_functions(self.linear, 2.0)
        qr = coupling_functions(self.qr, 1.0)
        self.assertEqual(gr["F"], 4.0)
        self.assertEqual(gr["f_gb"], 0.0)
        self.assertEqual(linear["f_gb"], -0.5)
        self.assertEqual(linear["f_gb_prime"], -0.25)
        self.assertEqual(linear["f_gb_second"], 0.0)
        self.assertEqual(qr["F"], 3.75)
        self.assertEqual(qr["F_prime"], -0.5)
        self.assertEqual(qr["F_second"], -0.5)
        self.assertEqual(qr["V"], 4.625)
        self.assertEqual(qr["V_prime"], 9.5)
        self.assertEqual(qr["V_second"], 10.5)
        self.assertEqual(qr["f_gb"], 0.0625)
        self.assertEqual(qr["f_gb_prime"], 0.125)
        self.assertEqual(qr["f_gb_second"], 0.125)
        self.assertEqual(effective_planck_mass_squared(self.qr), 3.75)
        self.assertEqual(
            scalar_equation_algebraic_residual(
                self.linear, gauss_bonnet=8.0, phi=0.0
            ),
            -2.0,
        )

    def test_linear_mass_and_schwarzschild_threshold_match_declared_action(self) -> None:
        self.assertEqual(
            linear_effective_mass_squared(
                self.qr, ricci_scalar=2.0, gauss_bonnet=8.0, phi=0.0
            ),
            8.5,
        )
        radius = schwarzschild_activation_radius(self.qr, schwarzschild_radius=8.0)
        expected = (3.0 * 0.5 * 8.0**2 / 3.0**2) ** (1.0 / 6.0)
        self.assertAlmostEqual(radius, expected)
        invariant = schwarzschild_gauss_bonnet(8.0, radius)
        self.assertAlmostEqual(0.5 * invariant / 4.0, 3.0**2)
        self.assertAlmostEqual(
            linear_effective_mass_squared(
                self.qr, gauss_bonnet=invariant, phi=0.0
            ),
            0.0,
            places=12,
        )
        self.assertEqual(
            scalar_equation_algebraic_residual(
                self.qr, gauss_bonnet=invariant, phi=0.0
            ),
            0.0,
        )
        doubled_parent = schwarzschild_activation_radius(self.qr, 16.0)
        self.assertAlmostEqual(doubled_parent / radius, 2.0 ** (1.0 / 3.0))

    def test_dimensionless_coordinates_and_invalid_parameters(self) -> None:
        scan = dimensionless_scan_coordinates(
            self.qr, phi=1.0, ricci_scalar=2.0, gauss_bonnet=3.0
        )
        self.assertEqual(scan["phi_over_planck_mass"], 0.5)
        self.assertEqual(scan["mu_times_L0"], 12.0)
        self.assertEqual(scan["eta_over_L0_squared"], 0.03125)
        self.assertEqual(scan["R_times_L0_squared"], 32.0)
        self.assertEqual(scan["G_times_L0_fourth"], 768.0)

        shared = dict(
            planck_mass=1.0,
            scalar_mass=1.0,
            quartic_coupling=1.0,
            pulse_width=1.0,
        )
        invalid = [
            (ModelID.GR_0, {**shared, "planck_mass": 0.0}),
            (ModelID.GR_0, {**shared, "scalar_mass": -1.0}),
            (ModelID.GR_0, {**shared, "quartic_coupling": 0.0}),
            (ModelID.GR_0, {**shared, "pulse_width": math.inf}),
            (ModelID.GR_0, {**shared, "ricci_coupling": 0.1}),
            (ModelID.SGB_L, shared),
            (ModelID.FGC_QR, shared),
            (ModelID.FGC_QR, {**shared, "ricci_coupling": 0.1}),
            (
                ModelID.FGC_QR,
                {
                    **shared,
                    "ricci_coupling": 0.1,
                    "quadratic_gb_coupling": -0.1,
                },
            ),
        ]
        for model_id, kwargs in invalid:
            with self.subTest(model_id=model_id, kwargs=kwargs):
                with self.assertRaises(ValueError):
                    ActionParameters(model_id, **kwargs)
        with self.assertRaises(ValueError):
            ActionParameters("GR-0", **shared)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            schwarzschild_gauss_bonnet(1.0, 0.0)
        with self.assertRaises(ValueError):
            schwarzschild_activation_radius(self.gr, 1.0)
        with self.assertRaises(ValueError):
            dimensionless_scan_coordinates(self.gr, areal_radius=1.0)
        with self.assertRaises(ValueError):
            coupling_functions(self.qr, 1.0e308)
        with self.assertRaises(ValueError):
            dimensionless_scan_coordinates(self.qr, gauss_bonnet=1.0e308)
        failed_planck = ActionParameters(
            ModelID.FGC_QR,
            **shared,
            scalar_field=3.0,
            ricci_coupling=-0.5,
            quadratic_gb_coupling=0.5,
        )
        with self.assertRaisesRegex(ValueError, "effective Planck mass squared"):
            action_certificate(failed_planck)

    def test_certificate_is_json_safe_and_fail_closed(self) -> None:
        radius = schwarzschild_activation_radius(self.qr, 8.0)
        payload = action_certificate(
            self.qr,
            ricci_scalar=0.0,
            schwarzschild_radius=8.0,
            areal_radius=radius,
        )
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(payload["model_id"], "FGC-QR")
        self.assertEqual(
            payload["classification"],
            "covariant_action_and_linear_background_certificate_not_collapse_solution",
        )
        self.assertTrue(payload["action_gate"]["candidate_covariant_action_declared"])
        self.assertTrue(payload["action_gate"]["separate_matter_and_regulator_fields"])
        self.assertEqual(
            payload["low_gradient_reference_gate"]["effective_planck_ratio"],
            0.9375,
        )
        self.assertFalse(
            payload["low_gradient_reference_gate"][
                "solution_level_gr_recovery_demonstrated"
            ]
        )
        self.assertFalse(
            payload["low_gradient_reference_gate"][
                "zero_regulator_solves_declared_scalar_equation"
            ]
        )
        self.assertFalse(
            payload["low_gradient_reference_gate"][
                "coupled_metric_chi_background_verified"
            ]
        )
        self.assertEqual(payload["nonclaims"], required_nonclaims())
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))
        self.assertFalse(
            payload["weak_field_kinetic_reference"][
                "full_coupled_quadratic_action_diagonalized"
            ]
        )
        activation = payload["schwarzschild_linear_activation_control"]
        self.assertAlmostEqual(
            activation["linear_effective_mass_squared_at_activation"],
            0.0,
            places=12,
        )
        self.assertEqual(
            activation["interpretation"],
            "linear_instability_scale_not_transition_or_eft_cutoff",
        )
        json.dumps(payload, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
