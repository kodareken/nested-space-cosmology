from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.evidence import reproduce_controlled_model  # noqa: E402


class ControlledModelEvidenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.record = reproduce_controlled_model()

    def test_record_is_finite_and_json_serializable(self) -> None:
        payload = json.dumps(self.record, allow_nan=False)
        self.assertIn("controlled_covariant_core_not_global_black_hole_transition", payload)
        self.assertIn("chosen_symmetric_benchmark_not_derived", payload)
        self.assertIn("classical_core_transfer_not_primordial_spectrum", payload)
        self.assertIn("not_dark_sector_or_recursive_transition", payload)

    def test_all_six_required_artifacts_have_bounded_status(self) -> None:
        self.assertEqual(self.record["project_version"], "0.11.0")
        self.assertEqual(self.record["schema_version"], 3)
        surfaces = self.record["covariant_core"]["surface_examples"]
        self.assertEqual(
            surfaces["contracting_equator_t_over_L_minus_2"], "trapped"
        )
        self.assertEqual(surfaces["bounce_quarter_sphere"], "normal")
        self.assertEqual(
            surfaces["expanding_equator_t_over_L_plus_2"], "anti_trapped"
        )

        audit = self.record["information_map"]["audit"]
        self.assertTrue(audit["is_isometry"])
        self.assertTrue(audit["is_trace_preserving"])
        self.assertTrue(audit["is_recoverable_on_code_range"])

        perturbations = self.record["linear_perturbations"]
        self.assertLess(
            perturbations["scalar_exact_regression"]["max_abs_error_from_exact_rotation"],
            1e-10,
        )
        self.assertLess(
            perturbations["tensor_refinement"]["final_determinant_error"], 1e-10
        )

        domain = self.record["domain_synthesis"]
        self.assertEqual(domain["model_id"], "DST-1")
        self.assertEqual(
            domain["quadratic_cycle_average"]["w_oscillatory"], 0.0
        )
        self.assertEqual(
            domain["quadratic_monomial_limit"]["density_scale_exponent"], 3.0
        )
        self.assertLess(
            abs(domain["fixed_h_spectator_endpoint"]["endpoint_residual"]),
            1e-10,
        )
        self.assertLess(
            abs(domain["exchange_ledger"]["total_continuity_residual"]),
            1e-12,
        )
        self.assertLess(
            abs(domain["kottler_balance"]["acceleration_at_balance_m_s^-2"]),
            1e-18,
        )
        self.assertGreater(
            domain["kottler_balance"]["radial_slope_at_balance_s^-2"], 0.0
        )
        self.assertEqual(
            domain["kottler_balance"]["stability"],
            "unstable_radial_separator",
        )

        discriminator = self.record["independent_discriminator"]
        self.assertEqual(discriminator["prediction"], "Omega_K<0_when_H_is_nonzero")
        self.assertEqual(
            discriminator["required_branch_postulate"],
            "observable_child_remains_globally_closed_with_K=+1",
        )
        self.assertLess(
            discriminator["illustrative_omega_K_not_present_day_prediction"], 0.0
        )

        gates = self.record["literal_model_gates"]
        self.assertEqual(
            gates["classification"],
            "conditional_literal_subclass_rejection_gates_not_model_confirmation",
        )
        self.assertEqual(gates["model_id"], "DSF-2_DEB-1")
        self.assertEqual(
            gates["external_benchmarks"]["GW170817_GRB170817A"]["source"],
            "https://arxiv.org/abs/1710.05834",
        )
        self.assertEqual(
            gates["external_benchmarks"]["DESI_DR2_flat_wCDM"]["source"],
            "https://arxiv.org/abs/2503.14738v3",
        )
        self.assertEqual(
            gates["external_benchmarks"]["DESI_DR2_flat_wCDM"]["published_doi"],
            "10.1103/tr6y-kpc6",
        )
        cone = gates["relative_cone"]
        self.assertFalse(cone["null_control"]["rejected"])
        self.assertTrue(cone["null_control"]["within_inclusive_bound"])
        self.assertTrue(cone["positive_1e_minus_14"]["rejected"])
        self.assertTrue(cone["negative_1e_minus_14"]["rejected"])
        self.assertEqual(
            cone["disformal_stationary_null_control"]["delta_tensor_over_photon"],
            0.0,
        )

        boundary = gates["boundary_scalings"]
        constant = boundary["constant_density"]
        self.assertEqual(constant["effective_w"], -1.0)
        self.assertTrue(constant["accelerates"])
        self.assertTrue(
            all(
                gate["within_approximate_interval"]
                for gate in constant["desi_dr2_flat_wcdm_summary_gates"].values()
            )
        )
        surface = boundary["constant_tension_surface"]
        self.assertAlmostEqual(surface["effective_w"], -2.0 / 3.0)
        self.assertTrue(surface["accelerates"])
        self.assertTrue(
            all(
                not gate["within_approximate_interval"]
                for gate in surface["desi_dr2_flat_wcdm_summary_gates"].values()
            )
        )
        inverse_area = boundary["inverse_area_density"]
        self.assertAlmostEqual(inverse_area["effective_w"], -1.0 / 3.0)
        self.assertFalse(inverse_area["accelerates"])
        self.assertTrue(
            all(
                not gate["within_approximate_interval"]
                for gate in inverse_area["desi_dr2_flat_wcdm_summary_gates"].values()
            )
        )

    def test_numeric_tree_contains_no_nonfinite_float(self) -> None:
        def visit(value: object) -> None:
            if isinstance(value, float):
                self.assertTrue(math.isfinite(value))
            elif isinstance(value, dict):
                for child in value.values():
                    visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value:
                    visit(child)

        visit(self.record)


if __name__ == "__main__":
    unittest.main()
