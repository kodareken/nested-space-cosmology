from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg5_stage_complete_refinement import (
    TDG5_METHOD_ORDERS,
    butcher_order_certificate,
    endpoint_only_alias_control,
    evaluate_cubic,
    hermite_coefficients,
    hermite_data_map_certificate,
    stage_complete_temporal_refinement_preflight,
    step_doubling_certificate,
)


Q = Fraction


class TDG5StageCompleteRefinementTests(unittest.TestCase):
    def test_Hermite_data_map_is_exact_injective_and_quantitatively_stable(self) -> None:
        evidence = hermite_data_map_certificate()
        self.assertEqual(evidence["determinant"], "1")
        self.assertTrue(evidence["left_inverse_exact"])
        self.assertTrue(evidence["right_inverse_exact"])
        self.assertTrue(evidence["data_map_is_injective"])
        self.assertTrue(evidence["no_between_sample_kernel_inside_declared_class"])
        self.assertEqual(evidence["sampling_infinity_norm"], "6")
        self.assertEqual(evidence["inverse_infinity_norm"], "9")
        self.assertEqual(evidence["condition_number_infinity_upper_bound"], "54")

    def test_exact_Hermite_round_trip_recovers_values_slopes_and_midpoint(self) -> None:
        data = (Q(2), Q(-1), Q(5), Q(3))
        coefficients = hermite_coefficients(data)
        self.assertEqual(evaluate_cubic(coefficients, Q(0)), data[0])
        self.assertEqual(evaluate_cubic(coefficients, Q(1)), data[2])
        derivative_left = coefficients[1]
        derivative_right = coefficients[1] + 2 * coefficients[2] + 3 * coefficients[3]
        self.assertEqual(derivative_left, data[1])
        self.assertEqual(derivative_right, data[3])
        self.assertEqual(evaluate_cubic(coefficients, Q(1, 2)), Q(3))

    def test_endpoint_only_alias_is_detected_by_stage_complete_data(self) -> None:
        evidence = endpoint_only_alias_control()
        self.assertTrue(evidence["endpoint_values_equal"])
        self.assertTrue(evidence["endpoint_slopes_differ"])
        self.assertEqual(evidence["adversarial_midpoint_value"], "1")
        self.assertTrue(evidence["endpoint_only_sampling_misses_nonzero_interior"])
        self.assertTrue(evidence["stage_complete_Hermite_data_detects_it"])

    def test_both_inherited_methods_satisfy_their_exact_order_conditions(self) -> None:
        expected = dict(TDG5_METHOD_ORDERS)
        for method, order in expected.items():
            evidence = butcher_order_certificate(method)
            self.assertEqual(evidence["formal_order"], order)
            self.assertTrue(evidence["all_required_order_conditions_pass"])
            self.assertEqual(
                set(item["passed"] for item in evidence["conditions"].values()),
                {True},
            )

    def test_step_doubling_factors_and_exact_scalar_contraction(self) -> None:
        rk4 = step_doubling_certificate("RK4")
        ssprk3 = step_doubling_certificate("SSPRK3")
        self.assertEqual(rk4["Richardson_denominator"], 15)
        self.assertEqual(ssprk3["Richardson_denominator"], 7)
        for evidence in (rk4, ssprk3):
            self.assertTrue(
                evidence["contraction_is_stricter_than_two_to_minus_formal_order"]
            )
            self.assertTrue(evidence["same_initial_state_and_same_spatial_operator_required"])
            self.assertTrue(evidence["fine_two_half_step_path_is_the_only_committable_path"])
            self.assertTrue(evidence["Richardson_factor_is_not_a_rigorous_PDE_error_bound"])

    def test_preflight_selects_one_sufficient_route_without_promoting_a_gate(self) -> None:
        evidence = stage_complete_temporal_refinement_preflight()
        self.assertTrue(evidence["all_exact_controls_pass"])
        self.assertIn("finite_dimensional", evidence["selected_sufficient_route"])
        self.assertEqual(set(evidence["method_controls"]), {"RK4", "SSPRK3"})
        self.assertFalse(evidence["actual_terminal_histories_consumed"])
        self.assertFalse(evidence["campaign_checkpoint_loaded"])
        self.assertFalse(evidence["state_advanced"])
        self.assertFalse(evidence["runtime_refinement_pair_implemented"])
        self.assertFalse(evidence["runtime_thresholds_frozen"])
        self.assertFalse(evidence["replacement_temporal_admission_defined"])
        self.assertFalse(evidence["PROTO14_frozen"])
        self.assertFalse(evidence["GR0_case_eligible"])

    def test_invalid_methods_data_and_intervals_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "unknown"):
            butcher_order_certificate("Euler")
        with self.assertRaisesRegex(ValueError, "positive"):
            step_doubling_certificate("RK4", probe_step=0)
        with self.assertRaisesRegex(TypeError, "exact rational"):
            step_doubling_certificate("RK4", probe_step=0.125)
        with self.assertRaisesRegex(ValueError, "exactly four"):
            hermite_coefficients((0, 1, 2))
        with self.assertRaisesRegex(TypeError, "exact rational"):
            hermite_coefficients((0, 1, 2, 3.0))
        with self.assertRaisesRegex(ValueError, "closed unit interval"):
            evaluate_cubic((0, 1, 2, 3), Q(3, 2))


if __name__ == "__main__":
    unittest.main()
