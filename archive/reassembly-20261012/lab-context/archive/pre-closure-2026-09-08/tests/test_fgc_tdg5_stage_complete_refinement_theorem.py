from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE = REPOSITORY / "src"
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from recursive_horizons.fgc.evolution.tdg5_stage_complete_refinement_theorem import (
    cubic_extremum_theorem_certificate,
    exact_cubic_absolute_maximum,
    exact_derivative_root_candidates,
    half_interval_restriction_certificate,
    independent_hermite_coefficients,
    independent_hermite_map_certificate,
    numerical_engine_shape_certificate,
    restrict_cubic_to_half,
    richardson_debit_certificate,
    stage_complete_refinement_theorem_certificate,
)


Q = Fraction
ENGINE = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/numerical_engine.py"
)


class TDG5StageCompleteRefinementTheoremTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.engine_source = ENGINE.read_text(encoding="utf-8")

    def test_independent_Hermite_expansion_is_injective_and_stable(self) -> None:
        evidence = independent_hermite_map_certificate()
        self.assertEqual(
            independent_hermite_coefficients((2, -1, 5, 3)),
            (Q(2), Q(-1), Q(8), Q(-4)),
        )
        self.assertEqual(
            evidence["sampling_determinant_by_fraction_Gaussian_elimination"],
            "1",
        )
        self.assertEqual(
            evidence["basis_determinant_by_fraction_Gaussian_elimination"],
            "1",
        )
        self.assertTrue(evidence["all_four_unit_data_round_trip_exactly"])
        self.assertEqual(evidence["sampling_infinity_norm"], "6")
        self.assertEqual(evidence["basis_infinity_norm"], "9")
        self.assertEqual(evidence["condition_number_infinity_upper_bound"], "54")
        self.assertTrue(evidence["finite_dimensional_data_map_is_stably_injective"])

    def test_both_half_interval_restrictions_are_exact_and_injective(self) -> None:
        evidence = half_interval_restriction_certificate()
        self.assertEqual(evidence["left_restriction_determinant"], "1/64")
        self.assertEqual(evidence["right_restriction_determinant"], "1/64")
        self.assertEqual(evidence["control_count"], 8)
        self.assertTrue(evidence["all_basis_and_probe_composition_controls_pass"])
        self.assertEqual(
            restrict_cubic_to_half((1, 2, 3, 4), half=0),
            (Q(1), Q(1), Q(3, 4), Q(1, 2)),
        )
        self.assertEqual(
            restrict_cubic_to_half((1, 2, 3, 4), half=1),
            (Q(13, 4), Q(4), Q(9, 4), Q(1, 2)),
        )

    def test_continuous_extrema_use_endpoints_and_derivative_roots(self) -> None:
        evidence = cubic_extremum_theorem_certificate()
        self.assertTrue(evidence["continuous_piecewise_cubic_extremum_theorem_completed"])
        self.assertEqual(evidence["maximum_candidate_count_per_half_interval"], 4)
        self.assertEqual(
            evidence["endpoint_only_bump_control"]["absolute_maximum"],
            "1",
        )
        self.assertEqual(
            evidence["endpoint_only_bump_control"]["maximizers"],
            ["1/2"],
        )
        self.assertEqual(
            exact_derivative_root_candidates((0, 4, -4, 0)),
            (Q(1, 2),),
        )
        self.assertEqual(
            exact_cubic_absolute_maximum((0, 0, 0, 0))["candidate_points"],
            ["0", "1"],
        )
        self.assertTrue(evidence["binary64_root_and_value_evaluation_requires_outward_debit"])
        self.assertFalse(evidence["binary64_extremum_evaluator_implemented"])

    def test_nonrational_extremum_and_invalid_half_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "rational square"):
            exact_derivative_root_candidates((0, Q(-1, 2), 0, 1))
        with self.assertRaisesRegex(ValueError, "half"):
            restrict_cubic_to_half((0, 1, 2, 3), half=2)
        with self.assertRaisesRegex(TypeError, "exact rational"):
            independent_hermite_coefficients((0, 1, 2, 3.0))

    def test_method_owned_Richardson_relation_keeps_its_assumptions_public(self) -> None:
        rk4 = richardson_debit_certificate("RK4", 4)
        ssprk3 = richardson_debit_certificate("SSPRK3", 3)
        self.assertEqual(rk4["Richardson_denominator"], 15)
        self.assertEqual(rk4["fine_error_over_difference_factor"], "1/15")
        self.assertEqual(ssprk3["Richardson_denominator"], 7)
        self.assertEqual(ssprk3["fine_error_over_difference_factor"], "1/7")
        for evidence in (rk4, ssprk3):
            self.assertTrue(
                evidence["derivation_assumes_the_tested_step_is_in_the_asymptotic_regime"]
            )
            self.assertFalse(evidence["trajectory_asymptotic_regime_proved"])
            self.assertFalse(evidence["rigorous_global_PDE_error_bound_proved"])
        with self.assertRaisesRegex(ValueError, "differs"):
            richardson_debit_certificate("RK4", 3)

    def test_immutable_engine_exposes_fresh_endpoint_RHS_record(self) -> None:
        evidence = numerical_engine_shape_certificate(self.engine_source)
        self.assertTrue(evidence["stage_and_proposal_fields_match_frozen_shape"])
        self.assertTrue(evidence["all_required_internal_and_endpoint_stage_names_present"])
        self.assertTrue(
            evidence["candidate_endpoint_call_uses_final_time_and_candidate_state"]
        )
        self.assertTrue(evidence["evaluate_function_calls_RHS"])
        self.assertTrue(evidence["evaluate_function_appends_StageRecord"])
        self.assertTrue(evidence["fresh_candidate_endpoint_RHS_record_is_available"])
        self.assertFalse(evidence["atomic_coarse_fine_transaction_implemented"])

        attacked = self.engine_source.replace(
            'evaluate("candidate_endpoint", final_time, candidate)',
            'evaluate("candidate_shadow", final_time, candidate)',
        )
        rejected = numerical_engine_shape_certificate(attacked)
        self.assertFalse(
            rejected["all_required_internal_and_endpoint_stage_names_present"]
        )
        self.assertFalse(rejected["fresh_candidate_endpoint_RHS_record_is_available"])

    def test_combined_certificate_authorizes_only_runtime_implementation(self) -> None:
        evidence = stage_complete_refinement_theorem_certificate(
            self.engine_source
        )
        self.assertTrue(evidence["all_independent_theorem_controls_pass"])
        self.assertTrue(
            evidence["TDG5_stage_complete_refinement_theorem_completed"]
        )
        self.assertTrue(
            evidence["runtime_refinement_pair_implementation_authorized"]
        )
        for name in (
            "runtime_refinement_pair_implemented",
            "runtime_thresholds_frozen",
            "replacement_temporal_admission_defined",
            "PROTO14_frozen",
            "fresh_GR0_dynamic_calibration_completed",
            "GR0_case_eligible",
            "state_advanced",
            "SGBL_trajectory_read",
            "FGCQR_trajectory_read",
            "physical_or_candidate_question_answered",
        ):
            self.assertFalse(evidence[name], name)


if __name__ == "__main__":
    unittest.main()
