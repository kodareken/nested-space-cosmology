from __future__ import annotations

from copy import deepcopy
import json
import unittest

from scripts import reproduce_fgc_tdg5_pref20 as pref20


class FGCTDG5PREF20ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = pref20.verify_canonical()

    def test_canonical_theorem_binder_reproduces(self) -> None:
        self.assertEqual(self.record["artifact_id"], pref20.ARTIFACT_ID)
        self.assertEqual(
            self.record["classification"],
            "completed_TDG5_stage_complete_refinement_theorem_runtime_pair_implementation_authorized",
        )
        status = self.record["gate_status"]
        self.assertTrue(
            status["TDG5_stage_complete_refinement_theorem_completed"]
        )
        self.assertTrue(
            status["runtime_refinement_pair_implementation_authorized"]
        )
        self.assertFalse(status["runtime_refinement_pair_implemented"])

    def test_exact_theorem_rederives_map_restriction_extrema_and_debits(self) -> None:
        theorem = self.record["artifact_payload"]["theorem_certificate"]
        self.assertTrue(theorem["all_independent_theorem_controls_pass"])
        hermite = theorem["independent_Hermite_map"]
        self.assertEqual(
            hermite["sampling_determinant_by_fraction_Gaussian_elimination"],
            "1",
        )
        self.assertEqual(hermite["condition_number_infinity_upper_bound"], "54")
        restriction = theorem["half_interval_restriction"]
        self.assertEqual(restriction["left_restriction_determinant"], "1/64")
        self.assertEqual(restriction["right_restriction_determinant"], "1/64")
        extrema = theorem["continuous_extremum_theorem"]
        self.assertTrue(
            extrema["continuous_piecewise_cubic_extremum_theorem_completed"]
        )
        self.assertEqual(
            extrema["endpoint_only_bump_control"]["absolute_maximum"], "1"
        )
        self.assertEqual(
            theorem["method_owned_Richardson_debits"]["RK4"][
                "Richardson_denominator"
            ],
            15,
        )
        self.assertEqual(
            theorem["method_owned_Richardson_debits"]["SSPRK3"][
                "Richardson_denominator"
            ],
            7,
        )

    def test_engine_shape_exposes_fresh_endpoint_RHS_without_runtime_pair(self) -> None:
        engine = self.record["artifact_payload"]["theorem_certificate"][
            "immutable_numerical_engine_shape"
        ]
        self.assertTrue(engine["stage_and_proposal_fields_match_frozen_shape"])
        self.assertTrue(
            engine["all_required_internal_and_endpoint_stage_names_present"]
        )
        self.assertTrue(
            engine["candidate_endpoint_call_uses_final_time_and_candidate_state"]
        )
        self.assertTrue(engine["fresh_candidate_endpoint_RHS_record_is_available"])
        self.assertFalse(engine["atomic_coarse_fine_transaction_implemented"])

    def test_authorization_is_limited_to_implementation_and_synthetic_validation(self) -> None:
        payload = self.record["artifact_payload"]
        authorization = payload["runtime_authorization"]
        boundary = payload["claim_boundary"]
        self.assertEqual(
            authorization["authorized_object"],
            "implementation_and_synthetic_validation_of_the_frozen_same_grid_refinement_pair",
        )
        self.assertFalse(authorization["production_trajectory_authorized"])
        self.assertTrue(boundary["historical_PROTO13_stop_preserved"])
        self.assertFalse(boundary["historical_CAL10_result_reclassified"])
        self.assertFalse(boundary["declared_cubic_is_exact_PDE_history"])
        self.assertFalse(boundary["trajectory_asymptotic_regime_proved"])
        self.assertFalse(boundary["rigorous_global_PDE_error_bound_proved"])
        self.assertTrue(
            boundary["runtime_refinement_pair_implementation_authorized"]
        )
        for name in (
            "runtime_refinement_pair_implemented",
            "runtime_thresholds_frozen",
            "replacement_temporal_admission_defined",
            "new_production_trajectory_authorized",
            "GR0_case_eligible",
            "candidate_branch_opened",
            "physical_question_answered",
        ):
            self.assertFalse(boundary[name], name)

    def test_tracked_result_is_canonical(self) -> None:
        raw = pref20.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, pref20._canonical_bytes(stored))

    def test_config_mutations_fail_closed(self) -> None:
        config = pref20.load_config()
        attacks = []

        value = deepcopy(config)
        value["immutable_lineage"]["freeze_result_sha256"] = "0" * 64
        attacks.append(value)

        value = deepcopy(config)
        value["theorem_result"]["sampling_determinant"] = "0"
        attacks.append(value)

        value = deepcopy(config)
        value["theorem_result"]["fresh_candidate_endpoint_RHS_record_is_available"] = (
            False
        )
        attacks.append(value)

        value = deepcopy(config)
        value["claims"]["runtime_refinement_pair_implemented"] = True
        attacks.append(value)

        value = deepcopy(config)
        value["claims"]["runtime_thresholds_frozen"] = True
        attacks.append(value)

        value = deepcopy(config)
        value["claims"]["replacement_temporal_admission_defined"] = True
        attacks.append(value)

        value = deepcopy(config)
        value["claims"]["GR0_case_eligible"] = True
        attacks.append(value)

        value = deepcopy(config)
        value["claims"]["FGCQR_holdout_execution_authorized"] = True
        attacks.append(value)

        for attacked in attacks:
            with self.assertRaises(ValueError):
                pref20.validate_config_data(attacked)


if __name__ == "__main__":
    unittest.main()
