from __future__ import annotations

from copy import deepcopy
import json
import unittest

from scripts import reproduce_fgc_tdg5_frz1 as freeze


class FGCTDG5FRZ1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = freeze.verify_canonical()

    def test_canonical_prospective_freeze_reproduces(self) -> None:
        self.assertEqual(self.record["artifact_id"], freeze.ARTIFACT_ID)
        self.assertEqual(
            self.record["classification"],
            "prospective_TDG5_stage_complete_temporal_refinement_design_freeze",
        )
        status = self.record["gate_status"]
        self.assertTrue(status["TDG5_replacement_design_frozen"])
        self.assertTrue(status["TDG5_exact_no_history_controls_pass"])
        self.assertFalse(
            status["TDG5_stage_complete_refinement_theorem_completed"]
        )

    def test_source_theorem_and_historical_stop_are_preserved(self) -> None:
        lineage = self.record["artifact_payload"]["immutable_lineage"]
        self.assertTrue(lineage["source_compact_result_verified"])
        self.assertTrue(lineage["historical_PROTO13_stop_preserved"])
        self.assertFalse(lineage["historical_CAL10_result_reclassified"])
        self.assertFalse(lineage["actual_terminal_history_arrays_consumed"])
        self.assertFalse(lineage["campaign_checkpoint_loaded"])
        self.assertFalse(lineage["state_advanced"])

    def test_exact_preflight_covers_data_map_alias_and_both_methods(self) -> None:
        preflight = self.record["artifact_payload"]["exact_no_history_preflight"]
        self.assertTrue(preflight["all_exact_controls_pass"])
        hermite = preflight["Hermite_data_map"]
        self.assertEqual(hermite["determinant"], "1")
        self.assertTrue(hermite["left_inverse_exact"])
        self.assertTrue(hermite["right_inverse_exact"])
        self.assertEqual(hermite["condition_number_infinity_upper_bound"], "54")
        alias = preflight["endpoint_only_alias_control"]
        self.assertTrue(alias["endpoint_values_equal"])
        self.assertEqual(alias["adversarial_midpoint_value"], "1")
        methods = preflight["method_controls"]
        self.assertEqual(set(methods), {"RK4", "SSPRK3"})
        self.assertEqual(
            methods["RK4"]["step_doubling"]["Richardson_denominator"],
            15,
        )
        self.assertEqual(
            methods["SSPRK3"]["step_doubling"]["Richardson_denominator"],
            7,
        )

    def test_design_freeze_promotes_no_runtime_or_physical_claim(self) -> None:
        boundary = self.record["artifact_payload"]["claim_boundary"]
        self.assertTrue(boundary["replacement_design_frozen"])
        for name in (
            "stage_complete_refinement_theorem_completed",
            "runtime_refinement_pair_implementation_authorized",
            "runtime_refinement_pair_implemented",
            "runtime_thresholds_frozen",
            "replacement_temporal_admission_defined",
            "new_trajectory_authorized",
            "GR0_case_eligible",
            "candidate_branch_opened",
            "physical_question_answered",
        ):
            self.assertFalse(boundary[name], name)

    def test_tracked_result_is_canonical(self) -> None:
        raw = freeze.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, freeze._canonical_bytes(stored))

    def test_config_mutations_fail_closed(self) -> None:
        config = freeze.load_config()
        attacks = []

        value = deepcopy(config)
        value["immutable_lineage"]["source_result_sha256"] = "0" * 64
        attacks.append(value)

        value = deepcopy(config)
        value["frozen_replacement_design"]["monitored_state"] = (
            "64_normal_flow_tracers"
        )
        attacks.append(value)

        value = deepcopy(config)
        value["frozen_replacement_design"]["Hermite_continuous_extension"][
            "sampling_matrix_determinant"
        ] = "0"
        attacks.append(value)

        value = deepcopy(config)
        value["frozen_replacement_design"]["same_grid_step_refinement"][
            "SSPRK3_Richardson_denominator"
        ] = 8
        attacks.append(value)

        value = deepcopy(config)
        value["claims"]["runtime_refinement_pair_implementation_authorized"] = (
            True
        )
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
                freeze.validate_config_data(attacked)


if __name__ == "__main__":
    unittest.main()
