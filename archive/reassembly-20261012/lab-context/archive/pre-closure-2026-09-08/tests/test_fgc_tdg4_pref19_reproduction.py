from __future__ import annotations

from copy import deepcopy
import json
import unittest

from scripts import reproduce_fgc_tdg4_pref19 as pref19


class FGCTDG4PREF19ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = pref19.verify_canonical()

    def test_canonical_theorem_binder_reproduces(self) -> None:
        self.assertEqual(self.record["artifact_id"], pref19.ARTIFACT_ID)
        self.assertEqual(
            self.record["classification"],
            "completed_TDG4_sampling_identifiability_theorem_current_sampled_rule_retired_as_continuum_admission",
        )
        status = self.record["gate_status"]
        self.assertTrue(
            status["TDG4_sampling_identifiability_theorem_completed"]
        )
        self.assertTrue(
            status[
                "finite_samples_alone_proved_insufficient_for_continuum_top_band_bound"
            ]
        )
        self.assertTrue(
            status[
                "current_temporal_spectral_rule_retired_as_continuum_admission"
            ]
        )
        self.assertTrue(status["replacement_temporal_gate_design_authorized"])

    def test_retirement_is_narrow_and_history_is_not_reclassified(self) -> None:
        payload = self.record["artifact_payload"]
        retirement = payload["retirement_scope"]
        boundary = payload["claim_boundary"]
        self.assertTrue(retirement["current_rule_retained_as_sampled_diagnostic"])
        self.assertTrue(retirement["historical_PROTO13_stop_preserved"])
        self.assertFalse(retirement["historical_CAL10_result_reclassified"])
        self.assertFalse(retirement["actual_history_classifications_changed"])
        self.assertTrue(boundary["historical_PROTO13_stop_preserved"])
        self.assertFalse(boundary["historical_CAL10_result_reclassified"])
        self.assertFalse(boundary["actual_unsampled_PROTO13_history_power_inferred"])
        self.assertFalse(boundary["PDE_solution_manifold_restriction_proved"])

    def test_no_successor_gate_or_candidate_is_promoted(self) -> None:
        status = self.record["gate_status"]
        boundary = self.record["artifact_payload"]["claim_boundary"]
        self.assertFalse(status["replacement_temporal_admission_defined"])
        self.assertFalse(status["PROTO14_frozen"])
        self.assertFalse(status["fresh_GR0_dynamic_calibration_completed"])
        self.assertFalse(status["GR0_case_eligible"])
        self.assertFalse(status["SGBL_execution_authorized"])
        self.assertFalse(status["FGCQR_holdout_execution_authorized"])
        self.assertFalse(boundary["new_trajectory_authorized"])
        self.assertFalse(boundary["candidate_branch_opened"])
        self.assertFalse(boundary["physical_question_answered"])

    def test_exact_theorem_and_freeze_crosscheck_cover_every_bin(self) -> None:
        payload = self.record["artifact_payload"]
        theorem = payload["theorem_certificate"]
        crosscheck = payload["freeze_crosscheck"]
        self.assertTrue(theorem["sampling_identifiability_theorem_completed"])
        self.assertEqual(
            set(theorem["smooth_kernel_witnesses"]),
            {"28", "29", "30", "31", "32"},
        )
        self.assertEqual(
            set(theorem["uniform_alias_crosschecks"]),
            {"28", "29", "30", "31", "32"},
        )
        self.assertTrue(crosscheck["all_five_declared_top_bins_crosschecked"])
        self.assertEqual(set(payload["mutation_controls"].values()), {True})

    def test_tracked_result_is_canonical(self) -> None:
        raw = pref19.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, pref19._canonical_bytes(stored))

    def test_config_mutations_fail_closed(self) -> None:
        config = pref19.load_config()
        attacks = []
        value = deepcopy(config)
        value["immutable_lineage"]["freeze_result_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(config)
        value["theorem_result"]["sample_count"] = 63
        attacks.append(value)
        value = deepcopy(config)
        value["theorem_result"][
            "actual_unsampled_power_of_PROTO13_histories_inferred"
        ] = True
        attacks.append(value)
        value = deepcopy(config)
        value["retirement_scope"]["historical_CAL10_result_reclassified"] = True
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
                pref19.validate_config_data(attacked)


if __name__ == "__main__":
    unittest.main()
