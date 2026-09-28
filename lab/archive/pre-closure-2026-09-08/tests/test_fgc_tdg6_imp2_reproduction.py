from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_tdg6_imp2 as imp2  # noqa: E402


class FGCTDG6IMP2ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = imp2.verify_canonical()

    def test_canonical_implementation_record_reproduces(self) -> None:
        self.assertEqual(self.record["artifact_id"], "FGC-1-TDG6-IMP2")
        self.assertEqual(
            self.record["classification"],
            "production_TDG6_three_level_compositor_implemented_and_synthetically_qualified_PROTO14_design_only_authorized",
        )
        status = self.record["gate_status"]
        self.assertTrue(status["production_compositor_implemented"])
        self.assertTrue(
            status["production_compositor_synthetic_qualification_passed"]
        )
        self.assertTrue(status["PROTO14_freeze_design_authorized"])
        self.assertFalse(status["production_trajectory_authorized"])
        self.assertFalse(status["PROTO14_frozen"])
        self.assertFalse(status["GR0_case_eligible"])
        self.assertFalse(status["FGCQR_holdout_execution_authorized"])

    def test_both_methods_bind_seven_paths_and_fine_only_commit(self) -> None:
        methods = self.record["artifact_payload"]["synthetic_qualification"][
            "method_controls"
        ]
        self.assertEqual(methods["RK4"]["total_shadow_stage_records"], 35)
        self.assertEqual(methods["RK4"]["fine_stage_records"], 20)
        self.assertEqual(methods["SSPRK3"]["total_shadow_stage_records"], 28)
        self.assertEqual(methods["SSPRK3"]["fine_stage_records"], 16)
        for evidence in methods.values():
            self.assertEqual(evidence["shadow_proposal_count"], 7)
            self.assertTrue(evidence["complete_admission_passed"])
            self.assertTrue(evidence["real_boundary_unchanged_during_prepare"])
            self.assertFalse(evidence["outer_and_medium_committed"])
            self.assertTrue(evidence["fine_quarter_steps_committed_atomically"])
            self.assertEqual(
                evidence["committed_state_sha256"],
                evidence["fine_shadow_state_sha256"],
            )
            self.assertEqual(
                evidence["committed_tracer_snapshot"],
                evidence["fine_shadow_tracer_snapshot"],
            )

    def test_interval_retry_checkpoint_and_failure_controls_are_bound(self) -> None:
        synthetic = self.record["artifact_payload"]["synthetic_qualification"]
        interval = synthetic["continuous_interval_controls"]
        self.assertEqual(interval["exact_zero_classification"], "exact_zero")
        self.assertFalse(
            interval["tiny_nonconvergent_complete_admission_passed"]
        )
        self.assertEqual(
            interval["tiny_nonconvergent_classification"],
            "resolved_order_failure",
        )
        self.assertTrue(interval["large_convergent_complete_admission_passed"])
        self.assertGreater(
            interval["large_convergent_minimum_retained_debit"], 9.0
        )
        self.assertTrue(interval["dense_audit_enclosed"])

        failure = synthetic["failure_and_atomicity_controls"]
        self.assertEqual(failure["late_non_source_stop_path"], "fine_2")
        self.assertTrue(
            failure["late_non_source_stop_preserved_all_real_boundaries"]
        )
        self.assertTrue(failure["source_only_retry_owned_by_PROTO7"])
        self.assertTrue(failure["tracer_boundary_drift_rejected_before_commit"])
        self.assertTrue(
            failure[
                "late_tracer_commit_exception_rolled_back_all_mutable_ledgers"
            ]
        )

        retry = synthetic["retry_and_checkpoint_controls"]
        self.assertEqual(retry["first_retry_channel_decision_count"], 18)
        self.assertFalse(retry["first_retry_physical_classification"])
        self.assertTrue(retry["failed_sink_preserved_input_ledger"])
        self.assertTrue(retry["accepted_step_reset_only_current_retry_counter"])
        self.assertTrue(retry["checkpoint_round_trip_passed"])
        self.assertTrue(retry["checkpoint_debit_mutation_rejected"])
        self.assertTrue(retry["checkpoint_claim_promotion_rejected"])
        self.assertEqual(
            retry["retry_exhaustion_reason"], "maximum_temporal_retries"
        )
        self.assertFalse(retry["retry_exhaustion_is_physical_classification"])

    def test_protocol_and_claim_mutations_fail_closed(self) -> None:
        controls = self.record["artifact_payload"]["mutation_controls"]
        self.assertTrue(controls)
        self.assertEqual(set(controls.values()), {True})
        config = imp2.load_config()
        for mutation in (
            lambda value: value["continuous_interval"].__setitem__(
                "exact_squared_admission", "7*U12^2<=L01^2"
            ),
            lambda value: value["claims"].__setitem__(
                "production_trajectory_authorized", True
            ),
            lambda value: value["claims"].__setitem__("PROTO14_frozen", True),
            lambda value: value["claims"].__setitem__("GR0_case_eligible", True),
            lambda value: value["claims"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
        ):
            attacked = deepcopy(config)
            mutation(attacked)
            with self.assertRaises(ValueError):
                imp2.validate_config_data(attacked)

    def test_noncanonical_or_promoted_result_is_rejected(self) -> None:
        promoted = deepcopy(self.record)
        promoted["gate_status"]["production_trajectory_authorized"] = True
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "promoted.json"
            path.write_text(
                json.dumps(promoted, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                imp2.verify_canonical(output_path=path)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "noncanonical.json"
            path.write_text(json.dumps(self.record), encoding="utf-8")
            with self.assertRaises(ValueError):
                imp2.verify_canonical(output_path=path)


if __name__ == "__main__":
    unittest.main()
