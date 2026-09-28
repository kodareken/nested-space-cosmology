from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import reproduce_fgc_tdg7_imp3 as imp3  # noqa: E402


class TDG7IMP3ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = imp3.verify_canonical()

    def test_canonical_build_is_the_tracked_hash_bound_result(self) -> None:
        self.assertEqual(imp3.build(imp3.load_config()), self.record)
        self.assertEqual(self.record["artifact_id"], "FGC-1-TDG7-IMP3")
        self.assertEqual(
            self.record["classification"],
            "authorized_TDG7_stage_safe_runtime_repair_implemented_and_"
            "synthetically_qualified_without_history_execution",
        )
        status = self.record["gate_status"]
        for name in (
            "CAL11_terminal_invalid_runtime_result_preserved",
            "TDG7_independent_binder_completed",
            "TDG7_runtime_repair_implementation_authorized",
            "TDG7_runtime_repair_implemented",
            "TDG7_runtime_repair_synthetic_qualification_passed",
            "one_shared_stage_safe_plan_owns_all_paths",
            "typed_lattice_stops_precede_work",
            "TDG6_threshold_and_all_of_admission_unchanged",
            "stateless_adapter_cannot_detect_discarded_external_retry",
            "future_protocol_cursor_required_before_calibration",
            "replacement_protocol_freeze_design_authorized",
        ):
            self.assertTrue(status[name], name)
        self.assertTrue(all(status[name] is False for name in imp3.CLAIM_FALSE))
        payload = self.record["artifact_payload"]
        self.assertTrue(payload["synthetic_qualification"]["synthetic_states_only"])
        self.assertFalse(payload["synthetic_qualification"]["CAL11_history_arrays_loaded"])
        self.assertFalse(payload["synthetic_qualification"]["PROTO14_history_arrays_loaded"])
        self.assertFalse(payload["synthetic_qualification"]["runs_namespace_created"])
        discriminator = payload["synthetic_qualification"][
            "protocol_discriminator_limitation"
        ]
        self.assertFalse(discriminator["fresh_vs_retry_protocol_discriminator_implemented"])
        self.assertTrue(
            discriminator["stateless_adapter_cannot_detect_discarded_external_retry"]
        )
        self.assertTrue(discriminator["future_protocol_cursor_required_before_calibration"])

    def test_both_methods_use_one_shared_plan_and_commit_only_fine(self) -> None:
        methods = self.record["artifact_payload"]["synthetic_qualification"][
            "method_controls"
        ]
        expected = {"RK4": (35, 20), "SSPRK3": (28, 16)}
        for method, (total, fine) in expected.items():
            with self.subTest(method=method):
                evidence = methods[method]
                self.assertEqual(evidence["shadow_proposal_count"], 7)
                self.assertEqual(evidence["total_shadow_stage_records"], total)
                self.assertEqual(evidence["fine_stage_records"], fine)
                self.assertTrue(evidence["one_shared_plan_owns_all_paths"])
                self.assertTrue(evidence["complete_admission_passed"])
                self.assertTrue(evidence["real_boundary_unchanged_during_prepare"])
                self.assertFalse(evidence["outer_and_medium_committed"])
                self.assertTrue(evidence["fine_quarter_steps_committed_atomically"])
                self.assertEqual(
                    evidence["committed_state_sha256"], evidence["fine_shadow_state_sha256"]
                )
                self.assertEqual(
                    evidence["committed_tracer_snapshot"],
                    evidence["fine_shadow_tracer_snapshot"],
                )
                plan = evidence["plan"]
                self.assertEqual(len(plan["boundaries_hex"]), 5)
                self.assertEqual(len(plan["one_full_stage_times_hex"]), 1)
                self.assertEqual(len(plan["two_half_stage_times_hex"]), 2)
                self.assertEqual(len(plan["four_quarter_stage_times_hex"]), 4)
                self.assertEqual(
                    sum(
                        len(stages)
                        for paths in evidence["actual_shadow_stage_records"].values()
                        for stages in paths
                    ),
                    total,
                )

    def test_stops_intervals_late_shadow_rollback_and_retry_evidence(self) -> None:
        synthetic = self.record["artifact_payload"]["synthetic_qualification"]
        stops = synthetic["typed_lattice_stop_controls"]
        self.assertEqual(len(stops), 11)
        for case, evidence in stops.items():
            self.assertEqual(evidence["case"], case)
            self.assertEqual(evidence["typed_stop"], "TDG7CoordinateLatticeStop")
            self.assertEqual(evidence["hook_count"], 0)
        frozen_minimum = synthetic["frozen_minimum_width_control"]
        self.assertTrue(frozen_minimum["caller_minimum_width_override_rejected"])
        self.assertEqual(frozen_minimum["runtime_hook_count"], 0)

        intervals = synthetic["TDG6_all_of_interval_controls"]
        self.assertEqual(intervals["channel_count"], 18)
        self.assertEqual(intervals["exact_zero_classification"], "exact_zero")
        self.assertTrue(intervals["exact_zero_complete_admission_passed"])
        self.assertEqual(
            intervals["one_channel_veto_classification"], "resolved_order_failure"
        )
        self.assertFalse(intervals["one_channel_veto_complete_admission_passed"])
        self.assertTrue(intervals["large_convergent_complete_admission_passed"])
        self.assertGreater(intervals["large_convergent_minimum_retained_debit"], 9.0)
        self.assertTrue(intervals["complete_admission_is_all_of"])

        failures = synthetic["late_shadow_and_tracer_controls"]
        self.assertEqual(failures["non_source"]["path"], "fine_2")
        self.assertFalse(failures["non_source"]["source_retry_owned_by_PROTO7"])
        self.assertEqual(failures["source_only"]["path"], "fine_2")
        self.assertTrue(failures["source_only"]["source_retry_owned_by_PROTO7"])
        self.assertTrue(failures["precommit_tracer_drift_rejected"])
        self.assertTrue(failures["late_tracer_commit_rollback_passed"])

        retry = synthetic["retry_replan_controls"]
        self.assertEqual(retry["first_retry_channel_decision_count"], 18)
        self.assertFalse(retry["retry_is_physical_classification"])
        self.assertTrue(retry["replan_strictly_reduced"])
        self.assertTrue(
            retry["retry_successor_binds_predecessor_plan_and_exact_updated_ledger"]
        )
        self.assertTrue(retry["canonical_durable_TDG6_evidence_retained"])
        self.assertTrue(retry["active_retry_ledger_rejected_by_fresh_entrypoint"])
        self.assertTrue(retry["retry_successor_resume_uses_exact_updated_ledger"])
        self.assertTrue(retry["over_limit_retry_lineage_rejected"])
        self.assertTrue(retry["mismatched_method_retry_lineage_rejected"])
        self.assertTrue(retry["forged_canonical_retry_evidence_rejected"])
        self.assertTrue(retry["full_ledger_prefix_retained"])
        self.assertTrue(retry["last_accepted_retry_count_retained"])
        self.assertTrue(retry["forged_full_ledger_prefix_rejected"])
        self.assertTrue(retry["two_minimum_to_minimum_successor"])
        self.assertEqual(
            retry["TDG6_minimum_macro_step_exhaustion_reason"], "minimum_macro_step"
        )
        self.assertFalse(retry["TDG6_minimum_macro_step_exhaustion_is_physical"])
        self.assertTrue(
            retry["TDG7_lattice_exhaustion_is_defensive_forged_wrapper_only"]
        )
        self.assertEqual(retry["defensive_lattice_stop"], "TDG7RetryLatticeExhausted")
        self.assertTrue(retry["defensive_lattice_stop_retains_TDG6_evidence"])
        self.assertFalse(retry["defensive_lattice_stop_is_physical_classification"])

    def test_lineage_hash_contract_and_nonpromotion_mutations_fail_closed(self) -> None:
        config = imp3.load_config()
        attacks = []
        value = deepcopy(config)
        value["immutable_lineage"]["predecessor_commit"] = "0" * 40
        attacks.append(value)
        value = deepcopy(config)
        value["immutable_lineage"]["PREF23_result_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(config)
        value["immutable_lineage"]["core_runtime_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(config)
        value["stage_safe_runtime_contract"][
            "one_immutable_shared_TDG7_plan_owns_one_two_four_paths"
        ] = False
        attacks.append(value)
        for claim in (
            "PROTO14_runtime_mutation_authorized",
            "fresh_replacement_protocol_frozen",
            "new_protocol_authorized",
            "new_output_namespace_authorized",
            "fresh_GR0_calibration_authorized",
            "FGCQR_holdout_execution_authorized",
            "physical_transition_claim_authorized",
        ):
            value = deepcopy(config)
            value["claims"][claim] = True
            attacks.append(value)
        for attacked in attacks:
            with self.subTest(attacked=attacked):
                with self.assertRaises(ValueError):
                    imp3.validate_config_data(attacked)

    def test_noncanonical_result_hash_mutation_and_claim_promotion_are_rejected(self) -> None:
        changed_hash = deepcopy(self.record)
        changed_hash["implementation_sha256"][
            "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_runtime.py"
        ] = "0" * 64
        promoted = deepcopy(self.record)
        promoted["gate_status"]["FGCQR_holdout_execution_authorized"] = True
        malformed = json.dumps(self.record).encode("utf-8")
        for name, value in (("hash.json", changed_hash), ("promoted.json", promoted)):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / name
                path.write_bytes(imp3.canonical(value))
                with self.assertRaises(ValueError):
                    imp3.verify_canonical(output_path=path)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "noncanonical.json"
            path.write_bytes(malformed)
            with self.assertRaises(ValueError):
                imp3.verify_canonical(output_path=path)


if __name__ == "__main__":
    unittest.main()
