from __future__ import annotations

from copy import deepcopy
import unittest

from scripts import reproduce_fgc_pro14_frz1 as pro14


class PRO14FreezeReproductionTests(unittest.TestCase):
    def test_canonical_record_reproduces(self) -> None:
        record = pro14.load_canonical_result()
        pro14.verify_canonical()
        self.assertEqual(record, pro14.record())
        self.assertEqual(record["artifact_id"], pro14.ARTIFACT_ID)
        self.assertTrue(record["gate_status"]["PROTO14_frozen"])
        self.assertFalse(
            record["gate_status"]["PROTO14_fresh_GR0_dynamic_calibration_authorized"]
        )
        self.assertFalse(record["gate_status"]["FGCQR_holdout_execution_authorized"])

    def test_exact_restart_manifest_is_preserved(self) -> None:
        members = pro14.load_canonical_result()["artifact_payload"]["restart_members"]
        self.assertEqual(
            [item["key"] for item in members],
            [
                "RK4-2049",
                "RK4-4097",
                "RK4-8193",
                "SSPRK3-4097",
                "SSPRK3-8193",
                "SSPRK3-16385",
            ],
        )
        self.assertTrue(all(item["coordinate_time"] == 23 / 16 for item in members))
        self.assertTrue(all(item["event_sample_count"] == 24 for item in members))
        self.assertTrue(all(item["tracer_count"] == 48 for item in members))

    def test_historical_stop_and_replacement_boundary_are_both_public(self) -> None:
        payload = pro14.load_canonical_result()["artifact_payload"]
        evidence = payload["evidence_boundary"]
        epistemic = payload["epistemic_boundary"]
        self.assertTrue(evidence["historical_PROTO13_temporal_stop_preserved"])
        self.assertFalse(evidence["historical_terminal_checkpoint_may_resume"])
        self.assertFalse(evidence["historical_stop_is_physical_obstruction"])
        self.assertTrue(evidence["production_compositor_implemented"])
        self.assertTrue(evidence["production_compositor_synthetic_qualification_passed"])
        self.assertFalse(epistemic["restart_state_has_pre_restart_temporal_error_bound"])
        self.assertFalse(epistemic["restart_state_may_seed_candidate_or_physical_claim"])
        self.assertTrue(epistemic["sampled_64_history_rule_is_nonveto_diagnostic_only"])
        self.assertFalse(epistemic["TDG6_debit_is_global_PDE_or_Raychaudhuri_error_bound"])

    def test_embedded_mutation_controls_all_pass(self) -> None:
        controls = pro14.load_canonical_result()["artifact_payload"]["mutation_controls"]
        self.assertEqual(
            set(controls),
            {
                "threshold",
                "sampled_diagnostic_veto",
                "restart_error_promotion",
                "terminal_candidate_promotion",
                "claim",
                "restart_state_hash",
                "namespace",
            },
        )
        self.assertTrue(all(controls.values()))

    def test_restart_and_namespace_mutations_differ_from_frozen_contract(self) -> None:
        config = pro14.load_config()
        for name, mutate in {
            "state": lambda value: value["restart_members"][0].__setitem__(
                "state_sha256", "0" * 64
            ),
            "namespace": lambda value: value["namespace_precondition"].__setitem__(
                "calibration_output_root", "runs/fgc-2-sf1/proto13/calibration"
            ),
            "claim": lambda value: value["claims"].__setitem__(
                "GR0_case_eligible", True
            ),
        }.items():
            with self.subTest(name=name):
                attacked = deepcopy(config)
                mutate(attacked)
                self.assertNotEqual(attacked, config)


if __name__ == "__main__":
    unittest.main()
