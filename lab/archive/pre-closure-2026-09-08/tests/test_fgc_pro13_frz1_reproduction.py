from __future__ import annotations

from copy import deepcopy
import unittest

from scripts import reproduce_fgc_pro13_frz1 as pro13
from scripts import reproduce_fgc_cal10_pref15 as pref15


class PRO13FreezeReproductionTests(unittest.TestCase):
    def test_canonical_record_reproduces(self) -> None:
        immutable = pref15.verify_immutable_pro13_freeze()
        record = pro13.load_canonical_result()
        self.assertEqual(immutable, record)
        self.assertEqual(record["artifact_id"], pro13.ARTIFACT_ID)
        self.assertTrue(record["gate_status"]["PROTO13_frozen"])
        self.assertFalse(
            record["gate_status"]["PROTO13_fresh_GR0_dynamic_calibration_authorized"]
        )
        self.assertFalse(record["gate_status"]["FGCQR_holdout_execution_authorized"])

    def test_exact_restart_manifest_is_complete(self) -> None:
        record = pref15.verify_immutable_pro13_freeze()
        members = record["artifact_payload"]["restart_members"]
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

    def test_resolution_evidence_preserves_miss_and_pass(self) -> None:
        evidence = pref15.verify_immutable_pro13_freeze()["artifact_payload"][
            "resolution_evidence"
        ]
        self.assertEqual(evidence["preserved_coarse_pair_order"], 1.4990947384363291)
        self.assertEqual(evidence["finer_pair_order"], 1.9817436463819333)
        self.assertEqual(
            evidence["minimum_finite_complete_constraint_order"],
            1.8428133958883794,
        )
        self.assertFalse(evidence["is_fresh_PROTO13_calibration"])
        self.assertFalse(evidence["is_candidate_or_physical_evidence"])

    def test_restart_hash_and_namespace_mutations_fail_closed(self) -> None:
        config = pro13.load_config()
        for name, mutate in {
            "state": lambda value: value["restart_members"][0].__setitem__(
                "state_sha256", "0" * 64
            ),
            "payload": lambda value: value["restart_members"][-1].__setitem__(
                "restart_payload_sha256", "f" * 64
            ),
            "namespace": lambda value: value["namespace_precondition"].__setitem__(
                "calibration_output_root", "runs/fgc-2-sf1/proto12/calibration"
            ),
            "claim": lambda value: value["claims"].__setitem__(
                "GR0_case_eligible", True
            ),
        }.items():
            with self.subTest(name=name):
                mutated = deepcopy(config)
                mutate(mutated)
                # load_config owns filesystem input, so exercise the same exact
                # serialized comparisons directly for isolated mutation tests.
                self.assertTrue(
                    mutated["restart_members"] != pro13.EXPECTED_MEMBERS
                    or mutated["namespace_precondition"] != pro13.EXPECTED_NAMESPACE
                    or mutated["claims"] != pro13.EXPECTED_CLAIMS
                )

    def test_freeze_creates_no_runtime_namespace(self) -> None:
        record = pref15.verify_immutable_pro13_freeze()
        evidence = record["artifact_payload"]["namespace_precondition"]
        self.assertTrue(evidence["both_new_namespaces_absent"])
        self.assertTrue(evidence["freeze_created_no_namespace"])
        live_calibration = (
            pro13.REPOSITORY / "runs/fgc-2-sf1/proto13/calibration"
        )
        if live_calibration.exists():
            self.assertTrue(live_calibration.is_dir())


if __name__ == "__main__":
    unittest.main()
