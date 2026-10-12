from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from scripts import reproduce_fgc_tdg1_pref16 as pref16


class FGCTDG1PREF16ReproductionTests(unittest.TestCase):
    def test_canonical_certificate_verifies(self) -> None:
        record = pref16.verify_canonical()
        self.assertEqual(record["artifact_id"], pref16.ARTIFACT_ID)
        self.assertEqual(
            record["classification"],
            "completed_TDG1_history_diagnosis_TDG2_required",
        )
        self.assertTrue(
            record["gate_status"]["TDG1_terminal_histories_diagnosed"]
        )
        self.assertFalse(
            record["gate_status"][
                "normalized_raw_temporal_tail_ratio_is_valid_general_spatial_convergence_test"
            ]
        )
        self.assertFalse(
            record["gate_status"]["replacement_temporal_admission_defined"]
        )
        self.assertFalse(record["gate_status"]["GR0_case_eligible"])

    def test_compact_findings_are_exact_and_claim_boundary_is_closed(self) -> None:
        record = pref16.verify_canonical()
        payload = record["artifact_payload"]
        compact = payload["compact_findings"]
        self.assertEqual(compact, pref16._expected_observed())
        self.assertEqual(compact["inherited_individual_failure_count"], 544)
        self.assertEqual(
            compact["failed_below_declared_absolute_amplitude_floor"], 149
        )
        self.assertEqual(
            compact["failed_above_declared_absolute_amplitude_floor"], 395
        )
        self.assertTrue(
            payload["conclusions"][
                "TDG2_absolute_tail_discriminator_design_may_begin"
            ]
        )
        self.assertFalse(
            payload["claim_boundary"]["replacement_temporal_admission_defined"]
        )
        self.assertFalse(payload["claim_boundary"]["new_trajectory_authorized"])
        self.assertFalse(payload["claim_boundary"]["candidate_branch_opened"])
        self.assertFalse(payload["claim_boundary"]["physical_question_answered"])

    def test_tracked_result_is_canonical(self) -> None:
        raw = pref16.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, pref16._canonical_bytes(stored))

    def test_complete_absence_is_allowed_but_partial_raw_bundle_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = tuple(root / name for name in ("a", "b", "c", "d"))
            self.assertFalse(pref16._raw_bundle_is_complete_or_absent(paths))
            paths[0].write_bytes(b"present")
            with self.assertRaisesRegex(ValueError, "bundle is partial"):
                pref16._raw_bundle_is_complete_or_absent(paths)
            for path in paths[1:]:
                path.write_bytes(b"present")
            self.assertTrue(pref16._raw_bundle_is_complete_or_absent(paths))

    def test_config_mutations_fail_closed(self) -> None:
        config = pref16.load_config()
        attacks = []
        value = deepcopy(config)
        value["immutable_lineage"]["checkpoint_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(config)
        value["observed_diagnosis"]["inherited_individual_failure_count"] = 543
        attacks.append(value)
        value = deepcopy(config)
        value["conclusions"][
            "replacement_gate_design_may_not_yet_begin"
        ] = False
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["GR0_case_eligible"] = True
        attacks.append(value)
        for attacked in attacks:
            with self.assertRaises(ValueError):
                pref16.validate_config_data(attacked)

    def test_all_mutation_controls_pass(self) -> None:
        record = pref16.verify_canonical()
        controls = record["artifact_payload"]["mutation_controls"]
        self.assertTrue(controls)
        self.assertEqual(set(controls.values()), {True})


if __name__ == "__main__":
    unittest.main()
