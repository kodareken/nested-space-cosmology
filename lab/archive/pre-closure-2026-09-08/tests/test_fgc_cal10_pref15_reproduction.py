from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts import reproduce_fgc_cal10_pref15 as pref15


class FGCCAL10PREF15ReproductionTests(unittest.TestCase):
    def test_canonical_certificate_verifies(self) -> None:
        record = pref15.verify_canonical()
        self.assertEqual(record["artifact_id"], pref15.ARTIFACT_ID)
        self.assertTrue(record["gate_status"]["PROTO13_campaign_terminated_normally"])
        self.assertTrue(record["gate_status"]["PROTO13_temporal_admission_failed"])
        self.assertFalse(record["gate_status"]["PROTO13_temporal_failure_cause_derived"])
        self.assertFalse(record["gate_status"]["GR0_case_eligible"])
        self.assertFalse(record["gate_status"]["FGCQR_holdout_execution_authorized"])

    def test_tracked_result_is_canonical(self) -> None:
        stored = json.loads(pref15.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(
            pref15.DEFAULT_OUTPUT.read_text(encoding="utf-8"),
            pref15._canonical(stored),
        )

    def test_complete_absence_is_allowed_but_partial_raw_bundle_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = tuple(root / name for name in ("a", "b", "c", "d"))
            self.assertFalse(pref15._raw_bundle_is_complete_or_absent(paths))
            paths[0].write_bytes(b"present")
            with self.assertRaises(ValueError):
                pref15._raw_bundle_is_complete_or_absent(paths)
            for path in paths[1:]:
                path.write_bytes(b"present")
            self.assertTrue(pref15._raw_bundle_is_complete_or_absent(paths))

    def test_claim_promotion_in_config_fails_closed(self) -> None:
        text = pref15.DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            "FGCQR_holdout_execution_authorized = false",
            "FGCQR_holdout_execution_authorized = true",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutated.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "configuration contract"):
                pref15.load_config(path)

    def test_temporal_threshold_change_in_config_fails_closed(self) -> None:
        text = pref15.DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            'maximum_nested_tail_ratio = "1/4"',
            'maximum_nested_tail_ratio = "1/3"',
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutated.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "configuration contract"):
                pref15.load_config(path)


if __name__ == "__main__":
    unittest.main()
