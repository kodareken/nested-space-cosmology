from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts import reproduce_fgc_rsp2_pref14 as pref14


class RSP2PREF14ReproductionTests(unittest.TestCase):
    def test_canonical_result_reproduces_from_raw_bundle(self) -> None:
        record = pref14.verify_canonical()
        self.assertEqual(record["artifact_id"], pref14.ARTIFACT_ID)
        self.assertTrue(record["gate_status"]["RSP2_target_order_cleared"])
        self.assertFalse(record["gate_status"]["PROTO13_frozen"])

    def test_complete_absence_is_allowed_but_partial_raw_bundle_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = tuple(root / name for name in ("a", "b", "c", "d"))
            self.assertFalse(pref14._raw_bundle_is_complete_or_absent(paths))
            paths[0].write_bytes(b"present")
            with self.assertRaises(ValueError):
                pref14._raw_bundle_is_complete_or_absent(paths)
            for path in paths[1:]:
                path.write_bytes(b"present")
            self.assertTrue(pref14._raw_bundle_is_complete_or_absent(paths))


if __name__ == "__main__":
    unittest.main()
