from __future__ import annotations

import copy
from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(REPOSITORY))

from scripts import reproduce_fgc_pro11_frz1 as reproduction  # noqa: E402
from scripts import reproduce_fgc_hlt9_mon9 as successor  # noqa: E402


class FGCPro11Frz1ReproductionTests(unittest.TestCase):
    def test_canonical_freeze_is_proved_from_immutable_successor(self) -> None:
        # PRO11's absent-namespace observation is pre-launch evidence. HLT9 is
        # its post-launch-safe verifier: it reads the exact canonical freeze
        # from immutable commit b53c0f1 and hash-matches the current artifact.
        lineage = successor._immutable_lineage(successor.load_config())
        result = reproduction.load_canonical_result()
        self.assertEqual(lineage["freeze"]["artifact_id"], "FGC-1-PRO11-FRZ1")
        self.assertEqual(
            lineage["freeze"]["result_sha256"],
            successor.EXPECTED_LINEAGE["protocol_freeze_result_sha256"],
        )
        self.assertTrue(lineage["freeze"]["claims_all_false"])
        self.assertEqual(result["artifact_id"], "FGC-1-PRO11-FRZ1")
        self.assertTrue(all(value is False for value in result["gate_status"].values()))

    def test_freeze_binds_only_the_numerical_map_and_nonclaims(self) -> None:
        result = reproduction.load_canonical_result()
        payload = result["artifact_payload"]
        validation = payload["protocol_validation"]
        revision = payload["numerical_map_revision"]
        diagnosis = payload["immutable_lineage"]["diagnosis"]
        self.assertEqual(validation["point_counts"], [2049, 4097, 8193])
        self.assertTrue(validation["only_reference_state_differentiation_changes"])
        self.assertEqual(revision["initial_q_rule"], "D_h(u-u_ref)+q_ref")
        self.assertEqual(revision["source_q_r_rule"], "D_h(q-q_ref)")
        self.assertFalse(revision["interior_q_reprojection_added"])
        self.assertFalse(revision["continuum_equation_or_source_row_changed"])
        self.assertEqual(revision["raw_source_residual_maximum"], "1/1000000000000")
        self.assertEqual(diagnosis["campaign_event_count"], 264)
        self.assertEqual(diagnosis["source_only_rejection_count"], 262)
        self.assertTrue(all(value is False for value in result["gate_status"].values()))
        self.assertTrue(all(value is False for value in result["nonclaims"].values()))

    def test_config_broadening_fails_closed(self) -> None:
        text = reproduction.DEFAULT_CONFIG.read_text()
        attacked = text.replace(
            "interior_q_reprojection_damping_and_constraint_cleaning_are_forbidden = true",
            "interior_q_reprojection_damping_and_constraint_cleaning_are_forbidden = false",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.toml"
            path.write_text(attacked)
            with self.assertRaises(ValueError):
                reproduction.load_config(path)

    def test_result_promotion_fails_reproduction(self) -> None:
        attacked = copy.deepcopy(reproduction.load_canonical_result())
        attacked["gate_status"]["PROTO11_fresh_GR0_dynamic_calibration_authorized"] = (
            True
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.json"
            path.write_text(reproduction._canonical(attacked))
            with self.assertRaises(ValueError):
                reproduction.verify_canonical(path)


if __name__ == "__main__":
    unittest.main()
