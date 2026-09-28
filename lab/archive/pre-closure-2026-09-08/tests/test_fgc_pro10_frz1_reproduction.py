from __future__ import annotations

import copy
from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(REPOSITORY))

from scripts import reproduce_fgc_pro10_frz1 as reproduction  # noqa: E402
from scripts import reproduce_fgc_hlt8_mon8 as successor  # noqa: E402


class FGCPro10Frz1ReproductionTests(unittest.TestCase):
    def test_canonical_freeze_is_proved_from_immutable_successor(self) -> None:
        # PRO10's empty-namespace observation was temporal pre-launch evidence.
        # Once the authorized campaign has populated that namespace, calling
        # PRO10.record() again must fail rather than rewrite history.  HLT8 is
        # the post-campaign verifier: it reads the canonical PRO10 blob from
        # its immutable checkpoint and independently hash-matches the current
        # artifact and every implementation premise inherited by the freeze.
        lineage = successor._immutable_lineage(successor.load_config())
        result = reproduction.load_canonical_result()
        self.assertEqual(lineage["freeze"]["artifact_id"], "FGC-1-PRO10-FRZ1")
        self.assertEqual(
            lineage["freeze"]["result_sha256"],
            successor.EXPECTED_LINEAGE["protocol_freeze_result_sha256"],
        )
        self.assertTrue(lineage["freeze"]["claims_all_false"])
        self.assertEqual(result["artifact_id"], "FGC-1-PRO10-FRZ1")
        self.assertTrue(all(value is False for value in result["gate_status"].values()))

    def test_freeze_binds_design_evidence_without_runtime_or_candidate_claim(self) -> None:
        result = reproduction.load_canonical_result()
        payload = result["artifact_payload"]
        validation = payload["protocol_validation"]
        diagnosis = payload["immutable_lineage"]["diagnosis"]
        revision = payload["spectral_revision"]
        self.assertEqual(validation["point_counts"], [2049, 4097, 8193])
        self.assertEqual(validation["maximum_nested_tail_ratio"], "1/4")
        self.assertTrue(validation["coarse_to_medium_direct_ratio_veto"])
        self.assertEqual(diagnosis["raw_failed_finest_pair_metric_count"], 16)
        self.assertEqual(diagnosis["diagnostically_saturated_metric_count"], 16)
        self.assertFalse(diagnosis["current_campaign_authorized"])
        self.assertFalse(revision["round_trip_is_a_continuum_error_bound"])
        self.assertFalse(
            revision["diagnostic_saturation_is_physical_resolution_evidence"]
        )
        self.assertTrue(all(value is False for value in result["gate_status"].values()))
        self.assertTrue(all(value is False for value in result["nonclaims"].values()))

    def test_config_broadening_fails_closed(self) -> None:
        text = reproduction.DEFAULT_CONFIG.read_text()
        attacked = text.replace(
            "coarse_to_medium_field_and_derivative_ratios_must_remain_direct_vetoes = true",
            "coarse_to_medium_field_and_derivative_ratios_must_remain_direct_vetoes = false",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.toml"
            path.write_text(attacked)
            with self.assertRaises(ValueError):
                reproduction.load_config(path)

    def test_result_promotion_fails_reproduction(self) -> None:
        attacked = copy.deepcopy(reproduction.load_canonical_result())
        attacked["gate_status"]["PROTO10_fresh_GR0_dynamic_calibration_authorized"] = (
            True
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.json"
            path.write_text(reproduction._canonical(attacked))
            with self.assertRaises(ValueError):
                reproduction.verify_canonical(path)


if __name__ == "__main__":
    unittest.main()
