from copy import deepcopy
import json
import unittest

from scripts import reproduce_fgc_tdg6_frz1 as reproduce


class TDG6FRZ1ReproductionTests(unittest.TestCase):
    def test_canonical_record_matches_tracked_result(self) -> None:
        record = reproduce.record()
        tracked = json.loads(reproduce.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(record, tracked)
        self.assertEqual(record["artifact_id"], "FGC-1-TDG6-FRZ1")
        self.assertEqual(
            record["classification"],
            "prospective_three_level_temporal_admission_threshold_freeze",
        )
        boundary = record["artifact_payload"]["claim_boundary"]
        self.assertTrue(boundary["runtime_thresholds_frozen"])
        self.assertTrue(boundary["replacement_temporal_admission_defined"])
        self.assertFalse(boundary["production_compositor_implementation_authorized"])
        self.assertFalse(boundary["PROTO14_frozen"])
        self.assertFalse(boundary["GR0_case_eligible"])
        self.assertFalse(boundary["candidate_branch_opened"])

    def test_config_mutations_fail_closed(self) -> None:
        config = reproduce._load_config(reproduce.DEFAULT_CONFIG)
        minimum = deepcopy(config)
        minimum["thresholds"]["minimum_observed_temporal_refinement_order"] = "149/100"
        with self.assertRaises(ValueError):
            reproduce._validate_config(minimum)

        absolute = deepcopy(config)
        absolute["thresholds"]["absolute_state_error_threshold_defined"] = True
        with self.assertRaises(ValueError):
            reproduce._validate_config(absolute)

        promoted = deepcopy(config)
        promoted["claims"]["FGCQR_holdout_execution_authorized"] = True
        with self.assertRaises(ValueError):
            reproduce._validate_config(promoted)


if __name__ == "__main__":
    unittest.main()
