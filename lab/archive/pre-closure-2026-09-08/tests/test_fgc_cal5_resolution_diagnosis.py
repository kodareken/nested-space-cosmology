from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.evolution.cal5_resolution_diagnosis import (  # noqa: E402
    diagnose_proto8_campaign_events,
)


class Cal5ResolutionDiagnosisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        path = REPOSITORY / "runs/fgc-2-sf1/proto8/calibration/events.jsonl"
        if not path.is_file():
            raise unittest.SkipTest("immutable PROTO8 campaign is unavailable")
        cls.events = [json.loads(line) for line in path.read_text().splitlines()]

    def test_public_campaign_reduces_to_one_level_resolution_test(self) -> None:
        result = diagnose_proto8_campaign_events(self.events)
        self.assertEqual(result["current_point_counts"], [1025, 2049, 4097])
        self.assertEqual(
            result["smallest_prospective_point_counts"], [2049, 4097, 8193]
        )
        self.assertFalse(result["threshold_changes_required"])
        self.assertTrue(result["new_8193_data_required_before_admission"])
        for key in ("amplitude_5over2_stop", "amplitude_3_stop"):
            for method in ("RK4", "SSPRK3"):
                item = result[key]["methods"][method]
                self.assertTrue(item["prospective_2049_coarsest_guard_passed"])
                self.assertTrue(item["conditional_8193_finest_guard_passed"])
                self.assertFalse(item["conditional_projection_is_admission_evidence"])

    def test_amplitude_three_tail_failure_is_not_hidden(self) -> None:
        result = diagnose_proto8_campaign_events(self.events)
        for method in ("RK4", "SSPRK3"):
            failures = result["amplitude_3_stop"]["methods"][method][
                "current_nested_tail_failures"
            ]
            self.assertEqual(len(failures), 1)
            self.assertEqual(failures[0]["family"], "derivative")
            self.assertEqual(failures[0]["field"], "phi_over_Lambda")
            self.assertEqual(failures[0]["pair"], "2049_to_4097")
            self.assertGreater(failures[0]["ratio"], failures[0]["maximum"])

    def test_non_source_retry_cannot_be_relabelled(self) -> None:
        attacked = copy.deepcopy(self.events)
        retry = next(
            item
            for item in attacked
            if item.get("event_type") == "rejected_unaccepted_source_only_proposal"
        )
        retry["non_source_failure_vetoed_retry"] = True
        with self.assertRaises(ValueError):
            diagnose_proto8_campaign_events(attacked)

    def test_relaxing_a_public_tail_does_not_change_the_frozen_threshold(self) -> None:
        attacked = copy.deepcopy(self.events)
        stop = next(
            item
            for item in attacked
            if item.get("event_type") == "common_event"
            and item.get("amplitude") == "3"
            and item.get("event_index") == 1
        )
        stop["assessment"]["primary_common_event"]["spatial_spectral_admission"][
            "derivative_power_tail_ratios"
        ]["phi_over_Lambda"][1] = 0.2
        with self.assertRaises(ValueError):
            diagnose_proto8_campaign_events(attacked)


if __name__ == "__main__":
    unittest.main()
