from __future__ import annotations

import copy
import json
from pathlib import Path
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.rsp1_campaign_diagnosis import (  # noqa: E402
    EXPECTED_MEMBER_KEYS,
    diagnose_rsp1_campaign,
)


EVENTS = REPOSITORY / "runs/fgc-2-sf1/rsp1/amplitude3-spectrum/events.jsonl"
RESULT = REPOSITORY / "runs/fgc-2-sf1/rsp1/amplitude3-spectrum/study-result.json"
STUDY_ID = "0de794592ecafb7b93aff5b316d5b393a4007008c877535baf826a3c8d45f53f"


def _load_events() -> list[dict]:
    return [json.loads(line) for line in EVENTS.read_text(encoding="utf-8").splitlines()]


class FGCRSP1CampaignDiagnosisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.events = _load_events()
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))
        cls.diagnosis = diagnose_rsp1_campaign(
            cls.events,
            cls.result,
            study_id=STUDY_ID,
        )

    def test_target_pass_and_generic_raw_failure_are_both_public(self) -> None:
        self.assertTrue(
            self.diagnosis[
                "amplitude_three_spectral_veto_cleared_for_successor_design"
            ]
        )
        self.assertFalse(
            self.diagnosis["generic_all_field_raw_spectral_admission_passed"]
        )
        self.assertFalse(self.diagnosis["GR0_case_eligible"])
        for method in ("RK4", "SSPRK3"):
            initial = self.diagnosis["initial_premise"]["methods"][method]
            endpoint = self.diagnosis["evolved_endpoint"]["methods"][method]
            self.assertEqual(initial["target_direct_pair_passed"], [False, False])
            self.assertEqual(endpoint["target_direct_pair_passed"], [True, True])
            self.assertFalse(
                endpoint["generic_raw_all_field_spectral_admission_passed"]
            )

    def test_all_six_members_reach_endpoint_without_source_retry(self) -> None:
        self.assertEqual(set(self.diagnosis["members"]), set(EXPECTED_MEMBER_KEYS))
        self.assertEqual(self.diagnosis["total_source_retry_count"], 0)
        self.assertTrue(self.diagnosis["all_six_members_reached_endpoint"])
        self.assertGreater(
            self.diagnosis["minimum_remaining_boundary_margin"], 80.0
        )
        for item in self.diagnosis["members"].values():
            self.assertEqual(item["source_retry_count"], 0)
            self.assertGreater(item["step_index"], 0)

    def test_ratio_reclassification_and_claim_promotion_fail_closed(self) -> None:
        events = copy.deepcopy(self.events)
        result = copy.deepcopy(self.result)
        events[-1]["assessment"]["primary_common_event"][
            "target_derivative_tail_ratios"
        ][0] = 0.25
        result["endpoint_assessment"] = events[-1]["assessment"]
        with self.assertRaisesRegex(ValueError, "target pair classification"):
            diagnose_rsp1_campaign(events, result, study_id=STUDY_ID)

        result = copy.deepcopy(self.result)
        result["mechanism_question_answered"] = True
        with self.assertRaisesRegex(ValueError, "terminal nonclaim"):
            diagnose_rsp1_campaign(self.events, result, study_id=STUDY_ID)

    def test_missing_member_and_hidden_generic_failure_fail_closed(self) -> None:
        events = copy.deepcopy(self.events)
        result = copy.deepcopy(self.result)
        del events[-1]["assessment"]["member_diagnostics"]["RK4-4097"]
        result["endpoint_assessment"] = events[-1]["assessment"]
        with self.assertRaisesRegex(ValueError, "member set"):
            diagnose_rsp1_campaign(events, result, study_id=STUDY_ID)

        events = copy.deepcopy(self.events)
        result = copy.deepcopy(self.result)
        events[-1]["assessment"]["comparator_common_event"][
            "raw_spatial_spectral_admission"
        ]["admission_passed"] = True
        result["endpoint_assessment"] = events[-1]["assessment"]
        with self.assertRaisesRegex(ValueError, "generic raw all-field"):
            diagnose_rsp1_campaign(events, result, study_id=STUDY_ID)


if __name__ == "__main__":
    unittest.main()
