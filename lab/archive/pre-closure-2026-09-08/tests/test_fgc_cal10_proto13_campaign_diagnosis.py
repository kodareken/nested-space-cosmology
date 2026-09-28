from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal10_proto13_campaign_diagnosis import (  # noqa: E402
    diagnose_proto13_campaign,
)


RAW_ROOT = REPOSITORY / "runs/fgc-2-sf1/proto13/calibration"
RAW_AVAILABLE = all(
    (RAW_ROOT / name).is_file()
    for name in ("events.jsonl", "campaign-result.json")
)


@unittest.skipUnless(RAW_AVAILABLE, "ignored PROTO13 raw campaign is absent")
class FGCCAL10PROTO13CampaignDiagnosisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.events = [
            json.loads(line)
            for line in (RAW_ROOT / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        cls.campaign = json.loads(
            (RAW_ROOT / "campaign-result.json").read_text(encoding="utf-8")
        )

    def test_campaign_localizes_only_the_first_temporal_gate(self) -> None:
        result = diagnose_proto13_campaign(self.events, self.campaign)
        self.assertEqual(result["event_log_line_count"], 41)
        self.assertEqual(result["terminal_event_index"], 63)
        self.assertTrue(
            result["terminal_common_spatial_and_constraint_admissions_passed"]
        )
        self.assertFalse(result["terminal_trapped_sign_passed"])
        self.assertTrue(result["both_temporal_method_admissions_failed"])
        self.assertTrue(result["both_temporal_individual_budget_layers_failed"])
        self.assertTrue(result["both_temporal_nested_ratio_layers_failed"])
        self.assertFalse(result["temporal_failure_cause_derived"])
        self.assertFalse(result["GR0_case_eligible"])
        self.assertFalse(result["candidate_or_mechanism_outcome_read"])

    def test_exact_failed_ratio_counts_are_preserved(self) -> None:
        result = diagnose_proto13_campaign(self.events, self.campaign)
        rk4 = result["temporal_methods"]["RK4"]["ratio_summaries"]
        ssprk3 = result["temporal_methods"]["SSPRK3"]["ratio_summaries"]
        self.assertEqual(
            rk4["field_power_tail_ratios"]["failed_comparison_count"], 510
        )
        self.assertEqual(
            rk4["derivative_power_tail_ratios"]["failed_comparison_count"],
            511,
        )
        self.assertEqual(
            ssprk3["field_power_tail_ratios"]["failed_comparison_count"], 453
        )
        self.assertEqual(
            ssprk3["derivative_power_tail_ratios"][
                "failed_comparison_count"
            ],
            453,
        )

    def test_post_outcome_threshold_change_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "threshold changed"):
            diagnose_proto13_campaign(
                self.events,
                self.campaign,
                maximum_nested_tail_ratio=1.0 / 3.0,
            )

    def test_temporal_pass_relabel_fails_closed(self) -> None:
        mutated = deepcopy(self.events)
        mutated[-1]["assessment"]["temporal_spectral_admission"][
            "admission_passed"
        ] = True
        with self.assertRaisesRegex(ValueError, "terminal temporal"):
            diagnose_proto13_campaign(mutated, self.campaign)

    def test_source_retry_relabel_fails_closed(self) -> None:
        mutated = deepcopy(self.campaign)
        mutated["amplitude_records"][0]["member_source_retry_counts"][
            "RK4-2049"
        ] = 1
        with self.assertRaisesRegex(ValueError, "retry"):
            diagnose_proto13_campaign(self.events, mutated)

    def test_candidate_outcome_relabel_fails_closed(self) -> None:
        mutated = deepcopy(self.campaign)
        mutated["FGCQR_outcome_read"] = True
        with self.assertRaisesRegex(ValueError, "campaign result"):
            diagnose_proto13_campaign(self.events, mutated)


if __name__ == "__main__":
    unittest.main()
