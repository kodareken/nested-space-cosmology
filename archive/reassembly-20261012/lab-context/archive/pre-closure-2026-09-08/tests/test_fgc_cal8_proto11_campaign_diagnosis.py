from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal8_proto11_campaign_diagnosis import (  # noqa: E402
    diagnose_proto11_campaign,
)


class FGCCAL8PROTO11CampaignDiagnosisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        root = REPOSITORY / "runs/fgc-2-sf1/proto11/calibration"
        cls.events = [
            json.loads(line)
            for line in (root / "events.jsonl").read_text(encoding="utf-8").splitlines()
        ]
        cls.campaign = json.loads(
            (root / "campaign-result.json").read_text(encoding="utf-8")
        )

    def test_campaign_localizes_two_independent_obstructions(self) -> None:
        result = diagnose_proto11_campaign(self.events, self.campaign)
        self.assertEqual(result["event_log_line_count"], 168)
        self.assertEqual(result["source_only_rejection_count"], 164)
        self.assertTrue(result["both_amplitudes_reached_one_evolved_common_event"])
        source = result["amplitude_five_halves_source_obstruction"]
        self.assertEqual(source["accepted_time_group_count"], 16)
        self.assertTrue(source["all_retry_sequences_are_exact_binary_halvings"])
        self.assertTrue(source["cross_member_rollback_verified"])
        spectral = result["amplitude_three_direct_spectral_obstruction"]
        self.assertEqual(
            spectral["sole_direct_coarse_to_medium_veto"],
            "phi_over_Lambda derivative tail",
        )
        for method in ("RK4", "SSPRK3"):
            ratios = spectral["methods"][method]["phi_derivative_tail_ratios"]
            self.assertGreaterEqual(ratios[0], 0.25)
            self.assertLess(ratios[1], 0.25)

    def test_non_source_relabel_fails_closed(self) -> None:
        mutated = deepcopy(self.events)
        rejected = next(
            item
            for item in mutated
            if item["event_type"] == "rejected_unaccepted_source_only_proposal"
        )
        rejected["complete_failure_set"].append("kinetic_condition_limit")
        with self.assertRaisesRegex(ValueError, "source-only rollback contract"):
            diagnose_proto11_campaign(mutated, self.campaign)

    def test_rollback_mutation_fails_closed(self) -> None:
        mutated = deepcopy(self.events)
        rejected = next(
            item
            for item in mutated
            if item["event_type"] == "rejected_unaccepted_source_only_proposal"
        )
        rejected["fields_bitwise_preserved"] = False
        with self.assertRaisesRegex(ValueError, "source-only rollback contract"):
            diagnose_proto11_campaign(mutated, self.campaign)

    def test_raw_threshold_mutation_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "failed source evaluation"):
            diagnose_proto11_campaign(
                self.events,
                self.campaign,
                raw_source_limit=2.0e-12,
            )

    def test_direct_coarse_veto_mutation_fails_closed(self) -> None:
        mutated = deepcopy(self.events)
        event = next(
            item
            for item in mutated
            if item.get("event_type") == "common_event"
            and item.get("amplitude") == "3"
            and item.get("event_index") == 1
        )
        event["assessment"]["primary_common_event"]["spatial_spectral_admission"][
            "coarse_to_medium_direct_passed_by_derivative"
        ]["phi_over_Lambda"] = True
        with self.assertRaisesRegex(ValueError, "direct spectral veto"):
            diagnose_proto11_campaign(mutated, self.campaign)


if __name__ == "__main__":
    unittest.main()
