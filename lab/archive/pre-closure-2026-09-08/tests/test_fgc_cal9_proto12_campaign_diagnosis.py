from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal9_proto12_campaign_diagnosis import (  # noqa: E402
    diagnose_proto12_campaign,
)


class FGCCAL9PROTO12CampaignDiagnosisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        root = REPOSITORY / "runs/fgc-2-sf1/proto12/calibration"
        cls.events = [
            json.loads(line)
            for line in (root / "events.jsonl").read_text(encoding="utf-8").splitlines()
        ]
        cls.campaign = json.loads(
            (root / "campaign-result.json").read_text(encoding="utf-8")
        )

    def _mutated_event(self, amplitude: str, event_index: int) -> list[dict]:
        records = list(self.events)
        position = next(
            index
            for index, item in enumerate(records)
            if item["amplitude"] == amplitude
            and item["event_index"] == event_index
        )
        records[position] = deepcopy(records[position])
        return records

    def test_campaign_localizes_one_declared_gate_at_two_amplitudes(self) -> None:
        result = diagnose_proto12_campaign(self.events, self.campaign)
        self.assertEqual(result["event_log_line_count"], 49)
        self.assertEqual(result["serialized_source_retry_count"], 0)
        self.assertTrue(result["some_adaptive_CFL_step_reductions_occurred"])
        self.assertTrue(result["SRC4_arithmetic_wall_cleared_for_both_amplitudes"])
        self.assertFalse(result["candidate_or_mechanism_outcome_read"])
        shared = result["shared_declared_veto"]
        self.assertEqual(shared["method"], "SSPRK3")
        self.assertEqual(shared["constraint_component"], "radial_momentum")
        self.assertEqual(shared["frozen_minimum_order"], 1.5)
        self.assertFalse(shared["common_numerical_cause_derived"])
        self.assertFalse(shared["common_physical_cause_derived"])

        five = result["amplitudes"]["5/2"]
        three = result["amplitudes"]["3"]
        self.assertEqual(five["terminal_event_index"], 24)
        self.assertEqual(three["terminal_event_index"], 23)
        self.assertAlmostEqual(
            five["comparator_terminal_radial_momentum_order"],
            1.494957481780717,
        )
        self.assertAlmostEqual(
            three["comparator_terminal_radial_momentum_order"],
            1.4990947384363291,
        )
        self.assertGreater(
            five["comparator_previous_radial_momentum_order"], 1.5
        )
        self.assertGreater(
            three["comparator_previous_radial_momentum_order"], 1.5
        )
        self.assertTrue(
            all(
                value > 0
                for value in three["terminal_CFL_retry_count_by_member"].values()
            )
        )

    def test_roundoff_enclosure_is_part_of_the_recomputed_order(self) -> None:
        result = diagnose_proto12_campaign(self.events, self.campaign)
        terminal = result["amplitudes"]["3"]
        raw = terminal["terminal_comparator_radial_momentum_raw_owned_norms"]
        effective = terminal[
            "terminal_comparator_radial_momentum_effective_norms_for_order_only"
        ]
        self.assertNotEqual(raw, effective)
        self.assertEqual(
            terminal["comparator_terminal_radial_momentum_order"],
            terminal["terminal_comparator_component_orders"]["radial_momentum"],
        )

    def test_source_retry_relabel_fails_closed(self) -> None:
        mutated = self._mutated_event("3", 23)
        target = next(
            item
            for item in mutated
            if item["amplitude"] == "3" and item["event_index"] == 23
        )
        target["assessment"]["member_diagnostics"]["RK4-2049"][
            "source_retry_count"
        ] = 1
        with self.assertRaisesRegex(ValueError, "retry provenance"):
            diagnose_proto12_campaign(mutated, self.campaign)

    def test_terminal_spatial_failure_relabel_fails_closed(self) -> None:
        mutated = self._mutated_event("5/2", 24)
        target = next(
            item
            for item in mutated
            if item["amplitude"] == "5/2" and item["event_index"] == 24
        )
        target["assessment"]["comparator_common_event"][
            "spatial_spectral_admission"
        ]["admission_passed"] = False
        with self.assertRaisesRegex(ValueError, "common-event admission"):
            diagnose_proto12_campaign(mutated, self.campaign)

    def test_serialized_order_mutation_fails_closed(self) -> None:
        mutated = self._mutated_event("3", 23)
        target = next(
            item
            for item in mutated
            if item["amplitude"] == "3" and item["event_index"] == 23
        )
        target["assessment"]["comparator_common_event"]["constraint_admission"][
            "component_finest_pair_orders"
        ]["radial_momentum"] = 1.5
        with self.assertRaisesRegex(ValueError, "serialized constraint order"):
            diagnose_proto12_campaign(mutated, self.campaign)

    def test_lowering_the_frozen_threshold_does_not_reclassify_the_run(self) -> None:
        with self.assertRaisesRegex(ValueError, "constraint-order decision"):
            diagnose_proto12_campaign(
                self.events,
                self.campaign,
                minimum_constraint_order=1.49,
            )


if __name__ == "__main__":
    unittest.main()
