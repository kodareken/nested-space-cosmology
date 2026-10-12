from __future__ import annotations

from copy import deepcopy
import unittest

from scripts import reproduce_fgc_rsp2_pref14 as pref14


class RSP2PREF14DiagnosisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = pref14.reproduce()

    def test_raw_checkpoint_reduction_clears_target_and_complete_admission(self) -> None:
        reduction = self.record["artifact_payload"][
            "independent_constraint_reduction"
        ]
        self.assertEqual(
            reduction["target_adjacent_pair_orders"],
            [1.4990947384363291, 1.9817436463819333],
        )
        self.assertEqual(
            reduction["target_adjacent_pair_margins"],
            [-0.0009052615636708783, 0.48174364638193334],
        )
        self.assertTrue(reduction["target_previous_pair_failed"])
        self.assertTrue(reduction["target_finest_pair_passed"])
        self.assertTrue(reduction["target_preasymptotic_on_tested_ladder"])
        self.assertTrue(reduction["complete_constraint_admission_passed"])
        self.assertEqual(
            reduction["classification"],
            "completed_target_and_complete_constraint_pass",
        )

    def test_every_frozen_mutation_control_is_visible(self) -> None:
        controls = self.record["artifact_payload"]["mutation_controls"]
        self.assertEqual(
            set(controls),
            {
                "threshold",
                "grid",
                "time",
                "component",
                "claim",
                "one_bit_state",
                "accepted_stage_count",
            },
        )
        self.assertTrue(all(controls.values()))

    def test_claim_and_threshold_config_mutations_fail_closed(self) -> None:
        for key, value in (
            ("FGCQR_holdout_execution_authorized", True),
            ("PROTO13_frozen", True),
        ):
            with self.subTest(key=key):
                config = pref14.load_config()
                config["claims"][key] = value
                with self.assertRaises(ValueError):
                    pref14.validate_config_data(config)
        config = pref14.load_config()
        config["constraint_reduction"]["minimum_finest_pair_order"] = "149/100"
        with self.assertRaises(ValueError):
            pref14.validate_config_data(config)

    def test_tracked_result_cannot_promote_PROTO13_or_hide_the_old_miss(self) -> None:
        promoted = deepcopy(self.record)
        promoted["gate_status"]["PROTO13_frozen"] = True
        with self.assertRaises(ValueError):
            pref14._validate_tracked_record(promoted)
        hidden = deepcopy(self.record)
        hidden["artifact_payload"]["independent_constraint_reduction"][
            "target_adjacent_pair_orders"
        ][0] = 1.5
        with self.assertRaises(ValueError):
            pref14._validate_tracked_record(hidden)


if __name__ == "__main__":
    unittest.main()
