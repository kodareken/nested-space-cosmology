from __future__ import annotations

from copy import deepcopy
import json
import unittest

from scripts import reproduce_fgc_tdg4_frz1 as freeze


class FGCTDG4FRZ1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = freeze.verify_canonical()

    def test_canonical_prospective_freeze_reproduces(self) -> None:
        self.assertEqual(self.record["artifact_id"], freeze.ARTIFACT_ID)
        self.assertEqual(
            self.record["classification"],
            "prospective_TDG4_sampling_identifiability_theorem_freeze",
        )
        self.assertTrue(
            self.record["gate_status"][
                "TDG4_sampling_identifiability_theorem_design_frozen"
            ]
        )
        self.assertFalse(
            self.record["gate_status"][
                "TDG4_sampling_identifiability_theorem_completed"
            ]
        )

    def test_freeze_reads_no_history_and_authorizes_no_successor_gate(self) -> None:
        payload = self.record["artifact_payload"]
        lineage = payload["immutable_lineage"]
        boundary = payload["claim_boundary"]
        self.assertFalse(lineage["actual_terminal_history_arrays_consumed"])
        self.assertFalse(lineage["campaign_checkpoint_loaded"])
        self.assertFalse(lineage["state_advanced"])
        self.assertFalse(boundary["sampling_identifiability_theorem_completed"])
        self.assertFalse(boundary["current_temporal_spectral_rule_retired"])
        self.assertFalse(boundary["replacement_temporal_gate_design_authorized"])
        self.assertFalse(boundary["replacement_temporal_admission_defined"])
        self.assertFalse(boundary["new_trajectory_authorized"])
        self.assertFalse(boundary["candidate_branch_opened"])

    def test_exact_preflight_covers_every_bin_and_assumption_class(self) -> None:
        payload = self.record["artifact_payload"]
        preflight = payload["exact_no_history_preflight"]
        self.assertTrue(preflight["all_exact_controls_pass"])
        self.assertEqual(
            set(preflight["generic_smooth_kernel_witnesses"]),
            {"28", "29", "30", "31", "32"},
        )
        self.assertEqual(
            set(preflight["exact_uniform_alias_witnesses"]),
            {"28", "29", "30", "31", "32"},
        )
        self.assertEqual(
            preflight["quantitative_regular_bound_control"][
                "coefficient_error_bound"
            ],
            "5/3",
        )
        self.assertEqual(set(payload["mutation_controls"].values()), {True})

    def test_tracked_result_is_canonical(self) -> None:
        raw = freeze.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, freeze._canonical_bytes(stored))

    def test_config_mutations_fail_closed(self) -> None:
        config = freeze.load_config()
        attacks = []
        value = deepcopy(config)
        value["immutable_lineage"]["source_result_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(config)
        value["frozen_theorem"]["sample_count"] = 63
        attacks.append(value)
        value = deepcopy(config)
        value["frozen_theorem"]["top_bins"] = [27, 28, 29, 30, 31]
        attacks.append(value)
        value = deepcopy(config)
        value["frozen_theorem"]["function_space"] = "piecewise_linear_only"
        attacks.append(value)
        value = deepcopy(config)
        value["claims"][
            "TDG4_sampling_identifiability_theorem_completed"
        ] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["replacement_temporal_gate_design_authorized"] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["GR0_case_eligible"] = True
        attacks.append(value)
        for attacked in attacks:
            with self.assertRaises(ValueError):
                freeze.validate_config_data(attacked)


if __name__ == "__main__":
    unittest.main()
