from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "scripts/reproduce_fgc_id2_all1.py"
CONFIG = REPOSITORY / "configs/fgc/fgc-1-id2-all1.toml"
RESULT = REPOSITORY / "results/fgc-1-id2-all1.json"
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

SPEC = importlib.util.spec_from_file_location("reproduce_fgc_id2_all1", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FGCID2ALL1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.canonical = MODULE.load_canonical_result(RESULT)
        cls.recomputed = MODULE.record(CONFIG)

    def test_canonical_result_reproduces_exactly(self) -> None:
        self.assertEqual(self.recomputed, self.canonical)

    def test_cartesian_product_and_hash_coverage_are_complete(self) -> None:
        payload = self.canonical["artifact_payload"]
        partition = payload["protocol_partition"]
        self.assertEqual(partition["GR0_candidate_count"], 7)
        self.assertEqual(partition["FGCQR_static_slice_count"], 35)
        self.assertEqual(partition["complete_grid_state_hash_count"], 252)
        self.assertTrue(partition["every_complete_grid_state_hash_unique"])
        self.assertEqual(len(payload["GR0_candidate_ledger"]), 7)
        self.assertEqual(len(payload["FGCQR_all_amplitude_all_case_ledger"]), 35)

    def test_static_eligibility_is_ordered_and_outcome_neutral(self) -> None:
        payload = self.canonical["artifact_payload"]
        self.assertEqual(
            payload["static_decision_summary"][
                "eligible_amplitudes_in_frozen_order"
            ],
            ["5/2", "3"],
        )
        ledger = payload["dynamic_calibration_eligibility"]
        self.assertEqual(
            [item["amplitude"] for item in ledger],
            ["2", "5/2", "3", "7/2", "4", "9/2", "5"],
        )
        self.assertTrue(all(item["dynamic_outcome_used"] is False for item in ledger))
        self.assertEqual(
            ledger[0]["ordered_ineligibility_reasons"],
            ["FGCQR-WIDTH-NINE-EIGHTHS:peak_compactness_below_protocol_minimum"],
        )
        self.assertIn(
            "FGCQR-WIDTH-SEVEN-EIGHTHS:peak_compactness_above_protocol_maximum",
            ledger[-1]["ordered_ineligibility_reasons"],
        )

    def test_static_and_nested_diagnostic_counts_are_distinguished(self) -> None:
        summary = self.canonical["artifact_payload"]["static_decision_summary"]
        self.assertEqual(summary["FGCQR_static_slice_pass_count"], 31)
        self.assertEqual(summary["FGCQR_static_slice_fail_count"], 4)
        self.assertEqual(summary["FGCQR_nested_spectral_diagnostic_pass_count"], 17)
        self.assertEqual(summary["FGCQR_nested_spectral_diagnostic_fail_count"], 18)
        self.assertTrue(
            summary[
                "nested_FGCQR_spectral_result_is_disclosed_but_not_a_PROTO3_static_five_boolean_decision"
            ]
        )

    def test_no_dynamic_or_physical_claim_is_promoted(self) -> None:
        boundary = self.canonical["artifact_payload"]["epistemic_boundary"]
        self.assertFalse(boundary["dynamic_trajectory_used"])
        self.assertFalse(boundary["fresh_GR0_calibration_completed"])
        self.assertFalse(boundary["PRO4_HLD1_resolved"])
        self.assertFalse(boundary["mechanism_question_answered"])
        self.assertTrue(all(value is False for value in self.canonical["nonclaims"].values()))
        claims = self.canonical["gate_status"]
        for name in (
            "fresh_GR0_dynamic_calibration_authorized",
            "fresh_GR0_dynamic_calibration_completed",
            "PROTO4_resolved_holdout_manifest_authorized",
            "FGCQR_holdout_execution_authorized",
            "classical_spherical_diagnostic_authorized",
            "retained_EFT_evolution_authorized",
            "physical_transition_claim_authorized",
            "FGCQR_mechanism_rejected",
            "general_gradient_route_rejected",
        ):
            self.assertFalse(claims[name])

    def _mutated_config(self, old: str, new: str) -> Path:
        text = CONFIG.read_text(encoding="utf-8")
        self.assertEqual(text.count(old), 1)
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "mutated.toml"
        path.write_text(text.replace(old, new), encoding="utf-8")
        return path

    def test_amplitude_order_case_and_threshold_mutations_fail(self) -> None:
        changed_order = self._mutated_config(
            'amplitude_candidates = ["2", "5/2", "3", "7/2", "4", "9/2", "5"]',
            'amplitude_candidates = ["5/2", "2", "3", "7/2", "4", "9/2", "5"]',
        )
        with self.assertRaisesRegex(ValueError, "physical inputs"):
            MODULE.load_config(changed_order)

        changed_case = self._mutated_config(
            'case_id = "FGCQR-WIDTH-NINE-EIGHTHS"',
            'case_id = "FGCQR-WIDTH-TEN-EIGHTHS"',
        )
        with self.assertRaisesRegex(ValueError, "case partition"):
            MODULE.load_config(changed_case)

        changed_threshold = self._mutated_config(
            'initial_compactness_maximum = "3/4"',
            'initial_compactness_maximum = "4/5"',
        )
        with self.assertRaisesRegex(ValueError, "threshold contract"):
            MODULE.load_config(changed_threshold)

    def test_claim_promotion_and_predecessor_substitution_fail(self) -> None:
        promoted = self._mutated_config(
            "FGCQR_holdout_execution_authorized = false",
            "FGCQR_holdout_execution_authorized = true",
        )
        with self.assertRaisesRegex(ValueError, "claims"):
            MODULE.load_config(promoted)

        substituted = self._mutated_config(
            'successor_monitor_result = "results/fgc-1-hlt2-mon2.json"',
            'successor_monitor_result = "results/fgc-1-hlt1-mon1.json"',
        )
        with self.assertRaisesRegex(ValueError, "must name"):
            MODULE.load_config(substituted)


if __name__ == "__main__":
    unittest.main()
