from __future__ import annotations

from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_rsp1_pref12 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_config,
    reproduce,
)


class FGCRSP1PREF12ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = reproduce()

    def test_stored_result_reproduces_from_terminal_checkpoint(self) -> None:
        self.assertEqual(DEFAULT_OUTPUT.read_text(encoding="utf-8"), _canonical(self.record))
        recomputation = self.record["artifact_payload"]["checkpoint_recomputation"]
        self.assertTrue(recomputation["all_six_terminal_states_restored"])
        self.assertTrue(recomputation["recomputed_endpoint_equals_terminal_event"])
        self.assertTrue(recomputation["recomputed_endpoint_equals_study_result"])
        self.assertTrue(recomputation["initial_event_equals_frozen_authorization"])

    def test_scoped_clearance_does_not_promote_calibration_or_physics(self) -> None:
        gates = self.record["gate_status"]
        self.assertTrue(gates["RSP1_study_terminated_normally"])
        self.assertTrue(gates["amplitude_three_spectral_veto_cleared_for_successor_design"])
        self.assertTrue(gates["PROTO12_design_may_begin"])
        self.assertFalse(gates["RSP1_generic_all_field_raw_spectral_admission_passed"])
        self.assertFalse(gates["PROTO12_frozen"])
        self.assertFalse(gates["fresh_GR0_dynamic_calibration_completed"])
        self.assertFalse(gates["GR0_case_eligible"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["retained_EFT_evolution_authorized"])
        self.assertFalse(gates["physical_transition_claim_authorized"])

    def test_claim_promotion_and_ratio_weakening_fail_closed(self) -> None:
        text = DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            "FGCQR_holdout_execution_authorized = false",
            "FGCQR_holdout_execution_authorized = true",
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "promoted.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "configuration contract"):
                load_config(path)

        text = DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            'RK4_ratios = ["0.22812049822366257", "0.18656614804836605"]',
            'RK4_ratios = ["0.251", "0.18656614804836605"]',
            1,
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "weakened.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "configured target ratios"):
                load_config(path)

    def test_raw_hash_mutation_fails_before_interpretation(self) -> None:
        text = DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            "0ab89d0ad869df6c5a3f28b1bbcbf0ae099ebf9b03807775b2ddd628a3dc699d",
            "1ab89d0ad869df6c5a3f28b1bbcbf0ae099ebf9b03807775b2ddd628a3dc699d",
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "hash.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "immutable RSP1 study hash"):
                reproduce(path)


if __name__ == "__main__":
    unittest.main()
