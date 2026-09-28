from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_cal9_pref13 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_config,
    reproduce,
)


class FGCCAL9PREF13ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = reproduce()
        cls.stored = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_canonical_result_reproduces_exactly(self) -> None:
        self.assertEqual(self.record, self.stored)
        self.assertEqual(
            DEFAULT_OUTPUT.read_text(encoding="utf-8"),
            _canonical(self.record),
        )

    def test_result_localizes_only_the_numerical_calibration_veto(self) -> None:
        gates = self.record["gate_status"]
        self.assertTrue(gates["PROTO12_campaign_terminated_normally"])
        self.assertTrue(
            gates["PROTO12_both_amplitudes_reached_late_evolved_common_events"]
        )
        self.assertTrue(gates["SRC4_arithmetic_wall_cleared_for_both_amplitudes"])
        self.assertTrue(
            gates[
                "PROTO12_dual_amplitude_SSPRK3_radial_momentum_order_veto_localized"
            ]
        )
        self.assertTrue(gates["RSP2_high_ladder_constraint_preflight_required"])
        self.assertFalse(gates["PROTO12_GR0_case_eligible"])
        self.assertFalse(gates["PROTO13_frozen"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["FGCQR_mechanism_rejected"])

    def test_retry_provenance_is_public_and_typed(self) -> None:
        diagnosis = self.record["artifact_payload"]["campaign_diagnosis"]
        self.assertEqual(diagnosis["serialized_source_retry_count"], 0)
        self.assertTrue(diagnosis["some_adaptive_CFL_step_reductions_occurred"])
        self.assertTrue(diagnosis["campaign_stop_was_not_CFL_retry_exhaustion"])
        self.assertGreater(
            diagnosis["terminal_CFL_retry_count_by_amplitude_and_member"]["3"]
            ["RK4-4097"],
            0,
        )

    def test_claim_promotion_in_config_fails_closed(self) -> None:
        text = DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            "FGCQR_holdout_execution_authorized = false",
            "FGCQR_holdout_execution_authorized = true",
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "configuration contract"):
                load_config(path)

    def test_threshold_reduction_in_config_fails_closed(self) -> None:
        text = DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            'minimum_constraint_finest_pair_order = "3/2"',
            'minimum_constraint_finest_pair_order = "149/100"',
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "configuration contract"):
                load_config(path)

    def test_CFL_counter_mutation_fails_checkpoint_binding(self) -> None:
        text = DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            "RK4_2049 = 217",
            "RK4_2049 = 218",
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "checkpoint runtime"):
                reproduce(path)


if __name__ == "__main__":
    unittest.main()
