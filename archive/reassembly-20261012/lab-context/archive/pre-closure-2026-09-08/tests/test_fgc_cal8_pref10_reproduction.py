from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_cal8_pref10 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_config,
    reproduce,
)


class FGCCAL8PREF10ReproductionTests(unittest.TestCase):
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

    def test_result_localizes_only_numerical_contract_obstructions(self) -> None:
        gates = self.record["gate_status"]
        self.assertTrue(gates["PROTO11_campaign_terminated_normally"])
        self.assertTrue(gates["PROTO11_both_amplitudes_reached_evolved_common_event"])
        self.assertTrue(gates["PROTO11_evolved_affine_source_arithmetic_obstruction_localized"])
        self.assertTrue(gates["PROTO11_direct_coarse_phi_derivative_veto_localized"])
        self.assertTrue(gates["SRC2_affine_source_arithmetic_preflight_required"])
        self.assertFalse(gates["PROTO11_GR0_case_eligible"])
        self.assertFalse(gates["PROTO12_frozen"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["FGCQR_mechanism_rejected"])

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

    def test_threshold_mutation_in_config_fails_closed(self) -> None:
        text = DEFAULT_CONFIG.read_text(encoding="utf-8").replace(
            'raw_source_residual_maximum = "1/1000000000000"',
            'raw_source_residual_maximum = "1/500000000000"',
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mutated.toml"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "failed source evaluation"):
                reproduce(path)


if __name__ == "__main__":
    unittest.main()
