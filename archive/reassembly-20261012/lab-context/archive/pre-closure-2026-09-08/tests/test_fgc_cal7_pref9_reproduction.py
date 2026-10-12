from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_cal7_pref9 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_config,
    reproduce,
)


class FGCCAL7PREF9ReproductionTests(unittest.TestCase):
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

    def test_result_localizes_only_a_numerical_source_floor(self) -> None:
        claims = self.record["gate_status"]
        self.assertTrue(claims["PROTO10_campaign_terminated_normally"])
        self.assertTrue(claims["PROTO10_center_adjacent_binary64_source_floor_localized"])
        self.assertTrue(claims["PROTO11_well_balanced_reference_map_required"])
        self.assertFalse(claims["PROTO11_frozen"])
        self.assertFalse(claims["fresh_GR0_dynamic_calibration_completed"])
        self.assertFalse(claims["FGCQR_holdout_execution_authorized"])
        self.assertFalse(claims["FGCQR_mechanism_rejected"])
        comparison = self.record["artifact_payload"]["well_balanced_control_summary"]["comparison"]
        self.assertTrue(comparison["all_well_balanced_source_residuals_pass_unchanged_raw_tolerance"])
        self.assertTrue(comparison["well_balanced_control_is_not_a_runtime_or_campaign_result"])

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

    def test_result_promotion_differs_from_canonical_bytes(self) -> None:
        mutated = json.loads(json.dumps(self.record))
        mutated["gate_status"]["FGCQR_mechanism_rejected"] = True
        self.assertNotEqual(_canonical(mutated), DEFAULT_OUTPUT.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
