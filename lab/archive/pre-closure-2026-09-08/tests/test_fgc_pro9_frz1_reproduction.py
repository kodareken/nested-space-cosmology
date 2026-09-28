from __future__ import annotations

import copy
from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(REPOSITORY))

from scripts import reproduce_fgc_pro9_frz1 as reproduction  # noqa: E402


class FGCPro9Frz1ReproductionTests(unittest.TestCase):
    def test_canonical_freeze_reproduction(self) -> None:
        reproduction.verify_canonical()

    def test_freeze_contains_no_runtime_or_candidate_claim(self) -> None:
        result = reproduction.load_canonical_result()
        validation = result["artifact_payload"]["protocol_validation"]
        self.assertEqual(validation["point_counts"], [2049, 4097, 8193])
        self.assertTrue(validation["new_8193_data_required"])
        self.assertFalse(validation["conditional_projection_is_admission"])
        self.assertTrue(all(value is False for value in result["gate_status"].values()))
        self.assertTrue(all(value is False for value in result["nonclaims"].values()))

    def test_config_broadening_fails_closed(self) -> None:
        text = reproduction.DEFAULT_CONFIG.read_text()
        attacked = text.replace(
            "only_the_three_grid_sizes_may_change = true",
            "only_the_three_grid_sizes_may_change = false",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.toml"
            path.write_text(attacked)
            with self.assertRaises(ValueError):
                reproduction.load_config(path)

    def test_result_promotion_fails_reproduction(self) -> None:
        attacked = copy.deepcopy(reproduction.load_canonical_result())
        attacked["gate_status"]["PROTO9_fresh_GR0_dynamic_calibration_authorized"] = True
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.json"
            path.write_text(reproduction._canonical(attacked))
            with self.assertRaises(ValueError):
                reproduction.verify_canonical(path)


if __name__ == "__main__":
    unittest.main()
