from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(REPOSITORY))

from scripts import reproduce_fgc_cal5_pref7 as reproduction  # noqa: E402


class Cal5Pref7ReproductionTests(unittest.TestCase):
    def test_stored_result_matches_full_reconstruction(self) -> None:
        reproduction.verify_canonical()

    def test_result_keeps_every_candidate_and_physical_gate_closed(self) -> None:
        result = reproduction.load_canonical_result()
        gates = result["gate_status"]
        self.assertTrue(gates["PROTO8_campaign_terminated_normally"])
        self.assertTrue(gates["PROTO9_resolution_ladder_revision_required"])
        self.assertFalse(gates["PROTO8_GR0_case_eligible"])
        self.assertFalse(gates["PROTO9_frozen"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["retained_EFT_evolution_authorized"])
        self.assertTrue(all(value is False for value in result["nonclaims"].values()))

    def test_mutated_result_is_rejected(self) -> None:
        attacked = copy.deepcopy(reproduction.load_canonical_result())
        attacked["gate_status"]["PROTO8_GR0_case_eligible"] = True
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.json"
            path.write_text(reproduction._canonical(attacked))
            with self.assertRaises(ValueError):
                reproduction.verify_canonical(path)

    def test_resolution_change_cannot_be_broadened_in_config(self) -> None:
        text = reproduction.DEFAULT_CONFIG.read_text()
        attacked = text.replace(
            "only_resolution_ladder_may_change = true",
            "only_resolution_ladder_may_change = false",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.toml"
            path.write_text(attacked)
            with self.assertRaises(ValueError):
                reproduction.load_config(path)


if __name__ == "__main__":
    unittest.main()
