from __future__ import annotations

import copy
from pathlib import Path
import sys
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.evolution.protocol_v9 import (  # noqa: E402
    validate_sf1_protocol_v9,
)


class FGCProtocolV9Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v9.toml").open("rb") as handle:
            cls.protocol = tomllib.load(handle)

    def test_overlay_validates_and_changes_only_the_resolution_ladder(self) -> None:
        record = validate_sf1_protocol_v9(self.protocol)
        self.assertTrue(record["frozen"])
        self.assertTrue(record["only_resolution_ladder_changes"])
        self.assertEqual(record["old_point_counts"], [1025, 2049, 4097])
        self.assertEqual(record["point_counts"], [2049, 4097, 8193])
        self.assertTrue(record["new_8193_data_required"])
        self.assertFalse(record["conditional_projection_is_admission"])
        self.assertTrue(all(value is False for value in record["claims"].values()))

    def test_resolution_or_threshold_broadening_fails_closed(self) -> None:
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["resolution_ladder"]["point_counts"] = [2049, 4097, 6145]
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v9(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["resolution_ladder"][
            "coarsest_and_finest_constraint_magnitude_guards_unchanged"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v9(attacked)

    def test_claim_promotion_and_old_output_relabelling_fail_closed(self) -> None:
        attacked = copy.deepcopy(self.protocol)
        attacked["claims"]["PROTO9_fresh_GR0_dynamic_calibration_authorized"] = True
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v9(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["resolution_ladder"][
            "no_old_state_may_be_relabelled_as_a_PROTO9_outcome"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v9(attacked)


if __name__ == "__main__":
    unittest.main()
