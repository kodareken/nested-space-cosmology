from __future__ import annotations

import copy
from pathlib import Path
import sys
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.evolution.protocol_v11 import (  # noqa: E402
    validate_sf1_protocol_v11,
)


class FGCProtocolV11Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v11.toml").open(
            "rb"
        ) as handle:
            cls.protocol = tomllib.load(handle)

    def test_overlay_freezes_only_the_well_balanced_derivative_map(self) -> None:
        record = validate_sf1_protocol_v11(self.protocol)
        self.assertTrue(record["frozen"])
        self.assertTrue(record["resolution_ladder_unchanged"])
        self.assertTrue(record["only_reference_state_differentiation_changes"])
        self.assertEqual(record["point_counts"], [2049, 4097, 8193])
        self.assertEqual(record["initial_q_rule"], "D_h(u-u_ref)+q_ref")
        self.assertEqual(record["source_q_r_rule"], "D_h(q-q_ref)")
        self.assertFalse(record["interior_q_reprojection_permitted"])
        self.assertEqual(record["raw_source_residual_maximum"], "1/1000000000000")
        self.assertFalse(record["new_numeric_tolerance_added"])
        self.assertTrue(all(value is False for value in record["claims"].values()))

    def test_threshold_or_interior_reprojection_broadening_fails_closed(self) -> None:
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["reference_state_differentiation"][
            "raw_source_residual_maximum"
        ] = "1/100000000000"
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v11(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["reference_state_differentiation"][
            "interior_q_is_not_reprojected_after_initialization_or_Runge_Kutta_stages"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v11(attacked)

    def test_reference_fitting_damping_or_incomplete_source_fails_closed(self) -> None:
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["reference_state_differentiation"][
            "reference_is_fixed_and_has_no_fitted_or_evolved_parameters"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v11(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["reference_state_differentiation"][
            "reference_subtraction_is_not_damping_or_constraint_cleaning"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v11(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["reference_state_differentiation"][
            "raw_source_residual_is_computed_from_the_complete_unredefined_equations"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v11(attacked)

    def test_claim_promotion_and_old_output_relabelling_fail_closed(self) -> None:
        attacked = copy.deepcopy(self.protocol)
        attacked["claims"]["PROTO11_fresh_GR0_dynamic_calibration_authorized"] = True
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v11(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["provenance"][
            "PROTO10_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO11_outcomes"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v11(attacked)


if __name__ == "__main__":
    unittest.main()
