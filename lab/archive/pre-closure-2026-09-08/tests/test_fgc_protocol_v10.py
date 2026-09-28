from __future__ import annotations

import copy
from pathlib import Path
import sys
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.evolution.protocol_v10 import (  # noqa: E402
    validate_sf1_protocol_v10,
)


class FGCProtocolV10Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v10.toml").open(
            "rb"
        ) as handle:
            cls.protocol = tomllib.load(handle)

    def test_overlay_validates_and_changes_only_finest_pair_interpretation(self) -> None:
        record = validate_sf1_protocol_v10(self.protocol)
        self.assertTrue(record["frozen"])
        self.assertTrue(record["resolution_ladder_unchanged"])
        self.assertTrue(record["only_finest_pair_spectral_interpretation_changes"])
        self.assertEqual(record["point_counts"], [2049, 4097, 8193])
        self.assertEqual(record["maximum_nested_tail_ratio"], "1/4")
        self.assertTrue(record["coarse_to_medium_direct_ratio_veto"])
        self.assertTrue(record["finest_pair_direct_ratio_remains_primary"])
        self.assertFalse(record["round_trip_is_a_continuum_error_bound"])
        self.assertFalse(record["diagnostic_saturation_is_physical_resolution_evidence"])
        self.assertTrue(all(value is False for value in record["claims"].values()))

    def test_threshold_resolution_or_coarse_veto_broadening_fails_closed(self) -> None:
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["spectral_interpretation"]["point_counts"] = [
            2049,
            4097,
            16385,
        ]
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v10(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["spectral_interpretation"][
            "maximum_nested_tail_ratio"
        ] = "1/2"
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v10(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["spectral_interpretation"][
            "coarse_to_medium_field_and_derivative_ratios_must_pass_directly"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v10(attacked)

    def test_omitted_guard_or_promoted_semantics_fail_closed(self) -> None:
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["spectral_interpretation"][
            "whole_profile_infinity_and_RMS_differences_must_contract_strictly_or_be_exact_zero"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v10(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["spectral_interpretation"][
            "round_trip_interpolation_diagnostic_is_not_a_continuum_error_bound"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v10(attacked)

    def test_claim_promotion_and_old_output_relabelling_fail_closed(self) -> None:
        attacked = copy.deepcopy(self.protocol)
        attacked["claims"]["PROTO10_fresh_GR0_dynamic_calibration_authorized"] = True
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v10(attacked)
        attacked = copy.deepcopy(self.protocol)
        attacked["replacement"]["provenance"][
            "PROTO9_outputs_are_immutable_diagnostic_history_and_forbidden_as_PROTO10_outcomes"
        ] = False
        with self.assertRaises(ValueError):
            validate_sf1_protocol_v10(attacked)


if __name__ == "__main__":
    unittest.main()
