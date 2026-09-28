from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.protocol_v12 import (  # noqa: E402
    validate_sf1_protocol_v12,
)


PROTOCOL = tomllib.loads(
    (REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v12.toml").read_text(
        encoding="utf-8"
    )
)


class FGCProtocolV12Tests(unittest.TestCase):
    def test_canonical_overlay_is_frozen_and_outcome_neutral(self) -> None:
        result = validate_sf1_protocol_v12(PROTOCOL)
        self.assertTrue(result["frozen"])
        self.assertTrue(result["outcome_neutral_contract_validated"])
        self.assertEqual(result["eligible_calibration_amplitudes"], ["5/2", "3"])
        self.assertEqual(result["point_counts"], [2049, 4097, 8193])
        self.assertEqual(
            result["RSP1_qualification_point_counts"], [4097, 8193, 16385]
        )
        self.assertFalse(result["FGCQR_outcomes_inspected_before_revision"])
        self.assertTrue(all(value is False for value in result["claims"].values()))

    def test_threshold_ladder_and_amplitude_mutations_fail_closed(self) -> None:
        for path, value in (
            (("replacement", "source_evaluator", "raw_source_residual_maximum"), "1/100000000000"),
            (("replacement", "spatial_spectral_interpretation", "maximum_nested_tail_ratio"), "1/3"),
            (("replacement", "spatial_spectral_interpretation", "calibration_point_counts"), [4097, 8193, 16385]),
            (("amendment", "eligible_calibration_amplitudes_in_frozen_order"), ["3"]),
        ):
            mutated = deepcopy(PROTOCOL)
            parent = mutated
            for key in path[:-1]:
                parent = parent[key]
            parent[path[-1]] = value
            with self.assertRaisesRegex(ValueError, "differs"):
                validate_sf1_protocol_v12(mutated)

    def test_guard_removal_and_claim_promotion_fail_closed(self) -> None:
        for key in (
            "a_saturated_pair_requires_both_pair_member_absolute_budgets",
            "a_saturated_pair_requires_both_pair_member_top_band_erasures_within_their_own_measured_map_scales",
            "all_three_round_trip_residuals_must_contract_strictly_or_be_exact_zero",
            "whole_profile_infinity_and_RMS_differences_must_contract_strictly_or_be_exact_zero",
        ):
            mutated = deepcopy(PROTOCOL)
            mutated["replacement"]["spatial_spectral_interpretation"][key] = False
            with self.assertRaisesRegex(ValueError, "differs"):
                validate_sf1_protocol_v12(mutated)

        promoted = deepcopy(PROTOCOL)
        promoted["claims"]["FGCQR_holdout_execution_authorized"] = True
        with self.assertRaisesRegex(ValueError, "differs"):
            validate_sf1_protocol_v12(promoted)

    def test_source_backend_may_not_mutate_equations_or_branch(self) -> None:
        for key, value in (
            ("branch", "FGC-QR"),
            ("continuum_equations", "altered_equations"),
            ("maximum_complete_refinement_iterations", 32),
            ("wall_clock_speed_is_not_a_scientific_gate", False),
        ):
            mutated = deepcopy(PROTOCOL)
            mutated["replacement"]["source_evaluator"][key] = value
            with self.assertRaisesRegex(ValueError, "differs"):
                validate_sf1_protocol_v12(mutated)


if __name__ == "__main__":
    unittest.main()
