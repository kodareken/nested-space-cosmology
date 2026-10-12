from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.protocol_v13 import (  # noqa: E402
    validate_sf1_protocol_v13,
)


PROTOCOL = tomllib.loads(
    (REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v13.toml").read_text(
        encoding="utf-8"
    )
)


class FGCProtocolV13Tests(unittest.TestCase):
    def test_canonical_overlay_is_frozen_and_outcome_neutral(self) -> None:
        result = validate_sf1_protocol_v13(PROTOCOL)
        self.assertTrue(result["frozen"])
        self.assertTrue(result["outcome_neutral_contract_validated"])
        self.assertTrue(result["premise_revision_only"])
        self.assertEqual(result["eligible_calibration_amplitudes"], ["3"])
        self.assertEqual(result["primary_point_counts"], [2049, 4097, 8193])
        self.assertEqual(result["comparator_point_counts"], [4097, 8193, 16385])
        self.assertEqual(result["restart_coordinate_time"], "23/16")
        self.assertEqual(result["final_coordinate_time"], "32")
        self.assertTrue(all(value is False for value in result["claims"].values()))

    def test_ladder_threshold_restart_and_partition_mutations_fail_closed(self) -> None:
        for path, value in (
            (("replacement", "method_owned_ladders", "primary_point_counts"), [4097, 8193, 16385]),
            (("replacement", "method_owned_ladders", "comparator_point_counts"), [2049, 4097, 8193]),
            (("replacement", "method_owned_ladders", "minimum_constraint_finest_pair_order"), "149/100"),
            (("replacement", "restart", "restart_coordinate_time"), "3/2"),
            (("replacement", "restart", "final_coordinate_time"), "31"),
            (("replacement", "partition", "calibration_amplitude"), "5/2"),
        ):
            with self.subTest(path=path):
                mutated = deepcopy(PROTOCOL)
                parent = mutated
                for key in path[:-1]:
                    parent = parent[key]
                parent[path[-1]] = value
                with self.assertRaisesRegex(ValueError, "differs"):
                    validate_sf1_protocol_v13(mutated)

    def test_inherited_equations_and_claims_cannot_be_promoted(self) -> None:
        inherited = deepcopy(PROTOCOL)
        inherited["inheritance"][
            "PROTO12_unredefined_GR0_REF1_equations_and_SRC4_source_backend_inherited_unchanged"
        ] = False
        with self.assertRaisesRegex(ValueError, "differs"):
            validate_sf1_protocol_v13(inherited)

        promoted = deepcopy(PROTOCOL)
        promoted["claims"]["FGCQR_holdout_execution_authorized"] = True
        with self.assertRaisesRegex(ValueError, "differs"):
            validate_sf1_protocol_v13(promoted)

    def test_cross_method_collocation_shortcut_fails_closed(self) -> None:
        mutated = deepcopy(PROTOCOL)
        mutated["replacement"]["cross_method_observables"][
            "raw_finest_grid_values_may_not_be_compared_as_if_collocated"
        ] = False
        with self.assertRaisesRegex(ValueError, "differs"):
            validate_sf1_protocol_v13(mutated)


if __name__ == "__main__":
    unittest.main()
