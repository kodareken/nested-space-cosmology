from __future__ import annotations

from copy import deepcopy
import tomllib
import unittest

from scripts import reproduce_fgc_pro14_frz1 as pro14
from recursive_horizons.fgc.evolution.protocol_v14 import (
    PROTO14_RESTART,
    PROTO14_TEMPORAL_ADMISSION,
    PROTO14_TERMINAL,
    validate_sf1_protocol_v14,
)


class FGCProtocolV14Tests(unittest.TestCase):
    def setUp(self) -> None:
        with pro14.PROTOCOL_CONFIG.open("rb") as handle:
            self.protocol = tomllib.load(handle)

    def test_overlay_validates_and_is_prospective_only(self) -> None:
        result = validate_sf1_protocol_v14(self.protocol)
        self.assertEqual(result["artifact_id"], "FGC-2-SF1-PROTO14")
        self.assertEqual(result["protocol_version"], 14)
        self.assertTrue(result["replacement_temporal_admission_defined"])
        self.assertFalse(result["historical_PROTO13_result_reclassified"])
        self.assertFalse(result["production_trajectory_authorized"])
        self.assertTrue(all(value is False for value in result["claims"].values()))

    def test_restart_boundary_is_calibration_only(self) -> None:
        restart = self.protocol["replacement"]["restart"]
        self.assertEqual(restart, PROTO14_RESTART)
        self.assertEqual(restart["restart_coordinate_time"], "23/16")
        self.assertTrue(restart["restart_arrays_are_fixed_discrete_calibration_inputs"])
        self.assertFalse(restart["pre_restart_temporal_error_bound_proved"])
        self.assertFalse(restart["restart_may_seed_candidate_or_physical_claim"])
        self.assertEqual(restart["TDG6_accumulated_debit_initial_values"], [0.0] * 18)

    def test_temporal_replacement_and_terminal_classes_are_exact(self) -> None:
        temporal = self.protocol["replacement"]["temporal_admission"]
        terminal = self.protocol["replacement"]["terminal_classification"]
        self.assertEqual(temporal, PROTO14_TEMPORAL_ADMISSION)
        self.assertEqual(terminal, PROTO14_TERMINAL)
        self.assertEqual(temporal["level_step_counts"], [1, 2, 4])
        self.assertEqual(temporal["exact_outward_order_test"], "8*U12^2<=L01^2")
        self.assertTrue(temporal["sampled_64_history_rule_retained_as_nonveto_diagnostic"])
        self.assertFalse(temporal["sampled_64_history_rule_may_admit_or_reject_PROTO14"])
        self.assertEqual(
            terminal["TDG6_retry_exhaustion_classification"],
            "invalid_implementation_or_nonconverged_run",
        )
        self.assertFalse(terminal["candidate_branch_may_open_from_PROTO14_result_alone"])

    def test_threshold_history_restart_and_claim_mutations_fail_closed(self) -> None:
        for name, mutate in {
            "threshold": lambda value: value["replacement"]["temporal_admission"].__setitem__(
                "minimum_observed_order", "1"
            ),
            "history_veto": lambda value: value["replacement"]["temporal_admission"].__setitem__(
                "sampled_64_history_rule_may_admit_or_reject_PROTO14", True
            ),
            "restart_promotion": lambda value: value["replacement"]["restart"].__setitem__(
                "pre_restart_temporal_error_bound_proved", True
            ),
            "candidate_promotion": lambda value: value["replacement"]["terminal_classification"].__setitem__(
                "candidate_branch_may_open_from_PROTO14_result_alone", True
            ),
            "claim": lambda value: value["claims"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
        }.items():
            with self.subTest(name=name):
                attacked = deepcopy(self.protocol)
                mutate(attacked)
                with self.assertRaises(ValueError):
                    validate_sf1_protocol_v14(attacked)


if __name__ == "__main__":
    unittest.main()
