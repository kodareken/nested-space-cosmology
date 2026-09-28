from __future__ import annotations

from math import pi
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.numerical_controls import (  # noqa: E402
    WaveControl,
    boundary_isolation_control,
    checkpoint_restart_control,
    convergence_control,
    minkowski_preservation_control,
    runtime_transaction_control,
    sbp_and_dissipation_controls,
)
from scripts.reproduce_fgc_hlt1_mon1 import load_config as load_hlt1_config  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
)


class FGCNumericalControlsTests(unittest.TestCase):
    def test_outcome_blind_controls_pass_both_independent_methods(self) -> None:
        manufactured = WaveControl("MMS-CENTER", 1, 1.0, True)
        gr0 = WaveControl("GR0-SCALAR", -1, 1.0, False)
        thresholds = {
            PRIMARY_METHOD: 2.5,
            COMPARATOR_METHOD: 1.9,
        }
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            with self.subTest(method=method):
                self.assertTrue(minkowski_preservation_control(method=method)["bitwise_preserved"])
                for control in (manufactured, gr0):
                    convergence = convergence_control(
                        control,
                        method=method,
                        point_counts=(33, 65, 129),
                        final_time=0.2,
                    )
                    self.assertGreaterEqual(
                        convergence["minimum_solution_order"], thresholds[method]
                    )
                    self.assertGreater(
                        convergence["minimum_constraint_order"], 1.5
                    )
                restart = checkpoint_restart_control(method=method)
                self.assertTrue(restart["canonical_roundtrip_bitwise"])
                self.assertTrue(restart["continuous_and_restart_endpoint_bitwise_equal"])
                self.assertLessEqual(
                    boundary_isolation_control(method=method)["maximum_measurement_difference"],
                    1.0e-13,
                )

        for record in sbp_and_dissipation_controls()["records"]:
            self.assertLessEqual(record["sbp_identity_residual"], 1.0e-14)
            self.assertTrue(record["dissipation_nonpositive"])

        transaction = runtime_transaction_control(load_hlt1_config()["thresholds"])
        self.assertTrue(transaction["all_transactions_passed"])
        self.assertEqual(len(transaction["records"]), 2)


if __name__ == "__main__":
    unittest.main()
