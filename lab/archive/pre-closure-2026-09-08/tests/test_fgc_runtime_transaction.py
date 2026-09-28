from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt1_mon1 import (  # noqa: E402
    _baseline_snapshot,
    load_config,
)
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    ClassicalHealthMonitor,
    ClassicalHealthStop,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    propose_step,
)
from recursive_horizons.fgc.evolution.runtime_transaction import (  # noqa: E402
    RuntimeStageTransaction,
)


class FGCRuntimeTransactionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.thresholds = load_config()["thresholds"]

    @staticmethod
    def proposal():
        initial = EvolutionState(
            np.zeros((9, 6)), np.zeros((9, 6)), np.zeros((9, 6))
        )

        def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
            zero = np.zeros(state.shape)
            return EvolutionRHS(zero, zero, zero, {})

        return propose_step(
            method=PRIMARY_METHOD,
            time=0.0,
            step_size=0.125,
            state=initial,
            rhs=rhs,
        )

    def transaction(self, snapshot_factory):
        return RuntimeStageTransaction(
            health_monitor=ClassicalHealthMonitor(self.thresholds),
            causal_state=CausalBudgetState(previous_speed_upper=1.0),
            boundary_geometry=BoundaryGeometry(96.0, 24.0, 16.0, 0.375),
            speed_evaluator=lambda _stage: 1.25,
            snapshot_factory=snapshot_factory,
        )

    def test_all_stages_and_endpoint_commit_atomically(self) -> None:
        def snapshots(stage, serial, margin):
            return replace(
                _baseline_snapshot(),
                time=stage.time,
                stage_index=serial,
                boundary_margin_over_required_buffer=margin,
            )

        transaction = self.transaction(snapshots)
        receipt = transaction(self.proposal())
        self.assertEqual(receipt["accepted_stage_count"], 5)
        self.assertEqual(transaction.health_monitor.state.accepted_stage_count, 5)
        self.assertEqual(transaction.causal_state.accepted_time, 0.125)
        self.assertEqual(transaction.causal_state.previous_speed_upper, 1.25)

    def test_late_stage_failure_rolls_back_all_stage_accepts_and_causal_debit(self) -> None:
        def snapshots(stage, serial, margin):
            snapshot = replace(
                _baseline_snapshot(),
                time=stage.time,
                stage_index=serial,
                boundary_margin_over_required_buffer=margin,
            )
            if stage.stage_name == "candidate_endpoint":
                snapshot = replace(snapshot, minimum_lapse=0.0)
            return snapshot

        transaction = self.transaction(snapshots)
        original_causal = transaction.causal_state
        with self.assertRaises(ClassicalHealthStop) as caught:
            transaction(self.proposal())
        self.assertEqual(caught.exception.reason, "nonpositive_lapse")
        self.assertEqual(transaction.health_monitor.state.accepted_stage_count, 0)
        self.assertEqual(transaction.health_monitor.state.last_accepted_stage_index, -1)
        self.assertEqual(transaction.causal_state, original_causal)
        with self.assertRaises(ClassicalHealthStop) as repeated:
            transaction(self.proposal())
        self.assertEqual(repeated.exception.reason, "nonpositive_lapse")


if __name__ == "__main__":
    unittest.main()
