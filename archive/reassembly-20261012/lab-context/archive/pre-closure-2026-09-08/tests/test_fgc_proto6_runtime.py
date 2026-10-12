from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    PRIMARY_METHOD,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.proto6_runtime import (  # noqa: E402
    attempt_proto6_step,
)


def _state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _rhs_with_failure(*, call_index: int | None, non_source: bool = False):
    calls = 0

    def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
        nonlocal calls
        fail = calls == call_index
        calls += 1
        zero = np.zeros(state.shape)
        return EvolutionRHS(
            zero,
            zero,
            zero,
            {
                "source_residual_infinity": 2.0e-12 if fail else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if fail and non_source else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


def _rhs_with_failures(
    *, source_calls: set[int], lapse_calls: set[int] | None = None
):
    calls = 0
    lapse_failures = set() if lapse_calls is None else lapse_calls

    def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
        nonlocal calls
        source_failed = calls in source_calls
        lapse_failed = calls in lapse_failures
        calls += 1
        zero = np.zeros(state.shape)
        return EvolutionRHS(
            zero,
            zero,
            zero,
            {
                "source_residual_infinity": 2.0e-12 if source_failed else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if lapse_failed else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


class FGCProto6RuntimeTests(unittest.TestCase):
    def transaction(self) -> GR0RuntimeStageTransaction:
        return GR0RuntimeStageTransaction(
            thresholds=GR0UniversalThresholds(),
            causal_state=CausalBudgetState(previous_speed_upper=1.0),
            boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
            grid_spacing=0.125,
            cfl_maximum=0.125,
        )

    def attempt(self, rhs, transaction, *, preaccept=None):
        return attempt_proto6_step(
            method=PRIMARY_METHOD,
            time=0.0,
            step_size=0.01,
            state=_state(),
            rhs=rhs,
            projector=None,
            transaction=transaction,
            previous_step_index=0,
            previous_transaction_serial=0,
            preaccept=preaccept,
        )

    def test_accepted_state_source_miss_is_terminal_and_latched(self) -> None:
        transaction = self.transaction()
        causal_before = transaction.causal_state
        with self.assertRaises(GR0RuntimeStop) as caught:
            self.attempt(_rhs_with_failure(call_index=0), transaction)
        self.assertEqual(caught.exception.reason, "newton_residual_limit")
        self.assertEqual(transaction.state.accepted_stage_count, 0)
        self.assertEqual(transaction.causal_state, causal_before)
        with self.assertRaises(GR0RuntimeStop) as repeated:
            self.attempt(_rhs_with_failure(call_index=None), transaction)
        self.assertEqual(repeated.exception.reason, "newton_residual_limit")

    def test_internal_source_only_miss_is_retryable_with_exact_rollback(self) -> None:
        transaction = self.transaction()
        state_before = transaction.state
        causal_before = transaction.causal_state
        preaccept_calls = 0

        def preaccept(_proposal):
            nonlocal preaccept_calls
            preaccept_calls += 1

        # Call zero is the accepted-state precheck; call four is RK4 k4.
        result = self.attempt(
            _rhs_with_failure(call_index=4),
            transaction,
            preaccept=preaccept,
        )
        self.assertIsNone(result.accepted)
        self.assertEqual(result.retry.failed_stage_name, "rk4_k4")
        self.assertEqual(result.retry.source_residual_infinity, 2.0e-12)
        self.assertEqual(preaccept_calls, 0)
        self.assertEqual(transaction.state, state_before)
        self.assertEqual(transaction.causal_state, causal_before)

    def test_complete_preview_serializes_every_source_only_failed_stage(self) -> None:
        transaction = self.transaction()
        result = self.attempt(
            _rhs_with_failures(source_calls={2, 4}),
            transaction,
        )
        self.assertIsNotNone(result.retry)
        self.assertEqual(
            [item.stage_name for item in result.retry.failed_stages],
            ["rk4_k2", "rk4_k4"],
        )
        self.assertEqual(transaction.state.accepted_stage_count, 0)

    def test_later_non_source_failure_prevents_earlier_source_retry(self) -> None:
        transaction = self.transaction()
        with self.assertRaises(GR0RuntimeStop) as caught:
            self.attempt(
                _rhs_with_failures(source_calls={2}, lapse_calls={4}),
                transaction,
            )
        # Original chronological priority still reports the earlier source
        # miss, but the complete preview prevents it from becoming retryable.
        self.assertEqual(caught.exception.reason, "newton_residual_limit")

    def test_candidate_endpoint_source_miss_retains_terminal_ownership(self) -> None:
        transaction = self.transaction()
        # Call five is the explicit candidate endpoint, not an internal stage.
        with self.assertRaises(GR0RuntimeStop) as caught:
            self.attempt(_rhs_with_failure(call_index=5), transaction)
        self.assertEqual(caught.exception.reason, "newton_residual_limit")
        self.assertIsNotNone(transaction.state.first_failed_transaction_serial)

    def test_simultaneous_non_source_failure_cannot_be_retried(self) -> None:
        transaction = self.transaction()
        with self.assertRaises(GR0RuntimeStop) as caught:
            self.attempt(
                _rhs_with_failure(call_index=4, non_source=True),
                transaction,
            )
        self.assertEqual(caught.exception.reason, "nonpositive_lapse")
        self.assertIn("newton_residual_limit", caught.exception.failures)

    def test_healthy_attempt_commits_once_after_preaccept(self) -> None:
        transaction = self.transaction()
        ordering = []

        def preaccept(_proposal):
            ordering.append(transaction.state.accepted_stage_count)
            return "tracer-preview"

        result = self.attempt(
            _rhs_with_failure(call_index=None),
            transaction,
            preaccept=preaccept,
        )
        self.assertIsNone(result.retry)
        self.assertEqual(result.preaccept_payload, "tracer-preview")
        self.assertEqual(ordering, [0])
        self.assertEqual(transaction.state.accepted_stage_count, 5)
        self.assertEqual(result.accepted.transaction_serial, 5)


if __name__ == "__main__":
    unittest.main()
