"""Focused controls for the recovery-only TDG8 persisted-retry overlay."""
from __future__ import annotations

from dataclasses import replace
import inspect
import unittest

import numpy as np

from recursive_horizons.fgc.evolution import proto15_runtime as p15
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6
from recursive_horizons.fgc.evolution import tdg7_stage_safe_runtime as tdg7
from recursive_horizons.fgc.evolution import tdg8_persisted_retry_replay as replay
from recursive_horizons.fgc.evolution.numerical_engine import (
    PRIMARY_METHOD,
    array_content_sha256,
)
from tests.test_fgc_tdg7_stage_safe_runtime import (
    START,
    TARGET,
    _objects,
    _rhs,
)


class TDG8PersistedRetryReplayTests(unittest.TestCase):
    def _first_successor(self):
        state, transaction, tracers, ledger = _objects()
        coordinates = np.linspace(0.0, 1.0, 9)
        prepared = tdg7.prepare_tdg7_initial_stage_safe_runtime(
            method=PRIMARY_METHOD,
            time=START,
            event_target=TARGET,
            requested_cap=0.125,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=coordinates,
            temporal_ledger=ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        failed = tdg7.TDG7PreparedStageSafeRuntime(
            prepared.plan,
            replace(
                prepared.prepared_tdg6,
                continuous_admission=p15._synthetic_failed_admission(
                    PRIMARY_METHOD
                ),
            ),
        )
        with self.assertRaises(tdg6.TDG6TemporalRetryRequired) as caught:
            tdg7.require_tdg7_stage_safe_admission(
                failed,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=START,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda _evidence: None,
            )
        successor = tdg7.replan_tdg7_after_tdg6_retry(
            failed, caught.exception
        )
        return (
            state,
            transaction,
            tracers,
            coordinates,
            prepared,
            successor,
        )

    def test_overlay_exposes_no_commit_route(self) -> None:
        self.assertFalse(any("commit" in name for name in replay.__all__))
        self.assertNotIn("def commit", inspect.getsource(replay))
        self.assertFalse(hasattr(replay, "commit_tdg8_persisted_retry_replay"))

    def test_exact_active_ledger_replay_is_nonmutating_and_noncommittable(self) -> None:
        state, transaction, tracers, coordinates, _initial, successor = (
            self._first_successor()
        )
        state_before = array_content_sha256(state.u, state.p, state.q)
        tracer_before = tdg7._tracer_snapshot(tracers)
        monitor_before = transaction.state
        causal_before = transaction.causal_state

        prepared = replay.prepare_tdg8_persisted_retry_replay(
            replay_plan=successor.plan,
            expected_prior_retry_count=1,
            method=PRIMARY_METHOD,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=coordinates,
            temporal_ledger=successor.temporal_ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        self.assertEqual(prepared.plan, successor.plan)
        self.assertNotIsInstance(prepared, tdg7.TDG7PreparedStageSafeRuntime)
        with self.assertRaises(TypeError):
            tdg7.commit_tdg7_stage_safe_runtime(
                prepared,  # type: ignore[arg-type]
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=successor.temporal_ledger,
                current_time=START,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
            )
        self.assertEqual(array_content_sha256(state.u, state.p, state.q), state_before)
        self.assertEqual(tdg7._tracer_snapshot(tracers), tracer_before)
        self.assertEqual(transaction.state, monitor_before)
        self.assertEqual(transaction.causal_state, causal_before)

    def test_requested_cap_and_selected_width_remain_distinct(self) -> None:
        state, transaction, tracers, coordinates, initial, successor = (
            self._first_successor()
        )
        self.assertEqual(initial.plan.requested_cap, 0.125)
        self.assertEqual(initial.plan.macro_width, 0.0625)
        prepared = replay.prepare_tdg8_persisted_retry_replay(
            replay_plan=initial.plan,
            expected_prior_retry_count=1,
            method=PRIMARY_METHOD,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=coordinates,
            temporal_ledger=successor.temporal_ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        self.assertEqual(prepared.plan.requested_cap, 0.125)
        self.assertEqual(prepared.plan.macro_width, 0.0625)
        self.assertEqual(
            prepared._prepared.prepared_tdg6.final_time,
            initial.plan.endpoint,
        )

    def test_admission_and_replan_delegate_without_a_commit_surface(self) -> None:
        state, transaction, tracers, coordinates, _initial, successor = (
            self._first_successor()
        )
        prepared = replay.prepare_tdg8_persisted_retry_replay(
            replay_plan=successor.plan,
            expected_prior_retry_count=1,
            method=PRIMARY_METHOD,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=coordinates,
            temporal_ledger=successor.temporal_ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        failed = replay.TDG8PersistedRetryReplay(
            tdg7.TDG7PreparedStageSafeRuntime(
                prepared.plan,
                replace(
                    prepared._prepared.prepared_tdg6,
                    continuous_admission=p15._synthetic_failed_admission(
                        PRIMARY_METHOD
                    ),
                ),
            ),
            1,
        )
        with self.assertRaises(tdg6.TDG6TemporalRetryRequired) as caught:
            replay.require_tdg8_persisted_retry_replay_admission(
                failed,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=successor.temporal_ledger,
                current_time=START,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda _evidence: None,
            )
        next_successor = replay.replan_tdg8_after_persisted_retry_rejection(
            failed, caught.exception
        )
        self.assertLess(next_successor.plan.macro_width, failed.plan.macro_width)

    def test_depth_three_replay_uses_count_two_ledger_and_yields_fourth_plan(self) -> None:
        state, transaction, tracers, coordinates, _initial, second_plan = (
            self._first_successor()
        )
        state_before = array_content_sha256(state.u, state.p, state.q)
        tracer_before = tdg7._tracer_snapshot(tracers)
        monitor_before = transaction.state
        causal_before = transaction.causal_state

        second_prepared = replay.prepare_tdg8_persisted_retry_replay(
            replay_plan=second_plan.plan,
            expected_prior_retry_count=1,
            method=PRIMARY_METHOD,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=coordinates,
            temporal_ledger=second_plan.temporal_ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        second_failed = replay.TDG8PersistedRetryReplay(
            tdg7.TDG7PreparedStageSafeRuntime(
                second_prepared.plan,
                replace(
                    second_prepared._prepared.prepared_tdg6,
                    continuous_admission=p15._synthetic_failed_admission(
                        PRIMARY_METHOD
                    ),
                ),
            ),
            1,
        )
        with self.assertRaises(tdg6.TDG6TemporalRetryRequired) as second_caught:
            replay.require_tdg8_persisted_retry_replay_admission(
                second_failed,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=second_plan.temporal_ledger,
                current_time=START,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda _evidence: None,
            )
        self.assertEqual(
            second_caught.exception.updated_ledger.current_macro_step_temporal_retry_count,
            2,
        )
        third_plan = replay.replan_tdg8_after_persisted_retry_rejection(
            second_failed, second_caught.exception
        )

        third_prepared = replay.prepare_tdg8_persisted_retry_replay(
            replay_plan=third_plan.plan,
            expected_prior_retry_count=2,
            method=PRIMARY_METHOD,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=coordinates,
            temporal_ledger=third_plan.temporal_ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        third_failed = replay.TDG8PersistedRetryReplay(
            tdg7.TDG7PreparedStageSafeRuntime(
                third_prepared.plan,
                replace(
                    third_prepared._prepared.prepared_tdg6,
                    continuous_admission=p15._synthetic_failed_admission(
                        PRIMARY_METHOD
                    ),
                ),
            ),
            2,
        )
        with self.assertRaises(tdg6.TDG6TemporalRetryRequired) as third_caught:
            replay.require_tdg8_persisted_retry_replay_admission(
                third_failed,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=third_plan.temporal_ledger,
                current_time=START,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda _evidence: None,
            )
        self.assertEqual(
            third_caught.exception.updated_ledger.current_macro_step_temporal_retry_count,
            3,
        )
        fourth_plan = replay.replan_tdg8_after_persisted_retry_rejection(
            third_failed, third_caught.exception
        )
        self.assertLess(fourth_plan.plan.macro_width, third_plan.plan.macro_width)
        self.assertEqual(
            fourth_plan.temporal_ledger.current_macro_step_temporal_retry_count,
            3,
        )
        self.assertEqual(array_content_sha256(state.u, state.p, state.q), state_before)
        self.assertEqual(tdg7._tracer_snapshot(tracers), tracer_before)
        self.assertEqual(transaction.state, monitor_before)
        self.assertEqual(transaction.causal_state, causal_before)

    def test_retry_count_and_plan_inputs_fail_closed(self) -> None:
        state, transaction, tracers, coordinates, _initial, successor = (
            self._first_successor()
        )
        for expected in (0, 2, 32):
            with self.subTest(expected=expected), self.assertRaisesRegex(
                ValueError, "TDG8 persisted replay"
            ):
                replay.prepare_tdg8_persisted_retry_replay(
                    replay_plan=successor.plan,
                    expected_prior_retry_count=expected,
                    method=PRIMARY_METHOD,
                    state=state,
                    rhs=_rhs,
                    projector=None,
                    transaction=transaction,
                    tracers=tracers,
                    coordinates=coordinates,
                    temporal_ledger=successor.temporal_ledger,
                    previous_step_index=0,
                    previous_transaction_serial=0,
                )
        with self.assertRaises(TypeError):
            replay.prepare_tdg8_persisted_retry_replay(
                replay_plan=object(),  # type: ignore[arg-type]
                expected_prior_retry_count=1,
                method=PRIMARY_METHOD,
                state=state,
                rhs=_rhs,
                projector=None,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
                temporal_ledger=successor.temporal_ledger,
                previous_step_index=0,
                previous_transaction_serial=0,
            )


if __name__ == "__main__":
    unittest.main()
