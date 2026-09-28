from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.run_fgc_gr0_calibration import NormalFlowTracers  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
)
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg6_temporal_admission_runtime as runtime,
)


def _state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 1] = 0.025
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _rhs(
    *,
    failed_time: float | None = None,
    source_only: bool = False,
):
    def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
        failed = failed_time is not None and (
            np.float64(time).tobytes() == np.float64(failed_time).tobytes()
        )
        return EvolutionRHS(
            state.u,
            state.p,
            state.q,
            {
                "source_residual_infinity": 2.0 if failed and source_only else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if failed and not source_only else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return rhs


def _transaction() -> GR0RuntimeStageTransaction:
    return GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(previous_speed_upper=1.0),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=1.0,
        cfl_maximum=1.0,
    )


def _tracers(state: EvolutionState) -> NormalFlowTracers:
    return NormalFlowTracers.create(
        minimum=0.25,
        maximum=0.75,
        spacing=0.25,
        state=state,
        coordinates=np.linspace(0.0, 1.0, 9),
        cutoff=1.0,
        outer_radius=128.0,
    )


def _interval(value: float, *, subintervals: int):
    coefficients = np.tile(
        np.asarray((value, 0.0, 0.0, 0.0), dtype=np.float64),
        (subintervals, 1),
    )
    return runtime._magnitude_interval(
        coefficients,
        np.zeros_like(coefficients),
        subinterval_count=subintervals,
        owned_row_count=1,
    )


def _synthetic_admission(
    *, method: str = PRIMARY_METHOD, failed_channel: str | None = None
):
    outer = {
        channel: _interval(16.0, subintervals=2)
        for channel in TDG6_COMPLETE_STATE_CHANNELS
    }
    fine = {
        channel: _interval(1.0, subintervals=4)
        for channel in TDG6_COMPLETE_STATE_CHANNELS
    }
    if failed_channel is not None:
        outer[failed_channel] = _interval(1.0e-12, subintervals=2)
        fine[failed_channel] = _interval(9.0e-13, subintervals=4)
    return runtime.classify_tdg6_runtime_intervals(
        method=method, outer_intervals=outer, finest_intervals=fine
    )


class _FailOnceTracers(NormalFlowTracers):
    fail_next_commit: bool = False

    def commit_advance(self, positions, proper_times) -> None:
        if self.fail_next_commit:
            self.fail_next_commit = False
            raise RuntimeError("injected tracer commit failure")
        super().commit_advance(positions, proper_times)


class TDG6TemporalAdmissionRuntimeTests(unittest.TestCase):
    def prepare(
        self,
        method: str = PRIMARY_METHOD,
        *,
        rhs=None,
        transaction=None,
        state=None,
        tracers=None,
        ledger=None,
        step_size: float = 0.125,
    ):
        initial = _state() if state is None else state
        return runtime.prepare_tdg6_gr0_compositor(
            method=method,
            time=0.0,
            step_size=step_size,
            state=initial,
            rhs=_rhs() if rhs is None else rhs,
            projector=None,
            transaction=_transaction() if transaction is None else transaction,
            tracers=_tracers(initial) if tracers is None else tracers,
            coordinates=np.linspace(0.0, 1.0, 9),
            temporal_ledger=(
                runtime.TDG6TemporalLedger.zero(initial_time=0.0)
                if ledger is None
                else ledger
            ),
            previous_step_index=0,
            previous_transaction_serial=0,
        )

    def failed_prepared_context(self):
        state = _state()
        transaction = _transaction()
        tracers = _tracers(state)
        ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
        prepared = self.prepare(
            transaction=transaction,
            state=state,
            tracers=tracers,
            ledger=ledger,
        )
        prepared = replace(
            prepared,
            continuous_admission=_synthetic_admission(
                failed_channel=TDG6_COMPLETE_STATE_CHANNELS[7]
            ),
        )
        return prepared, state, transaction, tracers, ledger

    def test_both_methods_prepare_seven_independent_complete_paths(self) -> None:
        for method, shadow_stages, fine_stages in (
            (PRIMARY_METHOD, 35, 20),
            (COMPARATOR_METHOD, 28, 16),
        ):
            with self.subTest(method=method):
                state = _state()
                transaction = _transaction()
                tracers = _tracers(state)
                ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
                monitor_before = transaction.state
                causal_before = transaction.causal_state
                tracer_before = runtime._tracer_snapshot(tracers)
                prepared = self.prepare(
                    method,
                    transaction=transaction,
                    state=state,
                    tracers=tracers,
                    ledger=ledger,
                )
                self.assertEqual(
                    len(prepared.outer.attempts)
                    + len(prepared.medium.attempts)
                    + len(prepared.fine.attempts),
                    7,
                )
                self.assertEqual(
                    prepared.outer.stage_record_count
                    + prepared.medium.stage_record_count
                    + prepared.fine.stage_record_count,
                    shadow_stages,
                )
                self.assertEqual(prepared.fine.stage_record_count, fine_stages)
                self.assertTrue(prepared.continuous_admission.complete_admission_passed)
                self.assertEqual(transaction.state, monitor_before)
                self.assertEqual(transaction.causal_state, causal_before)
                self.assertEqual(runtime._tracer_snapshot(tracers), tracer_before)
                self.assertEqual(prepared.original_temporal_ledger, ledger)
                self.assertNotEqual(
                    prepared.outer.final_tracer_snapshot.positions_sha256,
                    prepared.medium.final_tracer_snapshot.positions_sha256,
                )
                self.assertNotEqual(
                    prepared.medium.final_tracer_snapshot.positions_sha256,
                    prepared.fine.final_tracer_snapshot.positions_sha256,
                )

    def test_only_complete_fine_path_commits_all_runtime_ledgers(self) -> None:
        for method, expected_stages, expected_serial in (
            (PRIMARY_METHOD, 20, 20),
            (COMPARATOR_METHOD, 16, 16),
        ):
            with self.subTest(method=method):
                state = _state()
                transaction = _transaction()
                tracers = _tracers(state)
                ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
                prepared = self.prepare(
                    method,
                    transaction=transaction,
                    state=state,
                    tracers=tracers,
                    ledger=ledger,
                )
                committed = runtime.commit_tdg6_gr0_compositor(
                    prepared,
                    transaction=transaction,
                    tracers=tracers,
                    temporal_ledger=ledger,
                    current_time=0.0,
                    current_state=state,
                    current_step_index=0,
                    current_transaction_serial=0,
                )
                self.assertEqual(committed.committed_path, "fine_four_quarter_steps")
                self.assertFalse(committed.outer_path_committed)
                self.assertFalse(committed.medium_path_committed)
                self.assertTrue(committed.fine_quarter_steps_committed_atomically)
                self.assertEqual(committed.time, 0.125)
                self.assertEqual(committed.step_index, 4)
                self.assertEqual(committed.transaction_serial, expected_serial)
                self.assertEqual(transaction.state.accepted_stage_count, expected_stages)
                self.assertEqual(
                    runtime._tracer_snapshot(tracers),
                    prepared.fine.final_tracer_snapshot,
                )
                self.assertEqual(committed.temporal_ledger.accepted_macro_step_count, 1)
                self.assertEqual(committed.temporal_ledger.last_accepted_time, 0.125)
                self.assertEqual(
                    committed.temporal_ledger.accumulated_debit_vector,
                    prepared.continuous_admission.finest_pair_debit_vector,
                )
                self.assertEqual(
                    array_content_sha256(
                        committed.state.u, committed.state.p, committed.state.q
                    ),
                    array_content_sha256(
                        prepared.fine.final_accepted.state.u,
                        prepared.fine.final_accepted.state.p,
                        prepared.fine.final_accepted.state.q,
                    ),
                )

    def test_late_fine_non_source_failure_preserves_every_real_ledger(self) -> None:
        state = _state()
        transaction = _transaction()
        tracers = _tracers(state)
        ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
        monitor_before = transaction.state
        causal_before = transaction.causal_state
        tracer_before = runtime._tracer_snapshot(tracers)
        state_hash = array_content_sha256(state.u, state.p, state.q)
        with self.assertRaises(runtime.TDG6RefinementPathStop) as caught:
            self.prepare(
                rhs=_rhs(failed_time=0.078125, source_only=False),
                transaction=transaction,
                state=state,
                tracers=tracers,
                ledger=ledger,
            )
        self.assertEqual(caught.exception.path, "fine_2")
        self.assertIsNotNone(caught.exception.cause)
        self.assertIsNone(caught.exception.source_retry)
        self.assertEqual(transaction.state, monitor_before)
        self.assertEqual(transaction.causal_state, causal_before)
        self.assertEqual(runtime._tracer_snapshot(tracers), tracer_before)
        self.assertEqual(array_content_sha256(state.u, state.p, state.q), state_hash)
        self.assertEqual(ledger, runtime.TDG6TemporalLedger.zero(initial_time=0.0))

    def test_source_only_shadow_failure_retains_PROTO7_ownership(self) -> None:
        state = _state()
        transaction = _transaction()
        tracers = _tracers(state)
        ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
        with self.assertRaises(runtime.TDG6RefinementPathStop) as caught:
            self.prepare(
                rhs=_rhs(failed_time=0.078125, source_only=True),
                transaction=transaction,
                state=state,
                tracers=tracers,
                ledger=ledger,
            )
        self.assertEqual(caught.exception.path, "fine_2")
        self.assertTrue(caught.exception.source_retry_owned_by_PROTO7)
        self.assertIsNotNone(caught.exception.source_retry)
        self.assertFalse(caught.exception.temporal_retry)
        self.assertEqual(transaction.state.accepted_stage_count, 0)
        self.assertEqual(ledger.cumulative_temporal_retry_count, 0)

    def test_all_of_channel_gate_separates_magnitude_from_convergence(self) -> None:
        passing = _synthetic_admission()
        self.assertTrue(passing.complete_admission_passed)
        self.assertEqual(passing.failed_channels, ())
        failed_channel = TDG6_COMPLETE_STATE_CHANNELS[7]
        failed = _synthetic_admission(failed_channel=failed_channel)
        self.assertFalse(failed.complete_admission_passed)
        self.assertEqual(failed.failed_channels, (failed_channel,))
        owner = failed.channel_admissions[7]
        self.assertEqual(owner.classification, "resolved_order_failure")
        self.assertLess(owner.finest_pair_debit, 1.0e-12)

        large_outer = {
            channel: _interval(80.0, subintervals=2)
            for channel in TDG6_COMPLETE_STATE_CHANNELS
        }
        large_fine = {
            channel: _interval(10.0, subintervals=4)
            for channel in TDG6_COMPLETE_STATE_CHANNELS
        }
        large = runtime.classify_tdg6_runtime_intervals(
            method=PRIMARY_METHOD,
            outer_intervals=large_outer,
            finest_intervals=large_fine,
        )
        self.assertTrue(large.complete_admission_passed)
        self.assertGreater(min(large.finest_pair_debit_vector), 9.0)
        self.assertFalse(large.absolute_state_tolerance_used)
        self.assertFalse(large.physical_signal_used_for_normalization)

    def test_temporal_retry_is_durable_separate_and_nonphysical(self) -> None:
        prepared, state, transaction, tracers, ledger = (
            self.failed_prepared_context()
        )
        sink: list[dict[str, object]] = []
        with self.assertRaises(runtime.TDG6TemporalRetryRequired) as caught:
            runtime.require_tdg6_temporal_admission(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda evidence: sink.append(dict(evidence)),
            )
        retry = caught.exception
        self.assertEqual(len(sink), 1)
        self.assertEqual(sink[0]["event_type"], "rejected_TDG6_temporal_admission")
        self.assertEqual(retry.retry_step_size, 0.0625)
        self.assertFalse(retry.physical_classification)
        self.assertEqual(retry.updated_ledger.cumulative_temporal_retry_count, 1)
        self.assertEqual(
            retry.updated_ledger.current_macro_step_temporal_retry_count, 1
        )
        self.assertEqual(retry.updated_ledger.accepted_macro_step_count, 0)
        self.assertEqual(
            retry.updated_ledger.accumulated_debit_vector,
            ledger.accumulated_debit_vector,
        )
        self.assertEqual(len(retry.updated_ledger.serialized_temporal_rejections), 1)

    def test_success_after_temporal_retry_resets_only_the_current_counter(self) -> None:
        prepared, state, transaction, tracers, ledger = (
            self.failed_prepared_context()
        )
        with self.assertRaises(runtime.TDG6TemporalRetryRequired) as caught:
            runtime.require_tdg6_temporal_admission(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda _item: None,
            )
        retried_ledger = caught.exception.updated_ledger
        retried = self.prepare(
            transaction=transaction,
            state=state,
            tracers=tracers,
            ledger=retried_ledger,
            step_size=caught.exception.retry_step_size,
        )
        committed = runtime.commit_tdg6_gr0_compositor(
            retried,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=retried_ledger,
            current_time=0.0,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
        )
        self.assertEqual(committed.temporal_ledger.cumulative_temporal_retry_count, 1)
        self.assertEqual(
            committed.temporal_ledger.last_accepted_macro_step_temporal_retry_count,
            1,
        )
        self.assertEqual(
            committed.temporal_ledger.current_macro_step_temporal_retry_count,
            0,
        )
        self.assertEqual(len(committed.temporal_ledger.serialized_temporal_rejections), 1)

    def test_failed_durable_sink_cannot_create_a_retry_ledger(self) -> None:
        prepared, state, transaction, tracers, ledger = (
            self.failed_prepared_context()
        )

        def failed_sink(_evidence):
            raise OSError("injected durable sink failure")

        with self.assertRaisesRegex(OSError, "durable sink"):
            runtime.require_tdg6_temporal_admission(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=failed_sink,
            )
        self.assertEqual(ledger.cumulative_temporal_retry_count, 0)
        self.assertEqual(ledger.serialized_temporal_rejections, ())

    def test_retry_boundary_drift_fails_before_durable_serialization(self) -> None:
        prepared, state, transaction, tracers, ledger = (
            self.failed_prepared_context()
        )
        tracers.positions[0] = np.nextafter(tracers.positions[0], np.inf)
        sink: list[dict[str, object]] = []
        with self.assertRaisesRegex(RuntimeError, "boundary drifted"):
            runtime.require_tdg6_temporal_admission(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda item: sink.append(dict(item)),
            )
        self.assertEqual(sink, [])
        self.assertEqual(ledger.cumulative_temporal_retry_count, 0)

    def test_retry_count_exhaustion_is_numerical_and_serialized(self) -> None:
        prepared, state, transaction, tracers, ledger = (
            self.failed_prepared_context()
        )
        sink: list[dict[str, object]] = []
        for _ in range(TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP):
            prepared = replace(prepared, original_temporal_ledger=ledger)
            with self.assertRaises(runtime.TDG6TemporalRetryRequired) as caught:
                runtime.require_tdg6_temporal_admission(
                    prepared,
                    transaction=transaction,
                    tracers=tracers,
                    temporal_ledger=ledger,
                    current_time=0.0,
                    current_state=state,
                    current_step_index=0,
                    current_transaction_serial=0,
                    durable_rejection_sink=lambda item: sink.append(dict(item)),
                )
            ledger = caught.exception.updated_ledger
        prepared = replace(prepared, original_temporal_ledger=ledger)
        with self.assertRaises(runtime.TDG6TemporalRetryExhausted) as caught:
            runtime.require_tdg6_temporal_admission(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda item: sink.append(dict(item)),
            )
        self.assertEqual(caught.exception.reason, "maximum_temporal_retries")
        self.assertFalse(caught.exception.physical_classification)
        self.assertEqual(len(sink), TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP + 1)
        self.assertEqual(
            caught.exception.updated_ledger.cumulative_temporal_retry_count,
            TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP + 1,
        )

    def test_commit_rejects_tracer_or_debit_drift_before_mutation(self) -> None:
        state = _state()
        transaction = _transaction()
        tracers = _tracers(state)
        ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
        prepared = self.prepare(
            transaction=transaction,
            state=state,
            tracers=tracers,
            ledger=ledger,
        )
        monitor_before = transaction.state
        tracers.positions[0] = np.nextafter(tracers.positions[0], np.inf)
        with self.assertRaisesRegex(RuntimeError, "drifted"):
            runtime.commit_tdg6_gr0_compositor(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
            )
        self.assertEqual(transaction.state, monitor_before)

        tracers = _tracers(state)
        prepared = self.prepare(
            transaction=transaction,
            state=state,
            tracers=tracers,
            ledger=ledger,
        )
        drifted = replace(ledger, accepted_macro_step_count=1)
        with self.assertRaisesRegex(RuntimeError, "drifted"):
            runtime.commit_tdg6_gr0_compositor(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=drifted,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
            )
        self.assertEqual(transaction.state, monitor_before)

    def test_commit_rolls_back_monitor_causal_and_tracer_on_late_exception(self) -> None:
        state = _state()
        base = _tracers(state)
        tracers = _FailOnceTracers(**base.__dict__)
        transaction = _transaction()
        ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
        prepared = self.prepare(
            transaction=transaction,
            state=state,
            tracers=tracers,
            ledger=ledger,
        )
        monitor_before = transaction.state
        causal_before = transaction.causal_state
        tracer_before = runtime._tracer_snapshot(tracers)
        tracers.fail_next_commit = True
        with self.assertRaisesRegex(RuntimeError, "injected tracer"):
            runtime.commit_tdg6_gr0_compositor(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
            )
        self.assertEqual(transaction.state, monitor_before)
        self.assertEqual(transaction.causal_state, causal_before)
        self.assertEqual(runtime._tracer_snapshot(tracers), tracer_before)

    def test_checkpoint_round_trip_covers_debit_retry_and_rejections(self) -> None:
        prepared, state, transaction, tracers, ledger = (
            self.failed_prepared_context()
        )
        with self.assertRaises(runtime.TDG6TemporalRetryRequired) as caught:
            runtime.require_tdg6_temporal_admission(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda _item: None,
            )
        ledger = caught.exception.updated_ledger
        metadata, debit = runtime.tdg6_checkpoint_extension(ledger)
        restored = runtime.restore_tdg6_checkpoint_extension(metadata, debit)
        self.assertEqual(restored, ledger)
        self.assertTrue(metadata["complete_rejection_evidence_present"])
        self.assertEqual(len(metadata["serialized_temporal_rejections"]), 1)

        mutation = debit.copy()
        mutation[0] = np.nextafter(mutation[0], np.inf)
        with self.assertRaisesRegex(ValueError, "hash differs"):
            runtime.restore_tdg6_checkpoint_extension(metadata, mutation)
        promoted = dict(metadata)
        promoted["physical_classification"] = True
        with self.assertRaisesRegex(ValueError, "numerical scope"):
            runtime.restore_tdg6_checkpoint_extension(promoted, debit)

    def test_interval_lower_and_upper_enclose_dense_polynomial_control(self) -> None:
        coefficients = np.asarray(
            (
                (0.0, 4.0, -4.0, 0.0),
                (0.0, 9.0 / 16.0, -1.5, 1.0),
            )
        )
        interval = runtime._magnitude_interval(
            coefficients,
            np.zeros_like(coefficients),
            subinterval_count=2,
            owned_row_count=1,
        )
        theta = np.linspace(0.0, 1.0, 131_073)
        dense = float(
            np.max(
                np.abs(
                    coefficients[:, 0, None]
                    + theta
                    * (
                        coefficients[:, 1, None]
                        + theta
                        * (
                            coefficients[:, 2, None]
                            + theta * coefficients[:, 3, None]
                        )
                    )
                )
            )
        )
        self.assertLessEqual(interval.lower_bound, dense)
        self.assertLessEqual(dense, interval.upper_bound)

    def test_malformed_channel_checkpoint_and_time_inputs_fail_closed(self) -> None:
        complete = _synthetic_admission()
        outer = {
            item.channel: item.outer_difference
            for item in complete.channel_admissions
        }
        fine = {
            item.channel: item.finest_difference
            for item in complete.channel_admissions
        }
        reordered = dict(reversed(tuple(outer.items())))
        with self.assertRaisesRegex(ValueError, "reordered"):
            runtime.classify_tdg6_runtime_intervals(
                method=PRIMARY_METHOD,
                outer_intervals=reordered,
                finest_intervals=fine,
            )
        ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
        metadata, debit = runtime.tdg6_checkpoint_extension(ledger)
        missing = dict(metadata)
        missing.pop("channel_order")
        with self.assertRaisesRegex(ValueError, "incomplete"):
            runtime.restore_tdg6_checkpoint_extension(missing, debit)
        transaction = _transaction()
        transaction.causal_state = CausalBudgetState(
            accepted_time=1.0e16, previous_speed_upper=1.0
        )
        with self.assertRaises(ValueError):
            runtime.prepare_tdg6_gr0_compositor(
                method=PRIMARY_METHOD,
                time=1.0e16,
                step_size=2.0,
                state=_state(),
                rhs=_rhs(),
                projector=None,
                transaction=transaction,
                tracers=_tracers(_state()),
                coordinates=np.linspace(0.0, 1.0, 9),
                temporal_ledger=runtime.TDG6TemporalLedger.zero(
                    initial_time=1.0e16
                ),
                previous_step_index=0,
                previous_transaction_serial=0,
            )


if __name__ == "__main__":
    unittest.main()
