"""Direct synthetic controls for the isolated TDG7 stage-safe adapter."""

# ruff: noqa: E402

from __future__ import annotations

from dataclasses import replace
import inspect
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts.run_fgc_gr0_calibration import NormalFlowTracers  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (
    BoundaryGeometry,
    CausalBudgetState,
)  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)  # noqa: E402
from recursive_horizons.fgc.evolution.proto5_runtime import (
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)  # noqa: E402
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg7_stage_safe_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_MINIMUM_MACRO_STEP,
)  # noqa: E402
from recursive_horizons.fgc.evolution.tdg7_binary64_subdivision_lattice import (
    CoordinateLatticeLimitReached,
    cal11_event24_witness,
)  # noqa: E402


START, TARGET = 23.0 / 16.0, 24.0 / 16.0


def _state() -> EvolutionState:
    r = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 1] = 0.025
    u[:, 2] = 1.0
    u[:, 3] = r
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, np.zeros_like(u), q)


def _rhs(time: float, state: EvolutionState) -> EvolutionRHS:
    return EvolutionRHS(
        state.u,
        state.p,
        state.q,
        {
            "source_residual_infinity": 0.0,
            "source_refinement_iterations": 0,
            "source_residual_decreased_monotonically": True,
            "kinetic_condition_infinity": 1.0,
            "coordinate_speed_upper": 1.0,
            "minimum_lapse": 1.0,
            "minimum_radial_metric": 1.0,
            "minimum_areal_radius_away_from_center": 0.125,
        },
    )


def _objects():
    state = _state()
    transaction = GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(accepted_time=START, previous_speed_upper=1.0),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=1.0,
        cfl_maximum=1.0,
    )
    tracers = NormalFlowTracers.create(
        minimum=0.25,
        maximum=0.75,
        spacing=0.25,
        state=state,
        coordinates=np.linspace(0.0, 1.0, 9),
        cutoff=1.0,
        outer_radius=128.0,
    )
    return state, transaction, tracers, tdg6.TDG6TemporalLedger.zero(initial_time=START)


class TDG7StageSafeRuntimeTests(unittest.TestCase):
    def prepare(self, method=PRIMARY_METHOD):
        state, transaction, tracers, ledger = _objects()
        prepared = runtime.prepare_tdg7_initial_stage_safe_runtime(
            method=method,
            time=START,
            event_target=TARGET,
            requested_cap=0.125,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=np.linspace(0.0, 1.0, 9),
            temporal_ledger=ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        return prepared, state, transaction, tracers, ledger

    def test_public_api_separates_fresh_attempts_from_bound_retries(self):
        self.assertFalse(runtime.FRESH_VS_RETRY_PROTOCOL_DISCRIMINATOR_IMPLEMENTED)
        self.assertTrue(
            runtime.STATELESS_ADAPTER_CANNOT_DETECT_DISCARDED_EXTERNAL_RETRY
        )
        self.assertNotIn("_prepare_tdg7_stage_safe_runtime", runtime.__all__)
        initial_parameters = inspect.signature(
            runtime.prepare_tdg7_initial_stage_safe_runtime
        ).parameters
        retry_parameters = inspect.signature(
            runtime.prepare_tdg7_retry_successor
        ).parameters
        self.assertNotIn("retry_successor", initial_parameters)
        self.assertNotIn("minimum_width", initial_parameters)
        self.assertIn("temporal_ledger", initial_parameters)
        self.assertIn("successor", retry_parameters)
        self.assertNotIn("temporal_ledger", retry_parameters)
        self.assertNotIn("time", retry_parameters)
        self.assertNotIn("requested_cap", retry_parameters)

    def test_shared_plan_actual_stages_and_fine_only_commit_for_both_methods(self):
        for method, shadow_total, fine_total in (
            (PRIMARY_METHOD, 35, 20),
            (COMPARATOR_METHOD, 28, 16),
        ):
            with self.subTest(method=method):
                prepared, state, transaction, tracers, ledger = self.prepare(method)
                self.assertEqual(prepared.plan.one_full[0], prepared.plan.two_half[0])
                self.assertEqual(
                    prepared.plan.two_half[-1], prepared.plan.four_quarter[-1]
                )
                paths = (
                    prepared.prepared_tdg6.outer,
                    prepared.prepared_tdg6.medium,
                    prepared.prepared_tdg6.fine,
                )
                self.assertEqual(
                    sum(path.stage_record_count for path in paths), shadow_total
                )
                self.assertEqual(paths[-1].stage_record_count, fine_total)
                runtime.validate_tdg7_stage_safe_prepared(prepared)
                committed = runtime.commit_tdg7_stage_safe_runtime(
                    prepared,
                    transaction=transaction,
                    tracers=tracers,
                    temporal_ledger=ledger,
                    current_time=START,
                    current_state=state,
                    current_step_index=0,
                    current_transaction_serial=0,
                )
                self.assertTrue(committed.fine_quarter_steps_committed_atomically)
                self.assertFalse(committed.outer_path_committed)
                self.assertFalse(committed.medium_path_committed)

    def test_cal11_cap_and_lattice_stop_precede_all_hooks(self):
        witness = cal11_event24_witness()
        self.assertGreater(witness.conservative_reduction, 0.0)
        calls = []

        class Poison:
            def __getattr__(self, name):
                calls.append(name)
                raise AssertionError(name)

        with self.assertRaises(CoordinateLatticeLimitReached):
            runtime.prepare_tdg7_initial_stage_safe_runtime(
                method=PRIMARY_METHOD,
                time=np.nextafter(2.0, -np.inf),
                event_target=2.0,
                requested_cap=0.125,
                state=Poison(),
                rhs=Poison(),
                projector=Poison(),
                transaction=Poison(),
                tracers=Poison(),
                coordinates=Poison(),
                temporal_ledger=Poison(),
                previous_step_index=0,
                previous_transaction_serial=0,
            )
        self.assertEqual(calls, [])

    def test_every_typed_plan_stop_precedes_runtime_hooks(self):
        class Poison:
            def __getattr__(self, name):
                raise AssertionError(f"runtime hook touched: {name}")

        q = math.ulp(START)
        controls = (
            (math.nan, TARGET, 0.125, "nonfinite_coordinate"),
            (START, START, 0.125, "target_not_ahead"),
            (1.0, TARGET, 0.125, "outside_frozen_time_envelope"),
            (np.nextafter(2.0, -math.inf), 2.0, 0.125, "endpoint_not_q_aligned"),
            (START, START + q, 0.125, "target_not_stage_lattice_aligned"),
            (START, TARGET, q, "no_positive_aligned_width"),
            (START, TARGET, 8.0 * q, "below_minimum_aligned_width"),
            (START, TARGET, 0.0, "nonpositive_requested_cap"),
            (START, TARGET, -1.0, "nonpositive_requested_cap"),
            (START, TARGET, math.nan, "nonfinite_coordinate"),
            (START, TARGET, math.inf, "nonfinite_coordinate"),
        )
        for time, target, cap, reason in controls:
            with (
                self.subTest(reason=reason),
                self.assertRaisesRegex(runtime.TDG7CoordinateLatticeStop, reason),
            ):
                runtime.prepare_tdg7_initial_stage_safe_runtime(
                    method=PRIMARY_METHOD,
                    time=time,
                    event_target=target,
                    requested_cap=cap,
                    state=Poison(),
                    rhs=Poison(),
                    projector=Poison(),
                    transaction=Poison(),
                    tracers=Poison(),
                    coordinates=Poison(),
                    temporal_ledger=Poison(),
                    previous_step_index=0,
                    previous_transaction_serial=0,
                )
        with self.assertRaises(TypeError):
            runtime.prepare_tdg7_initial_stage_safe_runtime(
                method=PRIMARY_METHOD,
                time=START,
                event_target=TARGET,
                requested_cap=0.125,
                state=Poison(),
                rhs=Poison(),
                projector=Poison(),
                transaction=Poison(),
                tracers=Poison(),
                coordinates=Poison(),
                temporal_ledger=Poison(),
                previous_step_index=0,
                previous_transaction_serial=0,
                minimum_width=0.0,
            )

    def test_tdg6_all_of_admission_controls_remain_unmodified(self):
        def interval(value, count):
            coefficients = np.tile(np.asarray((value, 0.0, 0.0, 0.0)), (count, 1))
            return tdg6._magnitude_interval(
                coefficients,
                np.zeros_like(coefficients),
                subinterval_count=count,
                owned_row_count=1,
            )

        exact_zero = {
            channel: interval(0.0, 2) for channel in TDG6_COMPLETE_STATE_CHANNELS
        }
        exact_zero_fine = {
            channel: interval(0.0, 4) for channel in TDG6_COMPLETE_STATE_CHANNELS
        }
        zero = tdg6.classify_tdg6_runtime_intervals(
            method=PRIMARY_METHOD,
            outer_intervals=exact_zero,
            finest_intervals=exact_zero_fine,
        )
        self.assertTrue(zero.complete_admission_passed)
        outer = {channel: interval(16.0, 2) for channel in TDG6_COMPLETE_STATE_CHANNELS}
        fine = {channel: interval(1.0, 4) for channel in TDG6_COMPLETE_STATE_CHANNELS}
        outer[TDG6_COMPLETE_STATE_CHANNELS[3]] = interval(1.0e-12, 2)
        fine[TDG6_COMPLETE_STATE_CHANNELS[3]] = interval(9.0e-13, 4)
        veto = tdg6.classify_tdg6_runtime_intervals(
            method=PRIMARY_METHOD, outer_intervals=outer, finest_intervals=fine
        )
        self.assertFalse(veto.complete_admission_passed)
        large = tdg6.classify_tdg6_runtime_intervals(
            method=PRIMARY_METHOD,
            outer_intervals={
                channel: interval(80.0, 2) for channel in TDG6_COMPLETE_STATE_CHANNELS
            },
            finest_intervals={
                channel: interval(10.0, 4) for channel in TDG6_COMPLETE_STATE_CHANNELS
            },
        )
        self.assertTrue(large.complete_admission_passed)
        self.assertGreater(min(large.finest_pair_debit_vector), 9.0)

    def test_replan_consumes_tdg6_cap_and_is_strict_or_typed_stop(self):
        prepared, state, transaction, tracers, ledger = self.prepare()

        def interval(value, count):
            coefficients = np.tile(np.asarray((value, 0.0, 0.0, 0.0)), (count, 1))
            return tdg6._magnitude_interval(
                coefficients,
                np.zeros_like(coefficients),
                subinterval_count=count,
                owned_row_count=1,
            )

        outer = {channel: interval(16.0, 2) for channel in TDG6_COMPLETE_STATE_CHANNELS}
        fine = {channel: interval(1.0, 4) for channel in TDG6_COMPLETE_STATE_CHANNELS}
        outer[TDG6_COMPLETE_STATE_CHANNELS[0]] = interval(1.0e-12, 2)
        fine[TDG6_COMPLETE_STATE_CHANNELS[0]] = interval(9.0e-13, 4)
        failed = replace(
            prepared.prepared_tdg6,
            continuous_admission=tdg6.classify_tdg6_runtime_intervals(
                method=PRIMARY_METHOD, outer_intervals=outer, finest_intervals=fine
            ),
        )
        wrapped = runtime.TDG7PreparedStageSafeRuntime(prepared.plan, failed)
        with self.assertRaises(tdg6.TDG6TemporalRetryRequired) as caught:
            runtime.require_tdg7_stage_safe_admission(
                wrapped,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=START,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda evidence: None,
            )
        new = runtime.replan_tdg7_after_tdg6_retry(wrapped, caught.exception)
        self.assertLess(new.plan.macro_width, wrapped.plan.macro_width)
        self.assertNotEqual(new.plan.macro_width, 4.0 * wrapped.plan.quantum)
        with self.assertRaisesRegex(ValueError, "initial entrypoint"):
            runtime.prepare_tdg7_initial_stage_safe_runtime(
                method=PRIMARY_METHOD,
                time=new.plan.current,
                event_target=new.plan.event_target,
                requested_cap=new.plan.requested_cap,
                state=state,
                rhs=_rhs,
                projector=None,
                transaction=transaction,
                tracers=tracers,
                coordinates=np.linspace(0.0, 1.0, 9),
                temporal_ledger=new.temporal_ledger,
                previous_step_index=0,
                previous_transaction_serial=0,
            )
        resumed = runtime.prepare_tdg7_retry_successor(
            successor=new,
            method=PRIMARY_METHOD,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=np.linspace(0.0, 1.0, 9),
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        self.assertEqual(resumed.plan, new.plan)

        second_wrapped = runtime.TDG7PreparedStageSafeRuntime(
            resumed.plan,
            replace(
                resumed.prepared_tdg6,
                continuous_admission=failed.continuous_admission,
            ),
        )
        with self.assertRaises(tdg6.TDG6TemporalRetryRequired) as second_caught:
            runtime.require_tdg7_stage_safe_admission(
                second_wrapped,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=new.temporal_ledger,
                current_time=START,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda evidence: None,
            )
        latest_rejection = (
            second_caught.exception.updated_ledger.serialized_temporal_rejections[-1]
        )
        forged_prefix_ledger = replace(
            second_caught.exception.updated_ledger,
            serialized_temporal_rejections=(latest_rejection, latest_rejection),
        )
        forged_prefix_retry = tdg6.TDG6TemporalRetryRequired(
            second_caught.exception.evidence, forged_prefix_ledger
        )
        with self.assertRaisesRegex(ValueError, "durable TDG6 evidence"):
            runtime.replan_tdg7_after_tdg6_retry(second_wrapped, forged_prefix_retry)

        small_state, small_transaction, small_tracers, small_ledger = _objects()
        twice_minimum = 2.0 * float(TDG6_MINIMUM_MACRO_STEP)
        small_prepared = runtime.prepare_tdg7_initial_stage_safe_runtime(
            method=PRIMARY_METHOD,
            time=START,
            event_target=TARGET,
            requested_cap=twice_minimum,
            state=small_state,
            rhs=_rhs,
            projector=None,
            transaction=small_transaction,
            tracers=small_tracers,
            coordinates=np.linspace(0.0, 1.0, 9),
            temporal_ledger=small_ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        self.assertEqual(small_prepared.plan.macro_width, twice_minimum)
        small_wrapped = runtime.TDG7PreparedStageSafeRuntime(
            small_prepared.plan,
            replace(
                small_prepared.prepared_tdg6,
                continuous_admission=failed.continuous_admission,
            ),
        )
        with self.assertRaises(tdg6.TDG6TemporalRetryRequired) as small_caught:
            runtime.require_tdg7_stage_safe_admission(
                small_wrapped,
                transaction=small_transaction,
                tracers=small_tracers,
                temporal_ledger=small_ledger,
                current_time=START,
                current_state=small_state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda evidence: None,
            )
        minimum_successor = runtime.replan_tdg7_after_tdg6_retry(
            small_wrapped, small_caught.exception
        )
        self.assertEqual(
            minimum_successor.plan.macro_width,
            float(TDG6_MINIMUM_MACRO_STEP),
        )
        minimum_prepared = runtime.prepare_tdg7_retry_successor(
            successor=minimum_successor,
            method=PRIMARY_METHOD,
            state=small_state,
            rhs=_rhs,
            projector=None,
            transaction=small_transaction,
            tracers=small_tracers,
            coordinates=np.linspace(0.0, 1.0, 9),
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        minimum_wrapped = runtime.TDG7PreparedStageSafeRuntime(
            minimum_prepared.plan,
            replace(
                minimum_prepared.prepared_tdg6,
                continuous_admission=failed.continuous_admission,
            ),
        )
        with self.assertRaises(tdg6.TDG6TemporalRetryExhausted) as minimum_stop:
            runtime.require_tdg7_stage_safe_admission(
                minimum_wrapped,
                transaction=small_transaction,
                tracers=small_tracers,
                temporal_ledger=minimum_successor.temporal_ledger,
                current_time=START,
                current_state=small_state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda evidence: None,
            )
        self.assertEqual(minimum_stop.exception.reason, "minimum_macro_step")
        self.assertFalse(minimum_stop.exception.physical_classification)
        self.assertEqual(
            len(minimum_stop.exception.updated_ledger.serialized_temporal_rejections),
            2,
        )

        # TDG6 normally owns the preceding minimum-step stop.  This forged
        # Required wrapper attacks TDG7's defensive response if an impossible
        # caller nevertheless presents the same evidence as retryable.
        impossible_required = tdg6.TDG6TemporalRetryRequired(
            minimum_stop.exception.evidence,
            minimum_stop.exception.updated_ledger,
        )
        with self.assertRaises(runtime.TDG7RetryLatticeExhausted) as exhausted:
            runtime.replan_tdg7_after_tdg6_retry(minimum_wrapped, impossible_required)
        self.assertEqual(exhausted.exception.evidence, minimum_stop.exception.evidence)
        self.assertEqual(
            exhausted.exception.updated_ledger,
            minimum_stop.exception.updated_ledger,
        )
        self.assertEqual(
            exhausted.exception.requested_half_cap,
            float(TDG6_MINIMUM_MACRO_STEP) / 2.0,
        )
        self.assertEqual(
            exhausted.exception.lattice_reason, "below_minimum_aligned_width"
        )
        self.assertFalse(exhausted.exception.physical_classification)

        over_limit = replace(
            caught.exception.evidence,
            retry_count_for_current_macro_step=33,
            cumulative_temporal_retry_count=33,
        )
        forged_retry = tdg6.TDG6TemporalRetryRequired(
            over_limit, caught.exception.updated_ledger
        )
        with self.assertRaisesRegex(ValueError, "exceeds TDG6 retry limit"):
            runtime.replan_tdg7_after_tdg6_retry(wrapped, forged_retry)
        mismatched_method = replace(caught.exception.evidence, method=COMPARATOR_METHOD)
        mismatched_retry = tdg6.TDG6TemporalRetryRequired(
            mismatched_method, caught.exception.updated_ledger
        )
        with self.assertRaisesRegex(ValueError, "differs"):
            runtime.replan_tdg7_after_tdg6_retry(wrapped, mismatched_retry)
        forged_mapping = json.loads(
            caught.exception.updated_ledger.serialized_temporal_rejections[-1]
        )
        forged_mapping["failed_channels"] = []
        wrong_ledger = replace(
            caught.exception.updated_ledger,
            serialized_temporal_rejections=(
                json.dumps(forged_mapping, sort_keys=True, separators=(",", ":")),
            ),
        )
        forged_ledger_retry = tdg6.TDG6TemporalRetryRequired(
            caught.exception.evidence, wrong_ledger
        )
        with self.assertRaisesRegex(ValueError, "durable TDG6 evidence"):
            runtime.replan_tdg7_after_tdg6_retry(wrapped, forged_ledger_retry)

    def test_mutated_boundary_or_stage_is_refused_before_commit(self):
        prepared, state, transaction, tracers, ledger = self.prepare()
        # The lattice model itself refuses malformed 4-boundary/old-endpoint attacks.
        with self.assertRaises((CoordinateLatticeLimitReached, ValueError)):
            bad_plan = replace(
                prepared.plan,
                boundaries=(prepared.plan.boundaries[0], *prepared.plan.boundaries[2:]),
            )
            runtime.TDG7PreparedStageSafeRuntime(bad_plan, prepared.prepared_tdg6)
        outer_attempt = prepared.prepared_tdg6.outer.attempts[0]
        bad_record = replace(outer_attempt.proposal.stages[1], stage_name="rk4_k1")
        bad_proposal = replace(
            outer_attempt.proposal,
            stages=(
                outer_attempt.proposal.stages[0],
                bad_record,
                *outer_attempt.proposal.stages[2:],
            ),
        )
        bad_attempt = replace(outer_attempt, proposal=bad_proposal)
        bad_outer = replace(prepared.prepared_tdg6.outer, attempts=(bad_attempt,))
        with self.assertRaisesRegex(RuntimeError, "stage record differs"):
            runtime.TDG7PreparedStageSafeRuntime(
                prepared.plan, replace(prepared.prepared_tdg6, outer=bad_outer)
            )
        witness = cal11_event24_witness()
        old_endpoint = witness.current + witness.requested_cap
        with self.assertRaises(CoordinateLatticeLimitReached):
            replace(witness, boundaries=(*witness.boundaries[:4], old_endpoint))
        independently_rounded = np.nextafter(witness.boundaries[2], math.inf)
        with self.assertRaises(CoordinateLatticeLimitReached):
            replace(
                witness,
                boundaries=(
                    witness.boundaries[0],
                    witness.boundaries[1],
                    independently_rounded,
                    witness.boundaries[3],
                    witness.boundaries[4],
                ),
            )
        self.assertEqual(
            array_content_sha256(state.u, state.p, state.q),
            array_content_sha256(_state().u, _state().p, _state().q),
        )

    def test_proposal_and_accepted_stage_time_mutations_fail_for_both_methods(self):
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            with self.subTest(method=method):
                (
                    prepared,
                    _state_value,
                    _transaction_value,
                    _tracers_value,
                    _ledger_value,
                ) = self.prepare(method)
                attempt = prepared.prepared_tdg6.outer.attempts[0]
                moved_proposal_record = replace(
                    attempt.proposal.stages[1],
                    time=np.nextafter(attempt.proposal.stages[1].time, math.inf),
                )
                moved_proposal = replace(
                    attempt.proposal,
                    stages=(
                        attempt.proposal.stages[0],
                        moved_proposal_record,
                        *attempt.proposal.stages[2:],
                    ),
                )
                bad_attempt = replace(attempt, proposal=moved_proposal)
                bad_outer = replace(
                    prepared.prepared_tdg6.outer, attempts=(bad_attempt,)
                )
                with self.assertRaisesRegex(RuntimeError, "stage record differs"):
                    runtime.TDG7PreparedStageSafeRuntime(
                        prepared.plan, replace(prepared.prepared_tdg6, outer=bad_outer)
                    )
                moved_accepted_record = replace(
                    attempt.accepted.stage_records[1],
                    time=np.nextafter(attempt.accepted.stage_records[1].time, math.inf),
                )
                moved_accepted = replace(
                    attempt.accepted,
                    stage_records=(
                        attempt.accepted.stage_records[0],
                        moved_accepted_record,
                        *attempt.accepted.stage_records[2:],
                    ),
                )
                bad_accepted_attempt = replace(attempt, accepted=moved_accepted)
                bad_accepted_outer = replace(
                    prepared.prepared_tdg6.outer, attempts=(bad_accepted_attempt,)
                )
                with self.assertRaisesRegex(
                    RuntimeError, "accepted stage record differs"
                ):
                    runtime.TDG7PreparedStageSafeRuntime(
                        prepared.plan,
                        replace(prepared.prepared_tdg6, outer=bad_accepted_outer),
                    )

    def test_source_and_non_source_late_shadow_stops_keep_real_boundaries(self):
        for source_only in (False, True):
            with self.subTest(source_only=source_only):
                state, transaction, tracers, ledger = _objects()
                state_before = array_content_sha256(state.u, state.p, state.q)
                tracer_before = runtime._tracer_snapshot(tracers)
                monitor_before, causal_before = (
                    transaction.state,
                    transaction.causal_state,
                )

                # The third fine proposal's midpoint is not an outer/medium stage.
                def failing_rhs(time, incoming):
                    rhs = _rhs(time, incoming)
                    if np.float64(time).tobytes() == np.float64(1.4765625).tobytes():
                        diagnostics = dict(rhs.diagnostics)
                        diagnostics["source_residual_infinity"] = (
                            2.0 if source_only else 0.0
                        )
                        diagnostics["minimum_lapse"] = 1.0 if source_only else 0.0
                        return EvolutionRHS(rhs.du, rhs.dp, rhs.dq, diagnostics)
                    return rhs

                with self.assertRaises(tdg6.TDG6RefinementPathStop) as caught:
                    runtime.prepare_tdg7_initial_stage_safe_runtime(
                        method=PRIMARY_METHOD,
                        time=START,
                        event_target=TARGET,
                        requested_cap=0.125,
                        state=state,
                        rhs=failing_rhs,
                        projector=None,
                        transaction=transaction,
                        tracers=tracers,
                        coordinates=np.linspace(0.0, 1.0, 9),
                        temporal_ledger=ledger,
                        previous_step_index=0,
                        previous_transaction_serial=0,
                    )
                self.assertEqual(caught.exception.path, "fine_2")
                self.assertEqual(
                    caught.exception.source_retry_owned_by_PROTO7, source_only
                )
                self.assertEqual(
                    array_content_sha256(state.u, state.p, state.q), state_before
                )
                self.assertEqual(runtime._tracer_snapshot(tracers), tracer_before)
                self.assertEqual(transaction.state, monitor_before)
                self.assertEqual(transaction.causal_state, causal_before)

    def test_late_tracer_commit_failure_rolls_back_the_real_tdq6_owned_boundary(self):
        class FailOnceTracers(NormalFlowTracers):
            fail = False

            def commit_advance(self, positions, proper_times):
                if self.fail:
                    self.fail = False
                    raise RuntimeError("injected late tracer failure")
                return super().commit_advance(positions, proper_times)

        state, transaction, base, ledger = _objects()
        tracers = FailOnceTracers(
            base.labels,
            base.positions,
            base.proper_times,
            base.event_proper_times,
            base.event_fields,
            base.cutoff,
            base.outer_radius,
        )
        prepared = runtime.prepare_tdg7_initial_stage_safe_runtime(
            method=PRIMARY_METHOD,
            time=START,
            event_target=TARGET,
            requested_cap=0.125,
            state=state,
            rhs=_rhs,
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=np.linspace(0.0, 1.0, 9),
            temporal_ledger=ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        before = runtime._tracer_snapshot(tracers)
        monitor, causal = transaction.state, transaction.causal_state
        tracers.fail = True
        with self.assertRaisesRegex(RuntimeError, "injected late tracer failure"):
            runtime.commit_tdg7_stage_safe_runtime(
                prepared,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=START,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
            )
        self.assertEqual(runtime._tracer_snapshot(tracers), before)
        self.assertEqual(transaction.state, monitor)
        self.assertEqual(transaction.causal_state, causal)
