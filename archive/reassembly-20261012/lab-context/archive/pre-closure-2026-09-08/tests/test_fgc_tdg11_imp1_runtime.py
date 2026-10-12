"""Synthetic IMP1 state, guard, debit and restart controls; no campaign data."""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
import inspect
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_imp1_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_enclosure import (  # noqa: E402
    IMP1EnclosureLimits,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    imp1_checkpoint_extension,
    restore_imp1_checkpoint_extension,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    make_synthetic_fixture as _objects,
    synthetic_diagnostics as _diagnostics,
    synthetic_exponential_rhs as _rhs,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
)


def _prepare(
    method: str = PRIMARY_METHOD, *, rhs=_rhs, saturated: bool = False, target=TARGET
):
    state, transaction, tracers, ledger, coordinates = _objects(
        method, saturated=saturated
    )
    prepared = runtime.prepare_imp1_initial_runtime(
        method=method,
        time=START,
        event_target=target,
        requested_cap=1.0 / 32,
        state=state,
        rhs=rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    return prepared, state, transaction, tracers, ledger, coordinates


def _current(prepared, state, transaction, tracers, ledger, coordinates):
    return dict(
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=ledger,
        current_time=prepared.initial_time,
        current_state=state,
        current_step_index=prepared.previous_step_index,
        current_transaction_serial=prepared.previous_transaction_serial,
        coordinates=coordinates,
    )


def _quarter_impulse(time, current):
    """An intentionally non-smooth synthetic nonpass, not a physical RHS."""

    du = np.zeros_like(current.u)
    if time == START + 1.0 / 128:
        du[:, 4] = 1.0
    return EvolutionRHS(du, np.zeros_like(du), np.zeros_like(du), _diagnostics())


class IMP1RuntimeTests(unittest.TestCase):
    def test_genuine_nonpass_requires_durable_sink_then_bound_retry(self):
        items = _prepare(rhs=_quarter_impulse)
        prepared, state, tx, tracers, ledger, coordinates = items
        self.assertFalse(prepared.assessment.admission_passed)
        before = (runtime._tracer_snapshot(tracers), tx.state, tx.causal_state)
        written = []
        with self.assertRaises(runtime.IMP1TemporalRetryRequired) as caught:
            runtime.require_imp1_admission(
                prepared,
                **_current(*items),
                durable_rejection_sink=lambda record: written.append(
                    canonical_json_bytes(record)
                ),
            )
        self.assertEqual(len(written), 1)
        self.assertEqual(
            before, (runtime._tracer_snapshot(tracers), tx.state, tx.causal_state)
        )
        self.assertEqual(ledger.current_macro_step_temporal_retry_count, 0)
        retry = caught.exception
        self.assertEqual(
            retry.updated_ledger.current_macro_step_temporal_retry_count, 1
        )
        self.assertEqual(
            retry.updated_ledger.accumulated_debit_vector,
            ledger.accumulated_debit_vector,
        )
        self.assertEqual(retry.evidence["preparation_sha256"], prepared.receipt_sha256)
        successor = runtime.replan_imp1_retry(prepared, retry)
        self.assertEqual(successor.plan.macro_width, prepared.plan.macro_width / 2)
        with self.assertRaisesRegex(ValueError, "cannot replace"):
            runtime.prepare_imp1_retry_runtime(
                successor,
                state=state,
                rhs=lambda time, value: _quarter_impulse(time, value),
                projector=None,
                transaction=tx,
                tracers=tracers,
                coordinates=coordinates,
            )
        with self.assertRaises(ValueError):
            runtime.prepare_imp1_initial_runtime(
                method=PRIMARY_METHOD,
                time=START,
                event_target=TARGET,
                requested_cap=1.0 / 64,
                state=state,
                rhs=_quarter_impulse,
                projector=None,
                transaction=tx,
                tracers=tracers,
                coordinates=coordinates,
                temporal_ledger=retry.updated_ledger,
                previous_step_index=0,
                previous_transaction_serial=0,
            )
        next_prepared = runtime.prepare_imp1_retry_runtime(
            successor,
            state=state,
            rhs=_quarter_impulse,
            projector=None,
            transaction=tx,
            tracers=tracers,
            coordinates=coordinates,
        )
        self.assertTrue(next_prepared.assessment.admission_passed)
        next_items = (
            next_prepared,
            state,
            tx,
            tracers,
            retry.updated_ledger,
            coordinates,
        )
        committed = runtime.commit_imp1_runtime(next_prepared, **_current(*next_items))
        self.assertEqual(
            committed.temporal_ledger.current_macro_step_temporal_retry_count, 0
        )
        self.assertEqual(committed.temporal_ledger.cumulative_temporal_retry_count, 1)
        self.assertEqual(
            committed.temporal_ledger.last_accepted_macro_step_temporal_retry_count, 1
        )
        self.assertEqual(
            committed.temporal_ledger.serialized_temporal_rejections,
            retry.updated_ledger.serialized_temporal_rejections,
        )

    def test_failed_or_mutating_durable_sink_cannot_return_retry(self):
        items = _prepare(rhs=_quarter_impulse)
        prepared, state, tx, tracers, ledger, coordinates = items
        before = (
            runtime._tracer_snapshot(tracers),
            tx.state,
            tx.causal_state,
            canonical_json_bytes(imp1_checkpoint_extension(ledger)),
        )

        def failed(record):
            raise OSError("injected durable sink failure")

        with self.assertRaises(OSError):
            runtime.require_imp1_admission(
                prepared, **_current(*items), durable_rejection_sink=failed
            )

        def changed(record):
            record["assessment_sha256"] = "2" * 64

        with self.assertRaisesRegex(RuntimeError, "sink altered"):
            runtime.require_imp1_admission(
                prepared, **_current(*items), durable_rejection_sink=changed
            )
        self.assertEqual(
            before,
            (
                runtime._tracer_snapshot(tracers),
                tx.state,
                tx.causal_state,
                canonical_json_bytes(imp1_checkpoint_extension(ledger)),
            ),
        )

    def test_rejection_cannot_be_rebound_to_another_genuine_preparation(self):
        first = _prepare(rhs=_quarter_impulse)
        other = _prepare(rhs=_quarter_impulse, target=25.0 / 16.0)
        self.assertEqual(
            first[0].recorded_family.family_sha256,
            other[0].recorded_family.family_sha256,
        )
        self.assertNotEqual(first[0].receipt_sha256, other[0].receipt_sha256)
        with self.assertRaises(runtime.IMP1TemporalRetryRequired) as caught:
            runtime.require_imp1_admission(
                first[0], **_current(*first), durable_rejection_sink=lambda record: None
            )
        with self.assertRaises(ValueError):
            runtime.replan_imp1_retry(other[0], caught.exception)

    def test_monitor_count_and_partial_failure_validation_precede_source(self):
        for changed in (
            {"accepted_stage_count": -20},
            {"accepted_stage_count": 2},
            {"accepted_stage_count": True},
            {"first_failed_time": START},
        ):
            with self.subTest(changed=changed):
                state, tx, tracers, ledger, coordinates = _objects()
                tx.state = replace(tx.state, **changed)
                calls = []
                with self.assertRaises(ValueError):
                    runtime.prepare_imp1_initial_runtime(
                        method=PRIMARY_METHOD,
                        time=START,
                        event_target=TARGET,
                        requested_cap=0.03125,
                        state=state,
                        rhs=lambda *args: calls.append(args),
                        projector=None,
                        transaction=tx,
                        tracers=tracers,
                        coordinates=coordinates,
                        temporal_ledger=ledger,
                        previous_step_index=0,
                        previous_transaction_serial=0,
                    )
                self.assertEqual(calls, [])

    def test_advanced_monitor_time_cannot_be_stale(self):
        items = _prepare()
        committed = runtime.commit_imp1_runtime(items[0], **_current(*items))
        _, _, tx, tracers, _, coordinates = items
        tx.state = replace(tx.state, last_accepted_time=START)
        with self.assertRaises(ValueError):
            runtime.prepare_imp1_initial_runtime(
                method=PRIMARY_METHOD,
                time=committed.time,
                event_target=TARGET,
                requested_cap=0.03125,
                state=committed.state,
                rhs=_rhs,
                projector=None,
                transaction=tx,
                tracers=tracers,
                coordinates=coordinates,
                temporal_ledger=committed.temporal_ledger,
                previous_step_index=committed.step_index,
                previous_transaction_serial=committed.transaction_serial,
            )

    def test_both_methods_commit_restore_and_prepare_next_step(self):
        for method, records, serial in (
            (PRIMARY_METHOD, 35, 20),
            (COMPARATOR_METHOD, 28, 16),
        ):
            with self.subTest(method=method):
                items = _prepare(method)
                prepared, state, tx, tracers, ledger, coordinates = items
                initial_state_hash = array_content_sha256(state.u, state.p, state.q)
                inherited = ledger.inherited_snapshot
                old_debit = ledger.accumulated_debit_vector
                self.assertEqual(prepared.stage_record_count, records)
                self.assertTrue(prepared.assessment.admission_passed)
                self.assertEqual(len(prepared.assessment.public_fine_debits), 18)
                self.assertEqual(tx.state.accepted_stage_count, 0)
                committed = runtime.commit_imp1_runtime(prepared, **_current(*items))
                self.assertEqual(committed.step_index, 4)
                self.assertEqual(committed.transaction_serial, serial)
                self.assertIs(committed.state, prepared.fine.final_accepted.state)
                self.assertEqual(tx.state.accepted_stage_count, serial)
                self.assertFalse(committed.corrected_endpoint_adopted)
                self.assertFalse(committed.outer_path_adopted)
                self.assertFalse(committed.medium_path_adopted)
                self.assertFalse(committed.campaign_execution_authorized)
                self.assertEqual(
                    committed.temporal_ledger.inherited_snapshot, inherited
                )
                self.assertEqual(
                    committed.temporal_ledger.accumulated_debit_vector,
                    tuple(
                        a + b
                        for a, b in zip(
                            old_debit,
                            prepared.assessment.public_fine_debits,
                            strict=True,
                        )
                    ),
                )
                self.assertEqual(
                    array_content_sha256(state.u, state.p, state.q), initial_state_hash
                )
                encoded = imp1_checkpoint_extension(committed.temporal_ledger)
                restored = restore_imp1_checkpoint_extension(encoded)
                self.assertEqual(
                    canonical_json_bytes(encoded),
                    canonical_json_bytes(imp1_checkpoint_extension(restored)),
                )
                next_prepared = runtime.prepare_imp1_initial_runtime(
                    method=method,
                    time=committed.time,
                    event_target=TARGET,
                    requested_cap=1.0 / 32,
                    state=committed.state,
                    rhs=_rhs,
                    projector=None,
                    transaction=tx,
                    tracers=tracers,
                    coordinates=coordinates,
                    temporal_ledger=restored,
                    previous_step_index=committed.step_index,
                    previous_transaction_serial=committed.transaction_serial,
                )
                self.assertEqual(next_prepared.initial_time, committed.time)
                self.assertTrue(next_prepared.assessment.admission_passed)

    def test_lattice_stop_is_before_poison_runtime_hooks(self):
        class Poison:
            def __getattr__(self, name):
                raise AssertionError(f"premature runtime hook: {name}")

        with self.assertRaises(runtime.IMP1CoordinateLatticeStop):
            runtime.prepare_imp1_initial_runtime(
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

    def test_row_resource_stop_precedes_rhs(self):
        state, tx, tracers, ledger, coordinates = _objects()
        calls = []
        with self.assertRaises(runtime.IMP1RuntimeResourceStop):
            runtime.prepare_imp1_initial_runtime(
                method=PRIMARY_METHOD,
                time=START,
                event_target=TARGET,
                requested_cap=0.03125,
                state=state,
                rhs=lambda *args: calls.append(args),
                projector=None,
                transaction=tx,
                tracers=tracers,
                coordinates=coordinates,
                temporal_ledger=ledger,
                previous_step_index=0,
                previous_transaction_serial=0,
                limits=replace(IMP1EnclosureLimits(), maximum_owned_rows=1),
            )
        self.assertEqual(calls, [])

    def test_saturated_actual_endpoint_is_not_replaced_by_corrected_endpoint(self):
        state, tx, tracers, ledger, coordinates = _objects(saturated=True)

        def constant_rhs(time, current):
            du = np.zeros_like(current.u)
            du[:, 4] = 1.0
            return EvolutionRHS(
                du, np.zeros_like(du), np.zeros_like(du), _diagnostics()
            )

        prepared = runtime.prepare_imp1_initial_runtime(
            method=PRIMARY_METHOD,
            time=START,
            event_target=START + 0.25,
            requested_cap=0.25,
            state=state,
            rhs=constant_rhs,
            projector=None,
            transaction=tx,
            tracers=tracers,
            coordinates=coordinates,
            temporal_ledger=ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
        items = (prepared, state, tx, tracers, ledger, coordinates)
        index = TDG6_COMPLETE_STATE_CHANNELS.index("u:phi")
        self.assertEqual(prepared.assessment.public_fine_debits[index], Fraction(1, 4))
        committed = runtime.commit_imp1_runtime(prepared, **_current(*items))
        self.assertTrue(np.array_equal(committed.state.u[:, 4], state.u[:, 4]))
        self.assertEqual(
            committed.temporal_ledger.accumulated_debit_vector[index], Fraction(1, 4)
        )
        self.assertFalse(committed.corrected_endpoint_adopted)

    def test_source_path_failure_preserves_real_boundary(self):
        state, tx, tracers, ledger, coordinates = _objects()
        before = (
            runtime._tracer_snapshot(tracers),
            tx.state,
            tx.causal_state,
            canonical_json_bytes(imp1_checkpoint_extension(ledger)),
        )

        def bad_rhs(time, current):
            if time > START:
                raise ArithmeticError("synthetic source failure")
            return _rhs(time, current)

        with self.assertRaises(runtime.IMP1RefinementPathStop) as caught:
            runtime.prepare_imp1_initial_runtime(
                method=PRIMARY_METHOD,
                time=START,
                event_target=TARGET,
                requested_cap=0.03125,
                state=state,
                rhs=bad_rhs,
                projector=None,
                transaction=tx,
                tracers=tracers,
                coordinates=coordinates,
                temporal_ledger=ledger,
                previous_step_index=0,
                previous_transaction_serial=0,
            )
        self.assertFalse(caught.exception.temporal_retry)
        self.assertEqual(
            before,
            (
                runtime._tracer_snapshot(tracers),
                tx.state,
                tx.causal_state,
                canonical_json_bytes(imp1_checkpoint_extension(ledger)),
            ),
        )

    def test_stale_coordinates_monitor_guard_or_tracer_history_refuse_commit(self):
        for mutation in ("coordinates", "monitor", "guard", "history", "causal"):
            with self.subTest(mutation=mutation):
                items = _prepare()
                prepared, state, tx, tracers, ledger, coordinates = items
                if mutation == "coordinates":
                    coordinates[2] += 1.0 / 1024
                elif mutation == "monitor":
                    tx.state = replace(tx.state, accepted_stage_count=1)
                elif mutation == "guard":
                    tx.thresholds = replace(
                        tx.thresholds, source_residual_maximum=1.0e-6
                    )
                elif mutation == "history":
                    tracers.event_fields[0][0, 0] += 1.0
                else:
                    tx.causal_state = replace(tx.causal_state, previous_speed_upper=2.0)
                with self.assertRaises(ValueError):
                    runtime.commit_imp1_runtime(prepared, **_current(*items))
                self.assertEqual(ledger.accepted_macro_step_count, 0)

    def test_recorded_candidate_mutation_is_detected(self):
        items = _prepare()
        prepared = items[0]
        data = prepared.fine.final_accepted.state.u
        data.setflags(write=True)
        data[1, 4] = np.nextafter(data[1, 4], np.inf)
        data.setflags(write=False)
        with self.assertRaises(ValueError):
            runtime.commit_imp1_runtime(prepared, **_current(*items))

    def test_late_tracer_failure_rolls_back_complete_history(self):
        items = _prepare()
        prepared, state, tx, tracers, ledger, coordinates = items
        before = (runtime._tracer_snapshot(tracers), tx.state, tx.causal_state)
        original_commit = tracers.commit_advance

        def fail_after_mutation(positions, proper):
            original_commit(positions, proper)
            tracers.event_fields.append(np.ones((3, 6)))
            tracers.event_proper_times.append(np.ones(3))
            tracers.labels = tracers.labels + 0.125
            raise RuntimeError("injected late tracer commit")

        with patch.object(tracers, "commit_advance", side_effect=fail_after_mutation):
            with self.assertRaisesRegex(RuntimeError, "injected late tracer"):
                runtime.commit_imp1_runtime(prepared, **_current(*items))
        self.assertEqual(
            before, (runtime._tracer_snapshot(tracers), tx.state, tx.causal_state)
        )
        self.assertEqual(ledger.accepted_macro_step_count, 0)

    def test_new_ledger_failure_precedes_any_real_mutation(self):
        items = _prepare()
        prepared, state, tx, tracers, ledger, coordinates = items
        before = (runtime._tracer_snapshot(tracers), tx.state, tx.causal_state)
        with patch.object(
            runtime, "accept_imp1_step", side_effect=RuntimeError("ledger capacity")
        ):
            with self.assertRaisesRegex(RuntimeError, "ledger capacity"):
                runtime.commit_imp1_runtime(prepared, **_current(*items))
        self.assertEqual(
            before, (runtime._tracer_snapshot(tracers), tx.state, tx.causal_state)
        )

    def test_prepared_is_factory_only_and_has_no_endpoint_codec(self):
        with self.assertRaises(TypeError):
            runtime.TDG11IMP1PreparedRuntime()
        self.assertFalse(hasattr(runtime.TDG11IMP1CommittedRuntime, "as_mapping"))
        self.assertFalse(runtime.CAMPAIGN_EXECUTION_AUTHORIZED)
        self.assertFalse(runtime.FRESH_VS_RETRY_DURABLE_PROTOCOL_IMPLEMENTED)
        self.assertNotIn(
            "retry_successor",
            inspect.signature(runtime.prepare_imp1_initial_runtime).parameters,
        )
        retry_parameters = inspect.signature(
            runtime.prepare_imp1_retry_runtime
        ).parameters
        for forbidden in (
            "time",
            "requested_cap",
            "temporal_ledger",
            "method",
            "previous_step_index",
        ):
            self.assertNotIn(forbidden, retry_parameters)


if __name__ == "__main__":
    unittest.main()
