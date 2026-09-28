"""Genuine small synthetic HLT17 replay, mixed-owner, and invalid-premise controls."""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution import tdg11_c1r1_runtime as c1r1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_imp1_runtime as imp1  # noqa: E402
from recursive_horizons.fgc.evolution import hlt17_imp1_bridge as bridge  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_imp1_bridge import (  # noqa: E402
    FUTURE_INTEGRATION_SEAMS,
    HLT17BridgeResult,
    HLT17InvalidPremise,
    HLT17LiveBoundary,
    HLT17UnsupportedCallableBinding,
    PINNED_C1R1_PLAN,
    PINNED_C1R1_PREPARE,
    PINNED_C1R1_PRIVATE_DEPENDENCIES,
    PINNED_HLT17_PRIVATE_DEPENDENCIES,
    PINNED_IMP1_PRIVATE_DEPENDENCIES,
    authenticate_live_boundary,
    commit_hlt17_fine,
    commit_hlt17_kernel,
    finalize_hlt17_fine,
    prepare_hlt17_fresh,
    prepare_hlt17_next,
    replay_hlt17_rejection,
    require_hlt17_admission,
    seed_cursor_from_boundary,
    synthetic_callable_binding_digest,
    with_ledger,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    FRESH_READY,
    HLT17IMP1Error,
    RETRY_PENDING,
    encode_hlt17_cursor,
    imp1_enclosure_limits_for_member,
    restore_hlt17_cursor,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START,
    SYNTHETIC_TARGET,
    make_synthetic_fixture,
    synthetic_diagnostics,
    synthetic_exponential_rhs,
)


MODULE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_imp1_bridge.py"
START = SYNTHETIC_START
TARGET = SYNTHETIC_TARGET
REQUESTED_CAP = 1.0 / 8.0
FINE_DESCRIPTOR = sha256(b"hlt17-synthetic-fine-accepted-descriptor").hexdigest()


class ControllableRHS:
    """Same callable object; mode is test injection, not source configuration."""

    def __init__(self) -> None:
        self.mode = "impulse"

    def hlt17_synthetic_callable_binding(self) -> dict[str, object]:
        code = type(self).__call__.__code__
        return {
            "kind": "synthetic_failure_injection_double",
            "type": f"{type(self).__module__}.{type(self).__qualname__}",
            "call_code_sha256": sha256(code.co_code).hexdigest(),
            "call_names": list(code.co_names),
            "mutable_injection_mode_excluded": True,
            "physical_source_authenticated": False,
        }

    def __call__(self, time: float, current) -> EvolutionRHS:
        diagnostics = synthetic_diagnostics()
        zeros = np.zeros_like(current.u)
        if self.mode == "source":
            # Residual on the already-accepted instant is a PROTO7 terminal.
            if time > START:
                diagnostics["source_residual_infinity"] = 1.0
            return EvolutionRHS(zeros, zeros, zeros, diagnostics)
        if self.mode == "cfl":
            if time > START:
                diagnostics["coordinate_speed_upper"] = 1.0e9
            return EvolutionRHS(zeros, zeros, zeros, diagnostics)
        if self.mode == "ok":
            return synthetic_exponential_rhs(time, current)
        du = np.zeros_like(current.u)
        if time == START + 1.0 / 128.0:
            du[:, 4] = 1.0
        return EvolutionRHS(du, np.zeros_like(du), np.zeros_like(du), diagnostics)


def _with_defaults(func, defaults: tuple[object, ...]):
    cloned = types.FunctionType(
        func.__code__,
        func.__globals__,
        func.__name__,
        defaults,
        func.__closure__,
    )
    cloned.__kwdefaults__ = func.__kwdefaults__
    return cloned


def _descriptor(physical: str) -> str:
    digest = sha256(b"hlt17-synthetic-descriptor:" + physical.encode("ascii")).hexdigest()
    if digest == physical:
        digest = sha256(b"hlt17-synthetic-descriptor-alt").hexdigest()
    return digest


def _boundary(rhs: ControllableRHS | None = None, *, method: str = PRIMARY_METHOD):
    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(method)
    rhs = ControllableRHS() if rhs is None else rhs
    physical = array_content_sha256(state.u, state.p, state.q)
    member_key = "RK4-9" if method == PRIMARY_METHOD else "SSPRK3-9"
    return HLT17LiveBoundary(
        member_key=member_key,
        state=state,
        rhs=rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
        descriptor_sha256=_descriptor(physical),
        event_target=TARGET,
        limits=imp1_enclosure_limits_for_member(member_key),
    )


def _restore_pair(cursor, boundary):
    restored = restore_hlt17_cursor(encode_hlt17_cursor(cursor))
    live = with_ledger(boundary, restored.ledger)
    return restored, live


class HLT17BridgeTests(unittest.TestCase):
    def test_pinned_private_imp1_prepare_is_not_a_new_public_api(self) -> None:
        self.assertIs(PINNED_C1R1_PREPARE, c1r1._prepare)
        self.assertIsNot(PINNED_C1R1_PREPARE, imp1._prepare)
        self.assertNotIn("_prepare", c1r1.__all__)
        self.assertNotIn("_prepare", imp1.__all__)
        self.assertNotIn("prepare_imp1_initial_runtime", dict(PINNED_IMP1_PRIVATE_DEPENDENCIES))
        self.assertNotIn("tdg11_imp1_runtime._prepare", dict(PINNED_HLT17_PRIVATE_DEPENDENCIES))
        names = [name for name, _ in PINNED_HLT17_PRIVATE_DEPENDENCIES]
        self.assertIn("tdg11_c1r1_runtime._prepare", names)
        self.assertIn("tdg11_imp1_runtime._check_current_boundary", names)
        self.assertIn("tdg11_imp1_runtime._ledger_digest", names)
        self.assertIn("tdg6_temporal_admission_runtime._tracer_snapshot", names)
        self.assertIn("tdg6_temporal_admission_runtime._clone_tracers", names)
        self.assertIn("tdg6_temporal_admission_runtime._restore_tracers", names)
        for name, value in PINNED_C1R1_PRIVATE_DEPENDENCIES:
            _module_name, attr = name.rsplit(".", 1)
            self.assertTrue(callable(value))
            self.assertTrue(attr.startswith("_"))
        self.assertIn("member_codec", FUTURE_INTEGRATION_SEAMS)
        self.assertIn("content_addressed_store", FUTURE_INTEGRATION_SEAMS)
        self.assertIn("checkpoint_journal", FUTURE_INTEGRATION_SEAMS)
        self.assertIn("authenticated_origin", FUTURE_INTEGRATION_SEAMS)
        self.assertIn("store_implementation_identity", FUTURE_INTEGRATION_SEAMS)
        source = MODULE.read_text(encoding="utf-8")
        self.assertNotIn("from .hlt16", source)
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertIn("does not read a store", source.lower())
        self.assertIn("source_manifest", FUTURE_INTEGRATION_SEAMS)
        self.assertNotIn("imp1._prepare(", source)
        self.assertNotIn("prepare_imp1_initial_runtime(", source)
        self.assertNotIn("commit_imp1_runtime(", source)

    def test_genuine_prepare_nonadmit_persist_restore_replay_requested_cap_not_width(
        self,
    ) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        rhs.mode = "impulse"
        prepared = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        self.assertEqual(prepared.disposition, "prepared_fresh")
        assert prepared.prepared is not None
        self.assertEqual(prepared.prepared.plan.requested_cap, REQUESTED_CAP)
        self.assertEqual(prepared.prepared.plan.macro_width, 1.0 / 16.0)
        self.assertNotEqual(
            prepared.prepared.plan.requested_cap.hex(),
            prepared.prepared.plan.macro_width.hex(),
        )
        self.assertFalse(prepared.prepared.assessment.admission_passed)
        self.assertEqual(prepared.sink_writes, 0)

        closed = require_hlt17_admission(cursor, boundary, prepared.prepared)
        self.assertEqual(closed.disposition, "temporal_retry_required")
        self.assertEqual(closed.sink_writes, 1)
        self.assertEqual(closed.cursor.mode, RETRY_PENDING)
        assert closed.cursor.temporal is not None
        self.assertEqual(
            closed.cursor.temporal.predecessor_plan.requested_cap, REQUESTED_CAP
        )
        self.assertEqual(
            closed.cursor.temporal.imp1_half_cap, prepared.prepared.plan.macro_width / 2.0
        )
        inherited = cursor.ledger.inherited_snapshot
        debit = cursor.ledger.accumulated_debit_vector
        self.assertEqual(closed.cursor.ledger.inherited_snapshot, inherited)
        self.assertEqual(closed.cursor.ledger.accumulated_debit_vector, debit)

        restored, live = _restore_pair(closed.cursor, boundary)
        with patch.object(
            c1r1, "require_c1r1_admission", side_effect=AssertionError("replay sink")
        ), patch.object(
            imp1, "require_imp1_admission", side_effect=AssertionError("imp1 replay sink")
        ):
            replay = replay_hlt17_rejection(restored, live)
        self.assertEqual(replay.sink_writes, 0)
        self.assertEqual(
            replay.reconstructed.receipt_sha256,
            restored.temporal.preparation_sha256,
        )
        self.assertEqual(
            replay.reconstructed.plan.requested_cap.hex(), REQUESTED_CAP.hex()
        )
        self.assertNotEqual(
            replay.reconstructed.plan.requested_cap.hex(),
            replay.reconstructed.plan.macro_width.hex(),
        )
        self.assertEqual(
            replay.temporal_successor.plan, restored.temporal.successor_plan
        )
        self.assertEqual(
            replay.temporal_successor.plan.requested_cap,
            restored.temporal.imp1_half_cap,
        )

        with self.assertRaises(HLT17InvalidPremise):
            prepare_hlt17_fresh(restored, live, requested_cap=REQUESTED_CAP)

        rhs.mode = "impulse"
        successor = prepare_hlt17_next(restored, live, replay=replay)
        self.assertEqual(successor.disposition, "prepared_successor")
        assert successor.prepared is not None
        self.assertEqual(successor.prepared.plan, restored.temporal.successor_plan)
        self.assertEqual(successor.sink_writes, 0)

    def test_temporal_source_cfl_retry_across_serialize_restore(self) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        rhs.mode = "impulse"
        prepared = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        closed = require_hlt17_admission(cursor, boundary, prepared.prepared)
        temporal_count = closed.cursor.ledger.current_macro_step_temporal_retry_count
        cumulative = closed.cursor.ledger.cumulative_temporal_retry_count
        inherited = closed.cursor.ledger.inherited_snapshot
        debit = closed.cursor.ledger.accumulated_debit_vector
        half_cap = closed.cursor.temporal.successor_plan  # type: ignore[union-attr]
        restored, live = _restore_pair(closed.cursor, boundary)

        rhs.mode = "impulse"
        replay = replay_hlt17_rejection(restored, live)
        rhs.mode = "source"
        source = prepare_hlt17_next(restored, live, replay=replay)
        self.assertEqual(source.disposition, "source_retry_required")
        self.assertEqual(source.cursor.mode, RETRY_PENDING)
        self.assertEqual(source.cursor.latest_owner, "source")
        assert source.cursor.temporal is not None and source.cursor.overlay is not None
        self.assertEqual(source.cursor.temporal.successor_plan, half_cap)
        self.assertNotEqual(source.cursor.overlay.next_plan, half_cap)
        self.assertEqual(
            source.cursor.ledger.current_macro_step_temporal_retry_count, temporal_count
        )
        self.assertEqual(
            source.cursor.ledger.cumulative_temporal_retry_count, cumulative
        )
        self.assertEqual(source.cursor.ledger.inherited_snapshot, inherited)
        self.assertEqual(source.cursor.ledger.accumulated_debit_vector, debit)
        self.assertEqual(source.cursor.source_current, 1)
        self.assertEqual(source.cursor.source_total, 1)
        self.assertEqual(source.cursor.cfl_current, 0)
        self.assertEqual(source.sink_writes, 0)

        restored, live = _restore_pair(source.cursor, boundary)
        rhs.mode = "impulse"
        replay = replay_hlt17_rejection(restored, live)
        self.assertEqual(replay.temporal_successor.plan, half_cap)
        rhs.mode = "cfl"
        cfl = prepare_hlt17_next(restored, live, replay=replay)
        self.assertEqual(cfl.disposition, "cfl_retry_required")
        self.assertEqual(cfl.cursor.mode, RETRY_PENDING)
        self.assertEqual(cfl.cursor.latest_owner, "cfl")
        assert cfl.cursor.temporal is not None
        self.assertEqual(cfl.cursor.temporal.successor_plan, half_cap)
        self.assertEqual(
            cfl.cursor.ledger.current_macro_step_temporal_retry_count, temporal_count
        )
        self.assertEqual(cfl.cursor.source_current, 1)
        self.assertEqual(cfl.cursor.cfl_current, 1)
        self.assertEqual(cfl.cursor.source_total, 1)
        self.assertEqual(cfl.cursor.cfl_total, 1)
        self.assertEqual(cfl.cursor.ledger.inherited_snapshot, inherited)
        self.assertEqual(cfl.cursor.ledger.accumulated_debit_vector, debit)

        restored, live = _restore_pair(cfl.cursor, boundary)
        rhs.mode = "impulse"
        replay = replay_hlt17_rejection(restored, live)
        rhs.mode = "ok"
        nxt = prepare_hlt17_next(restored, live, replay=replay)
        self.assertEqual(nxt.disposition, "prepared_overlay")
        assert nxt.prepared is not None
        self.assertTrue(nxt.prepared.assessment.admission_passed)
        self.assertEqual(nxt.temporal_successor.plan, half_cap)
        self.assertNotEqual(nxt.prepared.plan, half_cap)
        admitted = require_hlt17_admission(restored, live, nxt.prepared)
        self.assertTrue(admitted.evidence["admission_passed"])
        committed = commit_hlt17_fine(
            restored, live, nxt.prepared, descriptor_sha256=FINE_DESCRIPTOR
        )
        self.assertEqual(committed.disposition, "accepted_fine")
        self.assertEqual(committed.cursor.mode, FRESH_READY)
        self.assertEqual(committed.cursor.source_current, 0)
        self.assertEqual(committed.cursor.cfl_current, 0)
        self.assertEqual(committed.cursor.source_total, 1)
        self.assertEqual(committed.cursor.cfl_total, 1)
        self.assertEqual(committed.cursor.ledger.current_macro_step_temporal_retry_count, 0)
        self.assertEqual(committed.cursor.ledger.cumulative_temporal_retry_count, cumulative)
        self.assertEqual(
            committed.cursor.ledger.serialized_temporal_rejections,
            restored.ledger.serialized_temporal_rejections,
        )
        self.assertEqual(committed.cursor.ledger.inherited_snapshot, inherited)
        self.assertEqual(committed.cursor.descriptor_sha256, FINE_DESCRIPTOR)
        self.assertNotEqual(committed.cursor.descriptor_sha256, restored.descriptor_sha256)
        self.assertTrue(committed.accepted_state_advanced)
        self.assertFalse(committed.evidence["physical_classification"])
        self.assertFalse(committed.evidence["corrected_endpoint_adopted"])
        assert committed.committed is not None
        self.assertEqual(committed.committed.time.hex(), committed.cursor.accepted_time.hex())
        self.assertEqual(
            committed.committed.step_index, committed.cursor.accepted_step_index
        )
        self.assertEqual(
            committed.committed.transaction_serial,
            committed.cursor.accepted_transaction_serial,
        )
        self.assertEqual(
            committed.committed.temporal_ledger.current_state_sha256,
            committed.cursor.physical_state_sha256,
        )

    def test_replay_source_failure_is_invalid_premise_and_spends_no_counters(
        self,
    ) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        rhs.mode = "impulse"
        prepared = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        closed = require_hlt17_admission(cursor, boundary, prepared.prepared)
        live = with_ledger(boundary, closed.cursor.ledger)
        before = (
            closed.cursor.source_current,
            closed.cursor.cfl_current,
            closed.cursor.ledger.current_macro_step_temporal_retry_count,
            closed.cursor.ledger.accumulated_debit_vector,
        )
        rhs.mode = "source"
        with self.assertRaises(HLT17InvalidPremise):
            replay_hlt17_rejection(closed.cursor, live)
        self.assertEqual(
            (
                closed.cursor.source_current,
                closed.cursor.cfl_current,
                closed.cursor.ledger.current_macro_step_temporal_retry_count,
                closed.cursor.ledger.accumulated_debit_vector,
            ),
            before,
        )
        rhs.mode = "impulse"
        with patch.object(
            bridge,
            "PINNED_C1R1_PREPARE",
            side_effect=c1r1.C1R1RuntimeResourceStop("maximum_owned_rows"),
        ):
            with self.assertRaises(HLT17InvalidPremise) as caught:
                replay_hlt17_rejection(closed.cursor, live)
        self.assertIn("resource", str(caught.exception))
        self.assertEqual(
            closed.cursor.ledger.current_macro_step_temporal_retry_count, before[2]
        )

    def test_swapped_descriptor_physical_and_altered_callbacks_fail_before_prepare(
        self,
    ) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        swapped = HLT17LiveBoundary(
            member_key=boundary.member_key,
            state=boundary.state,
            rhs=boundary.rhs,
            projector=boundary.projector,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            coordinates=boundary.coordinates,
            temporal_ledger=boundary.temporal_ledger,
            previous_step_index=boundary.previous_step_index,
            previous_transaction_serial=boundary.previous_transaction_serial,
            descriptor_sha256=cursor.physical_state_sha256,
            event_target=boundary.event_target,
            limits=boundary.limits,
        )
        with self.assertRaises(HLT17InvalidPremise):
            authenticate_live_boundary(cursor, swapped)

        class OtherRHS(ControllableRHS):
            pass

        other = HLT17LiveBoundary(
            member_key=boundary.member_key,
            state=boundary.state,
            rhs=OtherRHS(),
            projector=boundary.projector,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            coordinates=boundary.coordinates,
            temporal_ledger=boundary.temporal_ledger,
            previous_step_index=boundary.previous_step_index,
            previous_transaction_serial=boundary.previous_transaction_serial,
            descriptor_sha256=boundary.descriptor_sha256,
            event_target=boundary.event_target,
            limits=boundary.limits,
        )
        with self.assertRaises(HLT17InvalidPremise):
            authenticate_live_boundary(cursor, other)
        self.assertNotEqual(
            synthetic_callable_binding_digest(boundary.rhs, None),
            synthetic_callable_binding_digest(other.rhs, None),
        )

        retargeted = HLT17LiveBoundary(
            member_key=boundary.member_key,
            state=boundary.state,
            rhs=boundary.rhs,
            projector=boundary.projector,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            coordinates=boundary.coordinates,
            temporal_ledger=boundary.temporal_ledger,
            previous_step_index=boundary.previous_step_index,
            previous_transaction_serial=boundary.previous_transaction_serial,
            descriptor_sha256=boundary.descriptor_sha256,
            event_target=TARGET + 0.25,
            limits=boundary.limits,
        )
        with self.assertRaises(HLT17InvalidPremise):
            prepare_hlt17_fresh(cursor, retargeted, requested_cap=REQUESTED_CAP)

        calls: list[object] = []

        def forbidden(*args, **kwargs):
            calls.append((args, kwargs))
            raise AssertionError("prepare must not run after authentication failure")

        with patch.object(c1r1, "_prepare", side_effect=forbidden), patch.object(
            c1r1, "prepare_c1r1_initial_runtime", side_effect=forbidden
        ), patch.object(imp1, "_prepare", side_effect=forbidden), patch.object(
            imp1, "prepare_imp1_initial_runtime", side_effect=forbidden
        ):
            with self.assertRaises(HLT17InvalidPremise):
                authenticate_live_boundary(cursor, swapped)
        self.assertEqual(calls, [])

    def test_altered_preparation_rejection_and_ledger_fail_replay(self) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        rhs.mode = "impulse"
        prepared = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        closed = require_hlt17_admission(cursor, boundary, prepared.prepared)
        live = with_ledger(boundary, closed.cursor.ledger)
        mapping = encode_hlt17_cursor(closed.cursor)
        temporal = dict(mapping["temporal"])
        temporal["preparation_sha256"] = sha256(b"tampered-prep").hexdigest()
        mapping["temporal"] = temporal
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)

        mapping = encode_hlt17_cursor(closed.cursor)
        temporal = dict(mapping["temporal"])
        rejection = dict(temporal["rejection"])
        rejection["failed_channels"] = ["u:lambda"]
        temporal["rejection"] = rejection
        mapping["temporal"] = temporal
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)

        mapping = encode_hlt17_cursor(closed.cursor)
        ledger = dict(mapping["imp1_ledger"])
        ledger["origin_receipt_sha256"] = sha256(b"tampered-origin").hexdigest()
        mapping["imp1_ledger"] = ledger
        with self.assertRaises(HLT17IMP1Error):
            restore_hlt17_cursor(mapping)

        rhs.mode = "impulse"
        replay = replay_hlt17_rejection(closed.cursor, live)
        self.assertEqual(
            replay.reconstructed.assessment.assessment_sha256,
            closed.cursor.temporal.assessment_sha256,  # type: ignore[union-attr]
        )

    def test_fresh_source_overlay_without_temporal_evidence(self) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        rhs.mode = "source"
        result = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        self.assertEqual(result.disposition, "source_retry_required")
        self.assertEqual(result.cursor.mode, FRESH_READY)
        self.assertIsNone(result.cursor.temporal)
        self.assertIsNotNone(result.cursor.overlay)
        self.assertEqual(result.cursor.ledger.current_macro_step_temporal_retry_count, 0)
        self.assertFalse(result.cursor.public_fresh_permitted())
        with self.assertRaises(HLT17InvalidPremise):
            prepare_hlt17_fresh(result.cursor, with_ledger(boundary, result.cursor.ledger), requested_cap=REQUESTED_CAP)

    def test_synthetic_member_uses_truthful_nine_point_limits(self) -> None:
        boundary = _boundary()
        cursor = seed_cursor_from_boundary(boundary)
        limits = boundary.resolved_limits()
        self.assertEqual(boundary.member_key, "RK4-9")
        self.assertEqual(cursor.identity_scope, "synthetic")
        self.assertEqual(cursor.member.point_count, 9)
        self.assertEqual(limits.maximum_owned_rows, 4)
        self.assertEqual(limits.fallback_maximum_candidates_D01, 32)
        self.assertEqual(limits.fallback_maximum_candidates_D12, 64)
        self.assertEqual(boundary.state.shape[0], 9)
        production = imp1_enclosure_limits_for_member("RK4-2049")
        self.assertEqual(production.maximum_owned_rows, 2044)
        self.assertEqual(production.fallback_maximum_candidates_D01, 16352)
        self.assertEqual(production.fallback_maximum_candidates_D12, 32704)

    def test_production_member_refuses_nine_point_fixture(self) -> None:
        state, transaction, tracers, ledger, coordinates = make_synthetic_fixture()
        physical = array_content_sha256(state.u, state.p, state.q)
        with self.assertRaises(HLT17InvalidPremise):
            HLT17LiveBoundary(
                member_key="RK4-2049",
                state=state,
                rhs=ControllableRHS(),
                projector=None,
                transaction=transaction,
                tracers=tracers,
                coordinates=coordinates,
                temporal_ledger=ledger,
                descriptor_sha256=_descriptor(physical),
                event_target=TARGET,
                limits=imp1_enclosure_limits_for_member("RK4-2049"),
            )

    def test_nested_rhs_captured_rate_refuses_substitution_before_prepare(self) -> None:
        def make_rhs(rate: float):
            def rhs(time: float, state) -> EvolutionRHS:
                return EvolutionRHS(
                    rate * state.u,
                    rate * state.p,
                    rate * state.q,
                    synthetic_diagnostics(),
                )

            return rhs

        first = make_rhs(1.0)
        second = make_rhs(2.0)
        self.assertNotEqual(
            first.__code__.co_code, b""
        )
        self.assertEqual(first.__code__.co_code, second.__code__.co_code)
        self.assertNotEqual(
            synthetic_callable_binding_digest(first, None),
            synthetic_callable_binding_digest(second, None),
        )
        boundary = _boundary(first)
        cursor = seed_cursor_from_boundary(boundary)
        swapped = HLT17LiveBoundary(
            member_key=boundary.member_key,
            state=boundary.state,
            rhs=second,
            projector=None,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            coordinates=boundary.coordinates,
            temporal_ledger=boundary.temporal_ledger,
            descriptor_sha256=boundary.descriptor_sha256,
            event_target=boundary.event_target,
            limits=boundary.limits,
        )
        with self.assertRaises(HLT17InvalidPremise):
            authenticate_live_boundary(cursor, swapped)
        calls: list[object] = []

        def forbidden(*args, **kwargs):
            calls.append(1)
            raise AssertionError("prepare must not run after rate substitution")

        with patch.object(
            c1r1, "prepare_c1r1_initial_runtime", side_effect=forbidden
        ), patch.object(imp1, "prepare_imp1_initial_runtime", side_effect=forbidden):
            with self.assertRaises(HLT17InvalidPremise):
                prepare_hlt17_fresh(cursor, swapped, requested_cap=REQUESTED_CAP)
        self.assertEqual(calls, [])

    def test_defaults_consts_and_frozen_callable_parameters_are_bound(self) -> None:
        def rhs_default(time: float, state, rate=1.0) -> EvolutionRHS:
            return EvolutionRHS(
                rate * state.u, rate * state.p, rate * state.q, synthetic_diagnostics()
            )

        self.assertNotEqual(
            synthetic_callable_binding_digest(rhs_default, None),
            synthetic_callable_binding_digest(
                _with_defaults(rhs_default, (2.0,)), None
            ),
        )

        def rhs_const_one(time: float, state) -> EvolutionRHS:
            rate = 1.0
            return EvolutionRHS(
                rate * state.u, rate * state.p, rate * state.q, synthetic_diagnostics()
            )

        def rhs_const_two(time: float, state) -> EvolutionRHS:
            rate = 2.0
            return EvolutionRHS(
                rate * state.u, rate * state.p, rate * state.q, synthetic_diagnostics()
            )

        self.assertEqual(rhs_const_one.__code__.co_code, rhs_const_two.__code__.co_code)
        self.assertNotEqual(
            synthetic_callable_binding_digest(rhs_const_one, None),
            synthetic_callable_binding_digest(rhs_const_two, None),
        )

        @dataclass(frozen=True, slots=True)
        class ScaledRHS:
            rate: float

            def __call__(self, time: float, state) -> EvolutionRHS:
                return EvolutionRHS(
                    self.rate * state.u,
                    self.rate * state.p,
                    self.rate * state.q,
                    synthetic_diagnostics(),
                )

        self.assertNotEqual(
            synthetic_callable_binding_digest(ScaledRHS(1.0), None),
            synthetic_callable_binding_digest(ScaledRHS(2.0), None),
        )

        class Opaque:
            def __call__(self, time: float, state) -> EvolutionRHS:
                return synthetic_exponential_rhs(time, state)

        with self.assertRaises(HLT17UnsupportedCallableBinding):
            synthetic_callable_binding_digest(Opaque(), None)

        class MutableConfigured:
            def __init__(self, rate: float) -> None:
                self.rate = rate

            def __call__(self, time: float, state) -> EvolutionRHS:
                return EvolutionRHS(
                    self.rate * state.u,
                    self.rate * state.p,
                    self.rate * state.q,
                    synthetic_diagnostics(),
                )

        with self.assertRaises(HLT17UnsupportedCallableBinding):
            synthetic_callable_binding_digest(MutableConfigured(1.0), None)

    def test_ssprk3_9_synthetic_retry_restore_and_mixed_overlay(self) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs, method=COMPARATOR_METHOD)
        cursor = seed_cursor_from_boundary(boundary)
        self.assertEqual(cursor.member_key, "SSPRK3-9")
        self.assertEqual(cursor.identity_scope, "synthetic")
        self.assertEqual(cursor.member.point_count, 9)
        rhs.mode = "impulse"
        prepared = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        closed = require_hlt17_admission(cursor, boundary, prepared.prepared)
        restored, live = _restore_pair(closed.cursor, boundary)
        self.assertEqual(restored.member_key, "SSPRK3-9")
        rhs.mode = "impulse"
        replay = replay_hlt17_rejection(restored, live)
        rhs.mode = "source"
        source = prepare_hlt17_next(restored, live, replay=replay)
        self.assertEqual(source.disposition, "source_retry_required")
        self.assertEqual(source.cursor.mode, RETRY_PENDING)
        self.assertEqual(source.cursor.member_key, "SSPRK3-9")
        restored, live = _restore_pair(source.cursor, boundary)
        self.assertEqual(restored.overlay.owner, "source")
        self.assertEqual(restored.identity_scope, "synthetic")

    def test_bridge_result_replace_cannot_claim_success_without_payloads(self) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        invalid = HLT17BridgeResult(
            disposition="invalid_premise",
            cursor=cursor,
            evidence={"kind": "invalid_premise"},
        )
        with self.assertRaises(HLT17IMP1Error):
            replace(invalid, accepted_state_advanced=True)
        with self.assertRaises(HLT17IMP1Error):
            replace(invalid, disposition="accepted_fine", accepted_state_advanced=True)
        with self.assertRaises(HLT17IMP1Error):
            replace(invalid, disposition="prepared_fresh")
        with self.assertRaises(HLT17IMP1Error):
            replace(invalid, disposition="replay_authenticated")
        exhausted = HLT17BridgeResult(
            disposition="temporal_retry_exhausted",
            cursor=cursor,
            evidence={
                "kind": "temporal_retry_exhausted",
                "updated_imp1_ledger": {"schema": "x"},
                "executable_retry_cursor": False,
            },
        )
        with self.assertRaises(HLT17IMP1Error):
            replace(exhausted, evidence={"kind": "temporal_retry_exhausted"})
        with self.assertRaises(HLT17IMP1Error):
            replace(invalid, disposition="accepted_fine", accepted_state_advanced=True, prepared=object())  # type: ignore[arg-type]


def _runtime_bits(boundary: HLT17LiveBoundary) -> tuple[object, ...]:
    return bridge._runtime_fingerprint(boundary.transaction, boundary.tracers)


def _doubled_rhs(time: float, state) -> EvolutionRHS:
    return EvolutionRHS(
        2.0 * state.u, 2.0 * state.p, 2.0 * state.q, synthetic_diagnostics()
    )


def _identity_projector(_time: float, state):
    return state


def _foreign_prepared(boundary: HLT17LiveBoundary, **changes):
    kwargs = dict(
        method=boundary.temporal_ledger.method,
        time=START,
        event_target=boundary.event_target,
        requested_cap=REQUESTED_CAP,
        state=boundary.state,
        rhs=boundary.rhs,
        projector=boundary.projector,
        transaction=boundary.transaction,
        tracers=boundary.tracers,
        coordinates=boundary.coordinates,
        temporal_ledger=boundary.temporal_ledger,
        previous_step_index=boundary.previous_step_index,
        previous_transaction_serial=boundary.previous_transaction_serial,
        limits=boundary.resolved_limits(),
    )
    kwargs.update(changes)
    return c1r1.prepare_c1r1_initial_runtime(**kwargs)


class HLT17PreparedBindingTests(unittest.TestCase):
    def test_wrong_prepared_rhs_is_rejected_before_sink_or_commit(self) -> None:
        boundary = _boundary(synthetic_exponential_rhs)
        cursor = seed_cursor_from_boundary(boundary)
        wrong = _foreign_prepared(boundary, rhs=_doubled_rhs)
        self.assertTrue(wrong.assessment.admission_passed)
        self.assertIsNot(wrong.original_rhs, boundary.rhs)
        with patch.object(c1r1, "require_c1r1_admission") as require, patch.object(
            c1r1, "commit_c1r1_runtime"
        ) as commit:
            with self.assertRaises(HLT17InvalidPremise):
                require_hlt17_admission(cursor, boundary, wrong)
            with self.assertRaises(HLT17InvalidPremise):
                commit_hlt17_fine(
                    cursor, boundary, wrong, descriptor_sha256=FINE_DESCRIPTOR
                )
            require.assert_not_called()
            commit.assert_not_called()
        self.assertEqual(boundary.transaction.state.accepted_stage_count, 0)

    def test_wrong_projector_and_event_target_are_rejected_before_kernel(self) -> None:
        boundary = _boundary(synthetic_exponential_rhs)
        cursor = seed_cursor_from_boundary(boundary)
        wrong_proj = _foreign_prepared(boundary, projector=_identity_projector)
        with patch.object(c1r1, "commit_c1r1_runtime") as commit:
            with self.assertRaises(HLT17InvalidPremise):
                commit_hlt17_fine(
                    cursor, boundary, wrong_proj, descriptor_sha256=FINE_DESCRIPTOR
                )
            commit.assert_not_called()
        wrong_target = _foreign_prepared(boundary, event_target=TARGET + 0.25)
        with patch.object(c1r1, "require_c1r1_admission") as require:
            with self.assertRaises(HLT17InvalidPremise):
                require_hlt17_admission(cursor, boundary, wrong_target)
            require.assert_not_called()

    def test_predecessor_ledger_and_foreign_retry_plan_are_rejected(self) -> None:
        rhs = ControllableRHS()
        boundary = _boundary(rhs)
        cursor = seed_cursor_from_boundary(boundary)
        rhs.mode = "impulse"
        prepared = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        closed = require_hlt17_admission(cursor, boundary, prepared.prepared)
        live = with_ledger(boundary, closed.cursor.ledger)
        predecessor = closed.cursor.temporal.predecessor_ledger  # type: ignore[union-attr]
        stale = PINNED_C1R1_PREPARE(
            plan=closed.cursor.temporal.predecessor_plan,  # type: ignore[union-attr]
            method=closed.cursor.method,
            state=live.state,
            rhs=rhs,
            projector=None,
            transaction=live.transaction,
            tracers=live.tracers,
            coordinates=live.coordinates,
            temporal_ledger=predecessor,
            previous_step_index=live.previous_step_index,
            previous_transaction_serial=live.previous_transaction_serial,
            limits=live.resolved_limits(),
        )
        rhs.mode = "ok"
        with patch.object(c1r1, "require_c1r1_admission") as require, patch.object(
            c1r1, "commit_c1r1_runtime"
        ) as commit:
            with self.assertRaises(HLT17InvalidPremise):
                require_hlt17_admission(closed.cursor, live, stale)
            with self.assertRaises(HLT17InvalidPremise):
                commit_hlt17_fine(
                    closed.cursor, live, stale, descriptor_sha256=FINE_DESCRIPTOR
                )
            require.assert_not_called()
            commit.assert_not_called()

        next_plan = closed.cursor.next_attempt_plan()
        assert next_plan is not None
        other_plan = PINNED_C1R1_PLAN(START, TARGET, REQUESTED_CAP)
        self.assertNotEqual(other_plan, next_plan)
        foreign = PINNED_C1R1_PREPARE(
            plan=other_plan,
            method=closed.cursor.method,
            state=live.state,
            rhs=rhs,
            projector=None,
            transaction=live.transaction,
            tracers=live.tracers,
            coordinates=live.coordinates,
            temporal_ledger=live.temporal_ledger,
            previous_step_index=live.previous_step_index,
            previous_transaction_serial=live.previous_transaction_serial,
            limits=live.resolved_limits(),
        )
        with patch.object(c1r1, "commit_c1r1_runtime") as commit:
            with self.assertRaises(HLT17InvalidPremise):
                commit_hlt17_fine(
                    closed.cursor, live, foreign, descriptor_sha256=FINE_DESCRIPTOR
                )
            commit.assert_not_called()

    def test_late_adoption_failures_roll_back_complete_runtime_boundary(self) -> None:
        boundary = _boundary(synthetic_exponential_rhs)
        cursor = seed_cursor_from_boundary(boundary)
        prepared = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        admitted = require_hlt17_admission(cursor, boundary, prepared.prepared)
        self.assertTrue(admitted.evidence["admission_passed"])
        family = admitted.prepared
        assert family is not None
        before = _runtime_bits(boundary)
        self.assertEqual(boundary.transaction.state.accepted_stage_count, 0)

        for target, error in (
            (bridge, "live_boundary_identities"),
            (bridge, "close_fine_adoption"),
            (bridge, "HLT17BridgeResult"),
        ):
            with self.subTest(failure=error):
                with patch.object(
                    target, error, side_effect=RuntimeError(f"late {error}")
                ):
                    with self.assertRaises(RuntimeError):
                        commit_hlt17_fine(
                            cursor,
                            boundary,
                            family,
                            descriptor_sha256=FINE_DESCRIPTOR,
                        )
                self.assertEqual(_runtime_bits(boundary), before)
                self.assertEqual(boundary.transaction.state.accepted_stage_count, 0)

    def test_committed_runtime_is_returned_and_matches_cursor(self) -> None:
        boundary = _boundary(synthetic_exponential_rhs)
        cursor = seed_cursor_from_boundary(boundary)
        prepared = prepare_hlt17_fresh(cursor, boundary, requested_cap=REQUESTED_CAP)
        admitted = require_hlt17_admission(cursor, boundary, prepared.prepared)
        committed_runtime = commit_hlt17_kernel(cursor, boundary, admitted.prepared)
        result = finalize_hlt17_fine(
            cursor,
            boundary,
            admitted.prepared,
            committed_runtime,
            descriptor_sha256=FINE_DESCRIPTOR,
        )
        self.assertEqual(result.disposition, "accepted_fine")
        self.assertEqual(result.cursor.descriptor_sha256, FINE_DESCRIPTOR)
        self.assertNotEqual(result.cursor.descriptor_sha256, cursor.descriptor_sha256)
        committed = result.committed
        self.assertIs(committed, committed_runtime)
        assert committed is not None
        self.assertIsInstance(committed, c1r1.TDG11C1R1CommittedRuntime)
        self.assertEqual(
            committed.implementation_id,
            result.cursor.implementation_id,
        )
        self.assertEqual(committed.time.hex(), result.cursor.accepted_time.hex())
        self.assertEqual(committed.step_index, result.cursor.accepted_step_index)
        self.assertEqual(
            committed.transaction_serial, result.cursor.accepted_transaction_serial
        )
        self.assertEqual(
            committed.temporal_ledger.current_state_sha256,
            result.cursor.physical_state_sha256,
        )
        self.assertEqual(
            committed.temporal_ledger.current_step_index,
            result.cursor.accepted_step_index,
        )
        self.assertGreater(boundary.transaction.state.accepted_stage_count, 0)
        self.assertEqual(
            boundary.transaction.state.accepted_stage_count,
            committed.transaction_serial,
        )


if __name__ == "__main__":
    unittest.main()
