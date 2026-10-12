"""Focused classification tests for the one-proposal HLT16 adapter."""
from __future__ import annotations

import ast
import copy
import importlib
import inspect
import json
import unittest
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import Mock, patch

from recursive_horizons.fgc.evolution import proto15_runtime as p15
from recursive_horizons.fgc.evolution import hlt16_progression_attempt as attempt
from recursive_horizons.fgc.evolution.proto14_runtime import Proto14RunMember
from tests import test_hlt16_lifecycle as lifecycle_fixture


DESCRIPTOR = "2" * 64


def ledger() -> p15.TDG6TemporalLedger:
    return p15.TDG6TemporalLedger.zero(initial_time=23 / 16)


def cursor(value: p15.TDG6TemporalLedger | None = None) -> p15.Proto15Cursor:
    value = ledger() if value is None else value
    return p15._cursor({
        "protocol_artifact_id": "FGC-2-SF1-PROTO18", "campaign_id": "A",
        "member_key": "RK4-2049", "method": "RK4", "point_count": 2049,
        "committed_common_event_index": 23,
        "accepted_boundary_time": {"rational": "23/16", "binary64_hex": (23 / 16).hex()},
        "accepted_state_sha256": DESCRIPTOR, "previous_step_index": 1,
        "previous_transaction_serial": 3, "TDG6_ledger_sha256": p15._hash(p15._ledger_mapping(value)),
        "cursor_generation": 0, "event_target_time": {"rational": "3/2", "binary64_hex": (3 / 2).hex()},
        "attempt_serial": 0, "attempt_id": "A:RK4-2049:23:0", "mode": "FRESH_READY",
        "retry_successor_payload_or_none": None, "journal_tip_sha256": p15.ROOT_SHA256,
        "cursor_chain_parent_sha256": p15.ROOT_SHA256,
    })


def member(value: p15.TDG6TemporalLedger | None = None) -> Proto14RunMember:
    value = ledger() if value is None else value
    result = object.__new__(Proto14RunMember)
    result.method_label = "RK4"; result.integrator_id = "RK4"; result.point_count = 2049
    result.time = 23 / 16; result.temporal_ledger = value; result.state = "old"
    result.step_index = 1; result.transaction_serial = 3
    result.operator = object(); result.projector = None; result.transaction = object(); result.tracers = object()
    result.initial = SimpleNamespace(grid=SimpleNamespace(coordinates=object()))
    result.snapshot = lambda: {"state": result.state, "time": result.time, "ledger": result.temporal_ledger}
    def restore(boundary):
        result.state = boundary["state"]; result.time = boundary["time"]; result.temporal_ledger = boundary["ledger"]
        result.restored = True
    result.restore = restore
    result._accepted_boundary_preserved = lambda boundary: (
        result.state == boundary["state"] and result.time == boundary["time"] and result.temporal_ledger == boundary["ledger"]
    )
    return result


class FakeTemporalRequired(Exception):
    def __init__(self):
        self.retry_step_size = 1 / 32
        self.evidence = SimpleNamespace(as_mapping=lambda: {"event_type": "retry"})


class FakeTemporalExhausted(Exception):
    def __init__(self):
        self.reason = "minimum_macro_step"
        self.evidence = SimpleNamespace(as_mapping=lambda: {"event_type": "retry"})


class FakePathStop(Exception):
    def __init__(self, *, source_retry=None, cause=None, path="fine"):
        self.source_retry, self.cause, self.path = source_retry, cause, path


class FakeCFL(Exception):
    def __init__(self):
        self.observed_ratio, self.maximum_ratio = 2.0, 1.0


class FakeScientific(Exception):
    pass


class TestHLT16ProgressionAttempt(unittest.TestCase):
    def test_no_loop_persistence_or_candidate_surface(self) -> None:
        module = importlib.import_module("recursive_horizons.fgc.evolution.hlt16_progression_attempt")
        source = ast.unparse(ast.parse(inspect.getsource(module)))
        for forbidden in ("while ", "Path(", "open(", "FGCQR", "SGBL", "numpy"):
            self.assertNotIn(forbidden, source)

    def test_accepted_fine_mutates_only_in_memory_and_returns_complete_evidence(self) -> None:
        item = member(); old = cursor(); prepared = SimpleNamespace(plan="PLAN")
        committed = SimpleNamespace(state="new", time=1.5, step_index=2, transaction_serial=4,
                                    temporal_ledger=ledger(), continuous_admission="admitted")
        with patch.object(attempt, "_prepare_initial", return_value=prepared), \
             patch.object(attempt.tdg7, "require_tdg7_stage_safe_admission"), \
             patch.object(attempt.tdg7, "commit_tdg7_stage_safe_runtime", return_value=committed):
            result = attempt.attempt_once(item, old, 1.5, descriptor_sha256=DESCRIPTOR, requested_cap=1 / 16)
        self.assertEqual(result.disposition, "accepted_fine")
        self.assertEqual(item.state, "new")
        self.assertEqual(result.accepted_time, 1.5)
        self.assertFalse(result.accepted_boundary_restored)

    def test_temporal_outcomes_restore_boundary_and_preserve_counters(self) -> None:
        for raised, expected in ((FakeTemporalRequired(), "temporal_retry_required"), (FakeTemporalExhausted(), "temporal_retry_exhausted")):
            with self.subTest(expected=expected):
                item = member(); old = cursor(); prepared = SimpleNamespace(plan="PLAN")
                with patch.object(attempt, "_prepare_initial", return_value=prepared), \
                     patch.object(attempt.tdg6, "TDG6TemporalRetryRequired", FakeTemporalRequired), \
                     patch.object(attempt.tdg6, "TDG6TemporalRetryExhausted", FakeTemporalExhausted), \
                     patch.object(attempt.tdg7, "require_tdg7_stage_safe_admission", side_effect=raised):
                    result = attempt.attempt_once(item, old, 1.5, descriptor_sha256=DESCRIPTOR, requested_cap=1 / 16)
                self.assertEqual(result.disposition, expected)
                self.assertTrue(result.accepted_boundary_restored)
                self.assertEqual(item.state, "old")

    def test_pending_cursor_uses_replayed_successor_then_retry_entrypoint(self) -> None:
        item = member()
        pending_cursor = SimpleNamespace(
            mode="RETRY_PENDING",
            payload={"retry_successor_payload_or_none": {"half_cap": (1 / 32).hex()}},
        )
        replayed = SimpleNamespace(plan="REPLAYED")
        prepared = SimpleNamespace(plan="RETRY-PLAN")
        committed = SimpleNamespace(state="retry-new", time=1.5, step_index=2, transaction_serial=4,
                                    temporal_ledger=ledger(), continuous_admission="admitted")
        with patch.object(attempt, "_cursor_and_member_preflight"), \
             patch.object(attempt, "_replay_pending_rejection", return_value=replayed) as replay, \
             patch.object(attempt.tdg7, "prepare_tdg7_retry_successor", return_value=prepared) as prepare_retry, \
             patch.object(attempt.tdg7, "require_tdg7_stage_safe_admission"), \
             patch.object(attempt.tdg7, "commit_tdg7_stage_safe_runtime", return_value=committed):
            with patch.object(attempt, "_same_binary64", return_value=True):
                result = attempt.attempt_once(item, pending_cursor, 1.5, descriptor_sha256=DESCRIPTOR, requested_cap=1 / 32)  # type: ignore[arg-type]
        self.assertEqual(result.disposition, "accepted_fine")
        replay.assert_called_once_with(item, pending_cursor, 1.5)
        self.assertEqual(prepare_retry.call_args.kwargs["successor"], replayed)

    def _pending_replay_fixture(self, *, cap=1 / 16):
        old, prior = lifecycle_fixture.proto18()
        predecessor = lifecycle_fixture.plan(old, cap)
        retry = lifecycle_fixture.rejected(old, prior, predecessor)
        self.assertIsInstance(retry, p15.TDG6TemporalRetryRequired)
        successor = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(
                retry.evidence.initial_time,
                3 / 2,
                retry.retry_step_size,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
        )
        pending = lifecycle_fixture.close_retry(
            old,
            prior,
            retry,
            predecessor,
            successor,
            lifecycle_fixture.REJECTION,
            p15.ROOT_SHA256,
        )
        item = member(retry.updated_ledger)
        item.state = SimpleNamespace(u=object(), p=object(), q=object())
        prepared = SimpleNamespace(plan=predecessor)
        replayed = SimpleNamespace(plan=successor)
        return item, pending, retry, prepared, replayed

    def _pending_replay_two_fixture(self):
        old, base = lifecycle_fixture.proto18()
        first_plan = lifecycle_fixture.plan(old)
        first = lifecycle_fixture.rejected(old, base, first_plan)
        self.assertIsInstance(first, p15.TDG6TemporalRetryRequired)
        first_successor = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(
                first.evidence.initial_time,
                3 / 2,
                first.retry_step_size,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
        )
        pending_one = lifecycle_fixture.close_retry(
            old,
            base,
            first,
            first_plan,
            first_successor,
            lifecycle_fixture.REJECTION,
            p15.ROOT_SHA256,
        )
        second = lifecycle_fixture.rejected(
            pending_one,
            first.updated_ledger,
            first_successor,
        )
        self.assertIsInstance(second, p15.TDG6TemporalRetryRequired)
        second_successor = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(
                second.evidence.initial_time,
                3 / 2,
                second.retry_step_size,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
        )
        pending_two = lifecycle_fixture.close_retry(
            pending_one,
            first.updated_ledger,
            second,
            first_successor,
            second_successor,
            "4" * 64,
            lifecycle_fixture.REJECTION,
        )
        item = member(second.updated_ledger)
        item.state = SimpleNamespace(u=object(), p=object(), q=object())
        prepared = SimpleNamespace(plan=first_successor)
        replayed = SimpleNamespace(plan=second_successor)
        return item, pending_two, first, second, prepared, replayed

    def _pending_replay_three_with_older_prefix_fixture(self):
        """Build one accepted old rejection plus a three-record active suffix."""
        old, zero = lifecycle_fixture.proto18()

        old_plan = lifecycle_fixture.plan(old, 1 / 128)
        old_retry = lifecycle_fixture.rejected(old, zero, old_plan)
        self.assertIsInstance(old_retry, p15.TDG6TemporalRetryRequired)
        old_successor = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(
                old_retry.evidence.initial_time,
                3 / 2,
                old_retry.retry_step_size,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
        )
        old_pending = lifecycle_fixture.close_retry(
            old,
            zero,
            old_retry,
            old_plan,
            old_successor,
            "4" * 64,
            p15.ROOT_SHA256,
        )
        old_endpoint = float.fromhex(old_successor["boundaries_hex"][-1])
        base = lifecycle_fixture.accepted(old_retry.updated_ledger, old_endpoint)
        accepted_descriptor = "5" * 64
        fresh = lifecycle_fixture.close_fine(
            old_pending,
            old_retry.updated_ledger,
            base,
            accepted_descriptor,
            old_successor,
            lifecycle_fixture.time_identity(old_endpoint),
            8,
            12,
            lifecycle_fixture.CHECKPOINT_TIP,
        )

        active: list[p15.TDG6TemporalRetryRequired] = []
        current = fresh
        current_ledger = base
        current_plan = lifecycle_fixture.plan(fresh)
        predecessor_tip = lifecycle_fixture.CHECKPOINT_TIP
        for index, durable in enumerate(("6" * 64, "7" * 64, "8" * 64)):
            outcome = lifecycle_fixture.rejected(
                current, current_ledger, current_plan
            )
            self.assertIsInstance(outcome, p15.TDG6TemporalRetryRequired)
            assert isinstance(outcome, p15.TDG6TemporalRetryRequired)
            successor = p15._plan_mapping(
                p15.plan_forward_proto14_subdivision(
                    outcome.evidence.initial_time,
                    3 / 2,
                    outcome.retry_step_size,
                    minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
                )
            )
            current = lifecycle_fixture.close_retry(
                current,
                current_ledger,
                outcome,
                current_plan,
                successor,
                durable,
                predecessor_tip,
            )
            current_ledger = outcome.updated_ledger
            current_plan = successor
            predecessor_tip = durable
            active.append(outcome)

        item = member(current_ledger)
        item.time = old_endpoint
        item.step_index = 8
        item.transaction_serial = 12
        item.state = SimpleNamespace(u=object(), p=object(), q=object())
        nested = current.payload["retry_successor_payload_or_none"]
        prepared = SimpleNamespace(plan=nested["predecessor_plan"])
        replayed = SimpleNamespace(plan=nested["successor_plan"])
        return (
            item,
            current,
            base,
            tuple(active),
            prepared,
            replayed,
            accepted_descriptor,
        )

    def test_pending_replay_normalizes_only_tuple_list_persistence_shape(self) -> None:
        item, pending, retry, prepared, replayed = self._pending_replay_fixture()

        def reject(_prepared, **kwargs):
            kwargs["durable_rejection_sink"](retry.evidence.as_mapping())
            raise retry

        with patch.object(attempt, "array_content_sha256", return_value=lifecycle_fixture.STATE), \
             patch.object(attempt, "_prepare_initial", return_value=prepared), \
             patch.object(attempt.tdg7, "require_tdg7_stage_safe_admission", side_effect=reject), \
             patch.object(attempt.tdg7, "replan_tdg7_after_tdg6_retry", return_value=replayed):
            first = attempt._replay_pending_rejection(item, pending, 1.5)
            second = attempt._replay_pending_rejection(item, pending, 1.5)
        self.assertIs(first, replayed)
        self.assertIs(second, replayed)
        self.assertEqual(item.temporal_ledger, retry.updated_ledger)
        self.assertEqual(item.time, 23 / 16)
        self.assertEqual(item.step_index, 1)
        self.assertEqual(item.transaction_serial, 3)

    def test_pending_replay_preserves_requested_cap_when_selected_width_is_smaller(self) -> None:
        item, pending, retry, prepared, replayed = self._pending_replay_fixture(
            cap=1 / 8
        )
        predecessor = pending.payload["retry_successor_payload_or_none"][
            "predecessor_plan"
        ]
        self.assertEqual(float.fromhex(predecessor["requested_cap_hex"]), 1 / 8)
        self.assertEqual(float.fromhex(predecessor["macro_width_hex"]), 1 / 16)

        def reject(_prepared, **kwargs):
            kwargs["durable_rejection_sink"](retry.evidence.as_mapping())
            raise retry

        with patch.object(attempt, "array_content_sha256", return_value=lifecycle_fixture.STATE), \
             patch.object(attempt, "_prepare_initial", return_value=prepared) as prepare, \
             patch.object(attempt.tdg7, "require_tdg7_stage_safe_admission", side_effect=reject), \
             patch.object(attempt.tdg7, "replan_tdg7_after_tdg6_retry", return_value=replayed):
            self.assertIs(
                attempt._replay_pending_rejection(item, pending, 1.5),
                replayed,
            )
        self.assertEqual(prepare.call_args.args[3], 1 / 8)

    def test_second_pending_replay_reruns_only_the_latest_exact_predecessor(self) -> None:
        item, pending, first, second, prepared, replayed = (
            self._pending_replay_two_fixture()
        )

        def reject(_prepared, **kwargs):
            kwargs["durable_rejection_sink"](second.evidence.as_mapping())
            raise second

        with patch.object(attempt, "array_content_sha256", return_value=lifecycle_fixture.STATE), \
             patch.object(attempt, "_prepare_initial") as fresh, \
             patch.object(
                 attempt.tdg8_replay,
                 "prepare_tdg8_persisted_retry_replay",
                 return_value=prepared,
             ) as replay_prepare, \
             patch.object(
                 attempt.tdg8_replay,
                 "require_tdg8_persisted_retry_replay_admission",
                 side_effect=reject,
             ) as admission, \
             patch.object(
                 attempt.tdg8_replay,
                 "replan_tdg8_after_persisted_retry_rejection",
                 return_value=replayed,
             ):
            result = attempt._replay_pending_rejection(item, pending, 1.5)

        self.assertIs(result, replayed)
        fresh.assert_not_called()
        admission.assert_called_once()
        self.assertEqual(
            replay_prepare.call_args.kwargs["expected_prior_retry_count"], 1
        )
        self.assertEqual(
            replay_prepare.call_args.kwargs["temporal_ledger"],
            first.updated_ledger,
        )
        self.assertEqual(item.temporal_ledger, second.updated_ledger)
        self.assertEqual(item.time, 23 / 16)
        self.assertEqual(item.step_index, 1)
        self.assertEqual(item.transaction_serial, 3)

    def test_second_pending_cursor_reaches_the_third_tdg6_disposition(self) -> None:
        item, pending, _first, second, replay_prepared, replayed = (
            self._pending_replay_two_fixture()
        )
        nested = pending.payload["retry_successor_payload_or_none"]
        third_plan = nested["successor_plan"]
        third = lifecycle_fixture.rejected(
            pending,
            second.updated_ledger,
            third_plan,
        )
        self.assertIsInstance(third, p15.TDG6TemporalRetryRequired)
        current_prepared = SimpleNamespace(plan=third_plan)

        def reject_second(_prepared, **kwargs):
            kwargs["durable_rejection_sink"](second.evidence.as_mapping())
            raise second

        def reject_third(_prepared, **kwargs):
            kwargs["durable_rejection_sink"](third.evidence.as_mapping())
            raise third

        with patch.object(attempt, "array_content_sha256", return_value=lifecycle_fixture.STATE), \
             patch.object(
                 attempt.tdg8_replay,
                 "prepare_tdg8_persisted_retry_replay",
                 return_value=replay_prepared,
             ), \
             patch.object(
                 attempt.tdg7,
                 "prepare_tdg7_retry_successor",
                 return_value=current_prepared,
             ), \
             patch.object(
                 attempt.tdg8_replay,
                 "require_tdg8_persisted_retry_replay_admission",
                 side_effect=reject_second,
             ) as replay_admission, \
             patch.object(
                 attempt.tdg8_replay,
                 "replan_tdg8_after_persisted_retry_rejection",
                 return_value=replayed,
             ), \
             patch.object(
                 attempt.tdg7,
                 "require_tdg7_stage_safe_admission",
                 side_effect=reject_third,
             ) as current_admission:
            result = attempt.attempt_once(
                item,
                pending,
                1.5,
                descriptor_sha256=DESCRIPTOR,
                requested_cap=float.fromhex(nested["half_cap"]),
            )

        self.assertEqual(result.disposition, "temporal_retry_required")
        self.assertIs(result.retry_outcome, third)
        self.assertEqual(result.next_cap, third.retry_step_size)
        self.assertTrue(result.accepted_boundary_restored)
        replay_admission.assert_called_once()
        current_admission.assert_called_once()
        self.assertEqual(item.temporal_ledger, second.updated_ledger)
        self.assertEqual(item.time, 23 / 16)

    def test_third_pending_cursor_replays_latest_only_and_reaches_fourth_retry(self) -> None:
        (
            item,
            pending,
            base,
            active,
            replay_prepared,
            replayed,
            accepted_descriptor,
        ) = self._pending_replay_three_with_older_prefix_fixture()
        first, second, third = active

        latest_prior, checked_active = attempt._active_rejection_chain(
            pending,
            third.updated_ledger,
            physical_state_sha256=lifecycle_fixture.STATE,
        )
        self.assertEqual(latest_prior, second.updated_ledger)
        self.assertEqual(
            checked_active,
            tuple(
                p15._json_safe(outcome.evidence.as_mapping())
                for outcome in active
            ),
        )
        self.assertEqual(len(base.serialized_temporal_rejections), 1)
        self.assertEqual(
            third.updated_ledger.serialized_temporal_rejections[:1],
            base.serialized_temporal_rejections,
        )
        self.assertEqual(
            len(third.updated_ledger.serialized_temporal_rejections), 4
        )

        nested = pending.payload["retry_successor_payload_or_none"]
        fourth_plan = nested["successor_plan"]
        fourth = lifecycle_fixture.rejected(
            pending, third.updated_ledger, fourth_plan
        )
        self.assertIsInstance(fourth, p15.TDG6TemporalRetryRequired)
        current_prepared = SimpleNamespace(plan=fourth_plan)
        replayed_records: list[dict[str, object]] = []

        def reject_third(_prepared, **kwargs):
            mapping = p15._json_safe(third.evidence.as_mapping())
            replayed_records.append(mapping)
            kwargs["durable_rejection_sink"](mapping)
            raise third

        def reject_fourth(_prepared, **kwargs):
            kwargs["durable_rejection_sink"](fourth.evidence.as_mapping())
            raise fourth

        with patch.object(
            attempt,
            "array_content_sha256",
            return_value=lifecycle_fixture.STATE,
        ), patch.object(attempt, "_prepare_initial") as fresh_prepare, patch.object(
            attempt.tdg8_replay,
            "prepare_tdg8_persisted_retry_replay",
            return_value=replay_prepared,
        ) as replay_prepare, patch.object(
            attempt.tdg8_replay,
            "require_tdg8_persisted_retry_replay_admission",
            side_effect=reject_third,
        ) as replay_admission, patch.object(
            attempt.tdg8_replay,
            "replan_tdg8_after_persisted_retry_rejection",
            return_value=replayed,
        ), patch.object(
            attempt.tdg7,
            "prepare_tdg7_retry_successor",
            return_value=current_prepared,
        ), patch.object(
            attempt.tdg7,
            "require_tdg7_stage_safe_admission",
            side_effect=reject_fourth,
        ) as current_admission:
            result = attempt.attempt_once(
                item,
                pending,
                3 / 2,
                descriptor_sha256=accepted_descriptor,
                requested_cap=float.fromhex(nested["half_cap"]),
            )

        self.assertEqual(result.disposition, "temporal_retry_required")
        self.assertIs(result.retry_outcome, fourth)
        self.assertEqual(result.next_cap, fourth.retry_step_size)
        self.assertTrue(result.accepted_boundary_restored)
        fresh_prepare.assert_not_called()
        replay_admission.assert_called_once()
        current_admission.assert_called_once()
        self.assertEqual(
            replay_prepare.call_args.kwargs["expected_prior_retry_count"], 2
        )
        self.assertEqual(
            replay_prepare.call_args.kwargs["temporal_ledger"],
            second.updated_ledger,
        )
        self.assertEqual(
            replayed_records,
            [p15._json_safe(third.evidence.as_mapping())],
        )
        self.assertEqual(item.temporal_ledger, third.updated_ledger)
        self.assertEqual(item.time, float.fromhex(
            pending.payload["accepted_boundary_time"]["binary64_hex"]
        ))

    def test_second_pending_replay_rejects_malformed_earlier_prefix(self) -> None:
        for scenario in ("noncontiguous_count", "broken_plan_adjacency"):
            with self.subTest(scenario=scenario):
                item, pending, _first, _second, _prepared, _replayed = (
                    self._pending_replay_two_fixture()
                )
                records = list(item.temporal_ledger.serialized_temporal_rejections)
                earlier = json.loads(records[0])
                if scenario == "noncontiguous_count":
                    earlier["retry_count_for_current_macro_step"] = 7
                else:
                    earlier["attempted_macro_step_size"] *= 0.5
                    earlier["retry_step_size"] *= 0.5
                records[0] = p15._canonical(earlier).decode("utf-8")
                attacked_ledger = replace(
                    item.temporal_ledger,
                    serialized_temporal_rejections=tuple(records),
                )
                payload = dict(pending.payload)
                nested = dict(payload["retry_successor_payload_or_none"])
                nested["complete_rejection_prefix"] = records
                payload["retry_successor_payload_or_none"] = nested
                payload["TDG6_ledger_sha256"] = p15._hash(
                    p15._ledger_mapping(attacked_ledger)
                )
                payload.pop("cursor_chain_sha256", None)
                attacked_cursor = p15._cursor(payload)
                attacked_member = member(attacked_ledger)
                attacked_member.state = SimpleNamespace(
                    u=object(), p=object(), q=object()
                )
                with patch.object(
                    attempt,
                    "array_content_sha256",
                    return_value=lifecycle_fixture.STATE,
                ), patch.object(attempt, "_prepare_initial") as fresh, patch.object(
                    attempt.tdg8_replay, "prepare_tdg8_persisted_retry_replay"
                ) as replay_prepare:
                    with self.assertRaises(attempt.HLT16AttemptError):
                        attempt._replay_pending_rejection(
                            attacked_member, attacked_cursor, 1.5
                        )
                fresh.assert_not_called()
                replay_prepare.assert_not_called()

    def test_active_retry_prefix_is_bounded_by_the_frozen_cap(self) -> None:
        _item, pending, retry, _prepared, _replayed = self._pending_replay_fixture()
        serialized = p15._canonical(
            p15._json_safe(retry.evidence.as_mapping())
        ).decode("utf-8")
        count = p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP + 1
        over_limit = replace(
            retry.updated_ledger,
            cumulative_temporal_retry_count=count,
            current_macro_step_temporal_retry_count=count,
            serialized_temporal_rejections=(serialized,) * count,
        )
        nested = dict(pending.payload["retry_successor_payload_or_none"])
        nested["complete_rejection_prefix"] = [serialized] * count
        nested["retry_count"] = count
        synthetic = SimpleNamespace(
            payload={"retry_successor_payload_or_none": nested}
        )
        with self.assertRaisesRegex(
            attempt.HLT16AttemptError, "rejection prefix differs"
        ):
            attempt._active_rejection_chain(
                synthetic,  # type: ignore[arg-type]
                over_limit,
                physical_state_sha256=lifecycle_fixture.STATE,
            )

    def test_pending_replay_rejects_semantic_order_cardinality_and_scalar_type_drift(self) -> None:
        scenarios = ("scalar", "scalar_type", "channel_order", "failed_channels", "zero", "two")
        for scenario in scenarios:
            with self.subTest(scenario=scenario):
                item, pending, retry, prepared, replayed = self._pending_replay_fixture()

                def reject(_prepared, **kwargs):
                    sink = kwargs["durable_rejection_sink"]
                    mapping = copy.deepcopy(retry.evidence.as_mapping())
                    if scenario == "scalar":
                        mapping["attempted_macro_step_size"] = float(mapping["attempted_macro_step_size"]) + 1e-12
                    elif scenario == "scalar_type":
                        mapping["previous_step_index"] = float(mapping["previous_step_index"])
                    elif scenario == "channel_order":
                        values = list(mapping["channel_admissions"])
                        values[0], values[1] = values[1], values[0]
                        mapping["channel_admissions"] = tuple(values)
                    elif scenario == "failed_channels":
                        mapping["failed_channels"] = ("u:alpha",)
                    if scenario != "zero":
                        sink(mapping)
                    if scenario == "two":
                        sink(mapping)
                    raise retry

                with patch.object(attempt, "array_content_sha256", return_value=lifecycle_fixture.STATE), \
                     patch.object(attempt, "_prepare_initial", return_value=prepared), \
                     patch.object(attempt.tdg7, "require_tdg7_stage_safe_admission", side_effect=reject), \
                     patch.object(attempt.tdg7, "replan_tdg7_after_tdg6_retry", return_value=replayed):
                    with self.assertRaisesRegex(
                        attempt.HLT16AttemptError,
                        "pending TDG6 replay evidence differs",
                    ):
                        attempt._replay_pending_rejection(item, pending, 1.5)

    def test_source_cfl_scientific_and_unknown_paths_are_separated(self) -> None:
        scenarios = (
            (FakePathStop(source_retry={"retryable_source_only": True}), "source_retry_required", {}, attempt.RetryCounters()),
            (FakePathStop(cause=FakeCFL()), "cfl_retry_required", {"CFLRetryRequired": FakeCFL}, attempt.RetryCounters()),
            (FakeScientific(), "scientific_stop", {"GR0RuntimeStop": FakeScientific}, attempt.RetryCounters()),
            (RuntimeError("boom"), "invalid", {}, attempt.RetryCounters()),
        )
        for raised, expected, replacements, counters in scenarios:
            with self.subTest(expected=expected):
                item = member(); old = cursor()
                patches = [patch.object(attempt, "_prepare_initial", side_effect=raised),
                           patch.object(attempt.tdg6, "TDG6RefinementPathStop", FakePathStop)]
                if "CFLRetryRequired" in replacements:
                    patches.append(patch.object(attempt, "CFLRetryRequired", replacements["CFLRetryRequired"]))
                if "GR0RuntimeStop" in replacements:
                    patches.append(patch.object(attempt, "GR0RuntimeStop", replacements["GR0RuntimeStop"]))
                with patches[0], patches[1]:
                    if len(patches) == 3:
                        with patches[2]:
                            result = attempt.attempt_once(item, old, 1.5, descriptor_sha256=DESCRIPTOR, requested_cap=1 / 16, retry_counters=counters)
                    else:
                        result = attempt.attempt_once(item, old, 1.5, descriptor_sha256=DESCRIPTOR, requested_cap=1 / 16, retry_counters=counters)
                self.assertEqual(result.disposition, expected)
                self.assertTrue(result.accepted_boundary_restored)

    def test_source_and_cfl_caps_exhaust_without_durable_counter_mutation(self) -> None:
        source = FakePathStop(source_retry={"retryable_source_only": True})
        cfl = FakePathStop(cause=FakeCFL())
        for raised, expected, counter, replacement in (
            (source, "source_retry_exhausted", attempt.RetryCounters(source_current=32, source_total=32), None),
            (cfl, "cfl_retry_exhausted", attempt.RetryCounters(cfl_current=32, cfl_total=32), FakeCFL),
        ):
            with self.subTest(expected=expected):
                item = member(); old = cursor()
                with patch.object(attempt, "_prepare_initial", side_effect=raised), \
                     patch.object(attempt.tdg6, "TDG6RefinementPathStop", FakePathStop):
                    if replacement is None:
                        result = attempt.attempt_once(item, old, 1.5, descriptor_sha256=DESCRIPTOR, requested_cap=1 / 16, retry_counters=counter)
                    else:
                        with patch.object(attempt, "CFLRetryRequired", replacement):
                            result = attempt.attempt_once(item, old, 1.5, descriptor_sha256=DESCRIPTOR, requested_cap=1 / 16, retry_counters=counter)
                self.assertEqual(result.disposition, expected)
                self.assertEqual(counter.source_current + counter.cfl_current, 32)

    def test_preflight_and_pending_replay_fail_closed(self) -> None:
        item = member(); old = cursor()
        # A wrong descriptor must return invalid rather than invoke a solver.
        result = attempt.attempt_once(item, old, 1.5, descriptor_sha256="f" * 64, requested_cap=1 / 16)
        self.assertEqual(result.disposition, "invalid")
        self.assertEqual(item.state, "old")

        # A refinement stop raised before a snapshot exists remains durable-JSON
        # evidence; a live exception object must never leak into the record.
        item = member(); old = cursor()
        item.snapshot = Mock(side_effect=FakePathStop(cause=RuntimeError("preflight")))
        with patch.object(attempt.tdg6, "TDG6RefinementPathStop", FakePathStop):
            result = attempt.attempt_once(
                item, old, 1.5, descriptor_sha256=DESCRIPTOR,
                requested_cap=1 / 16,
            )
        self.assertEqual(result.disposition, "invalid")
        self.assertEqual(result.evidence["kind"], "preflight_refinement_stop")

    def test_retry_counter_contract_is_frozen(self) -> None:
        with self.assertRaises(ValueError):
            attempt.RetryCounters(source_cap=31)
        with self.assertRaises(ValueError):
            attempt.RetryCounters(source_current=34, source_total=34)

    def test_requested_cap_is_explicit_and_bounded(self) -> None:
        item = member(); old = cursor()
        for cap in (0.0, float("nan"), 1 / 8):
            with self.subTest(cap=cap):
                result = attempt.attempt_once(
                    item, old, 1.5, descriptor_sha256=DESCRIPTOR,
                    requested_cap=cap,
                )
                self.assertEqual(result.disposition, "invalid")


if __name__ == "__main__":
    unittest.main()
