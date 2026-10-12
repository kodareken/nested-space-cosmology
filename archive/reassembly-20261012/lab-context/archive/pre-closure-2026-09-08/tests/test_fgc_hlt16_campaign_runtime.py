"""Coordinator dispatch plus sealed-GEN0 persistence integration tests."""
from __future__ import annotations

import copy
from dataclasses import dataclass, replace
import json
import os
from pathlib import Path
import shutil
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as runtime
from recursive_horizons.fgc.evolution import hlt16_campaign_recovery as recovery
from recursive_horizons.fgc.evolution import proto15_runtime as p15
from recursive_horizons.fgc.evolution.boundary_domain import CausalBudgetState
from recursive_horizons.fgc.evolution.hlt16_campaign_schema import (
    CampaignCheckpoint,
    validate_checkpoint_successor,
)
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (
    HLT16CampaignStore,
    _walk_tree,
)
from recursive_horizons.fgc.evolution.hlt16_campaign_recovery import reconcile
from recursive_horizons.fgc.evolution.hlt16_member_codec import HLT16MemberSnapshot
from recursive_horizons.fgc.evolution.hlt16_progression_attempt import MemberAttemptOutcome, RetryCounters
from recursive_horizons.fgc.evolution.hlt16_state_store import (
    HLT16StateStore,
    HLT16StateStoreError,
    canonical,
)
from recursive_horizons.fgc.evolution.numerical_engine import array_content_sha256
from recursive_horizons.fgc.evolution.proto5_runtime import GR0RuntimeMonitorState
from recursive_horizons.fgc.evolution.proto19_progression_contract import construct_first_event
from recursive_horizons.fgc.evolution import proto19_progression_inputs as progression_inputs
from recursive_horizons.fgc.evolution.proto19_progression_inputs import (
    build_static_gr0_shells, reconstruct_gr0_members,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import TDG6TemporalLedger
from tests.fgc_gen0_fixture import (
    copy_sealed_gen0_store,
    reconstruct_sealed_gen0_members,
)
from tests import test_hlt16_lifecycle as lifecycle_fixture


KEY = "RK4-2049"
PHYSICAL_STATE_SHA256 = "9" * 64


@dataclass(frozen=True)
class _State:
    descriptor_sha256: str = "a" * 64
    cursor: object = None
    ledger: object = None
    source_current: int = 0
    source_total: int = 0
    cfl_current: int = 0
    cfl_total: int = 0
    pending_owner: str | None = None
    pending_cap_hex: str | None = None


class _Cursor:
    def __init__(self) -> None:
        self.sha256 = "b" * 64
        self.mode = "FRESH_READY"
        self.payload = {
            "protocol_artifact_id": "FGC-2-SF1-PROTO18",
            "accepted_state_sha256": "a" * 64,
            "accepted_boundary_time": {"binary64_hex": (23 / 16).hex()},
            "TDG6_ledger_sha256": "f" * 64,
        }


class _Store:
    def __init__(self) -> None:
        self.order: list[str] = []
        self.states = SimpleNamespace(
            persist_snapshot=lambda snapshot: SimpleNamespace(descriptor_sha256="c" * 64)
        )

    def active_writer_capability(self, _checkpoint):
        return object()

    def persist_snapshot(self, snapshot, **_kwargs):
        return self.states.persist_snapshot(snapshot)

    def publish_journal(self, record, **_kwargs):
        self.order.append("journal")

    def publish_checkpoint(self, current, *, previous, records, **_kwargs):
        self.order.append("checkpoint")

    def publish_terminal_lock(self, checkpoint, **_kwargs):
        self.order.append("terminal")


def _dispatch_member(*, time: float = 23 / 16) -> SimpleNamespace:
    """Minimal finite member shape required by the physical-ID bridge.

    The dispatch tests deliberately patch the hasher to a fixed digest, but
    still provide genuine finite ``u,p,q`` arrays so an accidental removal of
    the restored-state access cannot hide behind an opaque test double.
    """
    state = SimpleNamespace(
        u=np.array([[1.0]], dtype=np.float64),
        p=np.array([[2.0]], dtype=np.float64),
        q=np.array([[3.0]], dtype=np.float64),
    )
    return SimpleNamespace(time=time, state=state)


def _plan():
    return construct_first_event(Path("."), authorization_commit="0" * 40)


def _outcome(disposition: str) -> MemberAttemptOutcome:
    retry_plan = {"macro_width_hex": (1.0).hex()}
    return MemberAttemptOutcome(
        disposition=disposition,
        executed_plan=(
            {"synthetic": "plan"}
            if disposition in {"accepted_fine", "temporal_retry_required"}
            else retry_plan
            if disposition.startswith(("source", "cfl"))
            else None
        ),
        accepted_time=1.5 if disposition == "accepted_fine" else None,
        accepted_state=object() if disposition == "accepted_fine" else None,
        accepted_step_index=1 if disposition == "accepted_fine" else None,
        accepted_transaction_serial=1 if disposition == "accepted_fine" else None,
        accepted_ledger=object() if disposition == "accepted_fine" else None,
        retry_outcome=None,
        retry_counters=RetryCounters(1 if disposition.startswith("source") else 0,
                                     1 if disposition.startswith("source") else 0,
                                     1 if disposition.startswith("cfl") else 0,
                                     1 if disposition.startswith("cfl") else 0),
        next_cap=0.5 if disposition.endswith("required") else None,
        evidence={"kind": disposition},
        accepted_boundary_restored=disposition != "accepted_fine",
    )


class HLT16CampaignRuntimeDispatchTests(unittest.TestCase):
    """The coordinator's branching and publication order are independently pinned."""

    def _checkpoint(self, plan):
        state = _State(cursor={"opaque": True}, ledger={"opaque": True})
        return SimpleNamespace(plan_sha256=plan.sha256, terminal=None, members={KEY: state},
                               target={"binary64_hex": (1.5).hex()}, campaign_id="campaign",
                               journal_tip_sha256="d" * 64, journal_sequence=0,
                               event=23, generation=1)

    def _patches(self, store, disposition):
        cursor = _Cursor()
        result_checkpoint = SimpleNamespace(marker=disposition)
        record = {"record_sha256": "e" * 64}
        terminal_disposition = (
            "scientific_terminal" if disposition == "scientific_stop"
            else "invalid_terminal"
        )
        durable = runtime.DurableAttempt(terminal_disposition, result_checkpoint, ("e" * 64,))
        patches = [patch.multiple(
            runtime,
            _cursor=lambda value: cursor,
            restore_member_with_overlay=lambda *args, **kwargs: args[2],
            choose_requested_cap=lambda *args, **kwargs: 1.0,
            _record=lambda *args, **kwargs: dict(record),
            _record_after=lambda *args, **kwargs: dict(record),
            _successor=lambda *args, **kwargs: result_checkpoint,
            _publish=lambda store, previous, current, records, **_kwargs: (store.order.append("journal"), store.order.append("checkpoint")),
            _terminal_records=lambda *args, **kwargs: (durable, (dict(record),)),
            _retry_state=lambda state, outcome, owner: replace(state, pending_owner=owner, pending_cap_hex="0x1.0000000000000p-1"),
            MemberCampaignState=lambda *args, **kwargs: _State(),
            encode_member=lambda *args, **kwargs: object(),
            _ledger=lambda value: object(),
        ), patch.object(runtime, "array_content_sha256", return_value=PHYSICAL_STATE_SHA256), \
            patch.object(runtime.lifecycle, "close_fine", return_value=_Cursor()),
            patch.object(runtime.p15, "_ledger_mapping", return_value={})]
        class _Many:
            def __enter__(self):
                for item in patches: item.__enter__()
            def __exit__(self, *args):
                for item in reversed(patches): item.__exit__(*args)
        return _Many()

    def test_typed_progression_time_preserves_its_exact_rational_identity(self):
        plan = _plan()
        self.assertEqual(
            runtime._time(plan.target_time),
            {"rational": "3/2", "binary64_hex": (3 / 2).hex()},
        )

    def test_all_non_temporal_dispositions_publish_one_complete_suffix(self):
        plan = _plan()
        for disposition in ("accepted_fine", "source_retry_required", "cfl_retry_required",
                            "source_retry_exhausted", "cfl_retry_exhausted", "scientific_stop", "invalid"):
            with self.subTest(disposition=disposition):
                store = _Store(); checkpoint = self._checkpoint(plan)
                with self._patches(store, disposition):
                    result = runtime.execute_one_attempt(
                        plan, store, checkpoint, _dispatch_member(), key=KEY,
                        runner=lambda *args, **kwargs: _outcome(disposition),
                    )
                expected = (
                    "scientific_terminal" if disposition == "scientific_stop"
                    else "invalid_terminal" if disposition.endswith("exhausted") or disposition == "invalid"
                    else disposition
                )
                self.assertEqual(result.disposition, expected)
                self.assertEqual(store.order[:2], ["journal", "checkpoint"])
                if disposition.endswith("exhausted") or disposition in {"scientific_stop", "invalid"}:
                    self.assertEqual(store.order[-1], "terminal")

    def test_unknown_disposition_is_not_published(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        # A fake result bypasses the immutable outcome constructor specifically
        # to attack the coordinator's final ownership switch.
        bogus = SimpleNamespace(disposition="unowned", retry_counters=RetryCounters(), evidence={}, next_cap=None)
        with self._patches(store, "unowned"):
            with self.assertRaises(runtime.HLT16CampaignRuntimeError):
                runtime.execute_one_attempt(plan, store, checkpoint, _dispatch_member(), key=KEY,
                                            runner=lambda *args, **kwargs: bogus)
        self.assertEqual(store.order, [])

    def test_temporal_rejection_requires_sink_then_publishes_two_records(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        outcome = _outcome("temporal_retry_required")
        outcome = replace(outcome, retry_outcome=SimpleNamespace(updated_ledger=object()))
        with self._patches(store, "temporal_retry_required"), \
             patch.object(runtime.lifecycle, "close_retry", return_value=_Cursor()):
            result = runtime.execute_one_attempt(
                plan, store, checkpoint, _dispatch_member(), key=KEY,
                runner=lambda *args, **kwargs: (
                    kwargs["temporal_rejection_sink"](
                        {"kind": "rejection", "initial_state_sha256": PHYSICAL_STATE_SHA256}
                    ) or outcome
                ),
            )
        self.assertEqual(result.disposition, "temporal_retry_required")
        self.assertEqual(len(result.record_sha256s), 2)
        self.assertEqual(store.order, ["journal", "journal", "checkpoint"])

    def test_temporal_retry_successor_preserves_frozen_minimum_width(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        cursor, ledger = lifecycle_fixture.proto18()
        predecessor = lifecycle_fixture.plan(cursor)
        retry = lifecycle_fixture.rejected(cursor, ledger, predecessor)
        self.assertIsInstance(retry, p15.TDG6TemporalRetryRequired)
        checkpoint.members = {
            KEY: replace(
                checkpoint.members[KEY],
                descriptor_sha256=lifecycle_fixture.DESCRIPTOR,
            )
        }
        outcome = replace(
            _outcome("temporal_retry_required"),
            executed_plan=predecessor,
            retry_outcome=retry,
            next_cap=retry.retry_step_size,
        )
        ledger_mapping = p15._ledger_mapping
        close_retry = runtime.lifecycle.close_retry
        with self._patches(store, "temporal_retry_required"), \
             patch.object(runtime, "_cursor", return_value=cursor), \
             patch.object(runtime, "_ledger", return_value=ledger), \
             patch.object(runtime, "array_content_sha256", return_value=lifecycle_fixture.STATE), \
             patch.object(runtime.p15, "_ledger_mapping", side_effect=ledger_mapping), \
             patch.object(runtime.lifecycle, "close_retry", wraps=close_retry) as close_retry_spy:
            result = runtime.execute_one_attempt(
                plan, store, checkpoint, _dispatch_member(), key=KEY,
                runner=lambda *args, **kwargs: (
                    kwargs["temporal_rejection_sink"]({
                        "kind": "rejection",
                        "initial_state_sha256": lifecycle_fixture.STATE,
                    }) or outcome
                ),
            )

        successor = close_retry_spy.call_args.args[4]
        successor_mapping = p15._plan_mapping(successor)
        expected = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            23 / 16,
            3 / 2,
            retry.retry_step_size,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ))
        self.assertEqual(result.disposition, "temporal_retry_required")
        self.assertEqual(successor_mapping, expected)
        self.assertEqual(
            successor_mapping["minimum_width_hex_or_none"],
            "0x1.0000000000000p-30",
        )

        omitted = p15.plan_forward_proto14_subdivision(
            23 / 16, 3 / 2, retry.retry_step_size
        )
        self.assertIsNone(p15._plan_mapping(omitted)["minimum_width_hex_or_none"])
        with self.assertRaisesRegex(
            runtime.lifecycle.HLT16LifecycleError,
            "TDG6 retry successor is not the exact TDG7 half-step plan",
        ):
            close_retry(
                cursor,
                ledger,
                retry,
                predecessor,
                omitted,
                "e" * 64,
                "d" * 64,
                evolution_state_sha256=lifecycle_fixture.STATE,
            )

    def test_cap_selection_rejects_no_remainder_before_runner(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        called = False
        with patch.object(runtime, "_cursor", return_value=_Cursor()), \
             patch.object(runtime, "restore_member_with_overlay", return_value=_dispatch_member(time=1.5)), \
             patch.object(runtime, "array_content_sha256", return_value=PHYSICAL_STATE_SHA256), \
             patch.object(runtime, "choose_requested_cap", side_effect=runtime.HLT16CampaignRuntimeError("no remainder")):
            with self.assertRaises(runtime.HLT16CampaignRuntimeError):
                runtime.execute_one_attempt(plan, store, checkpoint, _dispatch_member(time=1.5), key=KEY,
                                            runner=lambda *args, **kwargs: self.fail("runner must not be called"))
        self.assertFalse(called)

    def test_source_retry_payload_preserves_every_frozen_provenance_field(self):
        outcome = _outcome("source_retry_required")
        payload = runtime._source_or_cfl_rejection_payload(
            owner="source", key=KEY, state=_State(), cursor=_Cursor(), outcome=outcome,
            exhausted=False,
        )
        self.assertEqual(payload, {
            "member_key": KEY,
            "evidence": {"kind": "source_retry_required"},
            "predecessor_descriptor_sha256": "a" * 64,
            "predecessor_cursor_sha256": "b" * 64,
            "retry_current": 1,
            "retry_total": 1,
            "retry_limit": 32,
            "attempted_cap_hex": (1.0).hex(),
            "next_cap_hex": (0.5).hex(),
            "exhausted": False,
            "no_commit": True,
        })

    def test_exhausted_source_retry_cannot_carry_a_successor_cap(self):
        outcome = _outcome("source_retry_exhausted")
        outcome = replace(outcome, next_cap=0.5)
        with self.assertRaisesRegex(runtime.HLT16CampaignRuntimeError, "retained a successor cap"):
            runtime._source_or_cfl_rejection_payload(
                owner="source", key=KEY, state=_State(), cursor=_Cursor(), outcome=outcome,
                exhausted=True,
            )

    def test_cfl_rejection_serializes_actual_macro_width_not_outer_requested_cap(self):
        # TDG7 may choose a lattice-safe macro width smaller than the outer
        # requested cap.  The durable retry relation is defined by the step
        # actually attempted, so the next cap is its exact binary half.
        macro = 0.75
        outcome = replace(
            _outcome("cfl_retry_required"),
            executed_plan={"macro_width_hex": macro.hex()},
            next_cap=macro * 0.5,
        )
        payload = runtime._source_or_cfl_rejection_payload(
            owner="cfl", key=KEY, state=_State(), cursor=_Cursor(), outcome=outcome,
            exhausted=False,
        )
        self.assertEqual(payload["attempted_cap_hex"], macro.hex())
        self.assertEqual(payload["next_cap_hex"], (macro * 0.5).hex())

    def test_source_cfl_rejection_refuses_absent_or_malformed_macro_width(self):
        base = _outcome("cfl_retry_required")
        for plan in (None, {}, {"macro_width_hex": "not-a-binary64"}):
            with self.subTest(plan=plan):
                outcome = replace(base, executed_plan=plan)
                with self.assertRaisesRegex(runtime.HLT16CampaignRuntimeError, "executed plan"):
                    runtime._source_or_cfl_rejection_payload(
                        owner="cfl", key=KEY, state=_State(), cursor=_Cursor(), outcome=outcome,
                        exhausted=False,
                    )

    def test_durable_tdg6_rejection_cannot_be_relabelled_as_an_accepted_step(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        accepted = _outcome("accepted_fine")
        with self._patches(store, "accepted_fine"):
            with self.assertRaisesRegex(runtime.HLT16CampaignRuntimeError, "contradicts its durable TDG6 rejection"):
                runtime.execute_one_attempt(
                    plan, store, checkpoint, _dispatch_member(), key=KEY,
                    runner=lambda *args, **kwargs: (
                        kwargs["temporal_rejection_sink"](
                            {"kind": "rejection", "initial_state_sha256": PHYSICAL_STATE_SHA256}
                        ) or accepted
                    ),
                )
        # The append-only rejection is intentionally preserved; there is no
        # contradictory checkpoint successor to make the stale proposal look
        # accepted.
        self.assertEqual(store.order, ["journal"])

    def test_accepted_snapshot_is_owned_by_the_successor_generation(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        seen: dict[str, object] = {}
        with self._patches(store, "accepted_fine"), \
             patch.object(runtime, "encode_member", side_effect=lambda *args, **kwargs: (
                 seen.update(kwargs) or object()
             )):
            runtime.execute_one_attempt(
                plan, store, checkpoint, _dispatch_member(), key=KEY,
                runner=lambda *args, **kwargs: _outcome("accepted_fine"),
            )
        self.assertEqual(seen["generation"], checkpoint.generation + 1)

    def test_orphan_payload_replay_requires_the_exact_recomputed_semantic_state(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        expected = "7" * 64
        with self._patches(store, "accepted_fine"), \
             patch.object(runtime, "encode_member", return_value=SimpleNamespace(arrays={})), \
             patch.object(runtime, "semantic_state_sha256", return_value=expected):
            result = runtime.execute_one_attempt(
                plan, store, checkpoint, _dispatch_member(), key=KEY,
                runner=lambda *args, **kwargs: _outcome("accepted_fine"),
                expected_orphan_semantic_sha256=expected,
            )
        self.assertEqual(result.disposition, "accepted_fine")

    def test_orphan_payload_replay_rejects_a_different_state_before_persistence(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        expected = "7" * 64
        persisted = False

        def forbidden_persist(_snapshot):
            nonlocal persisted
            persisted = True
            raise AssertionError("mismatched replay must not publish another state")

        store.persist_snapshot = forbidden_persist
        with self._patches(store, "accepted_fine"), \
             patch.object(runtime, "encode_member", return_value=SimpleNamespace(arrays={})), \
             patch.object(runtime, "semantic_state_sha256", return_value="8" * 64):
            with self.assertRaisesRegex(
                runtime.HLT16CampaignRuntimeError,
                "differs from the authenticated orphan payload",
            ):
                runtime.execute_one_attempt(
                    plan, store, checkpoint, _dispatch_member(), key=KEY,
                    runner=lambda *args, **kwargs: _outcome("accepted_fine"),
                    expected_orphan_semantic_sha256=expected,
                )
        self.assertFalse(persisted)
        self.assertEqual(store.order, [])

    def test_orphan_payload_replay_cannot_be_relabelled_as_a_retry(self):
        plan = _plan(); store = _Store(); checkpoint = self._checkpoint(plan)
        with self._patches(store, "source_retry_required"):
            with self.assertRaisesRegex(
                runtime.HLT16CampaignRuntimeError,
                "did not reproduce an accepted fine state",
            ):
                runtime.execute_one_attempt(
                    plan, store, checkpoint, _dispatch_member(), key=KEY,
                    runner=lambda *args, **kwargs: _outcome("source_retry_required"),
                    expected_orphan_semantic_sha256="7" * 64,
                )
        self.assertEqual(store.order, [])

    def test_tdg6_sink_accepts_distinct_physical_identity_and_rejects_mutation(self):
        """The descriptor address and TDG6's u,p,q address are both binding."""
        plan = _plan()
        outcome = replace(
            _outcome("temporal_retry_required"),
            retry_outcome=SimpleNamespace(updated_ledger=object()),
        )
        store = _Store()
        checkpoint = self._checkpoint(plan)
        with self._patches(store, "temporal_retry_required"), \
             patch.object(runtime.lifecycle, "close_retry", return_value=_Cursor()):
            result = runtime.execute_one_attempt(
                plan,
                store,
                checkpoint,
                _dispatch_member(),
                key=KEY,
                runner=lambda *args, **kwargs: (
                    kwargs["temporal_rejection_sink"](
                        {"initial_state_sha256": PHYSICAL_STATE_SHA256}
                    ) or outcome
                ),
            )
        self.assertEqual(result.disposition, "temporal_retry_required")
        self.assertNotEqual(PHYSICAL_STATE_SHA256, _Cursor().payload["accepted_state_sha256"])

        store = _Store()
        checkpoint = self._checkpoint(plan)
        with self._patches(store, "temporal_retry_required"):
            with self.assertRaisesRegex(
                runtime.HLT16CampaignRuntimeError,
                "restored evolution state",
            ):
                runtime.execute_one_attempt(
                    plan,
                    store,
                    checkpoint,
                    _dispatch_member(),
                    key=KEY,
                    runner=lambda *args, **kwargs: (
                        kwargs["temporal_rejection_sink"](
                            {"initial_state_sha256": "a" * 64}
                        ) or outcome
                    ),
                )
        self.assertEqual(store.order, [])

    def test_common_event_is_not_reused_for_a_later_generation(self):
        checkpoint = self._checkpoint(_plan())
        checkpoint.event = 24
        with self.assertRaisesRegex(runtime.HLT16CampaignRuntimeError, "only the first common event"):
            runtime.commit_common_event(checkpoint, _Store())

    def _common_event_checkpoint(self, *, pending_key: str | None = None):
        plan = _plan()
        cursors: dict[int, _Cursor] = {}
        members = {}
        for key in runtime.MEMBER_KEYS:
            state = _State(
                cursor={"member": key}, ledger={"member": key},
                pending_owner="source" if key == pending_key else None,
                pending_cap_hex=(0.5).hex() if key == pending_key else None,
            )
            cursor = _Cursor()
            cursor.payload.update({
                "campaign_id": "campaign",
                "member_key": key,
                "committed_common_event_index": 23,
                "event_target_time": {"binary64_hex": (1.5).hex()},
                "accepted_boundary_time": {"rational": "3/2", "binary64_hex": (1.5).hex()},
                "cursor_generation": 7,
                "attempt_serial": 4,
                "attempt_id": f"campaign:{key}:23:4",
                "journal_tip_sha256": "d" * 64,
                "cursor_chain_parent_sha256": "c" * 64,
                "retry_successor_payload_or_none": None,
            })
            cursors[id(state.cursor)] = cursor
            members[key] = state
        checkpoint = SimpleNamespace(
            plan_sha256=plan.sha256, terminal=None, members=members,
            target={"rational": "3/2", "binary64_hex": (1.5).hex()}, campaign_id="campaign",
            journal_tip_sha256="d" * 64, journal_sequence=12,
            event=23, generation=9, sha256="e" * 64,
        )
        return checkpoint, cursors

    def test_common_event_rejects_any_pending_member_even_when_every_cursor_is_fresh(self):
        checkpoint, cursors = self._common_event_checkpoint(pending_key="SSPRK3-16385")
        with patch.object(runtime, "_cursor", side_effect=lambda value: cursors[id(value)]):
            with self.assertRaisesRegex(runtime.HLT16CampaignRuntimeError, "retain pending retry state"):
                runtime.commit_common_event(checkpoint, _Store())

    def test_common_event_resets_exact_proto15_attempt_identity_for_all_members(self):
        checkpoint, cursors = self._common_event_checkpoint()
        captured: dict[str, object] = {}

        class _NewCursor:
            def __init__(self, payload):
                self.payload = dict(payload)
                self.sha256 = "f" * 64

            def validate(self):
                return None

        def successor(previous, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace(marker="event")

        record = {"record_sha256": "a" * 64, "sequence": 13}
        with patch.object(runtime, "_cursor", side_effect=lambda value: cursors[id(value)]), \
             patch.object(runtime.p15, "_cursor", _NewCursor), \
             patch.object(runtime, "_record", return_value=record), \
             patch.object(runtime, "_successor", side_effect=successor), \
             patch.object(runtime, "_publish"):
            result = runtime.commit_common_event(checkpoint, _Store())
        self.assertEqual(result.disposition, "event_complete")
        self.assertEqual(captured["event"], 24)
        self.assertEqual(captured["target"], {"rational": "25/16", "binary64_hex": (25 / 16).hex()})
        for key, state in captured["members"].items():
            payload = state.cursor
            self.assertEqual(payload["attempt_serial"], 0)
            self.assertEqual(payload["attempt_id"], f"campaign:{key}:24:0")
            self.assertEqual(payload["committed_common_event_index"], 24)
            self.assertIsNone(payload["retry_successor_payload_or_none"])
            self.assertIsNone(state.pending_owner)
            self.assertIsNone(state.pending_cap_hex)

class HLT16CampaignRuntimeIntegrationTests(unittest.TestCase):
    """Exercise real GEN0, codec, schema, lifecycle, runtime, and store layers.

    The only replacement is a typed finite outcome in place of a PDE proposal.
    The stubs deliberately mutate a real restored member into a metadata-
    consistent finite boundary, so codec and checkpoint validation still own
    all persistence decisions.
    """

    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]
        cls.plan = construct_first_event(cls.root, authorization_commit="2" * 40)
        cls.reconstructed = reconstruct_sealed_gen0_members(cls.root)
        cls.seed_temporary = TemporaryDirectory()
        cls.seed_root = Path(cls.seed_temporary.name) / "calibration"
        copy_sealed_gen0_store(cls.root, cls.seed_root)
        cls.seed_checkpoint = runtime.initialize_or_recover_generation_one(
            cls.plan, HLT16CampaignStore(cls.seed_root), cls.reconstructed
        )

    @classmethod
    def tearDownClass(cls):
        cls.seed_temporary.cleanup()

    def bridge(self):
        temporary = TemporaryDirectory()
        path = Path(temporary.name) / "calibration"
        shutil.copytree(self.seed_root, path)
        store = HLT16CampaignStore(path)
        token = store.acquire_active_writer(
            self.seed_checkpoint, host=os.uname().nodename, pid=os.getpid(),
        )
        store._test_writer_token = token
        return temporary, store, self.seed_checkpoint

    @staticmethod
    def _accepted_stub(member, cursor, target, *, requested_cap, retry_counters, **_kwargs):
        plan = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            member.time, target, requested_cap,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ))
        endpoint = float.fromhex(plan["boundaries_hex"][-1])
        prior = member.temporal_ledger
        ledger = TDG6TemporalLedger(
            last_accepted_time=endpoint,
            accepted_macro_step_count=prior.accepted_macro_step_count + 1,
            cumulative_temporal_retry_count=prior.cumulative_temporal_retry_count,
            current_macro_step_temporal_retry_count=0,
            last_accepted_macro_step_temporal_retry_count=prior.current_macro_step_temporal_retry_count,
            accumulated_debit_vector=prior.accumulated_debit_vector,
            serialized_temporal_rejections=prior.serialized_temporal_rejections,
        )
        member.time = endpoint
        member.step_index += 1
        member.transaction_serial += 1
        member.temporal_ledger = ledger
        monitor = member.transaction.state
        member.transaction.state = GR0RuntimeMonitorState(
            accepted_stage_count=monitor.accepted_stage_count + 1,
            last_transaction_serial=member.transaction_serial - 1,
            last_accepted_time=endpoint,
        )
        causal = member.transaction.causal_state
        member.transaction.causal_state = CausalBudgetState(
            accepted_time=endpoint,
            accumulated_characteristic_distance=causal.accumulated_characteristic_distance,
            previous_speed_upper=causal.previous_speed_upper,
        )
        return MemberAttemptOutcome(
            "accepted_fine", plan, endpoint, member.state, member.step_index,
            member.transaction_serial, ledger, None, retry_counters, None,
            {"kind": "finite-typed-stub"}, False,
        )

    @staticmethod
    def _retry_stub(disposition, *, exhausted):
        def runner(_member, _cursor, _target, *, retry_counters, **_kwargs):
            attempted = _kwargs["requested_cap"]
            source_next = retry_counters.source_current + 1
            cfl_next = retry_counters.cfl_current + 1
            counters = RetryCounters(
                source_current=source_next if disposition.startswith("source") else retry_counters.source_current,
                source_total=retry_counters.source_total + 1 if disposition.startswith("source") else retry_counters.source_total,
                cfl_current=cfl_next if disposition.startswith("cfl") else retry_counters.cfl_current,
                cfl_total=retry_counters.cfl_total + 1 if disposition.startswith("cfl") else retry_counters.cfl_total,
            )
            return MemberAttemptOutcome(
                disposition, {"macro_width_hex": attempted.hex()}, None, None,
                None, None, None, None,
                counters, None if exhausted else attempted * 0.5,
                {"kind": disposition}, True,
            )
        return runner

    @staticmethod
    def _terminal_stub(disposition):
        def runner(_member, _cursor, _target, *, retry_counters, **_kwargs):
            return MemberAttemptOutcome(
                disposition, None, None, None, None, None, None, None,
                retry_counters, None, {"kind": disposition}, True,
            )
        return runner

    @staticmethod
    def _temporal_retry_stub(
        member,
        cursor,
        target,
        *,
        requested_cap,
        retry_counters,
        temporal_rejection_sink,
        **_kwargs,
    ):
        """Emit native tuple-bearing TDG6 evidence without running the PDE."""
        prior = member.temporal_ledger
        if prior is None:
            raise AssertionError("synthetic retry member lacks its TDG6 ledger")
        if cursor.mode == "RETRY_PENDING":
            pending = cursor.payload["retry_successor_payload_or_none"]
            if not isinstance(pending, dict):
                raise AssertionError("synthetic retry cursor lacks its durable plan")
            executed = dict(pending["successor_plan"])
            if requested_cap.hex() != pending["half_cap"]:
                raise AssertionError("synthetic retry cap differs from durable cursor")
        else:
            executed = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
                member.time,
                target,
                requested_cap,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            ))
        admission = p15._synthetic_failed_admission(
            p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]]
        )
        attempted = float.fromhex(executed["macro_width_hex"])
        evidence = p15.TDG6TemporalRetryEvidence(
            event_type="rejected_TDG6_temporal_admission",
            method=p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]],
            initial_time=member.time,
            attempted_macro_step_size=attempted,
            retry_step_size=attempted * 0.5,
            initial_state_sha256=array_content_sha256(
                member.state.u, member.state.p, member.state.q
            ),
            previous_step_index=cursor.payload["previous_step_index"],
            previous_transaction_serial=cursor.payload[
                "previous_transaction_serial"
            ],
            retry_count_for_current_macro_step=(
                prior.current_macro_step_temporal_retry_count + 1
            ),
            cumulative_temporal_retry_count=(
                prior.cumulative_temporal_retry_count + 1
            ),
            failed_channels=admission.failed_channels,
            channel_admissions=admission.channel_admissions,
        )
        native = evidence.as_mapping()
        if not isinstance(native["failed_channels"], tuple):
            raise AssertionError("synthetic evidence lost its native tuple")
        temporal_rejection_sink(native)
        updated = p15._updated_retry_ledger(prior, native)
        retry = p15.TDG6TemporalRetryRequired(evidence, updated)
        return MemberAttemptOutcome(
            "temporal_retry_required",
            executed,
            None,
            None,
            None,
            None,
            None,
            retry,
            retry_counters,
            evidence.retry_step_size,
            {"kind": "synthetic_TDG6_temporal_retry"},
            True,
        )

    def test_real_bridge_is_idempotent_and_restartable(self):
        temporary, store, first = self.bridge()
        with temporary:
            store.release_active_writer(
                owner_token=store._test_writer_token,
                host=os.uname().nodename,
                pid=os.getpid(),
            )
            second = runtime.initialize_or_recover_generation_one(
                self.plan, store, self.reconstructed
            )
            self.assertEqual(first.sha256, second.sha256)
            store._test_writer_token = store.acquire_active_writer(
                second, host=os.uname().nodename, pid=os.getpid(),
            )
            self.assertEqual(store.recover().state, "live_local_writer")
            self.assertFalse(store.recover().safe_to_restart)

    def test_one_real_member_accepted_transition_round_trips_exactly(self):
        temporary, store, checkpoint = self.bridge()
        with temporary:
            key = "RK4-2049"
            member = copy.deepcopy(self.reconstructed.members[key])
            # A finite typed test transition spans one synthetic full event;
            # production cap selection remains separately frozen and tested.
            with patch.object(runtime, "choose_requested_cap", side_effect=lambda _c, _k, item: 1.5 - item.time):
                result = runtime.execute_one_attempt(
                    self.plan, store, checkpoint, member, key=key,
                    runner=self._accepted_stub,
                )
            restored = copy.deepcopy(self.reconstructed.members[key])
            runtime.restore_member_with_overlay(store, result.checkpoint, restored, key=key)
            self.assertEqual(restored.time.hex(), member.time.hex())
            self.assertEqual(restored.step_index, member.step_index)
            self.assertEqual(restored.transaction_serial, member.transaction_serial)
            for name in ("u", "p", "q"):
                self.assertEqual(getattr(restored.state, name).tobytes(), getattr(member.state, name).tobytes())
            self.assertEqual(store.recover().state, "live_local_writer")

    def test_native_tdg6_depth_three_retry_matches_immediate_and_recovered_checkpoint(self):
        """The durable boundary owns one representation before and after I/O."""
        temporary, store, checkpoint = self.bridge()
        with temporary:
            key = "RK4-2049"
            current = checkpoint
            final_records = None
            final_previous = None
            final_result = None
            for depth in range(1, 4):
                final_previous = current
                member = copy.deepcopy(self.reconstructed.members[key])
                with patch.object(
                    store,
                    "publish_checkpoint",
                    wraps=store.publish_checkpoint,
                ) as publish_checkpoint:
                    final_result = runtime.execute_one_attempt(
                        self.plan,
                        store,
                        current,
                        member,
                        key=key,
                        runner=self._temporal_retry_stub,
                    )
                self.assertEqual(final_result.disposition, "temporal_retry_required")
                final_records = publish_checkpoint.call_args.kwargs["records"]
                rejection = final_records[0]["payload"]["evidence"]
                self.assertIsInstance(rejection["failed_channels"], list)
                self.assertIsInstance(rejection["channel_admissions"], list)
                self.assertEqual(rejection["retry_count_for_current_macro_step"], depth)
                # This is the exact immediate validation that failed before
                # the durable record boundary normalized native tuples.
                validate_checkpoint_successor(
                    current,
                    final_result.checkpoint,
                    final_records,
                )
                current = final_result.checkpoint

            assert final_previous is not None
            assert final_records is not None
            assert final_result is not None
            final_state = final_result.checkpoint.members[key]
            final_cursor = p15.Proto15Cursor(dict(final_state.cursor))
            self.assertEqual(final_state.pending_owner, "temporal")
            self.assertEqual(
                final_state.ledger["current_macro_step_temporal_retry_count"],
                3,
            )
            self.assertEqual(
                final_cursor.payload["retry_successor_payload_or_none"]["retry_count"],
                3,
            )

            persisted = store.load_state(
                final_previous.members[key].descriptor_sha256
            )
            rebuilt = recovery._resolve_tdg6_rejection(
                final_previous,
                final_records[0],
                persisted.arrays,
            )
            self.assertFalse(rebuilt.physical_state_advanced)
            self.assertEqual(
                rebuilt.checkpoint.sha256,
                final_result.checkpoint.sha256,
            )
            self.assertEqual(
                tuple(item["record_sha256"] for item in rebuilt.records),
                final_result.record_sha256s,
            )
            rebuilt_state = rebuilt.checkpoint.members[key]
            self.assertEqual(rebuilt_state.cursor, final_state.cursor)
            self.assertEqual(rebuilt_state.ledger, final_state.ledger)
            self.assertEqual(
                p15.Proto15Cursor(dict(rebuilt_state.cursor)).sha256,
                final_cursor.sha256,
            )
            self.assertEqual(
                p15._hash(dict(rebuilt_state.ledger)),
                p15._hash(dict(final_state.ledger)),
            )

            for record in final_records:
                path = store.root / "journal" / (
                    f"{record['sequence']:020d}-{record['record_sha256']}.journal"
                )
                raw = path.read_bytes()
                decoded = json.loads(raw.decode("ascii"))
                self.assertEqual(raw, canonical(record))
                self.assertEqual(decoded, record)
                self.assertEqual(
                    decoded["record_sha256"],
                    next(
                        item["record_sha256"]
                        for item in rebuilt.records
                        if item["sequence"] == decoded["sequence"]
                    ),
                )

            checkpoint_path = store.root / "checkpoints" / (
                f"{final_result.checkpoint.generation:020d}-"
                f"{final_result.checkpoint.sha256}.json"
            )
            checkpoint_raw = checkpoint_path.read_bytes()
            self.assertEqual(
                checkpoint_raw,
                canonical(final_result.checkpoint.mapping()),
            )
            restored_checkpoint = CampaignCheckpoint.from_mapping(
                json.loads(checkpoint_raw.decode("ascii"))
            )
            self.assertEqual(
                restored_checkpoint.sha256,
                rebuilt.checkpoint.sha256,
            )
            self.assertEqual(
                store.authenticated_snapshot().checkpoint.sha256,
                rebuilt.checkpoint.sha256,
            )

    def test_public_raw_store_cannot_persist_an_evolved_campaign_snapshot(self):
        temporary, store, checkpoint = self.bridge()
        with temporary:
            key = "RK4-2049"
            member = copy.deepcopy(self.reconstructed.members[key])
            with patch.object(
                runtime,
                "choose_requested_cap",
                side_effect=lambda _c, _k, item: 1.5 - item.time,
            ):
                accepted = runtime.execute_one_attempt(
                    self.plan,
                    store,
                    checkpoint,
                    member,
                    key=key,
                    runner=self._accepted_stub,
                )
            address = accepted.checkpoint.members[key].descriptor_sha256
            persisted = store.load_state(address)
            snapshot = HLT16MemberSnapshot(
                persisted.arrays,
                persisted.descriptor["metadata"],
            )
            before = _walk_tree(store.root)
            raw = HLT16StateStore(store.root)
            with self.assertRaisesRegex(
                HLT16StateStoreError,
                "evolved persistence requires campaign writer authority",
            ):
                raw.persist_snapshot(snapshot)
            self.assertEqual(_walk_tree(store.root), before)

    def test_persisted_successor_rehydrates_exactly_without_historical_raw_reopen(self):
        """A generation-one restart owns arrays through HLT16, never PREF26."""
        temporary, store, checkpoint = self.bridge()
        with temporary:
            key = "RK4-2049"
            original = copy.deepcopy(self.reconstructed.members[key])
            with patch.object(runtime, "choose_requested_cap", side_effect=lambda _c, _k, item: 1.5 - item.time):
                accepted = runtime.execute_one_attempt(
                    self.plan, store, checkpoint, original, key=key,
                    runner=self._accepted_stub,
                )
            # The static shell constructor may build equations/operators, but
            # it must not reach the historical raw container.  The exact
            # evolved fields below come solely from the persisted descriptor.
            with patch.object(
                progression_inputs,
                "reimport_matches_pref26_evidence",
                side_effect=AssertionError("historical raw reopen"),
            ):
                static = build_static_gr0_shells(self.root)[key]
                runtime.restore_member_with_overlay(
                    store, accepted.checkpoint, static, key=key,
                )
            self.assertEqual(static.time.hex(), original.time.hex())
            self.assertEqual(static.step_index, original.step_index)
            self.assertEqual(static.transaction_serial, original.transaction_serial)
            self.assertEqual(static.source_retry_count, original.source_retry_count)
            self.assertEqual(static.CFL_retry_count, original.CFL_retry_count)
            for name in ("u", "p", "q"):
                self.assertEqual(
                    getattr(static.state, name).tobytes(),
                    getattr(original.state, name).tobytes(),
                )

    def test_active_writer_bound_to_generation_one_remains_valid_after_its_accepted_successor(self):
        temporary, store, checkpoint = self.bridge()
        with temporary:
            token = store._test_writer_token
            key = "RK4-2049"
            with patch.object(runtime, "choose_requested_cap", side_effect=lambda _c, _k, item: 1.5 - item.time):
                advanced = runtime.execute_one_attempt(
                    self.plan, store, checkpoint,
                    copy.deepcopy(self.reconstructed.members[key]),
                    key=key, runner=self._accepted_stub,
                )
            # The writer began on the exact authenticated predecessor.  The
            # new checkpoint is its direct, validated successor, so status
            # must retain the live local lease rather than treating ordinary
            # progression as an authority mismatch.
            status = store.recover()
            self.assertEqual(status.checkpoint_sha256, advanced.checkpoint.sha256)
            self.assertEqual(status.state, "live_local_writer")
            self.assertFalse(status.safe_to_restart)
            store.release_active_writer(
                owner_token=token, host=os.uname().nodename, pid=os.getpid(),
            )
            self.assertEqual(store.recover().state, "clean_checkpoint")

    def test_source_and_cfl_retry_then_exhaustion_close_invalid_checkpoints(self):
        for owner in ("source", "cfl"):
            with self.subTest(owner=owner):
                temporary, store, checkpoint = self.bridge()
                with temporary:
                    key = "RK4-2049"
                    required = None
                    for _attempt in range(32):
                        required = runtime.execute_one_attempt(
                            self.plan, store, checkpoint,
                            copy.deepcopy(self.reconstructed.members[key]), key=key,
                            runner=self._retry_stub(f"{owner}_retry_required", exhausted=False),
                        )
                        checkpoint = required.checkpoint
                    assert required is not None
                    self.assertEqual(required.disposition, f"{owner}_retry_required")
                    self.assertEqual(store.recover().state, "live_local_writer")
                    exhausted = runtime.execute_one_attempt(
                        self.plan, store, checkpoint, copy.deepcopy(self.reconstructed.members[key]), key=key,
                        runner=self._retry_stub(f"{owner}_retry_exhausted", exhausted=True),
                    )
                    self.assertEqual(exhausted.disposition, "invalid_terminal")
                    self.assertTrue(store.recover().terminal)

    def test_scientific_and_invalid_stops_close_terminal_checkpoints(self):
        for disposition, expected in (("scientific_stop", "scientific_terminal"), ("invalid", "invalid_terminal")):
            with self.subTest(disposition=disposition):
                temporary, store, checkpoint = self.bridge()
                with temporary:
                    result = runtime.execute_one_attempt(
                        self.plan, store, checkpoint, copy.deepcopy(self.reconstructed.members["RK4-2049"]),
                        key="RK4-2049", runner=self._terminal_stub(disposition),
                    )
                    self.assertEqual(result.disposition, expected)
                    status = store.recover()
                    self.assertTrue(status.terminal)
                    self.assertEqual(status.state, "terminal")

    def test_all_six_synthetically_advanced_exact_states_close_the_common_event(self):
        temporary, store, checkpoint = self.bridge()
        with temporary, patch.object(runtime, "choose_requested_cap", side_effect=lambda _c, _k, item: 1.5 - item.time):
            current = checkpoint
            for key in runtime.MEMBER_KEYS:
                result = runtime.execute_one_attempt(
                    self.plan, store, current, copy.deepcopy(self.reconstructed.members[key]),
                    key=key, runner=self._accepted_stub,
                )
                current = result.checkpoint
            event = runtime.commit_common_event(current, store)
            self.assertEqual(event.disposition, "event_complete")
            self.assertEqual(event.checkpoint.event, 24)
            self.assertEqual(event.checkpoint.target["rational"], "25/16")
            self.assertEqual(store.recover().state, "live_local_writer")

    def test_accepted_publication_cuts_have_one_recovery_path(self):
        for phase in (
            "state:after_link",
            "journal:after_link",
            "checkpoint:after_link",
        ):
            with self.subTest(phase=phase):
                temporary, store, checkpoint = self.bridge()
                with temporary:
                    store.fault_hook = lambda observed: (
                        (_ for _ in ()).throw(RuntimeError(phase))
                        if observed == phase else None
                    )
                    with patch.object(runtime, "choose_requested_cap", side_effect=lambda _c, _k, item: 1.5 - item.time):
                        with self.assertRaises(RuntimeError):
                            runtime.execute_one_attempt(
                                self.plan, store, checkpoint,
                                copy.deepcopy(self.reconstructed.members["RK4-2049"]),
                                key="RK4-2049", runner=self._accepted_stub,
                            )
                    store.fault_hook = None
                    decision = reconcile(store)
                    self.assertEqual(decision.disposition, "publication_replay_required")
                    # Explicit store recovery removes only the proven linked
                    # stage, never a one-link predecessor.  What remains then
                    # determines the unique continuation.
                    status = store.recover()
                    snapshot = store.authenticated_snapshot()
                    if snapshot.suffix_classification == "payload_orphan":
                        expected_semantic = snapshot.orphan_payload_semantic_sha256
                        self.assertIsNotNone(expected_semantic)
                        with patch.object(
                            runtime,
                            "choose_requested_cap",
                            side_effect=lambda _c, _k, item: 1.5 - item.time,
                        ):
                            runtime.execute_one_attempt(
                                self.plan,
                                store,
                                checkpoint,
                                copy.deepcopy(self.reconstructed.members["RK4-2049"]),
                                key="RK4-2049",
                                runner=self._accepted_stub,
                                expected_orphan_semantic_sha256=expected_semantic,
                            )
                    elif snapshot.suffix_classification != "clean":
                        self.assertEqual(
                            reconcile(store).disposition,
                            "reconciled_checkpoint",
                        )
                    self.assertEqual(store.recover().state, "live_local_writer")


if __name__ == "__main__":
    unittest.main()
