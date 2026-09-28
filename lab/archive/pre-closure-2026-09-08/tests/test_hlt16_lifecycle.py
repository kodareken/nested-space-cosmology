"""Adversarial unit tests for the pure HLT16 lifecycle authority."""
from __future__ import annotations

import ast
import importlib
import inspect
import json
import unittest
from dataclasses import replace
from fractions import Fraction

from recursive_horizons.fgc.evolution import proto15_runtime as p15
from recursive_horizons.fgc.evolution.hlt16_lifecycle import (
    HLT16LifecycleError, PROTO17_ARTIFACT_ID, PROTO18_ARTIFACT_ID,
    bridge_gen0, canonical_hash, canonical_mapping, close_fine,
    close_retry as _close_retry, retry_exhaustion as _retry_exhaustion,
    typed_retry_from_evidence,
)


MEMBER, STATE, DESCRIPTOR, REJECTION = "RK4-2049", "1" * 64, "2" * 64, "3" * 64
CHECKPOINT_TIP = "e" * 64


def close_retry(*args, **kwargs):
    kwargs["evolution_state_sha256"] = STATE
    return _close_retry(*args, **kwargs)


def retry_exhaustion(*args, **kwargs):
    kwargs["evolution_state_sha256"] = STATE
    return _retry_exhaustion(*args, **kwargs)


def zero() -> p15.TDG6TemporalLedger:
    return p15.TDG6TemporalLedger.zero(initial_time=23 / 16)


def raw(ledger: p15.TDG6TemporalLedger | None = None) -> dict[str, object]:
    ledger = zero() if ledger is None else ledger
    return {
        "protocol_artifact_id": PROTO17_ARTIFACT_ID, "campaign_id": "HLT16-TEST",
        "member_key": MEMBER, "method": "RK4", "point_count": 2049,
        "committed_common_event_index": 23,
        "accepted_boundary_time": {"rational": "23/16", "binary64_hex": (23 / 16).hex()},
        "accepted_state_sha256": STATE, "previous_step_index": 7,
        "previous_transaction_serial": 11,
        "TDG6_ledger_sha256": p15._hash(p15._ledger_mapping(ledger)),
        "cursor_generation": 0,
        "event_target_time": {"rational": "3/2", "binary64_hex": (3 / 2).hex()},
        "attempt_serial": 0, "attempt_id": "HLT16-TEST:RK4-2049:23:0",
        "mode": "FRESH_READY", "retry_successor_payload_or_none": None,
        "journal_tip_sha256": p15.ROOT_SHA256, "cursor_chain_parent_sha256": p15.ROOT_SHA256,
    }


def proto18() -> tuple[p15.Proto15Cursor, p15.TDG6TemporalLedger]:
    ledger = zero()
    return bridge_gen0(p15._cursor(raw(ledger)), ledger, DESCRIPTOR), ledger


def plan(cursor: p15.Proto15Cursor, cap: float = 1 / 16) -> dict[str, object]:
    return p15._plan_mapping(p15.plan_forward_proto14_subdivision(
        float.fromhex(cursor.payload["accepted_boundary_time"]["binary64_hex"]),
        float.fromhex(cursor.payload["event_target_time"]["binary64_hex"]), cap,
        minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
    ))


def time_identity(value: float) -> dict[str, str]:
    return {"rational": str(Fraction.from_float(value)), "binary64_hex": value.hex()}


def accepted(prior: p15.TDG6TemporalLedger, endpoint: float) -> p15.TDG6TemporalLedger:
    return p15.TDG6TemporalLedger(
        last_accepted_time=endpoint,
        accepted_macro_step_count=prior.accepted_macro_step_count + 1,
        cumulative_temporal_retry_count=prior.cumulative_temporal_retry_count,
        current_macro_step_temporal_retry_count=0,
        last_accepted_macro_step_temporal_retry_count=prior.current_macro_step_temporal_retry_count,
        accumulated_debit_vector=prior.accumulated_debit_vector,
        serialized_temporal_rejections=prior.serialized_temporal_rejections,
    )


def rejected(cursor: p15.Proto15Cursor, ledger: p15.TDG6TemporalLedger, old_plan: dict[str, object], *, exhausted: bool = False) -> object:
    admission = p15._synthetic_failed_admission(p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]])
    minimum_exhaustion = exhausted and float.fromhex(old_plan["macro_width_hex"]) <= float(p15.TDG6_MINIMUM_MACRO_STEP)
    count = (ledger.current_macro_step_temporal_retry_count + 1 if minimum_exhaustion
             else (p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP + 1 if exhausted
                   else ledger.current_macro_step_temporal_retry_count + 1))
    evidence = p15.TDG6TemporalRetryEvidence(
        event_type="rejected_TDG6_temporal_admission", method=p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]],
        initial_time=float.fromhex(old_plan["current_hex"]), attempted_macro_step_size=float.fromhex(old_plan["macro_width_hex"]),
        retry_step_size=float.fromhex(old_plan["macro_width_hex"]) * .5, initial_state_sha256=STATE,
        previous_step_index=cursor.payload["previous_step_index"], previous_transaction_serial=cursor.payload["previous_transaction_serial"],
        retry_count_for_current_macro_step=count, cumulative_temporal_retry_count=ledger.cumulative_temporal_retry_count + 1,
        failed_channels=admission.failed_channels, channel_admissions=admission.channel_admissions,
    )
    updated = p15._updated_retry_ledger(ledger, evidence.as_mapping())
    if exhausted:
        return p15.TDG6TemporalRetryExhausted(
            reason="maximum_temporal_retries" if count > p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP else "minimum_macro_step",
            evidence=evidence,
            updated_ledger=updated,
        )
    return p15.TDG6TemporalRetryRequired(evidence, updated)


class TestHLT16Lifecycle(unittest.TestCase):
    def test_pure_owner_has_no_io_or_solver_import(self) -> None:
        module = importlib.import_module("recursive_horizons.fgc.evolution.hlt16_lifecycle")
        source = ast.unparse(ast.parse(inspect.getsource(module)))
        for forbidden in ("Path(", "open(", "numpy", "numerical_engine", "FGC-QR", "SGB-L"):
            self.assertNotIn(forbidden, source)

    def test_bridge_exact_and_rejects_rebridging(self) -> None:
        old = p15._cursor(raw())
        new, ledger = proto18()
        self.assertEqual(new.payload["protocol_artifact_id"], PROTO18_ARTIFACT_ID)
        self.assertEqual(new.payload["accepted_state_sha256"], DESCRIPTOR)
        self.assertEqual(new.payload["cursor_chain_parent_sha256"], old.sha256)
        self.assertEqual(new.mode, "FRESH_READY")
        self.assertEqual(canonical_mapping(new)["attempt_serial"], 0)
        with self.assertRaises(HLT16LifecycleError):
            bridge_gen0(new, ledger, "4" * 64)

    def test_bridge_can_rebase_campaign_without_rewriting_authenticated_parent(self) -> None:
        old = p15._cursor(raw())
        ledger = zero()
        new = bridge_gen0(
            old,
            ledger,
            DESCRIPTOR,
            successor_campaign_id="HLT16-SUCCESSOR",
        )
        self.assertEqual(new.payload["campaign_id"], "HLT16-SUCCESSOR")
        self.assertEqual(
            new.payload["attempt_id"],
            "HLT16-SUCCESSOR:RK4-2049:23:0",
        )
        self.assertEqual(new.payload["cursor_chain_parent_sha256"], old.sha256)
        self.assertEqual(old.payload["campaign_id"], "HLT16-TEST")
        self.assertEqual(
            bridge_gen0(
                old,
                ledger,
                DESCRIPTOR,
                successor_campaign_id="HLT16-TEST",
            ).payload,
            bridge_gen0(old, ledger, DESCRIPTOR).payload,
        )
        mutated = dict(old.payload)
        mutated["cursor_chain_sha256"] = "f" * 64
        with self.assertRaises(HLT16LifecycleError):
            bridge_gen0(
                mutated,
                ledger,
                DESCRIPTOR,
                successor_campaign_id="HLT16-SUCCESSOR",
            )
        for bad in ("", "bad:campaign", 1):
            with self.subTest(bad=bad):
                with self.assertRaises(HLT16LifecycleError):
                    bridge_gen0(
                        old,
                        ledger,
                        DESCRIPTOR,
                        successor_campaign_id=bad,
                    )

    def test_bridge_rejects_boundary_and_zero_ledger_mutations(self) -> None:
        ledger = zero()
        for field, value in (("committed_common_event_index", 24), ("journal_tip_sha256", "f" * 64), ("cursor_generation", 1)):
            with self.subTest(field=field):
                changed = raw(ledger); changed[field] = value
                with self.assertRaises(HLT16LifecycleError):
                    bridge_gen0(changed, ledger, DESCRIPTOR)
        dirty = replace(ledger, accepted_macro_step_count=1)
        with self.assertRaises(HLT16LifecycleError):
            bridge_gen0(raw(dirty), dirty, DESCRIPTOR)

    def test_direct_accepted_fine_and_mutations(self) -> None:
        cursor, prior = proto18(); executed = plan(cursor)
        endpoint = float.fromhex(executed["boundaries_hex"][-1]); ledger = accepted(prior, endpoint)
        new = close_fine(cursor, prior, ledger, "4" * 64, executed, time_identity(endpoint), 8, 12, CHECKPOINT_TIP)
        self.assertEqual(new.mode, "FRESH_READY")
        self.assertEqual(new.payload["accepted_state_sha256"], "4" * 64)
        self.assertEqual(new.payload["cursor_chain_parent_sha256"], cursor.sha256)
        self.assertEqual(new.payload["journal_tip_sha256"], CHECKPOINT_TIP)
        self.assertNotEqual(canonical_hash(new), canonical_hash(cursor))
        cases = (
            lambda: close_fine(cursor, prior, ledger, DESCRIPTOR, executed, time_identity(endpoint), 8, 12, CHECKPOINT_TIP),
            lambda: close_fine(cursor, prior, ledger, "4" * 64, executed, time_identity(endpoint), 7, 12, CHECKPOINT_TIP),
            lambda: close_fine(cursor, prior, ledger, "4" * 64, executed, time_identity(endpoint), 8, 11, CHECKPOINT_TIP),
            lambda: close_fine(cursor, prior, replace(ledger, accepted_macro_step_count=9), "4" * 64, executed, time_identity(endpoint), 8, 12, CHECKPOINT_TIP),
            lambda: close_fine(cursor, prior, ledger, "4" * 64, plan(cursor, 1 / 32), time_identity(endpoint), 8, 12, CHECKPOINT_TIP),
            lambda: close_fine(cursor, prior, ledger, "4" * 64, executed, time_identity(endpoint), 8, 12, p15.ROOT_SHA256),
        )
        for action in cases:
            with self.assertRaises(HLT16LifecycleError):
                action()

    def test_retry_and_retry_fine_with_plan_journal_and_ledger_attacks(self) -> None:
        cursor, prior = proto18(); old_plan = plan(cursor); outcome = rejected(cursor, prior, old_plan)
        assert isinstance(outcome, p15.TDG6TemporalRetryRequired)
        successor = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            outcome.evidence.initial_time, 3 / 2, outcome.retry_step_size, minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP)))
        pending = close_retry(cursor, prior, outcome, old_plan, successor, REJECTION, p15.ROOT_SHA256)
        self.assertEqual(pending.mode, "RETRY_PENDING")
        self.assertEqual(pending.payload["journal_tip_sha256"], REJECTION)
        endpoint = float.fromhex(successor["boundaries_hex"][-1]); post = accepted(outcome.updated_ledger, endpoint)
        fresh = close_fine(pending, outcome.updated_ledger, post, "5" * 64, successor, time_identity(endpoint), 8, 12, CHECKPOINT_TIP)
        self.assertEqual(fresh.mode, "FRESH_READY")
        self.assertEqual(fresh.payload["journal_tip_sha256"], CHECKPOINT_TIP)
        corrupt = dict(successor); corrupt["requested_cap_hex"] = float(outcome.retry_step_size * .5).hex()
        for args in ((old_plan, corrupt, REJECTION, p15.ROOT_SHA256), (old_plan, successor, p15.ROOT_SHA256, p15.ROOT_SHA256)):
            with self.assertRaises(HLT16LifecycleError):
                close_retry(cursor, prior, outcome, *args)
        # The checkpoint journal tip is authoritative and may be newer than
        # the cursor's embedded acyclic ancestor.
        newer = close_retry(
            cursor, prior, outcome, old_plan, successor, REJECTION, "f" * 64
        )
        self.assertEqual(
            newer.payload["retry_successor_payload_or_none"][
                "predecessor_journal_tip_sha256"
            ],
            "f" * 64,
        )
        bad_ledger = replace(outcome.updated_ledger, accumulated_debit_vector=(1.0,) * 18)
        bad = p15.TDG6TemporalRetryRequired(outcome.evidence, bad_ledger)
        with self.assertRaises(HLT16LifecycleError):
            close_retry(cursor, prior, bad, old_plan, successor, REJECTION, p15.ROOT_SHA256)

    def test_persisted_retry_evidence_recursively_decodes_exactly(self) -> None:
        """JSON round-trip preserves a descriptor-distinct TDG6 field identity."""
        cursor, prior = proto18()
        outcome = rejected(cursor, prior, plan(cursor))
        assert isinstance(outcome, p15.TDG6TemporalRetryRequired)
        persisted = json.loads(json.dumps(outcome.evidence.as_mapping()))

        rebuilt = typed_retry_from_evidence(
            cursor,
            prior,
            persisted,
            evolution_state_sha256=STATE,
        )
        self.assertIsInstance(rebuilt, p15.TDG6TemporalRetryRequired)
        self.assertEqual(rebuilt.evidence.as_mapping(), outcome.evidence.as_mapping())
        self.assertEqual(p15._json_safe(rebuilt.evidence.as_mapping()), persisted)
        self.assertEqual(
            p15._canonical(p15._json_safe(rebuilt.evidence.as_mapping())),
            p15._canonical(persisted),
        )
        self.assertEqual(rebuilt.updated_ledger, outcome.updated_ledger)
        self.assertNotEqual(cursor.payload["accepted_state_sha256"], STATE)

        for mutation in ("extra_envelope_field", "wrong_physical_identity"):
            with self.subTest(mutation=mutation):
                attacked = json.loads(json.dumps(persisted))
                physical = STATE
                if mutation == "extra_envelope_field":
                    attacked["channel_admissions"][0]["outer_difference"][
                        "envelope"
                    ]["foreign"] = 1
                else:
                    physical = "4" * 64
                with self.assertRaises(HLT16LifecycleError):
                    typed_retry_from_evidence(
                        cursor,
                        prior,
                        attacked,
                        evolution_state_sha256=physical,
                    )

    def test_retry_exhaustion_is_nonphysical_and_never_pending(self) -> None:
        cursor, prior = proto18()
        old_plan = plan(cursor, float(p15.TDG6_MINIMUM_MACRO_STEP))
        outcome = rejected(cursor, prior, old_plan, exhausted=True)
        assert isinstance(outcome, p15.TDG6TemporalRetryExhausted)
        terminal = retry_exhaustion(cursor, prior, outcome, REJECTION, p15.ROOT_SHA256)
        self.assertEqual(terminal.reason, "minimum_macro_step")
        with self.assertRaises(HLT16LifecycleError):
            close_retry(cursor, prior, outcome, old_plan, old_plan, REJECTION, p15.ROOT_SHA256)  # type: ignore[arg-type]

    def test_second_pending_rejection_preserves_prefix_and_requires_exact_plan(self) -> None:
        cursor, prior = proto18()
        first_plan = plan(cursor, 1 / 16)
        first = rejected(cursor, prior, first_plan)
        assert isinstance(first, p15.TDG6TemporalRetryRequired)
        first_successor = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            first.evidence.initial_time, 3 / 2, first.retry_step_size,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ))
        pending_one = close_retry(
            cursor, prior, first, first_plan, first_successor, REJECTION,
            p15.ROOT_SHA256,
        )
        second = rejected(pending_one, first.updated_ledger, first_successor)
        assert isinstance(second, p15.TDG6TemporalRetryRequired)
        second_successor = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            second.evidence.initial_time, 3 / 2, second.retry_step_size,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ))
        pending_two = close_retry(
            pending_one, first.updated_ledger, second, first_successor,
            second_successor, "4" * 64, REJECTION,
        )
        prefix = pending_two.payload["retry_successor_payload_or_none"]["complete_rejection_prefix"]
        self.assertEqual(tuple(prefix), second.updated_ledger.serialized_temporal_rejections)
        self.assertEqual(len(prefix), 2)
        self.assertEqual(pending_two.payload["journal_tip_sha256"], "4" * 64)
        with self.assertRaises(HLT16LifecycleError):
            close_retry(
                pending_one, first.updated_ledger, second, first_plan,
                second_successor, "5" * 64, REJECTION,
            )

    def test_pending_retry_exhaustion_replays_exact_successor_without_commit(self) -> None:
        """A second TDG6 rejection may stop from one authenticated pending edge.

        Use the exact minimum-width successor: the first rejection is still
        retryable at ``2 * minimum`` and the second one exhausts because its
        frozen half-step would fall below the minimum.  The terminal evidence
        must retain the original accepted descriptor and time.
        """
        cursor, prior = proto18()
        old_plan = plan(cursor, 2.0 * float(p15.TDG6_MINIMUM_MACRO_STEP))
        first = rejected(cursor, prior, old_plan)
        assert isinstance(first, p15.TDG6TemporalRetryRequired)
        successor = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            first.evidence.initial_time, 3 / 2, first.retry_step_size,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ))
        pending = close_retry(cursor, prior, first, old_plan, successor, REJECTION, p15.ROOT_SHA256)
        self.assertEqual(pending.mode, "RETRY_PENDING")
        exhausted = rejected(pending, first.updated_ledger, successor, exhausted=True)
        assert isinstance(exhausted, p15.TDG6TemporalRetryExhausted)
        terminal = retry_exhaustion(
            pending, first.updated_ledger, exhausted, "4" * 64, REJECTION,
        )
        self.assertEqual(terminal.reason, "minimum_macro_step")
        self.assertEqual(terminal.accepted_state_sha256, cursor.payload["accepted_state_sha256"])
        self.assertEqual(terminal.accepted_time, cursor.payload["accepted_boundary_time"])
        self.assertEqual(
            terminal.updated_ledger.accepted_macro_step_count,
            first.updated_ledger.accepted_macro_step_count,
        )
        self.assertEqual(
            terminal.updated_ledger.last_accepted_time,
            first.updated_ledger.last_accepted_time,
        )
        self.assertEqual(
            terminal.updated_ledger.accumulated_debit_vector,
            first.updated_ledger.accumulated_debit_vector,
        )
        self.assertEqual(
            terminal.updated_ledger.current_macro_step_temporal_retry_count,
            first.updated_ledger.current_macro_step_temporal_retry_count + 1,
        )

    def test_pending_retry_exhaustion_rejects_fabricated_plan_or_evidence(self) -> None:
        cursor, prior = proto18()
        old_plan = plan(cursor, 2.0 * float(p15.TDG6_MINIMUM_MACRO_STEP))
        first = rejected(cursor, prior, old_plan)
        assert isinstance(first, p15.TDG6TemporalRetryRequired)
        successor = p15._plan_mapping(p15.plan_forward_proto14_subdivision(
            first.evidence.initial_time, 3 / 2, first.retry_step_size,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        ))
        pending = close_retry(cursor, prior, first, old_plan, successor, REJECTION, p15.ROOT_SHA256)
        exhausted = rejected(pending, first.updated_ledger, successor, exhausted=True)
        assert isinstance(exhausted, p15.TDG6TemporalRetryExhausted)

        # A typed outcome with otherwise valid TDG6 arithmetic is still not a
        # valid pending exhaustion when it claims a different plan origin.
        fabricated_evidence = replace(
            exhausted.evidence,
            initial_time=exhausted.evidence.initial_time + float(p15.TDG6_MINIMUM_MACRO_STEP),
        )
        fabricated = p15.TDG6TemporalRetryExhausted(
            reason=exhausted.reason,
            evidence=fabricated_evidence,
            updated_ledger=p15._updated_retry_ledger(
                first.updated_ledger, fabricated_evidence.as_mapping(),
            ),
        )
        with self.assertRaises(HLT16LifecycleError):
            retry_exhaustion(pending, first.updated_ledger, fabricated, "4" * 64, REJECTION)

        # Proto15Cursor.validate intentionally checks only the cursor envelope;
        # the lifecycle must independently reject a semantically malformed
        # nested retry payload before accepting an exhaustion terminal.
        bad_payload = dict(pending.payload)
        bad_retry = dict(bad_payload["retry_successor_payload_or_none"])
        bad_successor = dict(bad_retry["successor_plan"])
        bad_successor["requested_cap_hex"] = float(
            float.fromhex(str(bad_successor["requested_cap_hex"])) * 0.5
        ).hex()
        bad_retry["successor_plan"] = bad_successor
        bad_payload["retry_successor_payload_or_none"] = bad_retry
        malformed_pending = p15._cursor(bad_payload)
        with self.assertRaises(HLT16LifecycleError):
            retry_exhaustion(malformed_pending, first.updated_ledger, exhausted, "4" * 64, REJECTION)


if __name__ == "__main__":
    unittest.main()
