"""Focused publication-only recovery decisions; no numerical member is run."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
import unittest

from recursive_horizons.fgc.evolution import hlt16_campaign_recovery as recovery
from recursive_horizons.fgc.evolution import hlt16_lifecycle as lifecycle
from recursive_horizons.fgc.evolution import proto15_runtime as p15
from recursive_horizons.fgc.evolution.hlt16_campaign_schema import (
    CampaignCheckpoint,
    MemberCampaignState,
    journal_envelope,
    validate_checkpoint_successor,
)
from recursive_horizons.fgc.evolution.numerical_engine import array_content_sha256
from recursive_horizons.fgc.evolution.proto17_pure_construction import MEMBER_KEYS
import numpy as np


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


class _Store:
    def __init__(self, snapshot):
        self.snapshot = snapshot
        self.published: list[object] = []

    def authenticated_snapshot(self):
        return self.snapshot

    def active_writer_capability(self, _checkpoint):
        return object()

    def publish_journal(self, record, **_kwargs):
        self.published.append(("journal", record))

    def publish_checkpoint(self, current, *, previous, records, **_kwargs):
        self.published.append(("checkpoint", current, previous, records))

    def publish_terminal_lock(self, checkpoint, **_kwargs):
        self.published.append(("terminal_lock", checkpoint))


class _TDG6SuffixStore(_Store):
    """Synthetic store exposing one persisted descriptor payload only."""

    def __init__(self, snapshot, arrays):
        super().__init__(snapshot)
        self.arrays = arrays
        self.load_calls: list[str] = []

    def load_state(self, descriptor_sha256):
        self.load_calls.append(descriptor_sha256)
        return SimpleNamespace(arrays=self.arrays)


class HLT16CampaignRecoveryDecisionTests(unittest.TestCase):
    def _checkpoint(self):
        return SimpleNamespace(
            sha256="b" * 64,
            members={"RK4-2049": _State()},
            event=23,
            target={"rational": "3/2", "binary64_hex": (1.5).hex()},
            plan_sha256="c" * 64, authorization_commit="d" * 40,
            protocol="FGC-2-SF1-PROTO18", campaign_id="campaign",
            journal_tip_sha256="e" * 64, generation=1, journal_sequence=0,
            disposition="nonterminal",
        )

    def _snapshot(self, **changes):
        value = dict(
            checkpoint=self._checkpoint(), suffix_records=(), suffix_kinds=(),
            suffix_classification="clean", orphan_descriptor_sha256=None,
            orphan_payload_semantic_sha256=None, staging_paths=(),
            terminal_lock_present=False,
        )
        value.update(changes)
        return SimpleNamespace(**value)

    def test_clean_snapshot_is_a_read_only_noop(self):
        store = _Store(self._snapshot())
        result = recovery.reconcile(store)
        self.assertEqual(result.disposition, "clean_checkpoint")
        self.assertFalse(result.published)
        self.assertEqual(store.published, [])

    def test_payload_without_descriptor_requires_rerun_not_guessing(self):
        store = _Store(self._snapshot(
            suffix_classification="payload_orphan",
            orphan_payload_semantic_sha256="f" * 64,
        ))
        result = recovery.reconcile(store)
        self.assertEqual(result.disposition, "rerun_from_checkpoint_required")
        self.assertFalse(result.published)
        self.assertEqual(store.published, [])

    def test_any_visible_stage_requires_explicit_publication_replay(self):
        store = _Store(self._snapshot(staging_paths=("journal/.cut.hlt16-stage-deadbeef",)))
        result = recovery.reconcile(store)
        self.assertEqual(result.disposition, "publication_replay_required")
        self.assertFalse(result.published)
        self.assertEqual(store.published, [])

    def test_unrecognized_authenticated_suffix_is_invalid_not_reinterpreted(self):
        store = _Store(self._snapshot(
            suffix_classification="source_rejection+accepted_fine",
            suffix_kinds=("source_rejection", "accepted_fine"),
            suffix_records=({"record_sha256": "f" * 64}, {"record_sha256": "0" * 64}),
        ))
        result = recovery.reconcile(store)
        self.assertEqual(result.disposition, "invalid_ambiguous_suffix")
        self.assertFalse(result.published)
        self.assertEqual(store.published, [])

    def test_clean_terminal_without_its_store_lock_publishes_only_that_lock(self):
        checkpoint = self._checkpoint()
        checkpoint.disposition = "invalid_terminal"
        store = _Store(self._snapshot(checkpoint=checkpoint))
        result = recovery.reconcile(store)
        self.assertEqual(result.disposition, "reconciled_checkpoint")
        self.assertTrue(result.published)
        self.assertEqual(store.published, [("terminal_lock", checkpoint)])

    @staticmethod
    def _tdg6_cursor(
        key: str,
        descriptor: str,
        ledger: p15.TDG6TemporalLedger,
        *,
        accepted_time: float = 23 / 16,
        accepted_rational: str = "23/16",
    ):
        method, points = key.split("-", 1)
        return p15._cursor(
            {
                "protocol_artifact_id": "FGC-2-SF1-PROTO18",
                "campaign_id": "campaign",
                "member_key": key,
                "method": method,
                "point_count": int(points),
                "committed_common_event_index": 23,
                "accepted_boundary_time": {
                    "rational": accepted_rational,
                    "binary64_hex": accepted_time.hex(),
                },
                "accepted_state_sha256": descriptor,
                "previous_step_index": 7,
                "previous_transaction_serial": 11,
                "TDG6_ledger_sha256": p15._hash(p15._ledger_mapping(ledger)),
                "cursor_generation": 0,
                "event_target_time": {
                    "rational": "3/2", "binary64_hex": (3 / 2).hex(),
                },
                "attempt_serial": 0,
                "attempt_id": f"campaign:{key}:23:0",
                "mode": "FRESH_READY",
                "retry_successor_payload_or_none": None,
                "journal_tip_sha256": "e" * 64,
                "cursor_chain_parent_sha256": "f" * 64,
            }
        )

    def _tdg6_rejection_fixture(
        self,
        *,
        mutate_physical_identity: bool = False,
        accepted_time: float = 23 / 16,
        accepted_rational: str = "23/16",
        requested_cap: float = 1 / 64,
    ):
        """One real six-member checkpoint plus an uncheckpointed TDG6 suffix."""
        descriptor = "a" * 64
        ledger = p15.TDG6TemporalLedger.zero(initial_time=accepted_time)
        members = {}
        for key in MEMBER_KEYS:
            cursor = self._tdg6_cursor(
                key,
                descriptor,
                ledger,
                accepted_time=accepted_time,
                accepted_rational=accepted_rational,
            )
            members[key] = MemberCampaignState(
                descriptor_sha256=descriptor,
                cursor=dict(cursor.payload),
                ledger=p15._ledger_mapping(ledger),
                cfl_current=1 if key == "RK4-2049" else 0,
                cfl_total=1 if key == "RK4-2049" else 0,
                pending_owner="cfl" if key == "RK4-2049" else None,
                pending_cap_hex=requested_cap.hex() if key == "RK4-2049" else None,
            )
        checkpoint = CampaignCheckpoint(
            plan_sha256="c" * 64,
            authorization_commit="d" * 40,
            protocol="FGC-2-SF1-PROTO18",
            campaign_id="campaign",
            event=23,
            target={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
            members=members,
            journal_tip_sha256="e" * 64,
            parent_sha256="f" * 64,
            generation=7,
            journal_sequence=6,
            disposition="nonterminal",
        )
        arrays = {
            "u": np.array([[1.0]], dtype=np.float64),
            "p": np.array([[2.0]], dtype=np.float64),
            "q": np.array([[3.0]], dtype=np.float64),
        }
        physical = array_content_sha256(arrays["u"], arrays["p"], arrays["q"])
        cursor = p15.Proto15Cursor(members["RK4-2049"].cursor)
        plan = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(
                accepted_time,
                3 / 2,
                requested_cap,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
        )
        admission = p15._synthetic_failed_admission(
            p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]]
        )
        evidence = p15._json_safe(p15.TDG6TemporalRetryEvidence(
            event_type="rejected_TDG6_temporal_admission",
            method=p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]],
            initial_time=float.fromhex(plan["current_hex"]),
            attempted_macro_step_size=float.fromhex(plan["macro_width_hex"]),
            retry_step_size=float.fromhex(plan["macro_width_hex"]) * 0.5,
            initial_state_sha256=("b" * 64 if mutate_physical_identity else physical),
            previous_step_index=cursor.payload["previous_step_index"],
            previous_transaction_serial=cursor.payload["previous_transaction_serial"],
            retry_count_for_current_macro_step=1,
            cumulative_temporal_retry_count=1,
            failed_channels=admission.failed_channels,
            channel_admissions=admission.channel_admissions,
        ).as_mapping())
        record = journal_envelope(
            sequence=7,
            campaign_id="campaign",
            generation=8,
            event=23,
            previous_record_sha256="e" * 64,
            kind="tdg6_rejection",
            payload={
                "member_key": "RK4-2049",
                "predecessor_descriptor_sha256": descriptor,
                "predecessor_cursor_sha256": cursor.sha256,
                "evidence": evidence,
            },
        )
        snapshot = SimpleNamespace(
            checkpoint=checkpoint,
            suffix_records=(record,),
            suffix_kinds=("tdg6_rejection",),
            suffix_classification="tdg6_rejection",
            orphan_descriptor_sha256=None,
            orphan_payload_semantic_sha256=None,
            staging_paths=(),
            terminal_lock_present=False,
        )
        return _TDG6SuffixStore(snapshot, arrays), checkpoint, descriptor, physical

    def test_pending_rejection_preserves_live_shape_conservative_plan_exactly(self):
        """A second rejection must retain the durable requested-cap history."""
        accepted = float.fromhex("0x1.78554de5a30e0p+0")
        store, _previous, descriptor, physical = self._tdg6_rejection_fixture(
            accepted_time=accepted,
            accepted_rational="206891375579527/140737488355328",
            requested_cap=float.fromhex("0x1.aaa9612df9000p-9"),
        )
        first = recovery.preview_tdg6_recovery(store)
        previous = first.checkpoint
        before = previous.members["RK4-2049"]
        cursor = p15.Proto15Cursor(dict(before.cursor))
        ledger = p15._ledger_from_mapping(before.ledger)
        pending = cursor.payload["retry_successor_payload_or_none"]
        self.assertIsInstance(pending, dict)
        durable_predecessor = pending["successor_plan"]
        self.assertEqual(
            durable_predecessor["requested_cap_hex"],
            "0x1.aaa9612df9000p-10",
        )
        self.assertEqual(
            durable_predecessor["macro_width_hex"],
            "0x1.aaa9612df8000p-10",
        )
        self.assertEqual(
            durable_predecessor["conservative_reduction_hex"],
            "0x1.0000000000000p-50",
        )

        # This is the lossy reconstruction used before this regression: the
        # selected width becomes a new requested cap and the one-quantum debit
        # disappears, so it is not the plan the pending cursor authorized.
        lossy = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(
                accepted,
                3 / 2,
                float.fromhex(durable_predecessor["macro_width_hex"]),
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
        )
        self.assertNotEqual(lossy, durable_predecessor)
        self.assertEqual(lossy["conservative_reduction_hex"], "0x0.0p+0")

        admission = p15._synthetic_failed_admission(
            p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]]
        )
        attempted = float.fromhex(durable_predecessor["macro_width_hex"])
        evidence = p15._json_safe(
            p15.TDG6TemporalRetryEvidence(
                event_type="rejected_TDG6_temporal_admission",
                method=p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]],
                initial_time=accepted,
                attempted_macro_step_size=attempted,
                retry_step_size=attempted * 0.5,
                initial_state_sha256=physical,
                previous_step_index=cursor.payload["previous_step_index"],
                previous_transaction_serial=cursor.payload[
                    "previous_transaction_serial"
                ],
                retry_count_for_current_macro_step=2,
                cumulative_temporal_retry_count=2,
                failed_channels=admission.failed_channels,
                channel_admissions=admission.channel_admissions,
            ).as_mapping()
        )
        rejection = journal_envelope(
            sequence=previous.journal_sequence + 1,
            campaign_id=previous.campaign_id,
            generation=previous.generation + 1,
            event=previous.event,
            previous_record_sha256=previous.journal_tip_sha256,
            kind="tdg6_rejection",
            payload={
                "member_key": "RK4-2049",
                "predecessor_descriptor_sha256": descriptor,
                "predecessor_cursor_sha256": cursor.sha256,
                "evidence": evidence,
            },
        )
        second_store = _TDG6SuffixStore(
            SimpleNamespace(
                checkpoint=previous,
                suffix_records=(rejection,),
                suffix_kinds=("tdg6_rejection",),
                suffix_classification="tdg6_rejection",
                orphan_descriptor_sha256=None,
                orphan_payload_semantic_sha256=None,
                staging_paths=(),
                terminal_lock_present=False,
            ),
            store.arrays,
        )
        preview = recovery.preview_tdg6_recovery(second_store)

        typed = lifecycle.typed_retry_from_evidence(
            cursor,
            ledger,
            evidence,
            evolution_state_sha256=physical,
        )
        exact_successor_plan = p15.plan_forward_proto14_subdivision(
            accepted,
            3 / 2,
            typed.evidence.retry_step_size,
            minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
        )
        exact_cursor = lifecycle.close_retry(
            cursor,
            ledger,
            typed,
            durable_predecessor,
            exact_successor_plan,
            rejection["record_sha256"],
            previous.journal_tip_sha256,
            evolution_state_sha256=physical,
        )
        transition = preview.records[1]
        self.assertEqual(transition["sequence"], 10)
        self.assertEqual(transition["payload"]["successor_cursor"], exact_cursor.payload)
        self.assertEqual(
            transition["payload"]["successor_cursor_sha256"], exact_cursor.sha256
        )
        self.assertEqual(
            transition["payload"]["successor_cursor"][
                "retry_successor_payload_or_none"
            ]["predecessor_plan"],
            durable_predecessor,
        )
        self.assertEqual(preview.checkpoint.generation, 9)
        self.assertEqual(preview.checkpoint.journal_sequence, 10)
        self.assertEqual(preview.evolution_state_sha256, physical)
        self.assertFalse(preview.physical_state_advanced)
        validate_checkpoint_successor(
            previous,
            preview.checkpoint,
            preview.records,
        )

    def test_uncheckpointed_tdg6_rejection_recovers_with_distinct_descriptor_and_physical_identity(self):
        store, previous, descriptor, physical = self._tdg6_rejection_fixture()
        preview = recovery.preview_tdg6_recovery(store)
        self.assertEqual(store.published, [])
        self.assertFalse(preview.physical_state_advanced)
        self.assertEqual(preview.predecessor_checkpoint_sha256, previous.sha256)
        self.assertEqual(preview.evolution_state_sha256, physical)
        self.assertEqual(preview.records[-1]["kind"], "cursor_transition")
        result = recovery.reconcile(store)
        self.assertEqual(result.disposition, "reconciled_checkpoint")
        self.assertTrue(result.published)
        self.assertEqual(store.load_calls, [descriptor, descriptor])
        self.assertEqual([item[0] for item in store.published], ["journal", "checkpoint"])
        transition = store.published[0][1]
        current = store.published[1][1]
        self.assertEqual(current.sha256, preview.checkpoint.sha256)
        self.assertEqual(transition["kind"], "cursor_transition")
        state = current.members["RK4-2049"]
        self.assertEqual(state.descriptor_sha256, descriptor)
        self.assertEqual(
            state.cursor["accepted_boundary_time"],
            previous.members["RK4-2049"].cursor["accepted_boundary_time"],
        )
        self.assertEqual((state.cfl_current, state.cfl_total), (1, 1))
        self.assertEqual(state.ledger["cumulative_temporal_retry_count"], 1)
        self.assertEqual(state.pending_owner, "temporal")
        self.assertEqual(state.pending_cap_hex, (1 / 128).hex())
        self.assertNotEqual(descriptor, physical)
        fresh_evidence = preview.records[0]["payload"]["evidence"]
        expected_fresh_predecessor = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(
                float(fresh_evidence["initial_time"]),
                3 / 2,
                float(fresh_evidence["attempted_macro_step_size"]),
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
        )
        self.assertEqual(
            state.cursor["retry_successor_payload_or_none"]["predecessor_plan"],
            expected_fresh_predecessor,
        )

    def test_uncheckpointed_tdg6_rejection_rejects_physical_hash_mutation_without_publication(self):
        store, _previous, _descriptor, _physical = self._tdg6_rejection_fixture(
            mutate_physical_identity=True
        )
        with self.assertRaises(recovery.HLT16RecoveryError):
            recovery.reconcile(store)
        self.assertEqual(store.published, [])


if __name__ == "__main__":
    unittest.main()
