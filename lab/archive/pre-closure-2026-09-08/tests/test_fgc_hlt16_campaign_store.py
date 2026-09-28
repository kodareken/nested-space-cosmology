from __future__ import annotations

"""Integration tests for the HLT16 authenticated campaign store.

These tests begin with a byte-for-byte copy of sealed PROTO17 generation zero.
A synthetic directory would not establish that the bridge is tied to the real
six authenticated source identities.
"""

import hashlib
import os
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory
import threading
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution import hlt16_lifecycle as lifecycle
from recursive_horizons.fgc.evolution import proto15_runtime as p15
from recursive_horizons.fgc.evolution.hlt16_campaign_schema import MemberCampaignState, journal_envelope
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (
    HLT16CampaignStore,
    HLT16CampaignStoreError,
    _safe_read_snapshot,
    _validated_stage,
    _walk_tree,
)
from recursive_horizons.fgc.evolution.hlt16_member_codec import GEN0_PROTOCOL, RUNTIME_PROTOCOL, encode_member
from recursive_horizons.fgc.evolution.hlt16_state_store import (
    HLT16StateStoreError,
    atomic_publish_noreplace,
    authorize_gen0_import,
    canonical,
    encode_payload,
    read_nofollow,
)
from recursive_horizons.fgc.evolution.numerical_engine import array_content_sha256
from recursive_horizons.fgc.evolution.proto17_pure_construction import MEMBER_KEYS
from tests.fgc_gen0_fixture import copy_sealed_gen0_store, reconstruct_sealed_gen0_members


class HLT16CampaignStoreTests(unittest.TestCase):
    """Exercise recovery against real authenticated GEN0 evidence."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.reconstructed = reconstruct_sealed_gen0_members(ROOT)
        cls.cursors = cls.reconstructed.evidence.checkpoint["cursors"]
        cls.campaign_id = str(cls.cursors[MEMBER_KEYS[0]]["campaign_id"])
        cls.snapshots = {
            key: encode_member(
                cls.reconstructed.members[key], cls.cursors[key],
                protocol_artifact_id=RUNTIME_PROTOCOL,
                campaign_id=cls.campaign_id, generation=0,
                provenance_source_protocol=GEN0_PROTOCOL,
            )
            for key in MEMBER_KEYS
        }
        cls.authorities = {key: authorize_gen0_import(cls.snapshots[key]) for key in MEMBER_KEYS}
        cls.seed_temporary = TemporaryDirectory()
        cls.seed_root = Path(cls.seed_temporary.name) / "calibration"
        copy_sealed_gen0_store(ROOT, cls.seed_root)
        store = HLT16CampaignStore(cls.seed_root)
        cls.seed_checkpoint = cls._adopt(store)
        if store.recover().state != "clean_checkpoint":
            raise AssertionError("HLT16 seed did not close cleanly")

    @classmethod
    def tearDownClass(cls) -> None:
        cls.seed_temporary.cleanup()

    @classmethod
    def _member_factory(cls, key: str, descriptor: str) -> MemberCampaignState:
        ledger = cls.reconstructed.members[key].temporal_ledger
        cursor = lifecycle.bridge_gen0(cls.cursors[key], ledger, descriptor)
        return MemberCampaignState(descriptor_sha256=descriptor, cursor=dict(cursor.payload), ledger=p15._ledger_mapping(ledger))

    @classmethod
    def _adopt(cls, store: HLT16CampaignStore):
        return store.adopt_generation_zero(
            plan_sha256="1" * 64, authorization_commit="2" * 40, campaign_id=cls.campaign_id,
            snapshots=cls.snapshots, authorities=cls.authorities, member_factory=cls._member_factory,
            target={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
        )

    def adopted(self) -> tuple[TemporaryDirectory, HLT16CampaignStore]:
        temporary = TemporaryDirectory()
        root = Path(temporary.name) / "calibration"
        shutil.copytree(self.seed_root, root)
        return temporary, HLT16CampaignStore(root)

    def test_adoption_uses_actual_six_source_identities_and_closes_generation_one(self) -> None:
        self.assertEqual(set(HLT16CampaignStore(self.seed_root)._legacy_inventory()), set(MEMBER_KEYS))
        temporary, adopted = self.adopted()
        with temporary:
            status = adopted.recover()
            self.assertEqual(status.state, "clean_checkpoint")
            self.assertTrue(status.safe_to_restart)
            self.assertEqual(status.checkpoint_generation, 1)
            self.assertEqual(status.checkpoint_sha256, self.seed_checkpoint.sha256)
            self.assertEqual(sorted(path.name for path in adopted.root.iterdir()), ["checkpoints", "journal", "locks", "payloads", "receipts", "states"])
            snapshot = adopted.authenticated_snapshot()
            self.assertEqual(snapshot.checkpoint.sha256, self.seed_checkpoint.sha256)
            self.assertEqual(snapshot.suffix_classification, "clean")
            self.assertEqual(snapshot.suffix_records, ())
            self.assertIsNone(snapshot.orphan_descriptor_sha256)
            self.assertIsNone(snapshot.orphan_payload_semantic_sha256)

    def test_authenticated_historical_checkpoint_accessor_validates_full_store(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            selected = store.authenticated_checkpoint_at_generation(1)
            self.assertEqual(selected.sha256, self.seed_checkpoint.sha256)
            self.assertEqual(selected.generation, 1)
            with self.assertRaisesRegex(
                HLT16CampaignStoreError, "historical checkpoint generation differs",
            ):
                store.authenticated_checkpoint_at_generation(0)
            with self.assertRaisesRegex(
                HLT16CampaignStoreError, "historical checkpoint is unavailable",
            ):
                store.authenticated_checkpoint_at_generation(2)

            # The accessor is not a single-leaf escape hatch: unrelated tree
            # drift invalidates selection even though generation one itself is
            # still byte-for-byte present.
            (store.root / "foreign-root").mkdir()
            with self.assertRaises(HLT16CampaignStoreError):
                store.authenticated_checkpoint_at_generation(1)

    def test_foreign_launcher_cannot_adopt_a_second_generation_one_bridge(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            # The first bridge is already a closed immutable predecessor.  A
            # launcher presenting any different genesis authority must fail;
            # it may not convert this into a second competing bridge.
            with self.assertRaisesRegex(
                HLT16CampaignStoreError, "generation-one bridge authority differs",
            ):
                store.adopt_generation_zero(
                    plan_sha256="a" * 64,
                    authorization_commit="2" * 40,
                    campaign_id=self.campaign_id,
                    snapshots=self.snapshots,
                    authorities=self.authorities,
                    member_factory=self._member_factory,
                    target={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
                )
            self.assertEqual(store.recover().checkpoint_sha256, self.seed_checkpoint.sha256)

    def test_campaign_object_exposes_no_raw_state_persistence_bypass(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            before = _walk_tree(store.root)[0]
            self.assertFalse(hasattr(store, "states"))
            with self.assertRaises(AttributeError):
                getattr(store, "states")
            self.assertEqual(_walk_tree(store.root)[0], before)
            self.assertEqual(store.recover().checkpoint_sha256, self.seed_checkpoint.sha256)

    def _append_suffix(self, store: HLT16CampaignStore, kind: str) -> dict[str, object]:
        key = MEMBER_KEYS[0]
        member = self.seed_checkpoint.members[key]
        cursor = p15.Proto15Cursor(dict(member.cursor))
        if kind in {"source_rejection", "cfl_rejection"}:
            payload: dict[str, object] = {
                "member_key": key, "evidence": {"owner": "test", "reason": kind},
                "predecessor_descriptor_sha256": member.descriptor_sha256,
                "predecessor_cursor_sha256": cursor.sha256,
                "retry_current": 1, "retry_total": 1, "retry_limit": 32,
                "attempted_cap_hex": float(1 / 16).hex(),
                "next_cap_hex": float(1 / 32).hex(), "exhausted": False, "no_commit": True,
            }
        elif kind == "tdg6_rejection":
            persisted = store.load_state(member.descriptor_sha256)
            physical = array_content_sha256(
                persisted.arrays["u"], persisted.arrays["p"], persisted.arrays["q"]
            )
            ledger = p15._ledger_from_mapping(member.ledger)
            admission = p15._synthetic_failed_admission(
                p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]]
            )
            evidence = p15.TDG6TemporalRetryEvidence(
                event_type="rejected_TDG6_temporal_admission",
                method=p15.RUNTIME_METHOD_BY_LABEL[cursor.payload["method"]],
                initial_time=23 / 16,
                attempted_macro_step_size=1 / 16,
                retry_step_size=1 / 32,
                initial_state_sha256=physical,
                previous_step_index=cursor.payload["previous_step_index"],
                previous_transaction_serial=cursor.payload[
                    "previous_transaction_serial"
                ],
                retry_count_for_current_macro_step=1,
                cumulative_temporal_retry_count=(
                    ledger.cumulative_temporal_retry_count + 1
                ),
                failed_channels=admission.failed_channels,
                channel_admissions=admission.channel_admissions,
            )
            payload = {
                "member_key": key, "predecessor_descriptor_sha256": member.descriptor_sha256,
                "predecessor_cursor_sha256": cursor.sha256,
                "evidence": p15._json_safe(evidence.as_mapping()),
            }
        elif kind == "terminal_lock":
            payload = {
                "member_key": key, "reason": "test-terminal", "classification": "invalid",
                "descriptor_sha256": member.descriptor_sha256, "cursor_sha256": cursor.sha256,
                "owner": "invalid", "rejection_record_sha256": None, "evidence": {"owner": "test"},
            }
        elif kind == "common_event_commit":
            payload = {
                "member_keys": MEMBER_KEYS, "predecessor_event": 23, "completed_event": 24,
                "completed_target": {"rational": "3/2", "binary64_hex": (3 / 2).hex()},
                "next_target": {"rational": "25/16", "binary64_hex": (25 / 16).hex()},
                "predecessor_checkpoint_sha256": self.seed_checkpoint.sha256,
            }
        else:
            raise AssertionError(kind)
        record = journal_envelope(
            sequence=self.seed_checkpoint.journal_sequence + 1,
            campaign_id=self.seed_checkpoint.campaign_id,
            generation=self.seed_checkpoint.generation + 1,
            event=self.seed_checkpoint.event,
            previous_record_sha256=self.seed_checkpoint.journal_tip_sha256,
            kind=kind, payload=payload,
        )
        token = store.acquire_active_writer(
            self.seed_checkpoint, host=os.uname().nodename, pid=os.getpid(),
        )
        try:
            store.publish_journal(
                record,
                capability=store.active_writer_capability(self.seed_checkpoint),
            )
        finally:
            store.release_active_writer(
                owner_token=token, host=os.uname().nodename, pid=os.getpid(),
            )
        return record

    @staticmethod
    def _link_after_publication_stage(target: Path) -> Path:
        raw = target.read_bytes()
        stage = target.parent / f".{target.name}.hlt16-stage-{hashlib.sha256(raw).hexdigest()}"
        os.link(target, stage)
        return stage

    @staticmethod
    def _tree_bytes(root: Path) -> tuple[tuple[str, bytes], ...]:
        """Capture all leaves so an inspection test can prove no mutation."""
        leaves, _directories = _walk_tree(root)
        return tuple((relative, (root / relative).read_bytes()) for relative in leaves)

    def test_inspect_recovery_is_read_only_for_after_link_and_partial_stages(self) -> None:
        """A certificate/status read may observe a suffix but never repair it."""
        temporary, store = self.adopted()
        with temporary:
            target = next((store.root / "journal").glob("*.journal"))
            stage = self._link_after_publication_stage(target)
            before = self._tree_bytes(store.root)
            status = store.inspect_recovery()
            self.assertEqual(status.state, "recognized_staging_suffix")
            self.assertFalse(status.safe_to_restart)
            self.assertTrue(stage.exists())
            self.assertEqual(self._tree_bytes(store.root), before)

        temporary, store = self.adopted()
        with temporary:
            raw = b"partial-stage-must-not-be-repaired"
            stage = store.root / "payloads" / (
                ".partial.npz.hlt16-stage-" + hashlib.sha256(raw).hexdigest()
            )
            stage.write_bytes(raw)
            before = self._tree_bytes(store.root)
            with self.assertRaises(HLT16CampaignStoreError):
                store.inspect_recovery()
            self.assertTrue(stage.exists())
            self.assertEqual(self._tree_bytes(store.root), before)

    def test_inspect_recovery_rejects_mismatched_writer_authority_without_mutation(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            body = {
                "schema": "FGC-1-HLT16-active-writer-v1",
                "authorization_commit": "f" * 40,
                "plan_sha256": self.seed_checkpoint.plan_sha256,
                "checkpoint_sha256": self.seed_checkpoint.sha256,
                "host": os.uname().nodename,
                "pid": 99999999,
                "owner_token": "e" * 64,
            }
            atomic_publish_noreplace(
                store.root, "locks/active-write.lock", canonical(body)
            )
            before = self._tree_bytes(store.root)
            with self.assertRaisesRegex(
                HLT16CampaignStoreError, "active writer lock authority differs"
            ):
                store.inspect_recovery()
            self.assertEqual(self._tree_bytes(store.root), before)

    def test_read_only_snapshot_reports_each_recognized_journal_suffix_without_recovery_decision(self) -> None:
        for kind in ("source_rejection", "cfl_rejection", "tdg6_rejection", "terminal_lock", "common_event_commit"):
            with self.subTest(kind=kind):
                temporary, store = self.adopted()
                with temporary:
                    record = self._append_suffix(store, kind)
                    snapshot = store.authenticated_snapshot()
                    self.assertEqual(snapshot.suffix_classification, kind)
                    self.assertEqual(snapshot.suffix_kinds, (kind,))
                    self.assertEqual(snapshot.suffix_records[0]["record_sha256"], record["record_sha256"])
                    self.assertFalse(store.recover().safe_to_restart)

    def test_read_only_snapshot_reports_payload_orphan_and_does_not_repair_staging(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            arrays = {name: value.copy(order="C") for name, value in self.snapshots[MEMBER_KEYS[0]].arrays.items()}
            arrays["u"][0, 0] += 0.125
            raw, semantic, _manifest = encode_payload(arrays)
            atomic_publish_noreplace(store.root, f"payloads/{semantic}.npz", raw)
            snapshot = store.authenticated_snapshot()
            self.assertEqual(snapshot.suffix_classification, "payload_orphan")
            self.assertEqual(snapshot.orphan_payload_semantic_sha256, semantic)
            self.assertFalse(store.recover().safe_to_restart)

        temporary, store = self.adopted()
        with temporary:
            target = "a" * 64 + ".npz"
            raw = b"read-only-stage"
            stage = store.root / "payloads" / f".{target}.hlt16-stage-{hashlib.sha256(raw).hexdigest()}"
            stage.write_bytes(raw)
            snapshot = store.authenticated_snapshot()
            self.assertEqual(snapshot.suffix_classification, "staging_suffix")
            self.assertEqual(snapshot.staging_paths, (f"payloads/{stage.name}",))
            self.assertTrue(stage.exists(), "read-only snapshot must not repair a publication stage")

    def test_snapshot_inspects_proven_after_link_stage_without_relaxing_normal_reader(self) -> None:
        """A two-link stage is readable only through the snapshot proof path."""
        for target_kind in (
            "payload", "descriptor", "journal", "checkpoint", "active_lock", "terminal_lock"
        ):
            with self.subTest(target_kind=target_kind):
                temporary, store = self.adopted()
                with temporary:
                    if target_kind == "payload":
                        arrays = {
                            name: value.copy(order="C")
                            for name, value in self.snapshots[MEMBER_KEYS[0]].arrays.items()
                        }
                        arrays["u"][0, 0] += 0.125
                        raw, semantic, _manifest = encode_payload(arrays)
                        target = store.root / "payloads" / f"{semantic}.npz"
                        target.write_bytes(raw)
                    elif target_kind == "journal":
                        target = next((store.root / "journal").glob("*.journal"))
                    elif target_kind == "descriptor":
                        target = store.root / "states" / (
                            self.seed_checkpoint.members[MEMBER_KEYS[0]].descriptor_sha256 + ".json"
                        )
                    elif target_kind == "checkpoint":
                        target = next((store.root / "checkpoints").glob("00000000000000000001-*.json"))
                    else:
                        if target_kind == "active_lock":
                            token = store.acquire_active_writer(
                                self.seed_checkpoint,
                                host=os.uname().nodename,
                                pid=os.getpid(),
                                owner_token="a" * 64,
                            )
                            self.assertEqual(token, "a" * 64)
                            target = store.root / "locks" / "active-write.lock"
                        else:
                            body = {
                                "schema": "FGC-1-HLT16-terminal-lock-v1",
                                "checkpoint_sha256": self.seed_checkpoint.sha256,
                                "journal_tip_sha256": self.seed_checkpoint.journal_tip_sha256,
                            }
                            target = store.root / "locks" / "terminal.lock"
                            target.write_bytes(canonical({
                                **body,
                                "lock_sha256": hashlib.sha256(canonical(body)).hexdigest(),
                            }))
                    stage = self._link_after_publication_stage(target)
                    relative = str(target.relative_to(store.root))
                    with self.assertRaises(HLT16StateStoreError):
                        read_nofollow(store.root, relative, "ordinary reader")
                    snapshot = store.authenticated_snapshot()
                    self.assertEqual(snapshot.suffix_classification, "staging_suffix")
                    self.assertEqual(snapshot.staging_paths, (str(stage.relative_to(store.root)),))
                    self.assertTrue(stage.exists(), "snapshot must never clean a stage")
                    self.assertTrue(target.exists(), "snapshot must never remove the final object")

    def test_linked_stage_reader_rejects_forged_pair_and_symlink_races(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            target = next((store.root / "journal").glob("*.journal"))
            raw = target.read_bytes()
            forged = target.parent / f".{target.name}.hlt16-stage-{'0' * 64}"
            os.link(target, forged)
            with self.assertRaises(HLT16CampaignStoreError):
                store.authenticated_snapshot()

        temporary, store = self.adopted()
        with temporary:
            target = next((store.root / "journal").glob("*.journal"))
            stage = self._link_after_publication_stage(target)
            relative = str(target.relative_to(store.root))
            validated = _validated_stage(store.root, str(stage.relative_to(store.root)))
            stage.unlink()
            os.symlink(target.name, stage)
            with self.assertRaises(HLT16CampaignStoreError):
                _safe_read_snapshot(store.root, relative, "forged stage reader", {relative: validated})

    def test_symlink_foreign_and_partial_trees_fail_closed(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            os.symlink(store.root / "journal", store.root / "locks" / "bad")
            with self.assertRaises(HLT16CampaignStoreError):
                store.recover()
        temporary, store = self.adopted()
        with temporary:
            (store.root / "foreign-root").mkdir()
            with self.assertRaises(HLT16CampaignStoreError):
                store.recover()
        temporary, store = self.adopted()
        with temporary:
            (store.root / "states" / ("f" * 64 + ".json")).write_bytes(b"{}")
            with self.assertRaises(HLT16CampaignStoreError):
                store.recover()

    def test_nonpublication_hardlink_and_checkpoint_fork_cannot_select_a_convenient_history(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            checkpoint = next((store.root / "checkpoints").glob("00000000000000000001-*.json"))
            os.link(checkpoint, store.root / "checkpoints" / "foreign-copy.json")
            with self.assertRaises(HLT16CampaignStoreError):
                store.recover()
        temporary, store = self.adopted()
        with temporary:
            checkpoint = next((store.root / "checkpoints").glob("00000000000000000001-*.json"))
            fork = store.root / "checkpoints" / ("00000000000000000002-" + "a" * 64 + ".json")
            fork.write_bytes(checkpoint.read_bytes())
            with self.assertRaises(HLT16CampaignStoreError):
                store.recover()

    def test_one_link_stage_is_recognized_but_never_restartable_or_runnable(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            target = "a" * 64 + ".npz"
            raw = b"interrupted-payload"
            stage = store.root / "payloads" / f".{target}.hlt16-stage-{hashlib.sha256(raw).hexdigest()}"
            stage.write_bytes(raw)
            status = store.recover()
            self.assertEqual(status.state, "recognized_staging_suffix")
            self.assertFalse(status.safe_to_restart)
            with self.assertRaises(HLT16CampaignStoreError):
                store.acquire_active_writer(self.seed_checkpoint, host=os.uname().nodename)

    def test_writer_lease_is_token_and_host_bound(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            token = store.acquire_active_writer(self.seed_checkpoint, host=os.uname().nodename, pid=os.getpid(), owner_token="a" * 64)
            self.assertEqual(token, "a" * 64)
            self.assertEqual(store.recover().state, "live_local_writer")
            with self.assertRaises(HLT16CampaignStoreError):
                store.release_active_writer(owner_token="b" * 64, host=os.uname().nodename, pid=os.getpid())
            store.release_active_writer(owner_token=token, host=os.uname().nodename, pid=os.getpid())
            self.assertEqual(store.recover().state, "clean_checkpoint")

    def test_interrupted_writer_acquisition_replays_complete_or_discards_partial_stage(self) -> None:
        host = os.uname().nodename
        body = {
            "schema": "FGC-1-HLT16-active-writer-v1",
            "authorization_commit": self.seed_checkpoint.authorization_commit,
            "plan_sha256": self.seed_checkpoint.plan_sha256,
            "checkpoint_sha256": self.seed_checkpoint.sha256,
            "host": host,
            "pid": 99999999,
            "owner_token": "a" * 64,
        }
        raw = canonical(body)
        stage_name = (
            ".active-write.lock.hlt16-stage-" + hashlib.sha256(raw).hexdigest()
        )

        temporary, store = self.adopted()
        with temporary:
            stage = store.root / "locks" / stage_name
            stage.write_bytes(raw)
            self.assertEqual(
                store.recover_interrupted_writer_acquisition(
                    self.seed_checkpoint, host=host,
                ),
                "replayed_complete_writer_stage",
            )
            self.assertFalse(stage.exists())
            self.assertEqual(store.recover().writer_state, "stale_verified_lock")
            token = store.take_over_stale_writer(
                self.seed_checkpoint,
                host=host,
                pid=os.getpid(),
                owner_token="b" * 64,
            )
            store.release_active_writer(
                owner_token=token, host=host, pid=os.getpid(),
            )

        temporary, store = self.adopted()
        with temporary:
            stage = store.root / "locks" / stage_name
            stage.write_bytes(raw[: len(raw) // 2])
            self.assertEqual(
                store.recover_interrupted_writer_acquisition(
                    self.seed_checkpoint, host=host,
                ),
                "discarded_partial_writer_stage",
            )
            self.assertFalse(stage.exists())
            self.assertEqual(store.recover().state, "clean_checkpoint")

    def test_writer_acquisition_recovery_refuses_an_unrelated_stage(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            raw = b"unrelated"
            stage = store.root / "journal" / (
                ".00000000000000000001-" + "a" * 64
                + ".journal.hlt16-stage-" + hashlib.sha256(raw).hexdigest()
            )
            stage.write_bytes(raw)
            with self.assertRaisesRegex(
                HLT16CampaignStoreError, "non-writer publication stage",
            ):
                store.recover_interrupted_writer_acquisition(
                    self.seed_checkpoint, host=os.uname().nodename,
                )

    def test_stale_takeover_requires_same_host_and_exact_authority(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            store.acquire_active_writer(self.seed_checkpoint, host=os.uname().nodename, pid=99999999, owner_token="a" * 64)
            self.assertEqual(store.recover().state, "stale_verified_lock")
            with self.assertRaises(HLT16CampaignStoreError):
                store.take_over_stale_writer(self.seed_checkpoint, host="foreign.example", owner_token="b" * 64)
            token = store.take_over_stale_writer(self.seed_checkpoint, host=os.uname().nodename, pid=os.getpid(), owner_token="b" * 64)
            self.assertEqual(token, "b" * 64)
            store.release_active_writer(owner_token=token, host=os.uname().nodename, pid=os.getpid())

    def test_foreign_host_lease_is_not_declared_stale(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            body = {"schema": "FGC-1-HLT16-active-writer-v1", "authorization_commit": self.seed_checkpoint.authorization_commit,
                    "plan_sha256": self.seed_checkpoint.plan_sha256, "checkpoint_sha256": self.seed_checkpoint.sha256,
                    "host": "foreign.example", "pid": 1, "owner_token": "a" * 64}
            atomic_publish_noreplace(store.root, "locks/active-write.lock", canonical(body))
            status = store.recover()
            self.assertEqual(status.state, "foreign_host_writer")
            self.assertFalse(status.safe_to_restart)
            with self.assertRaises(HLT16CampaignStoreError):
                store.take_over_stale_writer(self.seed_checkpoint, host=os.uname().nodename)

    def test_exact_rejection_journal_replay_is_idempotent_but_a_same_slot_variant_is_a_fork(self) -> None:
        temporary, store = self.adopted()
        with temporary:
            token = store.acquire_active_writer(
                self.seed_checkpoint, host=os.uname().nodename, pid=os.getpid(), owner_token="f" * 64,
            )
            key = MEMBER_KEYS[0]
            member = self.seed_checkpoint.members[key]
            cursor = p15.Proto15Cursor(dict(member.cursor))
            payload = {
                "member_key": key,
                "evidence": {"owner": "test", "reason": "synthetic-source-retry"},
                "predecessor_descriptor_sha256": member.descriptor_sha256,
                "predecessor_cursor_sha256": cursor.sha256,
                "retry_current": 1,
                "retry_total": 1,
                "retry_limit": 32,
                "attempted_cap_hex": float(1 / 16).hex(),
                "next_cap_hex": float(1 / 32).hex(),
                "exhausted": False,
                "no_commit": True,
            }
            record = journal_envelope(
                sequence=self.seed_checkpoint.journal_sequence + 1,
                campaign_id=self.seed_checkpoint.campaign_id,
                generation=self.seed_checkpoint.generation + 1,
                event=self.seed_checkpoint.event,
                previous_record_sha256=self.seed_checkpoint.journal_tip_sha256,
                kind="source_rejection", payload=payload,
            )
            capability = store.active_writer_capability(self.seed_checkpoint)
            self.assertEqual(store.publish_journal(record, capability=capability), record["record_sha256"])
            self.assertEqual(store.publish_journal(record, capability=capability), record["record_sha256"])
            changed = journal_envelope(
                sequence=record["sequence"], campaign_id=record["campaign_id"],
                generation=record["generation"], event=record["event"],
                previous_record_sha256=record["previous_record_sha256"], kind="source_rejection",
                payload={**payload, "evidence": {"owner": "test", "reason": "different"}},
            )
            with self.assertRaises(HLT16CampaignStoreError):
                store.publish_journal(changed, capability=capability)
            store.release_active_writer(owner_token=token, host=os.uname().nodename, pid=os.getpid())

    def test_publication_crash_cuts_have_one_fail_closed_or_replay_disposition(self) -> None:
        # Every low-level publication cut is exercised through each campaign
        # owner.  A pre-link one-link stage has no independently provable
        # owner and is retained fail-closed.  Once the final hard link exists,
        # repeating the exact bridge may repair only that proved stage/final
        # pair; it may never choose an alternative genesis.
        cuts = (
            "before_write",
            "after_write",
            "after_file_fsync",
            "after_link",
            "after_parent_fsync",
            "after_stage_cleanup",
        )
        phases = tuple(
            f"{owner}:{cut}"
            for owner in ("state", "journal", "checkpoint")
            for cut in cuts
        )
        for phase in phases:
            with self.subTest(phase=phase), TemporaryDirectory() as temporary:
                root = Path(temporary) / "calibration"
                copy_sealed_gen0_store(ROOT, root)

                def fault(observed: str) -> None:
                    if observed == phase:
                        raise RuntimeError(phase)

                interrupted = HLT16CampaignStore(root, fault_hook=fault)
                with self.assertRaises(RuntimeError):
                    self._adopt(interrupted)
                before = _walk_tree(root)[0]
                cut = phase.split(":", 1)[1]
                if cut in {"before_write", "after_write", "after_file_fsync"}:
                    recovered = HLT16CampaignStore(root)
                    with self.assertRaisesRegex(
                        HLT16CampaignStoreError,
                        "unlinked bootstrap stage ownership is unprovable",
                    ):
                        self._adopt(recovered)
                    # Fail-closed inspection may create no new bridge object
                    # and may not remove the ambiguous inode.
                    self.assertEqual(_walk_tree(root)[0], before)
                    continue

                # Once link(stage, final) has occurred, the inode pair and
                # expected bytes prove the only replayable publication.
                try:
                    status = interrupted.recover()
                except HLT16CampaignStoreError:
                    status = None
                if status is not None:
                    self.assertFalse(status.terminal)
                recovered = HLT16CampaignStore(root)
                checkpoint = self._adopt(recovered)
                self.assertEqual(checkpoint.sha256, self.seed_checkpoint.sha256)
                self.assertEqual(recovered.recover().state, "clean_checkpoint")

    def test_same_named_foreign_unlinked_bootstrap_stage_is_retained_and_rejected(self) -> None:
        for variant in ("partial", "exact"):
            with self.subTest(variant=variant), TemporaryDirectory() as temporary:
                root = Path(temporary) / "calibration"
                copy_sealed_gen0_store(ROOT, root)
                store = HLT16CampaignStore(root)
                expected = store._expected_bootstrap_bridge(
                    plan_sha256="1" * 64,
                    authorization_commit="2" * 40,
                    campaign_id=self.campaign_id,
                    snapshots=self.snapshots,
                    member_factory=self._member_factory,
                    target={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
                    event=23,
                )
                target_relative, target_raw = expected.publications[0]
                target = Path(target_relative)
                digest = hashlib.sha256(target_raw).hexdigest()
                stage = root / target.parent / f".{target.name}.hlt16-stage-{digest}"
                stage.parent.mkdir(exist_ok=True)
                foreign_raw = (
                    b"foreign-partial-stage" if variant == "partial" else target_raw
                )
                stage.write_bytes(foreign_raw)
                with store._bootstrap_guard():
                    pass
                before = _walk_tree(root)[0]
                with self.assertRaisesRegex(
                    HLT16CampaignStoreError,
                    "unlinked bootstrap stage ownership is unprovable",
                ):
                    self._adopt(store)
                self.assertEqual(_walk_tree(root)[0], before)
                self.assertEqual(stage.read_bytes(), foreign_raw)

    def test_bootstrap_precompute_deduplicates_shared_payload_paths_only_when_bytes_match(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary) / "calibration"
            copy_sealed_gen0_store(ROOT, root)
            store = HLT16CampaignStore(root)
            expected = store._expected_bootstrap_bridge(
                plan_sha256="1" * 64,
                authorization_commit="2" * 40,
                campaign_id=self.campaign_id,
                snapshots=self.snapshots,
                member_factory=self._member_factory,
                target={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
                event=23,
            )
            paths = [path for path, _raw in expected.publications]
            self.assertEqual(len(paths), len(set(paths)))
            raw_by_path: dict[str, set[bytes]] = {}
            for snapshot in self.snapshots.values():
                raw, semantic, _manifest = encode_payload(snapshot.arrays)
                raw_by_path.setdefault(f"payloads/{semantic}.npz", set()).add(raw)
            # The current six GEN0 payloads need not collide, but if a future
            # method pair shares a semantic payload this map proves the only
            # permitted duplicate relation: byte-identical content.
            self.assertTrue(all(len(raws) == 1 for raws in raw_by_path.values()))
            self.assertEqual(
                {path for path in paths if path.startswith("payloads/")},
                set(raw_by_path),
            )

    def test_foreign_valid_shaped_bootstrap_prefixes_fail_before_new_bridge_leaf(self) -> None:
        for kind in ("payloads", "states", "journal"):
            with self.subTest(kind=kind), TemporaryDirectory() as temporary:
                root = Path(temporary) / "calibration"
                copy_sealed_gen0_store(ROOT, root)
                store = HLT16CampaignStore(root)
                expected = store._expected_bootstrap_bridge(
                    plan_sha256="1" * 64,
                    authorization_commit="2" * 40,
                    campaign_id=self.campaign_id,
                    snapshots=self.snapshots,
                    member_factory=self._member_factory,
                    target={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
                    event=23,
                )
                source_path, source_raw = next(
                    (path, raw) for path, raw in expected.publications
                    if path.startswith(f"{kind}/")
                )
                if kind == "payloads":
                    foreign = f"payloads/{'f' * 64}.npz"
                elif kind == "states":
                    foreign = f"states/{'f' * 64}.json"
                else:
                    foreign = f"journal/{0:020d}-{'f' * 64}.journal"
                self.assertNotEqual(source_path, foreign)
                target = root / foreign
                target.parent.mkdir(exist_ok=True)
                target.write_bytes(source_raw)
                with store._bootstrap_guard():
                    pass
                before = _walk_tree(root)[0]
                with self.assertRaisesRegex(HLT16CampaignStoreError, "foreign bridge leaf"):
                    self._adopt(store)
                self.assertEqual(_walk_tree(root)[0], before)

    def _retired_path(self, store, body, *, digest=None):
        raw = canonical(body)
        observed = hashlib.sha256(raw).hexdigest() if digest is None else digest
        path = store.root / "locks" / f".active-write.lock.hlt16-quarantine-{observed}"
        path.write_bytes(raw)
        return path

    def _retired_body(self):
        return {
            "schema": "FGC-1-HLT16-active-writer-v1",
            "authorization_commit": self.seed_checkpoint.authorization_commit,
            "plan_sha256": self.seed_checkpoint.plan_sha256,
            "checkpoint_sha256": self.seed_checkpoint.sha256,
            "host": os.uname().nodename,
            "pid": os.getpid(),
            "owner_token": "9" * 64,
        }

    def test_retired_lease_rejects_malformed_body(self):
        temporary, store = self.adopted()
        with temporary:
            path = store.root / "locks" / (".active-write.lock.hlt16-quarantine-" + "a" * 64)
            path.write_bytes(b"not-canonical-json")
            with self.assertRaises(HLT16CampaignStoreError):
                store.recover()

    def test_retired_lease_rejects_filename_body_digest_mismatch(self):
        temporary, store = self.adopted()
        with temporary:
            self._retired_path(store, self._retired_body(), digest="a" * 64)
            with self.assertRaisesRegex(HLT16CampaignStoreError, "filename hash"):
                store.recover()

    def test_retired_lease_rejects_symlink_or_hardlink_leaf(self):
        temporary, store = self.adopted()
        with temporary:
            path = self._retired_path(store, self._retired_body())
            alias = path.with_name(path.name + ".alias")
            os.link(path, alias)
            with self.assertRaises(HLT16CampaignStoreError):
                store.recover()

    def test_retired_lease_rejects_foreign_authority_or_checkpoint(self):
        temporary, store = self.adopted()
        with temporary:
            body = self._retired_body()
            body["checkpoint_sha256"] = "a" * 64
            self._retired_path(store, body)
            with self.assertRaisesRegex(HLT16CampaignStoreError, "authority"):
                store.recover()


class HLT16BootstrapGuardTests(unittest.TestCase):
    """The first durable bridge is serialized before it has a checkpoint lease."""

    def test_second_store_cannot_enter_bootstrap_guard_until_first_releases(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary) / "calibration"
            (root / "locks").mkdir(parents=True)
            first = HLT16CampaignStore(root)
            second = HLT16CampaignStore(root)
            entered = threading.Event()
            release = threading.Event()
            second_entered = threading.Event()
            errors: list[BaseException] = []

            def hold_first() -> None:
                try:
                    with first._bootstrap_guard():
                        entered.set()
                        if not release.wait(2.0):
                            raise AssertionError("bootstrap release timed out")
                except BaseException as exc:  # surface thread failures below
                    errors.append(exc)

            def wait_second() -> None:
                try:
                    with second._bootstrap_guard():
                        second_entered.set()
                except BaseException as exc:
                    errors.append(exc)

            owner = threading.Thread(target=hold_first)
            contender = threading.Thread(target=wait_second)
            owner.start()
            self.assertTrue(entered.wait(1.0))
            contender.start()
            time.sleep(0.1)
            self.assertFalse(second_entered.is_set())
            release.set()
            owner.join(2.0)
            contender.join(2.0)
            self.assertFalse(owner.is_alive())
            self.assertFalse(contender.is_alive())
            self.assertEqual(errors, [])
            self.assertTrue(second_entered.is_set())

    def test_interrupted_bootstrap_guard_leaves_no_durable_bridge_object(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary) / "calibration"
            (root / "locks").mkdir(parents=True)
            store = HLT16CampaignStore(root)
            with self.assertRaisesRegex(RuntimeError, "interrupted bootstrap"):
                with store._bootstrap_guard():
                    raise RuntimeError("interrupted bootstrap")
            # The advisory inode can remain, but no payload, state, journal,
            # or checkpoint has been written before exclusive ownership was
            # released.  A later launcher therefore starts from the same
            # immutable GEN0 boundary rather than a half-bridge.
            self.assertEqual(sorted(path.name for path in (root / "locks").iterdir()), ["bootstrap.guard"])
            with store._bootstrap_guard():
                pass


if __name__ == "__main__":
    unittest.main()
