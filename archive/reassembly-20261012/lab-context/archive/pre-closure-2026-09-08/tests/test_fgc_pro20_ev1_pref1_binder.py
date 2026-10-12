"""Outcome-neutral hostile-store controls for the PRO20 PREF1 binder core."""

from __future__ import annotations

import ast
from dataclasses import replace
from hashlib import sha256
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import (  # noqa: E402
    DEFAULT_MAX_BYTES,
    canonical_json_bytes,
    load_canonical_json,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    encode_hlt17_cursor,
    restore_hlt17_cursor,
)
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (  # noqa: E402
    HLT17GenerationBundle,
)
from recursive_horizons.fgc.evolution.pro20_ev1_pref1_binder import (  # noqa: E402
    CONTAINER_STAGING_PREFIX,
    PRO20EV1PREF1Error,
    PRODUCTION_NAMESPACE,
    _MAX_LEAF_BYTES,
    _bundle_state,
    bind_independent_seed,
    bind_terminal_store,
    independent_outcome,
    independent_progression,
)
from recursive_horizons.fgc.evolution.pro20_ev1_store import (  # noqa: E402
    PRO20EV1CampaignStore,
)
from tests.test_fgc_pro20_ev1_store import _live_cohort, _receipt  # noqa: E402


def _terminal_repository(root: Path) -> Path:
    (root / "runs/fgc-2-sf1").mkdir(parents=True)
    _checkpoints, bundles = _live_cohort()
    store = PRO20EV1CampaignStore.publish_seed_store(
        root, receipt=_receipt(), bundles=bundles
    )
    writer = store.acquire_writer(owner_token="pref1-test", host="test-host")
    view = store.authenticate()
    first = view.generations[-1]
    member_key = "RK4-2049"
    store.publish_attempt(
        writer,
        kind="resource_exhausted",
        member_key=member_key,
        successor_bundle=first.bundles[member_key],
        terminal_evidence={
            "kind": "resource_exhausted",
            "reason": "synthetic_resource_ceiling",
            "executable_retry_cursor": False,
            "spends_retry": False,
        },
    )
    store.close_writer(writer)
    return root


def _transition_terminal_repository(root: Path, transition: str) -> Path:
    (root / "runs/fgc-2-sf1").mkdir(parents=True)
    checkpoints, bundles = _live_cohort(controllable=transition == "source_retry")
    store = PRO20EV1CampaignStore.publish_seed_store(
        root, receipt=_receipt(), bundles=bundles
    )
    writer = store.acquire_writer(owner_token="pref1-transition", host="test-host")
    checkpoint = checkpoints["RK4-2049"]
    if transition == "source_retry":
        checkpoint.member.operator.mode = "source"
    result = checkpoint.advance_once(requested_cap=1.0 / 32.0)
    expected = "source_retry_required" if transition == "source_retry" else "accepted_fine"
    if result.disposition != expected:
        raise AssertionError(
            f"synthetic transition differs: {result.disposition!r} != {expected!r}"
        )
    store.publish_attempt(
        writer,
        kind="source_retry" if transition == "source_retry" else "accepted_fine",
        member_key="RK4-2049",
        successor_bundle=result.bundle,
    )
    last = store.authenticate().generations[-1]
    store.publish_attempt(
        writer,
        kind="resource_exhausted",
        member_key="RK4-2049",
        successor_bundle=last.bundles["RK4-2049"],
        terminal_evidence={
            "kind": "resource_exhausted",
            "reason": "synthetic_after_transition",
            "executable_retry_cursor": False,
            "spends_retry": False,
        },
    )
    store.close_writer(writer)
    return root


def _store(root: Path) -> Path:
    return root.joinpath(*PRODUCTION_NAMESPACE.split("/"))


def _canonical_write(path: Path, value: object) -> None:
    path.write_bytes(canonical_json_bytes(value))


def _relabel_terminal(root: Path) -> None:
    """Make every stored label self-consistent while retaining resource evidence."""

    store = _store(root)
    generation = store / "generation-0000000000000001"
    journal_path = generation / "journal.json"
    checkpoint_path = generation / "checkpoint.json"
    manifest_path = generation / "generation_manifest.json"
    lock_path = store / "terminal.lock"

    journal = load_canonical_json(journal_path.read_bytes())
    journal["kind"] = "invalid_premise"
    _canonical_write(journal_path, journal)
    journal_hash = sha256(journal_path.read_bytes()).hexdigest()

    checkpoint = load_canonical_json(checkpoint_path.read_bytes())
    checkpoint["journal_record_sha256"] = journal_hash
    checkpoint["kind"] = "invalid_premise"
    checkpoint["disposition"] = "invalid_premise"
    _canonical_write(checkpoint_path, checkpoint)
    checkpoint_hash = sha256(checkpoint_path.read_bytes()).hexdigest()

    manifest = load_canonical_json(manifest_path.read_bytes())
    manifest["journal_record_sha256"] = journal_hash
    manifest["checkpoint_sha256"] = checkpoint_hash
    manifest["kind"] = "invalid_premise"
    manifest["leaves"]["journal.json"] = journal_hash
    manifest["leaves"]["checkpoint.json"] = checkpoint_hash
    _canonical_write(manifest_path, manifest)

    lock = load_canonical_json(lock_path.read_bytes())
    lock["checkpoint_sha256"] = checkpoint_hash
    lock["journal_record_sha256"] = journal_hash
    lock["disposition"] = "invalid_premise"
    _canonical_write(lock_path, lock)


def _generation_dir(root: Path, index: int) -> Path:
    return _store(root) / f"generation-{index:016d}"


def _blob_bytes(store: Path, digest: str) -> bytes:
    matches = list(store.glob(f"generation-*/blobs/{digest}"))
    if not matches:
        raise AssertionError(f"missing blob {digest}")
    return matches[0].read_bytes()


def _commit_generation(
    store: Path,
    index: int,
    *,
    journal: dict[str, object] | None = None,
    checkpoint: dict[str, object] | None = None,
) -> tuple[str, str]:
    gen = store / f"generation-{index:016d}"
    if journal is not None:
        _canonical_write(gen / "journal.json", journal)
    if checkpoint is not None:
        _canonical_write(gen / "checkpoint.json", checkpoint)
    journal_hash = sha256((gen / "journal.json").read_bytes()).hexdigest()
    checkpoint_obj = load_canonical_json((gen / "checkpoint.json").read_bytes())
    if checkpoint_obj["journal_record_sha256"] != journal_hash:
        checkpoint_obj["journal_record_sha256"] = journal_hash
        _canonical_write(gen / "checkpoint.json", checkpoint_obj)
    checkpoint_hash = sha256((gen / "checkpoint.json").read_bytes()).hexdigest()
    leaves = {
        path.relative_to(gen).as_posix(): sha256(path.read_bytes()).hexdigest()
        for path in sorted(gen.rglob("*"))
        if path.is_file() and path.name != "generation_manifest.json"
    }
    new_blobs = sorted(
        path.name for path in (gen / "blobs").glob("*") if path.is_file()
    ) if (gen / "blobs").is_dir() else []
    manifest = load_canonical_json((gen / "generation_manifest.json").read_bytes())
    manifest["journal_record_sha256"] = journal_hash
    manifest["checkpoint_sha256"] = checkpoint_hash
    manifest["leaves"] = leaves
    manifest["new_blob_sha256s"] = new_blobs
    _canonical_write(gen / "generation_manifest.json", manifest)
    return journal_hash, checkpoint_hash


def _relink_from(store: Path, start_index: int) -> None:
    gens = sorted(path for path in store.glob("generation-*") if path.is_dir())
    prev_journal = sha256((gens[start_index] / "journal.json").read_bytes()).hexdigest()
    prev_checkpoint = sha256(
        (gens[start_index] / "checkpoint.json").read_bytes()
    ).hexdigest()
    for gen in gens[start_index + 1 :]:
        index = int(gen.name.split("-")[1])
        journal = load_canonical_json((gen / "journal.json").read_bytes())
        checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
        journal["journal_parent_sha256"] = prev_journal
        journal["predecessor_checkpoint_sha256"] = prev_checkpoint
        checkpoint["journal_parent_sha256"] = prev_journal
        checkpoint["predecessor_checkpoint_sha256"] = prev_checkpoint
        prev_journal, prev_checkpoint = _commit_generation(
            store, index, journal=journal, checkpoint=checkpoint
        )
    lock_path = store / "terminal.lock"
    last = gens[-1]
    lock = load_canonical_json(lock_path.read_bytes())
    lock["journal_record_sha256"] = sha256((last / "journal.json").read_bytes()).hexdigest()
    lock["checkpoint_sha256"] = sha256((last / "checkpoint.json").read_bytes()).hexdigest()
    lock["store_generation"] = int(last.name.split("-")[1])
    _canonical_write(lock_path, lock)


def _replace_member_cursor(
    root: Path,
    index: int,
    member_key: str,
    new_cursor: bytes,
) -> None:
    store = _store(root)
    gen = _generation_dir(root, index)
    checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
    journal = load_canonical_json((gen / "journal.json").read_bytes())
    state = checkpoint["members"][member_key]
    bundle = HLT17GenerationBundle(
        _blob_bytes(store, str(state["descriptor_sha256"])),
        _blob_bytes(store, str(state["payload_sha256"])),
        new_cursor,
    )
    derived = dict(_bundle_state(bundle))
    old_digest = str(state["cursor_sha256"])
    new_digest = derived["cursor_sha256"]
    checkpoint["members"][member_key] = derived
    blob_dir = gen / "blobs"
    blob_dir.mkdir(exist_ok=True)
    (blob_dir / new_digest).write_bytes(new_cursor)
    stale = blob_dir / old_digest
    if stale.exists() and old_digest != new_digest:
        stale.unlink()
    if journal["kind"] == "seed":
        journal["payload"]["members"][member_key] = derived
    else:
        payload = journal["payload"]
        payload["successor_cursor_sha256"] = derived["cursor_sha256"]
        payload["successor_descriptor_sha256"] = derived["descriptor_sha256"]
        payload["successor_payload_sha256"] = derived["payload_sha256"]
        for key in (
            "kernel_transition_parent_sha256",
            "kernel_transition_digest",
            "accepted_generation",
            "cursor_generation",
            "attempt_serial",
            "kernel_checkpoint_generation",
            "source_current",
            "source_total",
            "cfl_current",
            "cfl_total",
            "mode",
            "latest_owner",
            "accepted_time_hex",
            "physical_state_sha256",
            "tdg6_snapshot_sha256",
            "parked",
            "public_fresh",
            "step_index",
            "transaction_serial",
        ):
            payload[key] = derived[key]
    _commit_generation(store, index, journal=journal, checkpoint=checkpoint)
    gens = sorted(path for path in store.glob("generation-*") if path.is_dir())
    for gen in gens[index + 1 :]:
        later = int(gen.name.split("-")[1])
        later_checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
        later_journal = load_canonical_json((gen / "journal.json").read_bytes())
        later_state = later_checkpoint["members"][member_key]
        if later_state["cursor_sha256"] != old_digest:
            continue
        later_checkpoint["members"][member_key] = derived
        payload = later_journal["payload"]
        if later_journal["kind"] != "seed" and payload.get("member_key") == member_key:
            if payload.get("predecessor_cursor_sha256") == old_digest:
                payload["predecessor_cursor_sha256"] = new_digest
            if payload.get("successor_cursor_sha256") == old_digest:
                payload["successor_cursor_sha256"] = new_digest
            for key in (
                "kernel_transition_parent_sha256",
                "kernel_transition_digest",
                "source_current",
                "source_total",
                "cfl_current",
                "cfl_total",
                "mode",
                "latest_owner",
                "accepted_time_hex",
                "physical_state_sha256",
                "tdg6_snapshot_sha256",
                "parked",
                "public_fresh",
                "step_index",
                "transaction_serial",
            ):
                payload[key] = derived[key]
        _commit_generation(
            store, later, journal=later_journal, checkpoint=later_checkpoint
        )
    _relink_from(store, index)


def _cursor_from_generation(root: Path, index: int, member_key: str):
    gen = _generation_dir(root, index)
    checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
    digest = str(checkpoint["members"][member_key]["cursor_sha256"])
    raw = _blob_bytes(_store(root), digest)
    return restore_hlt17_cursor(load_canonical_json(raw)), raw


def _write_raw_cursor(
    root: Path,
    index: int,
    member_key: str,
    raw: bytes,
) -> None:
    store = _store(root)
    gen = _generation_dir(root, index)
    checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
    old_digest = str(checkpoint["members"][member_key]["cursor_sha256"])
    new_digest = sha256(raw).hexdigest()
    checkpoint["members"][member_key]["cursor_sha256"] = new_digest
    blob_dir = gen / "blobs"
    blob_dir.mkdir(exist_ok=True)
    (blob_dir / new_digest).write_bytes(raw)
    stale = blob_dir / old_digest
    if stale.exists() and old_digest != new_digest:
        stale.unlink()
    _commit_generation(store, index, checkpoint=checkpoint)
    _relink_from(store, index)


class PRO20EV1PREF1BinderCoreTests(unittest.TestCase):
    def test_completed_typed_terminal_is_bound_without_promotion(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            result = bind_terminal_store(root)
        self.assertEqual(result.terminal_kind, "resource_exhausted")
        self.assertEqual(result.disposition, "resource_exhausted")
        self.assertFalse(result.first_event_complete)
        self.assertFalse(result.runner_label_sufficient)
        self.assertFalse(result.calibration_eligible)
        self.assertEqual(len(result.generations), 2)
        self.assertEqual(result.session_count, 1)
        outcome = independent_outcome(result)
        self.assertEqual(
            outcome["classification"],
            "independently_bound_resource_exhausted_no_common_event",
        )
        self.assertFalse(outcome["common_event_completed"])
        self.assertFalse(outcome["calibration_authorized"])
        progression = independent_progression(result)
        self.assertEqual(progression["accepted_fine_count"], 0)
        self.assertFalse(progression["accepted_state_advanced"])
        self.assertFalse(progression["binder_mutated_store"])

    def test_accepted_fine_is_derived_from_bundle_bytes(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "accepted_fine")
            result = bind_terminal_store(root)
        self.assertEqual(result.generations[1].kind, "accepted_fine")
        self.assertEqual(result.generations[1].derived_kind, "accepted_fine")
        self.assertEqual(result.generations[2].derived_kind, "resource_exhausted")
        progression = independent_progression(result)
        self.assertEqual(progression["accepted_fine_count"], 1)
        self.assertTrue(progression["accepted_state_advanced"])
        self.assertFalse(progression["common_event_completed"])

    def test_source_retry_is_derived_from_cursor_overlay(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "source_retry")
            result = bind_terminal_store(root)
        self.assertEqual(result.generations[1].kind, "source_retry")
        self.assertEqual(result.generations[1].derived_kind, "source_retry")
        self.assertEqual(
            result.generations[1].member_states["RK4-2049"]["latest_owner"],
            "source",
        )

    def test_seed_must_match_a_separately_reconstructed_origin(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            terminal = bind_terminal_store(root)
            _checkpoints, expected = _live_cohort()
            receipt = bind_independent_seed(terminal, expected)
            self.assertTrue(
                receipt["independently_reconstructed_origin_matches_raw_seed"]
            )
            self.assertFalse(receipt["root_manifest_origin_label_sufficient"])
            altered = dict(expected)
            altered["RK4-2049"] = expected["RK4-4097"]
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "raw seed differs"):
                bind_independent_seed(terminal, altered)

    def test_self_consistent_runner_relabel_is_rejected_by_independent_kind(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            _relabel_terminal(root)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "stored kind differs"):
                bind_terminal_store(root)

    def test_terminal_lock_and_closed_writer_are_mandatory(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            (_store(root) / "terminal.lock").unlink()
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "lacks root/lock"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            close = next(_store(root).glob("session-*-close.json"))
            close.unlink()
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "unclosed writer"):
                bind_terminal_store(root)

    def test_terminal_claim_promotion_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            lock_path = _store(root) / "terminal.lock"
            lock = load_canonical_json(lock_path.read_bytes())
            lock["calibration_eligible"] = True
            _canonical_write(lock_path, lock)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "must remain false"):
                bind_terminal_store(root)

    def test_symlink_hardlink_foreign_and_staging_objects_fail_closed(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            lock = _store(root) / "terminal.lock"
            original = lock.read_bytes()
            lock.unlink()
            (lock.parent / "outside").write_bytes(original)
            lock.symlink_to("outside")
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "symlink"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            lock = _store(root) / "terminal.lock"
            os.link(lock, lock.with_name("foreign-hardlink"))
            with self.assertRaises(PRO20EV1PREF1Error):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            (_store(root) / "foreign.json").write_bytes(b"{}")
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "foreign root"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            container = _store(root).parent
            (container / ".calibration.stage-attacker").mkdir()
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "staging sibling"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            (_store(root).parent / "foreign-sibling").mkdir()
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "foreign sibling"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            nested = (
                _store(root)
                / "generation-0000000000000001"
                / "foreign"
                / "nested"
                / "too-deep"
            )
            nested.mkdir(parents=True)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "nesting exceeds"):
                bind_terminal_store(root)

    def test_leaf_byte_bound_matches_neutral_evidence_io(self) -> None:
        self.assertEqual(_MAX_LEAF_BYTES, DEFAULT_MAX_BYTES)

    def test_forward_blob_from_a_later_generation_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "accepted_fine")
            store = _store(root)
            gen1 = store / "generation-0000000000000001"
            checkpoint = load_canonical_json((gen1 / "checkpoint.json").read_bytes())
            digest = str(checkpoint["members"]["RK4-2049"]["cursor_sha256"])
            src = gen1 / "blobs" / digest
            dst_dir = store / "generation-0000000000000002" / "blobs"
            dst_dir.mkdir(exist_ok=True)
            os.replace(src, dst_dir / digest)
            _commit_generation(store, 1)
            _commit_generation(store, 2)
            _relink_from(store, 1)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "future or rollback blob"):
                bind_terminal_store(root)

    def test_rollback_blob_from_a_non_predecessor_generation_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "accepted_fine")
            store = _store(root)
            seed = load_canonical_json(
                (store / "generation-0000000000000000" / "checkpoint.json").read_bytes()
            )
            last = load_canonical_json(
                (store / "generation-0000000000000002" / "checkpoint.json").read_bytes()
            )
            last["members"]["RK4-2049"]["cursor_sha256"] = seed["members"]["RK4-2049"][
                "cursor_sha256"
            ]
            _commit_generation(store, 2, checkpoint=last)
            _relink_from(store, 2)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "future or rollback blob"):
                bind_terminal_store(root)

    def test_reintroduced_and_non_addressed_blobs_are_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "accepted_fine")
            store = _store(root)
            gen1 = store / "generation-0000000000000001"
            digest = next(path.name for path in (gen1 / "blobs").iterdir() if path.is_file())
            payload = (gen1 / "blobs" / digest).read_bytes()
            dst = store / "generation-0000000000000002" / "blobs"
            dst.mkdir(exist_ok=True)
            (dst / digest).write_bytes(payload)
            _commit_generation(store, 2)
            _relink_from(store, 2)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "reintroduced"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            store = _store(root)
            extra = b"non-addressed-blob"
            digest = sha256(extra).hexdigest()
            blob_dir = store / "generation-0000000000000001" / "blobs"
            blob_dir.mkdir(exist_ok=True)
            (blob_dir / digest).write_bytes(extra)
            _commit_generation(store, 1)
            _relink_from(store, 1)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "non-addressed blob"):
                bind_terminal_store(root)

    def test_malformed_fine_and_retry_counters_are_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "accepted_fine")
            cursor, _raw = _cursor_from_generation(root, 1, "RK4-2049")
            mutated = replace(
                cursor,
                source_current=1,
                source_total=max(1, cursor.source_total),
            )
            _replace_member_cursor(
                root, 1, "RK4-2049", canonical_json_bytes(encode_hlt17_cursor(mutated))
            )
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "reset per-macro source/CFL"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "source_retry")
            cursor, _raw = _cursor_from_generation(root, 1, "RK4-2049")
            self.assertIsNotNone(cursor.overlay)
            overlay = replace(
                cursor.overlay,
                owner_current_count=cursor.overlay.owner_current_count + 1,
            )
            mutated = replace(
                cursor,
                overlay=overlay,
                source_current=cursor.source_current + 1,
                source_total=cursor.source_total + 1,
            )
            _replace_member_cursor(
                root, 1, "RK4-2049", canonical_json_bytes(encode_hlt17_cursor(mutated))
            )
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "source current counter"):
                bind_terminal_store(root)

    def test_overlay_next_plan_none_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "source_retry")
            _cursor, raw = _cursor_from_generation(root, 1, "RK4-2049")
            mapping = load_canonical_json(raw)
            mapping["overlay"]["next_plan"] = None
            _write_raw_cursor(root, 1, "RK4-2049", canonical_json_bytes(mapping))
            with self.assertRaisesRegex(
                PRO20EV1PREF1Error, "overlay|next plan|bundle bytes"
            ):
                bind_terminal_store(root)

    def test_temporal_exhaustion_without_ledger_successor_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            store = _store(root)
            gen = store / "generation-0000000000000001"
            journal = load_canonical_json((gen / "journal.json").read_bytes())
            checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
            manifest = load_canonical_json((gen / "generation_manifest.json").read_bytes())
            evidence = {
                "kind": "temporal_retry_exhausted",
                "reason": "coordinate_lattice",
                "executable_retry_cursor": False,
            }
            journal["kind"] = "temporal_retry_exhausted"
            journal["payload"]["terminal_evidence"] = evidence
            journal["payload"]["executable_retry_cursor"] = False
            checkpoint["kind"] = "temporal_retry_exhausted"
            checkpoint["disposition"] = "temporal_retry_exhausted"
            manifest["kind"] = "temporal_retry_exhausted"
            _commit_generation(store, 1, journal=journal, checkpoint=checkpoint)
            _canonical_write(gen / "generation_manifest.json", manifest)
            _commit_generation(store, 1)
            lock = load_canonical_json((store / "terminal.lock").read_bytes())
            lock["disposition"] = "temporal_retry_exhausted"
            _canonical_write(store / "terminal.lock", lock)
            _relink_from(store, 1)
            with self.assertRaisesRegex(
                PRO20EV1PREF1Error, "updated IMP1 ledger|temporal exhaustion"
            ):
                bind_terminal_store(root)

    def test_invalid_premise_spends_retry_omission_and_true_are_rejected(self) -> None:
        for spends in (True, "omit"):
            with TemporaryDirectory() as tmp:
                root = _terminal_repository(Path(tmp).resolve())
                store = _store(root)
                gen = store / "generation-0000000000000001"
                journal = load_canonical_json((gen / "journal.json").read_bytes())
                checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
                evidence = {
                    "kind": "invalid_premise",
                    "reason": "synthetic_invalid_premise",
                    "executable_retry_cursor": False,
                }
                if spends is True:
                    evidence["spends_retry"] = True
                journal["kind"] = "invalid_premise"
                journal["payload"]["terminal_evidence"] = evidence
                journal["payload"]["executable_retry_cursor"] = False
                checkpoint["kind"] = "invalid_premise"
                checkpoint["disposition"] = "invalid_premise"
                _commit_generation(store, 1, journal=journal, checkpoint=checkpoint)
                lock = load_canonical_json((store / "terminal.lock").read_bytes())
                lock["disposition"] = "invalid_premise"
                _canonical_write(store / "terminal.lock", lock)
                _relink_from(store, 1)
                with self.assertRaisesRegex(PRO20EV1PREF1Error, "must not spend a retry"):
                    bind_terminal_store(root)

    def test_parked_member_with_nonzero_retry_current_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _transition_terminal_repository(Path(tmp).resolve(), "accepted_fine")
            cursor, _raw = _cursor_from_generation(root, 1, "RK4-2049")
            mutated = replace(
                cursor,
                event_target=cursor.accepted_time,
                source_current=1,
                source_total=max(1, cursor.source_total),
            )
            _write_raw_cursor(
                root,
                1,
                "RK4-2049",
                canonical_json_bytes(encode_hlt17_cursor(mutated)),
            )
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "nonzero retry current"):
                bind_terminal_store(root)

    def test_real_container_staging_sibling_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            parent = root.joinpath(*"runs/fgc-2-sf1".split("/"))
            (parent / f"{CONTAINER_STAGING_PREFIX}deadbeef").mkdir()
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "staging sibling"):
                bind_terminal_store(root)

    def test_outside_store_hardlink_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            lock = _store(root) / "terminal.lock"
            os.link(lock, root / "outside-store-hardlink")
            with self.assertRaises(PRO20EV1PREF1Error):
                bind_terminal_store(root)

    def test_seed_target_origin_and_genesis_drift_are_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            cursor, _raw = _cursor_from_generation(root, 0, "RK4-2049")
            mutated = replace(cursor, event_target=2.0)
            _replace_member_cursor(
                root, 0, "RK4-2049", canonical_json_bytes(encode_hlt17_cursor(mutated))
            )
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "seed event target"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            cursor, _raw = _cursor_from_generation(root, 0, "RK4-2049")
            mutated = replace(cursor, kernel_transition_parent_sha256="ab" * 32)
            _replace_member_cursor(
                root, 0, "RK4-2049", canonical_json_bytes(encode_hlt17_cursor(mutated))
            )
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "genesis digest"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            store = _store(root)
            gen = store / "generation-0000000000000000"
            checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
            checkpoint["event_origin_time_hex"] = (1.0).hex()
            _commit_generation(store, 0, checkpoint=checkpoint)
            _relink_from(store, 0)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "event bounds|addressing differs"):
                bind_terminal_store(root)

    def test_seam_and_checkpoint_event_drift_are_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            path = _store(root) / "root_manifest.json"
            manifest = load_canonical_json(path.read_bytes())
            manifest["runner_seam"] = "drifted runner seam"
            _canonical_write(path, manifest)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "runner_seam differs"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            store = _store(root)
            gen = store / "generation-0000000000000001"
            checkpoint = load_canonical_json((gen / "checkpoint.json").read_bytes())
            checkpoint["event_target_time_hex"] = (2.0).hex()
            _commit_generation(store, 1, checkpoint=checkpoint)
            _relink_from(store, 1)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "event bounds|addressing differs"):
                bind_terminal_store(root)
        with TemporaryDirectory() as tmp:
            root = _terminal_repository(Path(tmp).resolve())
            path = _store(root) / "root_manifest.json"
            manifest = load_canonical_json(path.read_bytes())
            manifest["campaign_id"] = "../evil"
            _canonical_write(path, manifest)
            with self.assertRaisesRegex(PRO20EV1PREF1Error, "safe identifier"):
                bind_terminal_store(root)

    def test_import_graph_excludes_production_decision_and_writer_modules(self) -> None:
        path = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/pro20_ev1_pref1_binder.py"
        )
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.append(node.module or "")
        for forbidden in (
            "pro20_ev1_runtime",
            "pro20_ev1_protocol",
            "pro20_ev1_store",
            "pro20_ev1_authority",
            "protocol_v19",
        ):
            self.assertFalse(
                any(forbidden in module for module in imported),
                f"binder imports forbidden production owner {forbidden}",
            )


if __name__ == "__main__":
    unittest.main()
