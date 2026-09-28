"""Focused authentication, bootstrap, and no-advance tests for TDG8 successor."""
from __future__ import annotations

from hashlib import sha256
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as campaign_runtime
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore
from recursive_horizons.fgc.evolution.numerical_engine import array_content_sha256
from recursive_horizons.fgc.evolution.proto17_pure_construction import MEMBER_KEYS
from recursive_horizons.fgc.evolution import tdg8_successor_runtime as tdg8


ROOT = Path(__file__).resolve().parents[1]
AUTHORIZATION_COMMIT = "a" * 40


def _tree_metadata(root: Path) -> tuple[tuple[object, ...], ...]:
    rows: list[tuple[object, ...]] = []
    for directory, names, files in os.walk(root, followlinks=False):
        base = Path(directory)
        for name in sorted((*names, *files)):
            path = base / name
            info = path.lstat()
            rows.append(
                (
                    path.relative_to(root).as_posix(),
                    info.st_mode,
                    info.st_ino,
                    info.st_size,
                    info.st_mtime_ns,
                    info.st_ctime_ns,
                )
            )
    return tuple(sorted(rows))


def _tree_bytes(root: Path) -> tuple[tuple[str, str], ...]:
    return tuple(
        sorted(
            (
                path.relative_to(root).as_posix(),
                sha256(path.read_bytes()).hexdigest(),
            )
            for path in root.rglob("*")
            if path.is_file()
        )
    )


def _member_identity(source) -> tuple[tuple[object, ...], ...]:
    return tuple(
        (
            key,
            member.time.hex(),
            member.step_index,
            member.transaction_serial,
            member.CFL_retry_count,
            member.source_retry_count,
            array_content_sha256(member.state.u, member.state.p, member.state.q),
            member.temporal_ledger.accepted_macro_step_count,
            member.temporal_ledger.cumulative_temporal_retry_count,
        )
        for key, member in source.members.items()
    )


class TDG8SuccessorRuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source_root = ROOT / tdg8.SOURCE_RELATIVE
        cls.source_before = _tree_metadata(cls.source_root)
        cls.source = tdg8.authenticate_and_reconstruct(ROOT)
        cls.source_after = _tree_metadata(cls.source_root)
        cls.plan = tdg8.construct_successor_plan(ROOT, AUTHORIZATION_COMMIT)

    def test_terminal_restoration_is_read_only_and_exactly_generation_one(self) -> None:
        self.assertEqual(self.source_before, self.source_after)
        self.assertEqual(tuple(self.source.members), MEMBER_KEYS)
        self.assertEqual(self.source.branch, "GR-0")
        self.assertEqual(self.source.amplitude, "3")
        self.assertEqual(self.source.restart_time.hex(), (23 / 16).hex())
        self.assertEqual(self.source.target_time.hex(), (3 / 2).hex())
        self.assertEqual(
            self.source.evidence.checkpoint["checkpoint_sha256"],
            tdg8.GENERATION_ZERO_CHECKPOINT_SHA256,
        )
        self.assertEqual(
            self.source.evidence.receipt["receipt_sha256"],
            tdg8.GENERATION_ZERO_RECEIPT_SHA256,
        )
        for member in self.source.members.values():
            self.assertEqual(member.time.hex(), (23 / 16).hex())
            self.assertEqual(member.temporal_ledger.last_accepted_time.hex(), (23 / 16).hex())
            self.assertEqual(member.temporal_ledger.accepted_macro_step_count, 0)

    def test_successor_plan_changes_only_campaign_authority_scope(self) -> None:
        self.assertEqual(self.plan.authorization_commit, AUTHORIZATION_COMMIT)
        self.assertEqual(self.plan.protocol_artifact_id, "FGC-2-SF1-PROTO18")
        self.assertEqual(self.plan.campaign_id, tdg8.CAMPAIGN_ID)
        self.assertEqual((self.plan.branch, self.plan.amplitude), ("GR-0", "3"))
        self.assertEqual(self.plan.common_event_index, 23)
        self.assertEqual(self.plan.start_time.rational, "23/16")
        self.assertEqual(self.plan.target_time.rational, "3/2")
        self.assertEqual(tuple(item.member_key for item in self.plan.members), MEMBER_KEYS)

    def test_install_is_the_exact_eight_leaf_prefix(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / tdg8.DESTINATION_RELATIVE
            result = tdg8.install_generation_zero_prefix(
                root, destination, self.source, self.plan.sha256
            )
            self.assertEqual(result["checkpoint"], self.source.evidence.checkpoint)
            self.assertEqual(result["receipt"], self.source.evidence.receipt)

            old_prefix = tdg8.SOURCE_RELATIVE.as_posix() + "/"
            expected = tuple(
                sorted(
                    (
                        str(row["path"]).removeprefix(old_prefix),
                        str(row["raw_sha256"]),
                    )
                    for row in self.source.evidence.leaves
                )
            )
            self.assertEqual(_tree_bytes(destination), expected)
            self.assertEqual(
                {
                    path.relative_to(destination).as_posix()
                    for path in destination.rglob("*")
                    if path.is_dir()
                },
                {"checkpoints", "receipts", "states"},
            )
            self.assertFalse((destination / "journal").exists())
            self.assertFalse((destination / "locks").exists())
            self.assertFalse((destination / "payloads").exists())

    def test_existing_destination_is_foreign_and_never_adopted(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / tdg8.DESTINATION_RELATIVE
            destination.mkdir(parents=True)
            sentinel = destination / "foreign-sentinel"
            sentinel.write_bytes(b"foreign TDG8 destination")
            before = _tree_bytes(destination)
            with self.assertRaisesRegex(
                tdg8.TDG8SuccessorRuntimeError, "already exists"
            ):
                tdg8.install_generation_zero_prefix(
                    root, destination, self.source, self.plan.sha256
                )
            self.assertEqual(_tree_bytes(destination), before)
            self.assertEqual(
                list(destination.parent.glob(".tdg8-successor-*")), []
            )

    def test_destination_arriving_at_exclusive_rename_is_preserved(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / tdg8.DESTINATION_RELATIVE
            original_rename = tdg8.hlt15._rename_exclusive

            def arriving_foreign(source: Path, target: Path) -> None:
                if target == destination:
                    target.mkdir()
                    (target / "foreign-sentinel").write_bytes(b"arrived during staging")
                original_rename(source, target)

            with patch.object(
                tdg8.hlt15,
                "_rename_exclusive",
                side_effect=arriving_foreign,
            ):
                with self.assertRaisesRegex(
                    tdg8.TDG8SuccessorRuntimeError, "destination appeared"
                ):
                    tdg8.install_generation_zero_prefix(
                        root, destination, self.source, self.plan.sha256
                    )
            self.assertEqual(
                (destination / "foreign-sentinel").read_bytes(),
                b"arrived during staging",
            )
            self.assertEqual(
                list(destination.parent.glob(".tdg8-successor-*")), []
            )

    def test_wrong_generation_one_hash_and_unsafe_destination_fail_closed(self) -> None:
        wrong = SimpleNamespace(generation=1, sha256="0" * 64)
        fake_store = SimpleNamespace(
            authenticated_checkpoint_at_generation=lambda generation: wrong
        )
        with patch.object(tdg8.pref28, "bind_terminal_store"), patch.object(
            tdg8,
            "_pref27_evidence",
            return_value=(self.source.evidence, self.source.source_identity),
        ), patch.object(tdg8, "HLT16CampaignStore", return_value=fake_store):
            with self.assertRaisesRegex(
                tdg8.TDG8SuccessorRuntimeError,
                "historical generation-one anchor differs",
            ):
                tdg8.authenticate_and_reconstruct(ROOT)

        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaisesRegex(
                tdg8.TDG8SuccessorRuntimeError, "destination differs"
            ):
                tdg8.install_generation_zero_prefix(
                    root,
                    root / "runs/foreign/calibration",
                    self.source,
                    self.plan.sha256,
                )

        with TemporaryDirectory() as temporary, TemporaryDirectory() as outside:
            root = Path(temporary)
            os.symlink(outside, root / "runs")
            destination = root / tdg8.DESTINATION_RELATIVE
            with self.assertRaisesRegex(
                tdg8.TDG8SuccessorRuntimeError, "unsafe"
            ):
                tdg8.install_generation_zero_prefix(
                    root, destination, self.source, self.plan.sha256
                )
            self.assertEqual(list(Path(outside).iterdir()), [])

    def test_bridge_initializes_and_recovers_without_pde_advance(self) -> None:
        before_members = _member_identity(self.source)
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / tdg8.DESTINATION_RELATIVE
            tdg8.install_generation_zero_prefix(
                root, destination, self.source, self.plan.sha256
            )
            store = HLT16CampaignStore(destination)
            with patch.object(
                campaign_runtime,
                "attempt_once",
                side_effect=AssertionError("bridge must not execute a PDE attempt"),
            ), patch.object(
                campaign_runtime,
                "execute_one_attempt",
                side_effect=AssertionError("bridge must not advance a member"),
            ):
                first = campaign_runtime.initialize_or_recover_generation_one(
                    self.plan, store, self.source
                )
                second = campaign_runtime.initialize_or_recover_generation_one(
                    self.plan, store, self.source
                )
                before_recovery = _tree_bytes(destination)
                recovered = campaign_runtime.recover_persisted_generation_one(
                    self.plan, store
                )
                after_recovery = _tree_bytes(destination)

            self.assertEqual(first.sha256, second.sha256)
            self.assertEqual(first.sha256, recovered.sha256)
            self.assertEqual(first.generation, 1)
            self.assertEqual(first.parent_sha256, tdg8.GENERATION_ZERO_CHECKPOINT_SHA256)
            self.assertEqual(first.protocol, "FGC-2-SF1-PROTO18")
            self.assertEqual(first.campaign_id, tdg8.CAMPAIGN_ID)
            self.assertEqual(first.plan_sha256, self.plan.sha256)
            self.assertEqual(first.event, 23)
            self.assertEqual(first.target["rational"], "3/2")
            self.assertTrue(
                all(
                    state.cursor["accepted_boundary_time"]["binary64_hex"]
                    == (23 / 16).hex()
                    for state in first.members.values()
                )
            )
            self.assertEqual(before_recovery, after_recovery)
        self.assertEqual(_member_identity(self.source), before_members)


if __name__ == "__main__":
    unittest.main()
