from __future__ import annotations

import ast
from copy import deepcopy
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
import shutil
import stat
import sys
from tempfile import TemporaryDirectory
import tomllib
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from recursive_horizons.fgc.evolution import proto19_sid2_binder as sid2  # noqa: E402


CONFIG_RAW = (ROOT / sid2.CONFIG_PATH).read_bytes()
STORE = ROOT / sid2.STORE_PATH


def canonical(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def rehash(value: dict[str, object], field: str) -> str:
    body = dict(value)
    body.pop(field, None)
    digest = sha256(canonical(body)).hexdigest()
    value[field] = digest
    return digest


def expected_compact_result(config_raw: bytes) -> dict[str, object]:
    config = tomllib.loads(config_raw.decode("utf-8"))
    return {
        "schema_version": 1,
        "artifact_id": sid2.ARTIFACT_ID,
        "project_version": "0.11.0",
        "target_protocol": "FGC-2-SF1-PROTO18",
        "classification": "post_recovery_generation8_metadata_only_binder",
        "gate_status": "pass",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "artifact_payload": {
            "predecessor": config["predecessor"],
            "live_anchor": sid2.expected_anchor(config_raw),
            "scope": config["scope"],
            "claims": config["claims"],
            "nonclaims": config["nonclaims"],
        },
    }


def tree_snapshot(root: Path) -> tuple[tuple[object, ...], ...]:
    """Test-only lexical snapshot that never follows a directory entry."""
    if root.is_symlink():
        item = root.lstat()
        return ((".", "symlink", item.st_mode, os.readlink(root)),)
    if not root.exists():
        return ()
    rows: list[tuple[object, ...]] = []

    def descend(directory: Path, prefix: str) -> None:
        with os.scandir(directory) as scan:
            entries = sorted(scan, key=lambda entry: entry.name)
        for entry in entries:
            relative = f"{prefix}/{entry.name}" if prefix else entry.name
            info = entry.stat(follow_symlinks=False)
            if stat.S_ISLNK(info.st_mode):
                rows.append((relative, "symlink", info.st_mode, os.readlink(entry.path)))
            elif stat.S_ISDIR(info.st_mode):
                rows.append((relative, "directory", info.st_mode))
                descend(Path(entry.path), relative)
            elif stat.S_ISREG(info.st_mode):
                raw = Path(entry.path).read_bytes()
                rows.append(
                    (
                        relative,
                        "file",
                        info.st_mode,
                        info.st_nlink,
                        len(raw),
                        sha256(raw).hexdigest(),
                    )
                )
            else:
                rows.append((relative, "special", info.st_mode))

    descend(root, "")
    return tuple(rows)


class SID2CompactResultTests(unittest.TestCase):
    """These tests remain valid after the ignored campaign progresses."""

    def test_expected_anchor_and_compact_validation_are_store_free(self) -> None:
        expected = expected_compact_result(CONFIG_RAW)
        encoded = sid2.canonical_result(expected)
        with patch.object(
            sid2,
            "bind_live_store",
            side_effect=AssertionError("compact validation reopened the live store"),
        ):
            self.assertEqual(sid2.validate_compact_result(CONFIG_RAW, encoded), expected)
            self.assertEqual(
                sid2.expected_anchor(CONFIG_RAW)["generation8"]["checkpoint_sha256"],
                "6d2b6e37e7c7cfd4ce64504353d886099003f959ea2aaee1a14fc778ce14fa46",
            )

    def test_compact_claim_or_anchor_promotion_is_rejected_without_store_access(self) -> None:
        for name, mutate in (
            (
                "resume",
                lambda value: value["artifact_payload"]["claims"].__setitem__(
                    "safe_to_resume_trajectory", True
                ),
            ),
            (
                "physical",
                lambda value: value["artifact_payload"]["claims"].__setitem__(
                    "physical_state_advanced", True
                ),
            ),
            (
                "cursor",
                lambda value: value["artifact_payload"]["live_anchor"][
                    "retry_state"
                ].__setitem__("cursor_sha256", "0" * 64),
            ),
        ):
            with self.subTest(name=name):
                result = deepcopy(expected_compact_result(CONFIG_RAW))
                mutate(result)
                with patch.object(
                    sid2,
                    "bind_live_store",
                    side_effect=AssertionError("compact validation reopened the live store"),
                ):
                    with self.assertRaises(sid2.SID2BinderError) as raised:
                        sid2.validate_compact_result(CONFIG_RAW, sid2.canonical_result(result))
                self.assertEqual(raised.exception.stop_id, "SID2_COMPACT_DRIFT")

    def test_config_cap_hash_state_and_claim_mutations_fail_closed_store_free(self) -> None:
        replacements = {
            "cap": ("0x1.aaa9612df9000p-10", "0x1.0000000000000p-10"),
            "hash": (
                "6d2b6e37e7c7cfd4ce64504353d886099003f959ea2aaee1a14fc778ce14fa46",
                "0" * 64,
            ),
            "state": (
                "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a",
                "1" * 64,
            ),
            "claim": ("safe_to_resume_trajectory = false", "safe_to_resume_trajectory = true"),
        }
        for name, (old, new) in replacements.items():
            with self.subTest(name=name):
                mutated = CONFIG_RAW.replace(old.encode("ascii"), new.encode("ascii"), 1)
                self.assertNotEqual(mutated, CONFIG_RAW)
                with self.assertRaises(sid2.SID2BinderError) as raised:
                    sid2.expected_anchor(mutated)
                self.assertEqual(raised.exception.stop_id, "SID2_CONFIG_DRIFT")

    def test_result_schema_and_nonclaim_boundary_are_exact(self) -> None:
        result = expected_compact_result(CONFIG_RAW)
        self.assertEqual(
            set(result),
            {
                "schema_version",
                "artifact_id",
                "project_version",
                "target_protocol",
                "classification",
                "gate_status",
                "source_config_sha256",
                "artifact_payload",
            },
        )
        payload = result["artifact_payload"]
        self.assertEqual(
            set(payload), {"predecessor", "live_anchor", "scope", "claims", "nonclaims"}
        )
        claims = payload["claims"]
        self.assertTrue(claims["metadata_only_recovery_bound"])
        self.assertTrue(claims["recovery_boundary_frozen"])
        for name in (
            "physical_state_advanced",
            "rejected_proposal_replayed",
            "trajectory_resume_authorized",
            "safe_to_resume_trajectory",
            "committed_resume_authority_present",
            "committed_resume_image_authenticated",
            "resume_authorized_by_artifact_alone",
            "GR0_calibration_completed",
            "candidate_execution_authorized",
            "physical_result_earned",
        ):
            self.assertFalse(claims[name], name)

    def test_module_import_boundary_excludes_runtime_recovery_and_branches(self) -> None:
        source = Path(sid2.__file__).read_text("utf-8")
        tree = ast.parse(source)
        imports: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module)
        for forbidden in (
            "run_fgc_pro19_event1",
            "proto19_sid1_authority",
            "hlt16_campaign_recovery",
            "progression_attempt",
            "proposal",
            "candidate",
        ):
            self.assertFalse(any(forbidden in name for name in imports), forbidden)

    def test_tracked_result_validates_compactly_when_present(self) -> None:
        result_path = ROOT / sid2.RESULT_PATH
        if not result_path.exists():
            self.skipTest("compact SID2 result has not yet been materialized")
        with patch.object(
            sid2,
            "bind_live_store",
            side_effect=AssertionError("tracked compact validation reopened the live store"),
        ):
            sid2.validate_compact_result(CONFIG_RAW, result_path.read_bytes())


class SID2LiveStoreTests(unittest.TestCase):
    """Exact temporal-boundary tests skip only after lawful later progression."""

    @classmethod
    def setUpClass(cls) -> None:
        checkpoints = STORE / "checkpoints"
        journals = STORE / "journal"
        if not checkpoints.is_dir() or not journals.is_dir():
            raise unittest.SkipTest("SID2 ignored live store is unavailable")
        checkpoint_names = sorted(path.name for path in checkpoints.iterdir())
        journal_names = sorted(path.name for path in journals.iterdir())
        if (
            len(checkpoint_names) != 9
            or len(journal_names) != 9
            or not checkpoint_names[-1].startswith("00000000000000000008-")
            or not journal_names[-1].startswith("00000000000000000008-")
        ):
            raise unittest.SkipTest("campaign lawfully moved beyond the SID2 generation-eight boundary")
        # At the exact boundary, corruption or implementation drift must fail,
        # not silently turn into a skip.
        cls.anchor = sid2.bind_live_store(ROOT)

    def clone_store(self, temporary: Path) -> Path:
        repository = temporary / "repository"
        destination = repository / sid2.STORE_PATH
        destination.parent.mkdir(parents=True)
        shutil.copytree(STORE, destination, symlinks=True)
        return repository

    def assert_rejected_unchanged(self, repository: Path) -> sid2.SID2BinderError:
        store = repository / sid2.STORE_PATH
        before = tree_snapshot(store)
        with self.assertRaises(sid2.SID2BinderError) as raised:
            sid2.bind_live_store(repository)
        self.assertEqual(tree_snapshot(store), before)
        return raised.exception

    def generation8_path(self, repository: Path) -> Path:
        return next((repository / sid2.STORE_PATH / "checkpoints").glob("00000000000000000008-*.json"))

    def test_live_binding_and_result_build_are_exact_and_nonmutating(self) -> None:
        before = tree_snapshot(STORE)
        anchor = sid2.bind_live_store(ROOT)
        self.assertEqual(anchor, sid2.expected_anchor(CONFIG_RAW))
        self.assertEqual(len(anchor["member_descriptors"]), 6)
        self.assertEqual(anchor["store"]["regular_leaf_count"], 52)
        self.assertEqual(
            anchor["store"]["lexical_inventory_sha256"],
            "8196b088a19f70ae1016a3d9dddb508add1ce1ce5cd227d105c954bd0c151c27",
        )
        self.assertEqual(
            anchor["physical_identity"]["evolution_state_sha256"],
            "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a",
        )
        self.assertEqual(anchor["retry_state"]["mode"], "RETRY_PENDING")
        result = sid2.build_sid2_result(CONFIG_RAW, ROOT)
        self.assertEqual(result, expected_compact_result(CONFIG_RAW))
        self.assertEqual(tree_snapshot(STORE), before)

    def test_parent_cursor_cap_hash_and_physical_state_mutations_fail_closed(self) -> None:
        def parent(repository: Path) -> None:
            path = self.generation8_path(repository)
            value = json.loads(path.read_bytes())
            value["parent_sha256"] = "0" * 64
            address = rehash(value, "checkpoint_sha256")
            path.write_bytes(canonical(value))
            path.rename(path.with_name(f"{8:020d}-{address}.json"))

        def cursor(repository: Path) -> None:
            path = self.generation8_path(repository)
            value = json.loads(path.read_bytes())
            cursor_value = value["members"]["RK4-2049"]["cursor"]
            cursor_value["cursor_chain_parent_sha256"] = "0" * 64
            rehash(cursor_value, "cursor_chain_sha256")
            address = rehash(value, "checkpoint_sha256")
            path.write_bytes(canonical(value))
            path.rename(path.with_name(f"{8:020d}-{address}.json"))

        def cap(repository: Path) -> None:
            path = self.generation8_path(repository)
            value = json.loads(path.read_bytes())
            value["members"]["RK4-2049"]["pending_cap_hex"] = "0x1.0000000000000p-10"
            address = rehash(value, "checkpoint_sha256")
            path.write_bytes(canonical(value))
            path.rename(path.with_name(f"{8:020d}-{address}.json"))

        def raw_hash(repository: Path) -> None:
            journal = next(
                (repository / sid2.STORE_PATH / "journal").glob(
                    "00000000000000000008-*.journal"
                )
            )
            value = json.loads(journal.read_bytes())
            value["record_sha256"] = "0" * 64
            journal.write_bytes(canonical(value))

        def physical_u(repository: Path) -> None:
            payload = (
                repository
                / sid2.STORE_PATH
                / "payloads/0e48e0e2396ea7d9ef79d1641e6944a2dfdaa3716bb65e0b882b97156ff9372f.npz"
            )
            with np.load(BytesIO(payload.read_bytes()), allow_pickle=False) as archive:
                arrays = {name: archive[name].copy() for name in sid2.ARRAY_NAMES}
            arrays["u"][0, 0] = np.nextafter(arrays["u"][0, 0], np.inf)
            stream = BytesIO()
            np.savez_compressed(stream, **arrays)
            payload.write_bytes(stream.getvalue())

        for name, mutate in (
            ("parent", parent),
            ("cursor", cursor),
            ("cap", cap),
            ("raw_hash", raw_hash),
            ("physical_u", physical_u),
        ):
            with self.subTest(name=name), TemporaryDirectory(prefix="sid2-mutation-") as directory:
                repository = self.clone_store(Path(directory))
                mutate(repository)
                self.assert_rejected_unchanged(repository)

    def test_symlinked_repository_directory_and_leaf_fail_closed(self) -> None:
        with TemporaryDirectory(prefix="sid2-symlink-root-") as directory:
            repository = self.clone_store(Path(directory))
            link = Path(directory) / "repository-link"
            os.symlink(repository, link)
            with self.assertRaises(sid2.SID2BinderError) as raised:
                sid2.bind_live_store(link)
            self.assertEqual(raised.exception.stop_id, "SID2_PATH_UNSAFE")

        with TemporaryDirectory(prefix="sid2-symlink-dir-") as directory:
            repository = self.clone_store(Path(directory))
            states = repository / sid2.STORE_PATH / "states"
            relocated = repository / "relocated-states"
            states.rename(relocated)
            os.symlink(relocated, states)
            self.assert_rejected_unchanged(repository)

        with TemporaryDirectory(prefix="sid2-symlink-leaf-") as directory:
            repository = self.clone_store(Path(directory))
            checkpoint = self.generation8_path(repository)
            relocated = repository / "relocated-checkpoint.json"
            checkpoint.rename(relocated)
            os.symlink(relocated, checkpoint)
            self.assert_rejected_unchanged(repository)

    def test_foreign_fork_and_staging_leaves_fail_closed(self) -> None:
        def foreign(repository: Path) -> None:
            (repository / sid2.STORE_PATH / "foreign").write_bytes(b"foreign")

        def fork(repository: Path) -> None:
            checkpoint = self.generation8_path(repository)
            shutil.copy2(
                checkpoint,
                checkpoint.with_name(f"{8:020d}-{'0' * 64}.json"),
            )

        def stage(repository: Path) -> None:
            (repository / sid2.STORE_PATH / "journal" / f".foreign.hlt16-stage-{'0' * 64}").write_bytes(b"")

        def sibling_stage(repository: Path) -> None:
            (repository / sid2.STORE_PATH).parent.joinpath(
                ".calibration.hlt16-stage-foreign"
            ).mkdir()

        for name, mutate in (
            ("foreign", foreign),
            ("fork", fork),
            ("stage", stage),
            ("sibling_stage", sibling_stage),
        ):
            with self.subTest(name=name), TemporaryDirectory(prefix="sid2-tree-") as directory:
                repository = self.clone_store(Path(directory))
                mutate(repository)
                self.assert_rejected_unchanged(repository)

    def test_active_writer_and_terminal_lock_fail_closed(self) -> None:
        for name in ("active-write.lock", "terminal.lock"):
            with self.subTest(name=name), TemporaryDirectory(prefix="sid2-lock-") as directory:
                repository = self.clone_store(Path(directory))
                (repository / sid2.STORE_PATH / "locks" / name).write_bytes(b"{}")
                error = self.assert_rejected_unchanged(repository)
                self.assertEqual(error.stop_id, "SID2_TREE_DRIFT")

    def test_hard_link_and_special_file_fail_closed(self) -> None:
        with TemporaryDirectory(prefix="sid2-hardlink-") as directory:
            repository = self.clone_store(Path(directory))
            states = repository / sid2.STORE_PATH / "states"
            source = next(states.glob("*.json"))
            os.link(source, states / f"{'0' * 64}.json")
            error = self.assert_rejected_unchanged(repository)
            self.assertEqual(error.stop_id, "SID2_PATH_UNSAFE")

        if not hasattr(os, "mkfifo"):
            self.skipTest("platform has no FIFO constructor")
        with TemporaryDirectory(prefix="sid2-fifo-") as directory:
            repository = self.clone_store(Path(directory))
            fifo = repository / sid2.STORE_PATH / "states" / f"{'0' * 64}.json"
            try:
                os.mkfifo(fifo)
            except OSError as error:
                self.skipTest(f"FIFO unsupported by temporary filesystem: {error}")
            error = self.assert_rejected_unchanged(repository)
            self.assertEqual(error.stop_id, "SID2_PATH_UNSAFE")

    def test_descriptor_and_payload_orphans_or_omissions_fail_closed(self) -> None:
        for name, mutate in (
            (
                "descriptor_orphan",
                lambda repository: (
                    repository / sid2.STORE_PATH / "states" / f"{'0' * 64}.json"
                ).write_bytes(canonical({})),
            ),
            (
                "descriptor_missing",
                lambda repository: next(
                    (repository / sid2.STORE_PATH / "states").glob(
                        "6b0dd985ab95ac52950348df5079a6b8fc4f5e7a78cf95031237cb0155c49c56.json"
                    )
                ).unlink(),
            ),
            (
                "payload_orphan",
                lambda repository: (
                    repository / sid2.STORE_PATH / "payloads" / f"{'0' * 64}.npz"
                ).write_bytes(b"foreign"),
            ),
        ):
            with self.subTest(name=name), TemporaryDirectory(prefix="sid2-orphan-") as directory:
                repository = self.clone_store(Path(directory))
                mutate(repository)
                self.assert_rejected_unchanged(repository)


if __name__ == "__main__":
    unittest.main()
