"""Focused prospective authority controls for TDG11-MSEL1."""

from __future__ import annotations

import ast
from contextlib import contextmanager, ExitStack
from copy import deepcopy
from hashlib import sha256
import importlib.util
import inspect
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from recursive_horizons.evidence_io import (
    canonical_json_bytes,
    InspectTree,
    InspectWorktree,
    TreeEntry,
    WorktreeInspection,
    git_read,
)
from recursive_horizons.fgc.evolution import tdg11_msel1_authority as authority
from recursive_horizons.fgc.evolution import tdg11_msel1_contract as contract
from recursive_horizons.fgc.evolution.numerical_engine import PRIMARY_METHOD


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/tdg11_msel1_authority.py"
SCRIPT_PATH = ROOT / "scripts/reproduce_fgc_tdg11_msel1_frz1.py"
_WRITER_MARKERS = (
    "publish_journal",
    "persist_snapshot",
    "active_writer_capability",
    "acquire_writer",
    "repair_stages=True",
    "build_static_gr0_shells",
)


def _git_env() -> dict[str, str]:
    env = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    env.update(
        {
            "GIT_AUTHOR_EMAIL": "tdg11@example.test",
            "GIT_AUTHOR_NAME": "tdg11-test",
            "GIT_COMMITTER_EMAIL": "tdg11@example.test",
            "GIT_COMMITTER_NAME": "tdg11-test",
            "GIT_CONFIG_GLOBAL": os.devnull,
            "GIT_CONFIG_NOSYSTEM": "1",
        }
    )
    return env


def _git(repo: Path, args: tuple[str, ...]) -> None:
    subprocess.run(
        ["git", "-c", "core.hooksPath=/dev/null", *args],
        cwd=repo,
        check=True,
        env=_git_env(),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _git_head(repo: Path) -> str:
    return (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, env=_git_env())
        .decode("ascii")
        .strip()
    )


def _write(root: Path, relative: str, payload: bytes) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _fixture_environment() -> dict[str, str]:
    environment = dict(contract.MINIMUM_ENVIRONMENT)
    environment.update(dict.fromkeys(contract.ENVIRONMENT_IMAGE_KEYS, "synthetic"))
    environment["python_executable_sha256"] = "a" * 64
    environment["numpy_extension_sha256"] = "b" * 64
    return environment


def _compact_bundle(root: Path) -> dict[str, object]:
    implementation_path = "src/owned.py"
    reference_path = "ref.txt"
    _write(root, implementation_path, b"owned-source\n")
    _write(root, reference_path, b"pinned-ref\n")
    implementation = [
        {
            "path": implementation_path,
            "sha256": sha256(b"owned-source\n").hexdigest(),
        }
    ]
    references = {reference_path: sha256(b"pinned-ref\n").hexdigest()}
    with (
        patch.object(contract, "IMPLEMENTATION_PATHS", (implementation_path,)),
        patch.object(contract, "PINNED_REFERENCES", references),
    ):
        config = contract.expected_config(
            environment=_fixture_environment(), implementation=implementation
        )
        config_raw = contract.render_config(config)
        _write(root, contract.CONFIG_PATH, config_raw)
        result = contract.compact_freeze(config_raw)
        _write(root, contract.RESULT_PATH, canonical_json_bytes(result))
        return result


def _load_script():
    spec = importlib.util.spec_from_file_location(
        "reproduce_fgc_tdg11_msel1_frz1", SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _successor_repo(
    temporary: str, *, extra: str | None = None
) -> tuple[Path, str, str]:
    repo = Path(temporary).resolve()
    subprocess.run(
        ["git", "-c", "init.defaultBranch=main", "init"],
        cwd=repo,
        check=True,
        env=_git_env(),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    _write(repo, "seed.txt", b"seed\n")
    _git(repo, ("add", "seed.txt"))
    _git(repo, ("commit", "-m", "base"))
    base = _git_head(repo)
    _write(repo, "delta.txt", b"delta\n")
    added = ["delta.txt"]
    if extra is not None:
        _write(repo, extra, b"extra\n")
        added.append(extra)
    _git(repo, ("add", *added))
    _git(repo, ("commit", "-m", "successor"))
    return repo, base, _git_head(repo)


class MSEL1AuthorityTests(unittest.TestCase):
    def test_compact_accepts_one_emitter_newline_but_not_noncanonical_json(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            _compact_bundle(root)
            path = root / contract.RESULT_PATH
            original = path.read_bytes()
            path.write_bytes(original + b"\n")
            with (
                patch.object(contract, "IMPLEMENTATION_PATHS", ("src/owned.py",)),
                patch.object(
                    contract,
                    "PINNED_REFERENCES",
                    {"ref.txt": sha256(b"pinned-ref\n").hexdigest()},
                ),
            ):
                self.assertFalse(
                    authority.validate_compact_bundle(root)["diagnostic_executed"]
                )
                path.write_bytes(original + b"\n\n")
                with self.assertRaises(authority.MSEL1AuthorityError):
                    authority.validate_compact_bundle(root)

    def test_only_software_images_allow_existing_package_hardlinks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / "software.bin").write_bytes(b"software")
            os.link(root / "software.bin", root / "software-link.bin")
            self.assertEqual(
                authority._hash_software_file(str(root / "software.bin")),
                sha256(b"software").hexdigest(),
            )
            with self.assertRaises(authority.MSEL1AuthorityError):
                authority._store_leaf_bytes(root, "software.bin")
            with self.assertRaises(authority.MSEL1AuthorityError):
                authority._store_leaf_bytes(root, str(root / "software.bin"))

    def test_software_same_size_race_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary).resolve() / "software.bin"
            path.write_bytes(b"AAAA")
            before = path.stat()
            original_read = os.read
            changed = False

            def raced_read(fd, size):
                nonlocal changed
                data = original_read(fd, size)
                if not changed:
                    changed = True
                    path.write_bytes(b"BBBB")
                    os.utime(
                        path,
                        ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000),
                    )
                return data

            with patch.object(authority.os, "read", side_effect=raced_read):
                with self.assertRaises(authority.MSEL1AuthorityError):
                    authority._hash_software_file(str(path))

    def test_live_authority_requires_predecessor_inspection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            with (
                patch.object(authority, "_require_committed_image"),
                self._authorize_after_git(_fixture_environment()),
                patch.object(
                    authority,
                    "inspect_predecessors",
                    side_effect=authority.MSEL1AuthorityError("predecessor mismatch"),
                ) as inspected,
            ):
                with self.assertRaisesRegex(
                    authority.MSEL1AuthorityError, "predecessor mismatch"
                ):
                    authority.authorize_execution(root, authority_commit="a" * 40)
                inspected.assert_called_once_with(root)

    def test_exported_api_and_receipt_schema(self) -> None:
        self.assertEqual(
            set(authority.__all__),
            {
                "MSEL1Authority",
                "MSEL1AuthorityError",
                "PredecessorInspection",
                "RestoredPredecessor",
                "authorize_execution",
                "emit_config_bytes",
                "emit_result_object",
                "inspect_predecessors",
                "observe_environment",
                "require_output_absent",
                "restore_predecessor",
                "snapshot_store",
                "validate_compact_bundle",
            },
        )
        receipt = authority.MSEL1Authority(
            authority_commit="a" * 40,
            config_sha256="b" * 64,
            result_sha256="c" * 64,
            environment=_fixture_environment(),
            store_snapshot=(
                contract.STORE_LEAF_COUNT,
                contract.STORE_SHA256,
            ),
            config={"artifact_id": contract.ARTIFACT_ID},
        )
        self.assertEqual(
            receipt.__dataclass_fields__.keys(),
            {
                "authority_commit",
                "config_sha256",
                "result_sha256",
                "environment",
                "store_snapshot",
                "config",
            },
        )
        self.assertEqual(
            authority.RestoredPredecessor.__dataclass_fields__.keys(),
            {
                "retry",
                "member",
                "checkpoint_sha256",
                "descriptor_sha256",
                "fingerprint",
            },
        )
        self.assertTrue(receipt.__dataclass_params__.frozen)
        self.assertIn("config", receipt.__dataclass_fields__)

    def test_module_import_is_store_and_numpy_free(self) -> None:
        tree = ast.parse(MODULE_PATH.read_text(encoding="utf-8"))
        imported: set[str] = set()
        for node in tree.body:
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    imported.add(node.module)
        for name in imported:
            self.assertFalse(
                any(
                    part in name
                    for part in (
                        "numpy",
                        "hlt16",
                        "proto15",
                        "proto19",
                        "campaign_runtime",
                        "subprocess",
                    )
                ),
                name,
            )
        source = MODULE_PATH.read_text(encoding="utf-8")
        for marker in _WRITER_MARKERS:
            self.assertNotIn(marker, source)
        restore_source = inspect.getsource(authority.restore_predecessor)
        self.assertNotIn("build_static_gr0_shells", restore_source)
        self.assertIn("deepcopy", restore_source)
        self.assertIn("templates[contract.MEMBER_KEY]", restore_source)

    def test_compact_bundle_validates_without_live_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            _compact_bundle(root)
            with (
                patch.object(contract, "IMPLEMENTATION_PATHS", ("src/owned.py",)),
                patch.object(
                    contract,
                    "PINNED_REFERENCES",
                    {"ref.txt": sha256(b"pinned-ref\n").hexdigest()},
                ),
                patch.object(authority, "git_read", side_effect=AssertionError("git")),
                patch.object(
                    authority,
                    "snapshot_store",
                    side_effect=AssertionError("store"),
                ),
                patch.object(
                    authority,
                    "observe_environment",
                    side_effect=AssertionError("environment"),
                ),
                patch(
                    "subprocess.Popen",
                    side_effect=AssertionError("process"),
                ),
            ):
                bundle = authority.validate_compact_bundle(root)
            config_digest = sha256(
                (root / contract.CONFIG_PATH).read_bytes()
            ).hexdigest()
        self.assertIs(bundle["diagnostic_executed"], False)
        self.assertIs(bundle["selection_result_earned"], False)
        self.assertEqual(bundle["config_sha256"], config_digest)

    def test_compact_rejects_source_and_reference_hash_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            _compact_bundle(root)
            (root / "src/owned.py").write_bytes(b"mutated\n")
            with (
                patch.object(contract, "IMPLEMENTATION_PATHS", ("src/owned.py",)),
                patch.object(
                    contract,
                    "PINNED_REFERENCES",
                    {"ref.txt": sha256(b"pinned-ref\n").hexdigest()},
                ),
                self.assertRaisesRegex(
                    authority.MSEL1AuthorityError, "implementation binding"
                ),
            ):
                authority.validate_compact_bundle(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            _compact_bundle(root)
            (root / "ref.txt").write_bytes(b"other\n")
            with (
                patch.object(contract, "IMPLEMENTATION_PATHS", ("src/owned.py",)),
                patch.object(
                    contract,
                    "PINNED_REFERENCES",
                    {"ref.txt": sha256(b"pinned-ref\n").hexdigest()},
                ),
                self.assertRaisesRegex(
                    authority.MSEL1AuthorityError, "pinned reference"
                ),
            ):
                authority.validate_compact_bundle(root)

    def test_compact_source_has_no_live_dependencies(self) -> None:
        compact = inspect.getsource(authority.validate_compact_bundle)
        tracked = inspect.getsource(authority._tracked_payload)
        for source in (compact, tracked):
            for forbidden in (
                "git_read",
                "InspectTree",
                "snapshot_store",
                "observe_environment",
                "HLT16",
                "restore_member",
                "numpy",
                "/bin/ps",
            ):
                self.assertNotIn(forbidden, source)

    def test_output_absence_rejects_existing_staging_and_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / contract.OUTPUT_NAMESPACE).mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.MSEL1AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            parent = (root / contract.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            (parent / f"{authority._STAGING_PREFIX}dead").mkdir()
            with self.assertRaisesRegex(
                authority.MSEL1AuthorityError, "staging namespace already exists"
            ):
                authority.require_output_absent(root)
        with (
            tempfile.TemporaryDirectory() as temporary,
            tempfile.TemporaryDirectory() as foreign,
        ):
            root = Path(temporary).resolve()
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.MSEL1AuthorityError, "output parent is unsafe"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            parent = (root / contract.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            (root / contract.OUTPUT_NAMESPACE).symlink_to("missing-target")
            with self.assertRaisesRegex(
                authority.MSEL1AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            parent = (root / contract.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            authority.require_output_absent(root)
            self.assertFalse((root / contract.OUTPUT_NAMESPACE).exists())

    def test_store_snapshot_uses_historical_domain_and_rejects_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            store = root / "store"
            _write(store, "a.txt", b"A")
            _write(store, "z.txt", b"Z")
            _write(store, "sub/c.txt", b"C")
            digest = sha256(b"TDG9-LOC1-STORE-SNAPSHOT-v1\n")
            for relative, payload in (
                ("a.txt", b"A"),
                ("z.txt", b"Z"),
                ("sub/c.txt", b"C"),
            ):
                digest.update(
                    (relative + "\0" + sha256(payload).hexdigest() + "\n").encode(
                        "ascii"
                    )
                )
            with patch.object(contract, "STORE_PATH", "store"):
                pair = authority._hash_store_tree(store)
            self.assertEqual(pair, (3, digest.hexdigest()))
            with (
                patch.object(contract, "STORE_PATH", "store"),
                self.assertRaisesRegex(
                    authority.MSEL1AuthorityError, "store snapshot differs"
                ),
            ):
                authority.snapshot_store(root)
            with (
                patch.object(contract, "STORE_PATH", "store"),
                patch.object(contract, "STORE_LEAF_COUNT", 3),
                patch.object(contract, "STORE_SHA256", pair[1]),
            ):
                self.assertEqual(authority.snapshot_store(root), pair)
            (store / "link").symlink_to("a.txt")
            with self.assertRaisesRegex(authority.MSEL1AuthorityError, "symlink"):
                authority._hash_store_tree(store)

    def test_git_blob_id_matches_committed_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repo, _base, head = _successor_repo(temporary)
            tree = git_read(repo, InspectTree(commit=head))
            payload = (repo / "delta.txt").read_bytes()
            match = next(item for item in tree if item.path == "delta.txt")
            self.assertEqual(authority._git_blob_id(payload), match.object_id)

    @contextmanager
    def _authorize_after_git(self, environment: dict[str, str]):
        config = {"environment": dict(environment)}
        compact = {"config": config}
        with ExitStack() as stack:
            stack.enter_context(
                patch.object(
                    authority,
                    "_tracked_payload",
                    return_value=(b"cfg", b"res", compact),
                )
            )
            stack.enter_context(
                patch.object(
                    authority, "observe_environment", return_value=dict(environment)
                )
            )
            stack.enter_context(patch.object(authority, "_require_process_preflight"))
            stack.enter_context(patch.object(authority, "require_output_absent"))
            stack.enter_context(patch.object(authority, "inspect_predecessors"))
            stack.enter_context(
                patch.object(
                    authority,
                    "snapshot_store",
                    return_value=(contract.STORE_LEAF_COUNT, contract.STORE_SHA256),
                )
            )
            yield

    def test_authorize_accepts_direct_successor_and_matching_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            environment = _fixture_environment()
            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                self._authorize_after_git(environment),
            ):
                receipt = authority.authorize_execution(repo, authority_commit=head)
        self.assertEqual(receipt.authority_commit, head)
        self.assertEqual(receipt.config, {"environment": environment})
        self.assertEqual(receipt.config_sha256, sha256(b"cfg").hexdigest())
        self.assertEqual(receipt.store_snapshot[0], contract.STORE_LEAF_COUNT)

    def test_altered_parent_merge_delta_dirty_and_hidden_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            environment = _fixture_environment()
            with (
                patch.object(contract, "BASE_COMMIT", "0" * 40),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                self._authorize_after_git(environment),
                self.assertRaisesRegex(
                    authority.MSEL1AuthorityError, "direct successor"
                ),
            ):
                authority.authorize_execution(repo, authority_commit=head)
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            _git(repo, ("checkout", "-b", "side"))
            _write(repo, "seed.txt", b"side\n")
            _git(repo, ("commit", "-am", "side"))
            _git(repo, ("checkout", "main"))
            _git(repo, ("merge", "side", "--no-ff", "-m", "merge"))
            merge = _git_head(repo)
            environment = _fixture_environment()
            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                self._authorize_after_git(environment),
                self.assertRaisesRegex(
                    authority.MSEL1AuthorityError, "direct successor"
                ),
            ):
                authority.authorize_execution(repo, authority_commit=merge)
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary, extra="bonus.txt")
            environment = _fixture_environment()
            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                self._authorize_after_git(environment),
                self.assertRaisesRegex(authority.MSEL1AuthorityError, "delta"),
            ):
                authority.authorize_execution(repo, authority_commit=head)
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            _write(repo, "delta.txt", b"dirty\n")
            environment = _fixture_environment()
            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                self._authorize_after_git(environment),
                self.assertRaisesRegex(authority.MSEL1AuthorityError, "not clean"),
            ):
                authority.authorize_execution(repo, authority_commit=head)
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            _git(repo, ("update-index", "--skip-worktree", "delta.txt"))
            environment = _fixture_environment()
            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                self._authorize_after_git(environment),
                self.assertRaisesRegex(authority.MSEL1AuthorityError, "not clean"),
            ):
                authority.authorize_execution(repo, authority_commit=head)
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            environment = _fixture_environment()
            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                self._authorize_after_git(environment),
                self.assertRaisesRegex(authority.MSEL1AuthorityError, "40-hex"),
            ):
                authority.authorize_execution(repo, authority_commit="HEAD")

    def test_nonregular_tree_mode_and_byte_mismatch_are_rejected(self) -> None:
        environment = _fixture_environment()
        tree = (
            TreeEntry(
                mode="120000",
                kind="blob",
                object_id="a" * 40,
                path="delta.txt",
            ),
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            with (
                patch.object(authority, "git_read") as reader,
                patch.object(contract, "BASE_COMMIT", "b" * 40),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                self._authorize_after_git(environment),
                self.assertRaisesRegex(authority.MSEL1AuthorityError, "nonregular"),
            ):
                reader.side_effect = [
                    "c" * 40,
                    ("b" * 40,),
                    (SimpleNamespace(path="delta.txt"),),
                    WorktreeInspection((), ()),
                    tree,
                ]
                authority.authorize_execution(root, authority_commit="c" * 40)

        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            original_git = git_read
            original_read = authority.read_regular_file

            def git_clean(root, operation):
                if isinstance(operation, InspectWorktree):
                    return WorktreeInspection((), ())
                return original_git(root, operation)

            def read_mutated(root, relative, **kwargs):
                if relative == "delta.txt":
                    return b"mutated-bytes\n"
                return original_read(root, relative, **kwargs)

            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                patch.object(authority, "git_read", side_effect=git_clean),
                patch.object(authority, "read_regular_file", side_effect=read_mutated),
                self._authorize_after_git(environment),
                self.assertRaisesRegex(
                    authority.MSEL1AuthorityError, "committed bytes differ"
                ),
            ):
                authority.authorize_execution(repo, authority_commit=head)

    def test_environment_and_namespace_fail_closed_on_live_path(self) -> None:
        environment = _fixture_environment()
        other = dict(environment)
        other["python_version"] = "3.14.6"
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            config = {"environment": environment}
            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                patch.object(
                    authority,
                    "_tracked_payload",
                    return_value=(b"cfg", b"res", {"config": config}),
                ),
                patch.object(authority, "observe_environment", return_value=other),
                patch.object(authority, "_require_process_preflight"),
                patch.object(authority, "require_output_absent"),
                patch.object(
                    authority,
                    "snapshot_store",
                    return_value=(contract.STORE_LEAF_COUNT, contract.STORE_SHA256),
                ),
                self.assertRaisesRegex(authority.MSEL1AuthorityError, "environment"),
            ):
                authority.authorize_execution(repo, authority_commit=head)
        with tempfile.TemporaryDirectory() as temporary:
            repo, base, head = _successor_repo(temporary)
            parent = (repo / contract.OUTPUT_NAMESPACE).parent
            parent.mkdir(parents=True)
            (parent / f"{authority._STAGING_PREFIX}x").mkdir()
            with (
                patch.object(contract, "BASE_COMMIT", base),
                patch.object(contract, "DELTA_PATHS", ("delta.txt",)),
                patch.object(
                    authority,
                    "_tracked_payload",
                    return_value=(
                        b"cfg",
                        b"res",
                        {"config": {"environment": environment}},
                    ),
                ),
                patch.object(
                    authority, "observe_environment", return_value=environment
                ),
                patch.object(authority, "_require_process_preflight"),
                patch.object(
                    authority,
                    "snapshot_store",
                    return_value=(contract.STORE_LEAF_COUNT, contract.STORE_SHA256),
                ),
                self.assertRaisesRegex(
                    authority.MSEL1AuthorityError, "staging namespace"
                ),
            ):
                authority.authorize_execution(repo, authority_commit=head)

    def test_process_preflight_excludes_self_and_rejects_foreign_runners(self) -> None:
        self_pid = os.getpid()
        parent = os.getppid()
        rows = (
            (
                self_pid,
                parent,
                "python -m unittest tests/test_fgc_tdg11_msel1_authority",
            ),
            (parent, 1, "/sbin/launchd"),
            (4242, 1, "python scripts/run_fgc_tdg10_qa2.py --run"),
        )
        with patch.object(authority, "_process_rows", return_value=rows):
            with self.assertRaisesRegex(
                authority.MSEL1AuthorityError, "scientific runner"
            ):
                authority._require_process_preflight()
        safe = (
            (
                self_pid,
                parent,
                "python -m unittest tests/test_fgc_tdg11_msel1_authority",
            ),
            (parent, 1, "/sbin/launchd"),
        )
        with patch.object(authority, "_process_rows", return_value=safe):
            authority._require_process_preflight()
        review = safe + ((99, 1, "codex review --wait"),)
        with patch.object(authority, "_process_rows", return_value=review):
            with self.assertRaisesRegex(authority.MSEL1AuthorityError, "review"):
                authority._require_process_preflight()
        pytest_row = safe + ((100, 1, "python -m pytest tests"),)
        with patch.object(authority, "_process_rows", return_value=pytest_row):
            with self.assertRaisesRegex(
                authority.MSEL1AuthorityError, "test partition"
            ):
                authority._require_process_preflight()

    def test_inspect_predecessors_validates_generations_and_never_restores_rejections(
        self,
    ) -> None:
        specs = []
        leaves: dict[str, bytes] = {}
        checkpoints: dict[int, SimpleNamespace] = {}
        for retry, generation, journal_seq, rejection_seq, prior in (
            (3, 9, 10, 11, 2),
            (4, 10, 12, 13, 3),
            (5, 11, 14, 15, 4),
        ):
            checkpoint_sha = f"{retry:x}" * 64
            journal_sha = f"{retry + 3:x}" * 64
            rejection_sha = f"{retry + 6:x}" * 64
            checkpoint_raw = json.dumps(
                {
                    "generation": generation,
                    "checkpoint_sha256": checkpoint_sha,
                    "journal_sequence": journal_seq,
                    "journal_tip_sha256": journal_sha,
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
            journal_raw = json.dumps(
                {
                    "kind": "cursor_transition",
                    "sequence": journal_seq,
                    "record_sha256": journal_sha,
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
            rejection_raw = json.dumps(
                {
                    "kind": "tdg6_rejection",
                    "sequence": rejection_seq,
                    "record_sha256": rejection_sha,
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
            specs.append(
                {
                    "retry": retry,
                    "generation": generation,
                    "journal_tip_sequence": journal_seq,
                    "rejection_sequence": rejection_seq,
                    "width_hex": "0x1p-1",
                    "checkpoint_sha256": checkpoint_sha,
                    "checkpoint_raw_sha256": sha256(checkpoint_raw).hexdigest(),
                    "journal_tip_sha256": journal_sha,
                    "journal_tip_raw_sha256": sha256(journal_raw).hexdigest(),
                    "rejection_sha256": rejection_sha,
                    "rejection_raw_sha256": sha256(rejection_raw).hexdigest(),
                    "prior_retry_count": prior,
                }
            )
            leaves[
                f"{contract.STORE_PATH}/checkpoints/{generation:020d}-{checkpoint_sha}.json"
            ] = checkpoint_raw
            leaves[
                f"{contract.STORE_PATH}/journal/{journal_seq:020d}-{journal_sha}.journal"
            ] = journal_raw
            leaves[
                f"{contract.STORE_PATH}/journal/{rejection_seq:020d}-{rejection_sha}.journal"
            ] = rejection_raw
            member = SimpleNamespace(
                descriptor_sha256=contract.DESCRIPTOR_SHA256,
                pending_owner="temporal",
                cursor={
                    "accepted_boundary_time": {
                        "binary64_hex": contract.ACCEPTED_TIME_HEX
                    },
                    "mode": "RETRY_PENDING",
                },
                ledger={"current_macro_step_temporal_retry_count": prior},
            )
            checkpoints[generation] = SimpleNamespace(
                sha256=checkpoint_sha,
                generation=generation,
                journal_sequence=journal_seq,
                journal_tip_sha256=journal_sha,
                members={contract.MEMBER_KEY: member},
            )

        class FakeStore:
            def __init__(self, path):
                self.path = path

            def authenticated_checkpoint_at_generation(self, generation):
                return checkpoints[generation]

            def publish_journal(self, *args, **kwargs):
                raise AssertionError("writer")

        def fake_leaf(_root, relative):
            return leaves[relative]

        with (
            patch.object(contract, "replay_specs", return_value=tuple(specs)),
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_store.HLT16CampaignStore",
                FakeStore,
            ),
            patch.object(authority, "_store_leaf_bytes", side_effect=fake_leaf),
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_runtime.restore_member_with_overlay",
                side_effect=AssertionError("restore"),
            ),
        ):
            inspected = authority.inspect_predecessors(ROOT)
        self.assertEqual(tuple(item.generation for item in inspected), (9, 10, 11))
        self.assertEqual(tuple(item.retry for item in inspected), (3, 4, 5))
        self.assertTrue(all(item.pending_owner == "temporal" for item in inspected))

        broken = deepcopy(checkpoints)
        broken[10] = SimpleNamespace(**{**broken[10].__dict__, "generation": 11})

        class WrongStore(FakeStore):
            def authenticated_checkpoint_at_generation(self, generation):
                return broken[generation]

        with (
            patch.object(contract, "replay_specs", return_value=tuple(specs)),
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_store.HLT16CampaignStore",
                WrongStore,
            ),
            patch.object(authority, "_store_leaf_bytes", side_effect=fake_leaf),
            self.assertRaisesRegex(
                authority.MSEL1AuthorityError, "checkpoint identity"
            ),
        ):
            authority.inspect_predecessors(ROOT)

        def journal_swap(_root, relative):
            payload = leaves[relative]
            if relative.endswith(".journal") and "00000000000000000011-" in relative:
                return json.dumps(
                    {
                        "kind": "cursor_transition",
                        "sequence": 11,
                        "record_sha256": specs[0]["rejection_sha256"],
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            return payload

        with (
            patch.object(contract, "replay_specs", return_value=tuple(specs)),
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_store.HLT16CampaignStore",
                FakeStore,
            ),
            patch.object(authority, "_store_leaf_bytes", side_effect=journal_swap),
            self.assertRaisesRegex(
                authority.MSEL1AuthorityError, "leaf hash|rejection"
            ),
        ):
            authority.inspect_predecessors(ROOT)

    def test_restore_predecessor_uses_supplied_template_and_rejects_hash_drift(
        self,
    ) -> None:
        import numpy as np

        from recursive_horizons.fgc.evolution.numerical_engine import (
            array_content_sha256,
        )

        state_u = np.zeros((9, 1), dtype=np.float64)
        coordinates = np.linspace(0.0, 1.0, 9, dtype=np.float64)
        member = SimpleNamespace(
            key=contract.MEMBER_KEY,
            method_label="RK4",
            integrator_id=PRIMARY_METHOD,
            spatial_order=4,
            point_count=contract.POINT_COUNT,
            time=float.fromhex(contract.ACCEPTED_TIME_HEX),
            step_index=1,
            transaction_serial=2,
            source_retry_count=0,
            CFL_retry_count=0,
            state=SimpleNamespace(u=state_u, p=state_u, q=state_u),
            initial=SimpleNamespace(grid=SimpleNamespace(coordinates=coordinates)),
            operator=SimpleNamespace(),
            transaction=SimpleNamespace(
                state=SimpleNamespace(
                    accepted_stage_count=0,
                    last_transaction_serial=-1,
                    last_accepted_time=0.0,
                    first_failed_premise=None,
                    first_failed_transaction_serial=None,
                    first_failed_time=None,
                ),
                causal_state=SimpleNamespace(
                    accepted_time=0.0,
                    accumulated_characteristic_distance=0.0,
                    previous_speed_upper=0.0,
                ),
            ),
            tracers=SimpleNamespace(
                labels=np.zeros(2, dtype=np.float64),
                positions=np.zeros(2, dtype=np.float64),
                proper_times=np.zeros(2, dtype=np.float64),
                event_proper_times=[],
                event_fields=[],
            ),
            temporal_ledger=SimpleNamespace(
                current_macro_step_temporal_retry_count=2,
                last_accepted_time=float.fromhex(contract.ACCEPTED_TIME_HEX),
                accepted_macro_step_count=0,
                cumulative_temporal_retry_count=2,
                last_accepted_macro_step_temporal_retry_count=0,
                accumulated_debit_vector=(0.0,) * 18,
                serialized_temporal_rejections=("{}", "{}"),
            ),
        )
        # transaction.state/causal_state/ledger need asdict; use real dataclasses.
        from dataclasses import dataclass

        @dataclass
        class Monitor:
            accepted_stage_count: int = 0
            last_transaction_serial: int = -1
            last_accepted_time: float = 0.0
            first_failed_premise: str | None = None
            first_failed_transaction_serial: int | None = None
            first_failed_time: float | None = None

        @dataclass
        class Causal:
            accepted_time: float = 0.0
            accumulated_characteristic_distance: float = 0.0
            previous_speed_upper: float = 0.0

        @dataclass
        class Ledger:
            last_accepted_time: float = float.fromhex(contract.ACCEPTED_TIME_HEX)
            accepted_macro_step_count: int = 0
            cumulative_temporal_retry_count: int = 2
            current_macro_step_temporal_retry_count: int = 2
            last_accepted_macro_step_temporal_retry_count: int = 0
            accumulated_debit_vector: tuple = (0.0,) * 18
            serialized_temporal_rejections: tuple = ("{}", "{}")

        member.transaction = SimpleNamespace(state=Monitor(), causal_state=Causal())
        member.temporal_ledger = Ledger()
        templates = {
            "RK4-2049": member,
            "RK4-4097": SimpleNamespace(key="RK4-4097"),
            "RK4-8193": SimpleNamespace(key="RK4-8193"),
            "SSPRK3-4097": SimpleNamespace(key="SSPRK3-4097"),
            "SSPRK3-8193": SimpleNamespace(key="SSPRK3-8193"),
            "SSPRK3-16385": SimpleNamespace(key="SSPRK3-16385"),
        }
        checkpoint = SimpleNamespace(
            sha256=contract.replay_specs()[0]["checkpoint_sha256"],
            members={
                contract.MEMBER_KEY: SimpleNamespace(
                    descriptor_sha256=contract.DESCRIPTOR_SHA256
                )
            },
        )

        class FakeStore:
            def __init__(self, path):
                self.path = path

            def authenticated_checkpoint_at_generation(self, generation):
                self.generation = generation
                return checkpoint

        calls: list[str] = []

        def fake_restore(store, ckpt, restored_member, *, key):
            calls.append(key)
            return restored_member

        state_digest = array_content_sha256(state_u, state_u, state_u)
        coord_digest = array_content_sha256(coordinates)
        with (
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_store.HLT16CampaignStore",
                FakeStore,
            ),
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_runtime.restore_member_with_overlay",
                side_effect=fake_restore,
            ),
            patch.object(contract, "PHYSICAL_STATE_SHA256", state_digest),
            patch.object(contract, "COORDINATES_SHA256", coord_digest),
        ):
            restored = authority.restore_predecessor(ROOT, 3, templates=templates)
        self.assertEqual(calls, [contract.MEMBER_KEY])
        self.assertEqual(restored.retry, 3)
        self.assertEqual(restored.fingerprint["member_key"], contract.MEMBER_KEY)
        self.assertEqual(restored.fingerprint["integrator_id"], PRIMARY_METHOD)
        self.assertEqual(restored.fingerprint["spatial_order"], 4)
        self.assertEqual(
            restored.fingerprint["current_macro_step_temporal_retry_count"], 2
        )
        self.assertNotIn("u", restored.fingerprint)
        self.assertNotIn("endpoints", restored.fingerprint)

        member.spatial_order = 2
        with (
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_store.HLT16CampaignStore",
                FakeStore,
            ),
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_runtime.restore_member_with_overlay",
                side_effect=fake_restore,
            ),
            patch.object(contract, "PHYSICAL_STATE_SHA256", state_digest),
            patch.object(contract, "COORDINATES_SHA256", coord_digest),
            self.assertRaisesRegex(authority.MSEL1AuthorityError, "restored identity"),
        ):
            authority.restore_predecessor(ROOT, 3, templates=templates)

        member.spatial_order = 4
        member.state.u = state_u + 1.0
        with (
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_store.HLT16CampaignStore",
                FakeStore,
            ),
            patch(
                "recursive_horizons.fgc.evolution.hlt16_campaign_runtime.restore_member_with_overlay",
                side_effect=fake_restore,
            ),
            patch.object(contract, "PHYSICAL_STATE_SHA256", state_digest),
            patch.object(contract, "COORDINATES_SHA256", coord_digest),
            self.assertRaisesRegex(authority.MSEL1AuthorityError, "restored identity"),
        ):
            authority.restore_predecessor(ROOT, 3, templates=templates)

    def test_reproducer_default_is_compact_check_and_emit_does_not_authorize(
        self,
    ) -> None:
        script = _load_script()
        with (
            patch.object(
                script.authority,
                "validate_compact_bundle",
                return_value={
                    "config_sha256": "d" * 64,
                    "diagnostic_executed": False,
                    "selection_result_earned": False,
                },
            ) as check,
            patch.object(
                script.authority,
                "authorize_execution",
                side_effect=AssertionError("authorize"),
            ),
            patch.object(
                script.authority,
                "snapshot_store",
                side_effect=AssertionError("store"),
            ),
            patch.object(
                script.authority,
                "observe_environment",
                side_effect=AssertionError("environment"),
            ),
        ):
            code = script.main([])
        self.assertEqual(code, 0)
        check.assert_called_once()
        with (
            patch.object(
                script.authority, "emit_config_bytes", return_value=b"toml\n"
            ) as emit_config,
            patch.object(
                script.authority,
                "authorize_execution",
                side_effect=AssertionError("authorize"),
            ),
            patch.object(
                script.authority,
                "snapshot_store",
                side_effect=AssertionError("store"),
            ),
            patch.object(
                script.authority,
                "validate_compact_bundle",
                side_effect=AssertionError("check"),
            ),
        ):
            code = script.main(["--emit-config"])
        self.assertEqual(code, 0)
        emit_config.assert_called_once()
        with (
            patch.object(
                script.authority,
                "emit_result_object",
                return_value={"artifact_id": contract.ARTIFACT_ID},
            ) as emit_result,
            patch.object(
                script.authority,
                "authorize_execution",
                side_effect=AssertionError("authorize"),
            ),
            patch.object(
                script.authority,
                "snapshot_store",
                side_effect=AssertionError("store"),
            ),
        ):
            code = script.main(["--emit-result"])
        self.assertEqual(code, 0)
        emit_result.assert_called_once()


if __name__ == "__main__":
    unittest.main()
