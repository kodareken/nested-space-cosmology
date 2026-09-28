from __future__ import annotations

import ast
from copy import deepcopy
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
from tempfile import TemporaryDirectory
import tomllib
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from recursive_horizons.fgc.evolution import proto19_pref28_binder as pref28  # noqa: E402


CONFIG_RAW = (ROOT / pref28.CONFIG_PATH).read_bytes()
STORE = ROOT / pref28.STORE_PATH


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
        "artifact_id": pref28.ARTIFACT_ID,
        "project_version": "0.11.0",
        "target_protocol": "FGC-2-SF1-PROTO18",
        "classification": "outcome_neutral_generation9_invalid_terminal_binder",
        "gate_status": "pass",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "artifact_payload": {
            "predecessor": config["predecessor"],
            "campaign": config["campaign"],
            "live_anchor": pref28.expected_boundary(config_raw),
            "scope": config["scope"],
            "claims": config["claims"],
            "nonclaims": config["nonclaims"],
        },
    }


def tree_snapshot(root: Path) -> tuple[tuple[object, ...], ...]:
    """Test-only lexical snapshot which never follows a directory entry."""
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


class PREF28GitAuthorityTests(unittest.TestCase):
    """Immutable authority lookup is isolated from hostile ambient Git state."""

    def test_git_subprocess_strips_ambient_git_state_and_disables_replacements(self) -> None:
        completed = subprocess.CompletedProcess(
            args=(), returncode=0, stdout=b"answer\n", stderr=b""
        )
        hostile = {
            "GIT_DIR": "/attacker/repository",
            "GIT_WORK_TREE": "/attacker/worktree",
            "GIT_CONFIG_GLOBAL": "/attacker/config",
            "GIT_CONFIG_SYSTEM": "/attacker/system-config",
            "GIT_CONFIG_COUNT": "1",
            "GIT_CONFIG_KEY_0": "include.path",
            "GIT_CONFIG_VALUE_0": "/attacker/redirect",
            "GIT_OBJECT_DIRECTORY": "/attacker/objects",
            "GIT_ALTERNATE_OBJECT_DIRECTORIES": "/attacker/alternate",
            "GIT_REPLACE_REF_BASE": "refs/attacker/replace/",
            "GIT_OPTIONAL_LOCKS": "1",
        }
        with (
            patch.dict(os.environ, hostile),
            patch.object(pref28.subprocess, "run", return_value=completed) as run,
        ):
            self.assertEqual(
                pref28._git(Path("/authorized"), "rev-parse", "HEAD"),
                b"answer\n",
            )

        command = run.call_args.args[0]
        self.assertEqual(
            command[:8],
            [
                "git",
                "--no-replace-objects",
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                "rev-parse",
            ],
        )
        self.assertEqual(run.call_args.kwargs["cwd"], Path("/authorized"))
        self.assertIs(run.call_args.kwargs["stdin"], subprocess.DEVNULL)
        environment = run.call_args.kwargs["env"]
        self.assertEqual(environment["GIT_OPTIONAL_LOCKS"], "0")
        self.assertEqual(environment["LC_ALL"], "C")
        self.assertFalse(set(environment) & (set(hostile) - {"GIT_OPTIONAL_LOCKS"}))

    @staticmethod
    def _make_repository(root: Path, content: str) -> str:
        clean_environment = {
            key: value for key, value in os.environ.items() if not key.startswith("GIT_")
        }
        clean_environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"})
        subprocess.run(
            ["git", "init", "-q", str(root)],
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=clean_environment,
        )
        (root / "marker").write_text(content, encoding="utf-8")
        subprocess.run(
            ["git", "-C", str(root), "add", "marker"],
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=clean_environment,
        )
        subprocess.run(
            [
                "git",
                "-C",
                str(root),
                "-c",
                "user.name=PREF28 test",
                "-c",
                "user.email=pref28@example.invalid",
                "commit",
                "-qm",
                "fixture",
            ],
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=clean_environment,
        )
        return (
            subprocess.run(
                ["git", "-C", str(root), "rev-parse", "HEAD"],
                check=True,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=clean_environment,
            )
            .stdout.decode("ascii")
            .strip()
        )

    def test_hostile_git_dir_work_tree_config_and_objects_cannot_redirect_lookup(
        self,
    ) -> None:
        with TemporaryDirectory(prefix="pref28-hostile-git-") as directory:
            base = Path(directory)
            authorized = base / "authorized"
            attacker = base / "attacker"
            authorized_commit = self._make_repository(authorized, "authorized\n")
            self._make_repository(attacker, "attacker\n")
            attacker_config = base / "attacker.gitconfig"
            attacker_config.write_text("this is not valid git config\n", encoding="utf-8")
            hostile = {
                "GIT_DIR": str(attacker / ".git"),
                "GIT_WORK_TREE": str(attacker),
                "GIT_CONFIG_GLOBAL": str(attacker_config),
                "GIT_CONFIG_COUNT": "1",
                "GIT_CONFIG_KEY_0": "include.path",
                "GIT_CONFIG_VALUE_0": str(attacker_config),
                "GIT_OBJECT_DIRECTORY": str(attacker / ".git" / "objects"),
                "GIT_ALTERNATE_OBJECT_DIRECTORIES": str(attacker / ".git" / "objects"),
                "GIT_REPLACE_REF_BASE": "refs/attacker/replace/",
            }
            with patch.dict(os.environ, hostile):
                self.assertEqual(
                    pref28._git(authorized, "rev-parse", "HEAD"),
                    f"{authorized_commit}\n".encode("ascii"),
                )
                self.assertEqual(
                    pref28._git(authorized, "show", "HEAD:marker"),
                    b"authorized\n",
                )


class PREF28CompactResultTests(unittest.TestCase):
    """Compact verification remains valid and store-blind after sealing."""

    def test_expected_boundary_and_compact_validation_are_live_io_free(self) -> None:
        expected = expected_compact_result(CONFIG_RAW)
        encoded = pref28.canonical_result(expected)
        with (
            patch.object(
                pref28,
                "bind_terminal_store",
                side_effect=AssertionError("compact validation reopened the live store"),
            ),
            patch.object(
                pref28,
                "_git",
                side_effect=AssertionError("compact validation queried Git"),
            ),
        ):
            self.assertEqual(pref28.validate_compact_result(CONFIG_RAW, encoded), expected)
            self.assertEqual(
                pref28.expected_boundary(CONFIG_RAW)["generation9"]["checkpoint_sha256"],
                "770728faa3b0a7233e755ceadb2cc10a7109d231989e3824a40004da7ded5ea7",
            )

    def test_compact_claim_and_terminal_promotions_fail_closed_store_free(self) -> None:
        for name, mutate in (
            (
                "root_cause",
                lambda result: result["artifact_payload"]["claims"].__setitem__(
                    "root_cause_localized", True
                ),
            ),
            (
                "physical",
                lambda result: result["artifact_payload"]["claims"].__setitem__(
                    "physical_result_earned", True
                ),
            ),
            (
                "terminal_owner",
                lambda result: result["artifact_payload"]["live_anchor"]["terminal"].__setitem__(
                    "owner", "temporal"
                ),
            ),
        ):
            with self.subTest(name=name):
                result = deepcopy(expected_compact_result(CONFIG_RAW))
                mutate(result)
                with patch.object(
                    pref28,
                    "bind_terminal_store",
                    side_effect=AssertionError("compact validation reopened the live store"),
                ):
                    with self.assertRaises(pref28.PREF28BinderError) as raised:
                        pref28.validate_compact_result(CONFIG_RAW, pref28.canonical_result(result))
                self.assertEqual(raised.exception.stop_id, "PREF28_COMPACT_DRIFT")

    def test_config_authority_hash_state_and_claim_mutations_fail_closed(self) -> None:
        replacements = {
            "authority": (
                "b8fea44384be62918670ccbab5bd1a2d92057224",
                "0" * 40,
            ),
            "terminal": (
                "disposition = \"invalid_terminal\"",
                "disposition = \"scientific_terminal\"",
            ),
            "member_map": (
                "4dac7e098558313bc812832ec60abd8ec9b055fed6c70cc0e34aaaeb4c5541f3",
                "0" * 64,
            ),
            "claim": ("root_cause_localized = false", "root_cause_localized = true"),
        }
        for name, (old, new) in replacements.items():
            with self.subTest(name=name):
                mutated = CONFIG_RAW.replace(old.encode("ascii"), new.encode("ascii"), 1)
                self.assertNotEqual(mutated, CONFIG_RAW)
                with self.assertRaises(pref28.PREF28BinderError) as raised:
                    pref28.expected_boundary(mutated)
                self.assertEqual(raised.exception.stop_id, "PREF28_CONFIG_DRIFT")

    def test_result_schema_and_outcome_neutral_boundary_are_exact(self) -> None:
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
            set(payload),
            {"predecessor", "campaign", "live_anchor", "scope", "claims", "nonclaims"},
        )
        claims = payload["claims"]
        for name in (
            "generation9_terminal_boundary_bound",
            "invalid_runtime_terminal_bound",
            "accepted_member_map_unchanged",
            "terminal_lock_valid",
        ):
            self.assertTrue(claims[name], name)
        for name in (
            "root_cause_localized",
            "common_event_completed",
            "GR0_calibration_completed",
            "candidate_execution_authorized",
            "candidate_runtime_configuration_state_or_outcome_opened",
            "physical_result_earned",
        ):
            self.assertFalse(claims[name], name)
        self.assertIn("does not identify", payload["nonclaims"][0])
        self.assertEqual(payload["live_anchor"]["terminal"]["owner"], "invalid")

    def test_module_import_boundary_excludes_execution_and_diagnostic_runtimes(self) -> None:
        source = Path(pref28.__file__).read_text("utf-8")
        tree = ast.parse(source)
        imports: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module)
        for forbidden in (
            "run_fgc_pro19_event1",
            "hlt16_progression_attempt",
            "hlt16_campaign_recovery",
            "tdg6",
            "tdg7",
            "proposal",
            "candidate",
        ):
            self.assertFalse(any(forbidden in name for name in imports), forbidden)
        self.assertNotIn("_replay_pending_rejection", source)

    def test_tracked_result_validates_compactly_when_present(self) -> None:
        result_path = ROOT / pref28.RESULT_PATH
        if not result_path.exists():
            self.skipTest("compact PREF28 result has not yet been materialized")
        with patch.object(
            pref28,
            "bind_terminal_store",
            side_effect=AssertionError("tracked compact validation reopened the live store"),
        ):
            pref28.validate_compact_result(CONFIG_RAW, result_path.read_bytes())


class PREF28LiveStoreTests(unittest.TestCase):
    """The exact terminal store is inspected but never modified."""

    @classmethod
    def setUpClass(cls) -> None:
        checkpoints = STORE / "checkpoints"
        journals = STORE / "journal"
        if not checkpoints.is_dir() or not journals.is_dir():
            raise unittest.SkipTest("PREF28 ignored live store is unavailable")
        checkpoint_names = sorted(path.name for path in checkpoints.iterdir())
        journal_names = sorted(path.name for path in journals.iterdir())
        if (
            checkpoint_names
            and journal_names
            and (
                checkpoint_names[-1] > f"{9:020d}-{'f' * 64}.json"
                or journal_names[-1] > f"{9:020d}-{'f' * 64}.journal"
            )
        ):
            raise unittest.SkipTest("campaign moved beyond the sealed generation-nine boundary")
        cls.anchor = pref28.bind_terminal_store(CONFIG_RAW, ROOT)

    def clone_store(self, temporary: Path) -> Path:
        repository = temporary / "repository"
        destination = repository / pref28.STORE_PATH
        destination.parent.mkdir(parents=True)
        shutil.copytree(STORE, destination, symlinks=True)
        return repository

    def assert_rejected_unchanged(self, repository: Path) -> pref28.PREF28BinderError:
        store = repository / pref28.STORE_PATH
        before = tree_snapshot(store)
        with (
            patch.object(pref28, "_bind_authority", return_value=dict(pref28.EXPECTED_PREDECESSOR)),
            self.assertRaises(pref28.PREF28BinderError) as raised,
        ):
            pref28.bind_terminal_store(CONFIG_RAW, repository)
        self.assertEqual(tree_snapshot(store), before)
        return raised.exception

    def generation9_path(self, repository: Path) -> Path:
        return next(
            (repository / pref28.STORE_PATH / "checkpoints").glob(
                "00000000000000000009-*.json"
            )
        )

    def test_live_binding_build_and_authority_cross_binding_are_nonmutating(self) -> None:
        before = tree_snapshot(STORE)
        anchor = pref28.bind_terminal_store(CONFIG_RAW, ROOT)
        self.assertEqual(anchor, pref28.expected_boundary(CONFIG_RAW))
        self.assertEqual(anchor["store"]["regular_leaf_count"], 56)
        self.assertEqual(anchor["store"]["generation8_baseline_leaf_count"], 52)
        self.assertEqual(
            anchor["unchanged_state"]["member_map_sha256"],
            "4dac7e098558313bc812832ec60abd8ec9b055fed6c70cc0e34aaaeb4c5541f3",
        )
        self.assertEqual(anchor["terminal"]["exception_type"], "HLT16AttemptError")
        result = pref28.build_pref28_result(CONFIG_RAW, ROOT)
        self.assertEqual(result, expected_compact_result(CONFIG_RAW))
        self.assertEqual(tree_snapshot(STORE), before)

    def test_terminal_member_journal_lock_and_suffix_mutations_fail_closed(self) -> None:
        def member(repository: Path) -> None:
            path = self.generation9_path(repository)
            value = json.loads(path.read_bytes())
            value["members"]["RK4-2049"]["source_current"] = 1
            address = rehash(value, "checkpoint_sha256")
            path.write_bytes(canonical(value))
            path.rename(path.with_name(f"{9:020d}-{address}.json"))

        def journal(repository: Path) -> None:
            directory = repository / pref28.STORE_PATH / "journal"
            path = next(directory.glob("00000000000000000009-*.journal"))
            value = json.loads(path.read_bytes())
            value["payload"]["evidence"]["exception_type"] = "RuntimeError"
            address = rehash(value, "record_sha256")
            path.write_bytes(canonical(value))
            path.rename(path.with_name(f"{9:020d}-{address}.journal"))

        def lock(repository: Path) -> None:
            path = repository / pref28.STORE_PATH / "locks/terminal.lock"
            value = json.loads(path.read_bytes())
            value["checkpoint_sha256"] = "0" * 64
            rehash(value, "lock_sha256")
            path.write_bytes(canonical(value))

        def missing_suffix(repository: Path) -> None:
            (repository / pref28.STORE_PATH / "locks/terminal.lock").unlink()

        def new_payload(repository: Path) -> None:
            (repository / pref28.STORE_PATH / "payloads" / f"{'0' * 64}.npz").write_bytes(b"foreign")

        for name, mutate in (
            ("member", member),
            ("journal", journal),
            ("lock", lock),
            ("missing_suffix", missing_suffix),
            ("new_payload", new_payload),
        ):
            with self.subTest(name=name), TemporaryDirectory(prefix="pref28-mutation-") as directory:
                repository = self.clone_store(Path(directory))
                mutate(repository)
                self.assert_rejected_unchanged(repository)

    def test_symlink_hardlink_and_special_file_fail_closed(self) -> None:
        with TemporaryDirectory(prefix="pref28-symlink-root-") as directory:
            repository = self.clone_store(Path(directory))
            link = Path(directory) / "repository-link"
            os.symlink(repository, link)
            with (
                patch.object(pref28, "_bind_authority", return_value=dict(pref28.EXPECTED_PREDECESSOR)),
                self.assertRaises(pref28.PREF28BinderError) as raised,
            ):
                pref28.bind_terminal_store(CONFIG_RAW, link)
            self.assertEqual(raised.exception.stop_id, "PREF28_PATH_UNSAFE")

        with TemporaryDirectory(prefix="pref28-symlink-leaf-") as directory:
            repository = self.clone_store(Path(directory))
            lock = repository / pref28.STORE_PATH / "locks/terminal.lock"
            relocated = repository / "relocated-terminal.lock"
            lock.rename(relocated)
            os.symlink(relocated, lock)
            self.assert_rejected_unchanged(repository)

        with TemporaryDirectory(prefix="pref28-hardlink-") as directory:
            repository = self.clone_store(Path(directory))
            states = repository / pref28.STORE_PATH / "states"
            source = next(states.glob("*.json"))
            os.link(source, states / f"{'0' * 64}.json")
            error = self.assert_rejected_unchanged(repository)
            self.assertEqual(error.stop_id, "PREF28_PATH_UNSAFE")

        if not hasattr(os, "mkfifo"):
            return
        with TemporaryDirectory(prefix="pref28-fifo-") as directory:
            repository = self.clone_store(Path(directory))
            fifo = repository / pref28.STORE_PATH / "states" / f"{'0' * 64}.json"
            try:
                os.mkfifo(fifo)
            except OSError as error:
                self.skipTest(f"FIFO unsupported by temporary filesystem: {error}")
            special = self.assert_rejected_unchanged(repository)
            self.assertEqual(special.stop_id, "PREF28_PATH_UNSAFE")

    def test_active_writer_staging_and_foreign_fork_fail_closed(self) -> None:
        def active_writer(repository: Path) -> None:
            (repository / pref28.STORE_PATH / "locks/active-write.lock").write_bytes(b"{}")

        def staging(repository: Path) -> None:
            (repository / pref28.STORE_PATH).parent.joinpath(
                ".calibration.hlt16-stage-foreign"
            ).mkdir()

        def fork(repository: Path) -> None:
            source = self.generation9_path(repository)
            shutil.copy2(source, source.with_name(f"{9:020d}-{'0' * 64}.json"))

        for name, mutate in (
            ("active_writer", active_writer),
            ("staging", staging),
            ("fork", fork),
        ):
            with self.subTest(name=name), TemporaryDirectory(prefix="pref28-tree-") as directory:
                repository = self.clone_store(Path(directory))
                mutate(repository)
                self.assert_rejected_unchanged(repository)


if __name__ == "__main__":
    unittest.main()
