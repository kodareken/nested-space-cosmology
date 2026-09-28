"""Adversarial tests for the future-only evidence I/O utility."""

from __future__ import annotations

import ast
from pathlib import Path
import errno
import inspect
import os
import stat
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
HISTORICAL_IMPORT_BASELINE = "5ec2530c213845bcb2031860da08fdff0b30ef2d"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from recursive_horizons import evidence_io  # noqa: E402
from recursive_horizons.evidence_io import (  # noqa: E402
    CanonicalJSONError,
    DEFAULT_MAX_BYTES,
    DeltaEntry,
    GitQueryError,
    InspectDelta,
    InspectTree,
    MAX_JSON_DEPTH,
    PostpublicationUncertainty,
    PrepublicationError,
    PublicationReceipt,
    ReadBlob,
    ResolveCommit,
    TreeEntry,
    UnsafePathError,
    UnsupportedPublication,
    canonical_json_bytes,
    git_read,
    load_canonical_json,
    publish_exclusive_directory,
    publish_exclusive_file,
    read_regular_file,
)


def _init_git_repository(path: Path) -> str:
    env = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    env.update(
        {
            "GIT_AUTHOR_EMAIL": "evidence-io@example.test",
            "GIT_AUTHOR_NAME": "evidence-io-test",
            "GIT_COMMITTER_EMAIL": "evidence-io@example.test",
            "GIT_COMMITTER_NAME": "evidence-io-test",
            "GIT_CONFIG_GLOBAL": os.devnull,
            "GIT_CONFIG_NOSYSTEM": "1",
        }
    )
    subprocess.run(
        ["git", "-c", "init.defaultBranch=main", "init"],
        cwd=path,
        check=True,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    (path / "tracked.txt").write_bytes(b"hello\n")
    subprocess.run(
        ["git", "-c", "core.hooksPath=/dev/null", "add", "tracked.txt"],
        cwd=path,
        check=True,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        ["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "init"],
        cwd=path,
        check=True,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return (
        subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=path,
            env=env,
        )
        .decode("ascii")
        .strip()
    )


def _imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                names.add(node.module)
    return names


class CanonicalJSONTests(unittest.TestCase):
    def test_round_trip_sorted_compact_ascii(self) -> None:
        payload = {"b": [False, True, None], "a": 1, "z": "x"}
        raw = canonical_json_bytes(payload)
        self.assertEqual(raw, b'{"a":1,"b":[false,true,null],"z":"x"}')
        self.assertEqual(load_canonical_json(raw), payload)

    def test_duplicate_keys_are_rejected(self) -> None:
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b'{"a":1,"a":2}')
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b'{"outer":{"k":1,"k":1}}')

    def test_nonfinite_values_are_rejected(self) -> None:
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes(float("nan"))
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes(float("inf"))
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b"[NaN]")
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b"[Infinity]")
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b"[-Infinity]")

    def test_parse_depth_and_byte_bounds_are_typed(self) -> None:
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b"[" * 2000 + b"0" + b"]" * 2000)
        with patch.object(evidence_io, "DEFAULT_MAX_BYTES", 8):
            with self.assertRaises(CanonicalJSONError):
                load_canonical_json(b'"123456789"')
            with self.assertRaises(CanonicalJSONError):
                canonical_json_bytes("123456789")

    def test_unsupported_and_noncanonical_values(self) -> None:
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes({"a": b"bytes"})
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes({"a": {1, 2}})
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes({"a": (1, 2)})
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes({1: "a"})
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes(1 + 2j)
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b'{"a": 1}')
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b'{"a":1}\n')
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json(b'{"b":1,"a":2}')
        with self.assertRaises(CanonicalJSONError):
            load_canonical_json('{"a":1}')  # type: ignore[arg-type]

    def test_subclasses_cycles_and_overdepth_are_rejected(self) -> None:
        class Text(str):
            pass

        class Mapping(dict):
            pass

        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes(Text("a"))
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes(Mapping(a=1))
        cyclic: list[object] = []
        cyclic.append(cyclic)
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes(cyclic)
        nested: object = []
        for _ in range(MAX_JSON_DEPTH + 2):
            nested = [nested]
        with self.assertRaises(CanonicalJSONError):
            canonical_json_bytes(nested)


class PathReadTests(unittest.TestCase):
    def test_reads_a_regular_single_linked_file(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            nested = root / "a" / "b"
            nested.mkdir(parents=True)
            (nested / "leaf.txt").write_bytes(b"payload")
            self.assertEqual(
                read_regular_file(root, "a/b/leaf.txt"),
                b"payload",
            )

    def test_path_traversal_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            (root / "leaf.txt").write_bytes(b"x")
            for relative in (
                "../leaf.txt",
                "a/../../leaf.txt",
                "/etc/passwd",
                "a/./leaf.txt",
                "leaf.txt/",
                "a//leaf.txt",
                "",
                r"a\b",
            ):
                with self.subTest(relative=relative):
                    with self.assertRaises(UnsafePathError):
                        read_regular_file(root, relative)

    def test_leaf_and_ancestor_symlinks_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            target = root / "target"
            target.mkdir()
            (target / "leaf.txt").write_bytes(b"ok")
            os.symlink(target, root / "via")
            os.symlink("leaf.txt", target / "link.txt")
            with self.assertRaises(UnsafePathError):
                read_regular_file(root, "via/leaf.txt")
            with self.assertRaises(UnsafePathError):
                read_regular_file(target, "link.txt")

    def test_hardlinks_special_files_and_directories_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            (root / "one.txt").write_bytes(b"one")
            os.link(root / "one.txt", root / "two.txt")
            (root / "dir").mkdir()
            os.mkfifo(root / "fifo")
            with self.assertRaises(UnsafePathError):
                read_regular_file(root, "two.txt")
            with self.assertRaises(UnsafePathError):
                read_regular_file(root, "dir")
            with self.assertRaises(UnsafePathError):
                read_regular_file(root, "fifo")

    def test_root_symlink_and_byte_bound_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            real = root / "real"
            real.mkdir()
            (real / "leaf.txt").write_bytes(b"abcdef")
            os.symlink(real, root / "link")
            with self.assertRaises(UnsafePathError):
                read_regular_file(root / "link", "leaf.txt")
            with self.assertRaises(UnsafePathError):
                read_regular_file(real, "leaf.txt", max_bytes=3)

    def test_symlinked_ancestor_of_root_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            base = Path(raw).resolve()
            real = base / "real" / "root"
            real.mkdir(parents=True)
            (real / "leaf").write_bytes(b"x")
            (base / "alias").symlink_to(base / "real", target_is_directory=True)
            with self.assertRaises(UnsafePathError):
                read_regular_file(base / "alias/root", "leaf")

    def test_parent_replacement_during_final_read_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            parent = root / "parent"
            parent.mkdir()
            (parent / "leaf").write_bytes(b"old")
            original_read = os.read
            calls = 0

            def replace_on_second_read(descriptor: int, count: int) -> bytes:
                nonlocal calls
                calls += 1
                if calls == 2:
                    parent.rename(root / "moved")
                    parent.mkdir()
                    (parent / "leaf").write_bytes(b"new")
                return original_read(descriptor, count)

            with patch.object(
                evidence_io.os, "read", side_effect=replace_on_second_read
            ):
                with self.assertRaises(UnsafePathError):
                    read_regular_file(root, "parent/leaf")

    def test_leaf_and_ancestor_identity_races_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            nested = root / "keep"
            nested.mkdir()
            (nested / "leaf.txt").write_bytes(b"alpha")
            original = os.stat

            def swap_leaf(
                path: str | os.PathLike[str] | int,
                *args: object,
                dir_fd: int | None = None,
                follow_symlinks: bool = True,
                **kwargs: object,
            ) -> os.stat_result:
                result = original(
                    path,
                    *args,
                    dir_fd=dir_fd,
                    follow_symlinks=follow_symlinks,
                    **kwargs,
                )
                if (
                    follow_symlinks is False
                    and dir_fd is not None
                    and path == "leaf.txt"
                    and stat.S_ISREG(result.st_mode)
                ):
                    os.unlink("leaf.txt", dir_fd=dir_fd)
                    replacement = os.open(
                        "leaf.txt",
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                        0o600,
                        dir_fd=dir_fd,
                    )
                    try:
                        os.write(replacement, b"beta")
                    finally:
                        os.close(replacement)
                return result

            with patch("os.stat", swap_leaf):
                with self.assertRaises(UnsafePathError):
                    read_regular_file(root, "keep/leaf.txt")

            (nested / "other.txt").write_bytes(b"stable")
            swapped = {"done": False}

            def swap_parent(
                path: str | os.PathLike[str] | int,
                *args: object,
                dir_fd: int | None = None,
                follow_symlinks: bool = True,
                **kwargs: object,
            ) -> os.stat_result:
                result = original(
                    path,
                    *args,
                    dir_fd=dir_fd,
                    follow_symlinks=follow_symlinks,
                    **kwargs,
                )
                if (
                    not swapped["done"]
                    and follow_symlinks is False
                    and dir_fd is not None
                    and path == "keep"
                    and stat.S_ISDIR(result.st_mode)
                ):
                    swapped["done"] = True
                    os.rename(nested, root / "moved")
                    replacement = root / "keep"
                    replacement.mkdir()
                    (replacement / "other.txt").write_bytes(b"other")
                return result

            with patch("os.stat", swap_parent):
                with self.assertRaises(UnsafePathError):
                    read_regular_file(root, "keep/other.txt")

    def test_parent_replacement_during_read_does_not_return_stale_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            keep = root / "keep"
            keep.mkdir()
            (keep / "leaf.txt").write_bytes(b"old-bytes")
            original = os.read
            replaced = {"done": False}

            def swap_parent_during_read(fd: int, n: int) -> bytes:
                if not replaced["done"]:
                    replaced["done"] = True
                    os.rename(keep, root / "moved")
                    keep.mkdir()
                    (keep / "leaf.txt").write_bytes(b"new-bytes")
                return original(fd, n)

            with patch("os.read", swap_parent_during_read):
                with self.assertRaises(UnsafePathError):
                    read_regular_file(root, "keep/leaf.txt")
            self.assertEqual((root / "keep" / "leaf.txt").read_bytes(), b"new-bytes")
            self.assertEqual((root / "moved" / "leaf.txt").read_bytes(), b"old-bytes")

    def test_root_ancestor_traversal_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            nested = root / "keep"
            nested.mkdir()
            (nested / "leaf.txt").write_bytes(b"secret")
            with self.assertRaises(UnsafePathError):
                read_regular_file(nested / "..", "keep/leaf.txt")
            with self.assertRaises(UnsafePathError):
                read_regular_file(Path("keep") / ".." / "keep", "leaf.txt")


class GitQueryTests(unittest.TestCase):
    def test_malformed_tree_and_delta_identities_are_rejected(self) -> None:
        for raw in (
            b"100644 blob HEAD\tfile\0",
            b"100644 commit " + b"a" * 40 + b"\tfile\0",
            b"999999 blob " + b"a" * 40 + b"\tfile\0",
        ):
            with self.assertRaises(GitQueryError):
                evidence_io._parse_ls_tree(raw)
        for raw in (b"R100\0one\0two\0", b"MM\0file\0", b"PASS\0file\0"):
            with self.assertRaises(GitQueryError):
                evidence_io._parse_name_status(raw)

    def test_closed_read_only_operations_and_hostile_environment(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            repo = Path(raw).resolve() / "repo"
            decoy = Path(raw).resolve() / "decoy"
            repo.mkdir()
            decoy.mkdir()
            commit = _init_git_repository(repo)
            (repo / "tracked.txt").write_bytes(b"changed\n")
            env = {
                key: value
                for key, value in os.environ.items()
                if not key.startswith("GIT_")
            }
            env.update(
                {
                    "GIT_AUTHOR_EMAIL": "evidence-io@example.test",
                    "GIT_AUTHOR_NAME": "evidence-io-test",
                    "GIT_COMMITTER_EMAIL": "evidence-io@example.test",
                    "GIT_COMMITTER_NAME": "evidence-io-test",
                    "GIT_CONFIG_GLOBAL": os.devnull,
                    "GIT_CONFIG_NOSYSTEM": "1",
                }
            )
            subprocess.run(
                ["git", "-c", "core.hooksPath=/dev/null", "add", "tracked.txt"],
                cwd=repo,
                check=True,
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                ["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "second"],
                cwd=repo,
                check=True,
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            second = (
                subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, env=env)
                .decode("ascii")
                .strip()
            )
            evil_config = Path(raw).resolve() / "evil.gitconfig"
            pwned = Path(raw).resolve() / "pwned"
            evil_config.write_text(
                "[alias]\n"
                "    rev-parse = \"!touch '{pwned}' && false\"\n"
                "    cat-file = \"!touch '{pwned}' && false\"\n"
                "    ls-tree = \"!touch '{pwned}' && false\"\n"
                "    diff-tree = \"!touch '{pwned}' && false\"\n"
                "[core]\n"
                "    sshCommand = \"touch '{pwned}'\"\n".format(pwned=pwned)
            )
            subprocess.run(
                ["git", "init", "--bare"],
                cwd=decoy,
                check=True,
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            hostile = {
                "GIT_DIR": str(decoy),
                "GIT_WORK_TREE": str(decoy),
                "GIT_CONFIG_GLOBAL": str(evil_config),
                "GIT_CONFIG_SYSTEM": str(evil_config),
                "GIT_EXEC_PATH": str(Path(raw).resolve() / "bin"),
                "GIT_SSH_COMMAND": f"touch {pwned}",
                "GIT_ALIAS_REV_PARSE": "status",
            }
            with patch.dict(os.environ, hostile, clear=False):
                resolved = git_read(repo, ResolveCommit("HEAD"))
                self.assertEqual(resolved, second)
                blob = git_read(repo, ReadBlob(commit=commit, path="tracked.txt"))
                self.assertEqual(blob, b"hello\n")
                tree = git_read(repo, InspectTree(commit=commit))
                self.assertEqual(
                    tree,
                    (
                        TreeEntry(
                            mode="100644",
                            kind="blob",
                            object_id=tree[0].object_id,
                            path="tracked.txt",
                        ),
                    ),
                )
                delta = git_read(repo, InspectDelta(commit=second, against=commit))
                self.assertEqual(delta, (DeltaEntry(status="M", path="tracked.txt"),))
            self.assertFalse(pwned.exists())

    def test_read_only_command_constraints_reject_mutation_and_flags(self) -> None:
        self.assertNotIn("run_git", evidence_io.__all__)
        self.assertFalse(hasattr(evidence_io, "run_git"))
        self.assertFalse(
            inspect.signature(git_read).parameters["operation"].kind
            is inspect.Parameter.VAR_POSITIONAL
        )
        with tempfile.TemporaryDirectory() as raw:
            repo = Path(raw).resolve()
            _init_git_repository(repo)
            for revision in (
                "--output=/tmp/pwned",
                "-C",
                "HEAD~1",
                "main",
                "refs/heads/main",
            ):
                with self.subTest(revision=revision):
                    with self.assertRaises(GitQueryError):
                        git_read(repo, ResolveCommit(revision))
            with self.assertRaises(GitQueryError):
                git_read(repo, ReadBlob(commit="HEAD", path="../tracked.txt"))
            with self.assertRaises(GitQueryError):
                git_read(repo, "commit -am pwned")  # type: ignore[arg-type]
            with self.assertRaises(GitQueryError):
                git_read(repo, ReadBlob(commit="HEAD", path="missing.txt"))
        source = inspect.getsource(evidence_io.git_read)
        self.assertNotIn("shell=True", inspect.getsource(evidence_io))
        self.assertIn("closed read-only set", source)

    def test_minimal_environment_ignores_path_and_loader_variables(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            repo = Path(raw).resolve() / "repo"
            bindir = Path(raw).resolve() / "bin"
            repo.mkdir()
            bindir.mkdir()
            _init_git_repository(repo)
            pwned = Path(raw).resolve() / "pwned"
            fake = bindir / "git"
            fake.write_text("#!/bin/sh\ntouch '{pwned}'\nexit 1\n".format(pwned=pwned))
            fake.chmod(0o755)
            hostile = {
                "PATH": str(bindir),
                "LD_PRELOAD": str(fake),
                "LD_LIBRARY_PATH": str(bindir),
                "DYLD_INSERT_LIBRARIES": str(fake),
                "DYLD_LIBRARY_PATH": str(bindir),
            }
            with patch.dict(os.environ, hostile, clear=False):
                env = evidence_io._git_environment()
                self.assertEqual(env["PATH"], evidence_io._GIT_MINIMAL_PATH)
                self.assertNotIn("LD_PRELOAD", env)
                self.assertNotIn("LD_LIBRARY_PATH", env)
                self.assertNotIn("DYLD_INSERT_LIBRARIES", env)
                self.assertNotIn("DYLD_LIBRARY_PATH", env)
                executable = evidence_io._git_executable()
                self.assertTrue(os.path.isabs(executable))
                self.assertIn(executable, evidence_io._GIT_EXECUTABLE_CANDIDATES)
                resolved = git_read(repo, ResolveCommit("HEAD"))
            self.assertEqual(len(resolved), 40)
            self.assertFalse(pwned.exists())

    def test_git_output_is_capped_without_unbounded_allocation(self) -> None:
        proc = subprocess.Popen(
            [
                sys.executable,
                "-c",
                "import sys; sys.stdout.buffer.write(b'x' * (16 * 1024 * 1024 + 1))",
            ],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        with self.assertRaises(GitQueryError):
            evidence_io._communicate_bounded(
                proc, max_bytes=DEFAULT_MAX_BYTES, timeout=10
            )


class PublicationTests(unittest.TestCase):
    def test_exclusive_file_and_directory_publication(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            receipt = publish_exclusive_file(root, "out.bin", b"abc")
            self.assertIsInstance(receipt, PublicationReceipt)
            self.assertEqual(receipt.kind, "file")
            self.assertEqual((root / "out.bin").read_bytes(), b"abc")
            nested = publish_exclusive_directory(
                root,
                "terminal",
                {"a.txt": b"A", "sub/b.txt": b"B"},
            )
            self.assertEqual(nested.kind, "directory")
            self.assertEqual((root / "terminal" / "a.txt").read_bytes(), b"A")
            self.assertEqual((root / "terminal" / "sub" / "b.txt").read_bytes(), b"B")

    def test_existing_and_partial_destinations_are_not_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            (root / "exists.bin").write_bytes(b"old")
            with self.assertRaises(PrepublicationError) as existing:
                publish_exclusive_file(root, "exists.bin", b"new")
            self.assertFalse(existing.exception.published)
            self.assertEqual((root / "exists.bin").read_bytes(), b"old")
            (root / "partial.bin").write_bytes(b"")
            with self.assertRaises(PrepublicationError):
                publish_exclusive_file(root, "partial.bin", b"payload")
            self.assertEqual((root / "partial.bin").read_bytes(), b"")
            (root / "dir-exists").mkdir()
            with self.assertRaises(PrepublicationError):
                publish_exclusive_directory(root, "dir-exists", {"a.txt": b"x"})
            self.assertFalse((root / "dir-exists" / "a.txt").exists())

    def test_racing_destination_fails_closed_without_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()

            def race(cut: str) -> None:
                if cut == "before_publish":
                    (root / "raced.bin").write_bytes(b"other")

            with self.assertRaises(PrepublicationError) as failed:
                publish_exclusive_file(root, "raced.bin", b"payload", _fault_hook=race)
            self.assertFalse(failed.exception.published)
            self.assertEqual((root / "raced.bin").read_bytes(), b"other")
            leftovers = list(root.glob(".raced.bin.evidence-io-stage-*"))
            self.assertEqual(leftovers, [])

    def test_postpublication_uncertainty_does_not_delete_or_retry(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()

            def fail_after(cut: str) -> None:
                if cut == "after_publish":
                    raise OSError(errno.EIO, "injected durability failure")

            with self.assertRaises(PostpublicationUncertainty) as failed:
                publish_exclusive_file(
                    root, "done.bin", b"kept", _fault_hook=fail_after
                )
            self.assertTrue(failed.exception.published)
            self.assertEqual(failed.exception.relative_path, "done.bin")
            self.assertEqual((root / "done.bin").read_bytes(), b"kept")
            with self.assertRaises(PrepublicationError):
                publish_exclusive_file(root, "done.bin", b"kept")
            self.assertEqual((root / "done.bin").read_bytes(), b"kept")

    def test_concurrent_publishers_have_one_winner(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            barrier = threading.Barrier(2)
            outcomes: list[object] = []

            def worker() -> None:
                barrier.wait()
                try:
                    outcomes.append(publish_exclusive_file(root, "shared.bin", b"same"))
                except Exception as error:
                    outcomes.append(error)

            threads = [threading.Thread(target=worker) for _ in range(2)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
            successes = [
                item for item in outcomes if isinstance(item, PublicationReceipt)
            ]
            failures = [
                item for item in outcomes if isinstance(item, PrepublicationError)
            ]
            self.assertEqual(len(successes), 1)
            self.assertEqual(len(failures), 1)
            self.assertEqual((root / "shared.bin").read_bytes(), b"same")

    def test_unsupported_no_replace_semantics_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            with patch.object(evidence_io, "_RENAME_FUNC", None):
                with self.assertRaises(UnsupportedPublication):
                    publish_exclusive_file(root, "nope.bin", b"x")
            self.assertFalse((root / "nope.bin").exists())
            with self.assertRaises(PrepublicationError):
                publish_exclusive_directory(root, "overlap", {"a": b"1", "a/b": b"2"})

    def test_foreign_stage_collision_is_retained(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            (root / "foreign-stage").write_bytes(b"foreign")
            with patch.object(evidence_io, "_stage_name", return_value="foreign-stage"):
                with self.assertRaises(PrepublicationError) as failed:
                    publish_exclusive_file(root, "out.bin", b"payload")
            self.assertFalse(failed.exception.published)
            self.assertTrue(failed.exception.staging_retained)
            self.assertEqual((root / "foreign-stage").read_bytes(), b"foreign")
            self.assertFalse((root / "out.bin").exists())

            nested = root / "foreign-dir"
            nested.mkdir()
            (nested / "inside.txt").write_bytes(b"keep-me")
            with patch.object(evidence_io, "_stage_name", return_value="foreign-dir"):
                with self.assertRaises(PrepublicationError) as failed_dir:
                    publish_exclusive_directory(root, "terminal", {"a.txt": b"A"})
            self.assertTrue(failed_dir.exception.staging_retained)
            self.assertEqual((nested / "inside.txt").read_bytes(), b"keep-me")
            self.assertFalse((root / "terminal").exists())

    def test_stage_tamper_before_publish_is_not_a_success(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()

            def tamper(cut: str) -> None:
                if cut != "before_publish":
                    return
                stages = list(root.glob(".out.bin.evidence-io-stage-*"))
                self.assertEqual(len(stages), 1)
                stages[0].write_bytes(b"tampered")

            with self.assertRaises(PrepublicationError) as failed:
                publish_exclusive_file(root, "out.bin", b"expected", _fault_hook=tamper)
            self.assertFalse(failed.exception.published)
            self.assertTrue(failed.exception.staging_retained)
            self.assertFalse((root / "out.bin").exists())
            leftovers = list(root.glob(".out.bin.evidence-io-stage-*"))
            self.assertEqual(len(leftovers), 1)
            self.assertEqual(leftovers[0].read_bytes(), b"tampered")

    def test_stage_inode_swap_is_retained_and_not_published(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()

            def swap(cut: str) -> None:
                if cut != "before_publish":
                    return
                stages = list(root.glob(".out.bin.evidence-io-stage-*"))
                self.assertEqual(len(stages), 1)
                stages[0].unlink()
                stages[0].write_bytes(b"substituted")

            with self.assertRaises(PrepublicationError) as failed:
                publish_exclusive_file(root, "out.bin", b"expected", _fault_hook=swap)
            self.assertFalse(failed.exception.published)
            self.assertTrue(failed.exception.staging_retained)
            self.assertFalse((root / "out.bin").exists())
            leftovers = list(root.glob(".out.bin.evidence-io-stage-*"))
            self.assertEqual(len(leftovers), 1)
            self.assertEqual(leftovers[0].read_bytes(), b"substituted")

    def test_directory_extra_mutated_and_replaced_leaves_are_not_success(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()

            def extra(cut: str) -> None:
                if cut != "before_publish":
                    return
                stages = list(root.glob(".term.evidence-io-stage-*"))
                self.assertEqual(len(stages), 1)
                (stages[0] / "extra.txt").write_bytes(b"nope")

            with self.assertRaises(PrepublicationError) as extra_failure:
                publish_exclusive_directory(
                    root,
                    "term",
                    {"a.txt": b"A"},
                    _fault_hook=extra,
                )
            self.assertFalse((root / "term").exists())
            self.assertTrue(extra_failure.exception.staging_retained)
            retained = next(root.glob(".term.evidence-io-stage-*"))
            self.assertEqual((retained / "extra.txt").read_bytes(), b"nope")

            def mutate(cut: str) -> None:
                if cut != "before_publish":
                    return
                stages = list(root.glob(".term2.evidence-io-stage-*"))
                self.assertEqual(len(stages), 1)
                (stages[0] / "a.txt").write_bytes(b"mutated")

            with self.assertRaises(PrepublicationError):
                publish_exclusive_directory(
                    root,
                    "term2",
                    {"a.txt": b"A"},
                    _fault_hook=mutate,
                )
            self.assertFalse((root / "term2").exists())

            def replace_leaf(cut: str) -> None:
                if cut != "before_publish":
                    return
                stages = list(root.glob(".term3.evidence-io-stage-*"))
                self.assertEqual(len(stages), 1)
                (stages[0] / "a.txt").unlink()
                (stages[0] / "a.txt").write_bytes(b"replaced")

            with self.assertRaises(PrepublicationError):
                publish_exclusive_directory(
                    root,
                    "term3",
                    {"a.txt": b"A"},
                    _fault_hook=replace_leaf,
                )
            self.assertFalse((root / "term3").exists())

    def test_unexpected_empty_directory_is_retained_and_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()

            def extra_empty(cut: str) -> None:
                if cut == "before_publish":
                    stage = next(root.glob(".terminal.evidence-io-stage-*"))
                    (stage / "unexpected").mkdir()

            with self.assertRaises(PrepublicationError) as failed:
                publish_exclusive_directory(
                    root, "terminal", {"result.json": b"{}"}, _fault_hook=extra_empty
                )
            self.assertTrue(failed.exception.staging_retained)
            self.assertFalse((root / "terminal").exists())
            stage = next(root.glob(".terminal.evidence-io-stage-*"))
            self.assertTrue((stage / "unexpected").is_dir())

    def test_every_nested_directory_is_fsynced_before_publication(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            synced: set[tuple[int, int]] = set()
            original_fsync = os.fsync

            def record_fsync(descriptor: int) -> None:
                metadata = os.fstat(descriptor)
                if stat.S_ISDIR(metadata.st_mode):
                    synced.add((metadata.st_dev, metadata.st_ino))
                original_fsync(descriptor)

            with patch.object(evidence_io.os, "fsync", side_effect=record_fsync):
                publish_exclusive_directory(root, "terminal", {"a/b/value": b"x"})
            for directory in (
                root / "terminal",
                root / "terminal/a",
                root / "terminal/a/b",
            ):
                metadata = directory.stat()
                self.assertIn((metadata.st_dev, metadata.st_ino), synced)

    def test_payload_type_and_mutation_cannot_change_the_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            with self.assertRaises(PrepublicationError):
                publish_exclusive_file(root, "x.bin", 1024)  # type: ignore[arg-type]
            with self.assertRaises(PrepublicationError):
                publish_exclusive_directory(root, "d", {"a.txt": 10**6})  # type: ignore[dict-item]
            self.assertFalse((root / "x.bin").exists())
            self.assertFalse((root / "d").exists())

            mutable = bytearray(b"first")
            files: dict[str, bytes | bytearray] = {"a.txt": mutable}

            def mutate(cut: str) -> None:
                if cut == "before_write":
                    mutable[:] = b"XXXXX"
                    files["b.txt"] = b"added"

            receipt = publish_exclusive_directory(
                root,
                "snap",
                files,
                _fault_hook=mutate,  # type: ignore[arg-type]
            )
            self.assertEqual((root / "snap" / "a.txt").read_bytes(), b"first")
            self.assertFalse((root / "snap" / "b.txt").exists())
            self.assertEqual(
                receipt.payload_sha256,
                evidence_io._inventory_digest({"a.txt": b"first"}),
            )

    def test_post_rename_dest_mutation_is_uncertainty_not_success(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()

            def mutate_dest(cut: str) -> None:
                if cut == "after_publish":
                    (root / "done.bin").write_bytes(b"changed")

            with self.assertRaises(PostpublicationUncertainty) as failed:
                publish_exclusive_file(
                    root, "done.bin", b"kept", _fault_hook=mutate_dest
                )
            self.assertTrue(failed.exception.published)
            self.assertEqual((root / "done.bin").read_bytes(), b"changed")


class AuthorityGitQueryTests(unittest.TestCase):
    def test_worktree_changes_are_observed_without_index_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            repo = Path(raw).resolve()
            _init_git_repository(repo)
            self.assertTrue(git_read(repo, evidence_io.InspectWorktree()).clean)
            (repo / "tracked.txt").write_bytes(b"changed\n")
            (repo / "new.txt").write_bytes(b"new\n")
            index_before = (repo / ".git/index").read_bytes()
            observed = git_read(repo, evidence_io.InspectWorktree())
            self.assertFalse(observed.clean)
            self.assertEqual(
                {
                    (item.index_status, item.worktree_status, item.path)
                    for item in observed.changes
                },
                {(" ", "M", "tracked.txt"), ("?", "?", "new.txt")},
            )
            self.assertEqual(observed.hidden_index_paths, ())
            self.assertEqual((repo / ".git/index").read_bytes(), index_before)

    def test_assume_unchanged_and_skip_worktree_cannot_appear_clean(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            repo = Path(raw).resolve()
            _init_git_repository(repo)
            for flag, undo in (
                ("--assume-unchanged", "--no-assume-unchanged"),
                ("--skip-worktree", "--no-skip-worktree"),
            ):
                subprocess.run(
                    ["git", "update-index", flag, "tracked.txt"],
                    cwd=repo,
                    check=True,
                    env=evidence_io._git_environment(),
                )
                observed = git_read(repo, evidence_io.InspectWorktree())
                self.assertFalse(observed.clean)
                self.assertEqual(observed.hidden_index_paths, ("tracked.txt",))
                subprocess.run(
                    ["git", "update-index", undo, "tracked.txt"],
                    cwd=repo,
                    check=True,
                    env=evidence_io._git_environment(),
                )
            self.assertTrue(git_read(repo, evidence_io.InspectWorktree()).clean)

    def test_parent_query_distinguishes_root_successor_and_merge(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            repo = Path(raw).resolve()
            first = _init_git_repository(repo)
            self.assertEqual(git_read(repo, evidence_io.CommitParents(first)), ())
            env = evidence_io._git_environment() | {
                "GIT_AUTHOR_NAME": "evidence-io-test",
                "GIT_AUTHOR_EMAIL": "evidence-io@example.test",
                "GIT_COMMITTER_NAME": "evidence-io-test",
                "GIT_COMMITTER_EMAIL": "evidence-io@example.test",
            }
            tree = (
                subprocess.check_output(
                    ["git", "rev-parse", "HEAD^{tree}"], cwd=repo, env=env
                )
                .decode()
                .strip()
            )

            def commit(parents):
                args = ["git", "commit-tree", tree]
                for parent in parents:
                    args.extend(("-p", parent))
                return (
                    subprocess.check_output(
                        args, input=b"synthetic\n", cwd=repo, env=env
                    )
                    .decode()
                    .strip()
                )

            second = commit((first,))
            merge = commit((first, second))
            self.assertEqual(
                git_read(repo, evidence_io.CommitParents(second)), (first,)
            )
            self.assertEqual(
                git_read(repo, evidence_io.CommitParents(merge)), (first, second)
            )

    def test_new_git_parsers_reject_ambiguous_or_malformed_streams(self) -> None:
        for status, flags in (
            (b"?? file", b""),
            (b"R  new\0old\0", b""),
            (b"XX file\0", b""),
            (b"?? file\0?? file\0", b""),
            (b"", b"H file"),
            (b"", b"Z file\0"),
        ):
            with self.assertRaises(GitQueryError):
                evidence_io._parse_worktree(status, flags)
        tree = b"tree " + b"a" * 40 + b"\n"
        for raw in (
            b"not a commit",
            tree + b"parent broken\n\nmessage",
            tree + b"author x\nparent " + b"b" * 40 + b"\n\nmessage",
        ):
            with self.assertRaises(GitQueryError):
                evidence_io._parse_commit_parents(raw)

    def test_git_ceiling_is_separate_from_publication_ceiling(self) -> None:
        self.assertEqual(evidence_io.MAX_DIRECTORY_ENTRIES, 1024)
        self.assertEqual(evidence_io.MAX_GIT_ENTRIES, 16384)
        listing = b"".join(
            b"100644 blob " + b"a" * 40 + b"\tf" + str(index).encode() + b".txt\0"
            for index in range(1322)
        )
        self.assertEqual(len(evidence_io._parse_ls_tree(listing)), 1322)
        with patch.object(evidence_io, "MAX_GIT_ENTRIES", 2):
            with self.assertRaises(GitQueryError):
                evidence_io._parse_ls_tree(listing)


class ImportIsolationTests(unittest.TestCase):
    def test_utility_is_stdlib_only_and_not_a_scientific_owner(self) -> None:
        path = SRC / "recursive_horizons" / "evidence_io.py"
        imports = _imported_modules(path)
        stdlib = sys.stdlib_module_names
        for name in sorted(imports):
            root_name = name.split(".", 1)[0]
            self.assertIn(root_name, stdlib)
        joined = "\n".join(sorted(imports))
        self.assertNotIn("recursive_horizons", joined)
        self.assertNotIn("numpy", joined)
        source = path.read_text(encoding="utf-8")
        self.assertIn("This module is infrastructure", source)
        self.assertIn("does not classify scientific results", source)

    def test_historical_binders_and_runners_do_not_import_the_utility(self) -> None:
        listed = subprocess.check_output(
            [
                "git",
                "--no-replace-objects",
                "--no-optional-locks",
                "-C",
                str(ROOT),
                "ls-tree",
                "-r",
                "--name-only",
                "-z",
                HISTORICAL_IMPORT_BASELINE,
                "--",
                "src",
                "scripts",
            ],
            env=evidence_io._git_environment(),
        )
        historical_paths = {
            item.decode("utf-8")
            for item in listed.split(b"\0")
            if item and item.endswith(b".py")
        }
        self.assertIn(
            "src/recursive_horizons/fgc/evolution/tdg10_qa2_pref1_binder.py",
            historical_paths,
        )
        self.assertNotIn("src/recursive_horizons/evidence_io.py", historical_paths)
        offenders: list[str] = []
        for relative in sorted(historical_paths):
            path = ROOT / relative
            text = path.read_text(encoding="utf-8")
            if "evidence_io" not in text:
                continue
            for name in _imported_modules(path):
                if name == "recursive_horizons.evidence_io" or name.endswith(
                    ".evidence_io"
                ):
                    offenders.append(relative)
                    break
            tree = ast.parse(text, filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom) and node.module in {
                    "recursive_horizons",
                    "recursive_horizons.evidence_io",
                }:
                    if any(alias.name == "evidence_io" for alias in node.names):
                        offenders.append(relative)
        self.assertEqual(offenders, [])


if __name__ == "__main__":
    unittest.main()
