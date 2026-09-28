from __future__ import annotations

import ast
from contextlib import redirect_stderr
from hashlib import sha256
from io import StringIO
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

from scripts import bootstrap_fgc_pro19_sid3 as bootstrap


ROOT = Path(__file__).resolve().parents[1]


class _Guard:
    def __init__(self, verified: object, authority: object, runner: object) -> None:
        self.verified = verified
        self.authority = authority
        self.runner = runner
        self.audits = 0

    def __enter__(self) -> "_Guard":
        return self

    def __exit__(self, *args: object) -> bool:
        return False

    def import_module(self, name: str) -> object:
        if name == bootstrap.AUTHORITY_MODULE:
            return self.authority
        if name == bootstrap.RUNNER_MODULE:
            return self.runner
        raise AssertionError(name)

    def audit(self) -> None:
        self.audits += 1


class SID3BootstrapTests(unittest.TestCase):
    def test_bootstrap_has_only_standard_library_imports_before_guard(self) -> None:
        source = (ROOT / bootstrap.BOOTSTRAP_PATH).read_text("utf-8")
        tree = ast.parse(source)
        roots: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                roots.update(alias.name.split(".", 1)[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                roots.add(node.module.split(".", 1)[0])
        self.assertTrue(roots - {"__future__"} <= sys.stdlib_module_names)
        self.assertNotIn("sys.path.insert", source)
        self.assertNotIn("sys.path.append", source)

    def test_canonical_json_rejects_duplicates_and_format_drift(self) -> None:
        self.assertEqual(bootstrap._json(b"{}\n", "x"), {})
        with self.assertRaisesRegex(bootstrap.SID3BootstrapError, "noncanonical"):
            bootstrap._json(b"{ }\n", "x")
        with self.assertRaisesRegex(bootstrap.SID3BootstrapError, "malformed"):
            bootstrap._json(b'{"x":1,"x":2}\n', "x")

    def test_nofollow_reader_rejects_leaf_ancestor_and_hardlink(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "real").mkdir()
            (root / "real" / "leaf").write_bytes(b"value")
            self.assertEqual(
                bootstrap._read_nofollow(root, "real/leaf", "leaf"), b"value"
            )
            os.symlink(root / "real" / "leaf", root / "link")
            with self.assertRaises(bootstrap.SID3BootstrapError):
                bootstrap._read_nofollow(root, "link", "leaf symlink")
            os.symlink(root / "real", root / "directory-link")
            with self.assertRaises(bootstrap.SID3BootstrapError):
                bootstrap._read_nofollow(
                    root, "directory-link/leaf", "ancestor symlink"
                )
            os.link(root / "real" / "leaf", root / "hardlink")
            with self.assertRaises(bootstrap.SID3BootstrapError):
                bootstrap._read_nofollow(root, "hardlink", "hardlink")

    def test_git_subprocess_strips_every_ambient_git_variable(self) -> None:
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
            "GIT_OPTIONAL_LOCKS": "1",
        }
        with (
            patch.dict(os.environ, hostile),
            patch.object(bootstrap.subprocess, "run", return_value=completed) as run,
        ):
            self.assertEqual(
                bootstrap._git(Path("/authorized"), ("rev-parse", "HEAD"), "git"),
                b"answer\n",
            )
        command = run.call_args.args[0]
        self.assertEqual(
            command[:8],
            (
                "git",
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                "-C",
                "/authorized",
            ),
        )
        environment = run.call_args.kwargs["env"]
        self.assertEqual(environment["GIT_OPTIONAL_LOCKS"], "0")
        self.assertEqual(environment["LC_ALL"], "C")
        self.assertFalse(set(environment) & (set(hostile) - {"GIT_OPTIONAL_LOCKS"}))

    @staticmethod
    def _make_repository(root: Path, content: str) -> str:
        clean_environment = {
            key: value
            for key, value in os.environ.items()
            if not key.startswith("GIT_")
        }
        clean_environment["LC_ALL"] = "C"
        subprocess.run(
            ("git", "init", "-q", str(root)),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=clean_environment,
        )
        (root / "marker").write_text(content, encoding="utf-8")
        subprocess.run(
            ("git", "-C", str(root), "add", "marker"),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=clean_environment,
        )
        subprocess.run(
            (
                "git",
                "-C",
                str(root),
                "-c",
                "user.name=SID3 Test",
                "-c",
                "user.email=sid3@example.invalid",
                "commit",
                "-qm",
                "fixture",
            ),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=clean_environment,
        )
        return (
            subprocess.run(
                ("git", "-C", str(root), "rev-parse", "HEAD"),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=clean_environment,
            )
            .stdout.decode("ascii")
            .strip()
        )

    def test_git_dir_work_tree_and_config_environment_cannot_redirect_commit(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            authorized = base / "authorized"
            attacker = base / "attacker"
            authorized_commit = self._make_repository(authorized, "authorized\n")
            self._make_repository(attacker, "attacker\n")
            attacker_config = base / "attacker.gitconfig"
            attacker_config.write_text(
                "this is not valid git config\n", encoding="utf-8"
            )
            hostile = {
                "GIT_DIR": str(attacker / ".git"),
                "GIT_WORK_TREE": str(attacker),
                "GIT_CONFIG_GLOBAL": str(attacker_config),
                "GIT_CONFIG_COUNT": "1",
                "GIT_CONFIG_KEY_0": "include.path",
                "GIT_CONFIG_VALUE_0": str(attacker_config),
            }
            with patch.dict(os.environ, hostile):
                self.assertEqual(
                    bootstrap._resolve_commit(
                        authorized, authorized_commit, "authorized commit"
                    ),
                    authorized_commit,
                )
                self.assertEqual(
                    bootstrap._git_show(
                        authorized, authorized_commit, "marker", "authorized marker"
                    ),
                    b"authorized\n",
                )

    def test_execution_closure_is_loaded_from_supplied_commit_bytes(self) -> None:
        source = b"def parse_execution_closure_record(x): return x\n\ndef verify_execution_closure(*a, **k): return None\n\nclass ImportOriginGuard: pass\n"
        with patch.object(bootstrap, "_git_show", return_value=source) as show:
            module, observed = bootstrap._load_closure_module(Path("/repo"), "a" * 40)
        self.assertEqual(observed, source)
        self.assertTrue(module.__file__.startswith("/__fgc_committed__/"))
        self.assertTrue(callable(module.parse_execution_closure_record))
        self.assertEqual(show.call_args.args[1], "a" * 40)

    def test_stdin_pseudo_origin_is_rebound_outside_repository(self) -> None:
        main = SimpleNamespace(__file__="<stdin>")
        with patch.dict(sys.modules, {"__main__": main}):
            bootstrap._mark_committed_stdin_origin("a" * 40)
        self.assertEqual(
            main.__file__,
            f"/__fgc_committed__/{'a' * 40}/{bootstrap.BOOTSTRAP_PATH}",
        )

    def _execute_fixture(
        self, mode: str, *, store_root: Path | None = None
    ) -> tuple[dict[str, object], Mock, _Guard]:
        closure_source = b"closure-source"
        bootstrap_source = b"bootstrap-source"
        pins = (
            SimpleNamespace(
                path=bootstrap.CLOSURE_MODULE_PATH,
                kind="python",
                sha256=sha256(closure_source).hexdigest(),
            ),
            SimpleNamespace(
                path=bootstrap.BOOTSTRAP_PATH,
                kind="python",
                sha256=sha256(bootstrap_source).hexdigest(),
            ),
        )
        record = SimpleNamespace(implementation_commit="a" * 40, files=pins)
        verified = SimpleNamespace(
            files=tuple(
                SimpleNamespace(pin=pin, content=pin.path.encode("ascii"))
                for pin in pins
            )
        )
        receipt = SimpleNamespace(
            authority_commit="c" * 40,
            implementation_commit="a" * 40,
            store_path=bootstrap.DEFAULT_STORE_PATH,
        )
        authority_api = SimpleNamespace(
            authorize_resume_image=Mock(return_value=receipt)
        )
        runner = SimpleNamespace(
            sid3_resume_preflight=Mock(return_value={"safe_to_resume": True}),
            sid3_resume=Mock(return_value={"state": "advanced"}),
        )
        guard = _Guard(verified, authority_api, runner)
        closure_module = SimpleNamespace(
            parse_execution_closure_record=Mock(return_value=record),
            verify_execution_closure=Mock(return_value=verified),
            ImportOriginGuard=Mock(return_value=guard),
        )

        def git_show(_root: Path, _commit: str, path: str, _label: str) -> bytes:
            if path == bootstrap.BOOTSTRAP_PATH:
                return bootstrap_source
            raise AssertionError(path)

        canonical_empty = b"{}\n"
        with (
            patch.object(bootstrap, "_require_isolated_interpreter"),
            patch.object(
                bootstrap, "_resolve_commit", side_effect=lambda _r, value, _l: value
            ),
            patch.object(
                bootstrap,
                "_load_closure_module",
                return_value=(closure_module, closure_source),
            ),
            patch.object(bootstrap, "_authority_bytes", return_value=canonical_empty),
            patch.object(bootstrap, "_git_show", side_effect=git_show),
        ):
            outcome = bootstrap.execute(
                root=ROOT,
                implementation_commit="a" * 40,
                authority_commit="c" * 40,
                process_marker=bootstrap.PROCESS_MARKER,
                mode=mode,
                store_root=store_root,
            )
        operation = (
            runner.sid3_resume_preflight if mode == "preflight" else runner.sid3_resume
        )
        return dict(outcome), operation, guard

    def test_guarded_preflight_and_resume_dispatch_only_runner_owned_operation(
        self,
    ) -> None:
        preflight, preflight_call, preflight_guard = self._execute_fixture("preflight")
        self.assertEqual(preflight, {"safe_to_resume": True})
        preflight_call.assert_called_once()
        self.assertGreaterEqual(preflight_guard.audits, 2)

        resumed, resume_call, resume_guard = self._execute_fixture("resume")
        self.assertEqual(resumed, {"state": "advanced"})
        resume_call.assert_called_once()
        self.assertGreaterEqual(resume_guard.audits, 2)

    def test_store_override_cannot_escape_authorized_path(self) -> None:
        with self.assertRaisesRegex(bootstrap.SID3BootstrapError, "store path differs"):
            self._execute_fixture("preflight", store_root=Path("/tmp/not-authorized"))

    def test_process_marker_is_required_and_exact(self) -> None:
        with redirect_stderr(StringIO()), self.assertRaises(SystemExit):
            bootstrap._arguments(
                [
                    "--implementation-commit",
                    "a" * 40,
                    "--authority-commit",
                    "c" * 40,
                    "--preflight",
                ]
            )
        with self.assertRaisesRegex(bootstrap.SID3BootstrapError, "process marker"):
            with patch.object(bootstrap, "_require_isolated_interpreter"):
                bootstrap.execute(
                    root=ROOT,
                    implementation_commit="a" * 40,
                    authority_commit="c" * 40,
                    process_marker="wrong",
                    mode="preflight",
                )


if __name__ == "__main__":
    unittest.main()
