"""Repository-check boundaries for prospective compact TDG9 AR1 authority."""

from __future__ import annotations

import ast
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import inspect
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts import check_repo


class TDG9AR1AUTH1RepositoryCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.repository = Path(self.temporary.name)

    def _write_bundle(self, *relative_paths: str) -> None:
        for relative in relative_paths:
            path = self.repository / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"placeholder\n")

    def test_absent_or_partial_bundle_is_rejected_without_reproduction(self) -> None:
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError("absent AR1 bundle was reproduced"),
            ),
        ):
            self.assertEqual(
                check_repo._check_tdg9_ar1_auth1_result(),
                ["FGC-1-TDG9-AR1-AUTH1 compact bundle is absent"],
            )

        self._write_bundle(check_repo._TDG9_AR1_AUTH1_BUNDLE[0])
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError("partial AR1 bundle was reproduced"),
            ),
        ):
            findings = check_repo._check_tdg9_ar1_auth1_result()
        self.assertEqual(len(findings), 1)
        self.assertIn("compact bundle is partial", findings[0])

    def test_unsafe_bundle_path_is_rejected_without_reproduction(self) -> None:
        relative = check_repo._TDG9_AR1_AUTH1_BUNDLE[0]
        target = self.repository / "unsafe-ar1-target"
        target.write_bytes(b"not a bundle leaf\n")
        path = self.repository / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.symlink_to(target)
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError("unsafe AR1 bundle was reproduced"),
            ),
        ):
            findings = check_repo._check_tdg9_ar1_auth1_result()
        self.assertEqual(len(findings), 1)
        self.assertIn("bundle path is unsafe", findings[0])

    def test_complete_bundle_routes_only_to_compact_reproducer(self) -> None:
        self._write_bundle(*check_repo._TDG9_AR1_AUTH1_BUNDLE)
        completed = subprocess.CompletedProcess(args=[], returncode=0, stdout="", stderr="")
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(check_repo.subprocess, "run", return_value=completed) as run,
        ):
            self.assertEqual(check_repo._check_tdg9_ar1_auth1_result(), [])
        run.assert_called_once()
        command = run.call_args.args[0]
        self.assertEqual(
            command[-2:],
            ["scripts/reproduce_fgc_tdg9_ar1_auth1.py", "--verify-compact"],
        )
        self.assertNotIn("--write", command)
        self.assertNotIn("--run", command)
        self.assertNotIn("runs/", " ".join(command))

    def test_focused_flag_is_store_blind_and_mutually_exclusive(self) -> None:
        with (
            patch.object(check_repo, "check_required", return_value=[]) as required,
            patch.object(
                check_repo, "_check_tdg9_ar1_auth1_result", return_value=[]
            ) as focused,
            redirect_stdout(StringIO()),
        ):
            self.assertEqual(check_repo.main(["--only-tdg9-ar1-auth1"]), 0)
        required.assert_called_once_with()
        focused.assert_called_once_with()

        with (
            redirect_stdout(StringIO()),
            redirect_stderr(StringIO()),
            self.assertRaises(SystemExit) as raised,
        ):
            check_repo.main(
                ["--only-tdg8-rcv3-pref2", "--only-tdg9-ar1-auth1"]
            )
        self.assertEqual(raised.exception.code, 2)

    def test_ordinary_route_and_foundation_use_compact_authority_only(self) -> None:
        source = inspect.getsource(check_repo.check_results)
        calls = {
            node.func.id
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }
        self.assertIn("_check_tdg9_ar1_auth1_result", calls)

        makefile = (
            Path(__file__).resolve().parents[1]
            / "mk"
            / "historical-certificates.mk"
        ).read_text(encoding="utf-8")
        foundation = makefile.split("verify-fgc-sf1-foundation:", 1)[1].split(
            "\n", 1
        )[0]
        self.assertIn("fgc-tdg9-ar1-auth1", foundation)
        self.assertNotIn("run-fgc-tdg9-ar1", foundation)


if __name__ == "__main__":
    unittest.main()
