"""Store-blind repository-check controls for prospective LOC1."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts import check_repo


class LOC1RepositoryCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.repository = Path(self.temporary.name)

    def _write(self, *paths: str) -> None:
        for relative in paths:
            path = self.repository / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"placeholder\n")

    def test_partial_bundle_fails_without_reproducer(self) -> None:
        self._write(check_repo._TDG9_LOC1_FRZ1_BUNDLE[0])
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(check_repo.subprocess, "run", side_effect=AssertionError("must stay compact")),
        ):
            findings = check_repo._check_tdg9_loc1_frz1_result()
        self.assertEqual(len(findings), 1)
        self.assertIn("partial", findings[0])

    def test_complete_bundle_calls_only_compact_reproducer(self) -> None:
        self._write(*check_repo._TDG9_LOC1_FRZ1_BUNDLE)
        completed = subprocess.CompletedProcess(args=[], returncode=0, stdout="", stderr="")
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(check_repo.subprocess, "run", return_value=completed) as run,
        ):
            self.assertEqual(check_repo._check_tdg9_loc1_frz1_result(), [])
        command = run.call_args.args[0]
        self.assertEqual(command[-2:], ["scripts/reproduce_fgc_tdg9_loc1_frz1.py", "--verify-compact"])
        self.assertNotIn("runs/", " ".join(command))
        self.assertNotIn("run_fgc_tdg9_loc1.py", " ".join(command))

    def test_focused_flag_routes_only_compact_check(self) -> None:
        with (
            patch.object(check_repo, "check_required", return_value=[]),
            patch.object(check_repo, "_check_tdg9_loc1_frz1_result", return_value=[]) as focused,
            redirect_stdout(StringIO()),
        ):
            self.assertEqual(check_repo.main(["--only-tdg9-loc1-frz1"]), 0)
        focused.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
