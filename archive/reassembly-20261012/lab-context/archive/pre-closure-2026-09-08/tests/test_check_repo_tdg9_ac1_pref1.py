"""Routing controls for the compact AC1 PREF1 repository audit."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import check_repo  # noqa: E402


class AC1PREF1RepositoryCheckTests(unittest.TestCase):
    def test_ac1_pref1_calls_only_its_compact_reproducer(self) -> None:
        completed = SimpleNamespace(returncode=0, stdout="", stderr="")
        with patch.object(check_repo.subprocess, "run", return_value=completed) as run:
            self.assertEqual(check_repo._check_tdg9_ac1_pref1_result(), [])
        commands = [tuple(call.args[0]) for call in run.call_args_list]
        self.assertEqual(
            commands,
            [
                (
                    check_repo.sys.executable,
                    "scripts/reproduce_fgc_tdg9_ac1_pref1.py",
                )
            ],
        )
        self.assertFalse(
            any(
                "run_fgc_tdg9_ac1" in item or "--live" in item
                for command in commands
                for item in command
            )
        )

    def test_partial_ac1_pref1_bundle_fails_without_reproducer(self) -> None:
        missing = ROOT / check_repo._TDG9_AC1_PREF1_BUNDLE[-1]
        original = check_repo._probe_unique_regular_leaf

        def probe(relative: str, label: str) -> bool:
            if relative == str(missing.relative_to(ROOT)):
                return False
            return original(relative, label)

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", side_effect=probe),
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg9_ac1_pref1_result()
        self.assertEqual(len(failures), 1)
        self.assertIn("compact bundle is partial", failures[0])
        run.assert_not_called()

    def test_compact_audit_passes_without_raw_access(self) -> None:
        completed = subprocess.run(
            [sys.executable, "scripts/check_repo.py", "--only-tdg9-ac1-pref1"],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            completed.returncode, 0, completed.stdout + completed.stderr
        )


if __name__ == "__main__":
    unittest.main()
