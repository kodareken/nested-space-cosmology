"""Routing controls for the compact TI2 repository audit."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import check_repo  # noqa: E402


class TI2RepositoryCheckTests(unittest.TestCase):
    def test_ti1_and_ti2_call_only_their_compact_reproducers(self) -> None:
        completed = SimpleNamespace(returncode=0, stdout="", stderr="")
        with patch.object(check_repo.subprocess, "run", return_value=completed) as run:
            self.assertEqual(check_repo._check_tdg9_ti1_frz1_result(), [])
            self.assertEqual(check_repo._check_tdg9_ti2_frz1_result(), [])
        commands = [tuple(call.args[0]) for call in run.call_args_list]
        self.assertEqual(
            commands,
            [
                (
                    check_repo.sys.executable,
                    "scripts/reproduce_fgc_tdg9_ti1_frz1.py",
                    "--verify-compact",
                ),
                (
                    check_repo.sys.executable,
                    "scripts/reproduce_fgc_tdg9_ti2_frz1.py",
                    "--verify-compact",
                ),
            ],
        )
        self.assertFalse(
            any(
                "run_fgc_tdg9_ti" in item
                for command in commands
                for item in command
            )
        )

    def test_partial_ti2_bundle_fails_without_reproducer(self) -> None:
        missing = ROOT / check_repo._TDG9_TI2_FRZ1_BUNDLE[-1]
        original = check_repo._probe_unique_regular_leaf

        def probe(relative: str, label: str) -> bool:
            if relative == str(missing.relative_to(ROOT)):
                return False
            return original(relative, label)

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", side_effect=probe),
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg9_ti2_frz1_result()
        self.assertEqual(len(failures), 1)
        self.assertIn("compact bundle is partial", failures[0])
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
