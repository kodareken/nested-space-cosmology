"""Routing controls for the compact TI1 repository audit."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import check_repo  # noqa: E402


class TI1RepositoryCheckTests(unittest.TestCase):
    def test_complete_bundle_calls_only_compact_reproducer(self) -> None:
        completed = SimpleNamespace(returncode=0, stdout="", stderr="")
        with patch.object(check_repo.subprocess, "run", return_value=completed) as run:
            self.assertEqual(check_repo._check_tdg9_ti1_frz1_result(), [])
        commands = [tuple(call.args[0]) for call in run.call_args_list]
        self.assertEqual(
            commands,
            [
                (
                    check_repo.sys.executable,
                    "scripts/reproduce_fgc_tdg9_ti1_frz1.py",
                    "--verify-compact",
                )
            ],
        )
        self.assertFalse(any("run_fgc_tdg9_ti1.py" in item for command in commands for item in command))


if __name__ == "__main__":
    unittest.main()
