"""Routing controls for the compact LOC2 repository audit."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]

from scripts import check_repo  # noqa: E402


class LOC2RepositoryCheckTests(unittest.TestCase):
    def test_complete_bundle_calls_only_compact_reproducers(self) -> None:
        completed = SimpleNamespace(returncode=0, stdout="", stderr="")
        with patch.object(check_repo.subprocess, "run", return_value=completed) as run:
            self.assertEqual(check_repo._check_tdg9_loc2_frz1_result(), [])
        commands = [tuple(call.args[0]) for call in run.call_args_list]
        self.assertIn((check_repo.sys.executable, "scripts/reproduce_fgc_tdg9_loc1_err1.py", "--verify-compact"), commands)
        self.assertIn((check_repo.sys.executable, "scripts/reproduce_fgc_tdg9_loc2_frz1.py", "--verify-compact"), commands)
        self.assertFalse(any("run_fgc_tdg9_loc2.py" in command for command in commands))


if __name__ == "__main__":
    unittest.main()
