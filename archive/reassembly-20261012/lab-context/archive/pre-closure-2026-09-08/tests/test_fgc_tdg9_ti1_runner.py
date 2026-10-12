"""Fail-closed controls for the immutable invalid TI1 attempt."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]


class HistoricalTI1RunnerTests(unittest.TestCase):
    def test_historical_entry_fails_closed_with_or_without_old_arguments(self) -> None:
        for arguments in (
            (),
            ("--authority-commit", "a" * 40),
            ("--authority-commit", "a" * 40, "--status"),
        ):
            with self.subTest(arguments=arguments):
                completed = subprocess.run(
                    [sys.executable, "scripts/run_fgc_tdg9_ti1.py", *arguments],
                    cwd=ROOT,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False,
                )
                self.assertEqual(completed.returncode, 2)
                self.assertEqual(completed.stderr, "")
                record = json.loads(completed.stdout)
                self.assertEqual(
                    record["classification"],
                    "historical_invalid_authority_retired",
                )
                self.assertFalse(record["state_advance_authorized"])


if __name__ == "__main__":
    unittest.main()
