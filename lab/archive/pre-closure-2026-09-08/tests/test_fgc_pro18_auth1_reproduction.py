from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]


class Proto18Auth1ReproductionTests(unittest.TestCase):
    def test_fixed_canonical_result_reproduces(self) -> None:
        completed = subprocess.run(
            [sys.executable, "scripts/reproduce_fgc_pro18_auth1.py", "--verify"],
            cwd=ROOT, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn('"verified": true', completed.stdout)
