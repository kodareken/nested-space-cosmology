"""The default reading route must not launch scientific recomputation."""
import contextlib
import io
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import demonstrate


class DemonstrationReadingRoute(unittest.TestCase):
    def test_default_authenticates_and_displays_without_subprocesses(self):
        output = io.StringIO()
        with patch.object(sys, "argv", ["demonstrate.py"]), \
             patch.object(demonstrate.subprocess, "run",
                          side_effect=AssertionError("default demo launched a scientific subprocess")), \
             contextlib.redirect_stdout(output):
            self.assertEqual(demonstrate.main(), 0)
        rendered = output.getvalue()
        self.assertIn("Displaying recorded results", rendered)
        self.assertIn("Interactions between compact modes", rendered)
        self.assertIn("A consistent charged compact sector", rendered)
        self.assertIn("--recompute", rendered)


if __name__ == "__main__":
    unittest.main()
