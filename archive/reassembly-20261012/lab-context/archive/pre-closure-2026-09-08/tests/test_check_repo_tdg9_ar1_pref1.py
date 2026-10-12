from __future__ import annotations

from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]


class CheckRepoTDG9AR1PREF1Tests(unittest.TestCase):
    def test_compact_audit_mode_passes_without_raw_access(self) -> None:
        completed = subprocess.run(
            ["/opt/homebrew/Caskroom/miniconda/base/bin/python3", "scripts/check_repo.py", "--only-tdg9-ar1-pref1"],
            cwd=ROOT, capture_output=True, text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)


if __name__ == "__main__":
    unittest.main()
