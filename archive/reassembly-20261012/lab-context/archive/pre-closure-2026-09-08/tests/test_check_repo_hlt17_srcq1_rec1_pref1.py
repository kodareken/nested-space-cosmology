"""Repository-check boundary for compact REC1-PREF1 verification."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.repo_checks import temporal_selection as checks


ROOT = Path(__file__).resolve().parents[1]
CONFIG = "configs/fgc/fgc-1-hlt17-srcq1-rec1-pref1.toml"
RESULT = "results/fgc-1-hlt17-srcq1-rec1-pref1.json"


class REC1PREF1RepositoryCheckTests(unittest.TestCase):
    def _copy(self, root: Path, relative: str) -> None:
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((ROOT / relative).read_bytes())

    def test_compact_bundle_passes_without_raw_or_runs(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            self._copy(root, CONFIG)
            self._copy(root, RESULT)
            with patch.object(checks, "REPOSITORY", root):
                self.assertEqual(checks._check_hlt17_srcq1_rec1_pref1_result(), [])
            self.assertFalse((root / "runs").exists())

    def test_absent_or_changed_compact_fails_closed(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            with patch.object(checks, "REPOSITORY", root):
                failures = checks._check_hlt17_srcq1_rec1_pref1_result()
            self.assertTrue(failures)
            self._copy(root, CONFIG)
            self._copy(root, RESULT)
            result = root / RESULT
            result.write_bytes(result.read_bytes().replace(b"all_six", b"all_6", 1))
            with patch.object(checks, "REPOSITORY", root):
                failures = checks._check_hlt17_srcq1_rec1_pref1_result()
            self.assertTrue(failures)

    def test_check_route_names_only_compact_verifier(self) -> None:
        source = Path(checks.__file__).read_text(encoding="utf-8")
        block = source.split(
            "def _check_hlt17_srcq1_rec1_pref1_result", 1
        )[1].split("\n\ndef ", 1)[0]
        self.assertIn("verify_compact", block)
        for forbidden in ("bind_live", "runs/", "--live", "git_read"):
            self.assertNotIn(forbidden, block)


if __name__ == "__main__":
    unittest.main()
