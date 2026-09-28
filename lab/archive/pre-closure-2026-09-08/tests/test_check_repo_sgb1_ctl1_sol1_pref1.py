from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts import check_repo
from scripts.repo_checks import catalog as checks
from recursive_horizons.fgc import sgb1_ctl1_sol1_pref1_binder as binder


ROOT = Path(__file__).resolve().parents[1]


class SOL1PREF1RepositoryRouteTests(unittest.TestCase):
    def test_compact_route_calls_only_verify_compact(self) -> None:
        with patch.object(binder, "verify_compact", return_value={}) as compact:
            self.assertEqual(checks._check_sgb1_ctl1_sol1_pref1_result(), [])
        compact.assert_called_once_with(checks.REPOSITORY)

    def test_compact_failure_is_not_silenced(self) -> None:
        with patch.object(
            binder, "verify_compact", side_effect=binder.SGBLSOL1PREF1Error("changed")
        ):
            errors = checks._check_sgb1_ctl1_sol1_pref1_result()
        self.assertEqual(len(errors), 1)
        self.assertIn("changed", errors[0])

    def test_focused_cli_is_compact_only(self) -> None:
        with (
            patch.object(check_repo, "_check_sgb1_ctl1_sol1_pref1_result", return_value=[]) as selected,
            patch.object(check_repo, "check_artifact_catalog", side_effect=AssertionError()),
            patch.object(check_repo, "check_results", side_effect=AssertionError()),
        ):
            self.assertEqual(check_repo.main(["--only-sgb1-ctl1-sol1-pref1"]), 0)
        selected.assert_called_once_with()

    def test_compact_bundle_passes_without_git_source_or_runs(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            for relative in (binder.CONFIG_PATH, binder.RESULT_PATH):
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes((ROOT / relative).read_bytes())
            with patch.object(checks, "REPOSITORY", root):
                self.assertEqual(checks._check_sgb1_ctl1_sol1_pref1_result(), [])
            self.assertFalse((root / ".git").exists())
            self.assertFalse((root / "src").exists())
            self.assertFalse((root / "runs").exists())

    def test_make_targets_are_compact_only(self) -> None:
        source = (ROOT / "mk/current-foundation.mk").read_text(encoding="utf-8")
        self.assertRegex(source, r"(?m)^fgc-sgb1-ctl1-sol1-pref1:")
        self.assertRegex(source, r"(?m)^verify-fgc-sgb1-ctl1-sol1-pref1:")
        self.assertNotRegex(source, r"(?m)^(run|status)-fgc-sgb1-ctl1-sol1-pref1:")


if __name__ == "__main__":
    unittest.main()
