"""Compact-only IMP1 repository and Make routing; no numerical qualification."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import inspect
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import check_repo  # noqa: E402
from scripts import reproduce_fgc_tdg11_imp1 as certificate  # noqa: E402
from scripts.repo_checks import temporal_selection  # noqa: E402


class IMP1RepositoryRouteTests(unittest.TestCase):
    def test_repository_helper_uses_only_compact_check(self):
        with (
            patch.object(certificate, "verify_compact", return_value={}) as compact,
            patch.object(
                certificate, "qualify", side_effect=AssertionError("qualification")
            ),
        ):
            self.assertEqual(temporal_selection._check_tdg11_imp1_result(), [])
        compact.assert_called_once_with(temporal_selection.REPOSITORY)

    def test_partial_or_bad_compact_input_is_not_silenced(self):
        for error in (
            FileNotFoundError("partial bundle"),
            certificate.IMP1CertificateError("changed bytes"),
        ):
            with patch.object(certificate, "verify_compact", side_effect=error):
                errors = temporal_selection._check_tdg11_imp1_result()
            self.assertEqual(len(errors), 1)
            self.assertIn(str(error), errors[0])

    def test_focused_cli_avoids_git_catalog_and_other_results(self):
        with (
            patch.object(check_repo, "check_required", return_value=[]) as required,
            patch.object(
                check_repo, "_check_tdg11_imp1_result", return_value=[]
            ) as compact,
            patch.object(
                check_repo,
                "check_artifact_catalog",
                side_effect=AssertionError("catalog Git"),
            ),
            patch.object(
                check_repo, "check_results", side_effect=AssertionError("other results")
            ),
            redirect_stdout(StringIO()) as output,
        ):
            self.assertEqual(check_repo.main(["--only-tdg11-imp1"]), 0)
        required.assert_called_once_with()
        compact.assert_called_once_with()
        self.assertIn("PASS TDG11 compact qualified implementation", output.getvalue())

    def test_other_focused_flag_cannot_be_combined(self):
        with patch("sys.stderr"):
            with self.assertRaises(SystemExit) as caught:
                check_repo.main(["--only-tdg11-imp1", "--only-tdg11-msel1-pref1"])
        self.assertEqual(caught.exception.code, 2)

    def test_compatibility_facade_exports_helper_and_repository_override(self):
        self.assertIs(
            check_repo._check_tdg11_imp1_result,
            temporal_selection._check_tdg11_imp1_result,
        )
        self.assertIn("_check_tdg11_imp1_result", check_repo.__all__)
        original = temporal_selection.REPOSITORY
        with patch.object(check_repo, "REPOSITORY", Path("/synthetic-repository")):
            self.assertEqual(
                temporal_selection.REPOSITORY, Path("/synthetic-repository")
            )
        self.assertEqual(temporal_selection.REPOSITORY, original)

    def test_aggregate_includes_imp1_once(self):
        self.assertEqual(
            inspect.getsource(check_repo.check_results).count(
                "_check_tdg11_imp1_result()"
            ),
            1,
        )

    def test_make_compact_target_never_qualifies_or_starts_a_live_runner(self):
        process = subprocess.run(
            ["make", "-n", "fgc-tdg11-imp1"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("scripts/reproduce_fgc_tdg11_imp1.py --check", process.stdout)
        self.assertNotIn("--qualify", process.stdout)
        self.assertNotIn("--write-result", process.stdout)
        self.assertNotIn("run_fgc", process.stdout)
        current = (ROOT / "mk/current-foundation.mk").read_text()
        self.assertNotRegex(current, r"(?m)^(?:run|qualify|bind)-fgc-tdg11-imp1\s*:")
        self.assertIn("verify-fgc-tdg11-imp1: fgc-tdg11-imp1", current)


if __name__ == "__main__":
    unittest.main()
