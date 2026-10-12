"""TDG11 compact routing must not acquire live diagnostic authority."""

from __future__ import annotations

from contextlib import redirect_stdout
import io
from pathlib import Path
import unittest
from unittest.mock import patch

from scripts import check_repo
from scripts.repo_checks import temporal_selection
from recursive_horizons.fgc.evolution import tdg11_msel1_authority as authority


class TDG11RepositoryRouteTests(unittest.TestCase):
    def test_compact_result_routes_only_to_compact_owner(self):
        with (
            patch.object(
                authority, "validate_compact_bundle", return_value={}
            ) as compact,
            patch.object(
                authority,
                "authorize_execution",
                side_effect=AssertionError("live authority"),
            ),
            patch.object(
                authority, "snapshot_store", side_effect=AssertionError("store")
            ),
        ):
            self.assertEqual(temporal_selection._check_tdg11_msel1_frz1_result(), [])
            compact.assert_called_once_with(temporal_selection.REPOSITORY)

    def test_compact_failure_is_not_silenced(self):
        with patch.object(
            authority,
            "validate_compact_bundle",
            side_effect=authority.MSEL1AuthorityError("changed bytes"),
        ):
            errors = temporal_selection._check_tdg11_msel1_frz1_result()
        self.assertEqual(len(errors), 1)
        self.assertIn("changed bytes", errors[0])

    def test_focused_cli_does_not_run_git_catalog_or_other_certificates(self):
        with (
            patch.object(
                check_repo, "_check_tdg11_msel1_frz1_result", return_value=[]
            ) as selected,
            patch.object(
                check_repo,
                "check_artifact_catalog",
                side_effect=AssertionError("catalog Git"),
            ),
            patch.object(
                check_repo, "check_results", side_effect=AssertionError("all results")
            ),
            redirect_stdout(io.StringIO()),
        ):
            self.assertEqual(check_repo.main(["--only-tdg11-msel1-frz1"]), 0)
            selected.assert_called_once_with()

    def test_compatibility_facade_preserves_repository_override(self):
        original = temporal_selection.REPOSITORY
        with patch.object(check_repo, "REPOSITORY", Path("/synthetic-repository")):
            self.assertEqual(
                temporal_selection.REPOSITORY, Path("/synthetic-repository")
            )
        self.assertEqual(temporal_selection.REPOSITORY, original)


if __name__ == "__main__":
    unittest.main()
