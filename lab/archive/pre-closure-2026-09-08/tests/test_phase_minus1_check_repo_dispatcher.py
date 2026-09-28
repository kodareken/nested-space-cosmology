"""Phase -1 repository-check dispatcher compatibility invariants."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch
import unittest

from scripts import check_repo
from scripts.repo_checks import (
    _shared,
    calibration,
    catalog,
    core,
    progression,
    public,
    temporal_diagnostics,
)


ROOT = Path(__file__).resolve().parents[1]
DISPATCHER = ROOT / "scripts/check_repo.py"
DOMAINS = (
    _shared,
    public,
    calibration,
    progression,
    temporal_diagnostics,
    catalog,
    core,
)


class PhaseMinusOneCheckRepoDispatcherTests(unittest.TestCase):
    def test_dispatcher_stays_thin_and_domains_exist(self) -> None:
        source = DISPATCHER.read_text(encoding="utf-8")
        self.assertLessEqual(len(source.splitlines()), 100)
        self.assertIn("class _CompatibilityFacade", source)
        for name in (
            "core.py",
            "public.py",
            "calibration.py",
            "progression.py",
            "temporal_diagnostics.py",
            "catalog.py",
        ):
            self.assertTrue((ROOT / "scripts/repo_checks" / name).is_file())

    def test_legacy_repository_patch_propagates_to_every_domain(self) -> None:
        replacement = ROOT / "dispatcher-fixture-that-need-not-exist"
        originals = {module: module.REPOSITORY for module in DOMAINS}
        with patch.object(check_repo, "REPOSITORY", replacement):
            self.assertIs(check_repo.REPOSITORY, replacement)
            for module in DOMAINS:
                self.assertIs(module.REPOSITORY, replacement)
        for module, original in originals.items():
            self.assertIs(module.REPOSITORY, original)

    def test_legacy_private_helper_patch_reaches_core_or_owner(self) -> None:
        sentinel = object()
        with patch.object(check_repo, "check_required", sentinel):
            self.assertIs(core.check_required, sentinel)
            self.assertIs(public.check_required, sentinel)
        with patch.object(check_repo, "_check_tdg10_qa2_pref1_result", sentinel):
            self.assertIs(core._check_tdg10_qa2_pref1_result, sentinel)
            self.assertIs(
                temporal_diagnostics._check_tdg10_qa2_pref1_result,
                sentinel,
            )

    def test_facade_reexports_existing_private_check_names(self) -> None:
        for name in (
            "_check_pro19_pref28_compact_result",
            "_check_tdg8_rcv3_pref2_result",
            "_check_tdg9_ar1_auth1_result",
            "_check_tdg10_qa1_pref1_result",
            "_check_tdg10_qa2_pref1_result",
            "_read_unique_regular_bytes",
        ):
            self.assertTrue(hasattr(check_repo, name), name)

    def test_focused_qa2_mode_does_not_invoke_git_aware_catalog(self) -> None:
        with (
            patch.object(check_repo, "check_required", return_value=[]),
            patch.object(check_repo, "_check_tdg10_qa2_pref1_result", return_value=[]),
            patch.object(
                check_repo,
                "check_artifact_catalog",
                side_effect=AssertionError("focused compact mode inspected Git catalog"),
            ),
            redirect_stdout(StringIO()),
        ):
            self.assertEqual(check_repo.main(["--only-tdg10-qa2-pref1"]), 0)

    def test_full_audit_retains_catalog_coverage(self) -> None:
        names = (
            "check_required", "check_links", "check_public_text", "check_fgc_alignment",
            "check_results", "check_git_data_boundary", "check_reference_sequence",
        )
        patches = [patch.object(check_repo, name, return_value=[]) for name in names]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        with (
            patch.object(check_repo, "check_artifact_catalog", return_value=[]) as catalog,
            redirect_stdout(StringIO()),
        ):
            self.assertEqual(check_repo.main([]), 0)
        catalog.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
