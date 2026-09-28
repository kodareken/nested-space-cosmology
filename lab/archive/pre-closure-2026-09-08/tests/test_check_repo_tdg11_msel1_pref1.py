"""Compact-only repository routing for FGC-1-TDG11-MSEL1-PREF1."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import inspect
from pathlib import Path
import unittest
from unittest.mock import patch

from scripts import check_repo
from scripts.repo_checks import temporal_selection
from recursive_horizons.fgc.evolution import tdg11_msel1_pref1_binder as binder


ROOT = Path(__file__).resolve().parents[1]
_CURRENT = ROOT / "mk/current-foundation.mk"
_CLOSED = ROOT / "mk/closed-live-targets.mk"
_COMPACT = (
    "fgc-tdg11-msel1-pref1:\n"
    "\tPYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B "
    "scripts/reproduce_fgc_tdg11_msel1_pref1.py --check\n"
)
_VERIFICATION = (
    "verify-fgc-tdg11-msel1-pref1: fgc-tdg11-msel1-pref1\n"
    "\tPYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(EVOLUTION_PYTHON) -B -m unittest \\\n"
    "\t\ttests/test_fgc_tdg11_msel1_pref1_binder.py \\\n"
    "\t\ttests/test_fgc_tdg11_msel1_pref1_protocol.py \\\n"
    "\t\ttests/test_fgc_tdg11_msel1_pref1_reconstruction.py \\\n"
    "\t\ttests/test_fgc_tdg11_msel1_pref1_localization.py \\\n"
    "\t\ttests/test_check_repo_tdg11_msel1_pref1.py -v\n"
    "\tPYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py "
    "--only-tdg11-msel1-pref1\n"
)


class TDG11PREF1RepositoryRouteTests(unittest.TestCase):
    def test_compact_result_routes_only_to_verify_compact(self) -> None:
        with (
            patch.object(binder, "verify_compact", return_value={}) as compact,
            patch.object(
                binder,
                "bind_raw_result",
                side_effect=AssertionError("live bind"),
            ),
        ):
            self.assertEqual(temporal_selection._check_tdg11_msel1_pref1_result(), [])
            compact.assert_called_once_with(temporal_selection.REPOSITORY)

        source = inspect.getsource(temporal_selection._check_tdg11_msel1_pref1_result)
        self.assertIn("verify_compact", source)
        self.assertNotIn("bind_raw_result", source)
        self.assertNotIn("--bind", source)
        self.assertNotIn("--write-result", source)
        self.assertNotIn("runs/", source)
        self.assertNotIn("store", source)

    def test_compact_failure_is_not_silenced(self) -> None:
        with patch.object(
            binder,
            "verify_compact",
            side_effect=binder.TDG11PREF1Error("changed bytes"),
        ):
            errors = temporal_selection._check_tdg11_msel1_pref1_result()
        self.assertEqual(len(errors), 1)
        self.assertIn("changed bytes", errors[0])
        self.assertTrue(errors[0].startswith("TDG11-MSEL1-PREF1 compact result differs:"))

    def test_focused_cli_does_not_run_git_catalog_or_other_certificates(self) -> None:
        output = StringIO()
        with (
            patch.object(check_repo, "check_required", return_value=[]) as required,
            patch.object(
                check_repo, "_check_tdg11_msel1_pref1_result", return_value=[]
            ) as selected,
            patch.object(
                check_repo,
                "check_artifact_catalog",
                side_effect=AssertionError("catalog Git"),
            ),
            patch.object(
                check_repo, "check_results", side_effect=AssertionError("all results")
            ),
            redirect_stdout(output),
        ):
            self.assertEqual(check_repo.main(["--only-tdg11-msel1-pref1"]), 0)
            required.assert_called_once_with()
            selected.assert_called_once_with()
        text = output.getvalue()
        self.assertIn("PASS required files", text)
        self.assertIn("PASS TDG11 compact independent selection binder", text)

    def test_compatibility_facade_preserves_repository_override(self) -> None:
        original = temporal_selection.REPOSITORY
        with patch.object(check_repo, "REPOSITORY", Path("/synthetic-repository")):
            self.assertEqual(
                temporal_selection.REPOSITORY, Path("/synthetic-repository")
            )
        self.assertEqual(temporal_selection.REPOSITORY, original)

    def test_private_check_is_exported_on_the_dispatcher(self) -> None:
        self.assertIs(
            check_repo._check_tdg11_msel1_pref1_result,
            temporal_selection._check_tdg11_msel1_pref1_result,
        )
        self.assertIn("_check_tdg11_msel1_pref1_result", check_repo.__all__)
        self.assertIn("_check_tdg11_msel1_frz1_result", check_repo.__all__)

    def test_aggregate_and_make_targets_are_compact_only(self) -> None:
        aggregate_source = inspect.getsource(check_repo.check_results)
        self.assertEqual(aggregate_source.count("_check_tdg11_msel1_pref1_result()"), 1)
        self.assertEqual(aggregate_source.count("_check_tdg10_qa2_pref1_result()"), 1)

        current = _CURRENT.read_text(encoding="utf-8")
        closed = _CLOSED.read_text(encoding="utf-8")
        historical = (ROOT / "mk/historical-certificates.mk").read_text(encoding="utf-8")
        self.assertIn(_COMPACT, current)
        self.assertIn(_VERIFICATION, current)
        self.assertNotIn("--bind", _COMPACT + _VERIFICATION)
        self.assertNotIn("--write-result", _COMPACT + _VERIFICATION)
        self.assertNotIn("run-fgc-tdg11-msel1", _COMPACT + _VERIFICATION)
        self.assertNotRegex(current, r"(?m)^run-fgc-tdg11-msel1\s*:")
        self.assertIn("\trun-fgc-tdg11-msel1\n", closed.replace(" \\\n", "\n"))
        self.assertNotRegex(historical, r"(?m)^run-fgc-tdg11-msel1\s*:")
        self.assertNotIn("bind-fgc-tdg11-msel1-pref1", current)
        self.assertIn("scripts/reproduce_fgc_tdg11_msel1_frz1.py --check", current)

    def test_pref1_modules_use_the_owned_prefix(self) -> None:
        for name in ("binder", "protocol", "reconstruction", "localization"):
            relative = (
                f"src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_{name}.py"
            )
            self.assertTrue((ROOT / relative).is_file(), relative)


if __name__ == "__main__":
    unittest.main()
