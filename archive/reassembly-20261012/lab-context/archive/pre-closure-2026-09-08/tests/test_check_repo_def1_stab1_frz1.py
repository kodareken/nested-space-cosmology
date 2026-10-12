"""Compact-only repository routing for FGC-1-DEF1-STAB1-FRZ1."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import inspect
from pathlib import Path
import unittest
from unittest.mock import patch

from scripts import check_repo
from scripts.repo_checks import catalog
from recursive_horizons.fgc import def1_stab1_frz1_certificate as cert


ROOT = Path(__file__).resolve().parents[1]
_CURRENT = ROOT / "mk/current-foundation.mk"
_CLOSED = ROOT / "mk/closed-live-targets.mk"
_COMPACT = (
    "fgc-def1-stab1-frz1:\n"
    "\tPYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B "
    "scripts/reproduce_fgc_def1_stab1_frz1.py --check\n"
)
_VERIFICATION = (
    "verify-fgc-def1-stab1-frz1: fgc-def1-stab1-frz1\n"
    "\tPYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. $(PYTHON) -B -m unittest \\\n"
    "\t\ttests/test_fgc_def1_stab1_frz1_certificate.py \\\n"
    "\t\ttests/test_check_repo_def1_stab1_frz1.py -v\n"
    "\tPYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/check_repo.py "
    "--only-def1-stab1-frz1\n"
)


class Def1Stab1Frz1RepositoryRouteTests(unittest.TestCase):
    def test_compact_result_routes_only_to_verify_compact(self) -> None:
        with (
            patch.object(cert, "verify_compact", return_value={}) as compact,
            patch.object(
                cert,
                "compose_canonical_artifacts",
                side_effect=AssertionError("live reconstruction"),
            ),
        ):
            self.assertEqual(catalog._check_def1_stab1_frz1_result(), [])
            compact.assert_called_once_with(catalog.REPOSITORY)

        source = inspect.getsource(catalog._check_def1_stab1_frz1_result)
        self.assertIn("verify_compact", source)
        self.assertNotIn("compose_canonical_artifacts", source)
        self.assertNotIn("--write-result", source)
        self.assertNotIn("runs/", source)
        self.assertNotIn("git", source.lower())

    def test_compact_failure_is_not_silenced(self) -> None:
        with patch.object(
            cert,
            "verify_compact",
            side_effect=cert.Def1Stab1Frz1Error("changed bytes"),
        ):
            errors = catalog._check_def1_stab1_frz1_result()
        self.assertEqual(len(errors), 1)
        self.assertIn("changed bytes", errors[0])
        self.assertTrue(errors[0].startswith("DEF1-STAB1-FRZ1 compact freeze differs:"))

    def test_focused_cli_does_not_run_git_catalog_or_other_certificates(self) -> None:
        output = StringIO()
        with (
            patch.object(
                check_repo, "_check_def1_stab1_frz1_result", return_value=[]
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
            self.assertEqual(check_repo.main(["--only-def1-stab1-frz1"]), 0)
            selected.assert_called_once_with()
        text = output.getvalue()
        self.assertIn("PASS FGC-1-DEF1-STAB1-FRZ1 compact freeze", text)

    def test_private_check_is_exported_on_the_dispatcher(self) -> None:
        self.assertIs(
            check_repo._check_def1_stab1_frz1_result,
            catalog._check_def1_stab1_frz1_result,
        )
        self.assertIn("_check_def1_stab1_frz1_result", check_repo.__all__)

    def test_aggregate_and_make_targets_are_compact_only(self) -> None:
        aggregate_source = inspect.getsource(check_repo.check_results)
        self.assertEqual(aggregate_source.count("_check_def1_stab1_frz1_result()"), 1)

        current = _CURRENT.read_text(encoding="utf-8")
        closed = _CLOSED.read_text(encoding="utf-8")
        historical = (ROOT / "mk/historical-certificates.mk").read_text(encoding="utf-8")
        self.assertIn(_COMPACT, current)
        self.assertIn(_VERIFICATION, current)
        self.assertNotIn("--write-result", _COMPACT + _VERIFICATION)
        self.assertNotIn("run-fgc-def1", _COMPACT + _VERIFICATION)
        self.assertNotRegex(current, r"(?m)^run-fgc-def1-stab1-frz1\s*:")
        self.assertNotRegex(current, r"(?m)^status-fgc-def1-stab1-frz1\s*:")
        self.assertNotIn("run-fgc-def1-stab1-frz1", closed)
        self.assertNotRegex(historical, r"(?m)^fgc-def1-stab1-frz1\s*:")
        self.assertIn("scripts/reproduce_fgc_def1_stab1_frz1.py --check", current)


if __name__ == "__main__":
    unittest.main()
