"""Compact-only repository routing for FGC-1-TDG10-QA1-PREF1."""

from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import inspect
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import check_repo  # noqa: E402


_RESULT = "results/fgc-1-tdg10-qa1-pref1.json"
_HISTORICAL = "mk/historical-certificates.mk"
_CLASSIFICATION = (
    "independently_bound_retry3_ssprk3_sbp4_exact_complete_C_all_channel_pass_terminal"
)
_EXPECTED_BUNDLE = (
    "Makefile",
    "README.md",
    "configs/fgc/fgc-1-tdg10-qa1-pref1.toml",
    "docs/active-code-map.md",
    "docs/claim-ledger.md",
    "docs/fgc-runtime-matrix.md",
    "docs/fgc-tdg10-qa1-pref1.md",
    "docs/research-roadmap.md",
    _HISTORICAL,
    "results/README.md",
    _RESULT,
    "scripts/check_repo.py",
    "scripts/reproduce_fgc_tdg10_qa1_pref1.py",
    "src/recursive_horizons/fgc/evolution/tdg10_qa1_pref1_binder.py",
    "tests/test_fgc_tdg10_qa1_pref1_binder.py",
    "tests/test_check_repo_tdg10_qa1_pref1.py",
)
_ORDINARY = (
    "fgc-tdg10-qa1-pref1:\n"
    "\t$(EVOLUTION_PYTHON) scripts/reproduce_fgc_tdg10_qa1_pref1.py\n"
)
_VERIFICATION = (
    "verify-fgc-tdg10-qa1-pref1: fgc-tdg10-qa1-pref1\n"
    "\tPYTHONPATH=src:. $(EVOLUTION_PYTHON) -m unittest \\\n"
    "\t\ttests/test_fgc_tdg10_qa1_pref1_binder.py \\\n"
    "\t\ttests/test_check_repo_tdg10_qa1_pref1.py -v\n"
    "\t$(EVOLUTION_PYTHON) scripts/check_repo.py --only-tdg10-qa1-pref1\n"
)


class QA1PREF1RepositoryCheckTests(unittest.TestCase):
    def test_bundle_is_the_exact_compact_repository_surface(self) -> None:
        self.assertEqual(check_repo._TDG10_QA1_PREF1_BUNDLE, _EXPECTED_BUNDLE)
        self.assertEqual(len(_EXPECTED_BUNDLE), 16)
        self.assertIn(_HISTORICAL, _EXPECTED_BUNDLE)
        self.assertEqual(
            set(check_repo._TDG10_QA1_PREF1_PUBLIC_IDENTIFIERS),
            {
                "Makefile",
                "README.md",
                "docs/active-code-map.md",
                "docs/claim-ledger.md",
                "docs/fgc-runtime-matrix.md",
                "docs/fgc-tdg10-qa1-pref1.md",
                "docs/research-roadmap.md",
                _HISTORICAL,
                "results/README.md",
            },
        )
        for forbidden in (
            ".git",
            ".qdrant-initialized",
            "runs/",
            "raw/",
            "store/",
            "scripts/run_fgc_tdg10_qa1.py",
        ):
            self.assertFalse(
                any(
                    path == forbidden or path.startswith(forbidden)
                    for path in _EXPECTED_BUNDLE
                )
            )

    def test_compact_route_invokes_only_default_reproducer(self) -> None:
        completed = SimpleNamespace(returncode=0, stdout="", stderr="")
        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", return_value=True),
            patch.object(check_repo.subprocess, "run", return_value=completed) as run,
        ):
            self.assertEqual(check_repo._check_tdg10_qa1_pref1_result(), [])

        run.assert_called_once()
        self.assertEqual(
            run.call_args.args[0],
            [
                check_repo.sys.executable,
                "scripts/reproduce_fgc_tdg10_qa1_pref1.py",
            ],
        )
        self.assertIs(run.call_args.kwargs["stdin"], check_repo.subprocess.DEVNULL)
        flattened = " ".join(run.call_args.args[0])
        self.assertNotIn("--live", flattened)
        self.assertNotIn("run_fgc_tdg10_qa1", flattened)
        self.assertNotIn("runs/", flattened)

    def test_absent_result_is_the_only_partial_bundle_finding(self) -> None:
        def probe(relative: str, _label: str) -> bool:
            return relative != _RESULT

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", side_effect=probe),
            patch.object(check_repo, "_read_unique_regular_bytes") as read,
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg10_qa1_pref1_result()

        self.assertEqual(
            failures,
            ["FGC-1-TDG10-QA1-PREF1 compact bundle is partial; missing: " + _RESULT],
        )
        read.assert_not_called()
        run.assert_not_called()

    def test_unsafe_bundle_fails_closed_before_public_reads_or_reproducer(self) -> None:
        with (
            patch.object(
                check_repo,
                "_probe_unique_regular_leaf",
                side_effect=ValueError("nofollow rejected symlink"),
            ),
            patch.object(check_repo, "_read_unique_regular_bytes") as read,
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg10_qa1_pref1_result()

        self.assertEqual(len(failures), 1)
        self.assertIn("bundle path is unsafe: Makefile", failures[0])
        read.assert_not_called()
        run.assert_not_called()

    def test_public_identifier_drift_fails_before_reproducer(self) -> None:
        def read(relative: str, _label: str) -> bytes:
            if relative == "README.md":
                return b"public surface intentionally missing the artifact"
            return (ROOT / relative).read_bytes()

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", return_value=True),
            patch.object(check_repo, "_read_unique_regular_bytes", side_effect=read),
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg10_qa1_pref1_result()

        self.assertEqual(len(failures), 3)
        self.assertTrue(all("README.md" in failure for failure in failures))
        run.assert_not_called()

    def test_removing_historical_include_fails_before_reproducer(self) -> None:
        def read(relative: str, _label: str) -> bytes:
            if relative == "Makefile":
                return (
                    ".DEFAULT_GOAL := test\n"
                    "include mk/current-foundation.mk\n"
                    "include mk/closed-live-targets.mk\n"
                ).encode("utf-8")
            return (ROOT / relative).read_bytes()

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", return_value=True),
            patch.object(check_repo, "_read_unique_regular_bytes", side_effect=read),
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg10_qa1_pref1_result()

        self.assertEqual(
            failures,
            [
                "FGC-1-TDG10-QA1-PREF1 public-surface identifier is absent: "
                "Makefile: include mk/historical-certificates.mk"
            ],
        )
        run.assert_not_called()

    def test_removing_historical_target_fails_before_reproducer(self) -> None:
        def read(relative: str, _label: str) -> bytes:
            if relative == _HISTORICAL:
                return b"# historical certificates without the QA1 PREF1 targets\n"
            return (ROOT / relative).read_bytes()

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", return_value=True),
            patch.object(check_repo, "_read_unique_regular_bytes", side_effect=read),
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg10_qa1_pref1_result()

        self.assertEqual(len(failures), 4)
        self.assertTrue(all(_HISTORICAL in failure for failure in failures))
        self.assertTrue(any("fgc-tdg10-qa1-pref1" in failure for failure in failures))
        self.assertTrue(
            any("scripts/reproduce_fgc_tdg10_qa1_pref1.py" in failure for failure in failures)
        )
        self.assertTrue(
            any(
                "scripts/check_repo.py --only-tdg10-qa1-pref1" in failure
                for failure in failures
            )
        )
        run.assert_not_called()

    def test_owner_classification_drift_fails_before_reproducer(self) -> None:
        def read(relative: str, _label: str) -> bytes:
            if relative == "docs/fgc-tdg10-qa1-pref1.md":
                return (
                    (ROOT / relative)
                    .read_text(encoding="utf-8")
                    .replace(_CLASSIFICATION, "relabeled_classification")
                    .encode("utf-8")
                )
            return (ROOT / relative).read_bytes()

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", return_value=True),
            patch.object(check_repo, "_read_unique_regular_bytes", side_effect=read),
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg10_qa1_pref1_result()

        self.assertEqual(len(failures), 1)
        self.assertIn("docs/fgc-tdg10-qa1-pref1.md", failures[0])
        self.assertIn(_CLASSIFICATION, failures[0])
        run.assert_not_called()

    def test_focused_cli_routes_only_required_and_pref1_checks(self) -> None:
        output = StringIO()
        with (
            patch.object(check_repo, "check_required", return_value=[]) as required,
            patch.object(
                check_repo, "_check_tdg10_qa1_pref1_result", return_value=[]
            ) as pref1,
            patch.object(
                check_repo,
                "check_results",
                side_effect=AssertionError("focused route entered aggregate checks"),
            ),
            redirect_stdout(output),
        ):
            self.assertEqual(check_repo.main(["--only-tdg10-qa1-pref1"]), 0)

        required.assert_called_once_with()
        pref1.assert_called_once_with()
        self.assertIn("PASS required files", output.getvalue())
        self.assertIn(
            "PASS TDG10 compact retry-3 exact complete-C terminal binder",
            output.getvalue(),
        )

    def test_aggregate_and_make_targets_use_compact_verifier_only(self) -> None:
        aggregate_source = inspect.getsource(check_repo.check_results)
        self.assertEqual(aggregate_source.count("_check_tdg10_qa1_pref1_result()"), 1)

        makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
        historical = (ROOT / _HISTORICAL).read_text(encoding="utf-8")
        current = (ROOT / "mk/current-foundation.mk").read_text(encoding="utf-8")
        self.assertIn("include mk/historical-certificates.mk", makefile)
        self.assertNotIn(_ORDINARY, makefile)
        self.assertNotIn(_VERIFICATION, makefile)
        self.assertIn(_ORDINARY, historical)
        self.assertIn(_VERIFICATION, historical)
        self.assertNotIn("--live", _ORDINARY + _VERIFICATION)
        self.assertIn("fgc-tdg10-qa1-pref1", current)

        foundation = next(
            line
            for line in historical.splitlines()
            if line.startswith("verify-fgc-sf1-foundation:")
        )
        self.assertIn("fgc-tdg10-qa1-frz1 fgc-tdg10-qa1-pref1", foundation)


if __name__ == "__main__":
    unittest.main()
