"""Routing controls for the compact TDG10-QA2 repository audit."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import check_repo  # noqa: E402


_EXPECTED_BUNDLE = (
    "Makefile",
    "README.md",
    "configs/fgc/fgc-1-tdg10-qa2-frz1.toml",
    "docs/claim-ledger.md",
    "docs/fgc-runtime-matrix.md",
    "docs/fgc-tdg10-qa2-frz1.md",
    "docs/research-roadmap.md",
    "results/README.md",
    "results/fgc-1-tdg10-qa2-frz1.json",
    "scripts/check_repo.py",
    "scripts/reproduce_fgc_tdg10_qa2_frz1.py",
    "scripts/run_fgc_tdg10_qa2.py",
    "src/recursive_horizons/fgc/evolution/tdg10_qa2_authority.py",
    "tests/test_check_repo_tdg10_qa2_frz1.py",
    "tests/test_fgc_tdg10_qa2_authority.py",
    "tests/test_fgc_tdg10_qa2_runner.py",
)


class QA2RepositoryCheckTests(unittest.TestCase):
    def test_qa2_bundle_is_the_expected_sixteen_path_prelaunch_set(self) -> None:
        self.assertEqual(check_repo._TDG10_QA2_FRZ1_BUNDLE, _EXPECTED_BUNDLE)
        self.assertEqual(len(check_repo._TDG10_QA2_FRZ1_BUNDLE), 16)
        self.assertNotIn(".qdrant-initialized", check_repo._TDG10_QA2_FRZ1_BUNDLE)
        self.assertNotIn("CHANGELOG.md", check_repo._TDG10_QA2_FRZ1_BUNDLE)
        self.assertFalse(
            any(
                path.startswith("runs/") or path.startswith("paper/")
                for path in check_repo._TDG10_QA2_FRZ1_BUNDLE
            )
        )

    def test_qa2_calls_only_its_compact_reproducer(self) -> None:
        completed = SimpleNamespace(returncode=0, stdout="", stderr="")
        with patch.object(check_repo.subprocess, "run", return_value=completed) as run:
            self.assertEqual(check_repo._check_tdg10_qa2_frz1_result(), [])
        commands = [tuple(call.args[0]) for call in run.call_args_list]
        self.assertEqual(
            commands,
            [
                (
                    check_repo.sys.executable,
                    "scripts/reproduce_fgc_tdg10_qa2_frz1.py",
                    "--verify-compact",
                )
            ],
        )
        self.assertFalse(
            any(
                "run_fgc_tdg10_qa2" in item or "runs/" in item
                for command in commands
                for item in command
            )
        )

    def test_partial_qa2_bundle_fails_without_reproducer(self) -> None:
        missing = ROOT / check_repo._TDG10_QA2_FRZ1_BUNDLE[-1]
        original = check_repo._probe_unique_regular_leaf

        def probe(relative: str, label: str) -> bool:
            if relative == str(missing.relative_to(ROOT)):
                return False
            return original(relative, label)

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", side_effect=probe),
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg10_qa2_frz1_result()
        self.assertEqual(len(failures), 1)
        self.assertIn("compact bundle is partial", failures[0])
        run.assert_not_called()

    def test_unsafe_qa2_bundle_fails_without_reproducer(self) -> None:
        original = check_repo._probe_unique_regular_leaf

        def probe(relative: str, label: str) -> bool:
            if relative == check_repo._TDG10_QA2_FRZ1_BUNDLE[0]:
                raise ValueError("symlink")
            return original(relative, label)

        with (
            patch.object(check_repo, "_probe_unique_regular_leaf", side_effect=probe),
            patch.object(check_repo.subprocess, "run") as run,
        ):
            failures = check_repo._check_tdg10_qa2_frz1_result()
        self.assertEqual(len(failures), 1)
        self.assertIn("bundle path is unsafe", failures[0])
        run.assert_not_called()

    def test_compact_audit_passes_without_store_or_shadow_access(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                "scripts/check_repo.py",
                "--only-tdg10-qa2-frz1",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            completed.returncode, 0, completed.stdout + completed.stderr
        )


if __name__ == "__main__":
    unittest.main()
