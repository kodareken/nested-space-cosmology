"""Repository-checker boundaries for the one-time PRO19 SID1 anchor."""
from __future__ import annotations

from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

from scripts import check_repo


class PRO19SID1RepositoryCheckTests(unittest.TestCase):
    def test_compact_result_never_reopens_generation_seven_anchor(self) -> None:
        with patch.object(
            check_repo.subprocess,
            "run",
            side_effect=AssertionError("ordinary compact verification reopened SID1"),
        ):
            self.assertEqual(check_repo._check_pro19_sid1_compact_result(), [])

    def test_explicit_prelaunch_owns_the_live_reproducer(self) -> None:
        completed = subprocess.CompletedProcess(
            args=["sid1"], returncode=0, stdout="", stderr=""
        )
        with patch.object(
            check_repo,
            "_check_pro19_prelaunch_result",
            return_value=[],
        ), patch.object(
            check_repo,
            "_check_pro19_sid1_compact_result",
            return_value=[],
        ), patch.object(check_repo.subprocess, "run", return_value=completed) as run:
            self.assertEqual(check_repo._check_pro19_sid1_prelaunch_result(), [])
        command = run.call_args.args[0]
        self.assertEqual(
            command[1:],
            ["scripts/reproduce_fgc_pro19_sid1_frz1.py", "--verify"],
        )

    def test_launch_manifest_closure_requires_sid1_authority_paths(self) -> None:
        required = {
            "results/fgc-1-pro19-sid1-frz1.json",
            "configs/fgc/fgc-1-pro19-sid1-frz1.toml",
            "docs/fgc-pro19-sid1-frz1.md",
            "src/recursive_horizons/fgc/evolution/proto19_sid1_authority.py",
            "scripts/reproduce_fgc_pro19_sid1_frz1.py",
            "tests/test_fgc_proto19_sid1_authority.py",
            "tests/test_check_repo_pro19_sid1_prelaunch.py",
        }
        manifest_path = (
            Path(check_repo.REPOSITORY)
            / "configs/fgc/fgc-1-pro19-launch-authority.toml"
        )
        self.assertTrue(manifest_path.exists())
        # The dedicated checker enforces these through its required closure;
        # this assertion keeps the test's own fixture contract explicit.
        source = (
            Path(check_repo.REPOSITORY) / "scripts/repo_checks/_shared.py"
        ).read_text(encoding="utf-8")
        for path in required:
            with self.subTest(path=path):
                self.assertIn(f'"{path}"', source)


if __name__ == "__main__":
    unittest.main()
