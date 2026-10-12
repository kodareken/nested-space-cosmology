"""Repository-check routing for the two-phase PRO19 SID3 authority bundle."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts import check_repo


class PRO19SID3RepositoryCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.repository = Path(self.temporary.name)

    def _write_bundle(self, *relative_paths: str) -> None:
        for relative in relative_paths:
            path = self.repository / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"placeholder\n")

    def test_absent_bundle_is_accepted_by_ordinary_compact_route(self) -> None:
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError("absent bundle spawned compact reproducer"),
            ),
        ):
            self.assertEqual(
                check_repo._check_pro19_sid3_compact_result(required=False),
                [],
            )

    def test_focused_route_rejects_absent_bundle(self) -> None:
        with patch.object(check_repo, "REPOSITORY", self.repository):
            self.assertEqual(
                check_repo._check_pro19_sid3_compact_result(required=True),
                ["FGC-1-PRO19-SID3-AUTH1 compact bundle is absent"],
            )

    def test_partial_bundle_is_rejected(self) -> None:
        self._write_bundle(check_repo._PRO19_SID3_BUNDLE[0])
        with patch.object(check_repo, "REPOSITORY", self.repository):
            findings = check_repo._check_pro19_sid3_compact_result()
        self.assertEqual(len(findings), 1)
        self.assertIn("compact bundle is partial", findings[0])
        for relative in check_repo._PRO19_SID3_BUNDLE[1:]:
            self.assertIn(relative, findings[0])

    def test_complete_bundle_routes_to_compact_reproducer(self) -> None:
        self._write_bundle(*check_repo._PRO19_SID3_BUNDLE)
        completed = subprocess.CompletedProcess(
            args=[], returncode=0, stdout="", stderr=""
        )
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.subprocess,
                "run",
                return_value=completed,
            ) as run,
        ):
            self.assertEqual(
                check_repo._check_pro19_sid3_compact_result(required=True),
                [],
            )
        run.assert_called_once_with(
            [
                check_repo.sys.executable,
                "scripts/reproduce_fgc_pro19_sid3_auth1.py",
                "--verify",
            ],
            cwd=self.repository,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    def test_compact_route_never_opens_mutable_store(self) -> None:
        self._write_bundle(*check_repo._PRO19_SID3_BUNDLE)
        completed = subprocess.CompletedProcess(
            args=[], returncode=0, stdout="", stderr=""
        )
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.os,
                "scandir",
                side_effect=AssertionError("SID3 compact route scanned the run store"),
            ),
            patch.object(
                check_repo,
                "_sid2_live_store_inventory",
                side_effect=AssertionError(
                    "SID3 compact route opened generation eight"
                ),
            ),
            patch.object(check_repo.subprocess, "run", return_value=completed),
        ):
            self.assertEqual(check_repo._check_pro19_sid3_compact_result(), [])

    def test_explicit_sid3_flag_requires_the_bundle(self) -> None:
        standard_checks = (
            "check_required",
            "check_git_data_boundary",
        )
        patches = [
            patch.object(check_repo, name, return_value=[]) for name in standard_checks
        ]
        mocks = [item.start() for item in patches]
        self.addCleanup(lambda: [item.stop() for item in reversed(patches)])
        focused = patch.object(
            check_repo,
            "_check_pro19_sid3_compact_result",
            return_value=[],
        )
        focused_mock = focused.start()
        self.addCleanup(focused.stop)

        with redirect_stdout(StringIO()):
            self.assertEqual(
                check_repo.main(["--only-pro19-sid3-prelaunch"]),
                0,
            )
        for check in mocks:
            check.assert_called_once_with()
        focused_mock.assert_called_once_with(required=True)

    def test_sid3_flag_is_mutually_exclusive(self) -> None:
        with (
            redirect_stdout(StringIO()),
            redirect_stderr(StringIO()),
            self.assertRaises(SystemExit) as raised,
        ):
            check_repo.main(
                [
                    "--only-pro19-sid2-postrecovery",
                    "--only-pro19-sid3-prelaunch",
                ]
            )
        self.assertEqual(raised.exception.code, 2)

    def test_required_closure_contains_sid3_python_owners_only(self) -> None:
        required = set(check_repo.REQUIRED)
        self.assertIn(
            "src/recursive_horizons/fgc/evolution/proto19_resume_authority.py",
            required,
        )
        self.assertIn("scripts/reproduce_fgc_pro19_sid3_auth1.py", required)
        self.assertIn("tests/test_check_repo_pro19_sid3_prelaunch.py", required)
        self.assertTrue(set(check_repo._PRO19_SID3_BUNDLE).isdisjoint(required))

    def test_compact_reproducer_has_no_run_store_literal(self) -> None:
        source = (
            Path(__file__).resolve().parents[1]
            / "scripts/reproduce_fgc_pro19_sid3_auth1.py"
        ).read_text("utf-8")
        self.assertNotIn("runs/fgc-2-sf1", source)
        self.assertNotIn("STORE_PATH", source)


if __name__ == "__main__":
    unittest.main()
