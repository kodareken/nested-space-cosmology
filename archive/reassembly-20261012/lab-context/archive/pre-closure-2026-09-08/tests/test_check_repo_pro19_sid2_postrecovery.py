"""Repository-checker boundaries for the one-time PRO19 SID2 binder."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import unittest
from unittest.mock import patch

from scripts import check_repo


class PRO19SID2RepositoryCheckTests(unittest.TestCase):
    def test_compact_result_does_not_open_store_or_spawn_reproducer(self) -> None:
        with patch.object(
            check_repo,
            "_sid2_live_store_inventory",
            side_effect=AssertionError("compact verification opened the live store"),
        ), patch.object(
            check_repo.os,
            "scandir",
            side_effect=AssertionError("compact verification scanned the live store"),
        ), patch.object(
            check_repo.subprocess,
            "run",
            side_effect=AssertionError("compact verification spawned a process"),
        ):
            self.assertEqual(check_repo._check_pro19_sid2_compact_result(), [])

    def test_explicit_postrecovery_flag_routes_through_live_sid2_gate(self) -> None:
        standard_checks = (
            "check_required",
            "check_links",
            "check_public_text",
            "check_fgc_alignment",
            "check_git_data_boundary",
            "check_reference_sequence",
        )
        patches = [
            patch.object(check_repo, name, return_value=[])
            for name in standard_checks
        ]
        live = patch.object(
            check_repo,
            "_check_pro19_sid2_postrecovery_result",
            return_value=[],
        )
        default = patch.object(
            check_repo,
            "check_results",
            side_effect=AssertionError("focused SID2 routing used default results"),
        )
        mocks = [item.start() for item in patches]
        live_mock = live.start()
        default.start()
        self.addCleanup(default.stop)
        self.addCleanup(live.stop)
        for item in reversed(patches):
            self.addCleanup(item.stop)

        with redirect_stdout(StringIO()):
            self.assertEqual(
                check_repo.main(["--only-pro19-sid2-postrecovery"]),
                0,
            )
        for check in mocks:
            check.assert_called_once_with()
        live_mock.assert_called_once_with()

    def test_required_closure_contains_every_sid2_owner(self) -> None:
        required = {
            "results/fgc-1-pro19-sid2-pref1.json",
            "configs/fgc/fgc-1-pro19-sid2-pref1.toml",
            "docs/fgc-pro19-sid2-pref1.md",
            "src/recursive_horizons/fgc/evolution/proto19_sid2_binder.py",
            "scripts/reproduce_fgc_pro19_sid2_pref1.py",
            "tests/test_fgc_proto19_sid2_binder.py",
            "tests/test_check_repo_pro19_sid2_postrecovery.py",
        }
        self.assertLessEqual(required, set(check_repo.REQUIRED))

    def test_sid2_flag_is_mutually_exclusive(self) -> None:
        with redirect_stdout(StringIO()), redirect_stderr(StringIO()), self.assertRaises(
            SystemExit
        ) as raised:
            check_repo.main(
                [
                    "--only-pro19-sid1-prelaunch",
                    "--only-pro19-sid2-postrecovery",
                ]
            )
        self.assertEqual(raised.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
