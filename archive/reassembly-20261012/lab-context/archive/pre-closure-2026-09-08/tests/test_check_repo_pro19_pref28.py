"""Repository-check boundaries for the outcome-neutral PRO19 PREF28 binder."""

from __future__ import annotations

import ast
from contextlib import redirect_stderr, redirect_stdout
from hashlib import sha256
from io import StringIO
import inspect
import json
from pathlib import Path
import subprocess
import tempfile
import textwrap
import tomllib
import unittest
from unittest.mock import patch

from scripts import check_repo
from scripts import reproduce_fgc_pro19_pref28 as pref28_reproducer


class PRO19PREF28RepositoryCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.repository = Path(self.temporary.name)

    def _write_bundle(self, *relative_paths: str) -> None:
        for relative in relative_paths:
            path = self.repository / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"placeholder\n")

    def _copy_postterminal_inputs(self, repository: Path) -> None:
        source = Path(__file__).resolve().parents[1]
        historical = {
            "configs/fgc/fgc-1-pro19-sid3-auth1.toml",
            "configs/fgc/fgc-1-pro19-sid3-execution-closure.json",
            "results/fgc-1-pro19-sid3-auth1.json",
            "configs/fgc/fgc-1-hlt16-mon16.toml",
            "results/fgc-1-hlt16-mon16.json",
            "configs/fgc/fgc-1-pro19-launch-authority.toml",
        }
        for relative in sorted(set(check_repo._PRO19_PREF28_BUNDLE) | historical):
            destination = repository / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((source / relative).read_bytes())

    def test_complete_bundle_uses_store_blind_compact_reproducer(self) -> None:
        self._write_bundle(*check_repo._PRO19_PREF28_BUNDLE)
        completed = subprocess.CompletedProcess(
            args=[], returncode=0, stdout="", stderr=""
        )
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.os,
                "scandir",
                side_effect=AssertionError("PREF28 compact route scanned runs"),
            ),
            patch.object(
                check_repo,
                "_check_pro19_pref28_postattempt_result",
                side_effect=AssertionError(
                    "PREF28 compact route used live verification"
                ),
            ),
            patch.object(
                check_repo.subprocess,
                "run",
                return_value=completed,
            ) as run,
        ):
            self.assertEqual(check_repo._check_pro19_pref28_compact_result(), [])

        run.assert_called_once_with(
            [
                check_repo.sys.executable,
                "scripts/reproduce_fgc_pro19_pref28.py",
                "--verify-compact",
            ],
            cwd=self.repository,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    def test_default_result_route_registers_compact_but_not_live_check(self) -> None:
        source = textwrap.dedent(inspect.getsource(check_repo.check_results))
        tree = ast.parse(source)
        calls = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }
        self.assertIn("_check_pro19_pref28_compact_result", calls)
        self.assertNotIn("_check_pro19_pref28_postattempt_result", calls)

    def test_compact_make_target_has_no_historical_live_dependency(self) -> None:
        makefile = (
            Path(__file__).resolve().parents[1]
            / "mk"
            / "historical-certificates.mk"
        ).read_text("utf-8")
        self.assertIn("\nfgc-pro19-pref28:\n", makefile)
        self.assertNotIn("fgc-pro19-pref28: fgc-pro19-sid3-auth1", makefile)

    def test_partial_bundle_is_rejected_without_reproduction(self) -> None:
        self.assertGreater(len(check_repo._PRO19_PREF28_BUNDLE), 1)
        self._write_bundle(check_repo._PRO19_PREF28_BUNDLE[0])
        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError("partial PREF28 bundle was reproduced"),
            ),
        ):
            findings = check_repo._check_pro19_pref28_compact_result()

        self.assertEqual(len(findings), 1)
        self.assertIn("compact bundle is partial", findings[0])
        for relative in check_repo._PRO19_PREF28_BUNDLE[1:]:
            self.assertIn(relative, findings[0])

    def test_unsafe_bundle_path_is_rejected_without_reproduction(self) -> None:
        relative = check_repo._PRO19_PREF28_BUNDLE[0]
        target = self.repository / "unsafe-pref28-target"
        target.write_bytes(b"not a bundle leaf\n")
        path = self.repository / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.symlink_to(target)

        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError("unsafe PREF28 bundle was reproduced"),
            ),
        ):
            findings = check_repo._check_pro19_pref28_compact_result()

        self.assertEqual(len(findings), 1)
        self.assertIn("bundle path is unsafe", findings[0])
        self.assertIn(relative, findings[0])

    def test_bundle_probe_rejects_inner_directory_symlink(self) -> None:
        self._write_bundle(*check_repo._PRO19_PREF28_BUNDLE)
        configs = self.repository / "configs"
        redirected = self.repository / "redirected-configs"
        configs.rename(redirected)
        configs.symlink_to(redirected, target_is_directory=True)

        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError(
                    "inner-directory symlink reached PREF28 reproduction"
                ),
            ),
        ):
            findings = check_repo._check_pro19_pref28_compact_result()

        self.assertEqual(len(findings), 1)
        self.assertIn("bundle path is unsafe", findings[0])
        self.assertIn("configs/fgc/fgc-1-pro19-pref28.toml", findings[0])

    def test_historical_reader_rejects_inner_directory_symlink(self) -> None:
        redirected = self.repository / "redirected-configs"
        leaf = redirected / "fgc/authority.toml"
        leaf.parent.mkdir(parents=True)
        leaf.write_bytes(b"authority\n")
        (self.repository / "configs").symlink_to(redirected, target_is_directory=True)

        with patch.object(check_repo, "REPOSITORY", self.repository):
            with self.assertRaisesRegex(ValueError, "ancestor is unsafe"):
                check_repo._read_unique_regular_bytes(
                    "configs/fgc/authority.toml", "historical authority fixture"
                )

    def test_direct_reproducer_rejects_unsafe_config_and_result_paths(self) -> None:
        for case in ("inner_config_directory", "result_leaf_symlink"):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as root:
                repository = Path(root)
                config = repository / pref28_reproducer.CONFIG_PATH
                result = repository / pref28_reproducer.RESULT_PATH
                if case == "inner_config_directory":
                    expected_error = "ancestor is unsafe"
                    redirected = repository / "redirected-configs"
                    redirected_config = redirected / Path(
                        pref28_reproducer.CONFIG_PATH
                    ).relative_to("configs")
                    redirected_config.parent.mkdir(parents=True)
                    redirected_config.write_bytes(b"unsafe config target\n")
                    (repository / "configs").symlink_to(
                        redirected, target_is_directory=True
                    )
                    result.parent.mkdir(parents=True)
                    result.write_bytes(b"unused result\n")
                else:
                    expected_error = "not a unique bounded regular file"
                    config.parent.mkdir(parents=True)
                    config.write_bytes(b"unused config\n")
                    result.parent.mkdir(parents=True)
                    target = repository / "unsafe-result-target"
                    target.write_bytes(b"unsafe result target\n")
                    result.symlink_to(target)

                with (
                    patch.object(pref28_reproducer, "ROOT", repository),
                    patch.object(
                        pref28_reproducer,
                        "validate_compact_result",
                        side_effect=AssertionError(
                            "unsafe direct reproducer input reached validation"
                        ),
                    ),
                    self.assertRaisesRegex(ValueError, expected_error),
                ):
                    pref28_reproducer.main(["--verify-compact"])

    def test_each_single_historical_bundle_leaf_is_partial_without_fallback(
        self,
    ) -> None:
        for relative in check_repo._PRO19_PREF28_BUNDLE:
            with self.subTest(relative=relative), tempfile.TemporaryDirectory() as root:
                repository = Path(root)
                path = repository / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"one historical PREF28 leaf\n")
                with (
                    patch.object(check_repo, "REPOSITORY", repository),
                    patch.object(
                        check_repo,
                        "_check_pro19_pref28_compact_result",
                        side_effect=AssertionError(
                            "partial historical bundle fell through to compact checking"
                        ),
                    ),
                ):
                    anchor, findings = check_repo._pref28_historical_prelaunch_anchor()

                self.assertIsNone(anchor)
                self.assertEqual(len(findings), 1)
                self.assertIn("historical prelaunch anchor is partial", findings[0])
                for missing in check_repo._PRO19_PREF28_BUNDLE:
                    if missing != relative:
                        self.assertIn(missing, findings[0])

    def test_consumed_historical_leaf_byte_tampering_is_rejected(self) -> None:
        cases = (
            (
                "sid3_config",
                "configs/fgc/fgc-1-pro19-sid3-auth1.toml",
                check_repo._check_hlt16_mon16_result,
                "FGC-1-PRO19-PREF28 historical authority differs: "
                "configs/fgc/fgc-1-pro19-sid3-auth1.toml",
            ),
            (
                "sid3_closure",
                "configs/fgc/fgc-1-pro19-sid3-execution-closure.json",
                check_repo._check_hlt16_mon16_result,
                "FGC-1-PRO19-PREF28 historical authority differs: "
                "configs/fgc/fgc-1-pro19-sid3-execution-closure.json",
            ),
            (
                "sid3_result",
                "results/fgc-1-pro19-sid3-auth1.json",
                check_repo._check_hlt16_mon16_result,
                "FGC-1-PRO19-PREF28 historical authority differs: "
                "results/fgc-1-pro19-sid3-auth1.json",
            ),
            (
                "hlt16_config",
                "configs/fgc/fgc-1-hlt16-mon16.toml",
                check_repo._check_hlt16_mon16_result,
                "FGC-1-PRO19-PREF28 historical authority differs: "
                "configs/fgc/fgc-1-hlt16-mon16.toml",
            ),
            (
                "hlt16_result",
                "results/fgc-1-hlt16-mon16.json",
                check_repo._check_hlt16_mon16_result,
                "FGC-1-PRO19-PREF28 historical authority differs: "
                "results/fgc-1-hlt16-mon16.json",
            ),
            (
                "launch_manifest",
                "configs/fgc/fgc-1-pro19-launch-authority.toml",
                check_repo._check_pro19_prelaunch_result,
                "FGC-1-PRO19-PREF28 historical authority differs: "
                "configs/fgc/fgc-1-pro19-launch-authority.toml",
            ),
        )
        for name, relative, checker, expected in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as root:
                repository = Path(root)
                self._copy_postterminal_inputs(repository)
                path = repository / relative
                if name in {"sid3_closure", "sid3_result", "hlt16_result"}:
                    record = json.loads(path.read_text("utf-8"))
                    record["exact_byte_tamper_probe"] = True
                    path.write_text(
                        json.dumps(
                            record,
                            sort_keys=True,
                            indent=2,
                            ensure_ascii=True,
                            allow_nan=False,
                        )
                        + "\n",
                        encoding="utf-8",
                    )
                else:
                    path.write_bytes(path.read_bytes() + b"\n# exact-byte tamper\n")

                with (
                    patch.object(check_repo, "REPOSITORY", repository),
                    patch.object(
                        check_repo,
                        "_check_pro19_pref28_compact_result",
                        return_value=[],
                    ),
                    patch.object(
                        check_repo,
                        "_check_pro19_frz1_result",
                        return_value=[],
                    ),
                    patch.object(
                        check_repo,
                        "_check_hlt16_mon16_result",
                        return_value=[],
                    )
                    if name == "launch_manifest"
                    else patch.object(
                        check_repo.subprocess,
                        "run",
                        side_effect=AssertionError(
                            "tampered historical leaf fell back to a reproducer"
                        ),
                    ),
                ):
                    findings = checker()

                self.assertEqual(findings, [expected])

    def test_valid_postterminal_mode_skips_both_legacy_reproducers(self) -> None:
        self._copy_postterminal_inputs(self.repository)
        source = Path(__file__).resolve().parents[1]
        current_checker_sha = sha256(
            (source / "scripts/check_repo.py").read_bytes()
        ).hexdigest()
        hlt16 = json.loads(
            (self.repository / "results/fgc-1-hlt16-mon16.json").read_text("utf-8")
        )
        manifest = tomllib.loads(
            (
                self.repository / "configs/fgc/fgc-1-pro19-launch-authority.toml"
            ).read_text("utf-8")
        )
        hlt16_checker_sha = next(
            item["sha256"]
            for item in hlt16["artifact_payload"]["inventory"]
            if item["path"] == "scripts/check_repo.py"
        )
        manifest_checker_sha = next(
            item["sha256"]
            for item in manifest["authority_path"]
            if item["path"] == "scripts/check_repo.py"
        )
        self.assertNotEqual(current_checker_sha, hlt16_checker_sha)
        self.assertNotEqual(current_checker_sha, manifest_checker_sha)

        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo,
                "_check_pro19_pref28_compact_result",
                return_value=[],
            ),
            patch.object(check_repo, "_check_pro19_frz1_result", return_value=[]),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError(
                    "postterminal mode used a legacy reproducer"
                ),
            ) as run,
        ):
            self.assertEqual(check_repo._check_hlt16_mon16_result(), [])
            self.assertEqual(check_repo._check_pro19_prelaunch_result(), [])
        run.assert_not_called()

    def test_postterminal_ordinary_checks_are_git_and_runs_blind(self) -> None:
        self._copy_postterminal_inputs(self.repository)
        self.assertFalse((self.repository / ".git").exists())
        self.assertFalse((self.repository / "runs").exists())

        with (
            patch.object(check_repo, "REPOSITORY", self.repository),
            patch.object(
                check_repo,
                "_check_pro19_pref28_compact_result",
                return_value=[],
            ),
            patch.object(check_repo, "_check_pro19_frz1_result", return_value=[]),
            patch.object(
                check_repo.os,
                "scandir",
                side_effect=AssertionError("postterminal ordinary check scanned runs"),
            ),
            patch.object(
                check_repo.subprocess,
                "run",
                side_effect=AssertionError("postterminal ordinary check opened Git"),
            ),
        ):
            self.assertEqual(check_repo._check_hlt16_mon16_result(), [])
            self.assertEqual(check_repo._check_pro19_prelaunch_result(), [])

    def test_explicit_postattempt_flag_routes_through_live_gate(self) -> None:
        standard_checks = (
            "check_required",
            "check_links",
            "check_public_text",
            "check_fgc_alignment",
            "check_git_data_boundary",
            "check_reference_sequence",
        )
        patches = [
            patch.object(check_repo, name, return_value=[]) for name in standard_checks
        ]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)

        live = patch.object(
            check_repo,
            "_check_pro19_pref28_postattempt_result",
            return_value=[],
        )
        live_mock = live.start()
        self.addCleanup(live.stop)
        default = patch.object(
            check_repo,
            "check_results",
            side_effect=AssertionError("focused PREF28 routing used default results"),
        )
        default.start()
        self.addCleanup(default.stop)

        with redirect_stdout(StringIO()):
            self.assertEqual(
                check_repo.main(["--only-pro19-pref28-postattempt"]),
                0,
            )
        live_mock.assert_called_once_with()

    def test_pref28_flag_is_mutually_exclusive(self) -> None:
        with (
            redirect_stdout(StringIO()),
            redirect_stderr(StringIO()),
            self.assertRaises(SystemExit) as raised,
        ):
            check_repo.main(
                [
                    "--only-pro19-sid3-prelaunch",
                    "--only-pro19-pref28-postattempt",
                ]
            )
        self.assertEqual(raised.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
