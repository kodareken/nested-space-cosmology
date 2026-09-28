"""Phase -1 Make/search/test routing invariants."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import re
import subprocess
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[1]
MAKEFILE = ROOT / "Makefile"
HISTORICAL = ROOT / "mk/historical-certificates.mk"
CURRENT = ROOT / "mk/current-foundation.mk"
CLOSED = ROOT / "mk/closed-live-targets.mk"
ACTIVE_MAP = ROOT / "docs/active-code-map.md"
MOVED_TEST = ROOT / "archive/historical-tests/test_fgc_pro19_sid3_real_store_preflight.py"
MOVED_TEST_SHA256 = "09d3515947530384063ad83c11faaea25ef7d43d527f23b15b50597067e63b3e"

EXPECTED_CLOSED = (
    "run-fgc-pro13-calibration",
    "run-fgc-pro14-calibration",
    "recover-fgc-pro19-sid1",
    "resume-fgc-pro19-sid3-event1",
    "run-fgc-pro19-event1",
    "resume-fgc-pro19-event1",
    "run-fgc-tdg8-successor-event1",
    "recover-fgc-tdg8-rcv1",
    "resume-fgc-tdg8-rcv1",
    "run-fgc-tdg8-rcv2",
    "run-fgc-tdg8-rcv3",
    "run-fgc-tdg8-rcv3-rec1",
    "run-fgc-tdg9-ar1",
    "run-fgc-tdg9-loc1",
    "run-fgc-tdg9-loc2",
    "run-fgc-tdg9-ti1",
    "run-fgc-tdg9-ti2",
    "run-fgc-tdg9-ac1",
    "run-fgc-tdg9-ur1",
    "run-fgc-tdg10-qa1",
    "run-fgc-tdg10-qa2",
    "run-fgc-tdg11-msel1",
    "run-fgc-hlt17-srcq1-rec1",
    "run-fgc-pro20-ev1",
)


def _closed_targets(source: str) -> tuple[str, ...]:
    match = re.search(
        r"CLOSED_HISTORICAL_STATE_TARGETS := \\\n(?P<body>(?:\t[^\n]+(?: \\)?\n)+)",
        source,
    )
    if match is None:
        raise AssertionError("closed target variable is absent")
    return tuple(
        line.strip().removesuffix(" \\")
        for line in match.group("body").splitlines()
        if line.strip()
    )


class PhaseMinusOneMakeRoutingTests(unittest.TestCase):
    def test_root_makefile_is_only_the_three_include_router(self) -> None:
        source = MAKEFILE.read_text(encoding="utf-8")
        self.assertLessEqual(len(source.splitlines()), 20)
        self.assertIn(".DEFAULT_GOAL := test", source)
        self.assertEqual(
            [line for line in source.splitlines() if line.startswith("include ")],
            [
                "include mk/historical-certificates.mk",
                "include mk/current-foundation.mk",
                "include mk/closed-live-targets.mk",
            ],
        )
        for path in (HISTORICAL, CURRENT, CLOSED):
            self.assertTrue(path.is_file())

    def test_exact_closed_state_changing_targets_are_tombstoned(self) -> None:
        closed = CLOSED.read_text(encoding="utf-8")
        self.assertEqual(_closed_targets(closed), EXPECTED_CLOSED)
        recipe = closed.split("$(CLOSED_HISTORICAL_STATE_TARGETS):", 1)[1]
        self.assertIn("@printf", recipe)
        self.assertIn("@exit 2", recipe)
        self.assertNotIn("PYTHON", recipe)
        historical = HISTORICAL.read_text(encoding="utf-8")
        for target in EXPECTED_CLOSED:
            self.assertNotRegex(historical, rf"(?m)^{re.escape(target)}\s*:")

    def test_verify_fgc_routes_to_current_foundation(self) -> None:
        source = CURRENT.read_text(encoding="utf-8")
        self.assertIn("verify-fgc: verify-fgc-sf1-foundation", source)
        self.assertIn("verify-current-development:", source)
        for focused in (
            "test_phase_minus1_artifact_catalog.py",
            "test_phase_minus1_check_repo_dispatcher.py",
            "test_phase_minus1_make_routing.py",
            "test_phase_minus1_public_authority.py",
            "test_historical_runner_import_boundary.py",
        ):
            self.assertIn(focused, source)
        self.assertIn("scripts/build_artifact_catalog.py --check", source)
        self.assertIn("verify-fgc-pro20-ev1-components", source)
        self.assertIn("verify-fgc-pro20-rsrc1-components", source)
        self.assertIn("tests/test_fgc_pro20_ev1_protocol.py", source)
        self.assertIn("tests/test_fgc_pro20_ev1_store.py", source)
        self.assertIn("tests/test_fgc_pro20_ev1_runtime.py", source)
        self.assertIn("tests/test_fgc_pro20_ev1_authority.py", source)
        self.assertIn("tests/test_fgc_pro20_rsrc1_member.py", source)
        self.assertIn("tests/test_fgc_pro20_rsrc1_attempt.py", source)
        self.assertIn("tests/test_fgc_pro20_rsrc1_isolation.py", source)
        self.assertIn("tests/test_fgc_pro20_rsrc1_seed.py", source)
        self.assertNotRegex(source, r"(?m)^run-fgc-pro20-rsrc1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-pro20-rsrc1:")
        self.assertRegex(source, r"(?m)^fgc-pro20-ev1-pref1:")
        self.assertRegex(source, r"(?m)^verify-fgc-pro20-ev1-pref1:")
        self.assertIn("scripts/reproduce_fgc_pro20_ev1_pref1.py --check", source)
        self.assertIn("tests/test_fgc_pro20_ev1_pref1_certificate.py", source)
        self.assertIn("tests/test_check_repo_pro20_ev1_pref1.py", source)
        self.assertRegex(source, r"(?m)^fgc-def1-stab1-frz1:")
        self.assertRegex(source, r"(?m)^verify-fgc-def1-stab1-frz1:")
        self.assertIn("scripts/reproduce_fgc_def1_stab1_frz1.py --check", source)
        self.assertIn("tests/test_fgc_def1_stab1_frz1_certificate.py", source)
        self.assertIn("tests/test_check_repo_def1_stab1_frz1.py", source)
        self.assertRegex(source, r"(?m)^fgc-def1-stab1-pref1:")
        self.assertRegex(source, r"(?m)^verify-fgc-def1-stab1-pref1:")
        self.assertIn("scripts/reproduce_fgc_def1_stab1_pref1.py --check", source)
        self.assertIn("tests/test_fgc_def1_stab1_pref1_binder.py", source)
        self.assertIn("tests/test_check_repo_def1_stab1_pref1.py", source)
        self.assertRegex(source, r"(?m)^fgc-sgb1-ctl1-sol1-frz1:")
        self.assertRegex(source, r"(?m)^verify-fgc-sgb1-ctl1-sol1-frz1:")
        self.assertIn("scripts/reproduce_fgc_sgb1_ctl1_sol1_frz1.py --check", source)
        self.assertIn("tests/test_fgc_sgb1_ctl1_sol1_frz1_certificate.py", source)
        self.assertIn("tests/test_check_repo_sgb1_ctl1_sol1_frz1.py", source)
        self.assertRegex(source, r"(?m)^fgc-sgb1-ctl1-sol1-pref1:")
        self.assertRegex(source, r"(?m)^verify-fgc-sgb1-ctl1-sol1-pref1:")
        self.assertIn("scripts/reproduce_fgc_sgb1_ctl1_sol1_pref1.py --check", source)
        self.assertIn("tests/test_fgc_sgb1_ctl1_sol1_pref1_binder.py", source)
        self.assertIn("tests/test_check_repo_sgb1_ctl1_sol1_pref1.py", source)
        foundation = next(
            line
            for line in source.splitlines()
            if line.startswith("verify-fgc-sf1-foundation:")
        )
        self.assertIn("fgc-def1-stab1-frz1", foundation)
        self.assertIn("fgc-def1-stab1-pref1", foundation)
        self.assertIn("fgc-sgb1-ctl1-sol1-frz1", foundation)
        self.assertIn("fgc-sgb1-ctl1-sol1-pref1", foundation)
        self.assertNotRegex(source, r"(?m)^run-fgc-def1-stab1-frz1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-def1-stab1-frz1:")
        self.assertNotRegex(source, r"(?m)^run-fgc-def1-stab1-pref1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-def1-stab1-pref1:")
        self.assertNotRegex(source, r"(?m)^run-fgc-sgb1-ctl1-sol1-frz1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-sgb1-ctl1-sol1-frz1:")
        self.assertNotRegex(source, r"(?m)^run-fgc-sgb1-ctl1-sol1-pref1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-sgb1-ctl1-sol1-pref1:")
        self.assertNotRegex(source, r"(?m)^run-fgc-pro20-ev1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-pro20-ev1:")
        self.assertNotIn("scripts/run_fgc_pro20_ev1.py --run", source)
        self.assertIn("run-fgc-pro20-ev1", EXPECTED_CLOSED)
        historical = HISTORICAL.read_text(encoding="utf-8")
        self.assertNotRegex(historical, r"(?m)^verify-fgc\s*:")
        self.assertNotIn("__retired-verify-fgc", historical)

    def test_active_command_map_routes_only_to_closed_pro20_checks(self) -> None:
        source = ACTIVE_MAP.read_text(encoding="utf-8")
        self.assertIn("make fgc-pro20-ev1-pref1", source)
        self.assertIn("make verify-fgc-pro20-ev1-pref1", source)
        self.assertIn("make verify-current-development", source)
        self.assertNotIn("make verify-fgc-pro20-ev1-prelaunch", source)
        self.assertNotIn("make status-fgc-pro20-ev1", source)
        self.assertIn("make run-fgc-pro20-ev1` is tombstoned", source)

    def test_make_database_exposes_tombstones_without_recipe_override(self) -> None:
        completed = subprocess.run(
            ["make", "-npRrq", "-f", str(MAKEFILE), ":"],
            cwd=ROOT,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertIn(completed.returncode, {0, 2})
        self.assertNotIn("overriding recipe", completed.stderr)
        self.assertNotIn("ignoring old recipe", completed.stderr)
        for target in EXPECTED_CLOSED:
            self.assertIn(f"{target}:", completed.stdout)

    def test_every_tombstone_refuses_before_python_or_prerequisites(self) -> None:
        for target in EXPECTED_CLOSED:
            with self.subTest(target=target):
                completed = subprocess.run(
                    ["make", "--no-print-directory", target],
                    cwd=ROOT,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False,
                )
                combined = completed.stdout + completed.stderr
                self.assertNotEqual(completed.returncode, 0)
                self.assertIn(
                    f"CLOSED historical operation '{target}' completed or terminal",
                    combined,
                )
                self.assertIn("state-changing rerun is not authorized", combined)
                self.assertIn("use the named compact verifier", combined)
                self.assertNotIn("python", combined.lower())
                self.assertNotIn("reproduce_", combined)
                self.assertNotIn("run_fgc_", combined)

    def test_moved_transient_test_is_hash_preserved_and_not_discovered(self) -> None:
        self.assertFalse((ROOT / "tests/test_fgc_pro19_sid3_real_store_preflight.py").exists())
        self.assertTrue(MOVED_TEST.is_file())
        self.assertEqual(sha256(MOVED_TEST.read_bytes()).hexdigest(), MOVED_TEST_SHA256)
        discovered = {path.name for path in (ROOT / "tests").glob("test*.py")}
        self.assertNotIn(MOVED_TEST.name, discovered)

    def test_search_and_tool_configs_exclude_noncanonical_trees(self) -> None:
        ignore = (ROOT / ".rgignore").read_text(encoding="utf-8")
        self.assertIn("archive/**", ignore)
        self.assertIn("runs/**", ignore)
        with (ROOT / "pyproject.toml").open("rb") as handle:
            config = tomllib.load(handle)
        pytest = config["tool"]["pytest"]["ini_options"]
        self.assertEqual(pytest["testpaths"], ["tests"])
        self.assertIn("archive", pytest["norecursedirs"])
        ruff = config["tool"]["ruff"]
        self.assertIn("archive", ruff["extend-exclude"])

    def test_current_foundation_owns_pref1_compact_routes_without_live_bind(self) -> None:
        source = CURRENT.read_text(encoding="utf-8")
        self.assertIn("fgc-tdg11-msel1-pref1", source)
        self.assertIn("verify-fgc-tdg11-msel1-pref1", source)
        self.assertIn("scripts/reproduce_fgc_tdg11_msel1_pref1.py --check", source)
        self.assertIn("tests/test_fgc_tdg11_msel1_pref1_binder.py", source)
        self.assertIn("tests/test_fgc_tdg11_msel1_pref1_protocol.py", source)
        self.assertIn("tests/test_fgc_tdg11_msel1_pref1_reconstruction.py", source)
        self.assertIn("tests/test_fgc_tdg11_msel1_pref1_localization.py", source)
        self.assertIn("tests/test_check_repo_tdg11_msel1_pref1.py", source)
        self.assertIn("scripts/check_repo.py --only-tdg11-msel1-pref1", source)
        self.assertIn("verify-fgc-sf1-foundation: fgc-tdg11-msel1-frz1 fgc-tdg11-msel1-pref1", source)
        self.assertIn("fgc-tdg11-msel1-pref1", source.split("CURRENT_DEVELOPMENT_CERTIFICATES", 1)[1].split("verify-current-development", 1)[0])
        self.assertNotRegex(source, r"(?m)^run-fgc-tdg11-msel1\s*:")
        self.assertNotIn("reproduce_fgc_tdg11_msel1_pref1.py --bind", source)
        self.assertNotIn("reproduce_fgc_tdg11_msel1_pref1.py --write-result", source)
        self.assertNotIn("bind-fgc-tdg11-msel1-pref1", source)
        self.assertIn("scripts/reproduce_fgc_tdg11_msel1_frz1.py --check", source)
        self.assertIn("scripts/run_fgc_tdg11_msel1.py --status", source)
        self.assertNotIn("scripts/run_fgc_tdg11_msel1.py --run", source)
        self.assertIn("fgc-pro20-ev1-pref1", source)
        self.assertIn("scripts/reproduce_fgc_pro20_ev1_pref1.py --check", source)
        self.assertIn("tests/test_fgc_pro20_ev1_pref1_certificate.py", source)
        self.assertNotIn("reproduce_fgc_pro20_ev1_pref1.py --live", source)
        self.assertIn("fgc-pro20-ev1-pref1", source.split("CURRENT_DEVELOPMENT_CERTIFICATES", 1)[1].split("verify-current-development", 1)[0])


if __name__ == "__main__":
    unittest.main()
