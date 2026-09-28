from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt12_mon12 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    EXPECTED_CLAIMS,
    EXPECTED_RUNTIME,
    _canonical,
    load_canonical_result,
    load_config,
    validate_run_plan,
    verify_canonical,
)


class FGCHLT12MON12ReproductionTests(unittest.TestCase):
    def test_frozen_authorization_scope_is_narrow(self) -> None:
        config = load_config(DEFAULT_CONFIG)
        self.assertEqual(config["claims"], EXPECTED_CLAIMS)
        self.assertTrue(config["claims"]["PROTO14_successor_runtime_implemented"])
        self.assertTrue(config["claims"]["PROTO14_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(config["claims"]["GR0_case_eligible"])
        self.assertFalse(config["claims"]["SGBL_execution_authorized"])
        self.assertFalse(config["claims"]["FGCQR_holdout_execution_authorized"])
        self.assertEqual(config["numerical_runtime"], EXPECTED_RUNTIME)
        self.assertEqual(
            validate_run_plan(REPOSITORY / config["run_plan_config"])["artifact_id"],
            "FGC-1-CAL11-RUN1-PLAN",
        )
        self.assertEqual(config["restart_admission"]["expected_member_count"], 6)
        self.assertTrue(config["restart_admission"]["TDG6_ledger_must_start_zero_at_restart"])

    def test_claim_and_runtime_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "FGCQR_holdout_execution_authorized = false",
                    "FGCQR_holdout_execution_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "frozen scope"):
                load_config(promoted)
            wrong_runtime = root / "wrong-runtime.toml"
            wrong_runtime.write_text(
                source.replace('python_version = "3.14.3"', 'python_version = "3.14.6"'),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "frozen scope"):
                load_config(wrong_runtime)
            plan_source = (REPOSITORY / "configs/fgc/fgc-1-pro14-run1.toml").read_text(
                encoding="utf-8"
            )
            wrong_plan = root / "wrong-plan.toml"
            wrong_plan.write_text(
                plan_source.replace(
                    'minimum_observed_order = "3/2"',
                    'minimum_observed_order = "149/100"',
                    1,
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "frozen scope"):
                validate_run_plan(wrong_plan)

    def test_canonical_result_is_checked_when_sealed(self) -> None:
        """Keep this test valid before the coordinator emits the hash-bound result."""

        if not DEFAULT_OUTPUT.exists():
            self.assertFalse(DEFAULT_OUTPUT.exists())
            return
        observed = load_canonical_result(DEFAULT_OUTPUT)
        self.assertEqual(observed["gate_status"], EXPECTED_CLAIMS)
        self.assertFalse(observed["artifact_payload"]["epistemic_boundary"]["trajectory_advanced"])
        promoted = deepcopy(observed)
        promoted["gate_status"]["FGCQR_holdout_execution_authorized"] = True
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "promoted.json"
            path.write_text(_canonical(promoted), encoding="utf-8")
            self.assertTrue(
                load_canonical_result(path)["gate_status"][
                    "FGCQR_holdout_execution_authorized"
                ]
            )
            with self.assertRaisesRegex(ValueError, "fresh reproduction"):
                verify_canonical(path, DEFAULT_CONFIG)


if __name__ == "__main__":
    unittest.main()
