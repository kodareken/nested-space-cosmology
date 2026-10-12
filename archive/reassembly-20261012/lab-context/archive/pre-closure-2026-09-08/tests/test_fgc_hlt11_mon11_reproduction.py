from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt11_mon11 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    EXPECTED_CLAIMS,
    EXPECTED_NUMERICAL_RUNTIME,
    _canonical,
    load_canonical_result,
    load_config,
    record,
    validate_run_plan,
)


class FGCHLT11MON11ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.observed = load_canonical_result(DEFAULT_OUTPUT)
        historical = cls.observed["artifact_payload"]["namespace_precondition"]
        cls.reproduced = record(
            DEFAULT_CONFIG, historical_namespace_evidence=historical
        )

    def test_canonical_result_reproduces(self) -> None:
        self.assertEqual(self.observed, self.reproduced)

    def test_authorization_is_narrow(self) -> None:
        self.assertEqual(self.observed["gate_status"], EXPECTED_CLAIMS)
        self.assertTrue(
            self.observed["gate_status"][
                "PROTO13_fresh_GR0_dynamic_calibration_authorized"
            ]
        )
        self.assertFalse(
            self.observed["gate_status"]["FGCQR_holdout_execution_authorized"]
        )
        admission = self.observed["artifact_payload"]["restart_admission"]
        self.assertEqual(admission["restored_member_count"], 6)
        self.assertTrue(admission["source_checks_passed"])
        self.assertTrue(admission["primary_common_event"]["admission_passed"])
        self.assertTrue(admission["comparator_common_event"]["admission_passed"])
        self.assertFalse(admission["trajectory_advanced"])

    def test_method_owned_contract_is_explicit(self) -> None:
        trapped = self.observed["artifact_payload"]["restart_admission"][
            "trapped_assessment"
        ]
        self.assertEqual(trapped["primary_point_counts"], [2049, 4097, 8193])
        self.assertEqual(
            trapped["comparator_point_counts"], [4097, 8193, 16385]
        )
        self.assertEqual(trapped["common_physical_point_count"], 2049)
        self.assertIn("cross_method_interval_overlap", trapped)
        frozen = self.observed["artifact_payload"]["frozen_run_plan"]
        self.assertEqual(
            frozen["numerical_runtime_contract"], EXPECTED_NUMERICAL_RUNTIME
        )
        environment = self.observed["artifact_payload"][
            "runtime_environment_observed"
        ]
        self.assertEqual(environment["python_version"], "3.14.3")
        self.assertEqual(environment["numpy_version"], "2.5.1")
        self.assertEqual(environment["blas"]["name"], "accelerate")
        lineage = self.observed["artifact_payload"][
            "runtime_lineage_identification"
        ]
        self.assertTrue(lineage["RSP2_PREF14_exact_raw_reproduction_passed"])

    def test_claim_and_ladder_mutations_fail_closed(self) -> None:
        config = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            promoted = root / "promoted.toml"
            promoted.write_text(
                config.replace(
                    "FGCQR_holdout_execution_authorized = false",
                    "FGCQR_holdout_execution_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "frozen scope"):
                load_config(promoted)
            plan_source = (
                REPOSITORY / "configs/fgc/fgc-1-pro13-run1.toml"
            ).read_text(encoding="utf-8")
            wrong = root / "wrong-plan.toml"
            wrong.write_text(
                plan_source.replace(
                    "primary_point_counts = [2049, 4097, 8193]",
                    "primary_point_counts = [4097, 8193, 16385]",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "PROTO13 freeze"):
                validate_run_plan(wrong)
            wrong_runtime = root / "wrong-runtime.toml"
            wrong_runtime.write_text(
                plan_source.replace(
                    'python_version = "3.14.3"',
                    'python_version = "3.14.6"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "PROTO13 freeze"):
                validate_run_plan(wrong_runtime)

    def test_noncanonical_and_promoted_results_fail(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(
                _canonical(self.reproduced).rstrip(), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "sorted canonical"):
                load_canonical_result(noncanonical)
            promoted = deepcopy(self.reproduced)
            promoted["gate_status"]["FGCQR_holdout_execution_authorized"] = True
            target = root / "promoted.json"
            target.write_text(_canonical(promoted), encoding="utf-8")
            observed = load_canonical_result(target)
            self.assertTrue(
                observed["gate_status"]["FGCQR_holdout_execution_authorized"]
            )
            self.assertNotEqual(observed, self.reproduced)


if __name__ == "__main__":
    unittest.main()
