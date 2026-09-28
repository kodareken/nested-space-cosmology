from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt4_mon4 import _canonical  # noqa: E402
from scripts.reproduce_fgc_hlt10_mon10 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    EXPECTED_CLAIMS,
    load_canonical_result,
    load_config,
    record,
    validate_run_plan,
)


class FGCHLT10MON10ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.observed = load_canonical_result(DEFAULT_OUTPUT)
        historical = cls.observed["artifact_payload"]["namespace_precondition"]
        cls.reproduced = record(
            DEFAULT_CONFIG,
            historical_namespace_evidence=historical,
        )

    def test_canonical_result_reproduces_exactly(self) -> None:
        self.assertEqual(self.observed, self.reproduced)

    def test_authorization_is_narrow_and_all_evidence_is_fresh(self) -> None:
        self.assertEqual(self.observed["gate_status"], EXPECTED_CLAIMS)
        plan = self.observed["artifact_payload"]["frozen_run_plan"]
        self.assertEqual(plan["expanded_input_count"], 12)
        self.assertEqual(plan["u_p_q_match_count"], 12)
        self.assertEqual(plan["SRC4_source_precheck_passed_count"], 12)
        self.assertEqual(plan["PROTO12_t0_passed_count"], 4)
        self.assertEqual(plan["source_point_batch_size"], 2048)
        self.assertTrue(
            self.observed["gate_status"][
                "PROTO12_fresh_GR0_dynamic_calibration_authorized"
            ]
        )
        self.assertFalse(
            self.observed["gate_status"]["SGBL_execution_authorized"]
        )
        self.assertFalse(
            self.observed["gate_status"]["FGCQR_holdout_execution_authorized"]
        )

    def test_controls_exercise_source_and_pairwise_surfaces(self) -> None:
        controls = self.observed["artifact_payload"]["runtime_controls"]
        self.assertEqual(controls["exact_reference_control_count"], 6)
        self.assertEqual(
            controls["SRC4_SRC3_acceleration_difference_infinity"], 0.0
        )
        self.assertGreater(
            controls["diagnostically_saturated_pair_classification_count"], 0
        )
        self.assertGreater(
            controls["public_t0_raw_adjacent_pair_failure_count"], 0
        )
        self.assertTrue(
            controls["historical_CAL8_evolved_coarse_pair_failure_retained"]
        )
        self.assertFalse(controls["new_numeric_tolerance_floor_or_fit_added"])

    def test_config_and_plan_broadening_fail_closed(self) -> None:
        config = DEFAULT_CONFIG.read_text(encoding="utf-8")
        plan_path = REPOSITORY / "configs/fgc/fgc-1-cal9-run1.toml"
        plan = plan_path.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            promoted = root / "promoted.toml"
            promoted.write_text(
                config.replace(
                    "classical_spherical_diagnostic_authorized = false",
                    "classical_spherical_diagnostic_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims"):
                load_config(promoted)
            backend = root / "backend.toml"
            backend.write_text(
                plan.replace(
                    'source_evaluator = "SRC4_tensor_contracted_reference_covariant_GR0_REF1"',
                    'source_evaluator = "outcome_selected_fallback"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "source or pairwise"):
                validate_run_plan(backend)
            threshold = root / "threshold.toml"
            threshold.write_text(
                plan.replace(
                    'source_residual_infinity_max = "1/1000000000000"',
                    'source_residual_infinity_max = "1/100000000000"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "changes more"):
                validate_run_plan(threshold)

    def test_noncanonical_duplicate_and_promoted_results_fail(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(
                _canonical(self.reproduced).rstrip(), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)
            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a":1,"a":2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
                load_canonical_result(duplicate)
            promoted = self.observed.copy()
            promoted["gate_status"] = dict(promoted["gate_status"])
            promoted["gate_status"]["FGCQR_holdout_execution_authorized"] = True
            promoted_path = root / "promoted.json"
            promoted_path.write_text(_canonical(promoted), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fresh reproduction"):
                from scripts.reproduce_fgc_hlt10_mon10 import verify_canonical

                verify_canonical(promoted_path, DEFAULT_CONFIG)


if __name__ == "__main__":
    unittest.main()
