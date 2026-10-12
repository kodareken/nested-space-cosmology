from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt4_mon4 import _canonical  # noqa: E402
from scripts.reproduce_fgc_hlt9_mon9 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    EXPECTED_CLAIMS,
    load_canonical_result,
    load_config,
    record,
    validate_run_plan,
)


class FGCHLT9MON9ReproductionTests(unittest.TestCase):
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
        payload = self.observed["artifact_payload"]
        plan = payload["frozen_run_plan"]
        self.assertEqual(plan["expanded_input_count"], 12)
        self.assertEqual(plan["u_p_match_count"], 12)
        self.assertEqual(plan["reference_mapped_q_count"], 12)
        self.assertEqual(plan["source_precheck_passed_count"], 12)
        self.assertEqual(plan["PROTO11_t0_passed_count"], 4)
        self.assertEqual(len(payload["initial_common_events"]), 4)
        self.assertTrue(
            all(
                item["PROTO11_admission_passed"]
                for item in payload["initial_common_events"]
            )
        )
        self.assertTrue(
            all(
                item["complete_raw_source_gate_passed"]
                for item in payload["frozen_run_inputs"]
            )
        )
        self.assertFalse(
            self.observed["gate_status"]["classical_spherical_diagnostic_authorized"]
        )
        self.assertFalse(
            self.observed["gate_status"]["retained_EFT_evolution_authorized"]
        )

    def test_reference_map_controls_attack_special_cases(self) -> None:
        controls = self.observed["artifact_payload"]["reference_map_controls"]
        self.assertEqual(controls["exact_reference_control_count"], 6)
        self.assertTrue(controls["all_exact_reference_controls_bitwise_preserved"])
        self.assertFalse(
            controls["one_bit_perturbation_used_exact_equilibrium_branch"]
        )
        self.assertGreater(controls["generic_reduction_defect_infinity"], 0.0)
        self.assertGreater(controls["generic_source_q_r_defect_infinity"], 0.0)
        self.assertEqual(
            controls["common_event_reduction_map_difference_infinity"], 0.0
        )
        self.assertFalse(controls["interior_q_reprojected_by_any_control"])
        self.assertFalse(controls["new_numeric_tolerance_added"])

    def test_config_and_plan_broadening_fail_closed(self) -> None:
        config = DEFAULT_CONFIG.read_text(encoding="utf-8")
        plan_path = REPOSITORY / "configs/fgc/fgc-1-cal8-run1.toml"
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
            reprojection = root / "reprojection.toml"
            reprojection.write_text(
                plan.replace(
                    "interior_q_reprojection_after_initialization_or_Runge_Kutta_stage = false",
                    "interior_q_reprojection_after_initialization_or_Runge_Kutta_stage = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "reference-map contract"):
                validate_run_plan(reprojection)
            weakened = root / "weakened.toml"
            weakened.write_text(
                plan.replace(
                    'source_residual_infinity_max = "1/1000000000000"',
                    'source_residual_infinity_max = "1/100000000000"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "changes more"):
                validate_run_plan(weakened)

    def test_noncanonical_and_duplicate_results_fail(self) -> None:
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


if __name__ == "__main__":
    unittest.main()
