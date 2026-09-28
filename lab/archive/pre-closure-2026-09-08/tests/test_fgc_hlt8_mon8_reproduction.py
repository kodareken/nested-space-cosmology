from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt4_mon4 import _canonical, _serial  # noqa: E402
from scripts.reproduce_fgc_hlt8_mon8 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_canonical_result,
    load_config,
    record,
    validate_run_plan,
    verify_canonical,
)


class FGCHLT8MON8ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.observed = load_canonical_result(DEFAULT_OUTPUT)
        historical = cls.observed["artifact_payload"]["namespace_precondition"]
        cls.reproduced = _serial(
            record(
                DEFAULT_CONFIG,
                historical_namespace_evidence=historical,
            )
        )

    def test_canonical_result_reproduces_and_authorizes_only_GR0(self) -> None:
        verify_canonical(DEFAULT_OUTPUT, DEFAULT_CONFIG)
        self.assertEqual(self.observed, self.reproduced)
        gates = self.observed["gate_status"]
        self.assertTrue(gates["PROTO10_successor_runtime_implemented"])
        self.assertTrue(gates["PROTO10_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(gates["PROTO10_resolved_holdout_manifest_authorized"])
        self.assertFalse(gates["classical_spherical_diagnostic_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["retained_EFT_evolution_authorized"])
        self.assertFalse(gates["physical_transition_claim_authorized"])
        self.assertTrue(all(value is False for value in self.observed["nonclaims"].values()))

    def test_inputs_raw_failures_and_guarded_admissions_are_all_public(self) -> None:
        payload = self.observed["artifact_payload"]
        frozen = payload["frozen_run_plan"]
        self.assertEqual(frozen["expanded_input_count"], 12)
        self.assertEqual(frozen["projected_state_hash_match_count"], 12)
        self.assertEqual(frozen["raw_PROTO9_t0_passed_count"], 0)
        self.assertEqual(frozen["guarded_PROTO10_t0_passed_count"], 4)
        self.assertEqual(frozen["t0_common_event_evaluated_count"], 4)
        self.assertTrue(
            all(
                item["accepted_state_source_gate_passed"]
                for item in payload["accepted_state_source_prechecks"]
            )
        )
        for event in payload["initial_common_events"]:
            self.assertFalse(event["raw_PROTO9_admission_passed"])
            self.assertTrue(event["guarded_PROTO10_admission_passed"])
            self.assertTrue(event["diagnostic_saturation_used"])
            self.assertFalse(event["diagnostic_saturation_is_physical_resolution"])
            self.assertFalse(event["round_trip_is_a_continuum_error_bound"])

    def test_adversarial_controls_fail_closed(self) -> None:
        controls = self.observed["artifact_payload"]["adversarial_guard_controls"]
        self.assertTrue(controls["all_adversarial_controls_passed"])
        self.assertFalse(controls["coarse_direct_failure_admitted"])
        self.assertFalse(controls["missing_tail_erasure_map_witness_admitted"])
        self.assertFalse(controls["missing_complete_profile_contraction_admitted"])
        self.assertFalse(controls["missing_round_trip_contraction_admitted"])
        self.assertEqual(controls["maximum_nested_tail_ratio"], 0.25)
        self.assertFalse(controls["new_epsilon_floor_added"])

    def test_config_and_run_plan_broadening_fail_closed(self) -> None:
        config = DEFAULT_CONFIG.read_text(encoding="utf-8")
        plan_path = REPOSITORY / "configs/fgc/fgc-1-cal7-run1.toml"
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
            omitted = root / "omitted.toml"
            omitted.write_text(
                plan.replace(
                    "common_event_spatial_guarded_saturation_requires_two_pair_round_trip_contraction = true\n",
                    "",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "spectral interpretation"):
                validate_run_plan(omitted)

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
