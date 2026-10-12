from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt2_mon2 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCHLT2MON2ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_decision_and_boundaries_are_exact(self) -> None:
        gates = self.payload["gate_status"]
        self.assertTrue(gates["PROTO4_successor_admission_monitor_implemented"])
        self.assertTrue(gates["PROTO4_branch_applicability_executable"])
        self.assertTrue(gates["PROTO4_weighted_spectral_admission_executable"])
        self.assertTrue(gates["PROTO4_common_event_constraint_admission_executable"])
        self.assertTrue(gates["FGCQR_candidate_health_definition_bound"])
        self.assertFalse(gates["SGBL_candidate_health_definition_bound"])
        for name in (
            "fresh_GR0_calibration_authorized",
            "SGBL_comparison_execution_authorized",
            "FGCQR_holdout_execution_authorized",
            "classical_spherical_diagnostic_authorized",
            "retained_EFT_evolution_authorized",
            "physical_transition_claim_authorized",
        ):
            self.assertFalse(gates[name])

        payload = self.payload["artifact_payload"]
        branches = payload["branch_applicability"]
        self.assertTrue(branches["GR0"]["runtime_monitor_complete"])
        self.assertEqual(branches["GR0"]["candidate_branch_stop_ids"], [])
        self.assertTrue(branches["FGCQR_with_exact_owner"]["runtime_monitor_complete"])
        self.assertFalse(branches["SGBL_without_owner"]["runtime_monitor_complete"])
        self.assertTrue(branches["SGBL_rejected_FGCQR_cross_owner"])

        evidence = payload["quantitative_evidence"]
        compact = evidence["declared_compact_initial_profile_control"]
        self.assertTrue(compact["nested_admission"]["admission_passed"])
        self.assertTrue(compact["at_least_one_last_bin_diagnostic_true"])
        self.assertTrue(compact["last_bin_diagnostic_ignored_by_decision"])
        self.assertEqual(
            [record["point_count"] for record in compact["records"]],
            [1025, 2049, 4097],
        )
        self.assertTrue(
            all(record["every_individual_budget_passed"] for record in compact["records"])
        )

        constraints = evidence["constraint_and_error_controls"]
        self.assertTrue(
            constraints["manufactured_second_order_common_event"]["admission_passed"]
        )
        self.assertFalse(constraints["nonconvergent_injection"]["admission_passed"])
        self.assertFalse(
            constraints["common_event_misalignment_injection"]["admission_passed"]
        )
        self.assertFalse(
            constraints["DEF1_missing_affine_alignment_injection"]
            ["common_event_alignment_passed"]
        )
        self.assertTrue(
            constraints["DEF1_affine_alignment_control"]
            ["common_event_alignment_passed"]
        )
        self.assertFalse(
            constraints["Richardson_additive_control"]["subtraction_available"]
        )
        boundary = payload["epistemic_boundary"]
        self.assertTrue(boundary["SGBL_health_definition_missing"])
        self.assertTrue(boundary["constraint_to_Raychaudhuri_stability_map_missing"])
        self.assertFalse(boundary["dynamic_trajectory_used"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_config_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    'maximum_top_band_field_power_fraction = "1/1048576"',
                    'maximum_top_band_field_power_fraction = "1/1024"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "spectral admission differs"):
                load_config(weakened)

            cross_owned = root / "cross-owned.toml"
            cross_owned.write_text(
                source.replace(
                    'SGBL_definition_owner = ""',
                    'SGBL_definition_owner = "FGC-1-HYP2-MD1"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "branch applicability differs"):
                load_config(cross_owned)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "FGCQR_holdout_execution_authorized = false",
                    "FGCQR_holdout_execution_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims must remain fail-closed"):
                load_config(promoted)

            subtracted = root / "subtracted.toml"
            subtracted.write_text(
                source.replace(
                    "components_are_summed_without_cancellation = true",
                    "components_are_summed_without_cancellation = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "error-budget contract differs"):
                load_config(subtracted)

    def test_noncanonical_and_duplicate_results_fail(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)

            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a":1,"a":2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unique-key"):
                load_canonical_result(duplicate)


if __name__ == "__main__":
    unittest.main()
