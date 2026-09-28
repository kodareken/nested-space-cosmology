from __future__ import annotations

from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]

from scripts.reproduce_fgc_con4_phy1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)
from scripts.reproduce_fgc_run1_sym1 import record as run1_record  # noqa: E402


class FGCCON4PHY1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_structural_map_monitor_and_nonclaim_evidence(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "physical_gauge_reduction_constraint_system_closed"
            ]
        )
        evidence = self.payload["artifact_payload"]["quantitative_evidence"]
        self.assertTrue(
            evidence["acceleration_structure"][
                "coordinate_time_accelerations_absent_structurally"
            ]
        )
        self.assertTrue(
            evidence["exact_acceleration_controls"]["all_exact_controls_passed"]
        )
        self.assertEqual(len(evidence["state_controls"]), 3)
        for control in evidence["state_controls"]:
            self.assertTrue(control["physical_constraints"]["exact_recomposition"])
            self.assertTrue(control["normal_gauge_map"]["analytic_matrix_equal"])
            self.assertTrue(
                control["normal_gauge_map"]["determinant_strictly_negative"]
            )
            self.assertEqual(control["gauge_characteristics"]["root_signs"], [-1, 1])
            self.assertFalse(
                control["gauge_characteristics"][
                    "constraint_preserving_boundary_map_derived"
                ]
            )
        self.assertTrue(
            evidence["predecessor_composition"]["conditional_closure"][
                "complete_boundary_free_spherical_constraint_system_closed_conditionally"
            ]
        )
        self.assertTrue(
            evidence["monitor_contract"][
                "smallness_on_one_run_is_not_a_propagation_proof"
            ]
        )
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_RUN1_retains_CON4_at_seven_of_eight_without_authorization(self) -> None:
        run1 = run1_record()
        audit = run1["scoped_run_authorization_audit"]["authorization_audits"][
            "classical_spherical_diagnostic"
        ]
        self.assertEqual(audit["passed_predicate_count"], 7)
        self.assertEqual(audit["required_predicate_count"], 8)
        passed = {
            predicate["id"]
            for predicate in audit["required_predicates"]
            if predicate["passed"]
        }
        self.assertEqual(
            passed,
            {
                "classical_spherical_domain_branch_and_multidirectional_early_kill_envelope",
                "nonlinear_unredefined_REF1_source_and_acceleration_solver_verified",
                "physical_gauge_and_reduction_constraint_system_closed",
                "regular_center_finite_mass_constraint_compatible_initial_data_family",
                "spherical_boundary_or_domain_of_dependence_control",
                "classical_health_and_typed_fail_closed_monitoring",
                "independent_solver_validation_and_affine_measurement_contract",
            },
        )
        for gate in (
            "classical_spherical_diagnostic_authorized",
            "FGCQR_holdout_execution_authorized",
            "retained_EFT_evolution_authorized",
            "physical_transition_claim_authorized",
        ):
            self.assertFalse(run1["gate_status"][gate])

    def test_config_and_result_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    "independent_exact_acceleration_controls_required = true",
                    "independent_exact_acceleration_controls_required = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof-contract"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "physical_transition_claim_authorized = false",
                    "physical_transition_claim_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "fail-closed"):
                load_config(promoted)

            fake_scale = root / "fake-scale.toml"
            fake_scale.write_text(
                source.replace(
                    "absolute_error_scale_supplied_here = false",
                    "absolute_error_scale_supplied_here = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "monitor contract"):
                load_config(fake_scale)

            noncanonical = root / "result.json"
            noncanonical.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)

            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a": 1, "a": 2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unique-key"):
                load_canonical_result(duplicate)


if __name__ == "__main__":
    unittest.main()
