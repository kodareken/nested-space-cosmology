from __future__ import annotations

from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]

from scripts.reproduce_fgc_ctr1_reg1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)
from scripts.reproduce_fgc_run1_sym1 import record as run1_record  # noqa: E402


class FGCCTR1REG1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_regular_center_payload_and_nonclaim_boundary(self) -> None:
        self.assertTrue(self.payload["gate_status"]["regular_center_formulation_verified"])
        payload = self.payload["artifact_payload"]
        for field in (
            "R_equals_r_times_A_and_v_equals_r_times_V",
            "all_removable_singular_limits_derived_analytically",
            "negative_power_guard_edge_fails_closed",
            "parity_and_elementary_flatness_enforced",
            "finite_curvature_series_controls_passed",
            "first_grid_point_convergence_passed",
        ):
            self.assertTrue(payload[field])
        evidence = payload["quantitative_evidence"]
        self.assertTrue(evidence["Minkowski_control"]["exact_zero_residual_gauge_and_curvature"])
        self.assertTrue(
            evidence["smooth_nontrivial_control"]["first_grid_point_convergence"][
                "all_regular_equation_and_gauge_controls_passed"
            ]
        )
        self.assertTrue(
            all(
                control["rejected_before_certificate"]
                for control in evidence["negative_controls"].values()
            )
        )
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_RUN1_combines_CTR1_and_ID1_but_remains_closed(self) -> None:
        run1 = run1_record()
        audit = run1["scoped_run_authorization_audit"]["authorization_audits"][
            "classical_spherical_diagnostic"
        ]
        self.assertEqual(audit["passed_predicate_count"], 7)
        predicate = next(
            item
            for item in audit["required_predicates"]
            if item["id"] == "regular_center_finite_mass_constraint_compatible_initial_data_family"
        )
        self.assertTrue(predicate["passed"])
        self.assertIn("CTR1-REG1 status=scope_bound_hash_bound_gate_passed", predicate["evidence"])
        self.assertIn("ID1-FAM1 status=scope_bound_hash_bound_gate_passed", predicate["evidence"])
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
                    "first_grid_point_convergence_required = true",
                    "first_grid_point_convergence_required = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof-contract"):
                load_config(weakened)

            incomplete = root / "incomplete.toml"
            incomplete.write_text(
                source.replace(
                    "negative_power_guard_edge_fails_closed = true\n",
                    "",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof_contract keys differ"):
                load_config(incomplete)

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

            widened = root / "widened.toml"
            widened.write_text(
                source.replace(
                    "certified_minimum_power = -4",
                    "certified_minimum_power = -5",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "series contract"):
                load_config(widened)

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
