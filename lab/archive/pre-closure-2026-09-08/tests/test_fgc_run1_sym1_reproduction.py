from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_run1_sym1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    _load_future_records,
    _load_unique_json,
    load_canonical_result,
    load_config,
    record,
)


class RUN1SYM1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_exact_initial_authorization_boundary(self) -> None:
        audit = self.payload["scoped_run_authorization_audit"]
        scoped = audit["authorization_audits"]["classical_spherical_diagnostic"]
        self.assertEqual(scoped["passed_predicate_count"], 7)
        self.assertEqual(scoped["required_predicate_count"], 8)
        passed = {
            item["id"]
            for item in scoped["required_predicates"]
            if item["passed"] is True
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
        self.assertEqual(
            scoped["missing_predicate_ids"],
            [
                "frozen_outcome_neutral_protocol_and_resolved_holdout_hash",
            ],
        )
        self.assertFalse(self.payload["gate_status"]["classical_spherical_diagnostic_authorized"])
        self.assertFalse(self.payload["gate_status"]["retained_EFT_evolution_authorized"])
        self.assertFalse(self.payload["gate_status"]["physical_transition_claim_authorized"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_claim_and_future_gate_config_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "retained_EFT_evolution_authorized = false",
                    "retained_EFT_evolution_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "fail-closed"):
                load_config(promoted)

            substituted = root / "substituted.toml"
            substituted.write_text(
                source.replace(
                    'required_gate = "regular_center_formulation_verified"',
                    'required_gate = "nonzero_width_finite_mass_constraint_compatible_family_constructed"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "identity or gate differs"):
                load_config(substituted)

    def test_partial_future_pair_and_duplicate_json_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "future.toml"
            config.write_text("artifact_id = 'test'\n", encoding="utf-8")
            paths = {
                name: {"config": root / f"{name}.toml", "result": root / f"{name}.json"}
                for name in (
                    "protocol_holdout",
                    "run_domain",
                    "multidirectional_health",
                    "nonlinear_source",
                    "physical_constraints",
                    "regular_center",
                    "initial_data",
                    "boundary_control",
                    "health_monitor",
                    "numerical_validation",
                )
            }
            paths["regular_center"] = {
                "config": config,
                "result": root / "missing.json",
            }
            with self.assertRaisesRegex(ValueError, "both exist or both be absent"):
                _load_future_records(paths)

            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a": 1, "a": 2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unique-key JSON"):
                _load_unique_json(duplicate, "duplicate")

    def test_noncanonical_result_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            bad_result = Path(directory) / "bad.json"
            bad_result.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(bad_result)


if __name__ == "__main__":
    unittest.main()
