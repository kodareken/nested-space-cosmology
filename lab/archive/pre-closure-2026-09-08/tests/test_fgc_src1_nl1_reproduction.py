from __future__ import annotations

from pathlib import Path
import json
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_run1_sym1 import record as run1_record  # noqa: E402
from scripts.reproduce_fgc_src1_nl1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCSRC1NL1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_quantitative_controls_and_nonclaims(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "nonlinear_REF1_source_and_branch_solver_verified"
            ]
        )
        quantitative = self.payload["artifact_payload"]["quantitative_evidence"]
        self.assertTrue(
            all(item["passed"] for item in quantitative["rational_exact_float_controls"])
        )
        self.assertTrue(quantitative["nonrational_Jacobian_control"]["passed"])
        self.assertTrue(
            quantitative["continuation_path"]["all_points_converged_on_declared_branch"]
        )
        self.assertTrue(
            quantitative["typed_failure_injection"][
                "all_declared_failures_typed_before_crossing"
            ]
        )
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_RUN1_retains_SRC1_among_seven_of_eight_without_authorization(self) -> None:
        run1 = run1_record()
        audit = run1["scoped_run_authorization_audit"]["authorization_audits"][
            "classical_spherical_diagnostic"
        ]
        self.assertEqual(audit["passed_predicate_count"], 7)
        self.assertEqual(audit["required_predicate_count"], 8)
        self.assertTrue(
            next(
                predicate["passed"]
                for predicate in audit["required_predicates"]
                if predicate["id"]
                == "nonlinear_unredefined_REF1_source_and_acceleration_solver_verified"
            )
        )
        self.assertTrue(
            next(
                predicate["passed"]
                for predicate in audit["required_predicates"]
                if predicate["id"]
                == "physical_gauge_and_reduction_constraint_system_closed"
            )
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
            relaxed = root / "relaxed.toml"
            relaxed.write_text(
                source.replace(
                    'residual_infinity_tolerance = "1/1000000000000"',
                    'residual_infinity_tolerance = "1/1000000"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "differ from PROTO1"):
                load_config(relaxed)

            second_equations = root / "second.toml"
            second_equations.write_text(
                source.replace(
                    "no_second_field_equation_implementation = true",
                    "no_second_field_equation_implementation = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "backend contract"):
                load_config(second_equations)

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
