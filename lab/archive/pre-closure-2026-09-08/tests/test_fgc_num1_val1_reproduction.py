from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_num1_val1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCNUM1VAL1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_all_controls_pass_without_promoting_a_physical_claim(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "independent_solver_validation_and_measurement_contract_passed"
            ]
        )
        payload = self.payload["artifact_payload"]
        for key in (
            "Minkowski_preservation_passed",
            "regular_center_manufactured_solution_passed",
            "GR0_scalar_control_passed",
            "exact_vs_floating_backend_passed",
            "constraint_convergence_passed",
            "boundary_isolation_passed",
            "all_typed_stops_injected",
            "checkpoint_restart_equivalence_passed",
            "direct_and_Raychaudhuri_routes_agree",
            "primary_and_comparator_methods_passed",
        ):
            self.assertTrue(payload[key])
        evidence = payload["quantitative_evidence"]
        self.assertEqual(evidence["typed_stop_count"], 26)
        self.assertEqual(len(evidence["method_controls"]), 2)
        self.assertTrue(
            evidence["runtime_transaction"]["all_transactions_passed"]
        )
        self.assertFalse(
            evidence["epistemic_boundary"]["FGCQR_trajectory_evaluated"]
        )
        self.assertFalse(
            evidence["epistemic_boundary"]["defocusing_margin_measured"]
        )
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_config_and_result_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    'direct_Raychaudhuri_agreement_residual_max = "1/100000000"',
                    'direct_Raychaudhuri_agreement_residual_max = "1/1000"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "thresholds differ"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "FGCQR_trajectory_evaluated = false",
                    "FGCQR_trajectory_evaluated = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "scope differs"):
                load_config(promoted)

            noncanonical = root / "result.json"
            noncanonical.write_text(
                _canonical(self.payload).rstrip(), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)


if __name__ == "__main__":
    unittest.main()
