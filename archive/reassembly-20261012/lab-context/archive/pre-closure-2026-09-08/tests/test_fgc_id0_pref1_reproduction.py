from __future__ import annotations

from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]

from scripts.reproduce_fgc_id0_pref1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCID0PREF1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_exact_obstruction_and_epistemic_boundary(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "PROTO1_initial_data_preflight_obstruction_verified"
            ]
        )
        for gate in (
            "PROTO1_initial_data_calibration_authorized",
            "PROTO1_FGCQR_holdout_execution_authorized",
            "classical_spherical_diagnostic_authorized",
            "retained_EFT_evolution_authorized",
            "physical_transition_claim_authorized",
        ):
            self.assertFalse(self.payload["gate_status"][gate])
        self.assertTrue(self.payload["gate_status"]["PROTO2_premise_revision_required"])

        artifact = self.payload["artifact_payload"]
        exact = artifact["exact_constraint_specialization"]
        self.assertEqual(exact["fixture_count"], 2)
        self.assertTrue(exact["all_exact_pair_equalities_passed"])
        analytic = artifact["analytic_entire_box_obstruction"]
        self.assertTrue(analytic["comparisons"]["compactness_strictly_below_separator"])
        self.assertTrue(analytic["comparisons"]["separator_strictly_below_protocol_floor"])
        self.assertTrue(analytic["bound"]["inverse_metric_bootstrap_closed"])
        numerical = artifact["independent_numerical_preflight"]
        self.assertTrue(numerical["original_box"]["all_candidates_below_declared_upper_bound"])
        self.assertTrue(
            numerical["premise_repair_scan"][
                "all_candidates_inside_original_initial_compactness_window"
            ]
        )
        boundary = artifact["epistemic_boundary"]
        self.assertEqual(boundary["protocol_design_route_closed"], "PROTO1_as_written")
        self.assertFalse(boundary["FGCQR_action_or_mechanism_tested"])
        self.assertFalse(boundary["FGCQR_action_or_mechanism_rejected"])
        self.assertFalse(boundary["general_gradient_mechanism_rejected"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_source_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    "entire_original_amplitude_box_bounded_analytically = true",
                    "entire_original_amplitude_box_bounded_analytically = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof-contract"):
                load_config(weakened)

            incomplete = root / "incomplete.toml"
            incomplete.write_text(
                source.replace("inverse_metric_bootstrap_closed = true\n", ""),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof_contract keys differ"):
                load_config(incomplete)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "FGCQR_mechanism_rejected = false",
                    "FGCQR_mechanism_rejected = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims must remain fail-closed"):
                load_config(promoted)

            tuned = root / "tuned.toml"
            tuned.write_text(
                source.replace(
                    'candidate_repair_scan = ["2", "5/2", "3", "7/2", "4", "9/2", "5"]',
                    'candidate_repair_scan = ["2", "5/2", "3", "7/2", "4", "9/2", "11/2"]',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "repair-scan amplitude"):
                record(tuned)

            rewritten_history = root / "rewritten-history.toml"
            rewritten_history.write_text(
                source.replace(
                    "d4f0cc8f58408ef4ee12fb32fe3619916e231795",
                    "0000000000000000000000000000000000000000",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ValueError, "historical pre-PROTO2 checkpoint differs"
            ):
                load_config(rewritten_history)

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
