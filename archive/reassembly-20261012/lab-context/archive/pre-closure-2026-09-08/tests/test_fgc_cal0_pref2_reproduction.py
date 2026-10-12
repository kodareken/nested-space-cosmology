from __future__ import annotations

from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]

from scripts.reproduce_fgc_cal0_pref2 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCCAL0PREF2ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_pre_holdout_obstruction_and_claim_boundary(self) -> None:
        gates = self.payload["gate_status"]
        self.assertTrue(
            gates["PROTO3_pre_holdout_numerical_contract_obstruction_verified"]
        )
        self.assertTrue(gates["PROTO4_premise_revision_required"])
        for gate in (
            "PROTO3_resolved_holdout_manifest_authorized",
            "classical_spherical_diagnostic_authorized",
            "FGCQR_holdout_execution_authorized",
            "retained_EFT_evolution_authorized",
            "physical_transition_claim_authorized",
        ):
            self.assertFalse(gates[gate])
        spectral = self.payload["artifact_payload"]["spectral_input_obstruction"]
        self.assertTrue(spectral["all_frozen_resolutions_fail_at_initial_time"])
        self.assertEqual(len(spectral["records"]), 6)
        self.assertTrue(all(item["alias_band_occupied"] for item in spectral["records"]))
        applicability = self.payload["artifact_payload"][
            "branch_applicability_obstruction"
        ]
        self.assertTrue(applicability["scope_conflation_verified"])
        self.assertFalse(applicability["branch_applicability_map_present"])
        boundary = self.payload["artifact_payload"]["epistemic_boundary"]
        self.assertFalse(boundary["FGCQR_evolution_outcome_inspected"])
        self.assertFalse(boundary["FGCQR_action_or_mechanism_tested"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_config_and_result_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    "declared_profile_must_be_evaluated_by_HLT1_estimator = true",
                    "declared_profile_must_be_evaluated_by_HLT1_estimator = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof-contract"):
                load_config(weakened)

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

            noncanonical = root / "result.json"
            noncanonical.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)

            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a":1,"a":2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
                load_canonical_result(duplicate)


if __name__ == "__main__":
    unittest.main()
