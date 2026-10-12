from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "src"))

from scripts.reproduce_fgc_eft0_led1 import (  # noqa: E402
    ACT1_CONFIG,
    ACT1_RESULT,
    COMP1_CONFIG,
    COMP1_RESULT,
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical_json_text,
    _check_result,
    load_canonical_result,
    load_config,
    record,
)


class FGCEFT0LED1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = record()

    def test_frozen_record_is_conditional_exact_and_provenanced(self) -> None:
        config = load_config()
        self.assertEqual(config.ledger.artifact_id, "FGC-1-EFT0-LED1")
        self.assertEqual(config.ledger.project_version, "0.11.0")
        self.assertEqual(self.payload["project_version"], "0.11.0")
        self.assertEqual(self.payload["source_configs"], {
            "eft0_ledger": "configs/fgc/fgc-1-eft0-led1.toml",
            "act1_action_predecessor": "configs/fgc/fgc-1-action-gate.toml",
            "comp1_predecessor": "configs/fgc/fgc-1-hyp1-con1-comp1.toml",
        })
        self.assertEqual(
            self.payload["source_results_sha256"].keys(),
            {
                "results/fgc-1-action-gate.json",
                "results/fgc-1-hyp1-con1-comp1.json",
            },
        )
        certificate = self.payload["conditional_local_ledger_certificate"]
        self.assertTrue(all(item["passes"] for item in certificate["local_controls"].values()))
        self.assertEqual(certificate["local_controls"]["F_over_Mpl_squared_lower"]["comparison"], "greater_than")
        self.assertEqual(certificate["local_controls"]["F_over_Mpl_squared_lower"]["strict_lower_bound"], "99/100")
        self.assertFalse(certificate["retained_eft_validity"])
        self.assertFalse(certificate["frequency_control"]["available"])
        self.assertFalse(certificate["operator_basis"]["representative_omitted_remainders_locally_evaluated"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))
        self.assertTrue(all(value is False for key, value in self.payload["gate_status"].items() if key != "conditional_declared_local_component_controls_passed"))

    def test_generated_result_is_canonical_and_byte_identical(self) -> None:
        self.assertEqual(self.payload, load_canonical_result())
        self.assertEqual(DEFAULT_OUTPUT.read_text(encoding="utf-8"), _canonical_json_text(self.payload))
        _check_result(DEFAULT_OUTPUT, self.payload)

    def test_config_result_and_control_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory(dir=REPOSITORY) as directory:
            root = Path(directory)
            unknown = root / "unknown.toml"
            unknown.write_text(source + "\nunknown_key = true\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "keys differ"):
                load_config(unknown)

            inferred = root / "inferred.toml"
            inferred.write_text(source.replace("inferred_from_fixture = false", "inferred_from_fixture = true"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must be false"):
                load_config(inferred)

            ambiguous = root / "ambiguous.toml"
            ambiguous.write_text(source.replace('id = "F_over_Mpl_squared_lower"', 'id = "F_lower"', 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "noncanonical order"):
                load_config(ambiguous)

            crossing = root / "crossing.toml"
            crossing.write_text(source.replace('strict_upper_bound = "1/100"', 'strict_upper_bound = "1/100000000000000"'), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "crosses its strict threshold"):
                record(crossing)

            duplicate = root / "duplicate.json"
            canonical = _canonical_json_text(self.payload)
            duplicate.write_text(canonical.replace("{\n  \"artifact_id\"", "{\n  \"artifact_id\": \"FGC-1-EFT0-LED1\",\n  \"artifact_id\"", 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unique-key"):
                load_canonical_result(duplicate)

            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(canonical.rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)

            stale = root / "stale.json"
            stale.write_text(canonical.replace('"retained_eft_validity": false', '"retained_eft_validity": true', 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "differs"):
                _check_result(stale, self.payload)

    def test_predecessor_inputs_exist_as_separate_provenance_not_a_validity_promotion(self) -> None:
        self.assertTrue(ACT1_CONFIG.is_file())
        self.assertTrue(ACT1_RESULT.is_file())
        self.assertTrue(COMP1_CONFIG.is_file())
        self.assertTrue(COMP1_RESULT.is_file())
        self.assertFalse(self.payload["gate_status"]["retained_eft_validity_envelope_passed"])
        self.assertFalse(self.payload["gate_status"]["global_or_open_eft_validity_passed"])


if __name__ == "__main__":
    unittest.main()
