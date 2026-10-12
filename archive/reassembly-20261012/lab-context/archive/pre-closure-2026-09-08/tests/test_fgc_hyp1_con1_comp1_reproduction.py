from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_con1_comp1 import (
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical_json_text,
    _theorem_composition,
    load_canonical_result,
    load_config,
    record,
)


class FGCHYP1CON1COMP1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = record()

    def test_frozen_config_and_in_memory_certificate(self) -> None:
        config = load_config()
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-CON1-COMP1")
        self.assertEqual(config.coordinate_radius, Q(4))
        self.assertEqual(config.phi_value, Q(1, 131072))
        self.assertEqual(config.phi_radial_derivative, Q(1, 131072))
        certificate = self.payload["compatible_constraint_certificate"]
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        composition = self.payload["theorem_composition"]
        self.assertTrue(all(composition["premises"].values()))
        self.assertTrue(composition["all_machine_checked_premises_pass"])
        self.assertTrue(composition["full_mhg_root_implies_normal_nabla_C_zero"])
        self.assertTrue(composition["gauge_extension_zero_at_full_mhg_root"])
        self.assertTrue(composition["unredefined_full_equations_at_local_root"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    @unittest.skipUnless(DEFAULT_OUTPUT.is_file(), "result generation is owned by integration")
    def test_generated_result_is_canonical_and_byte_identical(self) -> None:
        self.assertEqual(self.payload, load_canonical_result())

    def test_config_and_result_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            unknown = root / "unknown.toml"
            unknown.write_text(source + "\nunknown_key = true\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "keys differ"):
                load_config(unknown)

            traversal = root / "traversal.toml"
            traversal.write_text(
                source.replace(
                    'quantified_domain_config = "configs/fgc/fgc-1-hyp1-dom1-qift1.toml"',
                    'quantified_domain_config = "../fgc-1-hyp1-dom1-qift1.toml"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "stay inside|canonical"):
                load_config(traversal)

            noncanonical = root / "noncanonical.toml"
            noncanonical.write_text(source.replace('phi_value = "1/131072"', 'phi_value = "2/262144"'), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical rational"):
                load_config(noncanonical)

            degree = root / "degree.toml"
            degree.write_text(source.replace("maximum_degree = 2", "maximum_degree = 3"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "polynomial certificate"):
                load_config(degree)

            promoted = root / "promoted.toml"
            promoted.write_text(source.replace("evolution_authorized = false", "evolution_authorized = true"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(promoted)

            canonical = _canonical_json_text(self.payload)
            duplicate = root / "duplicate.json"
            duplicate.write_text(canonical.replace("{\n  \"activated_local_parameter_datum\"", "{\n  \"artifact_id\": \"FGC-1-HYP1-CON1-COMP1\",\n  \"activated_local_parameter_datum\"", 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unique-key"):
                load_canonical_result(duplicate)

            noncanonical_result = root / "noncanonical-result.json"
            noncanonical_result.write_text(canonical.rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical_result)

    def test_each_composition_premise_is_required(self) -> None:
        premises = dict(self.payload["theorem_composition"]["premises"])
        for key in premises:
            mutation = dict(premises)
            mutation[key] = False
            composition = _theorem_composition(mutation)
            self.assertFalse(composition["all_machine_checked_premises_pass"])
            if key in {
                "qift_root_is_external_predecessor_not_solved_here",
                "fo1_reduction_rows_zero_exact",
            }:
                self.assertTrue(composition["full_mhg_root_implies_normal_nabla_C_zero"])
                self.assertTrue(composition["gauge_extension_zero_at_full_mhg_root"])
                self.assertTrue(composition["unredefined_full_equations_at_local_root"])
            else:
                self.assertFalse(composition["full_mhg_root_implies_normal_nabla_C_zero"])
                self.assertFalse(composition["gauge_extension_zero_at_full_mhg_root"])
                self.assertFalse(composition["unredefined_full_equations_at_local_root"])


if __name__ == "__main__":
    unittest.main()
