from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.eft_ledger import (  # noqa: E402
    DEFAULT_CONFIG,
    EFT0_LED1_ARTIFACT_ID,
    load_eft_ledger_config,
    local_eft_ledger_certificate,
)


class FGCEFTLedgerTests(unittest.TestCase):
    def test_frozen_conditional_ledger_is_exact_and_unpromoted(self) -> None:
        config = load_eft_ledger_config()
        certificate = local_eft_ledger_certificate(config)
        self.assertEqual(config.artifact_id, EFT0_LED1_ARTIFACT_ID)
        self.assertEqual(config.project_version, "0.11.0")
        self.assertEqual(config.cutoff, Q(16))
        self.assertEqual(config.scalar_mass, Q(3))
        self.assertEqual(config.quartic_coupling, Q(1, 2))
        self.assertEqual(certificate["component_box"]["phi_absolute_upper"], "3/262144")
        self.assertTrue(all(control["passes"] for control in certificate["local_controls"].values()))
        self.assertFalse(certificate["frequency_control"]["available"])
        self.assertFalse(certificate["retained_eft_validity"])
        self.assertTrue(all(value is False for value in certificate["nonclaims"].values()))
        self.assertEqual(certificate["action_parameters"], {"M_Pl": "2", "mu": "3", "g4": "1/2", "beta": "-1/4", "eta": "1/2"})
        self.assertEqual([item["id"] for item in certificate["operator_basis"]["retained_operators"]], ["einstein_hilbert", "phi2_ricci", "phi_kinetic", "chi_kinetic", "phi_mass", "phi_quartic", "phi2_gauss_bonnet"])
        self.assertEqual([item["id"] for item in certificate["operator_basis"]["representative_omitted_operators"]], ["x_phi_squared", "x_chi_squared", "x_phi_x_chi", "phi6"])
        self.assertFalse(certificate["operator_basis"]["representative_omitted_remainders_locally_evaluated"])

    def test_missing_or_promoted_assumptions_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inferred = root / "inferred.toml"
            inferred.write_text(source.replace("inferred_from_fixture = false", "inferred_from_fixture = true"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must be false"):
                load_eft_ledger_config(inferred)

            floating = root / "floating.toml"
            floating.write_text(source.replace('Lambda = "16"', "Lambda = 16.0"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "float"):
                load_eft_ledger_config(floating)

            wrong_version = root / "wrong-version.toml"
            wrong_version.write_text(source.replace('project_version = "0.11.0"', 'project_version = "0.12.0"'), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "project_version must equal"):
                load_eft_ledger_config(wrong_version)

            bad_dimension = root / "bad-dimension.toml"
            bad_dimension.write_text(source.replace('canonical_dimension = "8"', 'canonical_dimension = "7"', 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "inconsistent canonical dimensions"):
                load_eft_ledger_config(bad_dimension)

            missing_mu = root / "missing-mu.toml"
            missing_mu.write_text(source.replace('mu = "3"\n', ""), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "model keys differ"):
                load_eft_ledger_config(missing_mu)

            missing_g4 = root / "missing-g4.toml"
            missing_g4.write_text(source.replace('g4 = "1/2"\n', ""), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "model keys differ"):
                load_eft_ledger_config(missing_g4)

            missing_g4_reference = root / "missing-g4-reference.toml"
            missing_g4_reference.write_text(source.replace('parameter_references = ["g4"]\n', ""), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "keys differ"):
                load_eft_ledger_config(missing_g4_reference)

            wrong_comparator = root / "wrong-comparator.toml"
            wrong_comparator.write_text(source.replace('comparison = "greater_than"', 'comparison = "less_than"', 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "comparison must be greater_than"):
                load_eft_ledger_config(wrong_comparator)

            dimensionally_ambiguous = root / "dimensionally-ambiguous.toml"
            dimensionally_ambiguous.write_text(source.replace('id = "F_over_Mpl_squared_lower"', 'id = "F_lower"', 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "noncanonical order"):
                load_eft_ledger_config(dimensionally_ambiguous)

            promoted = root / "promoted.toml"
            promoted.write_text(source.replace("retained_eft_validity = false", "retained_eft_validity = true"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must be false"):
                load_eft_ledger_config(promoted)

            crossing = root / "crossing.toml"
            crossing.write_text(source.replace('"1/100"', '"1/100000000000000"'), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "crosses its strict threshold"):
                local_eft_ledger_certificate(load_eft_ledger_config(crossing))


if __name__ == "__main__":
    unittest.main()
