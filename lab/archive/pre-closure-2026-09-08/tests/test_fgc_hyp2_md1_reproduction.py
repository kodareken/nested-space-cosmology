from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hyp2_md1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCHYP2MD1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_all_covector_margins_exact_crosscheck_and_nonclaims(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "quantitative_all_covector_weak_coupling_health_envelope_passed"
            ]
        )
        payload = self.payload["artifact_payload"]
        for key in (
            "all_nonzero_spatial_covectors_covered",
            "strong_hyperbolicity_margin_positive",
            "independent_chi_sector_included",
            "classical_principal_health_only_not_EFT_validity",
        ):
            self.assertTrue(payload[key])
        quantitative = payload["quantitative_evidence"]
        radial = quantitative["independent_exact_radial_crosscheck"]
        self.assertTrue(radial["independent_exact_spherical_symbol_agreement"])
        self.assertLess(
            radial["maximum_absolute_root_difference"], radial["tolerance"]
        )
        extrema = quantitative["witness_extrema"]
        self.assertLess(
            extrema["maximum_companion_deformation_operator_2"], 1.0 / 4096.0
        )
        self.assertLess(
            extrema["maximum_action_energy_deformation_operator_2"],
            1.0 / 8192.0,
        )
        self.assertGreater(
            extrema["minimum_cluster_resolvent_singular_margin"], 0.0
        )
        self.assertGreater(
            extrema["minimum_physical_energy_coercivity_lower"], 0.0
        )
        self.assertLess(
            extrema["maximum_pure_gauge_identity_coefficient_bound"], 1.0e-10
        )
        self.assertEqual(len(quantitative["witness_records"]), 4)
        self.assertTrue(
            all(
                record["passed"] and all(record["predicates"].values())
                for record in quantitative["witness_records"]
            )
        )
        boundary = quantitative["epistemic_boundary"]
        self.assertTrue(boundary["trajectory_containment_is_deferred_to_HLT1"])
        self.assertFalse(boundary["retained_EFT_remainder_or_frequency_control_supplied"])
        self.assertFalse(boundary["time_evolution_performed"])
        self.assertFalse(boundary["FGCQR_holdout_outcome_inspected"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_config_and_result_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    'companion_deformation_operator_2_maximum = "1/4096"',
                    'companion_deformation_operator_2_maximum = "1/64"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "thresholds differ"):
                load_config(weakened)

            scan = root / "scan.toml"
            scan.write_text(
                source.replace(
                    "directional_scan_cannot_substitute_for_coefficient_norm_bound = true",
                    "directional_scan_cannot_substitute_for_coefficient_norm_bound = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "every proof_contract"):
                load_config(scan)

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


if __name__ == "__main__":
    unittest.main()
