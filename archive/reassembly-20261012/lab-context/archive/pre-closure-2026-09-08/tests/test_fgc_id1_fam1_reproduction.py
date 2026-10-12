from __future__ import annotations

from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]

from scripts.reproduce_fgc_id1_fam1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)
from scripts.reproduce_fgc_run1_sym1 import (  # noqa: E402
    DEFAULT_OUTPUT as RUN1_OUTPUT,
    load_canonical_result as load_run1_result,
    record as run1_record,
)


class FGCID1FAM1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_family_gate_and_quantitative_margins(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "nonzero_width_finite_mass_constraint_compatible_family_constructed"
            ]
        )
        artifact = self.payload["artifact_payload"]
        for key in (
            "nonzero_width_family",
            "finite_Misner_Sharp_mass",
            "independent_chi_pulse",
            "explicit_nonzero_phi_seed",
            "physical_and_gauge_constraints_solved",
            "exact_outer_vacuum_buffer",
            "continuation_failure_reasons_typed",
        ):
            self.assertTrue(artifact[key])
        quantitative = artifact["quantitative_evidence"]
        self.assertEqual(len(quantitative["exact_full_evaluator_controls"]), 2)
        for control in quantitative["exact_full_evaluator_controls"]:
            self.assertTrue(control["exact_constraint_pair_equality"])
            self.assertTrue(control["exact_metric_defined_gauge_zero"])
            self.assertTrue(control["source_is_unredefined_ACT1_VAR1"])
        self.assertTrue(
            quantitative["exact_zero_seed_GR0_reduction"]["exact_pair_equality"]
        )
        convergence = quantitative["central_common_grid_convergence"]
        self.assertGreater(convergence["RK4"]["observed_common_grid_profile_order"], 3.0)
        self.assertGreater(
            convergence["SSPRK3"]["observed_common_grid_profile_order"], 2.5
        )
        family = quantitative["family_box"]
        self.assertEqual(family["exact_positive_parameter_volume"], "3/8388608")
        self.assertEqual(family["tensor_grid_member_count"], 27)
        self.assertTrue(all(family["required_inequalities"].values()))
        self.assertTrue(
            family[
                "smooth_local_open_family_follows_from_regular_ODE_continuous_dependence"
            ]
        )
        self.assertFalse(family["whole_declared_box_interval_enclosure_proven"])
        self.assertFalse(family["PROTO3_width_nine_eighths_case_qualified_here"])
        boundary = quantitative["epistemic_boundary"]
        self.assertTrue(boundary["initial_hypersurfaces_constructed"])
        for key in (
            "time_evolution_performed",
            "collapse_or_trapped_surface_derived",
            "Raychaudhuri_defocusing_tested",
            "FGCQR_holdout_outcomes_inspected",
            "retained_EFT_validity_established",
        ):
            self.assertFalse(boundary[key])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_RUN1_remains_closed_at_seven_of_eight(self) -> None:
        run1 = load_run1_result(RUN1_OUTPUT)
        self.assertEqual(run1, run1_record())
        audit = run1["scoped_run_authorization_audit"]["authorization_audits"][
            "classical_spherical_diagnostic"
        ]
        self.assertEqual(audit["passed_predicate_count"], 7)
        self.assertEqual(audit["required_predicate_count"], 8)
        initial = next(
            predicate
            for predicate in audit["required_predicates"]
            if predicate["id"]
            == "regular_center_finite_mass_constraint_compatible_initial_data_family"
        )
        self.assertTrue(initial["passed"])
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
                    "complete_unredefined_constraints_specialized = true",
                    "complete_unredefined_constraints_specialized = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof-contract"):
                load_config(weakened)

            widened = root / "widened.toml"
            widened.write_text(
                source.replace(
                    'chi_half_width_maximum = "2"',
                    'chi_half_width_maximum = "17/8"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "node endpoints"):
                load_config(widened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "physical_transition_claim_authorized = false",
                    "physical_transition_claim_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims must remain fail-closed"):
                load_config(promoted)

            incomplete = root / "incomplete.toml"
            incomplete.write_text(
                source.replace('chi_pi_r = "1/67"\n', "", 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "keys differ"):
                load_config(incomplete)

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
