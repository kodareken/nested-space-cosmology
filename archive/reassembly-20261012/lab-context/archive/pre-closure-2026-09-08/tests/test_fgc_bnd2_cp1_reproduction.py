from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_bnd2_cp1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCBND2CP1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_domain_control_passes_with_explicit_nonclaims(self) -> None:
        self.assertTrue(
            self.payload["gate_status"][
                "spherical_boundary_or_domain_of_dependence_control_passed"
            ]
        )
        payload = self.payload["artifact_payload"]
        for key in (
            "evolving_physical_and_auxiliary_cones_included",
            "SBP_stencil_reach_included",
            "measured_affine_interval_outside_boundary_domain_of_dependence",
            "outer_boundary_move_check_passed",
            "incoming_constraint_characteristics_controlled",
        ):
            self.assertTrue(payload[key])
        evidence = payload["quantitative_evidence"]
        self.assertEqual(
            evidence["outer_boundary_move_control"][
                "smallest_reference_margin_over_required_buffer"
            ],
            "689/40",
        )
        self.assertTrue(
            evidence["causal_budget_contract"]["runtime_controls"][
                "injected_failure"
            ]["trial_rejected"]
        )
        self.assertTrue(
            evidence["causal_budget_contract"][
                "candidate_endpoint_speed_is_evaluated_explicitly"
            ]
        )
        self.assertEqual(
            evidence["causal_budget_contract"]["runtime_controls"][
                "accepted_control"
            ]["accepted_endpoint_speed_upper"],
            1.3,
        )
        boundary = evidence["epistemic_boundary"]
        self.assertFalse(boundary["complete_nonlinear_constraint_preserving_IBVP_proved"])
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
                    'minimum_remaining_buffer = "16"',
                    'minimum_remaining_buffer = "8"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "causal budget differs"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "nonlinear_constraint_preserving_boundary_map_claimed = false",
                    "nonlinear_constraint_preserving_boundary_map_claimed = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "constraint-control"):
                load_config(promoted)

            extra = root / "extra.toml"
            extra.write_text(
                source.replace(
                    "no_holdout_execution = true\n\n[claims]",
                    "no_holdout_execution = true\nunexpected_contract_flag = true\n\n[claims]",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof_contract keys differ"):
                load_config(extra)

            noncanonical = root / "result.json"
            noncanonical.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)


if __name__ == "__main__":
    unittest.main()
