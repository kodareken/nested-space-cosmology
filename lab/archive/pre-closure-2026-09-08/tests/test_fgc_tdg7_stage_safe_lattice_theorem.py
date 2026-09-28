"""Independent theorem controls for the TDG7 stage-safe lattice."""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg7_stage_safe_lattice_theorem as theorem  # noqa: E402


class TDG7StageSafeLatticeTheoremTests(unittest.TestCase):
    def test_module_does_not_import_freeze_design_helper(self) -> None:
        source = (ROOT / "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_lattice_theorem.py").read_text()
        self.assertNotIn("tdg7_binary64_subdivision_lattice", source)

    def test_binary64_bits_and_ulp_are_exact(self) -> None:
        self.assertEqual(theorem.binary64_fraction(1.5), Fraction(3, 2))
        self.assertEqual(theorem.binary64_ulp_fraction(1.5), Fraction(1, 2**52))
        self.assertEqual(theorem.binary64_ulp_fraction(2.0), Fraction(1, 2**51))
        self.assertEqual(theorem.binary64_fraction(-0.5), Fraction(-1, 2))

    def test_cal11_witness_and_endpoint_only_four_q_counterexample(self) -> None:
        witness = theorem.cal11_stage_safe_witness_certificate()
        self.assertEqual(witness["historical_failed_counts"], [1, 2, 4])
        self.assertEqual(witness["macro_quantum_hex"], "0x1.0000000000000p-49")
        self.assertEqual(witness["fine_width_over_Q"], 7329968570070)
        self.assertTrue(
            witness[
                "endpoint_only_four_q_counterexample_derived_from_historical_one_step"
            ]
        )
        self.assertFalse(witness["endpoint_only_four_q_midpoint_exact"])
        self.assertEqual(witness["stage_triplet_count"], 7)

    def test_stage_safe_plan_matches_actual_stage_expression(self) -> None:
        plan = theorem.derive_stage_safe_lattice(23.0 / 16.0, 24.0 / 16.0, 0.0065)
        certificate = theorem.verify_stage_safe_plan(plan)
        self.assertEqual(certificate.stage_triplet_count, 7)
        self.assertTrue(certificate.all_stage_times_exact)
        self.assertTrue(certificate.all_midpoints_strictly_interior)
        for path in (plan.outer_stage_times, plan.medium_stage_times, plan.fine_stage_times):
            for left, midpoint, right in path:
                self.assertEqual(
                    theorem.binary64_fraction(midpoint),
                    theorem.binary64_fraction(left)
                    + (theorem.binary64_fraction(right) - theorem.binary64_fraction(left)) / 2,
                )

    def test_target_limited_path_lands_exactly(self) -> None:
        target = 24.0 / 16.0
        plan = theorem.derive_stage_safe_lattice(23.0 / 16.0, target, 1.0)
        self.assertTrue(plan.target_limited)
        self.assertEqual(plan.boundaries[-1].hex(), target.hex())

    def test_typed_failures_and_no_round_up(self) -> None:
        q = math.ulp(2.0)
        cases = [
            ((2.0, 2.0, 1.0), "target_not_ahead"),
            ((1.0, 2.0, 1.0), "outside_frozen_time_envelope"),
            ((2.0, 2.0 + 8 * q, 0.0), "nonpositive_requested_cap"),
            ((2.0, 2.0 + 8 * q, 7 * q), "no_positive_aligned_width"),
            ((2.0, 2.0 + 4 * q, 1.0), "target_not_stage_lattice_aligned"),
            ((math.nextafter(2.0, -math.inf), 2.0, 1.0), "endpoint_not_q_aligned"),
        ]
        for arguments, reason in cases:
            with self.subTest(reason=reason), self.assertRaisesRegex(
                theorem.IndependentTDG7TheoremError, reason
            ):
                theorem.derive_stage_safe_lattice(*arguments)
        with self.assertRaisesRegex(theorem.IndependentTDG7TheoremError, "below_minimum_aligned_width"):
            theorem.derive_stage_safe_lattice(2.0, 2.0 + 16 * q, 16 * q, minimum_width=32 * q)
        for value in (math.nan, math.inf, -math.inf, 10**400, True, "not-a-number"):
            with self.subTest(nonfinite=value), self.assertRaisesRegex(
                theorem.IndependentTDG7TheoremError, "nonfinite_coordinate"
            ):
                theorem.derive_stage_safe_lattice(value, 2.0, 1.0)

        cap = 17 * q
        floored = theorem.derive_stage_safe_lattice(2.0, 2.0 + 32 * q, cap)
        self.assertEqual(floored.macro_width, 16 * Fraction(q))
        self.assertLessEqual(floored.macro_width, theorem.binary64_fraction(cap))

    def test_mutated_plan_cannot_reintroduce_odd_q_fine_width(self) -> None:
        plan = theorem.derive_stage_safe_lattice(2.0, 2.0 + 32 * math.ulp(2.0), 1.0)
        # Four Q ticks retain positive endpoint spacing but make the resulting
        # fine width an odd multiple of Q, precisely the forbidden 4Q route.
        odd_macro = plan.macro_width - 4 * plan.quantum
        odd_fine = odd_macro / 4
        broken = replace(
            plan,
            macro_width=odd_macro,
            fine_width=odd_fine,
            conservative_reduction=plan.selection_budget - odd_macro,
            boundaries=tuple(plan.current + index * float(odd_fine) for index in range(5)),
        )
        with self.assertRaisesRegex(theorem.IndependentTDG7TheoremError, "internal_exactness_failure"):
            theorem.verify_stage_safe_plan(broken)

    def test_forged_plan_outside_envelope_is_rejected(self) -> None:
        q = theorem.binary64_ulp_fraction(1.0)
        macro = 8 * q
        fine = macro / 4
        forged = theorem.IndependentBinary64LatticePlan(
            current=1.0,
            target=float(Fraction(1) + macro),
            requested_cap=float(macro),
            quantum=q,
            macro_quantum=8 * q,
            selection_budget=macro,
            macro_width=macro,
            fine_width=fine,
            conservative_reduction=Fraction(0),
            boundaries=tuple(float(Fraction(1) + index * fine) for index in range(5)),
            target_limited=True,
        )
        with self.assertRaisesRegex(
            theorem.IndependentTDG7TheoremError, "internal_exactness_failure"
        ):
            theorem.verify_stage_safe_plan(forged)

    def test_ast_contract_accepts_current_engine_and_rejects_midpoint_mutation(self) -> None:
        source = (ROOT / "src/recursive_horizons/fgc/evolution/numerical_engine.py").read_text()
        contract = theorem.parse_immutable_engine_stage_contract(source)
        self.assertTrue(contract.candidate_endpoint_uses_final_time)
        self.assertIn(("rk4_k2", "midpoint"), contract.stage_time_kinds)
        broken = source.replace('evaluate("rk4_k2", start + dt / 2, y2)', 'evaluate("rk4_k2", start + dt, y2)')
        with self.assertRaisesRegex(theorem.IndependentTDG7TheoremError, "engine_stage_abscissa_contract_changed"):
            theorem.parse_immutable_engine_stage_contract(broken)
        decoy = broken.replace(
            "    final_time = start + dt",
            '    if False:\n        evaluate("rk4_k2", start + dt / 2, y2)\n\n    final_time = start + dt',
        )
        with self.assertRaisesRegex(
            theorem.IndependentTDG7TheoremError,
            "engine_stage_(control_flow|abscissa_contract)_changed",
        ):
            theorem.parse_immutable_engine_stage_contract(decoy)
        duplicate = source.replace(
            '        k2 = evaluate("rk4_k2", start + dt / 2, y2)',
            '        k2 = evaluate("rk4_k2", start + dt / 2, y2)\n'
            '        evaluate("rk4_k2", start + dt / 2, y2)',
        )
        with self.assertRaisesRegex(
            theorem.IndependentTDG7TheoremError,
            "engine_stage_abscissa_contract_changed",
        ):
            theorem.parse_immutable_engine_stage_contract(duplicate)
        wrong_candidate = source.replace(
            'evaluate("candidate_endpoint", final_time, candidate)',
            'evaluate("candidate_endpoint", start + dt, candidate)',
        )
        with self.assertRaisesRegex(
            theorem.IndependentTDG7TheoremError,
            "engine_stage_abscissa_contract_changed",
        ):
            theorem.parse_immutable_engine_stage_contract(wrong_candidate)

    def test_historical_tdg6_source_contract_is_parsed_and_mutation_rejected(self) -> None:
        source = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py"
        ).read_text()
        contract = theorem.parse_historical_tdg6_subdivision_contract(source)
        self.assertTrue(contract.step_is_width_over_count)
        self.assertTrue(contract.boundaries_are_start_plus_index_times_step)
        self.assertTrue(contract.expected_final_is_start_plus_width)
        self.assertTrue(contract.macro_endpoint_guard_present)
        self.assertTrue(contract.adjacent_uniformity_guard_present)
        broken_step = source.replace("step = width / count", "step = width / (count + 1)")
        with self.assertRaisesRegex(
            theorem.IndependentTDG7TheoremError,
            "historical_tdg6_subdivision_contract_changed",
        ):
            theorem.parse_historical_tdg6_subdivision_contract(broken_step)
        broken_guard = source.replace(
            "not _same_binary64(right - left, step)",
            "not _same_binary64(left - right, step)",
        )
        with self.assertRaisesRegex(
            theorem.IndependentTDG7TheoremError,
            "historical_tdg6_subdivision_contract_changed",
        ):
            theorem.parse_historical_tdg6_subdivision_contract(broken_guard)
        broken_action = source.replace(
            "TDG6 subdivision is not bitwise uniform",
            "TDG6 nonuniformity ignored",
        )
        with self.assertRaisesRegex(
            theorem.IndependentTDG7TheoremError,
            "historical_tdg6_subdivision_contract_changed",
        ):
            theorem.parse_historical_tdg6_subdivision_contract(broken_action)

    def test_deterministic_event_and_binade_controls(self) -> None:
        certificate = theorem.deterministic_stage_safe_property_certificate()
        self.assertGreater(certificate["successful_plan_count"], 3000)
        self.assertGreater(certificate["stage_triplet_count"], 20000)
        self.assertTrue(certificate["all_exact_and_strictly_interior"])


if __name__ == "__main__":
    unittest.main()
