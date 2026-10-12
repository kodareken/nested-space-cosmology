"""Focused qualification controls for the prospective TDG7 lattice helper."""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
import math
from pathlib import Path
import struct
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.tdg7_binary64_subdivision_lattice import (  # noqa: E402
    CoordinateLatticeLimitReached,
    cal11_event24_witness,
    plan_forward_proto14_subdivision,
    reconstruct_historical_tdg6_guard,
    validate_tdg7_binary64_subdivision_plan,
)


def _bits(value: float) -> bytes:
    return struct.pack(">d", value)


class TDG7Binary64SubdivisionLatticeTests(unittest.TestCase):
    def assert_exact_plan(self, plan: object) -> None:
        boundaries = plan.boundaries  # type: ignore[attr-defined]
        fine = Fraction.from_float(plan.fine_width)  # type: ignore[attr-defined]
        macro = Fraction.from_float(plan.macro_width)  # type: ignore[attr-defined]
        exact = tuple(Fraction.from_float(item) for item in boundaries)
        self.assertEqual(exact[-1] - exact[0], macro)
        self.assertEqual([exact[i + 1] - exact[i] for i in range(4)], [fine] * 4)
        self.assertEqual(plan.one_full, (boundaries[0], boundaries[4]))  # type: ignore[attr-defined]
        self.assertEqual(plan.two_half, (boundaries[0], boundaries[2], boundaries[4]))  # type: ignore[attr-defined]
        self.assertIs(plan.four_quarter, plan.boundaries)  # type: ignore[attr-defined]
        self.assertEqual(_bits(plan.one_full[-1]), _bits(plan.two_half[-1]))  # type: ignore[attr-defined]
        self.assertEqual(_bits(plan.two_half[-1]), _bits(plan.four_quarter[-1]))  # type: ignore[attr-defined]
        for stage_path in (
            plan.one_full_stage_times,  # type: ignore[attr-defined]
            plan.two_half_stage_times,  # type: ignore[attr-defined]
            plan.four_quarter_stage_times,  # type: ignore[attr-defined]
        ):
            for left, midpoint, right in stage_path:
                self.assertLess(left, midpoint)
                self.assertLess(midpoint, right)
                self.assertEqual(
                    Fraction.from_float(midpoint),
                    Fraction.from_float(left)
                    + (Fraction.from_float(right) - Fraction.from_float(left)) / 2,
                )

    def test_cal11_event24_witness_floors_the_unsafe_tail(self) -> None:
        cap = float.fromhex("0x1.aaa90b0fb5c26p-8")
        start = 23.0 / 16.0
        historical = reconstruct_historical_tdg6_guard(start, cap)
        self.assertEqual(historical.failed_counts, (1, 2, 4))
        self.assertEqual(
            tuple(item.hex() for item in historical.nominal_steps),
            (
                "0x1.aaa90b0fb5c26p-8",
                "0x1.aaa90b0fb5c26p-9",
                "0x1.aaa90b0fb5c26p-10",
            ),
        )
        self.assertEqual(
            tuple(item[0].hex() for item in historical.rounded_differences),
            (
                "0x1.aaa90b0fb5c00p-8",
                "0x1.aaa90b0fb5c00p-9",
                "0x1.aaa90b0fb5c00p-10",
            ),
        )

        plan = cal11_event24_witness()
        self.assertEqual(plan.quantum.hex(), "0x1.0000000000000p-52")
        self.assertEqual(plan.macro_width.hex(), "0x1.aaa90b0fb5800p-8")
        self.assertEqual(plan.fine_width.hex(), "0x1.aaa90b0fb5800p-10")
        self.assertLess(plan.macro_width, cap)
        self.assertEqual(plan.conservative_reduction.hex(), "0x1.0980000000000p-50")
        self.assertEqual(
            Fraction.from_float(plan.conservative_reduction),
            Fraction(531, 576460752303423488),
        )
        self.assertEqual(
            tuple(item.hex() for item in plan.boundaries),
            (
                "0x1.7000000000000p+0",
                "0x1.706aaa42c3ed6p+0",
                "0x1.70d5548587dacp+0",
                "0x1.713ffec84bc82p+0",
                "0x1.71aaa90b0fb58p+0",
            ),
        )
        self.assertEqual(
            tuple(stage[1].hex() for stage in plan.four_quarter_stage_times),
            (
                "0x1.7035552161f6bp+0",
                "0x1.709fff6425e41p+0",
                "0x1.710aa9a6e9d17p+0",
                "0x1.717553e9adbedp+0",
            ),
        )
        self.assertEqual(
            Fraction.from_float(plan.fine_width) / Fraction.from_float(plan.quantum),
            7329968570070,
        )
        self.assertNotEqual(plan.four_quarter, historical.rounded_boundaries[2])
        self.assert_exact_plan(plan)

    def test_target_limited_step_lands_bitwise_on_event_target(self) -> None:
        target = 24.0 / 16.0
        plan = plan_forward_proto14_subdivision(23.0 / 16.0, target, 1.0)
        self.assertTrue(plan.target_limited)
        self.assertEqual(_bits(plan.endpoint), _bits(target))
        self.assert_exact_plan(plan)

    def test_multistep_schedule_preserves_the_event_lattice(self) -> None:
        current, target, cap = 23.0 / 16.0, 24.0 / 16.0, 0.0065
        count = 0
        while _bits(current) != _bits(target):
            plan = plan_forward_proto14_subdivision(current, target, cap)
            self.assertLessEqual(plan.macro_width, cap)
            self.assert_exact_plan(plan)
            current = plan.endpoint
            count += 1
            self.assertLess(count, 100)
        self.assertGreater(count, 1)

    def test_upward_binade_crossing_uses_the_larger_endpoint_ulp(self) -> None:
        start = 2.0 - 2.0**-49
        target = 2.0 + 2.0**-49
        plan = plan_forward_proto14_subdivision(start, target, 1.0)
        self.assertEqual(plan.quantum, math.ulp(2.0))
        self.assert_exact_plan(plan)
        self.assertEqual(_bits(plan.endpoint), _bits(target))

    def test_off_lattice_phase_is_typed_not_repaired_by_a_smaller_width(self) -> None:
        # The last float below 2 has a phase incompatible with ulp(2).  Its
        # difference from 2 is real but cannot support a four-way exact plan.
        with self.assertRaisesRegex(CoordinateLatticeLimitReached, "endpoint_not_q_aligned"):
            plan_forward_proto14_subdivision(math.nextafter(2.0, -math.inf), 2.0, 1.0)

    def test_too_small_cap_is_a_typed_coordinate_limit(self) -> None:
        start, target = 23.0 / 16.0, 24.0 / 16.0
        q = max(math.ulp(start), math.ulp(target))
        with self.assertRaisesRegex(CoordinateLatticeLimitReached, "no_positive_aligned_width"):
            plan_forward_proto14_subdivision(start, target, 2.0 * q)

    def test_target_with_non_stage_lattice_remainder_is_typed(self) -> None:
        # Both endpoints are Q-aligned but the requested event interval has
        # only one Q tick, so it cannot be completed by a stage-safe macro plan.
        start = 1.5
        q = math.ulp(start)
        target = start + q
        with self.assertRaisesRegex(CoordinateLatticeLimitReached, "target_not_stage_lattice_aligned"):
            plan_forward_proto14_subdivision(start, target, 1.0)

    def test_four_q_endpoint_only_lattice_can_alias_a_fine_midpoint(self) -> None:
        # This is the missing condition in the former 4Q proposal.  The
        # endpoints of a 4Q macro are exact, but its quarter step is Q and its
        # half-stage is Q/2, which binary64 cannot represent at this binade.
        start = 2.0
        q = math.ulp(start)
        target = start + 8.0 * q
        old_macro = 4.0 * q
        old_fine = old_macro / 4.0
        old_midpoint_exact = (
            Fraction.from_float(start) + Fraction.from_float(old_fine) / 2
        )
        old_midpoint_binary64 = float(old_midpoint_exact)
        self.assertNotEqual(Fraction.from_float(old_midpoint_binary64), old_midpoint_exact)
        self.assertEqual(_bits(old_midpoint_binary64), _bits(start))
        # The 8Q scheme refuses to round this cap up or silently use the 4Q
        # endpoint-only lattice.
        with self.assertRaisesRegex(CoordinateLatticeLimitReached, "no_positive_aligned_width"):
            plan_forward_proto14_subdivision(start, target, old_macro)

    def test_nonfinite_and_negative_caps_are_typed(self) -> None:
        for cap in (math.nan, math.inf, -1.0):
            with self.subTest(cap=cap), self.assertRaisesRegex(
                CoordinateLatticeLimitReached,
                "nonfinite_coordinate|nonpositive_requested_cap",
            ):
                plan_forward_proto14_subdivision(23.0 / 16.0, 24.0 / 16.0, cap)

    def test_minimum_width_is_an_exact_typed_limit(self) -> None:
        cap = float.fromhex("0x1.aaa90b0fb5c26p-8")
        with self.assertRaisesRegex(CoordinateLatticeLimitReached, "below_minimum_aligned_width"):
            plan_forward_proto14_subdivision(
                23.0 / 16.0, 24.0 / 16.0, cap, minimum_width=cap
            )
        plan = plan_forward_proto14_subdivision(
            23.0 / 16.0, 24.0 / 16.0, cap, minimum_width=float.fromhex("0x1.0p-8")
        )
        self.assert_exact_plan(plan)

    def test_mutation_sensitive_exact_oracle_rejects_upward_rounding(self) -> None:
        plan = cal11_event24_witness()
        q = Fraction.from_float(plan.quantum)
        cap = Fraction.from_float(plan.requested_cap)
        mutated = Fraction.from_float(plan.macro_width) + 8 * q
        self.assertGreater(mutated, cap)
        # A floor-to-grid implementation passes; replacing it with a ceiling
        # would violate this cap proof on the historical witness.
        self.assertLessEqual(Fraction.from_float(plan.macro_width), cap)

    def test_replace_or_independent_endpoint_mutation_fails_closed(self) -> None:
        plan = cal11_event24_witness()
        broken = list(plan.boundaries)
        broken[-1] = math.nextafter(broken[-1], math.inf)
        with self.assertRaisesRegex(CoordinateLatticeLimitReached, "internal_exactness_failure"):
            replace(plan, boundaries=tuple(broken))
        with self.assertRaisesRegex(CoordinateLatticeLimitReached, "internal_exactness_failure"):
            replace(plan, macro_width=math.nextafter(plan.macro_width, math.inf))
        smaller_macro = plan.macro_width - 4.0 * plan.quantum
        smaller_fine = smaller_macro / 4.0
        smaller_boundaries = tuple(
            plan.current + index * smaller_fine for index in range(5)
        )
        with self.assertRaisesRegex(
            CoordinateLatticeLimitReached, "internal_exactness_failure"
        ):
            replace(
                plan,
                macro_width=smaller_macro,
                fine_width=smaller_fine,
                conservative_reduction=plan.selection_budget - smaller_macro,
                boundaries=smaller_boundaries,
            )
        self.assertIsNone(validate_tdg7_binary64_subdivision_plan(plan))


if __name__ == "__main__":
    unittest.main()
