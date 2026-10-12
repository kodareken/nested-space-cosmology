"""Pure adversarial controls for the LOC2 independent exact localizer."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg9_local_extrema as primary  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_local_extrema_independent_v2 as v2  # noqa: E402


Q = Fraction


class TDG9LocalExtremaV2Tests(unittest.TestCase):
    def test_tiny_value_polynomial_deflates_endpoint_roots(self) -> None:
        delta = Q(1, 1 << 417)
        coefficients = (Q(0), Q(0), -3 * delta, 2 * delta)
        first = primary.localize_absolute_maximum(
            [primary.LocalCubic(coefficients, {"case": "tiny"})],
            maximum_candidates=4,
        )
        second = v2.localize_absolute_maximum_independently_v2(
            [v2.IndependentLocalCubicV2(coefficients, {"case": "tiny"})],
            maximum_candidates=4,
        )
        self.assertEqual(first.candidate_count, 2)
        self.assertEqual(second.candidate_count, 2)
        self.assertEqual(first.classification, second.classification)
        self.assertEqual(
            [(item.location, item.location_ordinal) for item in first.candidates],
            [(item.location, item.location_ordinal) for item in second.candidates],
        )

    def test_constant_offset_does_not_change_derivative_roots(self) -> None:
        roots = v2.isolate_open_unit_roots(Q(0), Q(-2), Q(2))
        shifted = v2.IndependentLocalCubicV2((Q(9), Q(0), Q(-1), Q(2, 3)), {})
        _, a1, a2, a3 = shifted.coefficients
        self.assertEqual(roots, v2.isolate_open_unit_roots(a1, 2 * a2, 3 * a3))
        self.assertEqual(roots, ())

    def test_endpoint_deflation_controls(self) -> None:
        one_third = Q(1, 3)
        cases = (
            ((Q(0), -one_third, Q(1)), (one_third,)),
            ((one_third, -Q(4, 3), Q(1)), (one_third,)),
            ((Q(0), Q(1), Q(-1)), ()),
        )
        for coefficients, expected in cases:
            with self.subTest(coefficients=coefficients):
                roots = v2.isolate_open_unit_roots(*coefficients)
                self.assertEqual(tuple(item.lower for item in roots), expected)
                self.assertTrue(all(item.lower == item.upper for item in roots))
                self.assertTrue(all(item.method == "exact_endpoint_deflation" for item in roots))

    def test_exact_tiny_rational_interior_roots(self) -> None:
        left, right = Q(1, 1 << 417), Q(1, 3)
        roots = v2.isolate_open_unit_roots(left * right, -(left + right), Q(1))
        self.assertEqual(tuple((item.lower, item.upper) for item in roots), ((left, left), (right, right)))
        self.assertTrue(all(item.method == "exact_rational_square" for item in roots))

    def test_two_genuine_irrational_interior_roots(self) -> None:
        roots = v2.isolate_open_unit_roots(Q(1, 10), Q(-1), Q(1))
        self.assertEqual(len(roots), 2)
        self.assertLess(Q(0), roots[0].lower)
        self.assertLess(roots[0].upper, roots[1].lower)
        self.assertLess(roots[1].upper, Q(1))
        self.assertTrue(all(item.method == "adaptive_nonsquare" for item in roots))

    def test_scaling_and_negative_lead_normalize_identically(self) -> None:
        expected = v2.isolate_open_unit_roots(Q(1, 10), Q(-1), Q(1))
        self.assertEqual(expected, v2.isolate_open_unit_roots(Q(7, 10), Q(-7), Q(7)))
        self.assertEqual(expected, v2.isolate_open_unit_roots(Q(-1, 10), Q(1), Q(-1)))

    def test_repeated_linear_zero_and_negative_discriminant_controls(self) -> None:
        repeated = v2.isolate_open_unit_roots(Q(1, 9), Q(-2, 3), Q(1))
        self.assertEqual(tuple(item.lower for item in repeated), (Q(1, 3),))
        self.assertEqual(repeated[0].method, "exact_repeated")
        self.assertEqual(v2.isolate_open_unit_roots(Q(-1, 3), Q(1), Q(0))[0].lower, Q(1, 3))
        self.assertEqual(v2.isolate_open_unit_roots(Q(1), Q(0), Q(0)), ())
        self.assertEqual(v2.isolate_open_unit_roots(Q(1), Q(0), Q(1)), ())

    def test_near_boundary_and_close_nonsquare_roots_resolve(self) -> None:
        near = v2.isolate_open_unit_roots(Q(-3, 1 << 400), Q(0), Q(1))
        self.assertEqual(len(near), 1)
        self.assertGreater(near[0].lower, 0)
        close = v2.isolate_open_unit_roots(Q(1, 4) - Q(2, 10**40), Q(-1), Q(1))
        self.assertEqual(len(close), 2)
        self.assertLess(close[0].upper, close[1].lower)

    def test_proof_schedule_has_exact_final_nondoubling_cap(self) -> None:
        self.assertEqual(v2._refinement_schedule(1000), (192, 384, 768, 1000))  # type: ignore[attr-defined]
        with self.assertRaises(v2.RootIsolationInconclusive) as captured:
            v2._refinement_schedule(v2.GLOBAL_PROOF_BIT_CEILING + 1)  # type: ignore[attr-defined]
        self.assertEqual(captured.exception.reason, "root_proof_cap_exceeded")

    def test_absolute_integer_sqrt_bracket_has_frozen_width(self) -> None:
        lower, upper = v2._sqrt_bounds_integer_absolute(2, 192)  # type: ignore[attr-defined]
        self.assertLessEqual(lower * lower, 2)
        self.assertGreaterEqual(upper * upper, 2)
        self.assertEqual(upper - lower, Q(1, 1 << 192))
        exact = v2._sqrt_bounds_integer_absolute(4, 192)  # type: ignore[attr-defined]
        self.assertEqual(exact, (Q(2), Q(2)))

    def test_sturm_count_and_public_ordinals_cover_two_roots(self) -> None:
        self.assertEqual(v2._sturm_open_unit_root_count(10, -10, 1), 2)  # type: ignore[attr-defined]
        roots = v2.isolate_open_unit_roots(Q(1, 10), Q(-1), Q(1))
        self.assertEqual(len(roots), 2)
        self.assertLess(roots[0].upper, roots[1].lower)
        cubic = v2.IndependentLocalCubicV2((Q(1, 30), Q(1, 10), Q(-1, 2), Q(1, 3)), {})
        evidence = v2.localize_absolute_maximum_independently_v2(
            [cubic], maximum_candidates=4
        )
        interiors = [item for item in evidence.candidates if item.location == "interior_stationary"]
        self.assertEqual(evidence.candidate_count, 4)
        self.assertEqual([item.location_ordinal for item in interiors], [0, 1])
        expected_digest = sha256(b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n0|2\n").hexdigest()
        self.assertEqual(evidence.stationary_count_stream_sha256, expected_digest)

    def test_nonsquare_boundary_straddling_cap_fails_typed(self) -> None:
        with (
            patch.object(v2, "_sqrt_bounds_integer_absolute", return_value=(Q(0), Q(2))),
            self.assertRaises(v2.RootIsolationInconclusive) as captured,
        ):
            v2.isolate_open_unit_roots(Q(-4, 5), Q(0), Q(1))
        self.assertEqual(captured.exception.reason, "root_location_inconclusive")

    def test_candidate_ceiling_and_invalid_controls_fail_closed(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "candidate_ceiling_exhausted"):
            v2.localize_absolute_maximum_independently_v2(
                [v2.IndependentLocalCubicV2((0, 1, 0, 0), {})],
                maximum_candidates=1,
            )
        with self.assertRaises(TypeError):
            v2.isolate_open_unit_roots(True, 0, 1)


if __name__ == "__main__":
    unittest.main()
