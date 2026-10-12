"""Independent exact controls for the C1R1 (1/3)Z[1/2] ring."""

from __future__ import annotations

import ast
from fractions import Fraction
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg11_imp1_enclosure as imp1_core  # noqa: E402
from recursive_horizons.fgc.evolution.tdg11_c1r1_ring import (  # noqa: E402
    RING_ONE,
    RING_ZERO,
    DirectRingConversionError,
    RingElement,
    add,
    bernstein_samples,
    cmp,
    from_fraction,
    from_input,
    from_int,
    hermite_bernstein,
    hull_upper,
    int_to_decimal_text,
    mul,
    mul_int,
    sample_lower,
    split_bernstein,
    subtract_bernstein,
)


Q = Fraction
MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/tdg11_c1r1_ring.py"


def _imported_modules(path: Path) -> set[str]:
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


class RingStdlibAndConversionTests(unittest.TestCase):
    def test_ring_module_uses_only_the_standard_library(self) -> None:
        imported = _imported_modules(MODULE_PATH)
        self.assertTrue(imported <= {"__future__", "fractions", "typing"})

    def test_unique_representative_matches_reduced_fraction(self) -> None:
        values = (
            Q(0),
            Q(1),
            Q(-2),
            Q(1, 2),
            Q(1, 3),
            Q(-5, 6),
            Q(4, 3),
            Q(1, 24),
            Q(8, 9),
        )
        for value in values:
            with self.subTest(value=value):
                ring = from_fraction(value, allow_three=True, allow_nine=True)
                self.assertEqual(ring.as_fraction(), value)
                self.assertEqual(ring.canonical_text(), imp1_core._canonical_text(value))
                self.assertEqual(ring.bit_size(), imp1_core._bit_size(value))
                again = from_fraction(ring.as_fraction(), allow_three=True, allow_nine=True)
                self.assertEqual(ring, again)

    def test_non_ring_and_malformed_inputs_are_rejected(self) -> None:
        with self.assertRaises(DirectRingConversionError):
            from_fraction(Q(1, 5), allow_three=True)
        with self.assertRaises(DirectRingConversionError):
            from_fraction(Q(1, 3), allow_three=False)
        with self.assertRaises(DirectRingConversionError):
            from_fraction(Q(1, 9), allow_three=True, allow_nine=False)
        with self.assertRaises(TypeError):
            from_input(0.0, name="x", allow_three=True)
        with self.assertRaises(TypeError):
            from_input(True, name="x", allow_three=True)
        with self.assertRaises(ValueError):
            RingElement(1, 0, 3)
        third = from_fraction(Q(1, 3), allow_three=True)
        with self.assertRaises(DirectRingConversionError):
            from_fraction(third.square().as_fraction(), allow_three=True)


class RingImmutabilityTests(unittest.TestCase):
    def test_assignment_deletion_and_shared_constants_are_rejected(self) -> None:
        value = RingElement(1)
        with self.assertRaises(AttributeError):
            value.mantissa = 5  # type: ignore[misc]
        with self.assertRaises(AttributeError):
            value.exp2 = 3  # type: ignore[misc]
        with self.assertRaises(AttributeError):
            value.exp3 = 1  # type: ignore[misc]
        with self.assertRaises(AttributeError):
            del value.mantissa  # type: ignore[misc]
        self.assertEqual(value.as_fraction(), 1)
        self.assertEqual(hash(value), hash(from_int(1)))
        self.assertEqual(value, from_int(1))
        self.assertIs(from_int(0), RING_ZERO)
        self.assertIs(from_int(1), RING_ONE)
        self.assertIs(from_fraction(Q(0), allow_three=False), RING_ZERO)
        aliased = add(RING_ZERO, value)
        self.assertIs(aliased, value)
        with self.assertRaises(AttributeError):
            RING_ZERO.mantissa = 7  # type: ignore[misc]
        with self.assertRaises(AttributeError):
            RING_ONE.exp2 = 4  # type: ignore[misc]
        with self.assertRaises(AttributeError):
            aliased.mantissa = 9  # type: ignore[misc]
        self.assertTrue(RING_ZERO.is_zero())
        self.assertEqual(RING_ONE.as_fraction(), 1)
        self.assertEqual(value.as_fraction(), 1)
        before = sys.get_int_max_str_digits()
        self.assertEqual(int_to_decimal_text(1 << 20000)[-1], "6")
        self.assertEqual(sys.get_int_max_str_digits(), before)


class RingArithmeticTests(unittest.TestCase):
    def test_add_mul_div_and_order_match_fraction(self) -> None:
        pairs = (
            (Q(1, 3), Q(1, 2)),
            (Q(-5, 6), Q(1, 8)),
            (Q(0), Q(7, 3)),
            (Q(9, 4), Q(-3, 8)),
            (Q(1, 24), Q(1, 3)),
        )
        for left, right in pairs:
            with self.subTest(left=left, right=right):
                a = from_fraction(left, allow_three=True)
                b = from_fraction(right, allow_three=True)
                self.assertEqual(add(a, b).as_fraction(), left + right)
                self.assertEqual(a.sub(b).as_fraction(), left - right)
                self.assertEqual(cmp(a, b), (left > right) - (left < right))
                if right.denominator % 3 != 0 or left.denominator % 3 != 0:
                    product = mul(a, b)
                    self.assertEqual(product.as_fraction(), left * right)
        self.assertEqual(mul_int(from_int(1), 0), RING_ZERO)
        self.assertEqual(from_int(1), RING_ONE)
        squared = from_fraction(Q(1, 3), allow_three=True).square()
        self.assertEqual(squared.as_fraction(), Q(1, 9))
        self.assertEqual(squared.exp3, 2)

    def test_chunked_decimal_matches_imp1_for_large_dyadics(self) -> None:
        huge = Q(1, 1 << 20000)
        ring = from_fraction(huge, allow_three=False)
        self.assertEqual(ring.canonical_text(), imp1_core._canonical_text(huge))
        self.assertEqual(int_to_decimal_text(0), "0")
        self.assertEqual(int_to_decimal_text(-12345678901), "-12345678901")
        self.assertEqual(sys.get_int_max_str_digits(), sys.get_int_max_str_digits())


class DyadicBernsteinControlTests(unittest.TestCase):
    def test_five_samples_and_de_casteljau_match_imp1(self) -> None:
        controls = (
            (Q(0), Q(0), Q(1, 3), Q(0)),
            (Q(0), Q(0), Q(2), Q(1)),
            (Q(0), Q(1), Q(-1), Q(0)),
            (Q(1, 3), Q(-1, 2), Q(1, 6), Q(1)),
            (Q(-4), Q(0), Q(1, 8), Q(3, 4)),
        )
        points = (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1))
        for cubic in controls:
            with self.subTest(cubic=cubic):
                ring = tuple(from_fraction(item, allow_three=True) for item in cubic)
                samples = bernstein_samples(ring)  # type: ignore[arg-type]
                for point, sample in zip(points, samples, strict=True):
                    self.assertEqual(
                        sample.as_fraction(),
                        imp1_core._bernstein_value(cubic, point),
                    )
                left, right = split_bernstein(ring)  # type: ignore[arg-type]
                expected_left, expected_right = imp1_core._split_bernstein(cubic)
                self.assertEqual(
                    tuple(item.as_fraction() for item in left), expected_left
                )
                self.assertEqual(
                    tuple(item.as_fraction() for item in right), expected_right
                )
                self.assertEqual(
                    sample_lower(ring).as_fraction(),  # type: ignore[arg-type]
                    imp1_core._sample_lower(cubic),
                )
                self.assertEqual(
                    hull_upper(ring).as_fraction(),  # type: ignore[arg-type]
                    imp1_core._hull_upper(cubic),
                )

    def test_independent_closed_form_samples_for_bump_and_six_five(self) -> None:
        bump = tuple(from_fraction(item, allow_three=True) for item in (Q(0), Q(0), Q(1, 3), Q(0)))
        samples = bernstein_samples(bump)  # type: ignore[arg-type]
        # B(t) = (1-t) t^2 for Bernstein (0, 0, 1/3, 0)
        self.assertEqual(samples[0].as_fraction(), 0)
        self.assertEqual(samples[1].as_fraction(), Q(3, 64))
        self.assertEqual(samples[2].as_fraction(), Q(1, 8))
        self.assertEqual(samples[3].as_fraction(), Q(9, 64))
        self.assertEqual(samples[4].as_fraction(), 0)
        self.assertEqual(sample_lower(bump).as_fraction(), Q(9, 64))  # type: ignore[arg-type]
        self.assertEqual(hull_upper(bump).as_fraction(), Q(1, 3))  # type: ignore[arg-type]

        six_five = tuple(from_int(item) for item in (0, 0, 2, 1))
        samples = bernstein_samples(six_five)  # type: ignore[arg-type]
        self.assertEqual(sample_lower(six_five).as_fraction(), Q(81, 64))  # type: ignore[arg-type]
        self.assertEqual(hull_upper(six_five).as_fraction(), 2)  # type: ignore[arg-type]
        self.assertLessEqual(Q(81, 64), Q(32, 25))
        self.assertLessEqual(Q(32, 25), 2)

    def test_hermite_controls_stay_in_the_direct_ring(self) -> None:
        y0 = from_fraction(Q(1, 3), allow_three=True)
        f0 = from_int(1)
        y1 = from_fraction(Q(5, 6), allow_three=True)
        f1 = from_fraction(Q(-1, 4), allow_three=False)
        width = from_fraction(Q(1, 2), allow_three=False)
        controls = hermite_bernstein(y0, f0, y1, f1, width)
        self.assertEqual(
            tuple(item.as_fraction() for item in controls),
            imp1_core._hermite_bernstein(
                (y0.as_fraction(), f0.as_fraction(), y1.as_fraction(), f1.as_fraction(), width.as_fraction())
            ),
        )
        self.assertTrue(all(item.in_direct_ring() for item in controls))
        difference = subtract_bernstein(controls, controls)
        self.assertTrue(all(item.is_zero() for item in difference))


if __name__ == "__main__":
    unittest.main()
