from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.third_jet import (
    JET3_COMPONENT_ORDER,
    THIRD_DERIVATIVE_ORDER,
    Jet3,
)


class ThirdJetTests(unittest.TestCase):
    def test_exact_truncation_and_directional_two_jets(self) -> None:
        jet = Jet3(1, 2, 3, 5, 7, 11, 13, 17, 19, 23)
        self.assertEqual(jet.as_jet2().value, Q(1))
        self.assertEqual(jet.as_jet2().drr, Q(11))
        component_names = ("value", "dt", "dr", "dtt", "dtr", "drr")
        self.assertEqual(
            tuple(getattr(jet.partial_t_jet2(), name) for name in component_names),
            (Q(2), Q(5), Q(7), Q(13), Q(17), Q(19)),
        )
        self.assertEqual(
            tuple(getattr(jet.partial_r_jet2(), name) for name in component_names),
            (Q(3), Q(7), Q(11), Q(17), Q(19), Q(23)),
        )
        # Equality of mixed coordinate partials is encoded once in the jet,
        # rather than inferred from a numerical finite-difference probe.
        self.assertEqual(jet.partial_t_jet2().dr, jet.partial_r_jet2().dt)
        self.assertEqual(jet.partial_t_jet2().drr, jet.partial_r_jet2().dtr)

    def test_mapping_roundtrip_constant_and_orders(self) -> None:
        parsed = Jet3.from_mapping({"value": Q(3, 2), "dtrr": Q(-5, 7)})
        self.assertEqual(parsed.value, Q(3, 2))
        self.assertEqual(parsed.dtrr, Q(-5, 7))
        self.assertEqual(parsed.dttt, 0)
        self.assertEqual(tuple(parsed.as_mapping()), JET3_COMPONENT_ORDER)
        self.assertEqual(THIRD_DERIVATIVE_ORDER, ("dttt", "dttr", "dtrr", "drrr"))
        self.assertEqual(Jet3.constant(4), Jet3(4))

    def test_failure_closed_exactness_mapping_and_direction(self) -> None:
        with self.assertRaisesRegex(TypeError, "exact rational"):
            Jet3(1.0)
        with self.assertRaisesRegex(TypeError, "exact rational"):
            Jet3(True)
        with self.assertRaisesRegex(ValueError, "requires value"):
            Jet3.from_mapping({"dt": 1})
        with self.assertRaisesRegex(ValueError, "unknown keys"):
            Jet3.from_mapping({"value": 0, "drrrr": 1})
        with self.assertRaisesRegex(ValueError, "t or r"):
            Jet3(0).directional_jet2("theta")


if __name__ == "__main__":
    unittest.main()
