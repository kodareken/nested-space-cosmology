from fractions import Fraction
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    classify_tdg6_channel,
    tdg6_design_preflight,
    tdg6_path_accounting,
    three_halves_order_passes_squared,
)


class TDG6TemporalAdmissionDesignTests(unittest.TestCase):
    def test_exact_squared_three_halves_boundary(self) -> None:
        self.assertTrue(
            three_halves_order_passes_squared(
                outer_lower_squared=Fraction(8),
                finest_upper_squared=Fraction(1),
            )
        )
        self.assertFalse(
            three_halves_order_passes_squared(
                outer_lower_squared=Fraction(8) - Fraction(1, 2**52),
                finest_upper_squared=Fraction(1),
            )
        )

    def test_zero_floor_resolved_and_inconclusive_classes(self) -> None:
        zero = classify_tdg6_channel(
            CertifiedMagnitudeInterval(Fraction(0), Fraction(0)),
            CertifiedMagnitudeInterval(Fraction(0), Fraction(0)),
        )
        self.assertEqual(zero.classification, "exact_zero")
        self.assertTrue(zero.admission_passed)
        self.assertEqual(zero.finest_pair_debit, 0)

        floor = classify_tdg6_channel(
            CertifiedMagnitudeInterval(Fraction(0), Fraction(1, 100)),
            CertifiedMagnitudeInterval(Fraction(0), Fraction(1, 200)),
        )
        self.assertEqual(floor.classification, "enclosure_dominated_debit_only")
        self.assertIsNone(floor.order_threshold_passed)
        self.assertEqual(floor.finest_pair_debit, Fraction(1, 200))

        resolved = classify_tdg6_channel(
            CertifiedMagnitudeInterval(Fraction(8), Fraction(8)),
            CertifiedMagnitudeInterval(Fraction(1), Fraction(1)),
        )
        self.assertEqual(resolved.classification, "resolved_order_pass")
        self.assertTrue(resolved.order_threshold_passed)

        inconclusive = classify_tdg6_channel(
            CertifiedMagnitudeInterval(Fraction(0), Fraction(1, 100)),
            CertifiedMagnitudeInterval(Fraction(1, 1000), Fraction(1, 50)),
        )
        self.assertEqual(inconclusive.classification, "order_inconclusive")
        self.assertTrue(inconclusive.temporal_retry_permitted)

    def test_magnitude_does_not_replace_convergence(self) -> None:
        tiny = classify_tdg6_channel(
            CertifiedMagnitudeInterval(Fraction(1, 10**12), Fraction(1, 10**12)),
            CertifiedMagnitudeInterval(Fraction(9, 10**13), Fraction(9, 10**13)),
        )
        self.assertEqual(tiny.classification, "resolved_order_failure")
        self.assertFalse(tiny.admission_passed)

        large = classify_tdg6_channel(
            CertifiedMagnitudeInterval(Fraction(80), Fraction(80)),
            CertifiedMagnitudeInterval(Fraction(10), Fraction(10)),
        )
        self.assertEqual(large.classification, "resolved_order_pass")
        self.assertEqual(large.finest_pair_debit, 10)

    def test_three_level_path_accounting(self) -> None:
        self.assertEqual(len(TDG6_COMPLETE_STATE_CHANNELS), 18)
        self.assertEqual(
            tdg6_path_accounting("RK4"),
            {
                "coarse_proposals": 1,
                "medium_proposals": 2,
                "fine_proposals": 4,
                "total_shadow_proposals": 7,
                "total_shadow_stage_records": 35,
                "committed_fine_proposals": 4,
                "committed_fine_stage_records": 20,
            },
        )
        self.assertEqual(tdg6_path_accounting("SSPRK3")["total_shadow_stage_records"], 28)
        self.assertEqual(tdg6_path_accounting("SSPRK3")["committed_fine_stage_records"], 16)
        with self.assertRaises(ValueError):
            tdg6_path_accounting("Euler")

    def test_invalid_intervals_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "lower bound exceeds"):
            CertifiedMagnitudeInterval(Fraction(2), Fraction(1))
        with self.assertRaisesRegex(ValueError, "must be nonnegative"):
            CertifiedMagnitudeInterval(Fraction(-1), Fraction(1))
        with self.assertRaisesRegex(TypeError, "must be rational"):
            CertifiedMagnitudeInterval(True, Fraction(1))
        with self.assertRaisesRegex(TypeError, "outer_difference"):
            classify_tdg6_channel(  # type: ignore[arg-type]
                (Fraction(0), Fraction(1)),
                CertifiedMagnitudeInterval(Fraction(0), Fraction(1)),
            )

    def test_preflight_passes_without_absolute_threshold(self) -> None:
        result = tdg6_design_preflight()
        self.assertTrue(result["controls_passed"])
        self.assertTrue(result["no_absolute_state_tolerance"])
        self.assertFalse(result["physical_signal_used_for_normalization"])
        self.assertEqual(result["minimum_observed_order"], "3/2")


if __name__ == "__main__":
    unittest.main()
