from __future__ import annotations

from dataclasses import replace
from math import pi
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg2_absolute_tail_diagnosis import (  # noqa: E402
    TDG2_CLASS_NAMES,
    TDG2_MEASURE_NAMES,
    TDG2_METHOD_POINT_COUNTS,
    TDG2_SAMPLE_COUNT,
    AbsolutePowerInterval,
    AbsoluteTailObservation,
    absolute_tail_observation,
    classify_absolute_power_ladder,
    diagnose_absolute_tail_ladder,
    synthetic_absolute_tail_controls,
)


def _time() -> np.ndarray:
    return np.linspace(0.0, 63.0 / 16.0, TDG2_SAMPLE_COUNT, dtype=np.float64)


def _carrier() -> np.ndarray:
    index = np.arange(TDG2_SAMPLE_COUNT, dtype=np.float64)
    return np.sin(2.0 * pi * 29.0 * index / TDG2_SAMPLE_COUNT)


class TDG2AbsoluteTailDiagnosisTests(unittest.TestCase):
    def test_synthetic_controls_separate_every_frozen_class(self) -> None:
        observed = synthetic_absolute_tail_controls()
        self.assertTrue(observed["all_control_classes_separated"])
        self.assertFalse(observed["actual_terminal_histories_consumed"])
        self.assertFalse(observed["replacement_temporal_admission_defined"])
        seen = {
            classification
            for control in observed["controls"].values()
            for classification in control["observed_classifications"].values()
        }
        self.assertEqual(seen, set(TDG2_CLASS_NAMES))

    def test_absolute_power_preserves_amplitude_contraction(self) -> None:
        time = _time()
        carrier = _carrier()
        amplitudes = (1.0e-3, 1.25e-4, 1.5625e-5)
        diagnosis = diagnose_absolute_tail_ladder(
            (time, time, time),
            tuple(amplitude * carrier for amplitude in amplitudes),
            method="RK4",
            point_counts=TDG2_METHOD_POINT_COUNTS["RK4"],
        )
        for measure in TDG2_MEASURE_NAMES:
            result = diagnosis["classifications"][measure]
            self.assertEqual(
                result["classification"],
                "contracts_to_zero_at_required_power_order",
            )
            for order in result["zero_power_order_lower_bounds"]:
                self.assertIsNotNone(order)
                self.assertGreaterEqual(order, 3.0)

    def test_decimal_reference_encloses_binary64_fft(self) -> None:
        time = _time()
        values = 1.0e-3 * _carrier() + 2.0e-7 * np.cos(
            2.0 * pi * 7.0 * np.arange(TDG2_SAMPLE_COUNT) / TDG2_SAMPLE_COUNT
        )
        observed = absolute_tail_observation(
            time,
            values,
            round_trip_interpolation_infinity=0.0,
        )
        for measure in TDG2_MEASURE_NAMES:
            interval = getattr(observed, measure)
            self.assertLessEqual(interval.lower, interval.binary64_fft_value)
            self.assertGreaterEqual(interval.upper, interval.binary64_fft_value)
            self.assertGreaterEqual(interval.arithmetic_debit, 0.0)
        self.assertLessEqual(
            abs(
                observed.peak_top_band_coefficient_amplitude
                - observed.binary64_peak_top_band_coefficient_amplitude
            ),
            observed.coefficient_arithmetic_debit,
        )

    def test_interpolation_radius_expands_the_power_interval(self) -> None:
        time = _time()
        values = 1.0e-8 * _carrier()
        exact = absolute_tail_observation(
            time,
            values,
            round_trip_interpolation_infinity=0.0,
        )
        debited = absolute_tail_observation(
            time,
            values,
            round_trip_interpolation_infinity=1.0e-8,
        )
        self.assertEqual(exact.field_tail_power.interpolation_debit, 0.0)
        self.assertGreater(debited.field_tail_power.interpolation_debit, 0.0)
        self.assertLessEqual(debited.field_tail_power.lower, exact.field_tail_power.lower)
        self.assertGreaterEqual(debited.field_tail_power.upper, exact.field_tail_power.upper)
        self.assertFalse(debited.top_band_is_below_declared_amplitude_floor)

    def test_one_bit_signal_mutation_changes_absolute_evidence(self) -> None:
        time = _time()
        values = 1.0e-3 * _carrier()
        original = absolute_tail_observation(
            time,
            values,
            round_trip_interpolation_infinity=0.0,
        )
        bits = values.view(np.uint64).copy()
        bits[17] ^= np.uint64(1 << 44)
        mutated = absolute_tail_observation(
            time,
            bits.view(np.float64),
            round_trip_interpolation_infinity=0.0,
        )
        self.assertNotEqual(
            original.field_tail_power.observed,
            mutated.field_tail_power.observed,
        )

    def test_classifier_rejects_wrong_method_ladder_and_measure(self) -> None:
        time = _time()
        observation = absolute_tail_observation(
            time,
            _carrier(),
            round_trip_interpolation_infinity=0.0,
        )
        with self.assertRaisesRegex(ValueError, "method ladder differs"):
            classify_absolute_power_ladder(
                (observation, observation, observation),
                method="RK4",
                point_counts=(1025, 2049, 4097),
                measure="field_tail_power",
            )
        with self.assertRaisesRegex(ValueError, "unknown TDG2 measure"):
            classify_absolute_power_ladder(
                (observation, observation, observation),
                method="RK4",
                point_counts=TDG2_METHOD_POINT_COUNTS["RK4"],
                measure="total_power",
            )

    def test_nonfinite_nonuniform_and_negative_debit_fail_closed(self) -> None:
        time = _time()
        values = _carrier()
        attacks = []
        nonfinite = values.copy()
        nonfinite[3] = np.nan
        attacks.append((time, nonfinite, 0.0))
        folded = time.copy()
        folded[5] = folded[4]
        attacks.append((folded, values, 0.0))
        attacks.append((time, values, -1.0))
        for proper, signal, debit in attacks:
            with self.assertRaises((TypeError, ValueError)):
                absolute_tail_observation(
                    proper,
                    signal,
                    round_trip_interpolation_infinity=debit,
                )

    def test_claim_promotion_cannot_be_constructed(self) -> None:
        time = _time()
        observation = absolute_tail_observation(
            time,
            _carrier(),
            round_trip_interpolation_infinity=0.0,
        )
        result = classify_absolute_power_ladder(
            (observation, observation, observation),
            method="RK4",
            point_counts=TDG2_METHOD_POINT_COUNTS["RK4"],
            measure="field_tail_power",
        )
        with self.assertRaisesRegex(ValueError, "may not authorize"):
            replace(result, replacement_admission_authorized=True)

    def test_interval_owner_cannot_omit_a_debit(self) -> None:
        with self.assertRaisesRegex(ValueError, "omits"):
            AbsolutePowerInterval(
                observed=1.0,
                binary64_fft_value=1.0,
                arithmetic_debit=0.25,
                interpolation_debit=0.5,
                total_debit=0.1,
                lower=0.0,
                upper=2.0,
            )

    def test_observation_rejects_wrong_floor_decision_type(self) -> None:
        time = _time()
        observed = absolute_tail_observation(
            time,
            _carrier(),
            round_trip_interpolation_infinity=0.0,
        )
        data = {field: getattr(observed, field) for field in observed.__dataclass_fields__}
        data["top_band_is_below_declared_amplitude_floor"] = 1
        with self.assertRaisesRegex(TypeError, "must be bool"):
            AbsoluteTailObservation(**data)


if __name__ == "__main__":
    unittest.main()
