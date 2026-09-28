from __future__ import annotations

from dataclasses import replace
from math import pi
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg3_native_tail_diagnosis import (
    TDG3_ESTIMATOR_NAMES,
    TDG3_MEASURE_NAMES,
    TDG3_METHOD_POINT_COUNTS,
    TDG3_SAMPLE_COUNT,
    TDG3_TOP_BINS,
    _binary_piecewise_coefficients,
    diagnose_native_tail_ladder,
    native_tail_observation,
    synthetic_native_tail_controls,
)


def _time(jitter: float = 0.0) -> np.ndarray:
    base = np.linspace(0.0, 63.0 / 16.0, TDG3_SAMPLE_COUNT)
    if jitter == 0.0:
        return base
    phase = np.linspace(0.0, 1.0, TDG3_SAMPLE_COUNT)
    return base + jitter * np.sin(2.0 * pi * phase) * np.sin(pi * phase)


class TDG3NativeTailDiagnosisTests(unittest.TestCase):
    def test_closed_form_piecewise_integral_matches_affine_analytic_result(self) -> None:
        phase = np.linspace(0.0, 1.0, TDG3_SAMPLE_COUNT)
        coefficients = _binary_piecewise_coefficients(phase, phase)
        for frequency_bin, coefficient in zip(
            TDG3_TOP_BINS, coefficients, strict=True
        ):
            expected = 1j / (2.0 * pi * frequency_bin)
            self.assertAlmostEqual(coefficient.real, expected.real, places=13)
            self.assertAlmostEqual(coefficient.imag, expected.imag, places=13)

    def test_both_native_estimators_enclose_their_binary64_audits(self) -> None:
        proper = _time(0.01)
        phase = (proper - proper[0]) / (proper[-1] - proper[0])
        signal = 1.0e-3 * np.sin(2.0 * pi * 29.0 * phase)
        for estimator in TDG3_ESTIMATOR_NAMES:
            observation = native_tail_observation(
                proper, signal, estimator=estimator
            )
            self.assertFalse(observation.common_grid_resampling_performed)
            self.assertFalse(observation.continuum_surrogate_enclosure_claimed)
            self.assertGreaterEqual(observation.normalized_spacing_ratio, 1.0)
            for measure in TDG3_MEASURE_NAMES:
                interval = getattr(observation, measure)
                self.assertLessEqual(interval.lower, interval.binary64_value)
                self.assertGreaterEqual(interval.upper, interval.binary64_value)
                self.assertGreaterEqual(interval.arithmetic_debit, 0.0)

    def test_synthetic_controls_include_the_interpolation_owner_discriminator(self) -> None:
        observed = synthetic_native_tail_controls()
        self.assertTrue(observed["all_control_classes_separated"])
        self.assertFalse(observed["actual_terminal_histories_consumed"])
        self.assertFalse(observed["common_grid_resampling_performed"])
        self.assertFalse(observed["continuum_surrogate_enclosure_claimed"])
        control = observed["common_grid_interpolation_owner_control"]
        self.assertTrue(control["TDG2_identifies_interpolation_owner"])
        self.assertTrue(control["TDG3_native_estimators_agree_nonzero"])

    def test_amplitude_contraction_survives_native_grid_jitter(self) -> None:
        times = tuple(_time(value) for value in (0.02, 0.01, 0.005))
        values = []
        for amplitude, proper in zip(
            (1.0e-3, 1.25e-4, 1.5625e-5), times, strict=True
        ):
            phase = (proper - proper[0]) / (proper[-1] - proper[0])
            values.append(amplitude * np.sin(2.0 * pi * 29.0 * phase))
        diagnosis = diagnose_native_tail_ladder(
            times,
            tuple(values),
            method="RK4",
            point_counts=TDG3_METHOD_POINT_COUNTS["RK4"],
        )
        for measure in TDG3_MEASURE_NAMES:
            self.assertEqual(
                diagnosis["combined_classifications"][measure],
                "native_estimators_agree_zero_contracting",
            )

    def test_malformed_native_grids_and_unknown_estimators_fail_closed(self) -> None:
        proper = _time()
        signal = np.zeros_like(proper)
        folded = proper.copy()
        folded[7] = folded[6]
        with self.assertRaisesRegex(ValueError, "must increase"):
            native_tail_observation(
                folded, signal, estimator=TDG3_ESTIMATOR_NAMES[0]
            )
        with self.assertRaisesRegex(ValueError, "unknown TDG3 native estimator"):
            native_tail_observation(proper, signal, estimator="spline")
        with self.assertRaisesRegex(ValueError, "method ladder differs"):
            diagnose_native_tail_ladder(
                (proper, proper, proper),
                (signal, signal, signal),
                method="RK4",
                point_counts=(1025, 2049, 4097),
            )

    def test_observation_cannot_promote_a_continuum_enclosure(self) -> None:
        proper = _time()
        signal = np.sin(2.0 * pi * 29.0 * proper / proper[-1])
        observed = native_tail_observation(
            proper, signal, estimator=TDG3_ESTIMATOR_NAMES[0]
        )
        with self.assertRaisesRegex(ValueError, "continuum enclosure"):
            replace(observed, continuum_surrogate_enclosure_claimed=True)


if __name__ == "__main__":
    unittest.main()
