from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg1_temporal_diagnosis import (  # noqa: E402
    TDG1_FIELD_NAMES,
    TDG1_SAMPLE_COUNT,
    TDG1_SYNTHETIC_CONTROL_NAMES,
    diagnose_terminal_method_histories,
    synthetic_temporal_gate_controls,
    temporal_history_difference_diagnosis,
    temporal_signal_views,
)


class TDG1TemporalDiagnosisTests(unittest.TestCase):
    def test_synthetic_controls_expose_normalization_cancellation(self) -> None:
        record = synthetic_temporal_gate_controls()
        self.assertTrue(record["structural_counterexample_confirmed"])
        self.assertEqual(
            tuple(record["controls"]), TDG1_SYNTHETIC_CONTROL_NAMES
        )
        self.assertTrue(record["controls"]["zero"]["admission_passed"])
        for name in TDG1_SYNTHETIC_CONTROL_NAMES[1:]:
            control = record["controls"][name]
            self.assertTrue(control["every_individual_budget_passed"])
            self.assertFalse(control["every_nested_tail_ratio_passed"])
            self.assertFalse(control["admission_passed"])
            self.assertEqual(control["all_field_and_derivative_ratios_equal"], 1.0)
            self.assertTrue(control["ratio_identity_verified"])
        scaled = record["controls"]["amplitude_scaled_low_sine"]
        self.assertTrue(scaled["amplitude_error_contracts_by_four_per_resolution"])
        self.assertTrue(scaled["power_error_contracts_by_sixteen_per_resolution"])

    def test_signal_views_keep_every_transform_and_debit_public(self) -> None:
        proper = np.linspace(0.0, 63.0 / 16.0, TDG1_SAMPLE_COUNT)
        values = 0.125 + 0.01 * proper
        record = temporal_signal_views(proper, values)
        self.assertEqual(
            tuple(record["views"]),
            (
                "inherited_compact",
                "rectangular",
                "mean_centered_compact",
                "affine_detrended_compact",
            ),
        )
        self.assertEqual(
            tuple(record["resampling_sensitivity_views"]), ("16", "32", "64")
        )
        self.assertEqual(record["round_trip_interpolation_infinity"], 0.0)
        self.assertTrue(
            record[
                "resampling_views_are_not_independent_temporal_convergence_evidence"
            ]
        )
        # Closed-form detrending leaves only binary64-scale residue.  The
        # inherited normalized budget still rejects that residue; TDG1 keeps
        # both facts public and does not convert this diagnostic view into a
        # replacement pass.
        detrended = record["views"]["affine_detrended_compact"]
        self.assertLess(
            detrended["windowed_FFT_peak_amplitude"],
            2.0 * np.finfo(np.float64).eps,
        )
        self.assertFalse(detrended["budget"]["individual_admission_passed"])

    def test_raw_history_difference_order_is_separate_from_interpolation(self) -> None:
        proper = np.linspace(0.0, 1.0, TDG1_SAMPLE_COUNT)
        shape = np.sin(2.0 * np.pi * proper)
        records = (shape, shape + 0.25 * shape, shape + 0.3125 * shape)
        diagnosis = temporal_history_difference_diagnosis(
            (proper, proper, proper),
            records,
            point_counts=(2049, 4097, 8193),
        )
        self.assertEqual(diagnosis["classification"], "raw_order_pass")
        self.assertAlmostEqual(diagnosis["raw_observed_order"], 2.0)
        self.assertTrue(diagnosis["raw_order_passed"])
        self.assertFalse(diagnosis["interpolation_limited"])
        self.assertTrue(diagnosis["diagnostic_is_not_a_replacement_admission"])

        identical = temporal_history_difference_diagnosis(
            (proper, proper, proper),
            (shape, shape, shape),
            point_counts=(2049, 4097, 8193),
        )
        self.assertEqual(
            identical["classification"], "exactly_resolution_identical"
        )
        self.assertTrue(identical["raw_order_passed"])

    def test_method_diagnosis_uses_the_complete_frozen_shape(self) -> None:
        proper = np.repeat(
            np.linspace(0.0, 63.0 / 16.0, TDG1_SAMPLE_COUNT)[:, None],
            48,
            axis=1,
        )
        normalized = proper[:, 0] / proper[-1, 0]
        base = np.broadcast_to(
            (0.125 * np.sin(2.0 * np.pi * normalized))[:, None, None],
            (TDG1_SAMPLE_COUNT, 48, len(TDG1_FIELD_NAMES)),
        ).copy()
        record = diagnose_terminal_method_histories(
            method="RK4",
            point_counts=(2049, 4097, 8193),
            proper_times_by_resolution=(proper, proper, proper),
            field_histories_by_resolution=(base, base, base),
        )
        self.assertEqual(record["method"], "RK4")
        self.assertEqual(record["sample_count"], 64)
        self.assertEqual(record["tracer_count"], 48)
        self.assertEqual(
            record["history_difference_classification_counts"]
            ["exactly_resolution_identical"],
            48 * len(TDG1_FIELD_NAMES),
        )
        self.assertTrue(record["diagnostic_does_not_authorize_a_replacement_gate"])

    def test_malformed_inputs_fail_closed(self) -> None:
        proper = np.linspace(0.0, 1.0, TDG1_SAMPLE_COUNT)
        values = np.zeros_like(proper)
        with self.assertRaisesRegex(ValueError, "64"):
            temporal_signal_views(proper[:-1], values[:-1])
        with self.assertRaisesRegex(ValueError, "point counts"):
            temporal_history_difference_diagnosis(
                (proper, proper, proper),
                (values, values, values),
                point_counts=(2049, 4097, 4097),
            )
        bad = proper.copy()
        bad[-1] = bad[-2]
        with self.assertRaisesRegex(ValueError, "increasing"):
            temporal_history_difference_diagnosis(
                (proper, bad, proper),
                (values, values, values),
                point_counts=(2049, 4097, 8193),
            )


if __name__ == "__main__":
    unittest.main()
