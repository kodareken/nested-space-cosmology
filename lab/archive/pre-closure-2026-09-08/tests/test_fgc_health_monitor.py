from __future__ import annotations

from dataclasses import replace
from math import pi
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt1_mon1 import _baseline_snapshot, load_config  # noqa: E402
from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    ClassicalHealthMonitor,
    ClassicalHealthStop,
    compact_vacuum_buffer_window,
    curvature_scale_and_ratio,
    dimensionless_state_components,
    induced_scaled_infinity_norm,
    normalized_componentwise_infinity,
    windowed_spectral_support,
)


class FGCHealthMonitorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.thresholds = load_config()["thresholds"]

    def test_canonical_norms_and_curvature_ratio_are_executable(self) -> None:
        u, p, q = dimensionless_state_components(
            np.ones(6), np.ones(6), 2.0 * np.ones(6), length_unit=4.0
        )
        np.testing.assert_array_equal(u, (1, 1, 1, 0.25, 4, 4))
        np.testing.assert_array_equal(p, (4, 4, 4, 1, 16, 16))
        np.testing.assert_array_equal(q, (8, 8, 8, 2, 32, 32))
        self.assertEqual(
            induced_scaled_infinity_norm(
                ((1, -2), (3, 4)), row_scales=(2, 0.5), column_scales=(1, 2)
            ),
            4.0,
        )
        self.assertEqual(normalized_componentwise_infinity((0, 2e-8), (0, 2)), 1e-8)
        kappa, ratio = curvature_scale_and_ratio((256,), (4,), cutoff=16)
        self.assertEqual(kappa, 4.0)
        self.assertEqual(ratio, 1.0 / 16.0)

    def test_proper_spectral_support_and_alias_band(self) -> None:
        n = 128
        x = np.arange(n) / n
        low = windowed_spectral_support(np.sin(2 * pi * 2 * x), x)
        high = windowed_spectral_support(np.sin(2 * pi * 60 * x), x)
        self.assertEqual(low.support_bin, 2)
        self.assertFalse(low.alias_band_occupied)
        self.assertEqual(high.support_bin, 60)
        self.assertTrue(high.alias_band_occupied)
        window = compact_vacuum_buffer_window(
            x, support_minimum=0, support_maximum=(n - 1) / n, taper_width=1 / 8
        )
        self.assertEqual(window[0], 0.0)
        self.assertEqual(window[-1], 0.0)
        self.assertEqual(np.max(window), 1.0)

    def test_transaction_stops_at_equality_and_preserves_first_failure(self) -> None:
        monitor = ClassicalHealthMonitor(self.thresholds)
        accepted = monitor.assess_and_accept(_baseline_snapshot())
        self.assertEqual(accepted.accepted_stage_count, 1)
        trial = replace(
            _baseline_snapshot(),
            time=0.25,
            stage_index=2,
            newton_residual_infinity=self.thresholds.newton_residual_max,
        )
        with self.assertRaises(ClassicalHealthStop) as caught:
            monitor.assess_and_accept(trial)
        self.assertEqual(caught.exception.reason, "newton_residual_limit")
        stopped = monitor.state
        self.assertEqual(stopped.accepted_stage_count, 1)
        self.assertEqual(stopped.last_accepted_stage_index, 1)
        with self.assertRaises(ClassicalHealthStop) as repeated:
            monitor.assess_and_accept(replace(_baseline_snapshot(), time=0.375, stage_index=3))
        self.assertEqual(repeated.exception.reason, "newton_residual_limit")
        self.assertEqual(repeated.exception.state, stopped)

    def test_global_transaction_order_supports_SSPRK3_stage_abscissae(self) -> None:
        monitor = ClassicalHealthMonitor(self.thresholds)
        for serial, stage_time in enumerate((0.0, 1.0, 0.5, 1.0), start=1):
            state = monitor.assess_and_accept(
                replace(_baseline_snapshot(), time=stage_time, stage_index=serial)
            )
        self.assertEqual(state.accepted_stage_count, 4)
        self.assertEqual(state.last_accepted_time, 1.0)
        before = monitor.state
        with self.assertRaisesRegex(ValueError, "strictly increasing global"):
            monitor.assess_and_accept(
                replace(_baseline_snapshot(), time=2.0, stage_index=4)
            )
        self.assertEqual(monitor.state, before)


if __name__ == "__main__":
    unittest.main()
