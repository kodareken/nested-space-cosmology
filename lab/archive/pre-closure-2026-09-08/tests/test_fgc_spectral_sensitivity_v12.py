from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
    ProperRadialProfile,
    SpectralPowerBudget,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity import (  # noqa: E402
    SpectralTailSensitivity,
    resolved_or_saturated_admission,
    three_grid_profile_convergence,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity_v12 import (  # noqa: E402
    pairwise_resolved_or_saturated_admission,
)


class FGCSpectralSensitivityV12Tests(unittest.TestCase):
    @staticmethod
    def _profiles() -> tuple[ProperRadialProfile, ...]:
        profiles = []
        for count, error in ((33, 0.04), (65, 0.02), (129, 0.01)):
            coordinate = np.linspace(0.0, 2.0, count)
            base = np.sin(coordinate)
            deviations = np.zeros((count, len(PROTO4_SPECTRAL_FIELD_ORDER)))
            deviations[:, 2] = base + error * np.cos(3.0 * coordinate)
            deviations[:, 4] = 0.5 * base + 0.25 * error * np.cos(3.0 * coordinate)
            deviations[:, 5] = 0.25 * base + 0.5 * error * np.cos(3.0 * coordinate)
            round_trip = np.array((0.0, 0.0, error, 0.0, error, error))
            profiles.append(ProperRadialProfile(coordinate, deviations, round_trip))
        return tuple(profiles)

    @staticmethod
    def _budget(fraction: float) -> SpectralPowerBudget:
        return SpectralPowerBudget(
            sample_count=129,
            top_band_first_bin=56,
            nyquist_bin=64,
            total_field_power=1.0,
            top_band_field_power_fraction=fraction,
            total_derivative_weighted_power=1.0,
            top_band_derivative_weighted_power_fraction=fraction,
            rms_angular_scale=1.0,
            rms_scale_over_cutoff=1.0 / 16.0,
            last_occupied_bin=56,
            last_bin_alias_diagnostic=True,
            individual_admission_passed=True,
        )

    @staticmethod
    def _sensitivity() -> SpectralTailSensitivity:
        return SpectralTailSensitivity(
            sample_count=129,
            top_band_first_bin=56,
            nyquist_angular_frequency=64.0,
            round_trip_interpolation_infinity=1.0e-4,
            top_band_field_rms_amplitude=1.0e-6,
            top_band_derivative_rms_amplitude=1.0e-5,
            top_band_erasure_perturbation_infinity=1.0e-6,
            erasure_over_round_trip=0.01,
            field_tail_within_round_trip_scale=True,
            derivative_tail_within_nyquist_round_trip_scale=True,
            erasure_within_round_trip_scale=True,
            diagnostically_saturated=True,
        )

    def test_both_pairs_may_saturate_without_hiding_raw_ratios(self) -> None:
        names = tuple(PROTO4_SPECTRAL_FIELD_ORDER)
        budgets = tuple(
            {name: self._budget(fraction) for name in names}
            for fraction in (1.0e-10, 2.0e-10, 4.0e-10)
        )
        sensitivities = tuple(
            {name: self._sensitivity() for name in names} for _ in range(3)
        )
        convergence = three_grid_profile_convergence(
            (2049, 4097, 8193), self._profiles()
        )
        old = resolved_or_saturated_admission(
            (2049, 4097, 8193), budgets, sensitivities, convergence
        )
        pairwise = pairwise_resolved_or_saturated_admission(
            (2049, 4097, 8193), budgets, sensitivities, convergence
        )
        self.assertFalse(old.admission_passed)
        self.assertTrue(pairwise.admission_passed)
        self.assertTrue(pairwise.saturation_used)
        self.assertEqual(
            pairwise.field_power_tail_ratios["alpha_minus_1"], (2.0, 2.0)
        )
        self.assertEqual(
            pairwise.field_pair_classification["alpha_minus_1"],
            ("diagnostically_saturated", "diagnostically_saturated"),
        )

    def test_each_pair_requires_its_budget_and_map_witnesses(self) -> None:
        names = tuple(PROTO4_SPECTRAL_FIELD_ORDER)
        budgets = [
            {name: self._budget(fraction) for name in names}
            for fraction in (1.0e-10, 2.0e-10, 4.0e-10)
        ]
        sensitivities = [
            {name: self._sensitivity() for name in names} for _ in range(3)
        ]
        convergence = three_grid_profile_convergence(
            (2049, 4097, 8193), self._profiles()
        )

        budgets[0]["alpha_minus_1"] = replace(
            budgets[0]["alpha_minus_1"], individual_admission_passed=False
        )
        rejected_budget = pairwise_resolved_or_saturated_admission(
            (2049, 4097, 8193), budgets, sensitivities, convergence
        )
        self.assertFalse(rejected_budget.admission_passed)
        self.assertEqual(
            rejected_budget.field_pair_classification["alpha_minus_1"][0],
            "failed",
        )

        budgets[0]["alpha_minus_1"] = self._budget(1.0e-10)
        sensitivities[1]["alpha_minus_1"] = replace(
            sensitivities[1]["alpha_minus_1"], diagnostically_saturated=False
        )
        rejected_map = pairwise_resolved_or_saturated_admission(
            (2049, 4097, 8193), budgets, sensitivities, convergence
        )
        self.assertFalse(rejected_map.admission_passed)
        self.assertEqual(
            rejected_map.field_pair_classification["alpha_minus_1"],
            ("failed", "failed"),
        )


if __name__ == "__main__":
    unittest.main()
