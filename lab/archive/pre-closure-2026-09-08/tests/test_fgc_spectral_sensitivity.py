from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
    ProperRadialProfile,
    SpectralPowerBudget,
    windowed_spectral_power_budget,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity import (  # noqa: E402
    SpectralTailSensitivity,
    resolved_or_saturated_admission,
    spectral_tail_sensitivity,
    three_grid_profile_convergence,
)


class FGCSpectralSensitivityTests(unittest.TestCase):
    def test_top_band_erasure_exposes_ratio_sensitivity_without_becoming_error_bound(
        self,
    ) -> None:
        coordinate = np.linspace(0.0, 2.0 * np.pi, 257, endpoint=False)
        samples = np.sin(2.0 * coordinate) + 1.0e-8 * np.sin(112.0 * coordinate)
        window = np.ones_like(samples)
        budget = windowed_spectral_power_budget(
            samples,
            coordinate,
            cutoff=16.0,
            window=window,
        )
        witness = spectral_tail_sensitivity(
            samples,
            coordinate,
            window=window,
            round_trip_interpolation_infinity=1.0e-5,
            budget=budget,
        )
        self.assertTrue(witness.diagnostically_saturated)
        self.assertLess(witness.top_band_erasure_perturbation_infinity, 1.0e-5)
        self.assertGreater(witness.top_band_erasure_perturbation_infinity, 0.0)

        too_small = spectral_tail_sensitivity(
            samples,
            coordinate,
            window=window,
            round_trip_interpolation_infinity=(
                0.5 * witness.top_band_erasure_perturbation_infinity
            ),
            budget=budget,
        )
        self.assertFalse(too_small.diagnostically_saturated)

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

    def test_complete_profiles_and_round_trip_residuals_must_both_contract(
        self,
    ) -> None:
        convergence = three_grid_profile_convergence(
            (2049, 4097, 8193), self._profiles()
        )
        self.assertTrue(convergence.every_field_contracted)
        self.assertTrue(convergence.exact_zero_by_field["alpha_minus_1"])
        self.assertTrue(convergence.infinity_contracted_by_field["lambda_minus_1"])
        self.assertTrue(convergence.rms_contracted_by_field["chi_over_Lambda"])
        self.assertTrue(convergence.round_trip_contracted_by_field["phi_over_Lambda"])

        broken = list(self._profiles())
        final = broken[-1]
        deviations = final.deviations.copy()
        deviations[:, 2] += 0.2 * np.cos(3.0 * final.proper_coordinates)
        broken[-1] = ProperRadialProfile(
            final.proper_coordinates,
            deviations,
            final.round_trip_interpolation_infinity,
        )
        rejected = three_grid_profile_convergence((2049, 4097, 8193), broken)
        self.assertFalse(rejected.every_field_contracted)

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

    def test_direct_coarse_ratio_cannot_be_replaced_by_finest_saturation(self) -> None:
        names = tuple(PROTO4_SPECTRAL_FIELD_ORDER)
        budgets = tuple(
            {name: self._budget(fraction) for name in names}
            for fraction in (1.0e-8, 1.0e-10, 2.0e-10)
        )
        sensitivities = tuple(
            {name: self._sensitivity() for name in names} for _ in range(3)
        )
        convergence = three_grid_profile_convergence(
            (2049, 4097, 8193), self._profiles()
        )
        admitted = resolved_or_saturated_admission(
            (2049, 4097, 8193), budgets, sensitivities, convergence
        )
        self.assertTrue(admitted.admission_passed)
        self.assertTrue(admitted.saturation_used)
        self.assertTrue(
            all(
                value == "diagnostically_saturated"
                for value in admitted.finest_pair_field_classification.values()
            )
        )

        broken = list(budgets)
        broken[1] = {name: self._budget(1.0e-7) for name in names}
        rejected = resolved_or_saturated_admission(
            (2049, 4097, 8193), broken, sensitivities, convergence
        )
        self.assertFalse(rejected.admission_passed)
        self.assertFalse(
            all(rejected.coarse_to_medium_direct_passed_by_field.values())
        )


if __name__ == "__main__":
    unittest.main()
