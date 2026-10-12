from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    compact_vacuum_buffer_window,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    FGCQR_BRANCH_DEFINITION_OWNER,
    PROTO4_CONSTRAINT_ORDER,
    PROTO4_SPECTRAL_FIELD_ORDER,
    BranchStopApplicability,
    ConstraintEventSample,
    SpectralPowerBudget,
    branch_stop_applicability,
    causal_past_window,
    common_event_constraint_admission,
    conservative_observable_error_sum,
    conservative_richardson_error,
    minkowski_reference_deviations,
    nested_spectral_admission,
    normalized_constraint_norms,
    proper_radial_profile,
    spectral_field_budgets,
    windowed_spectral_power_budget,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    PROTO4_CANDIDATE_BRANCH_STOP_IDS,
    PROTO4_REPLACED_STOP_IDS,
    PROTO4_UNIVERSAL_STOP_IDS,
)


class FGCProto4AdmissionTests(unittest.TestCase):
    def test_branch_applicability_is_fail_closed_and_not_cross_owned(self) -> None:
        gr0 = branch_stop_applicability("GR-0")
        self.assertIsInstance(gr0, BranchStopApplicability)
        self.assertEqual(gr0.universal_stop_ids, PROTO4_UNIVERSAL_STOP_IDS)
        self.assertEqual(gr0.candidate_branch_stop_ids, ())
        self.assertEqual(gr0.replaced_stop_ids, PROTO4_REPLACED_STOP_IDS)
        self.assertTrue(gr0.runtime_monitor_complete)

        fgc_missing = branch_stop_applicability("FGC-QR")
        self.assertFalse(fgc_missing.runtime_monitor_complete)
        self.assertEqual(fgc_missing.candidate_branch_stop_ids, ())
        fgc = branch_stop_applicability(
            "FGC-QR", candidate_definition_owner=FGCQR_BRANCH_DEFINITION_OWNER
        )
        self.assertTrue(fgc.runtime_monitor_complete)
        self.assertEqual(
            fgc.candidate_branch_stop_ids, PROTO4_CANDIDATE_BRANCH_STOP_IDS
        )

        sgbl = branch_stop_applicability("SGB-L")
        self.assertFalse(sgbl.runtime_monitor_complete)
        self.assertIn("SGB-L", sgbl.missing_premise or "")
        with self.assertRaisesRegex(ValueError, "no certified"):
            branch_stop_applicability(
                "SGB-L", candidate_definition_owner=FGCQR_BRANCH_DEFINITION_OWNER
            )
        with self.assertRaisesRegex(ValueError, "GR-0 cannot"):
            branch_stop_applicability(
                "GR-0", candidate_definition_owner=FGCQR_BRANCH_DEFINITION_OWNER
            )

    def test_reference_deviations_use_exact_center_limits(self) -> None:
        radii = np.linspace(0.0, 2.0, 33)
        u = np.zeros((radii.size, 6))
        q = np.zeros_like(u)
        u[:, 0] = 1.0
        u[:, 1] = 0.25 * radii
        u[:, 2] = 1.0
        u[:, 3] = radii
        u[:, 4] = 8.0
        u[:, 5] = -4.0
        q[:, 1] = 0.25
        q[:, 3] = 1.0
        deviations = minkowski_reference_deviations(u, q, radii, cutoff=16.0)
        np.testing.assert_allclose(deviations[:, 0], 0.0, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(deviations[:, 1], 0.25, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(deviations[:, 2], 0.0, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(deviations[:, 3], 0.0, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(deviations[:, 4], 0.5, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(deviations[:, 5], -0.25, rtol=0.0, atol=0.0)

        u[:, 2] = 1.0 + 0.05 * radii
        profile = proper_radial_profile(
            u,
            q,
            radii,
            cutoff=16.0,
            measurement_radius_maximum=2.0,
        )
        self.assertEqual(profile.deviations.shape, (33, 6))
        self.assertTrue(np.all(np.diff(profile.proper_coordinates) > 0.0))
        self.assertGreater(profile.proper_coordinates[-1], 2.0)
        self.assertTrue(np.all(profile.round_trip_interpolation_infinity >= 0.0))

    def test_weighted_spectrum_rejects_power_not_last_bin_noise(self) -> None:
        coordinate = np.arange(256, dtype=np.float64) * (2.0 * np.pi / 256.0)
        low = np.sin(2.0 * coordinate)
        high = np.sin(112.0 * coordinate)
        low_budget = windowed_spectral_power_budget(
            low, coordinate, cutoff=16.0
        )
        high_budget = windowed_spectral_power_budget(
            high, coordinate, cutoff=16.0
        )
        self.assertTrue(low_budget.individual_admission_passed)
        self.assertFalse(high_budget.individual_admission_passed)
        self.assertGreater(
            high_budget.top_band_field_power_fraction,
            low_budget.top_band_field_power_fraction,
        )
        self.assertGreater(high_budget.rms_scale_over_cutoff, 0.25)

        strict_boundary = windowed_spectral_power_budget(
            np.sin(4.0 * coordinate), coordinate, cutoff=16.0
        )
        self.assertAlmostEqual(strict_boundary.rms_scale_over_cutoff, 0.25)
        self.assertFalse(strict_boundary.individual_admission_passed)

        compact = np.zeros_like(coordinate)
        interior = np.abs(coordinate - np.pi) < 1.0
        x = coordinate[interior] - np.pi
        compact[interior] = np.exp(1.0 - 1.0 / (1.0 - x * x))
        compact_budget = windowed_spectral_power_budget(
            compact, coordinate, cutoff=16.0
        )
        self.assertTrue(compact_budget.last_bin_alias_diagnostic)
        self.assertTrue(compact_budget.individual_admission_passed)

    def test_causal_temporal_window_requires_only_past_and_64_samples(self) -> None:
        times = np.linspace(-4.0, 0.0, 64)
        window = causal_past_window(times)
        self.assertEqual(window[0], 0.0)
        self.assertEqual(window[-1], 0.0)
        self.assertGreater(float(np.max(window)), 0.99)
        with self.assertRaisesRegex(ValueError, "64"):
            causal_past_window(times[:-1])

    @staticmethod
    def _budget(field: float, derivative: float, *, pass_individual: bool = True) -> SpectralPowerBudget:
        return SpectralPowerBudget(
            sample_count=129,
            top_band_first_bin=56,
            nyquist_bin=64,
            total_field_power=1.0,
            top_band_field_power_fraction=field,
            total_derivative_weighted_power=1.0,
            top_band_derivative_weighted_power_fraction=derivative,
            rms_angular_scale=1.0,
            rms_scale_over_cutoff=1.0 / 16.0,
            last_occupied_bin=64,
            last_bin_alias_diagnostic=True,
            individual_admission_passed=pass_individual,
        )

    def test_nested_spectral_tail_requires_decay_for_every_field(self) -> None:
        records = []
        for scale in (1.0, 0.125, 0.015625):
            records.append(
                {
                    name: self._budget(1.0e-8 * scale, 1.0e-6 * scale)
                    for name in PROTO4_SPECTRAL_FIELD_ORDER
                }
            )
        admitted = nested_spectral_admission((129, 257, 513), records)
        self.assertTrue(admitted.admission_passed)
        self.assertTrue(admitted.every_nested_tail_ratio_passed)

        broken = [dict(record) for record in records]
        broken[-1]["chi_over_Lambda"] = self._budget(1.0e-8, 1.0e-6)
        rejected = nested_spectral_admission((129, 257, 513), broken)
        self.assertFalse(rejected.admission_passed)
        self.assertFalse(rejected.every_nested_tail_ratio_passed)

        data = np.zeros((129, 6))
        coordinate = np.linspace(0.0, 8.0, 129)
        budgets = spectral_field_budgets(
            data,
            coordinate,
            cutoff=16.0,
            window=compact_vacuum_buffer_window(
                coordinate,
                support_minimum=0.0,
                support_maximum=8.0,
                taper_width=1.0,
            ),
        )
        self.assertEqual(tuple(budgets), PROTO4_SPECTRAL_FIELD_ORDER)
        self.assertTrue(all(value.individual_admission_passed for value in budgets.values()))

    @staticmethod
    def _constraint_sample(point_count: int, residual_scale: float, *, time: float = 0.0) -> ConstraintEventSample:
        residuals = np.full((point_count, len(PROTO4_CONSTRAINT_ORDER)), residual_scale)
        terms = np.zeros((point_count, len(PROTO4_CONSTRAINT_ORDER), 3))
        norms = normalized_constraint_norms(
            residuals,
            terms,
            planck_mass=2.0,
            length_unit=4.0,
        )
        return ConstraintEventSample(point_count, time, norms)

    def test_constraint_normalization_and_common_event_convergence(self) -> None:
        residuals = np.asarray([[1.0, 2.0], [3.0, 4.0]])
        terms = np.asarray(
            [
                [[1.0, -2.0], [0.0, 0.0]],
                [[0.0, 0.0], [1.0, -3.0]],
            ]
        )
        norms = normalized_constraint_norms(
            residuals,
            terms,
            planck_mass=2.0,
            length_unit=4.0,
            component_names=("a", "b"),
        )
        # max term sums are 3 and 4; the physical floor is 1/4.
        np.testing.assert_allclose(norms.component_infinity, (12.0, 8.0))
        self.assertEqual(norms.minimum_denominator, 0.25)
        self.assertEqual(norms.maximum_denominator, 4.0)

        samples = tuple(
            self._constraint_sample(n, (1.0 / (n - 1)) ** 2)
            for n in (33, 65, 129)
        )
        admitted = common_event_constraint_admission(samples)
        self.assertTrue(admitted.common_event_alignment_passed)
        self.assertTrue(admitted.coarsest_guard_passed)
        self.assertTrue(admitted.finest_guard_passed)
        self.assertTrue(admitted.convergence_order_passed)
        self.assertAlmostEqual(admitted.minimum_finite_observed_order or 0.0, 2.0)
        self.assertTrue(admitted.admission_passed)

        misaligned = list(samples)
        misaligned[-1] = self._constraint_sample(129, (1.0 / 128.0) ** 2, time=0.125)
        self.assertFalse(common_event_constraint_admission(misaligned).admission_passed)

        signed_zero = list(samples)
        signed_zero[-1] = self._constraint_sample(
            129, (1.0 / 128.0) ** 2, time=-0.0
        )
        self.assertFalse(
            common_event_constraint_admission(signed_zero).common_event_alignment_passed
        )

        with self.assertRaisesRegex(ValueError, "point_count must match"):
            ConstraintEventSample(130, 0.0, samples[-1].norms)

        # DEF1 cannot silently fall back to coordinate-time-only alignment.
        self.assertFalse(
            common_event_constraint_admission(
                samples,
                event_alignment_contract="DEF1_coordinate_time_affine_label_and_origin",
            ).common_event_alignment_passed
        )
        affine_samples = tuple(
            ConstraintEventSample(
                item.point_count,
                item.coordinate_time,
                item.norms,
                affine_label=1.5,
                affine_origin=0.0,
            )
            for item in samples
        )
        self.assertTrue(
            common_event_constraint_admission(
                affine_samples,
                event_alignment_contract="DEF1_coordinate_time_affine_label_and_origin",
            ).common_event_alignment_passed
        )

        nonconvergent = tuple(
            self._constraint_sample(n, 1.0e-4) for n in (33, 65, 129)
        )
        self.assertFalse(
            common_event_constraint_admission(nonconvergent).convergence_order_passed
        )
        self.assertEqual(
            common_event_constraint_admission(
                nonconvergent
            ).minimum_finite_observed_order,
            0.0,
        )

        exact = tuple(self._constraint_sample(n, 0.0) for n in (33, 65, 129))
        exact_result = common_event_constraint_admission(exact)
        self.assertTrue(exact_result.admission_passed)
        self.assertIsNone(exact_result.minimum_finite_observed_order)

    def test_richardson_and_error_sum_are_additive_only(self) -> None:
        error = conservative_richardson_error(
            np.asarray((1.0, 2.0)),
            np.asarray((1.25, 1.5)),
            refinement_ratio=2.0,
            assumed_order=2.0,
            event_alignment=common_event_constraint_admission(
                tuple(
                    self._constraint_sample(n, (1.0 / (n - 1)) ** 2)
                    for n in (33, 65, 129)
                )
            ),
        )
        self.assertAlmostEqual(error, 1.0 / 6.0)
        self.assertAlmostEqual(
            conservative_observable_error_sum(
                {"constraint": error, "interpolation": 0.25, "affine": 0.125}
            ),
            error + 0.375,
        )
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            conservative_observable_error_sum(
                {"constraint": error, "interpolation": -0.25}
            )

        misaligned = tuple(
            self._constraint_sample(
                n,
                (1.0 / (n - 1)) ** 2,
                time=(0.125 if n == 129 else 0.0),
            )
            for n in (33, 65, 129)
        )
        with self.assertRaisesRegex(ValueError, "common-event"):
            conservative_richardson_error(
                np.asarray((1.0, 2.0)),
                np.asarray((1.25, 1.5)),
                refinement_ratio=2.0,
                assumed_order=2.0,
                event_alignment=common_event_constraint_admission(misaligned),
            )


if __name__ == "__main__":
    unittest.main()
