from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
    gr0_semidiscrete_constraint_snapshot,
    project_gr0_semidiscrete_state,
    proto5_semidiscrete_constraint_admission,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.initial_data_preflight import PulseParameters  # noqa: E402


class FGCCalibrationRuntimeTests(unittest.TestCase):
    def snapshots(self, amplitude: float, method: str, order: int):
        records = []
        for point_count in (1025, 2049, 4097):
            initial = construct_gr0_grid_initial_data(
                PulseParameters(chi_amplitude=amplitude),
                point_count=point_count,
                constraint_method=method,
                diagnostic_spatial_order=order,
            )
            projected = project_gr0_semidiscrete_state(
                initial,
                spatial_order=order,
            )
            np.testing.assert_array_equal(projected.u, initial.state.u)
            np.testing.assert_array_equal(projected.p, initial.state.p)
            records.append(
                gr0_semidiscrete_constraint_snapshot(
                    projected,
                    initial.grid,
                    coordinate_time=0.0,
                    diagnostic_spatial_order=order,
                )
            )
        return tuple(records)

    def test_primary_eligible_inputs_pass_prospective_t0_contract(self) -> None:
        for amplitude in (2.5, 3.0):
            with self.subTest(amplitude=amplitude):
                result = proto5_semidiscrete_constraint_admission(
                    self.snapshots(amplitude, "RK4", 4),
                    method="RK4",
                    coarsest_guard_maximum=1.0 / 50.0,
                    finest_guard_maximum=1.0 / 1000.0,
                )
                self.assertTrue(result.admission_passed)
                self.assertGreaterEqual(
                    result.minimum_finite_finest_pair_order or 0.0,
                    1.5,
                )

    def test_comparator_uses_declared_lower_order_guard_without_hiding_raw_norms(self) -> None:
        records = self.snapshots(3.0, "SSPRK3", 2)
        result = proto5_semidiscrete_constraint_admission(
            records,
            method="SSPRK3",
            coarsest_guard_maximum=1.0 / 10.0,
            finest_guard_maximum=1.0 / 200.0,
        )
        self.assertTrue(result.admission_passed)
        self.assertGreater(records[0].raw_norms.global_infinity, 1.0 / 50.0)
        self.assertLess(records[-1].raw_norms.global_infinity, 1.0 / 200.0)

    def test_roundoff_enclosure_only_changes_order_classification(self) -> None:
        records = self.snapshots(2.5, "RK4", 4)
        enclosure = (
            CONSTRAINT_ROUNDOFF_OPERATION_BUDGET * np.finfo(np.float64).eps
        )
        for record in records:
            self.assertTrue(
                np.all(record.normalized_roundoff_zero_enclosure == enclosure)
            )
            self.assertGreaterEqual(
                record.raw_norms.global_infinity,
                record.effective_norms_for_order_only.global_infinity,
            )
        result = proto5_semidiscrete_constraint_admission(
            records,
            method="RK4",
            coarsest_guard_maximum=1.0 / 50.0,
            finest_guard_maximum=1.0 / 1000.0,
        )
        self.assertEqual(
            result.component_status["gauge_t"],
            "roundoff_enclosed_zero_pair",
        )

    def test_nonconvergent_component_cannot_hide_inside_roundoff_rule(self) -> None:
        records = list(self.snapshots(2.5, "RK4", 4))
        fine = records[-1]
        effective = fine.effective_norms_for_order_only.component_infinity.copy()
        effective[2] = 1.0e-4
        from recursive_horizons.fgc.evolution.proto4_admission import ConstraintNorms

        replacement = ConstraintNorms(
            fine.point_count,
            fine.effective_norms_for_order_only.component_names,
            effective,
            float(np.max(effective)),
            fine.effective_norms_for_order_only.minimum_denominator,
            fine.effective_norms_for_order_only.maximum_denominator,
        )
        from dataclasses import replace

        records[-1] = replace(fine, effective_norms_for_order_only=replacement)
        result = proto5_semidiscrete_constraint_admission(
            records,
            method="RK4",
            coarsest_guard_maximum=1.0 / 50.0,
            finest_guard_maximum=1.0 / 1000.0,
        )
        self.assertFalse(result.admission_passed)
        self.assertEqual(
            result.component_status["gauge_t"],
            "nonzero_reappeared_after_roundoff_enclosure",
        )


if __name__ == "__main__":
    unittest.main()
