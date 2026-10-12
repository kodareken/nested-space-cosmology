from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt10_mon10 as hlt10  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    SBPFirstDerivative,
    UniformRadialGrid,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    exact_spherical_minkowski_reference,
    project_gr0_reference_balanced_state,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    PROTO12_POINT_COUNTS,
    Proto12GR0EvolutionOperator,
    diagnose_gr0_proto12_accelerations_chunked,
    proto12_gr0_common_event,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity_v12 import (  # noqa: E402
    pairwise_resolved_or_saturated_admission,
)
from recursive_horizons.fgc.evolution.src4_vectorized_reference_source import (  # noqa: E402
    diagnose_gr0_vectorized_reference_accelerations,
)


class FGCProto12RuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.plan = hlt10.validate_run_plan(
            REPOSITORY / "configs/fgc/fgc-1-cal9-run1.toml"
        )
        cls.initials = [
            hlt10._initial_data(cls.plan, "5/2", "RK4", count)
            for count in PROTO12_POINT_COUNTS
        ]
        cls.states = [
            project_gr0_reference_balanced_state(initial, spatial_order=4)
            for initial in cls.initials
        ]
        cls.event = proto12_gr0_common_event(
            cls.states,
            [initial.grid for initial in cls.initials],
            accepted_stage_counts=(0, 0, 0),
            method="RK4",
            coordinate_time=0.0,
            cutoff=16.0,
            measurement_radius_maximum=24.0,
            taper_fraction=1.0 / 8.0,
            fixed_outer_rows=4,
        )

    def test_exact_reference_is_bitwise_stationary_through_SRC4(self) -> None:
        grid = UniformRadialGrid(0.0, 128.0, 257)
        state = exact_spherical_minkowski_reference(grid)
        before = array_content_sha256(state.u, state.p, state.q)
        operator = Proto12GR0EvolutionOperator(
            grid,
            spatial_order=4,
            ko_dissipation=1.0 / 64.0,
            raw_tolerance=1.0e-12,
            kinetic_condition_maximum=1.0e10,
            point_batch_size=64,
        )
        rhs = operator(0.0, state)
        self.assertEqual(before, array_content_sha256(state.u, state.p, state.q))
        self.assertTrue(np.array_equal(rhs.du, np.zeros_like(rhs.du)))
        self.assertTrue(np.array_equal(rhs.dp, np.zeros_like(rhs.dp)))
        self.assertTrue(np.array_equal(rhs.dq, np.zeros_like(rhs.dq)))
        self.assertTrue(rhs.diagnostics["source_raw_gate_passed"])
        self.assertTrue(rhs.diagnostics["SRC4_tensor_contracted_reference_source"])

    def test_fixed_batching_is_identical_to_one_SRC4_call(self) -> None:
        initial = self.initials[0]
        state = self.states[0]
        derivative = SBPFirstDerivative(initial.grid, 4)
        p_r, q_r, _ = reference_balanced_spatial_derivatives(state, derivative)
        indices = np.arange(1, 38)
        direct = diagnose_gr0_vectorized_reference_accelerations(
            state.u[indices],
            state.p[indices],
            state.q[indices],
            p_r[indices],
            q_r[indices],
            initial.grid.coordinates[indices],
        )
        chunked = diagnose_gr0_proto12_accelerations_chunked(
            state.u[indices],
            state.p[indices],
            state.q[indices],
            p_r[indices],
            q_r[indices],
            initial.grid.coordinates[indices],
            point_batch_size=7,
        )
        self.assertTrue(np.array_equal(direct.accelerations, chunked.accelerations))
        self.assertTrue(np.array_equal(direct.residuals, chunked.residuals))
        self.assertEqual(direct.residual_infinity, chunked.residual_infinity)
        self.assertEqual(chunked.point_batch_count, 6)

    def test_pairwise_event_retains_raw_evidence_and_passes(self) -> None:
        event = self.event
        self.assertTrue(event.admission_passed)
        self.assertTrue(event.spatial_spectral_admission.admission_passed)
        self.assertTrue(event.spatial_spectral_admission.saturation_used)
        self.assertFalse(event.interior_q_reprojected)
        self.assertEqual(
            dict(event.raw_spatial_spectral_admission.field_power_tail_ratios),
            dict(event.spatial_spectral_admission.field_power_tail_ratios),
        )
        self.assertEqual(
            dict(
                event.raw_spatial_spectral_admission.derivative_power_tail_ratios
            ),
            dict(event.spatial_spectral_admission.derivative_power_tail_ratios),
        )

    def test_removing_one_used_pair_guard_fails_closed(self) -> None:
        event = self.event
        admission = event.spatial_spectral_admission
        selected: tuple[str, str, int] | None = None
        for name, classes in admission.field_pair_classification.items():
            for pair, value in enumerate(classes):
                if value == "diagnostically_saturated":
                    selected = ("field", name, pair)
                    break
            if selected is not None:
                break
        if selected is None:
            for name, classes in admission.derivative_pair_classification.items():
                for pair, value in enumerate(classes):
                    if value == "diagnostically_saturated":
                        selected = ("derivative", name, pair)
                        break
                if selected is not None:
                    break
        self.assertIsNotNone(selected)
        kind, name, pair = selected  # type: ignore[misc]
        witnesses = [dict(record) for record in event.spectral_tail_sensitivities]
        grid_index = pair
        witnesses[grid_index][name] = replace(
            witnesses[grid_index][name],
            field_tail_within_round_trip_scale=False,
            derivative_tail_within_nyquist_round_trip_scale=False,
            erasure_within_round_trip_scale=False,
            diagnostically_saturated=False,
        )
        attacked = pairwise_resolved_or_saturated_admission(
            PROTO12_POINT_COUNTS,
            event.spectral_budgets,
            witnesses,
            event.profile_convergence,
        )
        self.assertFalse(attacked.admission_passed)
        classifications = (
            attacked.field_pair_classification
            if kind == "field"
            else attacked.derivative_pair_classification
        )
        self.assertEqual(classifications[name][pair], "failed")

    def test_invalid_batch_and_wrong_ladder_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "point_batch_size"):
            Proto12GR0EvolutionOperator(
                UniformRadialGrid(0.0, 128.0, 65),
                spatial_order=4,
                ko_dissipation=0.0,
                raw_tolerance=1.0e-12,
                kinetic_condition_maximum=1.0e10,
                point_batch_size=2049,
            )
        with self.assertRaisesRegex(ValueError, "frozen ladder"):
            proto12_gr0_common_event(
                self.states[::-1],
                [initial.grid for initial in self.initials[::-1]],
                accepted_stage_counts=(0, 0, 0),
                method="RK4",
                coordinate_time=0.0,
            )


if __name__ == "__main__":
    unittest.main()
