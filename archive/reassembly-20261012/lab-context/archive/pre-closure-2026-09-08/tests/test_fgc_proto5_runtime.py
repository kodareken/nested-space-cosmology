from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    GR0_UNIVERSAL_RUNTIME_STOP_IDS,
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    EvolutionRHS,
    EvolutionState,
    PRIMARY_METHOD,
    UniformRadialGrid,
    propose_step,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    CFLRetryRequired,
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    GR0UniversalThresholds,
    GR0_RUNTIME_STOP_PRIORITY,
    gr0_hat_lorentzian_margin,
    proto5_gr0_common_event,
    proto5_temporal_spectral_admission,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


def _synthetic_state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _synthetic_rhs(*, bad_time: float | None = None):
    def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
        zero = np.zeros(state.shape)
        diagnostics = {
            "source_residual_infinity": 0.0,
            "source_refinement_iterations": 0,
            "source_residual_decreased_monotonically": True,
            "kinetic_condition_infinity": 1.0,
            "coordinate_speed_upper": 1.0,
            "minimum_lapse": 0.0 if bad_time is not None and time >= bad_time else 1.0,
            "minimum_radial_metric": 1.0,
            "minimum_areal_radius_away_from_center": 0.125,
        }
        return EvolutionRHS(zero, zero, zero, diagnostics)

    return rhs


class FGCProto5RuntimeTests(unittest.TestCase):
    def transaction(self) -> GR0RuntimeStageTransaction:
        return GR0RuntimeStageTransaction(
            thresholds=GR0UniversalThresholds(),
            causal_state=CausalBudgetState(previous_speed_upper=1.0),
            boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
            grid_spacing=0.125,
            cfl_maximum=0.125,
        )

    def test_GR0_partition_is_exactly_the_nine_universal_stops(self) -> None:
        self.assertEqual(
            GR0_RUNTIME_STOP_PRIORITY,
            tuple(GR0_UNIVERSAL_RUNTIME_STOP_IDS),
        )
        self.assertNotIn("branch_displacement_limit", GR0_RUNTIME_STOP_PRIORITY)
        self.assertNotIn("phi_cutoff_proxy_limit", GR0_RUNTIME_STOP_PRIORITY)

    def test_hat_margin_is_positive_on_flat_regular_state(self) -> None:
        self.assertEqual(gr0_hat_lorentzian_margin(_synthetic_state()), 1.0)

    def test_runtime_transaction_accepts_atomically_and_latches_first_failure(self) -> None:
        state = _synthetic_state()
        healthy = propose_step(
            method=PRIMARY_METHOD,
            time=0.0,
            step_size=0.01,
            state=state,
            rhs=_synthetic_rhs(),
        )
        transaction = self.transaction()
        receipt = transaction(healthy)
        self.assertEqual(receipt["accepted_stage_count"], 5)
        self.assertEqual(transaction.state.accepted_stage_count, 5)

        failing = propose_step(
            method=PRIMARY_METHOD,
            time=0.01,
            step_size=0.01,
            state=state,
            rhs=_synthetic_rhs(bad_time=0.02),
        )
        causal_before = transaction.causal_state
        with self.assertRaises(GR0RuntimeStop) as caught:
            transaction(failing)
        self.assertEqual(caught.exception.reason, "nonpositive_lapse")
        self.assertEqual(transaction.state.accepted_stage_count, 5)
        self.assertEqual(transaction.causal_state, causal_before)
        with self.assertRaises(GR0RuntimeStop) as repeated:
            transaction(failing)
        self.assertEqual(repeated.exception.reason, "nonpositive_lapse")

    def test_CFL_rejection_is_retryable_and_not_a_scientific_stop(self) -> None:
        proposal = propose_step(
            method=COMPARATOR_METHOD,
            time=0.0,
            step_size=0.02,
            state=_synthetic_state(),
            rhs=_synthetic_rhs(),
        )
        transaction = self.transaction()
        causal_before = transaction.causal_state
        with self.assertRaises(CFLRetryRequired):
            transaction(proposal)
        self.assertIsNone(transaction.state.first_failed_premise)
        self.assertEqual(transaction.state.accepted_stage_count, 0)
        self.assertEqual(transaction.causal_state, causal_before)

    def test_actual_projected_t0_common_event_passes_primary_contract(self) -> None:
        states = []
        grids = []
        for count in (1025, 2049, 4097):
            initial = construct_gr0_grid_initial_data(
                PulseParameters(chi_amplitude=2.5),
                point_count=count,
                constraint_method="RK4",
                diagnostic_spatial_order=4,
            )
            states.append(project_gr0_semidiscrete_state(initial, spatial_order=4))
            grids.append(initial.grid)
        assessed = proto5_gr0_common_event(
            states,
            grids,
            method="RK4",
            coordinate_time=0.0,
        )
        self.assertTrue(assessed.constraint_admission.admission_passed)
        self.assertTrue(assessed.spatial_spectral_admission.admission_passed)
        self.assertTrue(assessed.admission_passed)

    def test_zero_temporal_histories_pass_without_invented_infinite_order(self) -> None:
        samples = 64
        tracers = 2
        fields = len(
            (
                "alpha_minus_1",
                "shift_over_r",
                "lambda_minus_1",
                "R_over_r_minus_1",
                "phi_over_Lambda",
                "chi_over_Lambda",
            )
        )
        times = np.repeat(np.linspace(0.0, 1.0, samples)[:, None], tracers, axis=1)
        histories = np.zeros((samples, tracers, fields))
        assessed = proto5_temporal_spectral_admission(
            point_counts=(1025, 2049, 4097),
            proper_times_by_resolution=(times, times, times),
            field_histories_by_resolution=(histories, histories, histories),
        )
        self.assertTrue(assessed.admission_passed)
        self.assertEqual(assessed.maximum_round_trip_interpolation_infinity, 0.0)


if __name__ == "__main__":
    unittest.main()
