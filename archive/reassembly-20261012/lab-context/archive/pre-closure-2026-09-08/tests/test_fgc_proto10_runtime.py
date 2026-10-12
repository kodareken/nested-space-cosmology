from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto10_runtime import (  # noqa: E402
    proto10_gr0_common_event,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
)
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity import (  # noqa: E402
    resolved_or_saturated_admission,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = Fraction


def _physical_event(amplitude: str, method: str):
    order = 4 if method == "RK4" else 2
    states = []
    grids = []
    for count in PROTO9_POINT_COUNTS:
        initial = construct_gr0_grid_initial_data(
            PulseParameters(
                chi_amplitude=float(Q(amplitude)),
                center=12.0,
                half_width=2.0,
                phi_amplitude=float(Q(1, 131072)),
                planck_mass=2.0,
                scalar_mass=3.0,
                quartic_coupling=0.5,
            ),
            point_count=count,
            outer_radius=128.0,
            constraint_method=method,
            diagnostic_spatial_order=order,
        )
        states.append(project_gr0_semidiscrete_state(initial, spatial_order=order))
        grids.append(initial.grid)
    return proto10_gr0_common_event(
        states,
        grids,
        accepted_stage_counts=(0, 0, 0),
        method=method,
        coordinate_time=0.0,
    )


def _flat_state(grid: UniformRadialGrid, order: int) -> EvolutionState:
    u = np.zeros((grid.point_count, 6), dtype=np.float64)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = grid.coordinates
    p = np.zeros_like(u)
    q = SBPFirstDerivative(grid, order).differentiate(
        u, center_parities=ADM_CENTER_PARITIES
    )
    return EvolutionState(u, p, q)


class FGCProto10RuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.events = {
            (amplitude, method): _physical_event(amplitude, method)
            for amplitude in ("5/2", "3")
            for method in ("RK4", "SSPRK3")
        }

    def test_all_four_physical_t0_events_preserve_raw_veto_and_pass_guards(self) -> None:
        for event in self.events.values():
            self.assertEqual(event.point_counts, PROTO9_POINT_COUNTS)
            self.assertFalse(event.raw_PROTO9_admission_passed)
            self.assertFalse(event.raw_spatial_spectral_admission.admission_passed)
            self.assertTrue(event.constraint_admission.admission_passed)
            self.assertTrue(event.spatial_spectral_admission.admission_passed)
            self.assertTrue(event.spatial_spectral_admission.saturation_used)
            self.assertTrue(event.profile_convergence.every_field_contracted)
            self.assertTrue(event.admission_passed)
            self.assertEqual(
                event.raw_spatial_spectral_admission.field_power_tail_ratios,
                event.spatial_spectral_admission.field_power_tail_ratios,
            )
            self.assertEqual(
                event.raw_spatial_spectral_admission.derivative_power_tail_ratios,
                event.spatial_spectral_admission.derivative_power_tail_ratios,
            )

    def test_old_or_reordered_ladder_fails_before_interpretation(self) -> None:
        grids = tuple(
            UniformRadialGrid(0.0, 128.0, count)
            for count in (1025, 2049, 4097)
        )
        states = tuple(_flat_state(grid, 4) for grid in grids)
        with self.assertRaisesRegex(ValueError, "exact point counts"):
            proto10_gr0_common_event(
                states,
                grids,
                accepted_stage_counts=(0, 0, 0),
                method="RK4",
                coordinate_time=0.0,
            )
        with self.assertRaises(ValueError):
            proto10_gr0_common_event(
                tuple(reversed(states)),
                tuple(reversed(grids)),
                accepted_stage_counts=(0, 0, 0),
                method="RK4",
                coordinate_time=0.0,
            )

    def test_direct_finest_route_remains_primary_without_saturation(self) -> None:
        grids = tuple(
            UniformRadialGrid(0.0, 128.0, count)
            for count in PROTO9_POINT_COUNTS
        )
        states = tuple(_flat_state(grid, 4) for grid in grids)
        event = proto10_gr0_common_event(
            states,
            grids,
            accepted_stage_counts=(0, 0, 0),
            method="RK4",
            coordinate_time=0.0,
        )
        self.assertTrue(event.raw_spatial_spectral_admission.admission_passed)
        self.assertTrue(event.spatial_spectral_admission.admission_passed)
        # This synthetic exact-spectrum control is not asserted to satisfy the
        # independent physical-constraint initial-data gate.
        self.assertFalse(event.spatial_spectral_admission.saturation_used)
        self.assertTrue(
            all(
                value == "directly_resolved"
                for value in (
                    *event.spatial_spectral_admission.finest_pair_field_classification.values(),
                    *event.spatial_spectral_admission.finest_pair_derivative_classification.values(),
                )
            )
        )

    def test_coarse_direct_failure_cannot_use_saturation(self) -> None:
        event = self.events[("5/2", "RK4")]
        budgets = [dict(record) for record in event.spectral_budgets]
        name = next(
            field
            for field in PROTO4_SPECTRAL_FIELD_ORDER
            if budgets[0][field].top_band_field_power_fraction > 0.0
        )
        coarse = budgets[0][name].top_band_field_power_fraction
        self.assertGreater(coarse, 0.0)
        budgets[1][name] = replace(
            budgets[1][name],
            top_band_field_power_fraction=coarse / 2.0,
        )
        attacked = resolved_or_saturated_admission(
            PROTO9_POINT_COUNTS,
            budgets,
            event.spectral_tail_sensitivities,
            event.profile_convergence,
        )
        self.assertFalse(attacked.coarse_to_medium_direct_passed_by_field[name])
        self.assertFalse(attacked.admission_passed)

    def test_missing_sensitivity_or_profile_guard_vetoes_saturation(self) -> None:
        event = self.events[("5/2", "RK4")]
        guarded = event.spatial_spectral_admission
        name = next(
            field
            for field in PROTO4_SPECTRAL_FIELD_ORDER
            if (
                guarded.finest_pair_field_classification[field]
                == "diagnostically_saturated"
                or guarded.finest_pair_derivative_classification[field]
                == "diagnostically_saturated"
            )
        )
        sensitivities = [dict(record) for record in event.spectral_tail_sensitivities]
        sensitivities[2][name] = replace(
            sensitivities[2][name], diagnostically_saturated=False
        )
        attacked_sensitivity = resolved_or_saturated_admission(
            PROTO9_POINT_COUNTS,
            event.spectral_budgets,
            sensitivities,
            event.profile_convergence,
        )
        self.assertFalse(attacked_sensitivity.admission_passed)

        complete = dict(
            event.profile_convergence.complete_profile_contraction_by_field
        )
        complete[name] = False
        convergence = replace(
            event.profile_convergence,
            complete_profile_contraction_by_field=complete,
            every_field_contracted=False,
        )
        attacked_profile = resolved_or_saturated_admission(
            PROTO9_POINT_COUNTS,
            event.spectral_budgets,
            event.spectral_tail_sensitivities,
            convergence,
        )
        self.assertFalse(attacked_profile.admission_passed)


if __name__ == "__main__":
    unittest.main()
