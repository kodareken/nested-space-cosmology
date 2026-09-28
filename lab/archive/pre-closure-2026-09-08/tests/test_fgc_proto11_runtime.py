from __future__ import annotations

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
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    Proto11GR0EvolutionOperator,
    exact_spherical_minkowski_reference,
    project_gr0_reference_balanced_state,
    proto11_gr0_common_event,
    reference_balanced_reduction_constraint,
    reference_balanced_semidiscrete_constraint_snapshot,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = Fraction


def _initial(method: str, count: int, amplitude: str = "5/2"):
    order = 4 if method == "RK4" else 2
    return construct_gr0_grid_initial_data(
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


def _operator(grid: UniformRadialGrid, order: int):
    return Proto11GR0EvolutionOperator(
        grid,
        spatial_order=order,
        ko_dissipation=1.0 / 64.0,
        raw_tolerance=1.0e-12,
        kinetic_condition_maximum=1.0e10,
    )


class FGCProto11RuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.states = {}
        for method in ("RK4", "SSPRK3"):
            order = 4 if method == "RK4" else 2
            records = []
            for count in PROTO9_POINT_COUNTS:
                initial = _initial(method, count)
                state = project_gr0_reference_balanced_state(
                    initial, spatial_order=order
                )
                records.append((initial, state))
            cls.states[method] = records

    def test_projection_preserves_physical_state_and_uses_one_reference_map(self) -> None:
        for method, records in self.states.items():
            order = 4 if method == "RK4" else 2
            for initial, state in records:
                legacy = project_gr0_semidiscrete_state(
                    initial, spatial_order=order
                )
                reference = exact_spherical_minkowski_reference(initial.grid)
                derivative = SBPFirstDerivative(initial.grid, order)
                expected_q = derivative.differentiate(
                    state.u - reference.u,
                    center_parities=ADM_CENTER_PARITIES,
                ) + reference.q
                self.assertTrue(np.array_equal(state.u, legacy.u))
                self.assertTrue(np.array_equal(state.p, legacy.p))
                self.assertTrue(np.array_equal(state.q, expected_q))
                self.assertEqual(
                    np.max(
                        np.abs(
                            reference_balanced_reduction_constraint(
                                state, derivative
                            )
                        )
                    ),
                    0.0,
                )

    def test_exact_reference_is_bitwise_stationary_but_one_bit_is_not_special(self) -> None:
        for method in ("RK4", "SSPRK3"):
            order = 4 if method == "RK4" else 2
            for count in PROTO9_POINT_COUNTS:
                grid = UniformRadialGrid(0.0, 128.0, count)
                state = exact_spherical_minkowski_reference(grid)
                before = array_content_sha256(state.u, state.p, state.q)
                rhs = _operator(grid, order)(0.0, state)
                self.assertEqual(before, array_content_sha256(state.u, state.p, state.q))
                self.assertTrue(np.array_equal(rhs.du, np.zeros_like(rhs.du)))
                self.assertTrue(np.array_equal(rhs.dp, np.zeros_like(rhs.dp)))
                self.assertTrue(np.array_equal(rhs.dq, np.zeros_like(rhs.dq)))
                self.assertTrue(rhs.diagnostics["source_raw_gate_passed"])
                self.assertTrue(
                    rhs.diagnostics[
                        "PROTO11_exact_reference_equilibrium_applied"
                    ]
                )

        grid = UniformRadialGrid(0.0, 128.0, PROTO9_POINT_COUNTS[0])
        reference = exact_spherical_minkowski_reference(grid)
        q = reference.q.copy()
        q[1, 3] = np.nextafter(q[1, 3], np.inf)
        rhs = _operator(grid, 4)(0.0, EvolutionState(reference.u, reference.p, q))
        self.assertFalse(
            rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]
        )
        self.assertLess(rhs.diagnostics["source_residual_infinity"], 1.0e-12)

    def test_source_and_common_event_do_not_reproject_q(self) -> None:
        records = self.states["RK4"]
        states = [item[1] for item in records]
        grids = [item[0].grid for item in records]
        hashes = [array_content_sha256(x.u, x.p, x.q) for x in states]
        event = proto11_gr0_common_event(
            states,
            grids,
            accepted_stage_counts=(0, 0, 0),
            method="RK4",
            coordinate_time=0.0,
        )
        self.assertTrue(event.reference_balanced_constraint_admission_passed)
        self.assertTrue(event.legacy_PROTO10_admission_passed)
        self.assertTrue(event.admission_passed)
        self.assertFalse(event.interior_q_reprojected)
        self.assertEqual(
            hashes, [array_content_sha256(x.u, x.p, x.q) for x in states]
        )

        state = states[0]
        snapshot = reference_balanced_semidiscrete_constraint_snapshot(
            state,
            grids[0],
            coordinate_time=0.0,
            diagnostic_spatial_order=4,
        )
        expected = reference_balanced_reduction_constraint(
            state, SBPFirstDerivative(grids[0], 4)
        )
        self.assertTrue(np.array_equal(snapshot.residuals[:, 4:], expected))

    def test_generic_q_defect_remains_visible_to_both_maps(self) -> None:
        grid = UniformRadialGrid(0.0, 128.0, PROTO9_POINT_COUNTS[0])
        reference = exact_spherical_minkowski_reference(grid)
        q = reference.q.copy()
        q[grid.point_count // 4, 5] += 1.0e-6
        state = EvolutionState(reference.u, reference.p, q)
        derivative = SBPFirstDerivative(grid, 4)
        _, q_r, _ = reference_balanced_spatial_derivatives(state, derivative)
        reduction = reference_balanced_reduction_constraint(state, derivative)
        self.assertGreater(np.max(np.abs(q_r)), 0.0)
        self.assertGreater(np.max(np.abs(reduction)), 0.0)


if __name__ == "__main__":
    unittest.main()
