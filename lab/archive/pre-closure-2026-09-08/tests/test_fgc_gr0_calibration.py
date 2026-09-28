from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    GR0EvolutionOperator,
    Q_CENTER_PARITIES,
    construct_gr0_grid_initial_data,
    make_gr0_center_boundary_projector,
    radial_null_observables,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


class FGCGR0CalibrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.initial = construct_gr0_grid_initial_data(
            PulseParameters(chi_amplitude=2.0),
            point_count=1025,
        )

    def test_constraint_slice_is_embedded_with_exact_center_and_vacua(self) -> None:
        initial = self.initial
        self.assertEqual(initial.support_minimum_index, 80)
        self.assertEqual(initial.support_maximum_index, 112)
        self.assertTrue(initial.no_initial_trapped_sphere)
        self.assertGreaterEqual(initial.peak_compactness, 0.1)
        self.assertLessEqual(initial.peak_compactness, 0.75)
        self.assertTrue(np.all(initial.state.u[0, ADM_CENTER_PARITIES == -1] == 0.0))
        self.assertTrue(np.all(initial.state.p[0, ADM_CENTER_PARITIES == -1] == 0.0))
        self.assertTrue(np.all(initial.state.q[0, Q_CENTER_PARITIES == -1] == 0.0))
        self.assertEqual(initial.state.q[0, 3], initial.state.u[0, 2])
        self.assertTrue(np.all(initial.state.u[:80, 2] == 1.0))
        self.assertTrue(np.all(initial.state.u[:80, 4:] == 0.0))
        self.assertTrue(np.all(initial.state.u[113:, 4:] == 0.0))

    def test_initial_null_signs_and_mass_match_the_constraint_solution(self) -> None:
        observed = radial_null_observables(self.initial.state)
        trapped = (observed.theta_plus[1:] < 0.0) & (
            observed.theta_minus[1:] < 0.0
        )
        self.assertFalse(np.any(trapped))
        self.assertAlmostEqual(
            float(np.max(observed.compactness)),
            self.initial.peak_compactness,
            places=13,
        )
        self.assertAlmostEqual(
            observed.misner_sharp_mass[-1], self.initial.outer_mass, places=12
        )

    def test_REF1_right_hand_side_consumes_the_solved_slice(self) -> None:
        rhs = GR0EvolutionOperator(
            self.initial.grid,
            spatial_order=4,
        )(0.0, self.initial.state)
        np.testing.assert_array_equal(rhs.du, self.initial.state.p)
        self.assertLessEqual(rhs.diagnostics["source_residual_infinity"], 1.0e-12)
        self.assertLess(rhs.diagnostics["kinetic_condition_infinity"], 1.0e10)
        self.assertLess(rhs.diagnostics["acceleration_infinity"], 4.0)
        self.assertLessEqual(rhs.diagnostics["gauge_constraint_infinity"], 1.0e-12)
        self.assertTrue(np.all(np.isfinite(rhs.dp)))
        self.assertTrue(np.all(np.isfinite(rhs.dq)))

    def test_projector_restores_center_and_causally_excluded_outer_rows(self) -> None:
        projector = make_gr0_center_boundary_projector(
            self.initial, fixed_outer_rows=4
        )
        u = self.initial.state.u.copy()
        p = self.initial.state.p.copy()
        q = self.initial.state.q.copy()
        u[0] = 7.0
        p[0] = 8.0
        q[0] = 9.0
        u[-4:] = 10.0
        p[-4:] = 11.0
        q[-4:] = 12.0
        projected = projector(0.25, EvolutionState(u, p, q))
        self.assertTrue(np.all(projected.u[0, ADM_CENTER_PARITIES == -1] == 0.0))
        self.assertTrue(np.all(projected.p[0, ADM_CENTER_PARITIES == -1] == 0.0))
        self.assertTrue(np.all(projected.q[0, Q_CENTER_PARITIES == -1] == 0.0))
        self.assertEqual(projected.q[0, 3], projected.u[0, 2])
        np.testing.assert_array_equal(projected.u[-4:], self.initial.state.u[-4:])
        np.testing.assert_array_equal(projected.p[-4:], self.initial.state.p[-4:])
        np.testing.assert_array_equal(projected.q[-4:], self.initial.state.q[-4:])


if __name__ == "__main__":
    unittest.main()
