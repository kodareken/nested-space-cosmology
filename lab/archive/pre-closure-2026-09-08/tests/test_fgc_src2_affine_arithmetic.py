from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    SBPFirstDerivative,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    project_gr0_reference_balanced_state,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.src2_affine_arithmetic import (  # noqa: E402
    CapturedAffineSystem,
    assemble_forward_affine_system,
    assemble_symmetric_affine_system,
    compare_affine_arithmetic,
    complete_residual,
    exact_binary_solve,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


class FGCSRC2AffineArithmeticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        initial = construct_gr0_grid_initial_data(
            PulseParameters(
                chi_amplitude=2.5,
                center=12.0,
                half_width=2.0,
                phi_amplitude=float(Fraction(1, 131072)),
                planck_mass=2.0,
                scalar_mass=3.0,
                quartic_coupling=0.5,
            ),
            point_count=129,
            outer_radius=128.0,
            constraint_method="RK4",
            diagnostic_spatial_order=4,
        )
        state = project_gr0_reference_balanced_state(initial, spatial_order=4)
        derivative = SBPFirstDerivative(initial.grid, 4)
        p_r, q_r, _differentiated_u = reference_balanced_spatial_derivatives(
            state, derivative
        )
        cls.jet = (
            state.u[1:],
            state.p[1:],
            state.q[1:],
            p_r[1:],
            q_r[1:],
            initial.grid.coordinates[1:],
        )

    def test_forward_capture_reconstructs_the_complete_affine_rows(self) -> None:
        system = assemble_forward_affine_system(*self.jet)
        acceleration = np.zeros_like(system.constant)
        acceleration[:, 0] = 0.125
        acceleration[:, 3] = -0.25
        complete = complete_residual(*self.jet, acceleration)
        reconstructed = system.constant + np.einsum(
            "nij,nj->ni", system.jacobian, acceleration
        )
        scale = max(1.0, float(np.max(np.abs(complete), initial=0.0)))
        self.assertLessEqual(
            float(np.max(np.abs(complete - reconstructed), initial=0.0)),
            4096.0 * np.finfo(np.float64).eps * scale,
        )

    def test_stable_routes_preserve_a_nontrivial_PROTO11_control(self) -> None:
        comparison, forward, symmetric = compare_affine_arithmetic(
            *self.jet,
            raw_tolerance=1.0e-12,
            condition_number_maximum=1.0e10,
            symmetric_seed_scales=(1.0, 16.0, 256.0),
        )
        self.assertTrue(comparison.baseline.strict_raw_gate_passed)
        self.assertEqual(forward.label, "forward_zero_plus_unit")
        self.assertEqual(len(symmetric), 3)
        self.assertTrue(
            all(item.strict_raw_gate_passed for item in comparison.candidates)
        )
        self.assertLessEqual(
            comparison.best.complete_residual_infinity,
            comparison.baseline.complete_residual_infinity,
        )
        self.assertTrue(
            all(
                item.complete_residual_infinity < comparison.raw_tolerance
                for item in comparison.candidates
            )
        )

    def test_exact_binary_solver_uses_only_selected_point_systems(self) -> None:
        jacobian = np.stack(
            (
                np.diag([0.5, 2.0, 4.0, 8.0, 16.0, 32.0]),
                np.diag([64.0, 32.0, 16.0, 8.0, 4.0, 2.0]),
            )
        )
        constant = np.asarray(
            (
                (-0.5, -4.0, -12.0, -32.0, -80.0, -192.0),
                (-448.0, -256.0, -144.0, -80.0, -44.0, -24.0),
            ),
            dtype=np.float64,
        )
        system = CapturedAffineSystem("synthetic", 1.0, constant, jacobian)
        base = np.full((2, 6), -7.0, dtype=np.float64)
        solved = exact_binary_solve(system, (1,), base_solution=base)
        self.assertTrue(np.array_equal(solved[0], base[0]))
        self.assertTrue(
            np.array_equal(
                solved[1], np.asarray((7.0, 8.0, 9.0, 10.0, 11.0, 12.0))
            )
        )

    def test_non_power_two_probe_and_broadened_controls_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "power of two"):
            assemble_symmetric_affine_system(*self.jet, seed_scale=3.0)
        with self.assertRaisesRegex(ValueError, "exact_binary_point_limit"):
            compare_affine_arithmetic(
                *self.jet,
                raw_tolerance=1.0e-12,
                condition_number_maximum=1.0e10,
                exact_binary_point_limit=0,
            )
        with self.assertRaisesRegex(ValueError, "selection_ratio"):
            compare_affine_arithmetic(
                *self.jet,
                raw_tolerance=1.0e-12,
                condition_number_maximum=1.0e10,
                exact_binary_selection_ratio=2.0,
            )


if __name__ == "__main__":
    unittest.main()
