from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_src1_nl1 import load_config  # noqa: E402
from recursive_horizons.fgc.action import ActionParameters, ModelID  # noqa: E402
from recursive_horizons.fgc.evolution.nonlinear_source import (  # noqa: E402
    solve_accelerations,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
    adm_accelerations_from_base_metric,
    adm_lower_jets_from_base_state,
    regular_center_acceleration_limit,
    solve_grid_accelerations,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    state_from_generalized_adm_pg_fixture,
)


class FGCVectorizedSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.loaded = load_config()
        cls.qift = cls.loaded["qift"]
        parameters = cls.qift.flat_fixture["action_parameters"]
        cls.action = ActionParameters(
            model_id=ModelID.FGC_QR,
            planck_mass=float(parameters["planck_mass"]),
            scalar_mass=float(parameters["scalar_mass"]),
            quartic_coupling=float(parameters["quartic_coupling"]),
            pulse_width=1.0,
            ricci_coupling=float(parameters["ricci_coupling"]),
            quadratic_gb_coupling=float(parameters["quadratic_gb_coupling"]),
        )
        cls.flat = state_from_generalized_adm_pg_fixture(cls.qift.flat_fixture)
        cls.activated = activated_compatible_state(cls.qift.flat_fixture)["state"]

    def test_flat_annulus_grid_has_zero_root(self) -> None:
        radii = np.linspace(0.5, 4.0, 12)
        u = np.zeros((radii.size, 6))
        u[:, 0] = 1.0
        u[:, 2] = 1.0
        u[:, 3] = radii
        p = np.zeros_like(u)
        q = np.zeros_like(u)
        q[:, 3] = 1.0
        result = solve_grid_accelerations(
            u, p, q, np.zeros_like(u), np.zeros_like(u), radii,
            action=self.action,
        )
        self.assertEqual(result.iterations, 0)
        self.assertLessEqual(result.residual_infinity, 2.0e-14)
        np.testing.assert_allclose(result.accelerations, 0.0, atol=2.0e-14)
        self.assertTrue(result.affine_in_accelerations_verified)

    def test_activated_vectorized_root_matches_scalar_REF1_solver(self) -> None:
        u, p, q, p_r, q_r = adm_lower_jets_from_base_state(self.activated)
        radius = float(self.qift.coordinate_radius)
        grid = solve_grid_accelerations(
            u[None, :], p[None, :], q[None, :], p_r[None, :], q_r[None, :],
            np.asarray((radius,)), action=self.action,
        )
        scalar = solve_accelerations(
            self.activated,
            branch_center=self.flat,
            config=self.loaded["solver"],
            reference=flat_spherical_annulus_reference(
                radial_domain_minimum=self.qift.radial_domain_minimum
            ),
            coordinate_radius=self.qift.coordinate_radius,
            tilde_normal_factor=self.qift.tilde_normal_factor,
            hat_normal_factor=self.qift.hat_normal_factor,
        )
        scalar_adm = adm_accelerations_from_base_metric(
            u, p, scalar.accelerations
        )
        np.testing.assert_allclose(
            grid.accelerations[0], scalar_adm, rtol=0.0, atol=1.0e-12
        )
        self.assertLessEqual(grid.residual_infinity, 1.0e-12)

    def test_center_limit_is_exact_for_even_quartics_and_zeros_odd_fields(self) -> None:
        radii = np.arange(1.0, 5.0)
        values = np.empty((4, 6))
        for field in range(6):
            values[:, field] = (field + 1.0) + (field + 2.0) * radii**2 + radii**4
        limit = regular_center_acceleration_limit(values)
        expected = np.arange(1.0, 7.0)
        expected[ADM_CENTER_PARITIES == -1] = 0.0
        np.testing.assert_allclose(limit.acceleration, expected, atol=1.0e-12)
        np.testing.assert_allclose(
            limit.three_point_acceleration, expected, atol=1.0e-12
        )
        self.assertLessEqual(limit.estimator_infinity, 1.0e-12)


if __name__ == "__main__":
    unittest.main()
