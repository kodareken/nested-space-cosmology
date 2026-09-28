from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.action import ActionParameters, ModelID  # noqa: E402
from recursive_horizons.fgc.evolution.gr0_direct_source import (  # noqa: E402
    gr0_ref1_residual_batch,
    solve_gr0_grid_accelerations,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    _batch_state_from_adm,
    batch_ref1_residual,
    solve_grid_accelerations,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    physical_constraint_projections,
)


class FGCGR0DirectSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.action = ActionParameters(
            ModelID.GR_0,
            planck_mass=2.0,
            scalar_mass=3.0,
            quartic_coupling=0.5,
            pulse_width=1.0,
        )

    @staticmethod
    def _nontrivial_jets(scale: float = 0.02) -> tuple[np.ndarray, ...]:
        rng = np.random.default_rng(20260822)
        radii = np.asarray((0.7, 1.1, 1.8))
        points = radii.size
        u = np.zeros((points, 6))
        u[:, 0] = 1.0 + scale * rng.normal(size=points)
        u[:, 1] = scale * rng.normal(size=points)
        u[:, 2] = 1.0 + scale * rng.normal(size=points)
        u[:, 3] = radii * (1.0 + scale * rng.normal(size=points))
        u[:, 4:] = scale * rng.normal(size=(points, 2))
        p = scale * rng.normal(size=(points, 6))
        q = scale * rng.normal(size=(points, 6))
        q[:, 3] += 1.0
        p_r = scale * rng.normal(size=(points, 6))
        q_r = scale * rng.normal(size=(points, 6))
        accelerations = scale * rng.normal(size=(4, points, 6))
        return u, p, q, p_r, q_r, accelerations, radii

    def test_flat_annulus_has_the_zero_root(self) -> None:
        radii = np.linspace(0.5, 4.0, 12)
        u = np.zeros((radii.size, 6))
        u[:, 0] = 1.0
        u[:, 2] = 1.0
        u[:, 3] = radii
        p = np.zeros_like(u)
        q = np.zeros_like(u)
        q[:, 3] = 1.0
        result = solve_gr0_grid_accelerations(
            u, p, q, np.zeros_like(u), np.zeros_like(u), radii
        )
        self.assertLessEqual(result.residual_infinity, 2.0e-14)
        np.testing.assert_allclose(result.accelerations, 0.0, atol=2.0e-14)

    def test_complete_rows_match_the_generic_REF1_authority(self) -> None:
        u, p, q, p_r, q_r, accelerations, radii = self._nontrivial_jets()
        expected = batch_ref1_residual(
            u,
            p,
            q,
            accelerations,
            p_r,
            q_r,
            radii,
            action=self.action,
        )
        observed = gr0_ref1_residual_batch(
            u, p, q, accelerations, p_r, q_r, radii
        ).full_residual
        np.testing.assert_allclose(observed, expected, rtol=2.0e-14, atol=2.0e-14)

    def test_unredefined_constraint_rows_match_the_exact_projection(self) -> None:
        u, p, q, p_r, q_r, accelerations, radii = self._nontrivial_jets()
        direct = gr0_ref1_residual_batch(
            u, p, q, accelerations, p_r, q_r, radii
        )
        state = _batch_state_from_adm(
            u, p, q, accelerations, p_r, q_r, action=self.action
        )
        exact = physical_constraint_projections(state)
        np.testing.assert_allclose(
            direct.hamiltonian_constraint,
            exact["H"].data,
            rtol=2.0e-14,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            direct.momentum_constraint,
            exact["M"].data,
            rtol=2.0e-14,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            direct.hamiltonian_constraint,
            np.broadcast_to(
                direct.hamiltonian_constraint[0][None, :],
                direct.hamiltonian_constraint.shape,
            ),
            rtol=2.0e-14,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            direct.momentum_constraint,
            np.broadcast_to(
                direct.momentum_constraint[0][None, :],
                direct.momentum_constraint.shape,
            ),
            rtol=2.0e-14,
            atol=2.0e-14,
        )

    def test_acceleration_root_matches_the_generic_REF1_solver(self) -> None:
        u, p, q, p_r, q_r, _accelerations, radii = self._nontrivial_jets(
            scale=0.001
        )
        direct = solve_gr0_grid_accelerations(u, p, q, p_r, q_r, radii)
        expected = solve_grid_accelerations(
            u,
            p,
            q,
            p_r,
            q_r,
            radii,
            action=self.action,
            branch_displacement_maximum=1.0,
        )
        np.testing.assert_allclose(
            direct.accelerations,
            expected.accelerations,
            rtol=2.0e-13,
            atol=2.0e-13,
        )
        self.assertLessEqual(direct.residual_infinity, 1.0e-12)


if __name__ == "__main__":
    unittest.main()
