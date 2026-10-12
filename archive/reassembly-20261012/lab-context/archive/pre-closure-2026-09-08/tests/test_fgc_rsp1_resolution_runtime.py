from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    exact_spherical_minkowski_reference,
)
from recursive_horizons.fgc.evolution.rsp1_resolution_runtime import (  # noqa: E402
    RSP1GR0EvolutionOperator,
    diagnose_gr0_reference_balanced_accelerations_chunked,
)
from recursive_horizons.fgc.evolution.src3_reference_balanced_source import (  # noqa: E402
    diagnose_gr0_reference_balanced_accelerations,
)


class FGCRSP1ResolutionRuntimeTests(unittest.TestCase):
    def test_chunking_is_bitwise_identical_to_one_direct_source_call(self) -> None:
        radii = np.arange(1, 18, dtype=np.float64) / 16.0
        points = radii.size
        u = np.zeros((points, 6), dtype=np.float64)
        p = np.zeros_like(u)
        q = np.zeros_like(u)
        p_r = np.zeros_like(u)
        q_r = np.zeros_like(u)
        u[:, 0] = 1.0
        u[:, 2] = 1.0
        u[:, 3] = radii
        q[:, 3] = 1.0
        # A dyadic non-reference perturbation makes this more than a zero
        # fixture while retaining a well-conditioned affine source.
        u[5:9, 4] = 1.0 / 4096.0
        q[5, 4] = 1.0 / 8192.0

        direct = diagnose_gr0_reference_balanced_accelerations(
            u, p, q, p_r, q_r, radii
        )
        chunked = diagnose_gr0_reference_balanced_accelerations_chunked(
            u, p, q, p_r, q_r, radii, point_batch_size=5
        )
        self.assertTrue(np.array_equal(chunked.accelerations, direct.accelerations))
        self.assertTrue(np.array_equal(chunked.residuals, direct.residuals))
        self.assertEqual(chunked.residual_infinity, direct.residual_infinity)
        self.assertTrue(chunked.raw_gate_passed)
        self.assertEqual(chunked.point_batch_count, 4)

    def test_exact_reference_is_bitwise_stationary(self) -> None:
        grid = UniformRadialGrid(0.0, 8.0, 257)
        state = exact_spherical_minkowski_reference(grid)
        before = (state.u.copy(), state.p.copy(), state.q.copy())
        operator = RSP1GR0EvolutionOperator(
            grid,
            spatial_order=4,
            ko_dissipation=1.0 / 64.0,
            raw_tolerance=1.0e-12,
            kinetic_condition_maximum=1.0e10,
            point_batch_size=64,
        )
        rhs = operator(0.0, state)
        self.assertTrue(np.array_equal(rhs.du, np.zeros_like(rhs.du)))
        self.assertTrue(np.array_equal(rhs.dp, np.zeros_like(rhs.dp)))
        self.assertTrue(np.array_equal(rhs.dq, np.zeros_like(rhs.dq)))
        self.assertEqual(rhs.diagnostics["source_residual_infinity"], 0.0)
        self.assertTrue(
            rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]
        )
        for observed, expected in zip((state.u, state.p, state.q), before, strict=True):
            self.assertTrue(np.array_equal(observed, expected))

    def test_invalid_batch_and_noncentre_grid_fail_closed(self) -> None:
        radii = np.asarray([1.0], dtype=np.float64)
        values = np.zeros((1, 6), dtype=np.float64)
        with self.assertRaisesRegex(ValueError, "point_batch_size"):
            diagnose_gr0_reference_balanced_accelerations_chunked(
                values,
                values,
                values,
                values,
                values,
                radii,
                point_batch_size=0,
            )
        with self.assertRaisesRegex(ValueError, "regular-centre"):
            RSP1GR0EvolutionOperator(
                UniformRadialGrid(1.0, 2.0, 17),
                spatial_order=4,
                ko_dissipation=0.0,
                raw_tolerance=1.0e-12,
                kinetic_condition_maximum=1.0e10,
            )


if __name__ == "__main__":
    unittest.main()
