from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.covariant_principal_health import (  # noqa: E402
    background_from_spherical_state,
    point_spectrum_diagnostics,
)
from recursive_horizons.fgc.evolution.initial_state_bridge import (  # noqa: E402
    solve_initial_second_jet,
    unaccelerated_initial_second_jet,
)
from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    InitialDataParameters,
    solve_initial_data,
)


class FGCInitialStateBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.parameters = InitialDataParameters(
            chi_amplitude=2.0,
            chi_half_width=2.0,
            phi_amplitude=1.0 / 131072.0,
        )
        cls.solution = solve_initial_data(cls.parameters, step_count=1024)
        cls.index = min(
            range(2, len(cls.solution.points) - 2),
            key=lambda item: abs(cls.solution.points[item].radius - 12.0),
        )

    def test_complete_second_jet_solves_the_same_ref1_rows(self) -> None:
        result = solve_initial_second_jet(self.solution, self.index)
        self.assertTrue(result.acceleration_solve.converged)
        self.assertTrue(result.acceleration_solve.branch_continuity_preserved)
        self.assertLess(result.acceleration_solve.residual_infinity, 1.0e-12)
        self.assertLess(
            result.acceleration_solve.jacobian_condition_infinity, 1.0e10
        )
        self.assertAlmostEqual(result.coordinate_radius, 12.0)

        background = background_from_spherical_state(result.accelerated_state)
        diagnostics = point_spectrum_diagnostics(
            background, (2.0**-0.5, 2.0**-0.5, 0.0)
        )
        self.assertLess(diagnostics["maximum_abs_imaginary_part"], 1.0e-10)
        self.assertGreater(
            diagnostics["kinetic"]["kinetic_smallest_singular_value"], 0.0
        )

    def test_five_point_spatial_reconstruction_converges(self) -> None:
        coarse = solve_initial_data(self.parameters, step_count=512)
        fine = self.solution
        finer = solve_initial_data(self.parameters, step_count=2048)

        def reconstructed(solution):
            index = min(
                range(2, len(solution.points) - 2),
                key=lambda item: abs(solution.points[item].radius - 12.0),
            )
            _state, diagnostic = unaccelerated_initial_second_jet(solution, index)
            return diagnostic

        first = reconstructed(coarse)
        second = reconstructed(fine)
        third = reconstructed(finer)
        lambda_ratio = abs(
            first.lambda_rr_five_point - second.lambda_rr_five_point
        ) / abs(second.lambda_rr_five_point - third.lambda_rr_five_point)
        k_ratio = abs(first.k_rr_five_point - second.k_rr_five_point) / abs(
            second.k_rr_five_point - third.k_rr_five_point
        )
        self.assertGreater(lambda_ratio, 12.0)
        self.assertGreater(k_ratio, 12.0)

    def test_bridge_rejects_support_edges_without_hidden_extrapolation(self) -> None:
        with self.assertRaisesRegex(ValueError, "two support neighbors"):
            unaccelerated_initial_second_jet(self.solution, 0)


if __name__ == "__main__":
    unittest.main()
