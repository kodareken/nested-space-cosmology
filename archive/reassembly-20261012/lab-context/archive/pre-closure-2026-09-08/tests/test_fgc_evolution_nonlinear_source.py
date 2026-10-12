from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_src1_nl1 import load_config  # noqa: E402
from recursive_horizons.fgc.evolution.floating_jet import (  # noqa: E402
    FloatFirstTangent,
    FloatJet2,
    scalar_primal_and_tangent,
)
from recursive_horizons.fgc.evolution.nonlinear_source import (  # noqa: E402
    AccelerationSolveStop,
    acceleration_jacobian,
    finite_difference_acceleration_jacobian,
    parameter_vector,
    solve_accelerations,
    state_from_parameter_vector,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.modified_harmonic_modes import (  # noqa: E402
    exact_comp1_acceleration_root,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    Jet2,
    state_from_generalized_adm_pg_fixture,
)


class FGCNonlinearSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.loaded = load_config()
        cls.qift = cls.loaded["qift"]
        cls.center = state_from_generalized_adm_pg_fixture(cls.qift.flat_fixture)
        cls.compatible = activated_compatible_state(cls.qift.flat_fixture)["state"]
        cls.reference = flat_spherical_annulus_reference(
            radial_domain_minimum=cls.qift.radial_domain_minimum
        )
        cls.arguments = {
            "reference": cls.reference,
            "coordinate_radius": cls.qift.coordinate_radius,
            "tilde_normal_factor": cls.qift.tilde_normal_factor,
            "hat_normal_factor": cls.qift.hat_normal_factor,
        }

    def test_floating_tangent_and_jet_preserve_derivatives_and_validate_finiteness(self) -> None:
        value = FloatFirstTangent.seed(1.5, 0.25)
        expression = (value**3 + 2 * value - 1) / (value + 1)
        primal, tangent = scalar_primal_and_tangent(expression)
        expected_primal = (1.5**3 + 2 * 1.5 - 1) / 2.5
        expected_tangent = (
            ((3 * 1.5**2 + 2) * 2.5 - (1.5**3 + 2 * 1.5 - 1))
            / 2.5**2
            * 0.25
        )
        self.assertAlmostEqual(primal, expected_primal)
        self.assertAlmostEqual(tangent, expected_tangent)
        jet = FloatJet2(value=1.0, dt=2.0, dtt=FloatFirstTangent.seed(3.0))
        self.assertIsInstance(jet, Jet2)
        self.assertEqual(scalar_primal_and_tangent((jet * jet).dtt)[1], 2.0)
        with self.assertRaisesRegex(ValueError, "finite"):
            FloatJet2(value=float("nan"))

    def test_flat_residual_and_analytic_jacobian_match_the_exact_matrix(self) -> None:
        residual, jacobian, condition = acceleration_jacobian(
            self.center, [0.0] * 6, **self.arguments
        )
        expected = np.asarray(
            (
                (36, 0, 9, 9, 0, 0),
                (0, 72, 0, 0, 0, 0),
                (4, 0, 1, -1, 0, 0),
                (64, 0, -16, 0, 0, 0),
                (0, 0, 0, 0, -1, 0),
                (0, 0, 0, 0, 0, -1),
            ),
            dtype=np.float64,
        )
        np.testing.assert_array_equal(residual, np.zeros(6))
        np.testing.assert_array_equal(jacobian, expected)
        self.assertEqual(condition, 80.0)

    def test_comp1_solve_matches_exact_root_and_strictly_decreases_residual(self) -> None:
        exact = exact_comp1_acceleration_root(
            self.compatible,
            **self.arguments,
            qift_acceleration_half_width=self.qift.acceleration_half_width,
        )["root_acceleration"]
        solved = solve_accelerations(
            self.compatible,
            branch_center=self.center,
            config=self.loaded["solver"],
            **self.arguments,
        )
        self.assertTrue(solved.converged)
        self.assertLessEqual(solved.residual_infinity, 1.0e-12)
        self.assertLess(
            max(abs(float(expected) - observed) for expected, observed in zip(exact, solved.accelerations, strict=True)),
            1.0e-18,
        )
        self.assertTrue(
            all(
                item.residual_infinity_after < item.residual_infinity_before
                for item in solved.iteration_history
            )
        )

    def test_nonrational_first_tangent_matches_finite_difference(self) -> None:
        parameters = np.asarray(parameter_vector(self.center))
        parameters[4] += np.sqrt(2.0) * 1.0e-7
        parameters[7] -= np.sqrt(3.0) * 1.0e-7
        state = state_from_parameter_vector(self.center, parameters)
        acceleration = np.asarray([np.pi, -np.sqrt(2), 1, -1, 2, -2]) * 1.0e-6
        _, analytic, _ = acceleration_jacobian(
            state, acceleration, **self.arguments
        )
        finite_difference = finite_difference_acceleration_jacobian(
            state, acceleration, step=2.0**-20, **self.arguments
        )
        np.testing.assert_allclose(analytic, finite_difference, rtol=1.0e-6, atol=1.0e-6)

    def test_parameter_roundtrip_and_typed_stops_fail_closed(self) -> None:
        parameters = np.asarray(parameter_vector(self.center))
        rebuilt = state_from_parameter_vector(self.center, parameters)
        np.testing.assert_array_equal(parameter_vector(rebuilt), parameters)

        cases = (
            (
                "nonfinite_input",
                lambda: solve_accelerations(
                    self.center,
                    branch_center=self.center,
                    warm_start=[float("nan"), 0, 0, 0, 0, 0],
                    config=self.loaded["solver"],
                    **self.arguments,
                ),
            ),
            (
                "initial_acceleration_outside_authorized_box",
                lambda: solve_accelerations(
                    self.center,
                    branch_center=self.center,
                    warm_start=[self.loaded["solver"].acceleration_half_width, 0, 0, 0, 0, 0],
                    config=self.loaded["solver"],
                    **self.arguments,
                ),
            ),
            (
                "kinetic_condition_limit",
                lambda: solve_accelerations(
                    self.center,
                    branch_center=self.center,
                    config=replace(self.loaded["solver"], condition_number_maximum=1.0),
                    **self.arguments,
                ),
            ),
        )
        for expected, function in cases:
            with self.subTest(reason=expected):
                with self.assertRaises(AccelerationSolveStop) as caught:
                    function()
                self.assertEqual(caught.exception.reason, expected)


if __name__ == "__main__":
    unittest.main()
