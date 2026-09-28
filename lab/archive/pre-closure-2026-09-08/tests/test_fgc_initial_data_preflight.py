from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
    compact_bump_with_derivatives,
    exact_gr0_constraint_crosscheck,
    gr0_constraint_residuals,
    gr0_constraint_rhs,
    observed_convergence_order,
    proto1_compactness_upper_bound,
    pulse_fields,
    solve_gr0_initial_slice,
)


EXACT_FIXTURES = (
    {
        "radius": Q(5, 2),
        "radial_metric": Q(4, 3),
        "angular_extrinsic_curvature": Q(1, 17),
        "radial_metric_derivative": -Q(1, 19),
        "angular_extrinsic_curvature_derivative": Q(1, 23),
        "phi": Q(1, 29),
        "phi_r": -Q(1, 31),
        "phi_pi": Q(1, 37),
        "chi": -Q(1, 41),
        "chi_r": Q(1, 43),
        "chi_pi": -Q(1, 47),
    },
    {
        "radius": Q(11, 3),
        "radial_metric": Q(7, 6),
        "angular_extrinsic_curvature": -Q(2, 19),
        "radial_metric_derivative": Q(3, 29),
        "angular_extrinsic_curvature_derivative": -Q(1, 31),
        "phi": -Q(2, 37),
        "phi_r": Q(1, 41),
        "phi_pi": -Q(1, 43),
        "chi": Q(3, 47),
        "chi_r": -Q(2, 53),
        "chi_pi": Q(1, 59),
    },
)


class FGCInitialDataPreflightTests(unittest.TestCase):
    def test_compact_bump_and_frozen_profiles(self) -> None:
        self.assertEqual(compact_bump_with_derivatives(10, center=12, half_width=2), (0, 0, 0))
        value, first, second = compact_bump_with_derivatives(12, center=12, half_width=2)
        self.assertEqual(value, 1.0)
        self.assertEqual(first, 0.0)
        self.assertEqual(second, -0.5)

        parameters = PulseParameters(0.125)
        fields = pulse_fields(12, parameters)
        self.assertEqual(fields["chi"], 0.125 / 12)
        self.assertEqual(fields["chi_r"], -0.125 / 12**2)
        self.assertEqual(fields["chi_pi"], 0.0)
        self.assertEqual(fields["phi"], 1.0 / 131072.0)
        self.assertEqual(fields["phi_pi"], 0.0)

    def test_specialized_constraints_equal_full_unredefined_evaluator(self) -> None:
        for fixture in EXACT_FIXTURES:
            result = exact_gr0_constraint_crosscheck(**fixture)
            self.assertTrue(result["exact_pair_equality"])
            self.assertTrue(result["source_is_unredefined_ACT1_VAR1"])
            self.assertEqual(result["specialized_H"], result["full_unredefined_H"])
            self.assertEqual(result["specialized_M"], result["full_unredefined_M"])

    def test_constraint_rhs_solves_specialized_residuals(self) -> None:
        parameters = PulseParameters(2.0)
        radius = 12.5
        radial_metric = 1.03
        angular_curvature = -0.002
        lam_r, k_r = gr0_constraint_rhs(
            radius,
            (radial_metric, angular_curvature),
            parameters,
        )
        fields = pulse_fields(radius, parameters)
        hamiltonian, momentum = gr0_constraint_residuals(
            radius=Q(str(radius)),
            radial_metric=Q(str(radial_metric)),
            angular_extrinsic_curvature=Q(str(angular_curvature)),
            radial_metric_derivative=Q(str(lam_r)),
            angular_extrinsic_curvature_derivative=Q(str(k_r)),
            phi=Q(str(fields["phi"])),
            phi_r=Q(str(fields["phi_r"])),
            phi_pi=0,
            chi=Q(str(fields["chi"])),
            chi_r=Q(str(fields["chi_r"])),
            chi_pi=Q(str(fields["chi_pi"])),
        )
        self.assertLess(abs(float(hamiltonian)), 2e-15)
        self.assertLess(abs(float(momentum)), 2e-15)

    def test_entire_PROTO1_box_bound_closes_strictly(self) -> None:
        bound = proto1_compactness_upper_bound()
        self.assertEqual(bound["chi_pi_absolute_bound"], Q(3, 160))
        self.assertEqual(bound["chi_r_absolute_bound"], Q(1, 50))
        self.assertEqual(bound["momentum_density_absolute_bound"], Q(3, 8000))
        self.assertLess(bound["compactness_absolute_bound"], Q(1, 64))
        self.assertLess(Q(1, 64), Q(1, 10))
        self.assertTrue(bound["inverse_metric_bootstrap_closed"])
        self.assertLess(
            bound["improved_inverse_radial_metric_squared_upper_bound"],
            bound["assumed_inverse_radial_metric_squared_upper_bound"],
        )

    def test_independent_integrators_are_finite_and_convergent(self) -> None:
        parameters = PulseParameters(0.125)
        method_records = {}
        for method in ("RK4", "SSPRK3"):
            solutions = tuple(
                solve_gr0_initial_slice(parameters, step_count=count, method=method)
                for count in (32, 64, 128)
            )
            order = observed_convergence_order(*(item.outer_mass for item in solutions))
            self.assertIsNotNone(order)
            self.assertGreater(order, 2.5)
            self.assertTrue(all(item.finite_mass for item in solutions))
            self.assertTrue(all(item.no_initial_trapped_sphere for item in solutions))
            method_records[method] = solutions[-1]
        self.assertLess(
            abs(
                method_records["RK4"].peak_compactness
                - method_records["SSPRK3"].peak_compactness
            ),
            1e-8,
        )

    def test_invalid_inputs_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "centre buffer"):
            PulseParameters(1.0, center=1.0, half_width=2.0)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            solve_gr0_initial_slice(PulseParameters(1.0), step_count=0)
        with self.assertRaisesRegex(ValueError, "unknown"):
            solve_gr0_initial_slice(PulseParameters(1.0), step_count=16, method="Euler")
        with self.assertRaisesRegex(TypeError, "exact rational"):
            exact_gr0_constraint_crosscheck(
                **{**EXACT_FIXTURES[0], "radius": 2.5},
            )


if __name__ == "__main__":
    unittest.main()
