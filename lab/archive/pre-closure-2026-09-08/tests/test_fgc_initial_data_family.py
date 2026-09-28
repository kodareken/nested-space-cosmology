from __future__ import annotations

from fractions import Fraction
from itertools import product
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    InitialDataParameters,
    InitialDataSolveStop,
    compact_family_fields,
    exact_constraint_and_gauge_crosscheck,
    gauge_compatible_metric_time_derivatives,
    polar_areal_constraint_residuals,
    polar_areal_constraint_rhs,
    solve_initial_data,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
    solve_gr0_initial_slice,
)


Q = Fraction
PHI_SEED = 1.0 / 131072.0


class FGCInitialDataFamilyTests(unittest.TestCase):
    def test_exact_kernel_and_gauge_match_complete_unredefined_evaluator(self) -> None:
        fixtures = (
            dict(
                radius=Q(5, 2),
                radial_metric=Q(4, 3),
                angular_extrinsic_curvature=Q(1, 17),
                radial_metric_derivative=-Q(1, 19),
                angular_extrinsic_curvature_derivative=Q(1, 23),
                phi=Q(1, 29),
                phi_r=-Q(1, 31),
                phi_rr=Q(1, 37),
                phi_pi=Q(1, 41),
                phi_pi_r=-Q(1, 43),
                chi=-Q(1, 47),
                chi_r=Q(1, 53),
                chi_rr=-Q(1, 59),
                chi_pi=-Q(1, 61),
                chi_pi_r=Q(1, 67),
            ),
            dict(
                radius=Q(11, 3),
                radial_metric=Q(7, 6),
                angular_extrinsic_curvature=-Q(2, 19),
                radial_metric_derivative=Q(3, 29),
                angular_extrinsic_curvature_derivative=-Q(1, 31),
                phi=-Q(2, 37),
                phi_r=Q(1, 41),
                phi_rr=-Q(1, 43),
                phi_pi=-Q(1, 47),
                phi_pi_r=Q(1, 53),
                chi=Q(3, 59),
                chi_r=-Q(2, 61),
                chi_rr=Q(1, 67),
                chi_pi=Q(1, 71),
                chi_pi_r=-Q(1, 73),
            ),
        )
        for fixture in fixtures:
            certificate = exact_constraint_and_gauge_crosscheck(**fixture)
            self.assertTrue(certificate["exact_constraint_pair_equality"])
            self.assertTrue(certificate["exact_metric_defined_gauge_zero"])
            self.assertTrue(certificate["source_is_unredefined_ACT1_VAR1"])
            self.assertNotEqual(certificate["constraint_jacobian_diagonal"][0], 0)
            self.assertNotEqual(certificate["constraint_jacobian_diagonal"][1], 0)

    def test_solved_rhs_is_fail_closed_and_zeroes_both_constraints(self) -> None:
        parameters = InitialDataParameters(2.0, 2.0, PHI_SEED)
        derivative, coefficients, diagnostics = polar_areal_constraint_rhs(
            12.0, (1.03, 0.002), parameters
        )
        residual = polar_areal_constraint_residuals(
            coefficients,
            radial_metric_derivative=derivative[0],
            angular_extrinsic_curvature_derivative=derivative[1],
        )
        self.assertLess(max(map(abs, residual)), 1.0e-14)
        self.assertGreater(abs(diagnostics["jacobian_determinant"]), 1.0)
        self.assertLess(diagnostics["jacobian_condition_infinity"], 100.0)
        with self.assertRaisesRegex(InitialDataSolveStop, "Jacobian") as stopped:
            polar_areal_constraint_rhs(
                12.0,
                (1.03, 0.002),
                parameters,
                minimum_abs_jacobian_entry=10.0,
            )
        self.assertEqual(stopped.exception.reason, "singular_constraint_jacobian")
        self.assertEqual(stopped.exception.outcome_label, "stopped_branch_loss")

    def test_regular_buffers_profiles_and_gauge_rates(self) -> None:
        parameters = InitialDataParameters(2.0, 2.0, PHI_SEED)
        for radius in (0.0, 1.0, parameters.support_minimum):
            fields = compact_family_fields(radius, parameters)
            self.assertTrue(all(value == 0.0 for value in fields.values()))
        fields = compact_family_fields(parameters.center, parameters)
        self.assertNotEqual(fields["chi"], 0.0)
        self.assertNotEqual(fields["phi"], 0.0)
        self.assertEqual(fields["phi_pi"], 0.0)
        h_tt_dt, h_tr_dt = gauge_compatible_metric_time_derivatives(
            12.0, 1.03, -0.002
        )
        self.assertEqual(h_tt_dt, 0.0)
        self.assertNotEqual(h_tr_dt, 0.0)

    def test_central_solution_converges_and_matches_GR0_control(self) -> None:
        parameters = InitialDataParameters(2.0, 2.0, PHI_SEED)
        rk = [
            solve_initial_data(parameters, step_count=count, method="RK4")
            for count in (64, 128, 256)
        ]
        ssp = [
            solve_initial_data(parameters, step_count=count, method="SSPRK3")
            for count in (64, 128, 256)
        ]
        rk_ratio = abs(rk[0].outer_mass - rk[1].outer_mass) / abs(
            rk[1].outer_mass - rk[2].outer_mass
        )
        ssp_ratio = abs(ssp[0].outer_mass - ssp[1].outer_mass) / abs(
            ssp[1].outer_mass - ssp[2].outer_mass
        )
        self.assertGreater(rk_ratio, 2.0**3)
        self.assertGreater(ssp_ratio, 2.0**2.5)
        fine_rk = solve_initial_data(parameters, step_count=2048, method="RK4")
        fine_ssp = solve_initial_data(parameters, step_count=2048, method="SSPRK3")
        gr0 = solve_gr0_initial_slice(
            PulseParameters(chi_amplitude=2.0), step_count=4096, method="RK4"
        )
        self.assertLess(abs(fine_rk.outer_mass - fine_ssp.outer_mass), 1.0e-8)
        self.assertLess(abs(fine_rk.outer_mass - gr0.outer_mass), 1.0e-7)
        self.assertTrue(fine_rk.regular_center)
        self.assertTrue(fine_rk.exact_inner_vacuum_buffer)
        self.assertTrue(fine_rk.exact_outer_vacuum_buffer)
        self.assertTrue(fine_rk.finite_mass)
        self.assertTrue(fine_rk.no_initial_trapped_sphere)
        self.assertTrue(fine_rk.compactness_inside_protocol_window)
        self.assertGreater(fine_rk.minimum_vacuum_metric_denominator, 0.5)
        self.assertGreater(fine_rk.minimum_effective_planck_coefficient, 3.9)
        self.assertLess(fine_rk.maximum_constraint_residual_infinity, 1.0e-14)

    def test_declared_parameter_box_vertices_all_construct(self) -> None:
        for amplitude, width, seed_factor in product(
            (2.0, 2.25), (1.875, 2.0), (0.5, 2.0)
        ):
            parameters = InitialDataParameters(
                amplitude, width, PHI_SEED * seed_factor
            )
            primary = solve_initial_data(parameters, step_count=512, method="RK4")
            comparator = solve_initial_data(
                parameters, step_count=512, method="SSPRK3"
            )
            self.assertTrue(primary.compactness_inside_protocol_window)
            self.assertTrue(primary.no_initial_trapped_sphere)
            self.assertGreater(primary.minimum_abs_constraint_jacobian_determinant, 3.0)
            self.assertLess(primary.maximum_constraint_jacobian_condition_infinity, 20.0)
            self.assertLess(abs(primary.outer_mass - comparator.outer_mass), 1.0e-6)

    def test_invalid_family_inputs_and_methods_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "centre buffer"):
            InitialDataParameters(2.0, 13.0, PHI_SEED)
        parameters = InitialDataParameters(2.0, 2.0, PHI_SEED)
        with self.assertRaisesRegex(ValueError, "unknown ID1 integration method"):
            solve_initial_data(parameters, step_count=16, method="Euler")
        with self.assertRaisesRegex(ValueError, "positive integer"):
            solve_initial_data(parameters, step_count=0)


if __name__ == "__main__":
    unittest.main()
