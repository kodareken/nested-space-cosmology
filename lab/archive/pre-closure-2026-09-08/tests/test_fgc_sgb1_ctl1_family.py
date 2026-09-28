from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from inspect import signature
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    compact_family_fields as id1_compact_family_fields,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    compact_bump_with_derivatives,
)
from recursive_horizons.fgc.sgb1_ctl1_constraints import (  # noqa: E402
    SGBLAnnularInputs,
    sgbl_constraint_coefficients,
)
from recursive_horizons.fgc.sgb1_ctl1_family import (  # noqa: E402
    CONSTRAINT_INTEGRATOR_NAMES,
    DECLARED_ALPHA_GB,
    DECLARED_CHI_AMPLITUDE,
    DECLARED_PHI_AMPLITUDE,
    MAX_RADIAL_CONSTRAINT_STEPS,
    SGBLConstraintFamilyStop,
    SGBLFamilyParameters,
    SGBLFamilyPoint,
    SGBLFamilySolution,
    integrate_sgbl_radial_constraints,
    integrate_sgbl_vacuum_constraint_refinement,
    sgbl_affine_constraint_coefficients,
    sgbl_analytic_exterior_vacuum_denominator,
    sgbl_branch_health_gate,
    sgbl_compact_family_fields,
    sgbl_constraint_rhs,
    sgbl_floating_exterior_vacuum_evaluation,
    sgbl_future_normal_chi_momentum,
    sgbl_initial_vacuum_outer_continuation,
    sgbl_polar_areal_state,
    sgbl_reconstruct_constraints_by_finite_difference,
    sgbl_scalar_source_on_vacuum_geometry,
    sgbl_tensor_constraint_pair,
)
from recursive_horizons.fgc.sgb1_ctl1_family import (  # noqa: E402
    sgbl_radial_constraint_rk4_step,
    sgbl_radial_constraint_ssprk3_step,
)
from recursive_horizons.fgc.spherical_reduction import residuals  # noqa: E402


def _nominal() -> SGBLFamilyParameters:
    return SGBLFamilyParameters()


class SGBLFamilyMatchingTests(unittest.TestCase):
    def test_compact_family_fields_take_no_lambda(self) -> None:
        self.assertNotIn("lambda", signature(id1_compact_family_fields).parameters)
        self.assertNotIn("radial_metric", signature(id1_compact_family_fields).parameters)
        self.assertNotIn("lambda", signature(sgbl_compact_family_fields).parameters)
        self.assertNotIn("radial_metric", signature(sgbl_compact_family_fields).parameters)
        self.assertEqual(tuple(signature(sgbl_compact_family_fields).parameters), ("radius", "parameters"))

    def test_pi_chi_equals_ingoing_bump_derivative_over_radius(self) -> None:
        parameters = _nominal()
        self.assertEqual(parameters.chi_amplitude, DECLARED_CHI_AMPLITUDE)
        self.assertEqual(parameters.phi_amplitude, DECLARED_PHI_AMPLITUDE)
        self.assertEqual(parameters.alpha_gb, DECLARED_ALPHA_GB)
        self.assertEqual(parameters.beta, 0.0)
        self.assertEqual(parameters.eta, 0.0)
        for radius in (0.0, 1.0, parameters.support_minimum, 11.0, 12.0, 13.0, 14.0):
            fields = sgbl_compact_family_fields(radius, parameters)
            independent = sgbl_future_normal_chi_momentum(radius, parameters)
            self.assertEqual(fields["chi_pi"], independent)
            self.assertEqual(fields["phi_pi"], 0.0)
            if radius > 0.0:
                _bump, bump_r, _bump_rr = compact_bump_with_derivatives(
                    radius, center=parameters.center, half_width=parameters.chi_half_width
                )
                self.assertEqual(independent, parameters.chi_amplitude * bump_r / radius)

    def test_empty_buffer_vanishes_and_does_not_divide_at_zero(self) -> None:
        parameters = _nominal()
        for radius in (0.0, 1.0, parameters.support_minimum):
            fields = sgbl_compact_family_fields(radius, parameters)
            self.assertTrue(all(value == 0.0 for value in fields.values()))

    def test_module_does_not_invoke_id1_solver_or_qr_action(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_family as family

        source = Path(family.__file__).read_text()
        self.assertNotIn("FGCQRActionParameters", source)
        self.assertNotIn("solve_initial_data", source)
        self.assertNotIn("InitialDataParameters", source)


class SGBLFamilyExactConstraintTests(unittest.TestCase):
    def test_exact_affine_kernel_matches_annular_coefficients(self) -> None:
        point = SGBLAnnularInputs(
            radius=Q(5, 2), radial_metric=Q(4, 3),
            angular_extrinsic_curvature=Q(1, 17),
            phi=Q(1, 29), phi_r=-Q(1, 31), phi_rr=Q(1, 37),
            phi_pi=Q(1, 41), phi_pi_r=-Q(1, 43),
            chi_r=Q(1, 53), chi_pi=-Q(1, 61),
            planck_mass=2, scalar_mass=3, quartic_coupling=Q(1, 2),
            alpha_gb=-Q(1, 4),
        )
        expected = sgbl_constraint_coefficients(point)
        generic = sgbl_affine_constraint_coefficients(
            radius=point.radius, radial_metric=point.radial_metric,
            angular_extrinsic_curvature=point.angular_extrinsic_curvature,
            phi=point.phi, phi_r=point.phi_r, phi_rr=point.phi_rr,
            phi_pi=point.phi_pi, phi_pi_r=point.phi_pi_r,
            chi_r=point.chi_r, chi_pi=point.chi_pi,
            planck_mass=point.planck_mass, scalar_mass=point.scalar_mass,
            quartic_coupling=point.quartic_coupling, alpha_gb=point.alpha_gb,
        )
        self.assertEqual(generic.hamiltonian_constant, expected.hamiltonian_constant)
        self.assertEqual(generic.hamiltonian_lambda_r, expected.hamiltonian_lambda_r)
        self.assertEqual(generic.momentum_constant, expected.momentum_constant)
        self.assertEqual(generic.momentum_k_r, expected.momentum_k_r)

    def test_flat_vacuum_tensor_constraints_vanish(self) -> None:
        state = sgbl_polar_areal_state(
            radius=Q(5, 2), radial_metric=1, angular_extrinsic_curvature=0,
            radial_metric_derivative=0, angular_extrinsic_curvature_derivative=0,
            phi=0, phi_r=0, phi_rr=0, phi_pi=0, phi_pi_r=0,
            chi=0, chi_r=0, chi_rr=0, chi_pi=0, chi_pi_r=0,
        )
        self.assertEqual(sgbl_tensor_constraint_pair(state), (Q(0), Q(0)))
        affine = sgbl_affine_constraint_coefficients(
            radius=Q(5, 2), radial_metric=1, angular_extrinsic_curvature=0,
            phi=0, phi_r=0, phi_rr=0, phi_pi=0, phi_pi_r=0,
            chi_r=0, chi_pi=0, planck_mass=2, scalar_mass=3,
            quartic_coupling=Q(1, 2), alpha_gb=-Q(1, 4),
        )
        self.assertEqual(
            affine.evaluate(radial_metric_derivative=0, angular_extrinsic_curvature_derivative=0),
            (Q(0), Q(0)),
        )

    def test_floating_kernel_matches_exact_nonzero_coefficient_route(self) -> None:
        # Exercise the separately written floating formula, not the exact
        # dispatch back into the annular kernel. All matter/momentum inputs
        # are nonzero and dyadic, so both routes start from identical values.
        point = dict(
            radius=2.5, radial_metric=1.25, angular_extrinsic_curvature=0.0625,
            phi=0.03125, phi_r=-0.0625, phi_rr=0.125,
            phi_pi=0.015625, phi_pi_r=-0.03125,
            chi_r=0.25, chi_pi=-0.125,
            planck_mass=2.0, scalar_mass=3.0, quartic_coupling=0.5,
            alpha_gb=-0.25,
        )
        floating = sgbl_affine_constraint_coefficients(**point)
        exact = sgbl_affine_constraint_coefficients(**{key: Q(value) for key, value in point.items()})
        for name in (
            "hamiltonian_constant", "hamiltonian_lambda_r",
            "momentum_constant", "momentum_k_r",
        ):
            self.assertAlmostEqual(getattr(floating, name), float(getattr(exact, name)), places=13)

    def test_outer_schwarzschild_is_constraint_vacuum_not_static_sgbl(self) -> None:
        # Polar-areal Schwarzschild: r=4, 2M=3, lambda=2, k=0, lambda_r=-3/4.
        state = sgbl_polar_areal_state(
            radius=4, radial_metric=2, angular_extrinsic_curvature=0,
            radial_metric_derivative=-Q(3, 4), angular_extrinsic_curvature_derivative=0,
            phi=0, phi_r=0, phi_rr=0, phi_pi=0, phi_pi_r=0,
            chi=0, chi_r=0, chi_rr=0, chi_pi=0, chi_pi_r=0,
            radial_metric_second_derivative=Q(39, 32),
        )
        self.assertEqual(sgbl_tensor_constraint_pair(state), (Q(0), Q(0)))
        gauss_bonnet = residuals(state)["GB"]
        phi_source = sgbl_scalar_source_on_vacuum_geometry(state)
        self.assertNotEqual(gauss_bonnet, 0)
        self.assertEqual(phi_source, -Q(1, 4) * gauss_bonnet)
        self.assertNotEqual(phi_source, 0)

    def test_chi_stress_is_preserved_in_the_affine_kernel(self) -> None:
        point = dict(
            radius=Q(5, 2), radial_metric=Q(4, 3),
            angular_extrinsic_curvature=Q(1, 17),
            phi=Q(1, 29), phi_r=-Q(1, 31), phi_rr=Q(1, 37),
            phi_pi=Q(1, 41), phi_pi_r=-Q(1, 43),
            chi_r=Q(1, 53), chi_pi=-Q(1, 61),
            planck_mass=2, scalar_mass=3, quartic_coupling=Q(1, 2),
            alpha_gb=-Q(1, 4),
        )
        with_chi = sgbl_affine_constraint_coefficients(**point)
        without = sgbl_affine_constraint_coefficients(**(point | {"chi_r": 0, "chi_pi": 0}))
        self.assertEqual(
            with_chi.hamiltonian_constant - without.hamiltonian_constant,
            -(point["chi_pi"] ** 2 + point["chi_r"] ** 2 / point["radial_metric"] ** 2) / 2,
        )
        self.assertEqual(
            with_chi.momentum_constant - without.momentum_constant,
            -point["chi_pi"] * point["chi_r"],
        )


class SGBLFamilyIntegratorTests(unittest.TestCase):
    def test_chi_and_normal_momentum_are_independent_of_solved_lambda(self) -> None:
        parameters = _nominal()
        before = {
            radius: (
                sgbl_compact_family_fields(radius, parameters)["chi"],
                sgbl_compact_family_fields(radius, parameters)["chi_pi"],
            )
            for radius in (11.0, 12.0, 13.0)
        }
        fields = sgbl_compact_family_fields(12.0, parameters)
        for radial_metric in (1.0, 1.25, 2.0):
            _derivative, _coefficients, diagnostics = sgbl_constraint_rhs(
                12.0, (radial_metric, 0.0), parameters
            )
            self.assertEqual(diagnostics["chi_pi"], fields["chi_pi"])
            self.assertEqual(diagnostics["chi"], fields["chi"])
        solution = integrate_sgbl_radial_constraints(parameters, step_count=8, method="RK4")
        after = {
            radius: (
                sgbl_compact_family_fields(radius, parameters)["chi"],
                sgbl_compact_family_fields(radius, parameters)["chi_pi"],
            )
            for radius in (11.0, 12.0, 13.0)
        }
        self.assertEqual(before, after)
        for point in solution.points:
            expected = sgbl_future_normal_chi_momentum(point.radius, parameters)
            self.assertEqual(point.chi_pi, expected)
            self.assertEqual(point.phi_pi, 0.0)
        self.assertFalse(solution.branch_owned_and_healthy)
        self.assertTrue(solution.not_a_time_evolution_method)
        self.assertFalse(solution.is_full_static_sgbl_solution)
        self.assertTrue(solution.exact_inner_minkowski_buffer)
        self.assertFalse(solution.exact_outer_gr_vacuum_continuation)
        self.assertEqual(
            solution.analytic_exterior_vacuum_formula,
            "k=J/r^3, lambda^{-2}=1-2M/r+J^2/r^4",
        )
        self.assertTrue(solution.sampled_nodes_untrapped)
        self.assertFalse(solution.continuous_no_initial_trapped_sphere)
        self.assertFalse(hasattr(SGBLFamilySolution, "no_initial_trapped_sphere"))

    def test_constraint_integrators_refine_independent_fd_not_just_rhs_zeros(self) -> None:
        parameters = _nominal()
        self.assertEqual(CONSTRAINT_INTEGRATOR_NAMES, ("RK4", "SSPRK3"))
        self.assertIsNotNone(sgbl_radial_constraint_rk4_step)
        self.assertIsNotNone(sgbl_radial_constraint_ssprk3_step)
        fd_maxima = []
        for count in (8, 16, 32):
            solution = integrate_sgbl_radial_constraints(
                parameters, step_count=count, method="RK4"
            )
            self.assertLess(solution.maximum_defining_rhs_residual_infinity, 1.0e-10)
            self.assertGreater(solution.minimum_vacuum_metric_denominator, 0.0)
            self.assertTrue(all(point.radial_metric > 0.0 for point in solution.points))
            self.assertTrue(solution.sampled_nodes_untrapped)
            self.assertFalse(solution.continuous_no_initial_trapped_sphere)
            self.assertLess(solution.peak_compactness, 1.0)
            self.assertTrue(solution.finite_mass)
            records = sgbl_reconstruct_constraints_by_finite_difference(solution)
            fd_norm = max(
                max(abs(item["finite_difference_H"]), abs(item["finite_difference_M"]))
                for item in records
            )
            defining = max(
                max(abs(item["defining_H"]), abs(item["defining_M"]))
                for item in records
            )
            self.assertLess(defining, 1.0e-10)
            self.assertGreater(fd_norm, defining)
            for item in records:
                self.assertEqual(
                    item["chi_pi"],
                    sgbl_future_normal_chi_momentum(item["radius"], parameters),
                )
            fd_maxima.append(fd_norm)
        self.assertLess(fd_maxima[1], fd_maxima[0])
        self.assertLess(fd_maxima[2], fd_maxima[1])
        ssp = integrate_sgbl_radial_constraints(parameters, step_count=16, method="SSPRK3")
        rk = integrate_sgbl_radial_constraints(parameters, step_count=16, method="RK4")
        self.assertEqual(rk.method, "RK4")
        self.assertEqual(ssp.method, "SSPRK3")
        self.assertGreater(rk.outer_mass, 0.0)
        self.assertGreater(ssp.outer_mass, 0.0)
        self.assertTrue(rk.sampled_nodes_untrapped)
        self.assertTrue(ssp.sampled_nodes_untrapped)
        self.assertFalse(rk.continuous_no_initial_trapped_sphere)
        self.assertTrue(rk.not_a_time_evolution_method)
        self.assertTrue(ssp.not_a_time_evolution_method)
        self.assertEqual(
            sgbl_analytic_exterior_vacuum_denominator(4, Q(3, 2), 0),
            Q(1, 4),
        )
        outer_lambda, outer_k, denominator = sgbl_floating_exterior_vacuum_evaluation(
            parameters.outer_radius,
            mass=rk.outer_mass,
            momentum_constant=rk.vacuum_momentum_constant,
        )
        alias = sgbl_initial_vacuum_outer_continuation(
            parameters.outer_radius,
            mass=rk.outer_mass,
            momentum_constant=rk.vacuum_momentum_constant,
        )
        self.assertEqual(alias, (outer_lambda, outer_k, denominator))
        self.assertEqual(outer_lambda, rk.floating_outer_radial_metric)
        self.assertEqual(outer_k, rk.floating_outer_angular_extrinsic_curvature)
        self.assertGreater(denominator, 0.0)
        self.assertFalse(rk.exact_outer_gr_vacuum_continuation)

    def test_typed_singular_metric_and_resource_stops(self) -> None:
        parameters = _nominal()
        floating = sgbl_affine_constraint_coefficients(
            radius=2.5, radial_metric=1.0, angular_extrinsic_curvature=0.0,
            phi=0.0, phi_r=-5.0, phi_rr=0.0, phi_pi=0.0, phi_pi_r=0.0,
            chi_r=0.0, chi_pi=0.0, planck_mass=2.0, scalar_mass=1.0,
            quartic_coupling=1.0, alpha_gb=-0.25,
        )
        self.assertEqual(floating.jacobian_diagonal, (0.0, 0.0))
        with self.assertRaises(SGBLConstraintFamilyStop) as uncertified:
            floating.derivatives()
        self.assertEqual(uncertified.exception.reason, "uncertified_floating_jacobian")
        self.assertIs(uncertified.exception.diagnostics["exact_kernel"], False)

        exact = sgbl_affine_constraint_coefficients(
            radius=Q(5, 2), radial_metric=1, angular_extrinsic_curvature=0,
            phi=0, phi_r=-5, phi_rr=0, phi_pi=0, phi_pi_r=0,
            chi_r=0, chi_pi=0, planck_mass=2, scalar_mass=1,
            quartic_coupling=1, alpha_gb=-Q(1, 4),
        )
        with self.assertRaises(SGBLConstraintFamilyStop) as singular:
            exact.derivatives()
        self.assertEqual(singular.exception.reason, "singular_constraint_jacobian")
        self.assertIs(singular.exception.diagnostics["exact_kernel"], True)

        with self.assertRaises(SGBLConstraintFamilyStop) as stopped:
            sgbl_constraint_rhs(2.5, (0.0, 0.0), parameters)
        self.assertEqual(stopped.exception.reason, "nonpositive_metric")
        with self.assertRaises(SGBLConstraintFamilyStop) as resource:
            integrate_sgbl_radial_constraints(
                parameters, step_count=MAX_RADIAL_CONSTRAINT_STEPS + 1, method="RK4"
            )
        self.assertEqual(resource.exception.reason, "resource_limit")
        with self.assertRaisesRegex(ValueError, "unknown SGB-L radial constraint integrator"):
            integrate_sgbl_radial_constraints(parameters, step_count=8, method="Euler")
        with self.assertRaisesRegex(ValueError, "declared SGB-L fixture"):
            replace(parameters, center=13.0)

    def test_result_records_refuse_forged_scope_claims_and_inconsistent_data(self) -> None:
        solution = integrate_sgbl_radial_constraints(_nominal(), step_count=8, method="RK4")
        with self.assertRaises(TypeError):
            replace(solution, branch_owned_and_healthy=True)
        with self.assertRaises(TypeError):
            replace(solution, is_full_static_sgbl_solution=True)
        with self.assertRaises(TypeError):
            replace(solution, no_initial_trapped_sphere=True)
        with self.assertRaises(TypeError):
            replace(solution, exact_outer_gr_vacuum_continuation=True)
        with self.assertRaises(ValueError):
            replace(solution, step_count=3)
        with self.assertRaises(ValueError):
            replace(solution, points=solution.points[:-1])
        point = solution.points[1]
        with self.assertRaises(ValueError):
            replace(point, compactness=0.0)
        with self.assertRaises(ValueError):
            replace(point, misner_sharp_mass=0.0)
        with self.assertRaises(ValueError):
            replace(point, radius=float("nan"))
        for name in ("chi", "chi_pi", "phi", "phi_pi", "hamiltonian_residual", "radial_metric_derivative"):
            altered = replace(solution.points[0], **{name: 1.0})
            with self.subTest(field=name), self.assertRaisesRegex(ValueError, "declared family"):
                replace(solution, points=(altered,) + solution.points[1:])
        with self.assertRaisesRegex(ValueError, "declared family"):
            replace(solution, parameters=replace(solution.parameters, chi_amplitude=4.0))
        self.assertIs(type(point), SGBLFamilyPoint)
        self.assertFalse(solution.continuous_no_initial_trapped_sphere)
        self.assertTrue(solution.sampled_nodes_untrapped)

    def test_known_vacuum_constraint_solution_refines_without_a_mass_tolerance(self) -> None:
        errors = []
        for count in (8, 16, 32):
            record = integrate_sgbl_vacuum_constraint_refinement(
                mass=1.0,
                momentum_constant=0.5,
                start_radius=8.0,
                end_radius=12.0,
                step_count=count,
                method="RK4",
            )
            self.assertEqual(record.method, "RK4")
            self.assertTrue(record.not_a_time_evolution_method)
            self.assertFalse(record.is_full_static_sgbl_solution)
            errors.append(record.radial_metric_error + record.k_error)
        self.assertLess(errors[1], errors[0])
        self.assertLess(errors[2], errors[1])
        self.assertGreater(errors[0], 0.0)

    def test_health_gate_stays_closed_without_proving_nonhyperbolicity(self) -> None:
        gate = sgbl_branch_health_gate()
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["missing_health_closes_the_gate"], True)
        self.assertIs(gate["missing_health_proves_nonhyperbolicity"], False)
        self.assertEqual(gate["interval_invertibility"], "unqualified")
        self.assertEqual(gate["source_runtime"], "unqualified")
        self.assertEqual(gate["kinetic_cone_health"], "unqualified")
        self.assertEqual(gate["evolving_center_after_matter_arrives"], "unqualified")
        self.assertEqual(gate["continuous_no_initial_trapped_sphere"], "unqualified")
        self.assertIs(gate["sampled_nodes_untrapped_is_not_continuous_no_initial_trap"], True)


if __name__ == "__main__":
    unittest.main()
