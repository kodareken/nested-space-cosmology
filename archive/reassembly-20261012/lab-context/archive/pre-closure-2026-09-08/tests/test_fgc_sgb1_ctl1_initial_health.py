from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from inspect import getsource
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.exact_interval import Interval, interval  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_constraints import (  # noqa: E402
    SGBLAnnularInputs,
    sgbl_constraint_coefficients,
)
from recursive_horizons.fgc.sgb1_ctl1_family import SGBLFamilyParameters  # noqa: E402
from recursive_horizons.fgc.interval_tangent import IntervalFirstTangent  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_initial_health import (  # noqa: E402
    BUMP_X_DERIVATIVE_BOUND,
    DECLARED_BASE_CELLS,
    CERTIFICATE_CHART_KIND,
    DECLARED_C_DOMAIN,
    DECLARED_J_DOMAIN,
    DECLARED_ODE_POLICY,
    SGBLContinuousCompactnessRecord,
    SGBLEvolvingCenterMonitorContract,
    SGBLExactInitialSlice,
    SGBLInitialHealthStop,
    SGBLValidatedODEPolicy,
    sgbl_chart_versus_lambda_k_regression,
    sgbl_centered_mean_value_algebraic_compactness,
    sgbl_centered_mean_value_compactness_rhs,
    sgbl_centered_mean_value_rhs,
    sgbl_compact_bump_enclosure,
    sgbl_constraint_rhs_state_jacobian,
    sgbl_continuous_initial_compactness,
    sgbl_evolving_center_monitor_contract,
    sgbl_exact_minkowski_rhs_reduction,
    sgbl_execute_evolving_center_monitor,
    sgbl_generic_affine_coefficients,
    sgbl_initial_health_gate,
    sgbl_integrator_containment_regression,
    sgbl_interval_affine_coefficients,
    sgbl_interval_sqrt,
    sgbl_intersect_algebraic_compactness_enclosures,
    sgbl_misner_sharp_angular_momentum_interval,
    sgbl_misner_sharp_angular_momentum_radial_derivative,
    sgbl_misner_sharp_chart_roundtrip,
    sgbl_misner_sharp_chart_state,
    sgbl_misner_sharp_ck_chart_roundtrip,
    sgbl_misner_sharp_ck_chart_state,
    sgbl_misner_sharp_compactness_chain_rule_identity,
    sgbl_misner_sharp_compactness_interval,
    sgbl_misner_sharp_compactness_invariant_residual,
    sgbl_misner_sharp_compactness_radial_derivative,
    sgbl_misner_sharp_compactness_state_gradient,
    sgbl_validate_ode_cell_inventory,
    sgbl_validated_cj_constraint_ode,
    sgbl_validated_lambda_k_constraint_ode,
)


class SGBLBumpEnclosureTests(unittest.TestCase):
    def test_outside_support_the_bump_is_exactly_zero(self) -> None:
        value, first, second = sgbl_compact_bump_enclosure(
            interval(0, 1), center=Q(12), half_width=Q(2)
        )
        self.assertEqual(value, interval(0))
        self.assertEqual(first, interval(0))
        self.assertEqual(second, interval(0))

    def test_interior_bump_stays_inside_the_global_bounds(self) -> None:
        value, first, second = sgbl_compact_bump_enclosure(
            interval(Q(47, 4), Q(49, 4)), center=Q(12), half_width=Q(2)
        )
        self.assertGreaterEqual(value.lower, 0)
        self.assertLessEqual(value.upper, 1)
        self.assertLessEqual(first.abs_upper(), BUMP_X_DERIVATIVE_BOUND / 2)
        self.assertGreater(value.upper, 0)

    def test_near_boundary_interior_is_not_clipped_to_zero(self) -> None:
        spec = SGBLExactInitialSlice()
        radius = interval(
            spec.support_minimum + Q(1, 64),
            spec.support_minimum + Q(1, 32),
        )
        value, first, second = sgbl_compact_bump_enclosure(
            radius,
            center=spec.center,
            half_width=spec.half_width,
        )
        self.assertGreater(value.upper, 0)
        self.assertGreater(first.width(), 0)
        self.assertGreater(second.width(), 0)

    def test_boundary_crossing_uses_the_global_enclosure(self) -> None:
        spec = SGBLExactInitialSlice()
        radius = interval(
            spec.support_minimum - Q(1, 64),
            spec.support_minimum + Q(1, 64),
        )
        value, first, second = sgbl_compact_bump_enclosure(
            radius,
            center=spec.center,
            half_width=spec.half_width,
        )
        self.assertEqual(value, interval(0, 1))
        self.assertLess(first.lower, 0)
        self.assertGreater(first.upper, 0)
        self.assertLess(second.lower, 0)
        self.assertGreater(second.upper, 0)

    def test_wholly_interior_near_edge_is_not_clipped_to_seven_eighths(self) -> None:
        radius = interval(Q(101, 10), Q(51, 5))
        scaled = (radius - Q(12)) / Q(2)
        self.assertGreater(scaled.lower, -1)
        self.assertLess(scaled.upper, -Q(7, 8))
        value, first, _second = sgbl_compact_bump_enclosure(
            radius, center=Q(12), half_width=Q(2)
        )
        self.assertLess(value.upper, 1)
        self.assertGreater(value.upper, 0)
        self.assertNotEqual(value, interval(0, 1))
        self.assertGreater(first.upper, 0)


class SGBLValidatedODETests(unittest.TestCase):
    def test_minkowski_buffer_rhs_matches_the_exact_affine_kernel(self) -> None:
        spec = SGBLExactInitialSlice()
        reduction = sgbl_exact_minkowski_rhs_reduction(spec)
        self.assertEqual(reduction["lambda_r"], 0)
        self.assertEqual(reduction["k_r"], 0)
        self.assertEqual(reduction["compactness"], 0)
        self.assertEqual(reduction["compactness_r"], 0)
        self.assertEqual(reduction["angular_momentum"], 0)
        self.assertEqual(reduction["angular_momentum_r"], 0)
        self.assertEqual(reduction["denominator"], 1)
        point = SGBLAnnularInputs(
            radius=spec.support_minimum,
            radial_metric=1,
            angular_extrinsic_curvature=0,
            phi=0,
            phi_r=0,
            phi_rr=0,
            phi_pi=0,
            phi_pi_r=0,
            chi_r=0,
            chi_pi=0,
            planck_mass=2,
            scalar_mass=3,
            quartic_coupling=Q(1, 2),
            alpha_gb=-Q(1, 4),
        )
        exact = sgbl_constraint_coefficients(point)
        boxed = sgbl_interval_affine_coefficients(
            radius=Interval.singleton(spec.support_minimum),
            radial_metric=Interval.singleton(1),
            angular_extrinsic_curvature=Interval.singleton(0),
            spec=spec,
        )
        self.assertEqual(boxed["hamiltonian_lambda_r"].lower, exact.hamiltonian_lambda_r)
        self.assertEqual(boxed["momentum_k_r"].lower, exact.momentum_k_r)
        self.assertTrue(boxed["hamiltonian_lambda_r"].is_singleton())
        self.assertTrue(boxed["momentum_k_r"].is_singleton())
        self.assertFalse(boxed["hamiltonian_lambda_r"].contains_zero())
        self.assertFalse(boxed["momentum_k_r"].contains_zero())

    def test_nominal_chi_three_is_a_validated_ode_result_not_a_mass_proxy(self) -> None:
        record = sgbl_continuous_initial_compactness(SGBLExactInitialSlice())
        self.assertEqual(record.minkowski_buffer_compactness, 0)
        self.assertTrue(record.sampled_nodes_are_not_the_certificate)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertNotIn("GB_EFFECTIVE_DENSITY_BOUND", record.theorem)
        self.assertNotIn("1/256", record.theorem)
        self.assertIn("Picard", record.theorem)
        self.assertGreaterEqual(record.ode_cell_count, 1)
        self.assertGreaterEqual(record.rhs_evaluations, 1)
        self.assertIsNotNone(record.ode)
        self.assertTrue(record.ode.chart_coordinates)
        self.assertTrue(record.chart_coordinates)
        self.assertEqual(record.chart_kind, CERTIFICATE_CHART_KIND)
        self.assertEqual(record.ode.chart_kind, CERTIFICATE_CHART_KIND)
        self.assertIsInstance(record.ode.mean_value_strictly_tighter, bool)
        self.assertIsInstance(record.ode.compactness_mean_value_strictly_tighter, bool)
        self.assertIsInstance(record.ode.algebraic_mean_value_strictly_tighter, bool)
        self.assertIsInstance(record.correlated_strictly_tighter_than_direct, bool)
        self.assertGreaterEqual(record.compactness_rhs_evaluations, 1)
        self.assertEqual(record.policy.compactness_domain, DECLARED_C_DOMAIN)
        self.assertGreater(DECLARED_C_DOMAIN.upper, 1)
        self.assertLess(DECLARED_C_DOMAIN.lower, 0)
        if record.continuous_no_initial_trapped_sphere:
            self.assertEqual(record.classification, "continuous_no_initial_trapped_sphere")
            self.assertLess(record.compact_support_compactness_upper, 1)
            self.assertLess(record.exterior_compactness_upper, 1)
            self.assertGreater(record.compactness_margin, 0)
            self.assertIsNone(record.missing_theorem)
            self.assertEqual(record.picard_strict_inclusions, record.ode_cell_count)
            self.assertIsNotNone(record.ode)
            self.assertTrue(all(cell.picard_strict_self_map for cell in record.ode.cells))
            self.assertTrue(all(cell.residual_contains_origin for cell in record.ode.cells))
            self.assertTrue(
                all(cell.compactness_invariant_contains_zero for cell in record.ode.cells)
            )
            self.assertIsNotNone(record.compactness_invariant_residual)
            self.assertTrue(record.compactness_invariant_residual.contains_zero())
        else:
            self.assertEqual(record.classification, "interval_inconclusive")
            self.assertIsNotNone(record.inconclusive_reason)
            self.assertIsNotNone(record.missing_theorem)
            self.assertNotEqual(
                record.missing_theorem,
                "validated_sgbl_constraint_ode_enclosure_of_lambda_and_k_on_compact_support_cells",
            )
            self.assertFalse(record.covers_complete_support)
            self.assertEqual(record.coverage_left, record.slice.support_minimum)
            self.assertLessEqual(record.coverage_right, record.slice.support_maximum)
            if record.ode is not None and record.ode.cells:
                self.assertEqual(
                    record.ode.coverage_right, record.ode.cells[-1].radius.upper
                )
                self.assertGreater(
                    record.ode_cell_count,
                    1,
                    "retained left siblings must not vanish from a failed later interval",
                )
            self.assertNotIn("GB_EFFECTIVE_DENSITY_BOUND", record.theorem)
            if record.ode is not None and record.ode.cells:
                self.assertTrue(
                    all(
                        cell.compactness.subset_of(cell.direct_compactness)
                        and cell.compactness.subset_of(cell.centered_compactness)
                        for cell in record.ode.cells
                    )
                )
                self.assertLessEqual(
                    record.compact_support_compactness_upper,
                    record.direct_compact_support_compactness_upper,
                )

    def test_small_amplitude_still_uses_the_ode_enclosure(self) -> None:
        record = sgbl_continuous_initial_compactness(
            SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        )
        self.assertEqual(record.classification, "continuous_no_initial_trapped_sphere")
        self.assertTrue(record.continuous_no_initial_trapped_sphere)
        self.assertEqual(record.minkowski_buffer_compactness, 0)
        self.assertLess(record.compactness_upper_bound, 1)
        self.assertGreater(record.compactness_margin, 0)
        self.assertEqual(record.ode_cell_count, DECLARED_BASE_CELLS)
        self.assertTrue(record.covers_complete_support)
        self.assertEqual(record.coverage_left, record.slice.support_minimum)
        self.assertEqual(record.coverage_right, record.slice.support_maximum)
        self.assertTrue(record.sampled_nodes_are_not_the_certificate)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertTrue(record.chart_coordinates)
        self.assertEqual(record.chart_kind, CERTIFICATE_CHART_KIND)
        self.assertTrue(all(cell.chart_coordinates for cell in record.ode.cells))
        self.assertTrue(all(cell.chart_kind == CERTIFICATE_CHART_KIND for cell in record.ode.cells))
        self.assertGreater(record.denominator_margin, 0)
        self.assertTrue(all(cell.compactness_invariant_contains_zero for cell in record.ode.cells))
        self.assertTrue(record.compactness_invariant_residual.contains_zero())
        self.assertTrue(
            all(cell.compactness.subset_of(cell.direct_compactness) for cell in record.ode.cells)
        )
        self.assertLessEqual(
            record.compact_support_compactness_upper,
            record.direct_compact_support_compactness_upper,
        )
        comparison = sgbl_chart_versus_lambda_k_regression(record)
        self.assertTrue(comparison["not_the_continuum_certificate"])
        self.assertEqual(comparison["chart_kind"], CERTIFICATE_CHART_KIND)
        self.assertEqual(comparison["chart_cells"], DECLARED_BASE_CELLS)
        self.assertIsNone(comparison["lambda_k_obstruction"])
        self.assertTrue(comparison["lambda_k_consistent"])
        self.assertTrue(comparison["all_overlapping_cells_consistent"])
        self.assertIn("cj_obstruction", comparison)
        self.assertGreaterEqual(comparison["cj_cells"], 1)
        if comparison["cj_overlapping_radial_cells"]:
            self.assertTrue(comparison["cj_consistent"])

    def test_family_float_parameters_round_trip_to_the_exact_slice(self) -> None:
        spec = SGBLExactInitialSlice.from_family_parameters(SGBLFamilyParameters())
        self.assertEqual(spec.chi_amplitude, 3)
        self.assertEqual(spec.phi_amplitude, Q(1, 131072))
        self.assertEqual(spec.center, 12)
        self.assertEqual(spec.alpha_gb, -Q(1, 4))


class SGBLValidatedODERegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = sgbl_continuous_initial_compactness(SGBLExactInitialSlice())

    def test_rk4_and_ssprk3_nodes_are_regression_containment_only(self) -> None:
        if self.record.ode is None or not self.record.ode.cells:
            self.skipTest("nominal ODE produced no graph cells")
        for method in ("RK4", "SSPRK3"):
            with self.subTest(method=method):
                report = sgbl_integrator_containment_regression(
                    self.record, method=method, step_count=16
                )
                self.assertTrue(report["not_the_continuum_certificate"])
                self.assertEqual(report["nodes"], 17)
                self.assertGreater(report["eligible_nodes"], 0)
                self.assertTrue(report["all_eligible_nodes_contained"])
                self.assertEqual(
                    report["covers_complete_support"],
                    self.record.continuous_no_initial_trapped_sphere,
                )


class SGBLValidatedODEAttackTests(unittest.TestCase):
    def test_lambda_domain_containing_zero_is_a_domain_stop(self) -> None:
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            SGBLValidatedODEPolicy(lambda_domain=interval(Q(-1, 8), Q(2)))
        self.assertEqual(stopped.exception.reason, "domain_error")

    def test_rhs_budget_is_a_resource_stop_not_a_fitted_floor(self) -> None:
        policy = SGBLValidatedODEPolicy(max_rhs_evaluations=1)
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_continuous_initial_compactness(SGBLExactInitialSlice(), policy=policy)
        self.assertEqual(stopped.exception.reason, "resource_limit")

    def test_declared_one_cell_without_bisection_does_not_refit_the_tree(self) -> None:
        policy = SGBLValidatedODEPolicy(base_cells=1, max_bisection_depth=0, max_cells=1)
        record = sgbl_continuous_initial_compactness(SGBLExactInitialSlice(), policy=policy)
        self.assertEqual(record.policy.base_cells, 1)
        self.assertEqual(record.policy.max_bisection_depth, 0)
        if not record.continuous_no_initial_trapped_sphere:
            self.assertEqual(record.classification, "interval_inconclusive")
            self.assertIn(
                record.inconclusive_reason,
                {
                    "picard_strict_self_map_failed",
                    "jacobian_diagonal_contains_zero",
                    "physical_domain_escape",
                    "compactness_enclosure_not_below_one",
                    "residual_enclosure_misses_origin",
                    "exterior_compactness_not_below_one",
                    "zero_in_interval_reciprocal",
                    "mean_value_inconsistent",
                    "compactness_mean_value_inconsistent",
                    "algebraic_compactness_mean_value_inconsistent",
                    "compactness_invariant_misses_origin",
                    "chart_denominator_not_strictly_positive",
                    "sqrt_branch_not_positive",
                    "chart_roundtrip_misses_origin",
                },
            )

    def test_tampered_pass_with_c_not_below_one_is_refused(self) -> None:
        record = sgbl_continuous_initial_compactness(
            SGBLExactInitialSlice(),
            policy=SGBLValidatedODEPolicy(base_cells=1, max_bisection_depth=0, max_cells=1),
        )
        with self.assertRaises(ValueError):
            replace(
                record,
                classification="continuous_no_initial_trapped_sphere",
                compact_support_compactness_upper=Q(2),
                exterior_compactness_upper=Q(2),
                compactness_margin=Q(-1),
                inconclusive_reason=None,
                missing_theorem=None,
            )

    def test_direct_compactness_formula_on_minkowski_is_zero(self) -> None:
        compactness = sgbl_misner_sharp_compactness_interval(
            Interval.singleton(10), Interval.singleton(1), Interval.singleton(0)
        )
        self.assertTrue(compactness.is_singleton())
        self.assertEqual(compactness.lower, 0)

    def test_compactness_domain_assuming_c_at_most_one_is_a_domain_stop(self) -> None:
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            SGBLValidatedODEPolicy(compactness_domain=interval(Q(-1), Q(1)))
        self.assertEqual(stopped.exception.reason, "domain_error")

    def test_compactness_domain_missing_zero_is_a_domain_stop(self) -> None:
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            SGBLValidatedODEPolicy(compactness_domain=interval(Q(2), Q(8)))
        self.assertEqual(stopped.exception.reason, "domain_error")

    def test_angular_momentum_domain_missing_zero_is_a_domain_stop(self) -> None:
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            SGBLValidatedODEPolicy(angular_momentum_domain=interval(Q(1), Q(8)))
        self.assertEqual(stopped.exception.reason, "domain_error")


def _natural_rhs(
    radius: Interval,
    radial_metric: Interval,
    angular_extrinsic_curvature: Interval,
    spec: SGBLExactInitialSlice,
) -> tuple[Interval, Interval]:
    coefficients = sgbl_interval_affine_coefficients(
        radius=radius,
        radial_metric=radial_metric,
        angular_extrinsic_curvature=angular_extrinsic_curvature,
        spec=spec,
    )
    return (
        -coefficients["hamiltonian_constant"] / coefficients["hamiltonian_lambda_r"],
        -coefficients["momentum_constant"] / coefficients["momentum_k_r"],
    )


class SGBLMeanValueRHSTests(unittest.TestCase):
    def test_singleton_tangent_is_exact_and_contains_a_finite_difference(self) -> None:
        spec = SGBLExactInitialSlice()
        radius = Interval.singleton(spec.support_minimum)
        lam = Interval.singleton(1)
        k = Interval.singleton(0)
        jacobian = sgbl_constraint_rhs_state_jacobian(radius, lam, k, spec)
        self.assertFalse(jacobian["d_lambda_r_d_lambda"].contains_zero())
        self.assertLessEqual(jacobian["d_lambda_r_d_lambda"].width(), Q(1, 2**32))
        step = Q(1, 64)
        plus = _natural_rhs(radius, Interval.singleton(1 + step), k, spec)
        minus = _natural_rhs(radius, Interval.singleton(1 - step), k, spec)
        fd_lambda = (plus[0] - minus[0]) / (2 * step)
        fd_k = (plus[1] - minus[1]) / (2 * step)
        remainder = step * step
        self.assertLessEqual(
            abs(fd_lambda.midpoint() - jacobian["d_lambda_r_d_lambda"].midpoint()),
            remainder,
        )
        self.assertLessEqual(
            abs(fd_k.midpoint() - jacobian["d_k_r_d_lambda"].midpoint()),
            remainder,
        )
        seeded = sgbl_generic_affine_coefficients(
            radius=radius,
            radial_metric=IntervalFirstTangent.seed(lam),
            angular_extrinsic_curvature=IntervalFirstTangent.constant(k),
            spec=spec,
        )
        self.assertIsInstance(seeded["hamiltonian_lambda_r"], IntervalFirstTangent)

    def test_centered_mean_value_intersects_the_natural_extension(self) -> None:
        spec = SGBLExactInitialSlice()
        radius = Interval.singleton(Q(21, 2))
        lam = interval(Q(15, 16), Q(17, 16))
        k = interval(Q(-1, 32), Q(1, 32))
        natural = _natural_rhs(radius, lam, k, spec)
        centered = sgbl_centered_mean_value_rhs(radius, lam, k, spec)
        sample = _natural_rhs(
            radius, Interval.singleton(lam.midpoint()), Interval.singleton(k.midpoint()), spec
        )
        for index in range(2):
            overlap = Interval(
                max(natural[index].lower, centered[index].lower),
                min(natural[index].upper, centered[index].upper),
            )
            self.assertLessEqual(overlap.lower, overlap.upper)
            self.assertTrue(sample[index].subset_of(natural[index]))
            self.assertTrue(sample[index].subset_of(centered[index]))
            self.assertLessEqual(overlap.width(), natural[index].width())

    def test_sign_flipped_lambda_derivative_is_not_the_tangent_enclosure(self) -> None:
        spec = SGBLExactInitialSlice()
        radius = Interval.singleton(Q(21, 2))
        jacobian = sgbl_constraint_rhs_state_jacobian(
            radius, Interval.singleton(1), Interval.singleton(0), spec
        )
        derivative = jacobian["d_lambda_r_d_lambda"]
        self.assertFalse(derivative.contains_zero())
        flipped = Interval(-derivative.upper, -derivative.lower)
        self.assertFalse(derivative.subset_of(flipped))
        self.assertFalse(flipped.subset_of(derivative))

    def test_jacobian_seed_on_a_nonpositive_lambda_is_a_domain_stop(self) -> None:
        spec = SGBLExactInitialSlice()
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_constraint_rhs_state_jacobian(
                Interval.singleton(Q(21, 2)),
                interval(Q(-1, 8), Q(1, 8)),
                Interval.singleton(0),
                spec,
            )
        self.assertEqual(stopped.exception.reason, "domain_error")

    def test_centered_compactness_rhs_intersects_the_natural_chain_rule(self) -> None:
        spec = SGBLExactInitialSlice()
        radius = Interval.singleton(Q(21, 2))
        lam = interval(Q(15, 16), Q(17, 16))
        k = interval(Q(-1, 32), Q(1, 32))
        natural_lambda_r, natural_k_r = _natural_rhs(radius, lam, k, spec)
        natural_c_r = sgbl_misner_sharp_compactness_radial_derivative(
            radius, lam, k, natural_lambda_r, natural_k_r
        )
        centered = sgbl_centered_mean_value_compactness_rhs(radius, lam, k, spec)
        sample_lambda_r, sample_k_r = _natural_rhs(
            radius, Interval.singleton(lam.midpoint()), Interval.singleton(k.midpoint()), spec
        )
        sample = sgbl_misner_sharp_compactness_radial_derivative(
            radius,
            Interval.singleton(lam.midpoint()),
            Interval.singleton(k.midpoint()),
            sample_lambda_r,
            sample_k_r,
        )
        overlap = Interval(
            max(natural_c_r.lower, centered.lower),
            min(natural_c_r.upper, centered.upper),
        )
        self.assertLessEqual(overlap.lower, overlap.upper)
        self.assertTrue(sample.subset_of(natural_c_r))
        self.assertTrue(sample.subset_of(centered))
        self.assertLessEqual(overlap.width(), natural_c_r.width())


class SGBLCompactnessAugmentationTests(unittest.TestCase):
    def test_singleton_chain_rule_identity_is_exact_zero(self) -> None:
        spec = SGBLExactInitialSlice()
        radius = Interval.singleton(spec.support_minimum)
        lam = Interval.singleton(1)
        k = Interval.singleton(0)
        lambda_r, k_r = _natural_rhs(radius, lam, k, spec)
        identity = sgbl_misner_sharp_compactness_chain_rule_identity(
            radius, lam, k, lambda_r, k_r
        )
        self.assertEqual(lambda_r.lower, 0)
        self.assertEqual(k_r.lower, 0)
        self.assertTrue(identity["chain_rule_residual"].is_singleton())
        self.assertEqual(identity["chain_rule_residual"].lower, 0)
        self.assertTrue(identity["invariant_residual"].is_singleton())
        self.assertEqual(identity["invariant_residual"].lower, 0)
        self.assertEqual(identity["algebraic_compactness"].lower, 0)
        self.assertEqual(identity["formula"].lower, 0)

    def test_interior_singleton_chain_rule_matches_the_closed_form(self) -> None:
        radius = Interval.singleton(Q(21, 2))
        lam = Interval.singleton(1)
        k = Interval.singleton(Q(-1, 8))
        lambda_r = Interval.singleton(Q(1, 4))
        k_r = Interval.singleton(Q(-1, 16))
        identity = sgbl_misner_sharp_compactness_chain_rule_identity(
            radius, lam, k, lambda_r, k_r
        )
        self.assertTrue(identity["chain_rule_residual"].is_singleton())
        self.assertEqual(identity["chain_rule_residual"].lower, 0)
        self.assertTrue(identity["invariant_residual"].is_singleton())
        self.assertEqual(identity["invariant_residual"].lower, 0)
        expected = (
            2 * radius.lower * (k.lower ** 2)
            + 2 * (radius.lower ** 2) * k.lower * k_r.lower
            + 2 * lambda_r.lower / (lam.lower ** 3)
        )
        self.assertEqual(identity["formula"].lower, expected)
        self.assertEqual(identity["jet_tangent"].lower, expected)
        self.assertNotEqual(expected, 0)
        affine_lambda_r, affine_k_r = _natural_rhs(
            radius, lam, k, SGBLExactInitialSlice()
        )
        affine_identity = sgbl_misner_sharp_compactness_chain_rule_identity(
            radius, lam, k, affine_lambda_r, affine_k_r
        )
        self.assertTrue(affine_identity["chain_rule_residual"].contains_zero())
        self.assertTrue(affine_identity["invariant_residual"].contains_zero())

    def test_broken_sign_on_lambda_term_is_not_the_chain_rule(self) -> None:
        radius = Interval.singleton(Q(21, 2))
        lam = Interval.singleton(1)
        k = Interval.singleton(Q(-1, 8))
        lambda_r = Interval.singleton(Q(1, 4))
        k_r = Interval.singleton(Q(-1, 16))
        identity = sgbl_misner_sharp_compactness_chain_rule_identity(
            radius, lam, k, lambda_r, k_r
        )
        flipped = (
            2 * radius.lower * (k.lower ** 2)
            + 2 * (radius.lower ** 2) * k.lower * k_r.lower
            - 2 * lambda_r.lower / (lam.lower ** 3)
        )
        self.assertNotEqual(flipped, identity["formula"].lower)
        self.assertNotEqual(flipped, identity["jet_tangent"].lower)
        residual = identity["jet_tangent"] - Interval.singleton(flipped)
        self.assertFalse(residual.is_singleton() and residual.lower == 0)
        self.assertFalse(residual.contains_zero())

    def test_broken_k_term_sign_is_not_the_chain_rule(self) -> None:
        radius = Interval.singleton(Q(21, 2))
        lam = Interval.singleton(1)
        k = Interval.singleton(Q(-1, 8))
        lambda_r = Interval.singleton(Q(1, 4))
        k_r = Interval.singleton(Q(-1, 16))
        identity = sgbl_misner_sharp_compactness_chain_rule_identity(
            radius, lam, k, lambda_r, k_r
        )
        flipped = (
            2 * radius.lower * (k.lower ** 2)
            - 2 * (radius.lower ** 2) * k.lower * k_r.lower
            + 2 * lambda_r.lower / (lam.lower ** 3)
        )
        self.assertNotEqual(flipped, identity["formula"].lower)
        residual = identity["jet_tangent"] - Interval.singleton(flipped)
        self.assertFalse(residual.is_singleton() and residual.lower == 0)
        self.assertFalse(residual.contains_zero())

    def test_direct_enclosure_contains_the_correlated_graph(self) -> None:
        record = sgbl_continuous_initial_compactness(
            SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        )
        self.assertIsNotNone(record.ode)
        self.assertTrue(record.ode.cells)
        for cell in record.ode.cells:
            self.assertTrue(cell.compactness.subset_of(cell.direct_compactness))
            self.assertTrue(cell.compactness.subset_of(cell.centered_compactness))
            self.assertTrue(cell.compactness.subset_of(cell.propagated_compactness))
            self.assertTrue(cell.invariant_residual.contains_zero())
            fresh_residual = sgbl_misner_sharp_compactness_invariant_residual(
                cell.radius, cell.lambda_box, cell.k_box, cell.compactness
            )
            self.assertTrue(fresh_residual.contains_zero())
        self.assertLessEqual(
            record.compact_support_compactness_upper,
            record.direct_compact_support_compactness_upper,
        )

    def test_tampered_pass_with_invariant_missing_zero_is_refused(self) -> None:
        record = sgbl_continuous_initial_compactness(
            SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        )
        with self.assertRaises(ValueError):
            replace(
                record,
                compactness_invariant_residual=interval(Q(1, 8), Q(1, 4)),
            )

    def test_jacobian_carries_compactness_derivatives(self) -> None:
        spec = SGBLExactInitialSlice()
        jacobian = sgbl_constraint_rhs_state_jacobian(
            Interval.singleton(Q(21, 2)),
            Interval.singleton(1),
            Interval.singleton(Q(-1, 8)),
            spec,
        )
        self.assertIn("d_compactness_r_d_lambda", jacobian)
        self.assertIn("d_compactness_r_d_k", jacobian)
        self.assertGreaterEqual(jacobian["d_compactness_r_d_lambda"].width(), 0)
        flipped = Interval(
            -jacobian["d_compactness_r_d_lambda"].upper,
            -jacobian["d_compactness_r_d_lambda"].lower,
        )
        if not jacobian["d_compactness_r_d_lambda"].contains_zero():
            self.assertFalse(jacobian["d_compactness_r_d_lambda"].subset_of(flipped))


class SGBLAlgebraicCompactnessMeanValueTests(unittest.TestCase):
    def test_singleton_centered_algebraic_c_matches_the_closed_form(self) -> None:
        radius = Interval.singleton(Q(21, 2))
        lam = Interval.singleton(Q(5, 4))
        k = Interval.singleton(Q(-1, 8))
        natural = sgbl_misner_sharp_compactness_interval(radius, lam, k)
        centered = sgbl_centered_mean_value_algebraic_compactness(radius, lam, k)
        gradient = sgbl_misner_sharp_compactness_state_gradient(radius, lam, k)
        expected = 1 + (radius.lower ** 2) * (k.lower ** 2) - 1 / (lam.lower ** 2)
        self.assertTrue(natural.is_singleton())
        self.assertTrue(centered.is_singleton())
        self.assertEqual(natural.lower, expected)
        self.assertEqual(centered.lower, expected)
        self.assertEqual(gradient["d_c_d_r"].lower, 2 * radius.lower * (k.lower ** 2))
        self.assertEqual(gradient["d_c_d_lambda"].lower, 2 / (lam.lower ** 3))
        self.assertEqual(gradient["d_c_d_k"].lower, 2 * (radius.lower ** 2) * k.lower)
        r_jet = IntervalFirstTangent.seed(radius)
        lam_const = IntervalFirstTangent.constant(lam)
        k_const = IntervalFirstTangent.constant(k)
        compactness_r = (
            IntervalFirstTangent.constant(1)
            + (r_jet ** 2) * (k_const ** 2)
            - (lam_const ** 2).reciprocal()
        )
        self.assertEqual(compactness_r.tangent.lower, gradient["d_c_d_r"].lower)
        lam_jet = IntervalFirstTangent.seed(lam)
        compactness_lambda = (
            IntervalFirstTangent.constant(1)
            + (IntervalFirstTangent.constant(radius) ** 2) * (k_const ** 2)
            - (lam_jet ** 2).reciprocal()
        )
        self.assertEqual(compactness_lambda.tangent.lower, gradient["d_c_d_lambda"].lower)
        k_jet = IntervalFirstTangent.seed(k)
        compactness_k = (
            IntervalFirstTangent.constant(1)
            + (IntervalFirstTangent.constant(radius) ** 2) * (k_jet ** 2)
            - (lam_const ** 2).reciprocal()
        )
        self.assertEqual(compactness_k.tangent.lower, gradient["d_c_d_k"].lower)

    def test_dense_rational_samples_sit_in_natural_and_centered_enclosures(self) -> None:
        radius = interval(Q(21, 2), Q(11))
        lam = interval(Q(7, 8), Q(9, 8))
        k = interval(Q(-1, 8), Q(1, 4))
        natural = sgbl_misner_sharp_compactness_interval(radius, lam, k)
        centered = sgbl_centered_mean_value_algebraic_compactness(radius, lam, k)
        overlap = Interval(
            max(natural.lower, centered.lower),
            min(natural.upper, centered.upper),
        )
        self.assertLessEqual(overlap.lower, overlap.upper)
        radii = (radius.lower, radius.midpoint(), radius.upper)
        lambdas = (lam.lower, lam.midpoint(), lam.upper)
        ks = (k.lower, k.midpoint(), k.upper)
        for sample_r in radii:
            for sample_lam in lambdas:
                for sample_k in ks:
                    point = sgbl_misner_sharp_compactness_interval(
                        Interval.singleton(sample_r),
                        Interval.singleton(sample_lam),
                        Interval.singleton(sample_k),
                    )
                    self.assertTrue(point.subset_of(natural))
                    self.assertTrue(point.subset_of(centered))
                    self.assertTrue(point.subset_of(overlap))

    def test_centered_algebraic_enclosure_is_not_wider_after_intersection(self) -> None:
        radius = interval(Q(45, 4), Q(23, 2))
        lam = interval(Q(3, 4), Q(5, 4))
        k = interval(Q(-1, 4), Q(1, 4))
        enclosures = sgbl_intersect_algebraic_compactness_enclosures(radius, lam, k)
        self.assertTrue(enclosures["algebraic"].subset_of(enclosures["natural"]))
        self.assertTrue(enclosures["algebraic"].subset_of(enclosures["centered"]))
        self.assertLessEqual(enclosures["algebraic"].width(), enclosures["natural"].width())
        self.assertLessEqual(enclosures["algebraic"].width(), enclosures["centered"].width())

    def test_wrong_lambda_partial_is_not_the_mean_value_enclosure(self) -> None:
        radius = interval(Q(21, 2), Q(11))
        lam = interval(Q(1, 2), Q(2))
        k = interval(Q(1, 8), Q(1, 4))
        centered = sgbl_centered_mean_value_algebraic_compactness(radius, lam, k)
        r_mid = Interval.singleton(radius.midpoint())
        lam_mid = Interval.singleton(lam.midpoint())
        k_mid = Interval.singleton(k.midpoint())
        center = sgbl_misner_sharp_compactness_interval(r_mid, lam_mid, k_mid)
        box_gradient = sgbl_misner_sharp_compactness_state_gradient(radius, lam, k)
        mid_gradient = sgbl_misner_sharp_compactness_state_gradient(r_mid, lam_mid, k_mid)
        flipped_lambda = Interval(
            -box_gradient["d_c_d_lambda"].upper, -box_gradient["d_c_d_lambda"].lower
        )
        edge = sgbl_misner_sharp_compactness_interval(r_mid, Interval.singleton(lam.upper), k_mid)
        self.assertTrue(edge.subset_of(centered))
        self.assertGreater(edge.lower, center.lower)
        wrong_taylor = (
            center.lower
            + (-mid_gradient["d_c_d_lambda"].lower) * (lam.upper - lam.midpoint())
        )
        self.assertLess(wrong_taylor, center.lower)
        self.assertFalse(box_gradient["d_c_d_lambda"].contains_zero())
        self.assertFalse(flipped_lambda.contains_zero())
        self.assertFalse(box_gradient["d_c_d_lambda"].subset_of(flipped_lambda))
        self.assertFalse(flipped_lambda.subset_of(box_gradient["d_c_d_lambda"]))

    def test_disjoint_propagated_interval_is_an_inconsistent_intersection(self) -> None:
        radius = Interval.singleton(Q(21, 2))
        lam = Interval.singleton(1)
        k = Interval.singleton(0)
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_intersect_algebraic_compactness_enclosures(
                radius, lam, k, propagated=interval(Q(1, 2), Q(1))
            )
        self.assertEqual(stopped.exception.reason, "interval_inconclusive")
        self.assertEqual(
            stopped.exception.payload["obstruction"],
            "compactness_invariant_misses_origin",
        )

    def test_wrong_gradient_can_miss_the_natural_enclosure(self) -> None:
        radius = Interval.singleton(Q(11))
        lam = interval(Q(1, 2), Q(2))
        k = Interval.singleton(Q(1, 4))
        natural = sgbl_misner_sharp_compactness_interval(radius, lam, k)
        r_mid = Interval.singleton(radius.midpoint())
        lam_mid = Interval.singleton(lam.midpoint())
        k_mid = Interval.singleton(k.midpoint())
        center = sgbl_misner_sharp_compactness_interval(r_mid, lam_mid, k_mid)
        mid_gradient = sgbl_misner_sharp_compactness_state_gradient(r_mid, lam_mid, k_mid)
        sample = sgbl_misner_sharp_compactness_interval(
            r_mid, Interval.singleton(lam.lower), k_mid
        )
        self.assertTrue(sample.subset_of(natural))
        self.assertLess(sample.lower, center.lower)
        wrong_taylor = (
            center.lower
            + (-mid_gradient["d_c_d_lambda"].lower) * (lam.lower - lam.midpoint())
        )
        self.assertGreater(wrong_taylor, center.lower)
        self.assertNotEqual(wrong_taylor, sample.lower)

    def test_tampered_pass_cannot_drop_the_centered_intersection(self) -> None:
        record = sgbl_continuous_initial_compactness(
            SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        )
        self.assertTrue(record.continuous_no_initial_trapped_sphere)
        widest = max(cell.direct_compactness.upper for cell in record.ode.cells)
        with self.assertRaises(ValueError):
            replace(
                record,
                compact_support_compactness_upper=max(widest, Q(2)),
                exterior_compactness_upper=max(widest, Q(2)),
                compactness_margin=Q(-1),
            )


class SGBLMisnerSharpChartTests(unittest.TestCase):
    def test_interval_sqrt_is_exact_on_perfect_squares_and_outward_otherwise(self) -> None:
        self.assertEqual(sgbl_interval_sqrt(Interval.singleton(1)), Interval.singleton(1))
        self.assertEqual(sgbl_interval_sqrt(Interval.singleton(4)), Interval.singleton(2))
        self.assertEqual(
            sgbl_interval_sqrt(Interval.singleton(Q(4, 9))), Interval.singleton(Q(2, 3))
        )
        root_two = sgbl_interval_sqrt(Interval.singleton(2))
        self.assertLessEqual(root_two.lower ** 2, 2)
        self.assertGreaterEqual(root_two.upper ** 2, 2)
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_interval_sqrt(interval(Q(-1), Q(1)))
        self.assertEqual(stopped.exception.reason, "interval_inconclusive")
        self.assertEqual(
            stopped.exception.payload["obstruction"],
            "chart_denominator_not_strictly_positive",
        )

    def test_ck_singleton_roundtrip_is_exact_on_minkowski(self) -> None:
        radius = Interval.singleton(10)
        lam = Interval.singleton(1)
        k = Interval.singleton(0)
        roundtrip = sgbl_misner_sharp_ck_chart_roundtrip(radius, lam, k)
        self.assertEqual(roundtrip["compactness"].lower, 0)
        self.assertEqual(roundtrip["denominator"].lower, 1)
        self.assertTrue(roundtrip["lambda_residual"].is_singleton())
        self.assertEqual(roundtrip["lambda_residual"].lower, 0)
        self.assertTrue(roundtrip["k_residual"].is_singleton())
        self.assertEqual(roundtrip["k_residual"].lower, 0)
        state = sgbl_misner_sharp_ck_chart_state(radius, Interval.singleton(0), k)
        self.assertEqual(state["radial_metric"].lower, 1)
        self.assertEqual(state["negative_sqrt_branch"].upper, -1)

    def test_ck_singleton_roundtrip_and_wrong_c_r_sign(self) -> None:
        radius = Interval.singleton(Q(21, 2))
        lam = Interval.singleton(1)
        k = Interval.singleton(Q(-1, 8))
        roundtrip = sgbl_misner_sharp_ck_chart_roundtrip(radius, lam, k)
        self.assertTrue(roundtrip["k_residual"].is_singleton())
        self.assertEqual(roundtrip["k_residual"].lower, 0)
        self.assertTrue(roundtrip["lambda_residual"].contains_zero())
        self.assertTrue(roundtrip["denominator_residual"].contains_zero())
        k_r = Interval.singleton(Q(-1, 16))
        lambda_r = Interval.singleton(Q(1, 4))
        formula = sgbl_misner_sharp_compactness_radial_derivative(
            radius, lam, k, lambda_r, k_r
        )
        flipped = (
            2 * radius.lower * (k.lower ** 2)
            + 2 * (radius.lower ** 2) * k.lower * k_r.lower
            - 2 * lambda_r.lower / (lam.lower ** 3)
        )
        self.assertNotEqual(flipped, formula.lower)

    def test_ck_lost_denominator_is_fail_closed(self) -> None:
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_misner_sharp_ck_chart_state(
                Interval.singleton(10), Interval.singleton(2), Interval.singleton(0)
            )
        self.assertEqual(
            stopped.exception.payload["obstruction"],
            "chart_denominator_not_strictly_positive",
        )

    def test_ck_direct_containment_on_a_box(self) -> None:
        radius = Interval.singleton(Q(21, 2))
        k = interval(Q(-1, 16), Q(1, 16))
        compactness = interval(Q(-1, 8), Q(1, 8))
        state = sgbl_misner_sharp_ck_chart_state(radius, compactness, k)
        algebraic = sgbl_misner_sharp_compactness_interval(
            radius, state["radial_metric"], k
        )
        residual = compactness - algebraic
        self.assertTrue(residual.contains_zero())

    def test_singleton_chart_roundtrip_is_exact_on_minkowski(self) -> None:
        radius = Interval.singleton(10)
        lam = Interval.singleton(1)
        k = Interval.singleton(0)
        roundtrip = sgbl_misner_sharp_chart_roundtrip(radius, lam, k)
        self.assertEqual(roundtrip["compactness"].lower, 0)
        self.assertEqual(roundtrip["angular_momentum"].lower, 0)
        self.assertEqual(roundtrip["denominator"].lower, 1)
        self.assertEqual(roundtrip["lambda_residual"].lower, 0)
        self.assertEqual(roundtrip["k_residual"].lower, 0)
        self.assertTrue(roundtrip["lambda_residual"].is_singleton())
        self.assertTrue(roundtrip["k_residual"].is_singleton())

    def test_singleton_chart_roundtrip_and_j_r_identity(self) -> None:
        radius = Interval.singleton(Q(21, 2))
        lam = Interval.singleton(1)
        k = Interval.singleton(Q(-1, 8))
        roundtrip = sgbl_misner_sharp_chart_roundtrip(radius, lam, k)
        self.assertTrue(roundtrip["k_residual"].is_singleton())
        self.assertEqual(roundtrip["k_residual"].lower, 0)
        self.assertTrue(roundtrip["lambda_residual"].contains_zero())
        self.assertTrue(roundtrip["denominator_residual"].contains_zero())
        angular_momentum = sgbl_misner_sharp_angular_momentum_interval(radius, k)
        self.assertEqual(angular_momentum.lower, (radius.lower ** 3) * k.lower)
        k_r = Interval.singleton(Q(-1, 16))
        formula = sgbl_misner_sharp_angular_momentum_radial_derivative(
            radius, angular_momentum, k_r
        )
        expanded = 3 * (radius.lower ** 2) * k.lower + (radius.lower ** 3) * k_r.lower
        self.assertTrue(formula.is_singleton())
        self.assertEqual(formula.lower, expanded)
        flipped = (
            -3 * angular_momentum.lower / radius.lower + (radius.lower ** 3) * k_r.lower
        )
        self.assertNotEqual(flipped, formula.lower)

    def test_lost_denominator_is_fail_closed(self) -> None:
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_misner_sharp_chart_state(
                Interval.singleton(10), Interval.singleton(2), Interval.singleton(0)
            )
        self.assertEqual(stopped.exception.reason, "interval_inconclusive")
        self.assertEqual(
            stopped.exception.payload["obstruction"],
            "chart_denominator_not_strictly_positive",
        )

    def test_negative_sqrt_branch_is_not_lambda(self) -> None:
        state = sgbl_misner_sharp_chart_state(
            Interval.singleton(10), Interval.singleton(0), Interval.singleton(0)
        )
        self.assertEqual(state["radial_metric"].lower, 1)
        self.assertTrue(state["radial_metric"].strictly_positive())
        self.assertEqual(state["negative_sqrt_branch"].upper, -1)
        self.assertFalse(state["negative_sqrt_branch"].strictly_positive())
        self.assertNotEqual(state["radial_metric"], state["negative_sqrt_branch"])

    def test_zero_containing_d_cannot_be_hidden_by_a_negative_sqrt(self) -> None:
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_interval_sqrt(interval(Q(-1, 8), Q(1, 8)))
        self.assertEqual(
            stopped.exception.payload["obstruction"],
            "chart_denominator_not_strictly_positive",
        )

    def test_tampered_pass_from_product_box_graph_is_refused(self) -> None:
        record = sgbl_continuous_initial_compactness(
            SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        )
        with self.assertRaises(ValueError):
            replace(record, chart_coordinates=False)
        with self.assertRaises(ValueError):
            replace(record, chart_kind="c_j")
        lambda_k = sgbl_validated_lambda_k_constraint_ode(
            SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        )
        self.assertFalse(lambda_k.chart_coordinates)
        self.assertEqual(lambda_k.chart_kind, "lambda_k")
        cj = sgbl_validated_cj_constraint_ode(
            SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        )
        self.assertEqual(cj.chart_kind, "c_j")


class SGBLODEInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        cls.record = sgbl_continuous_initial_compactness(cls.spec)
        cls.cells = cls.record.ode.cells

    def test_partition_functions_concatenate_left_and_right_siblings(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_initial_health as owner

        combiner = getsource(owner._combine_partition)
        self.assertIn("left_cells + right_cells", combiner)
        for name in ("_partition_cell", "_partition_chart_cell", "_partition_ck_cell"):
            source = getsource(getattr(owner, name))
            self.assertIn("_combine_partition", source)
            self.assertIn("right_cells, right_reason", source)
            self.assertNotIn("return _partition_cell(\n        mid", source)
            self.assertNotIn("return _partition_chart_cell(\n        mid", source)
            self.assertNotIn("return _partition_ck_cell(\n        mid", source)

    def test_validator_rejects_the_historical_dropped_left_inventory(self) -> None:
        origin = self.spec.support_minimum
        sgbl_validate_ode_cell_inventory(
            self.cells, origin=origin, terminus=self.spec.support_maximum
        )
        self.assertGreaterEqual(len(self.cells), 2)
        dropped_left = self.cells[1:]
        with self.assertRaises(ValueError):
            sgbl_validate_ode_cell_inventory(dropped_left, origin=origin)

    def test_validator_rejects_a_gapped_inventory(self) -> None:
        gapped = (self.cells[0], self.cells[2])
        with self.assertRaises(ValueError):
            sgbl_validate_ode_cell_inventory(
                gapped, origin=self.cells[0].radius.lower
            )

    def test_complete_ode_record_cannot_omit_a_prefix_or_suffix(self) -> None:
        with self.assertRaises(ValueError):
            replace(self.record.ode, cells=self.cells[:-1], obstruction=None)
        prefix = replace(
            self.record.ode,
            cells=self.cells[:4],
            obstruction="picard_strict_self_map_failed",
        )
        self.assertFalse(prefix.covers_complete_support)
        self.assertEqual(prefix.coverage_left, self.spec.support_minimum)
        self.assertEqual(prefix.coverage_right, prefix.cells[-1].radius.upper)

    def test_pass_cannot_claim_complete_coverage_on_a_prefix(self) -> None:
        prefix = replace(
            self.record.ode,
            cells=self.cells[:4],
            obstruction="picard_strict_self_map_failed",
        )
        with self.assertRaises(ValueError):
            replace(
                self.record,
                ode=prefix,
            )

    def test_cross_chart_complete_tiling_on_small_amplitude(self) -> None:
        lambda_k = sgbl_validated_lambda_k_constraint_ode(self.spec)
        cj = sgbl_validated_cj_constraint_ode(self.spec)
        for ode in (self.record.ode, lambda_k, cj):
            self.assertTrue(ode.covers_complete_support)
            self.assertTrue(ode.tiles_compact_support)
            self.assertEqual(ode.coverage_left, self.spec.support_minimum)
            self.assertEqual(ode.coverage_right, self.spec.support_maximum)
            sgbl_validate_ode_cell_inventory(
                ode.cells,
                origin=self.spec.support_minimum,
                terminus=self.spec.support_maximum,
            )


class SGBLCenterMonitorTests(unittest.TestCase):
    def test_unexecuted_monitor_is_not_a_pre_holdout_pass(self) -> None:
        contract = sgbl_evolving_center_monitor_contract()
        self.assertIs(contract.executed, False)
        self.assertIs(contract.pre_holdout_pass, False)
        self.assertTrue(contract.unexecuted_monitor_is_not_a_pre_holdout_pass)
        self.assertEqual(contract.classification, "monitor_not_executed")
        self.assertFalse(contract.evolving_center_after_matter_arrives)
        self.assertIs(contract.SGBL_branch_owned_and_healthy, False)
        self.assertIn("sgb1_ctl1_center.sgbl_initial_center_series", contract.grounded_in)
        self.assertIn("regular_center.LaurentSeries", contract.grounded_in)
        self.assertTrue(contract.initial_empty_center_series["initial_empty_center"])
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_execute_evolving_center_monitor()
        self.assertEqual(stopped.exception.reason, "monitor_not_executed")
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_execute_evolving_center_monitor({"trajectory": "forged"})
        self.assertEqual(stopped.exception.reason, "missing_trajectory_owner")
        with self.assertRaises(TypeError):
            replace(contract, pre_holdout_pass=True)
        with self.assertRaises(ValueError):
            SGBLEvolvingCenterMonitorContract(
                executed=True,
                classification="monitor_not_executed",
                grounded_in=contract.grounded_in,
                required_checks=contract.required_checks,
                initial_empty_center_series=dict(contract.initial_empty_center_series),
                inconclusive_reason=None,
                missing_theorem=None,
            )

    def test_aggregate_health_stays_false(self) -> None:
        compactness = sgbl_continuous_initial_compactness(
            SGBLExactInitialSlice(),
            policy=SGBLValidatedODEPolicy(base_cells=1, max_bisection_depth=0, max_cells=1),
        )
        monitor = sgbl_evolving_center_monitor_contract()
        gate = sgbl_initial_health_gate(compactness, monitor)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["pre_holdout_pass"], False)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertTrue(gate["sampled_nodes_are_not_the_certificate"])
        self.assertTrue(gate["unexecuted_monitor_is_not_a_pre_holdout_pass"])
        self.assertIs(SGBLContinuousCompactnessRecord.__dataclass_params__.frozen, True)
        self.assertEqual(DECLARED_ODE_POLICY.base_cells, DECLARED_BASE_CELLS)
        self.assertEqual(DECLARED_ODE_POLICY.compactness_domain, DECLARED_C_DOMAIN)
        self.assertEqual(DECLARED_ODE_POLICY.angular_momentum_domain, DECLARED_J_DOMAIN)
        self.assertGreater(DECLARED_C_DOMAIN.upper, 1)
        self.assertLess(DECLARED_J_DOMAIN.lower, 0)
        self.assertGreater(DECLARED_J_DOMAIN.upper, 0)


class SGBLInitialHealthImportTests(unittest.TestCase):
    def test_module_does_not_import_qr_health_or_use_a_gb_mass_proxy(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_initial_health as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        forbidden = {
            "proto1_compactness_upper_bound",
            "FGCQRActionParameters",
            "solve_initial_data",
            "validate_regular_profile",
            "regular_center_series_certificate",
            "newton_residual_limit",
            "run_fgc_gr0_calibration_v1",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("HYP2", source)
        self.assertNotIn("1e-12", source)
        self.assertNotIn("GB_EFFECTIVE_DENSITY_BOUND", source)


if __name__ == "__main__":
    unittest.main()
