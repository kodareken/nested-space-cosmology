from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    METHODS,
    PRIMARY_METHOD,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from recursive_horizons.fgc.sgb1_ctl1_family import (  # noqa: E402
    CONSTRAINT_INTEGRATOR_NAMES,
)
from recursive_horizons.fgc.sgb1_ctl1_interval_health import (  # noqa: E402
    ADM_SOURCE_CHART as INTERVAL_CHART,
)
from recursive_horizons.fgc.sgb1_ctl1_runtime import (  # noqa: E402
    ADM_SOURCE_CHART,
    LOWER_JET_GROUPS,
    MAX_SYNTHETIC_GRID_POINTS,
    METRIC_DTT_CHART,
    MISSING_LOWER_TWO_JET_INTERFACE,
    TIME_METHOD_SPATIAL_ORDER,
    SGBLRuntime,
    SGBLRuntimeAction,
    SGBLRuntimeStop,
    complete_lower_jets_from_evolution_state,
    sgbl_declared_polynomial_grid,
    sgbl_exact_floating_source_comparison,
    sgbl_floating_source_coefficients,
    sgbl_floating_source_residual,
    sgbl_floating_source_solve,
    sgbl_grid_source_solve,
    sgbl_runtime_health_facts,
    sgbl_synthetic_member_from_source_point,
)
from recursive_horizons.fgc.exact_linear_algebra import inverse  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    SGBLSourceInputs,
    SGBLSourceJacobianSolveStop,
    sgbl_source_coefficients,
    sgbl_source_residual,
    sgbl_source_solve,
)
from recursive_horizons.fgc.spherical_reduction import BASE_FIELD_ORDER  # noqa: E402


def _fixture_a() -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(2, Q(1, 3), Q(-1, 2), Q(1, 5), Q(-1, 7)),
        shift=(Q(1, 3), Q(-2, 5), Q(1, 4), Q(-1, 6), Q(2, 9)),
        radial_metric=(Q(3, 2), Q(1, 8), Q(-1, 3), Q(1, 7), Q(-1, 4)),
        areal_radius=(4, Q(-1, 5), Q(3, 2), Q(1, 9), Q(-2, 5)),
        phi=(Q(1, 2), Q(2, 3), Q(-4, 5), Q(1, 4), Q(-1, 8)),
        chi=(Q(-3, 4), Q(1, 2), Q(3, 5), Q(-2, 7), Q(1, 6)),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=Q(-1, 4),
    )


def _fixture_b() -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=3,
        alpha=(Q(5, 2), Q(-2, 7), Q(3, 4), Q(-1, 8), Q(2, 11)),
        shift=(Q(-2, 5), Q(1, 6), Q(-3, 7), Q(2, 9), Q(-1, 5)),
        radial_metric=(Q(4, 3), Q(-1, 4), Q(2, 5), Q(-1, 9), Q(3, 8)),
        areal_radius=(Q(7, 2), Q(3, 8), Q(-5, 4), Q(2, 7), Q(1, 10)),
        phi=(Q(-2, 3), Q(-5, 8), Q(1, 2), Q(-3, 11), Q(2, 13)),
        chi=(Q(4, 5), Q(-1, 3), Q(-2, 9), Q(3, 10), Q(-1, 12)),
        planck_mass=Q(3, 2),
        scalar_mass=Q(2, 3),
        quartic_coupling=Q(5, 7),
        alpha_gb=Q(2, 5),
    )


def _singular() -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(1, 0, 0, 0, 0),
        shift=(0, 0, 0, 0, 0),
        radial_metric=(1, 0, 0, 0, 0),
        areal_radius=(Q(5, 2), 0, 1, 0, 0),
        phi=(0, 0, -5, 0, 0),
        chi=(0, 0, 0, 0, 0),
        planck_mass=2,
        scalar_mass=1,
        quartic_coupling=1,
        alpha_gb=-Q(1, 4),
    )


def _comparison_index(grid: UniformRadialGrid, radius: float) -> int:
    return int(np.argmin(np.abs(grid.coordinates - radius)))


class SGBLFloatingSourceTests(unittest.TestCase):
    def test_nonflat_nonzero_shift_two_scalar_matches_exact_source(self) -> None:
        for point in (_fixture_a(), _fixture_b()):
            with self.subTest(radius=point.coordinate_radius):
                witness = sgbl_exact_floating_source_comparison(point)
                self.assertLess(witness["residual_infinity_difference"], 1.0e-12)
                self.assertLess(witness["jacobian_infinity_difference"], 1.0e-11)
                self.assertLess(witness["solve_infinity_difference"], 1.0e-11)
                self.assertLess(witness["solved_residual_infinity"], 1.0e-11)
                self.assertEqual(witness["chart"], ADM_SOURCE_CHART)
                self.assertEqual(witness["chart"], INTERVAL_CHART)
                self.assertTrue(witness["independent_chi"])
                self.assertIs(witness["SGBL_branch_owned_and_healthy"], False)
                exact = sgbl_source_coefficients(point)
                floating = sgbl_floating_source_coefficients(point)
                self.assertEqual(floating.classification, "exact_inverse")
                self.assertIsNotNone(floating.inverse)
                self.assertEqual(
                    floating.exact_determinant, exact.jacobian_determinant
                )
                self.assertEqual(
                    floating.exact_inverse,
                    tuple(tuple(row) for row in inverse(exact.jacobian)),
                )
                self.assertIs(floating.SGBL_branch_owned_and_healthy, False)

    def test_solved_root_re_evaluates_the_actual_mhg_residual(self) -> None:
        point = _fixture_a()
        accelerations, residual, coefficients = sgbl_floating_source_solve(point)
        reconstructed = coefficients.evaluate(accelerations)
        self.assertLess(float(np.max(np.abs(reconstructed))), 1.0e-11)
        self.assertLess(float(np.max(np.abs(residual))), 1.0e-11)
        exact = np.asarray(
            [float(entry) for entry in sgbl_source_solve(point)], dtype=np.float64
        )
        self.assertLess(float(np.max(np.abs(accelerations - exact))), 1.0e-11)

    def test_chi_momenta_are_not_dropped(self) -> None:
        point = _fixture_a()
        chi = point.chi
        zeroed = replace(
            point,
            chi=replace(chi, dt=Q(0), dtr=Q(0)),
        )
        residual = sgbl_floating_source_residual(point)
        zeroed_residual = sgbl_floating_source_residual(zeroed)
        self.assertGreater(float(np.max(np.abs(residual - zeroed_residual))), 0.0)

    def test_injected_wrong_chart_is_refused(self) -> None:
        point = _fixture_a()
        with self.assertRaises(SGBLRuntimeStop) as stopped:
            sgbl_floating_source_residual(point, chart=METRIC_DTT_CHART)
        self.assertEqual(stopped.exception.reason, "source_chart_error")
        self.assertNotEqual(METRIC_DTT_CHART, ADM_SOURCE_CHART)
        self.assertEqual(BASE_FIELD_ORDER[0], "h_tt")

    def test_injected_wrong_action_is_refused(self) -> None:
        with self.assertRaises(SGBLRuntimeStop) as stopped:
            SGBLRuntimeAction(
                planck_mass=2,
                scalar_mass=3,
                quartic_coupling=Q(1, 2),
                alpha_gb=Q(-1, 4),
                beta=1,
            )
        self.assertEqual(stopped.exception.reason, "action_identity_error")
        with self.assertRaises(SGBLRuntimeStop) as stopped:
            SGBLRuntimeAction(
                planck_mass=2,
                scalar_mass=3,
                quartic_coupling=Q(1, 2),
                alpha_gb=Q(-1, 4),
                eta=1,
            )
        self.assertEqual(stopped.exception.reason, "action_identity_error")

    def test_injected_sign_flip_does_not_match_the_true_residual(self) -> None:
        point = _fixture_a()
        flipped = replace(point, alpha_gb=-point.alpha_gb)
        residual = sgbl_floating_source_residual(point)
        flipped_residual = sgbl_floating_source_residual(flipped)
        exact = np.asarray(
            [float(entry) for entry in sgbl_source_residual(point)], dtype=np.float64
        )
        self.assertLess(float(np.max(np.abs(residual - exact))), 1.0e-12)
        self.assertGreater(float(np.max(np.abs(flipped_residual - exact))), 1.0e-6)
        self.assertGreater(
            float(np.max(np.abs(flipped_residual - residual))), 1.0e-6
        )

    def test_exact_singular_point_keeps_source_jacobian_ownership(self) -> None:
        coefficients = sgbl_floating_source_coefficients(_singular())
        self.assertEqual(coefficients.classification, "singular_source_jacobian")
        self.assertIsNone(coefficients.inverse)
        self.assertEqual(coefficients.exact_determinant, 0)
        with self.assertRaises(SGBLSourceJacobianSolveStop) as stopped:
            sgbl_floating_source_solve(_singular())
        self.assertEqual(stopped.exception.reason, "singular_source_jacobian")


class SGBLGridCompletionTests(unittest.TestCase):
    def test_evolution_state_without_operator_or_jets_is_blocked(self) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        with self.assertRaises(SGBLRuntimeStop) as stopped:
            complete_lower_jets_from_evolution_state(member.state)
        self.assertEqual(stopped.exception.reason, "missing_lower_two_jets")
        self.assertIn("p_r", MISSING_LOWER_TWO_JET_INTERFACE)
        self.assertEqual(LOWER_JET_GROUPS, ("u", "p", "q", "p_r", "q_r"))

    def test_method_owned_reconstruction_matches_declared_polynomial_jets(self) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        operator = SBPFirstDerivative(member.grid, member.spatial_order)
        p_r, q_r = complete_lower_jets_from_evolution_state(
            member.state, spatial_operator=operator
        )
        self.assertTrue(np.allclose(p_r, member.declared_p_r, rtol=0.0, atol=1.0e-12))
        self.assertTrue(np.allclose(q_r, member.declared_q_r, rtol=0.0, atol=1.0e-12))
        comparator = SBPFirstDerivative(member.grid, 2)
        p_r2, q_r2 = complete_lower_jets_from_evolution_state(
            member.state, spatial_operator=comparator
        )
        self.assertTrue(np.allclose(p_r2, member.declared_p_r, rtol=0.0, atol=1.0e-12))
        self.assertTrue(np.allclose(q_r2, member.declared_q_r, rtol=0.0, atol=1.0e-12))

    def test_reconstructed_interior_source_matches_exact_solve(self) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        operator = SBPFirstDerivative(member.grid, member.spatial_order)
        p_r, q_r = complete_lower_jets_from_evolution_state(
            member.state, spatial_operator=operator
        )
        index = _comparison_index(member.grid, float(point.coordinate_radius))
        solved = sgbl_grid_source_solve(
            u=member.state.u[index : index + 1],
            p=member.state.p[index : index + 1],
            q=member.state.q[index : index + 1],
            p_r=p_r[index : index + 1],
            q_r=q_r[index : index + 1],
            radii=member.grid.coordinates[index : index + 1],
            action=member.action,
        )
        exact = np.asarray(
            [float(entry) for entry in sgbl_source_solve(point)], dtype=np.float64
        )
        self.assertLess(
            float(np.max(np.abs(solved["accelerations"][0] - exact))),
            1.0e-10,
        )
        self.assertLess(solved["solved_residual_infinity"], 1.0e-9)
        self.assertIs(solved["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(solved["source_solve_admission_qualified"], False)
        self.assertTrue(np.isfinite(solved["scaled_residual_witness"]))
        self.assertEqual(solved["chart"], ADM_SOURCE_CHART)

    def test_floating_grid_inversion_failure_is_inconclusive_not_exact_singular(
        self,
    ) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        operator = SBPFirstDerivative(member.grid, member.spatial_order)
        p_r, q_r = complete_lower_jets_from_evolution_state(
            member.state, spatial_operator=operator
        )
        with patch(
            "recursive_horizons.fgc.sgb1_ctl1_runtime.np.linalg.inv",
            side_effect=np.linalg.LinAlgError("injected floating inversion"),
        ), self.assertRaises(SGBLRuntimeStop) as stopped:
            sgbl_grid_source_solve(
                u=member.state.u,
                p=member.state.p,
                q=member.state.q,
                p_r=p_r,
                q_r=q_r,
                radii=member.grid.coordinates,
                action=member.action,
            )
        self.assertEqual(stopped.exception.reason, "source_inversion_inconclusive")

    def test_grid_evidence_arrays_are_immutable(self) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        operator = SBPFirstDerivative(member.grid, member.spatial_order)
        p_r, q_r = complete_lower_jets_from_evolution_state(
            member.state, spatial_operator=operator
        )
        solved = sgbl_grid_source_solve(
            u=member.state.u,
            p=member.state.p,
            q=member.state.q,
            p_r=p_r,
            q_r=q_r,
            radii=member.grid.coordinates,
            action=member.action,
        )
        for name in (
            "accelerations",
            "constant",
            "jacobian",
            "inverse",
            "solved_residual",
        ):
            self.assertFalse(solved[name].flags.writeable)
        coefficients = sgbl_floating_source_coefficients(point)
        self.assertFalse(coefficients.constant.flags.writeable)

    def test_zeros_are_not_a_hidden_completion(self) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        index = _comparison_index(member.grid, float(point.coordinate_radius))
        radius = member.grid.coordinates[index : index + 1]
        zeros = np.zeros((1, 6), dtype=np.float64)
        slice_u = member.state.u[index : index + 1]
        slice_p = member.state.p[index : index + 1]
        slice_q = member.state.q[index : index + 1]
        solved_zero = sgbl_grid_source_solve(
            u=slice_u,
            p=slice_p,
            q=slice_q,
            p_r=zeros,
            q_r=zeros,
            radii=radius,
            action=member.action,
        )
        solved_true = sgbl_grid_source_solve(
            u=slice_u,
            p=slice_p,
            q=slice_q,
            p_r=member.declared_p_r[index : index + 1],
            q_r=member.declared_q_r[index : index + 1],
            radii=radius,
            action=member.action,
        )
        self.assertGreater(
            float(
                np.max(
                    np.abs(
                        solved_zero["accelerations"][0]
                        - solved_true["accelerations"][0]
                    )
                )
            ),
            1.0e-6,
        )


class SGBLTimeMethodTests(unittest.TestCase):
    def test_constraint_integrators_are_not_time_methods(self) -> None:
        self.assertEqual(CONSTRAINT_INTEGRATOR_NAMES, ("RK4", "SSPRK3"))
        self.assertEqual(TIME_METHOD_SPATIAL_ORDER[PRIMARY_METHOD], 4)
        self.assertEqual(TIME_METHOD_SPATIAL_ORDER[COMPARATOR_METHOD], 2)
        self.assertTrue(set(CONSTRAINT_INTEGRATOR_NAMES).isdisjoint(METHODS))
        point = _fixture_a()
        with self.assertRaisesRegex(ValueError, "time method"):
            sgbl_synthetic_member_from_source_point(point, method="RK4")

    def test_health_facts_remain_closed(self) -> None:
        health = sgbl_runtime_health_facts(origin=_fixture_a())
        self.assertIs(health["SGBL_branch_owned_and_healthy"], False)
        self.assertFalse(health["physical_execution_authorized"])
        self.assertEqual(health["cone_certificate"], "unqualified")
        self.assertEqual(
            health["center"]["evolving_center_after_matter_arrives"],
            "unqualified",
        )
        self.assertIsNotNone(health["constraint"])
        self.assertIsNotNone(health["principal"])
        self.assertIs(health["principal"].SGBL_branch_owned_and_healthy, False)

    def test_both_time_methods_advance_on_small_declared_data(self) -> None:
        point = _fixture_a()
        for method in METHODS:
            member = sgbl_synthetic_member_from_source_point(point, method=method)
            runtime = SGBLRuntime(member)
            snapshot = runtime.capture()
            before = snapshot.physical_fingerprint
            receipt = runtime.attempt_step(1.0e-4)
            self.assertEqual(receipt.method, method)
            self.assertEqual(
                receipt.spatial_order, TIME_METHOD_SPATIAL_ORDER[method]
            )
            self.assertNotEqual(runtime.member.physical_fingerprint, before)
            self.assertEqual(runtime.member.counters.accepted_steps, 1)
            self.assertGreater(runtime.member.counters.accepted_stages, 1)
            self.assertEqual(runtime.member.time, 1.0e-4)
            self.assertTrue(runtime.member.not_a_radial_constraint_integrator)
            self.assertIs(runtime.member.SGBL_branch_owned_and_healthy, False)
            self.assertFalse(runtime.member.physical_execution_authorized)
            self.assertLess(receipt.source_residual_infinity, 1.0e-6)
            self.assertTrue(np.isfinite(receipt.scaled_residual_witness))
            self.assertIs(receipt.source_solve_admission_qualified, False)
            self.assertTrue(receipt.synthetic_transition_only)
            self.assertIsNone(runtime.member.declared_p_r)
            self.assertIsNone(runtime.member.declared_q_r)
            runtime.restore(snapshot)
            self.assertEqual(runtime.member.physical_fingerprint, before)
            self.assertEqual(runtime.member.counters.accepted_steps, 0)
            self.assertEqual(runtime.member.time, 0.0)
            self.assertIsNotNone(runtime.member.declared_p_r)
            self.assertIsNotNone(runtime.member.declared_q_r)


class SGBLTransactionTests(unittest.TestCase):
    def test_source_failure_rolls_back_state_monitor_causal_tracers_and_counters(
        self,
    ) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        damaged = np.array(member.state.u, copy=True)
        damaged[:, 3] = 0.0
        runtime = SGBLRuntime(replace(member, state=EvolutionState(
            damaged, member.state.p, member.state.q
        )))
        fingerprint = runtime.member.physical_fingerprint
        counters = runtime.member.counters
        monitor = runtime.member.monitor
        causal = runtime.member.causal
        tracers = runtime.member.tracers
        with self.assertRaises(SGBLRuntimeStop) as stopped:
            runtime.attempt_step(1.0e-4)
        self.assertEqual(stopped.exception.reason, "source_chart_domain")
        restored = runtime.member
        self.assertEqual(restored.physical_fingerprint, fingerprint)
        self.assertEqual(restored.counters, counters)
        self.assertEqual(restored.monitor, monitor)
        self.assertEqual(restored.causal, causal)
        self.assertEqual(restored.tracers, tracers)

    def test_cfl_failure_rolls_back_without_rhs_side_effects(self) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=COMPARATOR_METHOD
        )
        runtime = SGBLRuntime(member)
        fingerprint = runtime.member.physical_fingerprint
        with self.assertRaises(SGBLRuntimeStop) as stopped:
            runtime.attempt_step(10.0)
        self.assertEqual(stopped.exception.reason, "cfl_violation")
        self.assertEqual(runtime.member.physical_fingerprint, fingerprint)
        self.assertEqual(runtime.member.counters.accepted_steps, 0)
        self.assertEqual(runtime.member.counters.rhs_evaluations, 0)

    def test_health_failure_rolls_back(self) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        damaged = np.array(member.state.u, copy=True)
        damaged[:, 0] = -1.0
        runtime = SGBLRuntime(replace(member, state=EvolutionState(
            damaged, member.state.p, member.state.q
        )))
        fingerprint = runtime.member.physical_fingerprint
        with self.assertRaises(SGBLRuntimeStop) as stopped:
            runtime.attempt_step(1.0e-4)
        self.assertEqual(stopped.exception.reason, "health_premise_failed")
        self.assertEqual(runtime.member.physical_fingerprint, fingerprint)
        self.assertEqual(runtime.member.counters.accepted_steps, 0)

    def test_health_and_forbidden_scope_cannot_be_attached(self) -> None:
        point = _fixture_a()
        member = sgbl_synthetic_member_from_source_point(
            point, method=PRIMARY_METHOD
        )
        self.assertIs(member.SGBL_branch_owned_and_healthy, False)
        with self.assertRaises(TypeError):
            replace(member, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(FrozenInstanceError):
            member.time = 1.0  # type: ignore[misc]
        self.assertLessEqual(member.grid.point_count, MAX_SYNTHETIC_GRID_POINTS)
        with self.assertRaises(SGBLRuntimeStop) as stopped:
            sgbl_declared_polynomial_grid(point, UniformRadialGrid(0.0, 1.0, 9))
        self.assertEqual(stopped.exception.reason, "source_chart_domain")


class SGBLRuntimeBoundaryTests(unittest.TestCase):
    def test_module_does_not_import_qr_source_or_old_runners(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_runtime as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        forbidden = {
            "background_from_spherical_state",
            "weak_coupling_health_certificate",
            "canonical_health_monitor_values",
            "solve_accelerations",
            "quantified_implicit_branch_certificate",
            "interval_full_residual_acceleration_jacobian",
            "compact_ref1_interval_principal_certificate",
            "FGCQRActionParameters",
            "solve_initial_data",
            "run_fgc_gr0_calibration_v1",
            "run_fgc_gr0_calibration_v2",
            "diagnose_gr0_grid_accelerations",
            "diagnose_gr0_vectorized_reference_accelerations",
            "GR0DirectSource",
            "vectorized_source",
            "nonlinear_source",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertNotIn("solve_accelerations", source)
        self.assertNotIn("FGCQRActionParameters", source)
        self.assertNotIn("weak_coupling_health_certificate", source)


if __name__ == "__main__":
    unittest.main()
