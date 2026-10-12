from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.exact_interval import Interval  # noqa: E402
from recursive_horizons.fgc.exact_linear_algebra import matrix_multiply  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_reference import (  # noqa: E402
    modified_harmonic_full_residuals,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.sgb1_ctl1_interval_health import (  # noqa: E402
    ADM_SOURCE_CHART,
    SGBLIntervalBox,
    SGBLIntervalInconclusive,
    SGBLIntervalLimits,
    jacobian_contains_point,
    sgbl_enclose_adm_source_jacobian,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    HAT_NORMAL_FACTOR,
    REFERENCE_RADIAL_MINIMUM,
    SGBLSourceInputs,
    SGBLSourceJacobianSolveStop,
    SOURCE_ACCELERATION_ORDER,
    TILDE_NORMAL_FACTOR,
    sgbl_source_coefficients,
    sgbl_source_residual,
    sgbl_source_state,
)
from recursive_horizons.fgc.spherical_reduction import BASE_FIELD_ORDER  # noqa: E402


EXPECTED_DET_A = Q(6519803496633, 6553600000)
ZERO6 = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))


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


def _flat(*, alpha_gb: Q | int = Q(-1, 4)) -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(1, 0, 0, 0, 0),
        shift=(0, 0, 0, 0, 0),
        radial_metric=(1, 0, 0, 0, 0),
        areal_radius=(Q(5, 2), 0, 1, 0, 0),
        phi=(0, 0, 0, 0, 0),
        chi=(0, 0, 0, 0, 0),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=alpha_gb,
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


def _inf_row_sum(matrix: tuple[tuple[Q, ...], ...]) -> Q:
    return max(sum(abs(entry) for entry in row) for row in matrix)


def _metric_dtt_jacobian(point: SGBLSourceInputs) -> tuple[tuple[Q, ...], ...]:
    origin = sgbl_source_residual(point, ZERO6)
    state = sgbl_source_state(point, ZERO6)
    reference = flat_spherical_annulus_reference(
        radial_domain_minimum=REFERENCE_RADIAL_MINIMUM,
    )
    columns = []
    for field in BASE_FIELD_ORDER:
        perturbed = replace(state, **{field: replace(getattr(state, field), dtt=Q(1))})
        residual = tuple(
            modified_harmonic_full_residuals(
                perturbed,
                reference=reference,
                coordinate_radius=point.coordinate_radius,
                tilde_normal_factor=TILDE_NORMAL_FACTOR,
                hat_normal_factor=HAT_NORMAL_FACTOR,
            )["full_residual_vector"]
        )
        columns.append(tuple(entry - origin[row] for row, entry in enumerate(residual)))
    return tuple(tuple(columns[column][row] for column in range(6)) for row in range(6))


class SGBLIntervalInputTests(unittest.TestCase):
    def test_metric_chart_and_float_width_are_refused(self) -> None:
        point = _flat()
        with self.assertRaises(SGBLIntervalInconclusive) as stopped:
            SGBLIntervalBox(center=point, parameter_half_width=0, chart="metric_dtt")
        self.assertEqual(stopped.exception.reason, "interval_chart_error")
        with self.assertRaises(TypeError):
            SGBLIntervalBox(center=point, parameter_half_width=0.0)
        with self.assertRaises(ValueError):
            SGBLIntervalBox(center=_flat(alpha_gb=0), parameter_half_width=0)
        coupling = SGBLIntervalBox(
            center=_flat(alpha_gb=0),
            parameter_half_width=0,
            allow_zero_coupling_control=True,
        )
        self.assertEqual(coupling.center.alpha_gb, 0)
        self.assertEqual(coupling.chart, ADM_SOURCE_CHART)

    def test_width_cap_is_a_resource_stop_not_a_fitted_floor(self) -> None:
        with self.assertRaises(SGBLIntervalInconclusive) as stopped:
            SGBLIntervalBox(
                center=_flat(),
                parameter_half_width=2,
                limits=SGBLIntervalLimits(max_parameter_half_width=1),
            )
        self.assertEqual(stopped.exception.reason, "resource_limit")


class SGBLIntervalSingletonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.point = _fixture_a()
        cls.exact = sgbl_source_coefficients(cls.point)
        cls.record = sgbl_enclose_adm_source_jacobian(
            SGBLIntervalBox(center=cls.point, parameter_half_width=0)
        )

    def test_singleton_reproduces_the_frozen_adm_source_jacobian(self) -> None:
        record = self.record
        self.assertEqual(record.exact_center_det, EXPECTED_DET_A)
        self.assertEqual(record.exact_center_jacobian, self.exact.jacobian)
        self.assertTrue(record.independent_exact_singleton_reduction)
        self.assertTrue(record.source_invertible_over_box)
        self.assertEqual(record.classification, "exact_singleton_inverse")
        self.assertEqual(record.rho_infinity, 0)
        self.assertEqual(record.remainder_bound, 0)
        self.assertEqual(record.acceleration_order, SOURCE_ACCELERATION_ORDER)
        self.assertEqual(record.tilde_normal_factor, Q(4))
        self.assertEqual(record.hat_normal_factor, Q(9))
        self.assertEqual(record.action["beta"], 0)
        self.assertEqual(record.action["eta"], 0)
        self.assertEqual(record.action["alpha_gb"], Q(-1, 4))
        self.assertTrue(
            jacobian_contains_point(record.jacobian_box, self.exact.jacobian)
        )
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(record.strongly_hyperbolic, False)
        with self.assertRaisesRegex(ValueError, "proved interval inverse"):
            replace(record, source_invertible_over_box=False)
        with self.assertRaisesRegex(ValueError, "determinant"):
            replace(record, exact_center_det=record.exact_center_det + 1)
        self.assertIsNone(record.principal_cone_margins)
        self.assertEqual(record.cone_certificate, "unqualified")
        with self.assertRaises(TypeError):
            replace(record, SGBL_branch_owned_and_healthy=True)

    def test_nonzero_shift_two_scalar_center_is_not_a_metric_dtt_jacobian(self) -> None:
        metric = _metric_dtt_jacobian(self.point)
        difference = max(
            abs(metric[row][column] - self.exact.jacobian[row][column])
            for row in range(6)
            for column in range(6)
        )
        self.assertGreater(difference, 1)
        self.assertNotEqual(
            tuple(tuple(row) for row in metric),
            self.exact.jacobian,
        )


class SGBLIntervalProvedEnclosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.point = _flat()
        cls.width = Q(1, 2**20)
        cls.record = sgbl_enclose_adm_source_jacobian(
            SGBLIntervalBox(
                center=cls.point,
                parameter_half_width=cls.width,
                acceleration_half_width=Q(1, 2**10),
            )
        )

    def test_declared_flat_box_proves_rho_below_one_and_contains_samples(self) -> None:
        record = self.record
        self.assertTrue(record.source_invertible_over_box)
        self.assertIsNotNone(record.rho_infinity)
        self.assertLess(record.rho_infinity, 1)
        self.assertGreater(record.rho_infinity, 0)
        self.assertEqual(record.classification, "neumann_inverse_enclosure")
        self.assertIsNone(record.inconclusive_reason)
        self.assertFalse(record.independent_exact_singleton_reduction)
        sampled = sgbl_source_coefficients(
            replace(self.point, phi=replace(self.point.phi, dr=self.width / 2))
        )
        self.assertTrue(
            jacobian_contains_point(record.jacobian_box, sampled.jacobian)
        )
        defect = tuple(
            tuple(
                int(row == column) - entry
                for column, entry in enumerate(
                    matrix_multiply(record.exact_center_inverse, sampled.jacobian)[row]
                )
            )
            for row in range(6)
        )
        self.assertLessEqual(_inf_row_sum(defect), record.rho_infinity)
        inverse_box = record.source_inverse["inverse_enclosure"]
        for row in range(6):
            for column in range(6):
                center = record.exact_center_inverse[row][column]
                self.assertLessEqual(inverse_box[row][column].lower, center)
                self.assertLessEqual(center, inverse_box[row][column].upper)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(record.strongly_hyperbolic, False)


class SGBLIntervalWrappingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.point = _flat()
        cls.width = Q(1, 16)
        cls.record = sgbl_enclose_adm_source_jacobian(
            SGBLIntervalBox(center=cls.point, parameter_half_width=cls.width)
        )

    def test_wide_wrapping_is_inconclusive_and_still_contains_samples(self) -> None:
        record = self.record
        self.assertEqual(record.classification, "interval_inconclusive")
        self.assertEqual(record.inconclusive_reason, "neumann_rho_not_below_one")
        self.assertFalse(record.source_invertible_over_box)
        self.assertIsNone(record.source_inverse)
        sampled = sgbl_source_coefficients(
            replace(self.point, phi=replace(self.point.phi, dr=Q(1, 32)))
        )
        self.assertTrue(
            jacobian_contains_point(record.jacobian_box, sampled.jacobian)
        )
        self.assertNotEqual(sampled.jacobian_determinant, 0)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        with self.assertRaisesRegex(ValueError, "inconclusive interval"):
            replace(
                record,
                source_invertible_over_box=True,
                source_inverse={"forged": True},
                rho_infinity=Q(0),
            )


class SGBLIntervalFailureModeTests(unittest.TestCase):
    def test_exact_singular_center_keeps_source_jacobian_ownership(self) -> None:
        with self.assertRaises(SGBLSourceJacobianSolveStop) as stopped:
            sgbl_enclose_adm_source_jacobian(
                SGBLIntervalBox(center=_singular(), parameter_half_width=0)
            )
        self.assertEqual(stopped.exception.reason, "singular_source_jacobian")
        self.assertEqual(stopped.exception.coefficients.jacobian_determinant, 0)

    def test_resource_limit_stops_before_a_fitted_inverse(self) -> None:
        with self.assertRaises(SGBLIntervalInconclusive) as stopped:
            sgbl_enclose_adm_source_jacobian(
                SGBLIntervalBox(
                    center=_flat(),
                    parameter_half_width=Q(1, 2**20),
                    limits=SGBLIntervalLimits(max_residual_evaluations=2),
                )
            )
        self.assertEqual(stopped.exception.reason, "resource_limit")

    def test_injected_sign_error_is_not_contained_in_the_true_enclosure(self) -> None:
        exact = sgbl_source_coefficients(_fixture_a())
        true_box = tuple(
            tuple(Interval.singleton(entry) for entry in row)
            for row in exact.jacobian
        )
        flipped = (tuple(-entry for entry in true_box[0]),) + true_box[1:]
        self.assertTrue(jacobian_contains_point(true_box, exact.jacobian))
        self.assertFalse(jacobian_contains_point(flipped, exact.jacobian))

    def test_module_does_not_import_qr_builder_health_or_old_runners(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_interval_health as owner

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
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("run_fgc_gr0_calibration", source)


if __name__ == "__main__":
    unittest.main()
