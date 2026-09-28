from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.sgb1_ctl1_admission import (  # noqa: E402
    INCONCLUSIVE_REASONS,
    SGBLAdmissionInconclusive,
    SGBLSourceAdmissionRecord,
    sgbl_parametric_source_admission,
    sgbl_source_admission_health_gate,
)
from recursive_horizons.fgc.sgb1_ctl1_interval_health import (  # noqa: E402
    SGBLIntervalBox,
    SGBLIntervalLimits,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    SGBLSourceInputs,
    SGBLSourceJacobianSolveStop,
    sgbl_source_solve,
)


ZERO6 = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))


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


class SGBLAdmissionInputTests(unittest.TestCase):
    def test_zero_acceleration_width_is_a_typed_interior_stop(self) -> None:
        box = SGBLIntervalBox(center=_flat(), parameter_half_width=0, acceleration_half_width=0)
        with self.assertRaises(SGBLAdmissionInconclusive) as stopped:
            sgbl_parametric_source_admission(box)
        self.assertEqual(stopped.exception.reason, "acceleration_box_not_interior")

    def test_exact_singular_center_keeps_source_jacobian_ownership(self) -> None:
        box = SGBLIntervalBox(
            center=_singular(),
            parameter_half_width=0,
            acceleration_half_width=Q(1, 8),
        )
        with self.assertRaises(SGBLSourceJacobianSolveStop) as stopped:
            sgbl_parametric_source_admission(box)
        self.assertEqual(stopped.exception.reason, "singular_source_jacobian")

    def test_resource_limit_stops_before_a_fitted_box(self) -> None:
        declared = Q(1, 2**20)
        box = SGBLIntervalBox(
            center=_flat(),
            parameter_half_width=declared,
            acceleration_half_width=Q(1, 2**10),
            limits=SGBLIntervalLimits(max_residual_evaluations=2),
        )
        with self.assertRaises(SGBLAdmissionInconclusive) as stopped:
            sgbl_parametric_source_admission(box)
        self.assertEqual(stopped.exception.reason, "resource_limit")
        self.assertEqual(box.parameter_half_width, declared)
        self.assertEqual(box.acceleration_half_width, Q(1, 2**10))


class SGBLAdmissionCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.point = _flat()
        cls.parameter_width = Q(1, 2**20)
        cls.acceleration_width = Q(1, 2**10)
        cls.box = SGBLIntervalBox(
            center=cls.point,
            parameter_half_width=cls.parameter_width,
            acceleration_half_width=cls.acceleration_width,
        )
        cls.record = sgbl_parametric_source_admission(cls.box)

    def test_declared_flat_box_proves_unique_root_by_krawczyk_and_residual(self) -> None:
        record = self.record
        self.assertEqual(record.classification, "parametric_krawczyk_unique_root")
        self.assertTrue(record.unique_acceleration_root_for_every_declared_parameter_point)
        self.assertIsNone(record.inconclusive_reason)
        self.assertEqual(record.acceleration_center, ZERO6)
        self.assertEqual(sgbl_source_solve(self.point), ZERO6)
        self.assertIsNotNone(record.krawczyk)
        self.assertLess(record.krawczyk["rho_infinity_upper_bound"], 1)
        self.assertGreater(record.krawczyk["minimum_strict_componentwise_inclusion_margin"], 0)
        self.assertTrue(
            record.krawczyk["krawczyk_image_strictly_inside_displacement_box"]
        )
        self.assertTrue(all(entry.contains_zero() for entry in record.residual_on_krawczyk_image))
        self.assertTrue(record.interval_newton_strictly_inside_acceleration_box)
        self.assertEqual(record.box.parameter_half_width, self.parameter_width)
        self.assertEqual(record.box.acceleration_half_width, self.acceleration_width)

    def test_certificate_is_not_a_fitted_threshold_or_proto4_floor(self) -> None:
        record = self.record
        self.assertIs(record.uses_proto4_newton_residual_limit, False)
        self.assertIs(record.finite_residual_alone_is_admission, False)
        self.assertIs(record.box_was_fitted_after_the_outcome, False)
        self.assertIs(record.source_solve_admission_qualified, False)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertNotIn("newton_residual_limit", INCONCLUSIVE_REASONS)
        self.assertNotIn("newton_residual_limit", record.inconclusive_reason or "")
        gate = sgbl_source_admission_health_gate(record)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["source_solve_admission_qualified"], False)
        self.assertIs(gate["parametric_unique_root_certified"], True)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["uses_proto4_newton_residual_limit"], False)
        with self.assertRaises(TypeError):
            replace(record, SGBL_branch_owned_and_healthy=True)

    def test_wrong_declared_center_does_not_refit_the_box(self) -> None:
        shifted = (Q(1), Q(0), Q(0), Q(0), Q(0), Q(0))
        record = sgbl_parametric_source_admission(
            self.box, acceleration_center=shifted
        )
        self.assertEqual(record.classification, "interval_inconclusive")
        self.assertFalse(record.unique_acceleration_root_for_every_declared_parameter_point)
        self.assertEqual(record.box.acceleration_half_width, self.acceleration_width)
        self.assertEqual(record.acceleration_center, shifted)
        self.assertNotIn(record.inconclusive_reason, {None, "newton_residual_limit"})


class SGBLAdmissionWrappingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.width = Q(1, 16)
        cls.box = SGBLIntervalBox(
            center=_flat(),
            parameter_half_width=cls.width,
            acceleration_half_width=Q(1, 8),
        )
        cls.record = sgbl_parametric_source_admission(cls.box)

    def test_wide_wrapping_is_inconclusive_and_does_not_fit_a_smaller_box(self) -> None:
        record = self.record
        self.assertEqual(record.classification, "interval_inconclusive")
        self.assertFalse(record.unique_acceleration_root_for_every_declared_parameter_point)
        self.assertEqual(record.inconclusive_reason, "neumann_rho_not_below_one")
        self.assertEqual(record.box.parameter_half_width, self.width)
        self.assertIs(record.source_solve_admission_qualified, False)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        with self.assertRaisesRegex(ValueError, "Krawczyk payload"):
            replace(
                record,
                unique_acceleration_root_for_every_declared_parameter_point=True,
                classification="parametric_krawczyk_unique_root",
                krawczyk={"forged": True},
                inconclusive_reason=None,
            )


class SGBLAdmissionImportTests(unittest.TestCase):
    def test_module_does_not_import_proto4_qr_health_or_old_runners(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_admission as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        forbidden = {
            "newton_residual_limit",
            "branch_stop_applicability",
            "background_from_spherical_state",
            "weak_coupling_health_certificate",
            "canonical_health_monitor_values",
            "solve_accelerations",
            "FGCQRActionParameters",
            "solve_initial_data",
            "proto1_compactness_upper_bound",
            "run_fgc_gr0_calibration_v1",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertIs(SGBLSourceAdmissionRecord.__dataclass_params__.frozen, True)


if __name__ == "__main__":
    unittest.main()
