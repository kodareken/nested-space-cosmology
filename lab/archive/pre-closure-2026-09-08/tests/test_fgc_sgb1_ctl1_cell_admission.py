from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.exact_interval import Interval, interval  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_admission import (  # noqa: E402
    sgbl_parametric_source_admission,
)
from recursive_horizons.fgc.sgb1_ctl1_cell_admission import (  # noqa: E402
    CANONICAL_SLOT_NAMES,
    FORBIDDEN_HEALTH_IMPORTS,
    INCONCLUSIVE_REASONS,
    SGBLCellAdmissionGraphRecord,
    SGBLCellAdmissionInconclusive,
    SGBLCellAdmissionRecord,
    SGBLCellProductBox,
    SGBLCellSlot,
    WRAPPING_OBSTRUCTIONS,
    sgbl_cell_admission_health_gate,
    sgbl_cell_admission_slot_table,
    sgbl_cell_product_box_from_family_cell,
    sgbl_cell_product_box_from_point,
    sgbl_family_cell_whole_cell_admission,
    sgbl_flat_product_box,
    sgbl_whole_cell_admission_from_health,
    sgbl_whole_cell_source_admission,
)
from recursive_horizons.fgc.sgb1_ctl1_family_principal import (  # noqa: E402
    sgbl_family_cell_principal_box,
)
from recursive_horizons.fgc.sgb1_ctl1_initial_health import (  # noqa: E402
    SGBLExactInitialSlice,
    sgbl_continuous_initial_compactness,
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


EXPECTED_DET_A = Q(6519803496633, 6553600000)
ZERO6 = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))
FAMILY_OBSTRUCTIONS = WRAPPING_OBSTRUCTIONS | {
    "resource_limit",
    "interval_chart_domain",
    "zero_in_interval_reciprocal",
}


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


def _point_component(point: SGBLSourceInputs, slot_name: str):
    field, component = slot_name.split(".")
    attr = {
        "alpha": "alpha",
        "shift": "shift",
        "lambda": "radial_metric",
        "areal_radius": "areal_radius",
        "phi": "phi",
        "chi": "chi",
    }[field]
    return getattr(getattr(point, attr), component)


class SGBLCellAdmissionContractTests(unittest.TestCase):
    def test_slot_table_owns_every_lower_jet_and_acceleration_slot(self) -> None:
        table = sgbl_cell_admission_slot_table()
        slots = " ".join(str(entry["slot"]) for entry in table)
        for required in (
            "alpha.value",
            "lambda.dr",
            "chi.dtr",
            "areal_radius.value",
            "coordinate_radius",
            "alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt",
        ):
            self.assertIn(required, slots)
        self.assertEqual(len(CANONICAL_SLOT_NAMES), 30)
        self.assertTrue(all("owner" in entry and "formula" in entry for entry in table))
        self.assertNotIn("parameter_half_width", slots)

    def test_module_does_not_import_proto4_qr_health_or_old_runners(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_cell_admission as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        self.assertTrue(set(FORBIDDEN_HEALTH_IMPORTS).isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertFalse(hasattr(SGBLCellProductBox, "parameter_half_width"))
        self.assertNotIn("parameter_half_width", CANONICAL_SLOT_NAMES)
        self.assertNotIn("newton_residual_limit", INCONCLUSIVE_REASONS)
        self.assertIs(SGBLCellAdmissionRecord.__dataclass_params__.frozen, True)


class SGBLCellAdmissionFlatTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.box = sgbl_flat_product_box()
        cls.record = sgbl_whole_cell_source_admission(cls.box)

    def test_flat_nonuniform_product_box_proves_unique_root(self) -> None:
        box = self.box
        record = self.record
        radii = {slot.radius for slot in box.slots}
        self.assertGreater(len(radii), 1)
        self.assertFalse(box.parameter_radii_are_uniform)
        self.assertEqual(box.coverage, "flat_product")
        self.assertEqual(record.classification, "parametric_krawczyk_unique_root")
        self.assertTrue(record.unique_acceleration_root_for_every_declared_parameter_point)
        self.assertIsNone(record.inconclusive_reason)
        self.assertEqual(record.acceleration_center, ZERO6)
        self.assertEqual(sgbl_source_solve(box.center), ZERO6)
        self.assertIsNotNone(record.krawczyk)
        self.assertLess(record.krawczyk["rho_infinity_upper_bound"], 1)
        self.assertGreater(record.krawczyk["minimum_strict_componentwise_inclusion_margin"], 0)
        self.assertTrue(record.krawczyk["krawczyk_image_strictly_inside_displacement_box"])
        self.assertTrue(all(entry.contains_zero() for entry in record.residual_on_krawczyk_image))
        self.assertTrue(all(entry.contains_zero() for entry in record.residual_constant_box))
        self.assertTrue(record.interval_newton_strictly_inside_acceleration_box)
        slot_map = box.slot_map()
        self.assertEqual(len(slot_map), 30)
        self.assertTrue(all(entry["owner"] and entry["formula"] for entry in slot_map))

    def test_flat_certificate_keeps_aggregate_health_false(self) -> None:
        record = self.record
        self.assertIs(record.uses_proto4_newton_residual_limit, False)
        self.assertIs(record.finite_residual_alone_is_admission, False)
        self.assertIs(record.box_was_fitted_after_the_outcome, False)
        self.assertIs(record.source_solve_admission_qualified, False)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(record.FRZ1, False)
        self.assertIs(record.PREF1, False)
        self.assertIs(record.promoted_singleton_witness, False)
        self.assertIs(record.uniformized_nonuniform_radii, False)
        gate = sgbl_cell_admission_health_gate(record)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["source_solve_admission_qualified"], False)
        self.assertIs(gate["parametric_unique_root_certified"], True)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["promoted_singleton_witness"], False)
        with self.assertRaises(TypeError):
            replace(record, SGBL_branch_owned_and_healthy=True)


class SGBLCellAdmissionFixtureATests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.point = _fixture_a()
        cls.box = sgbl_cell_product_box_from_point(
            cls.point,
            acceleration_half_width=Q(1, 1 << 10),
            coverage="exact_singleton",
        )
        cls.record = sgbl_whole_cell_source_admission(cls.box)
        cls.uniform = sgbl_parametric_source_admission(
            SGBLIntervalBox(
                center=cls.point,
                parameter_half_width=0,
                acceleration_half_width=Q(1, 1 << 10),
            )
        )

    def test_exact_fixture_a_singleton_reproduces_the_source_jacobian(self) -> None:
        record = self.record
        self.assertTrue(self.box.is_singleton)
        self.assertEqual(self.box.coverage, "exact_singleton")
        self.assertEqual(record.exact_center_det, EXPECTED_DET_A)
        self.assertTrue(record.independent_exact_singleton_reduction)
        self.assertEqual(record.classification, "parametric_krawczyk_unique_root")
        self.assertTrue(record.unique_acceleration_root_for_every_declared_parameter_point)
        self.assertEqual(record.acceleration_center, sgbl_source_solve(self.point))
        self.assertTrue(self.uniform.unique_acceleration_root_for_every_declared_parameter_point)
        self.assertEqual(record.exact_center_jacobian, self.uniform.jacobian_record.exact_center_jacobian)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(record.FRZ1, False)
        self.assertIs(record.PREF1, False)


class SGBLCellAdmissionFamilyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        cls.health = sgbl_continuous_initial_compactness(cls.spec)
        interior = [
            cell
            for cell in cls.health.ode.cells
            if (
                cell.radius.lower > cls.spec.support_minimum
                and cell.radius.upper < cls.spec.support_maximum
            )
        ]
        cls.cells = (interior[0], interior[len(interior) // 2])
        cls.boxes = tuple(
            sgbl_cell_product_box_from_family_cell(cell, cls.spec) for cell in cls.cells
        )
        cls.records = tuple(sgbl_whole_cell_source_admission(box) for box in cls.boxes)
        cls.feeder = tuple(
            sgbl_family_cell_principal_box(cell, cls.spec) for cell in cls.cells
        )

    def test_completed_family_cells_pass_or_preserve_exact_wrapping(self) -> None:
        self.assertEqual(self.health.classification, "continuous_no_initial_trapped_sphere")
        self.assertEqual(len(self.records), 2)
        for box, record, feeder in zip(self.boxes, self.records, self.feeder, strict=True):
            radii = {slot.radius for slot in box.slots}
            self.assertGreater(len(radii), 1)
            self.assertFalse(box.parameter_radii_are_uniform)
            self.assertEqual(box.coverage, "family_cell_product")
            self.assertEqual(len(box.slot_map()), 30)
            self.assertTrue(feeder.admission.unique_acceleration_root_for_every_declared_parameter_point)
            self.assertEqual(feeder.missing_owner.slot, "dtt_over_cell_product_box")
            self.assertIs(record.promoted_singleton_witness, False)
            self.assertIs(record.uniformized_nonuniform_radii, False)
            self.assertIs(record.SGBL_branch_owned_and_healthy, False)
            self.assertIs(record.FRZ1, False)
            self.assertIs(record.PREF1, False)
            if record is self.records[0]:
                direct = sgbl_family_cell_whole_cell_admission(self.cells[0], self.spec)
                self.assertEqual(direct.classification, record.classification)
                self.assertEqual(
                    direct.unique_acceleration_root_for_every_declared_parameter_point,
                    record.unique_acceleration_root_for_every_declared_parameter_point,
                )
            if record.unique_acceleration_root_for_every_declared_parameter_point:
                self.assertEqual(record.classification, "parametric_krawczyk_unique_root")
                self.assertIsNotNone(record.krawczyk)
                self.assertLess(record.krawczyk["rho_infinity_upper_bound"], 1)
                self.assertGreater(
                    record.krawczyk["minimum_strict_componentwise_inclusion_margin"], 0
                )
                self.assertFalse(record.box.is_singleton)
            else:
                self.assertEqual(record.classification, "interval_inconclusive")
                self.assertIn(record.inconclusive_reason, FAMILY_OBSTRUCTIONS)
                self.assertFalse(record.unique_acceleration_root_for_every_declared_parameter_point)
                self.assertIsNone(record.krawczyk)
                self.assertNotEqual(record.inconclusive_reason, "newton_residual_limit")

    def test_aggregate_health_stays_false_on_completed_cells(self) -> None:
        graph = sgbl_whole_cell_admission_from_health(
            self.health,
            spec=self.spec,
            max_cells=1,
        )
        self.assertEqual(graph.classification, "whole_cell_family_records")
        self.assertEqual(len(graph.cell_records), 1)
        gate = sgbl_cell_admission_health_gate(graph)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["execution_authorized"], False)
        self.assertIs(gate["promoted_singleton_witness"], False)
        self.assertIs(graph.SGBL_branch_owned_and_healthy, False)

    def test_correlated_uniformization_of_family_radii_is_refused(self) -> None:
        box = self.boxes[0]
        max_radius = box.max_slot_radius
        self.assertGreater(max_radius, 0)
        tampered = tuple(
            SGBLCellSlot(
                name=slot.name,
                interval=interval(
                    _point_component(box.center, slot.name) - max_radius,
                    _point_component(box.center, slot.name) + max_radius,
                ),
                owner=slot.owner,
                formula=slot.formula,
            )
            for slot in box.slots
        )
        with self.assertRaises(SGBLCellAdmissionInconclusive) as stopped:
            replace(box, slots=tampered)
        self.assertEqual(stopped.exception.reason, "correlated_uniformization")

    def test_singleton_witness_promotion_is_refused(self) -> None:
        box = self.boxes[0]
        self.assertFalse(self.cells[0].radius.is_singleton())
        singletons = tuple(
            SGBLCellSlot(
                name=slot.name,
                interval=Interval.singleton(_point_component(box.center, slot.name)),
                owner=slot.owner,
                formula=slot.formula,
            )
            for slot in box.slots
        )
        with self.assertRaises(SGBLCellAdmissionInconclusive) as stopped:
            replace(box, slots=singletons)
        self.assertEqual(stopped.exception.reason, "promoted_singleton_witness")


class SGBLCellAdmissionIncompleteTests(unittest.TestCase):
    def test_nominal_chi_three_refuses_the_incomplete_graph(self) -> None:
        spec = SGBLExactInitialSlice()
        self.assertEqual(spec.chi_amplitude, 3)
        health = sgbl_continuous_initial_compactness(spec)
        self.assertEqual(health.classification, "interval_inconclusive")
        self.assertEqual(health.inconclusive_reason, "picard_strict_self_map_failed")
        record = sgbl_whole_cell_admission_from_health(health, spec=spec)
        self.assertEqual(record.classification, "family_geometry_incomplete")
        self.assertEqual(record.inconclusive_reason, "family_geometry_incomplete")
        self.assertEqual(record.cell_records, ())
        self.assertEqual(
            record.missing_theorem,
            "complete_validated_constraint_ODE_graph_covering_the_compact_support",
        )
        gate = sgbl_cell_admission_health_gate(record)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertTrue(gate["family_geometry_incomplete"])
        self.assertIs(gate["parametric_unique_root_certified"], False)
        with self.assertRaises(ValueError):
            SGBLCellAdmissionGraphRecord(
                compactness=health,
                cell_records=(sgbl_whole_cell_source_admission(sgbl_flat_product_box()),),
                classification="family_geometry_incomplete",
                inconclusive_reason="family_geometry_incomplete",
                missing_theorem="complete_validated_constraint_ODE_graph_covering_the_compact_support",
                theorem="probe",
            )


class SGBLCellAdmissionAttackTests(unittest.TestCase):
    def test_omission_of_a_required_slot_is_refused(self) -> None:
        box = sgbl_flat_product_box()
        with self.assertRaises(SGBLCellAdmissionInconclusive) as stopped:
            replace(box, slots=box.slots[:-1])
        self.assertEqual(stopped.exception.reason, "omitted_slot")
        self.assertIn("chi.drr", stopped.exception.payload["omitted"])

    def test_uniform_half_width_alias_is_refused(self) -> None:
        point = sgbl_flat_product_box().center
        with self.assertRaises(SGBLCellAdmissionInconclusive) as stopped:
            sgbl_cell_product_box_from_point(
                point,
                acceleration_half_width=Q(1, 8),
                uniform_half_width=Q(1, 8),
            )
        self.assertEqual(stopped.exception.reason, "uniformization_alias")
        with self.assertRaises(TypeError):
            sgbl_cell_product_box_from_point(
                point,
                slot_half_widths={"lambda.value": 0.5},
                acceleration_half_width=Q(1, 8),
            )
        with self.assertRaises(TypeError):
            sgbl_cell_product_box_from_point(
                point,
                acceleration_half_width=True,
            )

    def test_zero_acceleration_width_is_a_typed_interior_stop(self) -> None:
        with self.assertRaises(SGBLCellAdmissionInconclusive) as stopped:
            sgbl_flat_product_box(acceleration_half_width=0)
        self.assertEqual(stopped.exception.reason, "acceleration_box_not_interior")

    def test_exact_singular_center_keeps_source_jacobian_ownership(self) -> None:
        box = sgbl_cell_product_box_from_point(
            _singular(),
            acceleration_half_width=Q(1, 8),
            coverage="exact_singleton",
        )
        with self.assertRaises(SGBLSourceJacobianSolveStop) as stopped:
            sgbl_whole_cell_source_admission(box)
        self.assertEqual(stopped.exception.reason, "singular_source_jacobian")

    def test_zero_denominator_is_a_typed_reciprocal_stop(self) -> None:
        box = sgbl_flat_product_box(
            slot_half_widths={"lambda.value": 2},
            acceleration_half_width=Q(1, 8),
        )
        with self.assertRaises(SGBLCellAdmissionInconclusive) as stopped:
            sgbl_whole_cell_source_admission(box)
        self.assertIn(
            stopped.exception.reason,
            {"zero_in_interval_reciprocal", "interval_chart_domain"},
        )

    def test_resource_limit_stops_before_a_fitted_box(self) -> None:
        widths = {"lambda.value": Q(1, 1 << 22), "areal_radius.value": Q(1, 1 << 24)}
        box = sgbl_flat_product_box(
            slot_half_widths=widths,
            acceleration_half_width=Q(1, 1 << 10),
            limits=SGBLIntervalLimits(max_residual_evaluations=2),
        )
        with self.assertRaises(SGBLCellAdmissionInconclusive) as stopped:
            sgbl_whole_cell_source_admission(box)
        self.assertEqual(stopped.exception.reason, "resource_limit")
        self.assertEqual(
            box.slot_interval("lambda", "value").radius(),
            Q(1, 1 << 22),
        )
        self.assertEqual(
            box.slot_interval("areal_radius", "value").radius(),
            Q(1, 1 << 24),
        )

    def test_wrong_declared_center_does_not_refit_the_box(self) -> None:
        box = sgbl_flat_product_box()
        original_slots = box.slots
        original_acceleration = box.acceleration_box
        shifted = (Q(1), Q(0), Q(0), Q(0), Q(0), Q(0))
        with self.assertRaises(SGBLCellAdmissionInconclusive) as stopped:
            sgbl_whole_cell_source_admission(box, acceleration_center=shifted)
        self.assertEqual(stopped.exception.reason, "acceleration_box_not_interior")
        self.assertEqual(box.slots, original_slots)
        self.assertEqual(box.acceleration_box, original_acceleration)
        self.assertNotIn(stopped.exception.reason, {None, "newton_residual_limit"})

    def test_forged_uniqueness_without_krawczyk_is_refused(self) -> None:
        record = sgbl_whole_cell_source_admission(sgbl_flat_product_box())
        with self.assertRaisesRegex(ValueError, "Krawczyk"):
            replace(
                record,
                unique_acceleration_root_for_every_declared_parameter_point=True,
                classification="parametric_krawczyk_unique_root",
                krawczyk={"forged": True},
                inconclusive_reason=None,
            )


if __name__ == "__main__":
    unittest.main()
