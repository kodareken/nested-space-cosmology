from __future__ import annotations

import ast
from pathlib import Path
import sys
import unittest

from fractions import Fraction as Q

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.exact_interval import Interval  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_cone import (  # noqa: E402
    ORTHONORMAL_MINKOWSKI_CHART,
    _contains_float_matrix,
    sgbl_enclose_principal_coefficients,
    sgbl_enclose_principal_cone,
    sgbl_nonflat_source_principal_box,
)
from recursive_horizons.fgc.sgb1_ctl1_family_principal import (  # noqa: E402
    BUMP_XXX_X_DERIVATIVE_BOUND,
    FORBIDDEN_HEALTH_IMPORTS,
    SGBLFamilyPrincipalRecord,
    SGBLFamilyPrincipalStop,
    sgbl_compact_bump_third_derivative_enclosure,
    sgbl_complete_interval_family_fields,
    sgbl_family_affine_jets,
    sgbl_family_cell_principal_box,
    sgbl_family_principal_from_health,
    sgbl_family_principal_health_gate,
    sgbl_family_principal_slot_table,
    sgbl_minkowski_buffer_principal_box,
    sgbl_nonflat_fixture_a_interval_principal_box,
)
from recursive_horizons.fgc.sgb1_ctl1_initial_health import (  # noqa: E402
    SGBLExactInitialSlice,
    sgbl_compact_bump_enclosure,
    sgbl_continuous_initial_compactness,
)
from recursive_horizons.fgc.sgb1_ctl1_principal import (  # noqa: E402
    sgbl_covariant_principal_background,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    SGBLSourceInputs,
    sgbl_source_solve,
    sgbl_source_state,
)
from recursive_horizons.fgc.evolution.covariant_principal_health import (  # noqa: E402
    principal_coefficient_tensors,
)


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


def _contains_float_tensor4(enclosure, array: np.ndarray) -> bool:
    from recursive_horizons.fgc.sgb1_ctl1_cone import BINARY64_UNIT, DEFAULT_ROUND_UNITS

    for first in range(4):
        for second in range(4):
            for third in range(4):
                for fourth in range(4):
                    dyadic = Q(*float(array[first, second, third, fourth]).as_integer_ratio())
                    entry = enclosure[first][second][third][fourth]
                    scale = max(Q(1), abs(entry.lower), abs(entry.upper), abs(dyadic))
                    slack = DEFAULT_ROUND_UNITS * BINARY64_UNIT * scale
                    if dyadic < entry.lower - slack or dyadic > entry.upper + slack:
                        return False
    return True


class SGBLFamilyPrincipalContractTests(unittest.TestCase):
    def test_slot_table_owns_every_required_jet_and_curvature_slot(self) -> None:
        table = sgbl_family_principal_slot_table()
        slots = " ".join(entry["slot"] for entry in table)
        for required in (
            "chi",
            "chi_r",
            "chi_rr",
            "chi_pi",
            "chi_pi_r / chi_tr",
            "lambda_r",
            "k_r",
            "lambda_rr, k_rr",
            "phi_rrr",
            "alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt",
            "Riemann_coord",
            "Hess(phi)_coord",
            "orthonormal frame",
            "SGBLPrincipalBackgroundBox",
        ):
            self.assertIn(required, slots)
        self.assertTrue(all("owner" in entry and "formula" in entry for entry in table))
        self.assertGreaterEqual(len(table), 20)

    def test_float_inputs_are_refused_as_enclosures(self) -> None:
        from recursive_horizons.fgc.sgb1_ctl1_family_principal import _as_interval

        with self.assertRaises(SGBLFamilyPrincipalStop) as stopped:
            _as_interval(0.5)
        self.assertEqual(stopped.exception.reason, "float_enclosure_refused")
        fields = sgbl_complete_interval_family_fields(
            Interval.singleton(12),
            SGBLExactInitialSlice(chi_amplitude=Q(1, 8)),
        )
        self.assertIsInstance(fields["chi"], Interval)

    def test_module_does_not_import_hyp2_or_qr_health_bits(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_family_principal as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        self.assertTrue(set(FORBIDDEN_HEALTH_IMPORTS).isdisjoint(imported))
        self.assertNotIn("WeakCouplingThresholds", imported)
        self.assertNotIn("weak_coupling_health_certificate", imported)
        self.assertNotIn("background_from_spherical_state", imported)
        source = Path(owner.__file__).read_text()
        self.assertNotIn("run_fgc_gr0_calibration", source)
        self.assertNotIn("sample_count=512", source)


class SGBLFamilyPrincipalProfileTests(unittest.TestCase):
    def test_third_bump_derivative_vanishes_exactly_at_the_peak(self) -> None:
        radius = Interval.singleton(12)
        third = sgbl_compact_bump_third_derivative_enclosure(
            radius, center=Q(12), half_width=Q(2)
        )
        self.assertEqual(third, Interval.singleton(0))
        value, first, second = sgbl_compact_bump_enclosure(
            radius, center=Q(12), half_width=Q(2)
        )
        self.assertEqual(first, Interval.singleton(0))
        self.assertLess(second.upper, 0)

    def test_outside_support_the_third_derivative_is_exactly_zero(self) -> None:
        third = sgbl_compact_bump_third_derivative_enclosure(
            Interval(0, 1), center=Q(12), half_width=Q(2)
        )
        self.assertEqual(third, Interval.singleton(0))

    def test_complete_fields_extend_chi_two_jet_and_keep_phi_momenta_zero(self) -> None:
        spec = SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        fields = sgbl_complete_interval_family_fields(Interval.singleton(Q(97, 8)), spec)
        for name in (
            "chi",
            "chi_r",
            "chi_rr",
            "chi_pi",
            "chi_pi_r",
            "chi_tr",
            "phi_rrr",
        ):
            self.assertIn(name, fields)
        self.assertEqual(fields["chi_pi_r"], fields["chi_tr"])
        self.assertEqual(fields["phi_pi"], Interval.singleton(0))
        self.assertEqual(fields["phi_pi_r"], Interval.singleton(0))
        self.assertGreater(fields["chi"].abs_upper(), 0)
        self.assertLessEqual(fields["phi_rrr"].abs_upper(), BUMP_XXX_X_DERIVATIVE_BOUND / 8)

    def test_affine_jets_enclose_lambda_rr_without_substituting_zero(self) -> None:
        spec = SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        jets = sgbl_family_affine_jets(
            Interval.singleton(Q(97, 8)),
            Interval.singleton(1),
            Interval.singleton(0),
            spec,
        )
        self.assertIn("lambda_r", jets)
        self.assertIn("k_r", jets)
        self.assertIn("lambda_rr", jets)
        self.assertIn("k_rr", jets)
        self.assertFalse(
            jets["lambda_rr"].is_singleton()
            and jets["lambda_rr"].lower == 0
            and jets["chi"].abs_upper() == 0
        )


class SGBLFamilyPrincipalBufferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice()
        cls.box = sgbl_minkowski_buffer_principal_box(cls.spec)
        cls.record = sgbl_enclose_principal_cone(cls.box)

    def test_exact_buffer_is_curvature_free_and_proves_the_cone(self) -> None:
        self.assertTrue(self.box.curvature_free)
        self.assertTrue(self.box.is_singleton)
        self.assertEqual(self.box.chart, ORTHONORMAL_MINKOWSKI_CHART)
        self.assertEqual(self.record.classification, "continuum_cone_proved")
        self.assertTrue(self.record.real_complete_basis)
        self.assertTrue(self.record.positive_symmetrizer)
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertEqual(self.box.effective_planck_prime, Interval.singleton(0))


class SGBLFamilyPrincipalNonflatFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.box = sgbl_nonflat_fixture_a_interval_principal_box()
        cls.state = sgbl_source_state(_fixture_a(), sgbl_source_solve(_fixture_a()))
        cls.background = sgbl_covariant_principal_background(cls.state)
        cls.record = sgbl_enclose_principal_cone(cls.box)

    def test_exact_interval_box_contains_the_numpy_adapter_and_keeps_the_typed_cone(self) -> None:
        self.assertFalse(self.box.curvature_free)
        self.assertEqual(self.box.chart, ORTHONORMAL_MINKOWSKI_CHART)
        self.assertTrue(
            _contains_float_tensor4(self.box.riemann_lower, self.background.riemann_lower)
        )
        numpy_hess = self.background.hessian_gb_lower / self.background.gb_coupling_prime
        from recursive_horizons.fgc.sgb1_ctl1_cone import BINARY64_UNIT, DEFAULT_ROUND_UNITS

        for first in range(4):
            for second in range(4):
                dyadic = Q(*float(numpy_hess[first, second]).as_integer_ratio())
                entry = self.box.hessian_phi_lower[first][second]
                scale = max(Q(1), abs(entry.lower), abs(entry.upper), abs(dyadic))
                slack = DEFAULT_ROUND_UNITS * BINARY64_UNIT * scale
                self.assertGreaterEqual(dyadic, entry.lower - slack)
                self.assertLessEqual(dyadic, entry.upper + slack)
        enclosure = sgbl_enclose_principal_coefficients(self.box)
        numpy_tensors = principal_coefficient_tensors(self.background)
        self.assertTrue(
            _contains_float_matrix(enclosure.time_time, numpy_tensors.time_time)
        )
        self.assertEqual(self.record.classification, "interval_inconclusive")
        self.assertEqual(
            self.record.inconclusive_reason, "resolvent_margin_not_strictly_positive"
        )
        adapter_box = sgbl_nonflat_source_principal_box()
        self.assertFalse(adapter_box.curvature_free)
        self.assertEqual(adapter_box.chart, self.box.chart)
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertGreater(self.record.companion_deformation_upper, 1)


class SGBLFamilyPrincipalSmallAmplitudeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice(chi_amplitude=Q(1, 8))
        cls.health = sgbl_continuous_initial_compactness(cls.spec)
        interior = [
            cell
            for cell in cls.health.ode.cells
            if (cell.radius.lower > cls.spec.support_minimum and cell.radius.upper < cls.spec.support_maximum)
        ]
        cls.cells = (interior[0], interior[len(interior) // 2])
        cls.cell_records = tuple(
            sgbl_family_cell_principal_box(cell, cls.spec) for cell in cls.cells
        )

    def test_completed_family_cells_produce_valid_principal_boxes(self) -> None:
        self.assertEqual(self.health.classification, "continuous_no_initial_trapped_sphere")
        self.assertIsNone(self.health.ode.obstruction)
        self.assertEqual(self.health.ode.cells[0].radius.lower, self.spec.support_minimum)
        self.assertEqual(self.health.ode.cells[-1].radius.upper, self.spec.support_maximum)
        self.assertEqual(len(self.cell_records), 2)
        for record in self.cell_records:
            self.assertEqual(record.classification, "principal_box")
            self.assertIsNotNone(record.box)
            self.assertEqual(record.box.chart, ORTHONORMAL_MINKOWSKI_CHART)
            self.assertEqual(record.box.effective_planck_prime, Interval.singleton(0))
            self.assertTrue(record.admission.unique_acceleration_root_for_every_declared_parameter_point)
            self.assertIn("chi_tr", record.jets.fields)
            self.assertIn("lambda_rr", record.jets.fields)
            self.assertIn("k_rr", record.jets.fields)
            self.assertIsNotNone(record.missing_owner)
            self.assertEqual(record.missing_owner.slot, "dtt_over_cell_product_box")
            self.assertNotEqual(record.missing_owner.reason, "substituted_zero")

    def test_aggregate_health_stays_false_on_completed_cells(self) -> None:
        record = sgbl_family_principal_from_health(
            self.health,
            spec=self.spec,
            max_cells=1,
            include_buffer=True,
        )
        self.assertEqual(record.classification, "family_principal_boxes")
        self.assertIsNotNone(record.buffer_box)
        self.assertTrue(record.buffer_box.curvature_free)
        self.assertGreaterEqual(
            sum(1 for cell in record.cell_records if cell.box is not None), 1
        )
        gate = sgbl_family_principal_health_gate(record)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["execution_authorized"], False)
        self.assertIs(gate["holdout_authorized"], False)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)


class SGBLFamilyPrincipalIncompleteTests(unittest.TestCase):
    def test_nominal_chi_three_returns_one_typed_incomplete_without_fabricating_exterior(self) -> None:
        spec = SGBLExactInitialSlice()
        self.assertEqual(spec.chi_amplitude, 3)
        health = sgbl_continuous_initial_compactness(spec)
        self.assertEqual(health.classification, "interval_inconclusive")
        self.assertEqual(health.inconclusive_reason, "picard_strict_self_map_failed")
        record = sgbl_family_principal_from_health(health, spec=spec)
        self.assertEqual(record.classification, "family_geometry_incomplete")
        self.assertEqual(record.inconclusive_reason, "family_geometry_incomplete")
        self.assertIsNone(record.buffer_box)
        self.assertIsNone(record.exterior_box)
        self.assertEqual(record.cell_records, ())
        self.assertEqual(len(record.missing_owners), 1)
        self.assertEqual(record.missing_owners[0].reason, "family_geometry_incomplete")
        self.assertEqual(
            record.missing_theorem,
            "complete_validated_constraint_ODE_graph_covering_the_compact_support",
        )
        gate = sgbl_family_principal_health_gate(record)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertTrue(gate["family_geometry_incomplete"])
        self.assertFalse(gate["exterior_box_emitted"])
        with self.assertRaises(ValueError):
            SGBLFamilyPrincipalRecord(
                slice=spec,
                compactness=health,
                buffer_box=sgbl_minkowski_buffer_principal_box(spec),
                cell_records=(),
                exterior_box=sgbl_minkowski_buffer_principal_box(spec),
                classification="family_geometry_incomplete",
                inconclusive_reason="family_geometry_incomplete",
                missing_theorem="complete_validated_constraint_ODE_graph_covering_the_compact_support",
                missing_owners=record.missing_owners,
                theorem="probe",
            )


if __name__ == "__main__":
    unittest.main()
