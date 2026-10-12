from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config
from recursive_horizons.fgc.modified_harmonic_constraints import (
    ACTIVATED_AREAL_RADIUS_DRR,
    ACTIVATED_H_TT_DRR,
    ACTIVATED_PHI,
    ACTIVATED_PHI_RADIAL_DERIVATIVE,
    ACCELERATION_EVALUATION_NODES,
    ACCELERATION_POLYNOMIAL_DEGREE,
    PHYSICAL_PROJECTION_ORDER,
    QIFT1_ACCELERATION_HALF_WIDTH,
    QIFT1_PARAMETER_HALF_WIDTH,
    acceleration_polynomial_certificate,
    activated_compatible_state,
    compatible_constraint_certificate,
    normal_gauge_extension_map,
    activated_principal_cone_control,
    physical_constraint_projections,
    required_con1_comp1_nonclaims,
)
from recursive_horizons.fgc.modified_harmonic_reference import modified_harmonic_gauge_constraint


class CON1COMP1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        predecessor = load_config()
        cls.fixture = predecessor.flat_fixture
        cls.datum = activated_compatible_state(cls.fixture)
        cls.configuration = {
            "qift1_record": {"artifact_id": "FGC-1-HYP1-DOM1-QIFT1"},
            "qift1_configuration": {"flat_fixture": cls.fixture},
            "reference_radius": Q(4),
            "activated_phi": ACTIVATED_PHI,
            "activated_phi_radial_derivative": ACTIVATED_PHI_RADIAL_DERIVATIVE,
            "parameter_half_width": QIFT1_PARAMETER_HALF_WIDTH,
            "acceleration_half_width": QIFT1_ACCELERATION_HALF_WIDTH,
            "formula": "H=E_tt-2sE_tr+s^2E_rr; M=E_tr-sE_rr; s=h_tr/h_rr",
            "projection_order": PHYSICAL_PROJECTION_ORDER,
            "polynomial_degree": ACCELERATION_POLYNOMIAL_DEGREE,
            "evaluation_nodes": ACCELERATION_EVALUATION_NODES,
        }

    def test_activated_datum_exact_targets_and_spatial_compatibility(self) -> None:
        datum = self.datum
        state = datum["state"]
        self.assertEqual(state.phi.value, ACTIVATED_PHI)
        self.assertEqual(state.phi.dr, ACTIVATED_PHI_RADIAL_DERIVATIVE)
        self.assertEqual(state.areal_radius.drr, ACTIVATED_AREAL_RADIUS_DRR)
        self.assertEqual(state.h_tt.drr, ACTIVATED_H_TT_DRR)
        self.assertEqual(state.h_tr.drr, 0)
        self.assertEqual(state.h_rr.drr, 0)
        gauge = modified_harmonic_gauge_constraint(
            state, reference=datum["reference"], coordinate_radius=4,
            tilde_normal_factor=datum["tilde_normal_factor"],
        )
        self.assertEqual(gauge.constraint_up, (Q(0),) * 4)
        self.assertEqual(gauge.covariant_constraint_derivative[1], (Q(0),) * 4)
        projection = physical_constraint_projections(state)
        self.assertEqual(projection["H"], 0)
        self.assertEqual(projection["M"], 0)
        self.assertEqual(projection["metric_tr"], 0)

    def test_all_unredefined_constraint_acceleration_coefficients_vanish(self) -> None:
        polynomial = acceleration_polynomial_certificate(self.datum["state"])
        self.assertEqual(polynomial["coefficient_count_per_projection"], 28)
        self.assertEqual(len(polynomial["coefficients"]), 28)
        self.assertTrue(polynomial["all_56_coefficients_zero"])
        self.assertTrue(polynomial["all_H_coefficients_zero"])
        self.assertTrue(polynomial["all_M_coefficients_zero"])

    def test_normal_gauge_map_is_exact_diagonal_and_nonsingular(self) -> None:
        datum = self.datum
        gauge = modified_harmonic_gauge_constraint(
            datum["state"], reference=datum["reference"], coordinate_radius=4,
            tilde_normal_factor=datum["tilde_normal_factor"],
        )
        normal_map = normal_gauge_extension_map(datum["state"], gauge, hat_normal_factor=datum["hat_normal_factor"])
        self.assertTrue(normal_map["diagonal_entries_nonzero"])
        self.assertTrue(normal_map["off_diagonal_entries_zero"])
        self.assertTrue(normal_map["determinant_nonzero"])

    def test_activated_principal_cone_control_is_pointwise_and_separated(self) -> None:
        control = activated_principal_cone_control(self.fixture, self.datum["state"])
        self.assertEqual(
            control["fixture_id"],
            "CON1_COMP1_activated_compatible_principal_control",
        )
        self.assertTrue(control["compatible_state_pullback_exact"])
        self.assertEqual(
            control["regulator_factor"][2],
            Q(
                1298074214621900990924908803391487,
                1298074214614817441201214220926976,
            ),
        )
        self.assertTrue(control["metric_null_discriminant_positive"])
        self.assertTrue(control["regulator_discriminant_positive"])
        self.assertTrue(control["physical_vs_regulator_resultant_nonzero"])
        self.assertTrue(control["cones_share_no_root"])
        self.assertFalse(control["regulator_cone_proportional_to_metric"])

    def test_certificate_keeps_qift_root_external_and_nonclaims_false(self) -> None:
        certificate = compatible_constraint_certificate(self.configuration)
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(certificate["local_exact_checks"]["qift_root_is_external_predecessor_not_solved_here"])
        self.assertTrue(certificate["local_exact_checks"]["normal_shift_zero"])
        self.assertTrue(
            certificate["local_exact_checks"][
                "metric_defined_all_non_normal_nabla_C_components_zero"
            ]
        )
        self.assertTrue(
            certificate["local_exact_checks"][
                "unredefined_Ett_and_Etr_zero_for_every_acceleration"
            ]
        )
        self.assertTrue(certificate["local_exact_checks"]["ref1_scalar_equations_unmodified"])
        self.assertEqual(certificate["nonclaims"], required_con1_comp1_nonclaims())
        self.assertTrue(all(value is False for value in certificate["nonclaims"].values()))

    def test_mutations_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "activated_phi"):
            activated_compatible_state(self.fixture, activated_phi=Q(1, 65536))
        with self.assertRaisesRegex(ValueError, "activated_phi_radial_derivative"):
            activated_compatible_state(self.fixture, activated_phi_radial_derivative=Q(1, 65536))
        with self.assertRaisesRegex(ValueError, "parameter_half_width"):
            activated_compatible_state(self.fixture, parameter_half_width=Q(1, 32768))

        broken = replace(self.datum["state"], areal_radius=replace(self.datum["state"].areal_radius, drr=0))
        self.assertNotEqual(physical_constraint_projections(broken)["H"], 0)

        with patch(
            "recursive_horizons.fgc.modified_harmonic_constraints.first_order_kinematic_residuals",
            return_value={
                "u_time_definition_residual": (Q(1),) + (Q(0),) * 5,
                "radial_reduction_constraint": (Q(0),) * 6,
                "mixed_partial_compatibility_residual": (Q(0),) * 6,
            },
        ):
            with self.assertRaisesRegex(ValueError, "local compatibility checks failed"):
                compatible_constraint_certificate(self.configuration)

        def corrupted_source(candidate):
            source = dict(physical_constraint_projections(candidate))
            source["source_is_unredefined_residual"] = False
            return source

        with patch(
            "recursive_horizons.fgc.modified_harmonic_constraints.physical_constraint_projections",
            side_effect=corrupted_source,
        ):
            with self.assertRaisesRegex(ValueError, "local compatibility checks failed"):
                compatible_constraint_certificate(self.configuration)

        malformed = dict(self.configuration)
        malformed["formula"] = "wrong"
        with self.assertRaisesRegex(ValueError, "formula"):
            compatible_constraint_certificate(malformed)
        malformed = dict(self.configuration)
        malformed["qift1_record"] = {"artifact_id": "wrong"}
        with self.assertRaisesRegex(ValueError, "QIFT1"):
            compatible_constraint_certificate(malformed)

    def test_polynomial_and_normal_map_mutation_probes(self) -> None:
        state = self.datum["state"]
        # A deliberately corrupted evaluator gives a nonzero a_h_tt^2
        # coefficient, proving the coefficient result is reconstructed rather
        # than hardcoded as all zero.
        def corrupted(candidate):
            H, M = (
                physical_constraint_projections(candidate)["H"],
                physical_constraint_projections(candidate)["M"],
            )
            return H + candidate.h_tt.dtt**2, M

        self.assertFalse(
            acceleration_polynomial_certificate(state, evaluator=corrupted)["all_56_coefficients_zero"]
        )

        gauge = modified_harmonic_gauge_constraint(
            state, reference=self.datum["reference"], coordinate_radius=4,
            tilde_normal_factor=4,
        )
        # Supplying a zero FGC coupling makes the normal extension map singular
        # only by leaving the frozen datum/formulation; exercise the public
        # detector directly with a state whose beta/eta branch remains valid is
        # not possible, so mutation is asserted on the returned matrix status.
        normal = normal_gauge_extension_map(state, gauge, hat_normal_factor=9)
        self.assertNotEqual(normal["determinant"], 0)
        singular = dict(normal)
        singular["diagonal_entries_nonzero"] = False
        singular["determinant_nonzero"] = False
        with patch(
            "recursive_horizons.fgc.modified_harmonic_constraints.normal_gauge_extension_map",
            return_value=singular,
        ):
            with self.assertRaisesRegex(ValueError, "local compatibility checks failed"):
                compatible_constraint_certificate(self.configuration)


if __name__ == "__main__":
    unittest.main()
