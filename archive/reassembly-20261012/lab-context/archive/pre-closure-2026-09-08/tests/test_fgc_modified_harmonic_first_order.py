from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.reproduce_fgc_hyp1_mhg_implicit import load_config  # noqa: E402
from scripts.reproduce_fgc_hyp1_reduction import (  # noqa: E402
    load_config as load_reduction_config,
)
from recursive_horizons.fgc.exact_linear_algebra import matrix_mul  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_first_order import (  # noqa: E402
    FO1_EQUATION_ORDER,
    FO1_IMPLICIT_NONACCELERATION_GROUPS,
    FO1_STATE_ORDER,
    FirstOrderFieldJet,
    first_order_kinematic_residuals,
    first_order_state_from_spherical_state,
    modified_harmonic_first_order_certificate,
    modified_harmonic_first_order_dae_residuals,
    required_fo1_rc1_nonclaims,
    spherical_state_from_first_order_state,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    state_from_generalized_adm_pg_fixture,
)


class FGCModifiedHarmonicFirstOrderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.implicit = load_config()
        reduction = load_reduction_config()
        cls.activated_fixture = next(
            fixture
            for fixture in reduction.configuration["fixtures"]
            if fixture["fixture_id"] == "FGCQR_activated_generic"
        )
        cls.reference = flat_spherical_annulus_reference(
            radial_domain_minimum=cls.implicit.radial_domain_minimum
        )
        cls.adapter = {
            "flat_fixture": cls.implicit.fixture,
            "activated_fixture": cls.activated_fixture,
            "reference": cls.reference,
            "flat_coordinate_radius": cls.implicit.coordinate_radius,
            "activated_coordinate_radius": Q(9, 2),
            "tilde_normal_factor": cls.implicit.tilde_normal_factor,
            "hat_normal_factor": cls.implicit.hat_normal_factor,
        }
        cls.certificate = modified_harmonic_first_order_certificate(cls.adapter)

    def test_exact_field_jet_validation_and_lossless_roundtrip(self) -> None:
        with self.assertRaisesRegex(TypeError, "exact rational"):
            FirstOrderFieldJet(  # type: ignore[arg-type]
                0.5, Q(0), Q(0), Q(0), Q(0), Q(0), Q(0), Q(0), Q(0)
            )
        with self.assertRaisesRegex(TypeError, "exact rational"):
            FirstOrderFieldJet(  # type: ignore[arg-type]
                True, Q(0), Q(0), Q(0), Q(0), Q(0), Q(0), Q(0), Q(0)
            )

        for fixture in (self.implicit.fixture, self.activated_fixture):
            state = state_from_generalized_adm_pg_fixture(fixture)
            first_order = first_order_state_from_spherical_state(state)
            self.assertEqual(
                spherical_state_from_first_order_state(first_order), state
            )
            kinematic = first_order_kinematic_residuals(first_order)
            self.assertTrue(
                all(
                    all(entry == 0 for entry in kinematic[name])
                    for name in (
                        "u_time_definition_residual",
                        "radial_reduction_constraint",
                        "mixed_partial_compatibility_residual",
                    )
                )
            )

    def test_kinematic_violations_remain_visible_and_fail_closed(self) -> None:
        state = state_from_generalized_adm_pg_fixture(self.implicit.fixture)
        first_order = first_order_state_from_spherical_state(state)
        changed_h_tt = replace(first_order.h_tt, q_t=Q(1))
        inconsistent = replace(first_order, h_tt=changed_h_tt)
        kinematic = first_order_kinematic_residuals(inconsistent)
        self.assertEqual(
            kinematic["mixed_partial_compatibility_residual"][0], Q(1)
        )
        self.assertEqual(
            kinematic[
                "radial_reduction_constraint_time_derivative_on_definition_shell"
            ][0],
            Q(1),
        )
        with self.assertRaisesRegex(ValueError, "kinematic equations"):
            spherical_state_from_first_order_state(inconsistent)

        dae = modified_harmonic_first_order_dae_residuals(
            inconsistent,
            reference=self.reference,
            coordinate_radius=self.implicit.coordinate_radius,
            tilde_normal_factor=self.implicit.tilde_normal_factor,
            hat_normal_factor=self.implicit.hat_normal_factor,
        )
        self.assertEqual(dae["dae_residual_vector"][-6], Q(1))
        self.assertFalse(dae["nonlinear_p_t_solve_performed"])

    def test_exact_dae_lift_flat_linearization_and_activated_control(self) -> None:
        certificate = self.certificate
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(certificate["exact_local_first_order_dae_lift_derived"])
        self.assertTrue(
            certificate["complete_flat_root_linearized_first_order_map_derived"]
        )
        self.assertTrue(
            certificate["kinematic_radial_reduction_constraint_identity_derived"]
        )
        self.assertEqual(len(FO1_STATE_ORDER), 18)
        self.assertEqual(len(FO1_EQUATION_ORDER), 18)

        linearized = certificate["flat_root_linearized_first_order_system"]
        self.assertEqual(
            linearized["acceleration_jacobian_determinant"], Q(-165888)
        )
        self.assertEqual(linearized["acceleration_jacobian_rank"], 6)
        self.assertTrue(linearized["implicit_branch_derivative_identity_exact"])
        self.assertTrue(linearized["p_q_principal_matrix_matches_mhg1_exact"])
        self.assertEqual(len(linearized["linearized_radial_principal_matrix"]), 18)
        self.assertTrue(
            all(
                len(row) == 18
                for row in linearized["linearized_radial_principal_matrix"]
            )
        )

        activated = certificate["controls"]["activated_off_shell"]
        self.assertEqual(activated["solution_status"], "off_shell_control_not_a_solution")
        self.assertTrue(any(activated["complete_ref1_residual"]))
        self.assertEqual(activated["acceleration_block_rank"], 6)

    def test_complete_argument_jacobian_has_independent_assembly_checks(self) -> None:
        jacobian = self.certificate["controls"]["flat"]["argument_jacobian"]
        self.assertEqual(len(jacobian["jacobian"]), 6)
        self.assertTrue(all(len(row) == 36 for row in jacobian["jacobian"]))
        self.assertTrue(jacobian["seeded_primals_reproduce_base_residual"])
        self.assertTrue(jacobian["weighted_direction_superposition_exact"])

        linearized = self.certificate["flat_root_linearized_first_order_system"]
        branch = [list(row) for row in linearized["implicit_branch_derivative"]]
        branch[0][0] += 1
        acceleration = linearized["acceleration_jacobian"]
        blocks = jacobian["blocks"]
        nonacceleration = tuple(
            tuple(
                entry
                for group in FO1_IMPLICIT_NONACCELERATION_GROUPS
                for entry in blocks[group][row]
            )
            for row in range(6)
        )
        mutated_product = matrix_mul(
            acceleration, tuple(tuple(row) for row in branch)
        )
        self.assertTrue(
            any(
                mutated_product[row][column] + nonacceleration[row][column] != 0
                for row in range(6)
                for column in range(30)
            )
        )

    def test_nonclaims_and_configuration_adapter_fail_closed(self) -> None:
        self.assertEqual(
            self.certificate["nonclaims"], required_fo1_rc1_nonclaims()
        )
        self.assertTrue(
            all(value is False for value in required_fo1_rc1_nonclaims().values())
        )
        with self.assertRaisesRegex(ValueError, "unexpected keys"):
            modified_harmonic_first_order_certificate(
                {**self.adapter, "extra": True}
            )
        with self.assertRaisesRegex(ValueError, "inside the reference annulus"):
            modified_harmonic_first_order_certificate(
                {**self.adapter, "flat_coordinate_radius": Q(1, 2)}
            )


if __name__ == "__main__":
    unittest.main()
