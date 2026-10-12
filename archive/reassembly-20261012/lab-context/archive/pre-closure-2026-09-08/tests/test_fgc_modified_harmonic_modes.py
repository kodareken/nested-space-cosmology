from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path[:0] = [str(Path(__file__).resolve().parents[1]), str(Path(__file__).resolve().parents[1] / "src")]

from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config
from recursive_horizons.fgc.modified_harmonic_constraints import activated_compatible_state
from recursive_horizons.fgc.modified_harmonic import modified_harmonic_symbol, standard_first_order_principal_system
from recursive_horizons.fgc.modified_harmonic_modes import (
    auxiliary_polynomial_mode_bases,
    exact_comp1_acceleration_root,
    mode1_point_certificate,
    quadratic_quotient_parameters,
    quotient_inverse,
    quotient_solve_3x3,
    ref1_branch_principal_matrix,
)
from recursive_horizons.fgc.reference_connection import flat_spherical_annulus_reference
from recursive_horizons.fgc.spherical_reduction import state_from_generalized_adm_pg_fixture


class MODE1Tests(unittest.TestCase):
    def _assert_branch_equals_mhg(self, state, radius: Q) -> None:
        reference = flat_spherical_annulus_reference(radial_domain_minimum=Q(1, 2))
        branch = ref1_branch_principal_matrix(state, reference=reference, coordinate_radius=radius, tilde_normal_factor=Q(4), hat_normal_factor=Q(9))
        mhg = standard_first_order_principal_system(modified_harmonic_symbol(state, tilde_normal_factor=Q(4), hat_normal_factor=Q(9)))
        self.assertEqual(branch["normal_principal_matrix"], mhg["normal_principal_matrix"])

    def test_flat_and_activated_fixture_branch_regressions(self) -> None:
        cfg = load_config()
        self._assert_branch_equals_mhg(state_from_generalized_adm_pg_fixture(cfg.flat_fixture), Q(4))
        self._assert_branch_equals_mhg(
            state_from_generalized_adm_pg_fixture(cfg.activated_fixture),
            cfg.coordinate_radii["FGCQR_activated_generic"],
        )

    def test_comp1_point_has_exact_auxiliary_bases_and_mhg_agreement(self) -> None:
        cfg = load_config()
        datum = activated_compatible_state(cfg.flat_fixture)
        out = mode1_point_certificate(datum["state"], reference=datum["reference"], coordinate_radius=Q(4), tilde_normal_factor=Q(4), hat_normal_factor=Q(9))
        self.assertTrue(out["ref1_branch_equals_mhg1"])
        self.assertTrue(out["all_auxiliary_modes_pointwise_semisimple"])

    def test_comp1_has_exact_polynomial_auxiliary_atlas(self) -> None:
        cfg = load_config()
        datum = activated_compatible_state(cfg.flat_fixture)
        root = exact_comp1_acceleration_root(
            datum["state"], reference=datum["reference"], coordinate_radius=Q(4),
            tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
        )
        self.assertTrue(root["solved_full_residual_zero"])
        self.assertTrue(root["root_strictly_inside_qift_acceleration_box"])
        atlas = auxiliary_polynomial_mode_bases(root["solved_state"], tilde_normal_factor=Q(4), hat_normal_factor=Q(9))
        self.assertTrue(atlas["tilde"]["all_residues_zero"])
        self.assertTrue(atlas["hat"]["all_residues_zero"])
        self.assertNotEqual(atlas["hat"]["pivot_unit_norm"], 0)
        self.assertTrue(atlas["hat"]["free_coordinate_identity_exact"])
        self.assertTrue(atlas["hat"]["rank_two_at_each_hat_root_by_free_coordinate_identity"])

    def test_generic_activated_fixture_has_polynomial_auxiliary_atlas(self) -> None:
        cfg = load_config()
        state = state_from_generalized_adm_pg_fixture(cfg.activated_fixture)
        atlas = auxiliary_polynomial_mode_bases(state, tilde_normal_factor=Q(4), hat_normal_factor=Q(9))
        self.assertTrue(atlas["tilde"]["all_residues_zero"])
        self.assertTrue(atlas["hat"]["all_residues_zero"])

    def test_quotient_and_atlas_fail_closed_on_corruption(self) -> None:
        with self.assertRaises(ValueError):
            quadratic_quotient_parameters((Q(1), Q(2), Q(3)))
        with self.assertRaises(TypeError):
            quadratic_quotient_parameters((Q(1), Q(2), 1))
        with self.assertRaises(ValueError):
            quotient_inverse((Q(0), Q(0)), b=Q(1), d=Q(1))
        zero = ((Q(0), Q(0)),) * 3
        with self.assertRaises(ValueError):
            quotient_solve_3x3((zero, zero, zero), zero, b=Q(1), d=Q(1))

        cfg = load_config()
        datum = activated_compatible_state(cfg.flat_fixture)
        bad_generator = tuple((Q(1),) if row == 0 else (Q(0),) for row in range(6))
        with patch("recursive_horizons.fgc.modified_harmonic_modes.right_gauge_generators", return_value=tuple((bad_generator[row], bad_generator[row]) for row in range(6))):
            with self.assertRaises(ValueError):
                auxiliary_polynomial_mode_bases(datum["state"], tilde_normal_factor=Q(4), hat_normal_factor=Q(9))

    def test_generic_rational_nonflat_nonzero_acceleration_regression(self) -> None:
        cfg = load_config()
        state = state_from_generalized_adm_pg_fixture(cfg.activated_fixture)
        self.assertNotEqual(state.phi.dtt, 0)
        self._assert_branch_equals_mhg(state, cfg.coordinate_radii["FGCQR_activated_generic"])


if __name__ == "__main__":
    unittest.main()
