from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1]), str(Path(__file__).resolve().parents[1] / "src")]

from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config
from recursive_horizons.fgc.exact_interval import interval
from recursive_horizons.fgc.modified_harmonic_constraints import activated_compatible_state
from recursive_horizons.fgc.modified_harmonic_first_order import exact_full_residual_first_order_argument_jacobian
from recursive_horizons.fgc.modified_harmonic_interval_principal import (
    compact_ref1_interval_principal_certificate,
    interval_ref1_first_order_blocks,
    kinetic_neumann_inverse_enclosure,
)
from recursive_horizons.fgc.modified_harmonic_modes import (
    exact_comp1_acceleration_root,
    ref1_branch_principal_matrix,
)
from recursive_horizons.fgc.spherical_reduction import BASE_FIELD_ORDER


class IntervalPrincipalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cfg = load_config()
        datum = activated_compatible_state(cfg.flat_fixture)
        cls.reference = datum["reference"]
        cls.center = exact_comp1_acceleration_root(
            datum["state"], reference=cls.reference, coordinate_radius=Q(4),
            tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
        )["solved_state"]

    def test_point_collapse_recovers_complete_ref1_and_mode1_branch(self) -> None:
        boxed = interval_ref1_first_order_blocks(
            self.center, parameter_half_width=Q(0), acceleration_half_width=Q(0),
            reference=self.reference, coordinate_radius=Q(4),
            tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
        )
        exact = exact_full_residual_first_order_argument_jacobian(
            self.center, reference=self.reference, coordinate_radius=Q(4),
            tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
        )["blocks"]
        for name in ("p_t", "p_r", "q_r"):
            self.assertEqual(
                boxed["blocks"][name],
                tuple(tuple(interval(value) for value in row) for row in exact[name]),
            )
        certificate = compact_ref1_interval_principal_certificate(
            self.center, reference=self.reference, coordinate_radius=Q(4),
            tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
            parameter_half_width=Q(1, 2**60), acceleration_half_width=Q(1, 2**50),
        )
        mode1 = ref1_branch_principal_matrix(
            self.center, reference=self.reference, coordinate_radius=Q(4),
            tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
        )["normal_principal_matrix"]
        for row, exact_row in enumerate(mode1):
            for column, value in enumerate(exact_row):
                self.assertTrue(certificate["solved_branch_principal_box"][row][column].lower <= value <= certificate["solved_branch_principal_box"][row][column].upper)

    def test_compact_box_has_all_axes_and_a_contracting_kinetic_inverse(self) -> None:
        certificate = compact_ref1_interval_principal_certificate(
            self.center, reference=self.reference, coordinate_radius=Q(4),
            tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
        )
        self.assertTrue(certificate["box"]["all_30_parameter_axes_have_nonzero_width"])
        self.assertTrue(certificate["box"]["all_6_acceleration_axes_have_nonzero_width"])
        self.assertTrue(certificate["box"]["seeded_primals_reproduce_same_complete_residual_box"])
        self.assertLess(certificate["kinetic_neumann_inverse"]["neumann_rho_infinity_upper_bound"], 1)
        self.assertTrue(certificate["kinetic_neumann_inverse"]["kinetic_inverse_regular_over_entire_box"])
        self.assertTrue(certificate["complete_REF1_acceleration_Krawczyk"]["krawczyk_image_strictly_inside_acceleration_box"])
        self.assertTrue(certificate["complete_REF1_acceleration_Krawczyk"]["unique_acceleration_root_for_every_declared_parameter_point"])
        self.assertTrue(certificate["regularity"]["all_declared_denominators_and_sign_branches_regular"])
        self.assertFalse(certificate["nonclaims"]["uniform_strong_hyperbolicity_proven"])

    def test_kinetic_and_input_fail_closed(self) -> None:
        zero = tuple(tuple(Q(0) for _ in BASE_FIELD_ORDER) for _ in BASE_FIELD_ORDER)
        with self.assertRaises(ValueError):
            kinetic_neumann_inverse_enclosure(zero, tuple(tuple(interval(0) for _ in BASE_FIELD_ORDER) for _ in BASE_FIELD_ORDER))
        identity = tuple(tuple(Q(int(row == column)) for column in range(len(BASE_FIELD_ORDER))) for row in range(len(BASE_FIELD_ORDER)))
        expansive = tuple(tuple(interval(3 if row == column else 0) for column in range(len(BASE_FIELD_ORDER))) for row in range(len(BASE_FIELD_ORDER)))
        with self.assertRaises(ValueError):
            kinetic_neumann_inverse_enclosure(identity, expansive)
        with self.assertRaises(ValueError):
            interval_ref1_first_order_blocks(
                self.center, parameter_half_width=Q(1, 2**60), acceleration_half_width=Q(1, 2**50),
                reference=self.reference, coordinate_radius=Q(1, 2),
                tilde_normal_factor=Q(4), hat_normal_factor=Q(9),
            )


if __name__ == "__main__":
    unittest.main()
