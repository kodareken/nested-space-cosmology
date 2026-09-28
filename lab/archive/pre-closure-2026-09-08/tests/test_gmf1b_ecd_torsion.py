from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_ecd_torsion import (  # noqa: E402
    contact_tetrad_stress_interpretation,
    effective_contact_lagrangian_crosscheck,
    required_nonclaims,
    spin_tensor_components,
    torsion_quadratic_invariant,
)


class ECDTorsionTests(unittest.TestCase):
    def setUp(self) -> None:
        # C = 1/(2 sqrt(pi) r sqrt(a)) = 1 at r = 1/(2 sqrt(pi)), a = 1.
        self.radius = 1.0 / (2.0 * math.sqrt(math.pi))
        self.F, self.G = 0.75 + 0j, 0.5 + 0j
        self.metric = (
            (-1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 1.0),
        )

    def test_spin_tensor_is_hodge_dual_of_lowered_axial_current(self) -> None:
        spin = spin_tensor_components(self.F, self.G, self.radius, 1.0, 0.7, 0.2)
        self.assertTrue(spin["hodge_dual_verified"])
        self.assertLess(spin["hodge_dual_residual"], 1.0e-12)
        self.assertAlmostEqual(spin["spin_tensor_components"]["S_012"], 5.0 / 16.0)
        self.assertAlmostEqual(spin["spin_tensor_components"]["S_123"], 13.0 / 16.0)
        self.assertTrue(spin["totally_antisymmetric"])
        # The spin tensor is angle-independent.
        other = spin_tensor_components(self.F, self.G, self.radius, 1.0, 1.1, -0.9)
        self.assertAlmostEqual(other["spin_tensor_components"]["S_012"], 5.0 / 16.0)
        self.assertTrue(other["hodge_dual_verified"])

    def test_torsion_quadratic_relation(self) -> None:
        torsion = torsion_quadratic_invariant(self.F, self.G, self.radius, 1.0, 1.0, 0.7, 0.2)
        self.assertTrue(torsion["relation_verified"])
        self.assertAlmostEqual(torsion["torsion_squared"], 27.0 / 8.0)
        self.assertAlmostEqual(torsion["expected_minus_3_2_kappa2_A2"], 27.0 / 8.0)

    def test_full_effective_contact_lagrangian_crosscheck(self) -> None:
        contact = effective_contact_lagrangian_crosscheck(
            self.F, self.G, self.radius, 1.0, 1.0, 0.7, 0.2
        )
        self.assertTrue(contact["coefficient_verified"])
        self.assertAlmostEqual(contact["torsion_squared_consistency_input"], 27.0 / 8.0)
        self.assertAlmostEqual(contact["full_effective_contact_lagrangian"], -27.0 / 64.0)
        self.assertAlmostEqual(contact["expected_3_16_kappa_A2"], -27.0 / 64.0)
        self.assertTrue(contact["derived_after_eliminating_connection_everywhere"])
        self.assertFalse(contact["einstein_hilbert_piece_alone_claimed_to_equal_full_contact"])
        self.assertFalse(contact["separate_connection_pieces_added_to_reduced_action"])

    def test_contact_tetrad_stress_is_positive_cosmological_constant(self) -> None:
        spin = spin_tensor_components(self.F, self.G, self.radius, 1.0, 0.7, 0.2)
        stress = contact_tetrad_stress_interpretation(
            1.0, spin["axial_squared"], self.metric
        )
        self.assertAlmostEqual(stress["rho"], 27.0 / 64.0)
        self.assertAlmostEqual(stress["p"], -27.0 / 64.0)
        self.assertAlmostEqual(stress["w"], -1.0)
        self.assertAlmostEqual(stress["rho_plus_p"], 0.0)
        self.assertAlmostEqual(stress["rho_plus_3p"], -27.0 / 32.0)
        self.assertEqual(stress["interpretation"], "positive_cosmological_constant_w_minus_1")
        self.assertTrue(stress["contact_sector_tetrad_variation_complete"])
        self.assertFalse(stress["rho_plus_p_negative"])
        self.assertTrue(stress["accelerates"])

    def test_validation_and_nonclaims(self) -> None:
        with self.assertRaises(ValueError):
            spin_tensor_components(0.75, 0.5, 0.0, 1.0, 0.7, 0.2)
        with self.assertRaises(ValueError):
            spin_tensor_components(True, 0.5, 1.0, 1.0, 0.7, 0.2)
        with self.assertRaises(ValueError):
            torsion_quadratic_invariant(0.75, 0.5, 1.0, 1.0, 0.0, 0.7, 0.2)
        nonclaims = required_nonclaims()
        self.assertFalse(nonclaims["full_ecd_tetrad_dependence_varied"])
        self.assertFalse(nonclaims["rho_plus_p_negative_derived"])
        self.assertFalse(nonclaims["dark_energy_value_fixed"])
        self.assertFalse(nonclaims["bounce_constructed"])
        self.assertFalse(nonclaims["global_solution_constructed"])


if __name__ == "__main__":
    unittest.main()
