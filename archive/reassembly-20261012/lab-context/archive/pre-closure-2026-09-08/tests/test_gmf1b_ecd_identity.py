from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_ecd_identity import (  # noqa: E402
    contact_tetrad_stress,
    metric_null_contact_contraction,
    required_nonclaims,
    vc_action_projection_identity,
    vc_contact_action_sources,
    vc_hehl_datta_projection,
    vc_reduced_contact_density,
    vc_signature_bridge,
)
from recursive_horizons.gmf1b_ecd_symmetry import ventrella_chiral_axial_bilinears  # noqa: E402


class ECDIdentityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.angles = ((0.7, 0.2), (1.1, -0.9), (2.0, 1.4))

    def test_exact_contact_fixture(self) -> None:
        # C=1 at r=1/(2 sqrt(pi)), a=1.
        radius = 1.0 / (2.0 * math.sqrt(math.pi))
        self.assertAlmostEqual(vc_reduced_contact_density(0.75, 0.5, 1.0, 1.0, radius, 1.0), -27.0 / 64.0)
        sources = vc_contact_action_sources(0.75, 0.5, 1.0, 1.0, radius, 1.0)
        self.assertAlmostEqual(sources["source_F"].real, -9.0 / 16.0)
        self.assertAlmostEqual(sources["source_G"].real, -27.0 / 32.0)

    def test_action_projection_identity_and_one_amplitude_semantics(self) -> None:
        identity = vc_action_projection_identity(0.75 + 0.2j, 0.5 - 0.3j, 1.4, 1.1, 1.3, 0.9, self.angles)
        self.assertTrue(identity["weak_form_identity_closed"])
        self.assertTrue(identity["all_tested_angle_reconstructions_closed"])
        zero = vc_action_projection_identity(1.0, 0.0, 1.0, 1.0, 1.0, 1.0, self.angles)
        self.assertTrue(zero["weak_form_identity_closed"])
        axial = ventrella_chiral_axial_bilinears(1.0, 0.0, 1.0, 1.0)
        self.assertTrue(axial["axial_current_nonzero"])
        self.assertFalse(axial["axial_contact_invariant_nonzero"])

    def test_direct_physical_hehl_datta_phase_and_local_basis_reconstruction(self) -> None:
        radius = 1.0 / (2.0 * math.sqrt(math.pi))
        projection = vc_hehl_datta_projection(0.75, 0.5, 1.0, radius, 1.0, 0.7, 0.2)
        self.assertAlmostEqual(projection["expected_retained_F"].real, -9.0 / 32.0)
        self.assertAlmostEqual(projection["expected_retained_F"].imag, 0.0)
        self.assertAlmostEqual(projection["expected_retained_G"].real, -27.0 / 64.0)
        self.assertTrue(projection["local_retained_basis_reconstructed"])
        self.assertTrue(projection["coefficient_sign_phase_closed"])

    def test_full_signature_contact_and_equation_bridge(self) -> None:
        radius = 1.0 / (2.0 * math.sqrt(math.pi))
        bridge = vc_signature_bridge(0.75, 0.5, 1.0, radius, 1.0, 0.7, 0.2)
        self.assertTrue(bridge["signature_bridge_closed"])
        self.assertTrue(bridge["canonical_action_defined_before_translation"])
        self.assertEqual(bridge["vc_clifford_residual"], 0.0)
        self.assertEqual(bridge["canonical_clifford_residual"], 0.0)
        self.assertEqual(bridge["gamma5_translation_residual"], 0.0)
        self.assertEqual(bridge["adjoint_translation_residual"], 0.0)
        self.assertAlmostEqual(bridge["raw_vc_axial_bilinear"][0].real, 0.0)
        self.assertAlmostEqual(bridge["raw_vc_axial_bilinear"][0].imag, 13.0 / 8.0)
        self.assertAlmostEqual(bridge["J_squared_plus"], 9.0 / 4.0)
        self.assertAlmostEqual(bridge["A_squared_minus"], -9.0 / 4.0)
        self.assertAlmostEqual(bridge["contact_lagrangian_plus"], -27.0 / 64.0)
        self.assertAlmostEqual(bridge["contact_lagrangian_minus"], -27.0 / 64.0)
        self.assertEqual(bridge["canonical_lhs_to_vc_equation_residual"], 0.0)
        self.assertLess(bridge["maximum_signature_bridge_residual"], 1.0e-12)
        self.assertFalse(bridge["line_element_physically_changed"])

    def test_contact_stress_and_null_contraction(self) -> None:
        metric = ((-1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0), (0.0, 0.0, 0.0, 1.0))
        stress = contact_tetrad_stress(2.0, -3.0, metric)
        self.assertEqual(stress[0][0], 9.0 / 8.0)
        null = metric_null_contact_contraction(
            2.0, -3.0, metric, (1.0, 1.0, 0.0, 0.0)
        )
        self.assertEqual(null["contact_null_contraction"], 0.0)
        self.assertTrue(null["contact_sector_tetrad_variation_complete"])
        self.assertFalse(null["full_ecd_stress_derived"])
        with self.assertRaises(ValueError):
            metric_null_contact_contraction(
                1.0, 1.0, metric, (1.0, 0.0, 0.0, 0.0)
            )

    def test_validation_json_and_nonclaims(self) -> None:
        with self.assertRaises(ValueError):
            vc_reduced_contact_density(True, 1.0, 1.0, 1.0, 1.0, 1.0)
        with self.assertRaises(ValueError):
            vc_action_projection_identity(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, ((0.0, 0.0),))
        payload = {"identity": vc_action_projection_identity(1.0, 2.0, 1.0, 1.0, 1.0, 1.0, self.angles), "nonclaims": required_nonclaims()}
        json.dumps(payload, default=lambda value: {"real": value.real, "imag": value.imag}, allow_nan=False)
        self.assertFalse(required_nonclaims()["global_solution_constructed"])


if __name__ == "__main__":
    unittest.main()
