from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_ecd_kinetic import (  # noqa: E402
    anticommutator_identity,
    contact_stress_null_energy_gate,
    effective_action_bookkeeping,
    free_dirac_positive_frequency_plane_wave_control,
    required_nonclaims,
    tetrad_current_invariant_cancellation,
)


class ECDKineticTests(unittest.TestCase):
    def setUp(self) -> None:
        self.radius = 1.0 / (2.0 * math.sqrt(math.pi))
        self.F, self.G = 0.75 + 0j, 0.5 + 0j
        self.metric = (
            (-1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0, 1.0),
        )

    def test_anticommutator_identity_is_exact(self) -> None:
        identity = anticommutator_identity()
        self.assertTrue(identity["verified_over_all_64_index_triples"])
        self.assertLess(identity["maximum_residual"], 1.0e-12)
        self.assertTrue(identity["representation_independent"])

    def test_effective_action_uses_int1_sign_without_double_counting(self) -> None:
        action = effective_action_bookkeeping(self.F, self.G, self.radius, 1.0, 1.0)
        self.assertAlmostEqual(action["canonical_J_squared"], 9.0 / 4.0)
        self.assertAlmostEqual(action["vc_A_squared"], -9.0 / 4.0)
        self.assertAlmostEqual(
            action["canonical_full_contact_coefficient_over_kappa_J2"], -3.0 / 16.0
        )
        self.assertAlmostEqual(
            action["vc_full_contact_coefficient_over_kappa_A2"], 3.0 / 16.0
        )
        self.assertAlmostEqual(action["canonical_full_contact_lagrangian"], -27.0 / 64.0)
        self.assertAlmostEqual(action["vc_full_contact_lagrangian"], -27.0 / 64.0)
        self.assertTrue(action["signature_translation_verified"])
        self.assertAlmostEqual(
            action["hehl_datta_equation_coefficient_over_kappa"], 3.0 / 8.0
        )
        self.assertTrue(action["field_variation_factor_two_verified"])
        self.assertFalse(action["separate_connection_pieces_added_to_reduced_action"])

    def test_tetrad_current_response_cancels_fixed_current_partial(self) -> None:
        cancellation = tetrad_current_invariant_cancellation()
        self.assertEqual(cancellation["coframe_directions_checked"], 16)
        self.assertTrue(cancellation["fixed_coordinate_current_partial_is_nonzero"])
        self.assertGreater(
            cancellation["maximum_fixed_coordinate_current_metric_partial"], 1.0e-3
        )
        self.assertLess(cancellation["baseline_invariant_residual"], 1.0e-12)
        self.assertLess(cancellation["maximum_total_variation_residual"], 1.0e-12)
        self.assertTrue(cancellation["current_response_cancels_metric_partial"])
        self.assertTrue(cancellation["internal_current_invariant_is_tetrad_independent"])
        self.assertFalse(cancellation["fixed_coordinate_current_anisotropic_stress_is_physical"])

    def test_contact_stress_has_zero_null_contraction(self) -> None:
        contact = contact_stress_null_energy_gate(
            1.0,
            -9.0 / 4.0,
            self.metric,
            ((1.0, 1.0, 0.0, 0.0), (1.0, 0.0, 0.0, 1.0)),
        )
        self.assertAlmostEqual(contact["rho"], 27.0 / 64.0)
        self.assertAlmostEqual(contact["p"], -27.0 / 64.0)
        self.assertAlmostEqual(contact["w"], -1.0)
        self.assertAlmostEqual(contact["rho_plus_p"], 0.0)
        self.assertEqual(contact["null_vectors_checked"], 2)
        self.assertTrue(contact["all_contact_null_contractions_zero"])
        self.assertFalse(contact["contact_null_energy_condition_violated"])
        self.assertTrue(contact["contact_sector_tetrad_variation_complete"])

    def test_contact_gate_rejects_non_null_vector(self) -> None:
        with self.assertRaises(ValueError):
            contact_stress_null_energy_gate(
                1.0, -9.0 / 4.0, self.metric, ((1.0, 0.0, 0.0, 0.0),)
            )

    def test_free_dirac_plane_wave_control_is_bounded(self) -> None:
        control = free_dirac_positive_frequency_plane_wave_control()
        self.assertAlmostEqual(control["momentum_norm"], 0.0)
        self.assertLess(control["massless_dirac_residual"], 1.0e-12)
        self.assertLess(control["current_equals_two_momentum_residual"], 1.0e-12)
        self.assertAlmostEqual(control["null_contractions"]["parallel"], 0.0)
        self.assertAlmostEqual(control["null_contractions"]["opposite"], 8.0)
        self.assertTrue(control["positive_frequency_plane_wave_nec_nonnegative"])
        self.assertTrue(control["positive_frequency_plane_wave_control_only"])
        self.assertFalse(control["general_classical_dirac_nec_proven"])

    def test_validation_and_nonclaims(self) -> None:
        with self.assertRaises(ValueError):
            anticommutator_identity(0.0)
        with self.assertRaises(ValueError):
            effective_action_bookkeeping(True, 0.5, self.radius, 1.0, 1.0)
        with self.assertRaises(ValueError):
            effective_action_bookkeeping(0.75, 0.5, 0.0, 1.0, 1.0)
        nonclaims = required_nonclaims()
        self.assertFalse(nonclaims["full_ecd_tetrad_dependence_varied"])
        self.assertFalse(nonclaims["torsion_free_dirac_stress_fully_derived"])
        self.assertFalse(nonclaims["general_classical_dirac_nec_proven"])
        self.assertFalse(nonclaims["rho_plus_p_negative_derived"])
        self.assertFalse(nonclaims["metric_null_defocusing_derived"])
        self.assertFalse(nonclaims["bounce_constructed"])
        self.assertFalse(nonclaims["global_solution_constructed"])


if __name__ == "__main__":
    unittest.main()
