from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_ecd_symmetry import (  # noqa: E402
    audit_chiral_annular_packet,
    double_null_metric_null_cone,
    finster_equal_profile_axial_control,
    geometric_normalization,
    minimal_ecd_contact_scalar,
    normalized_annular_bump,
    required_nonclaims,
    ventrella_chiral_axial_bilinears,
    ventrella_gamma_matrix_certificate,
)


class GMF1BECDSymmetryTests(unittest.TestCase):
    fixture_radius = 1.0 / (2.0 * math.sqrt(math.pi))

    def test_exact_ecd_fixture_coefficients(self) -> None:
        self.assertAlmostEqual(geometric_normalization(self.fixture_radius, 1.0), 1.0, places=14)
        axial = ventrella_chiral_axial_bilinears(0.75, 0.5, self.fixture_radius, 1.0)
        self.assertAlmostEqual(axial["A_hat_0"], 13.0 / 8.0, places=14)
        self.assertAlmostEqual(axial["A_hat_r"], 5.0 / 8.0, places=14)
        self.assertAlmostEqual(axial["axial_squared"], -9.0 / 4.0, places=14)
        self.assertAlmostEqual(minimal_ecd_contact_scalar(1.0, axial["axial_squared"]), -27.0 / 64.0, places=14)

    def test_finster_control_cancels_and_ventrella_pair_survives(self) -> None:
        control = finster_equal_profile_axial_control(
            0.75 + 0.2j,
            0.5 - 0.3j,
            normalization=1.25,
            theta=0.7,
            phi=-0.4,
        )
        self.assertTrue(control["profiles_populated"])
        self.assertEqual(control["gamma5_definition"], "i*gamma0*gamma1*gamma2*gamma3")
        self.assertFalse(control["axial_current_nonzero"])
        self.assertFalse(control["algebraic_torsion_source_nonzero"])
        self.assertFalse(control["axial_contact_invariant_nonzero"])
        self.assertTrue(
            all(abs(float(control[key])) < 1.0e-12 for key in ("A_hat_0", "A_hat_x", "A_hat_y", "A_hat_z", "axial_squared"))
        )
        self.assertTrue(control["single_spinor_1_axial_current_nonzero"])
        self.assertGreater(control["single_spinor_1_max_abs_axial_component"], 1.0e-6)
        self.assertLess(control["single_spinor_1_max_imaginary_residual"], 1.0e-12)
        self.assertTrue(control["pauli_doublet_cancellation_certified_to_binary64_roundoff"])
        self.assertLess(control["maximum_direct_axial_residual"], 1.0e-12)
        axial = ventrella_chiral_axial_bilinears(0.75, 0.5, self.fixture_radius, 1.0)
        self.assertTrue(axial["axial_current_nonzero"])
        self.assertTrue(axial["algebraic_torsion_source_nonzero"])
        self.assertTrue(axial["axial_contact_invariant_nonzero"])
        self.assertTrue(axial["contact_action_density_nonzero"])
        self.assertTrue(axial["so3_axial_current_compatible"])
        self.assertFalse(axial["o3_parity_compatibility_proven"])
        self.assertFalse(axial["full_ecd_contact_stress_verified"])

    def test_zero_amplitude_controls_and_gamma_certificate(self) -> None:
        for F, G in ((0.0, 0.5), (0.75, 0.0)):
            state = ventrella_chiral_axial_bilinears(F, G, self.fixture_radius, 1.0)
            self.assertEqual(state["axial_squared"], 0.0)
            self.assertTrue(state["axial_current_nonzero"])
            self.assertTrue(state["algebraic_torsion_source_nonzero"])
            self.assertFalse(state["axial_contact_invariant_nonzero"])
            self.assertFalse(state["contact_action_density_nonzero"])
        zero = ventrella_chiral_axial_bilinears(0.0, 0.0, self.fixture_radius, 1.0)
        self.assertFalse(zero["axial_current_nonzero"])
        self.assertFalse(zero["algebraic_torsion_source_nonzero"])
        self.assertFalse(zero["axial_contact_invariant_nonzero"])
        self.assertFalse(zero["contact_action_density_nonzero"])
        cert = ventrella_gamma_matrix_certificate(0.75, 0.5, self.fixture_radius, 1.0, 0.7, -0.4)
        self.assertTrue(cert["formula_certified_to_binary64_roundoff"])
        self.assertLess(cert["maximum_formula_residual"], 1.0e-12)

    def test_metric_null_cone(self) -> None:
        for u, v in ((1.0, 0.0), (0.0, 1.0)):
            state = double_null_metric_null_cone(1.0, u, v)
            self.assertTrue(state["metric_null_characteristic"])
            self.assertEqual(state["dirac_principal_determinant_factor"], 0.0)
        state = double_null_metric_null_cone(1.0, 1.0, 1.0)
        self.assertFalse(state["metric_null_characteristic"])
        self.assertFalse(state["algebraic_torsion_changes_metric_cone"])

    def test_compact_annular_support_and_257_finite_samples(self) -> None:
        for radius in (0.0, 0.5, 1.0, 2.0, 2.5):
            self.assertEqual(normalized_annular_bump(radius, 1.0, 2.0), 0.0)
        self.assertEqual(normalized_annular_bump(1.5, 1.0, 2.0), 1.0)
        finite = []
        nonzero = []
        for index in range(257):
            audit = audit_chiral_annular_packet(0.5 + 2.0 * index / 256.0, 1.0, 2.0, 0.75, 0.5, 1.0, 1.0)
            axial = audit["axial"]
            finite.extend(value for value in axial.values() if isinstance(value, float))
            nonzero.append(abs(float(axial["axial_squared"])) > 0.0)
        self.assertTrue(all(math.isfinite(value) for value in finite))
        self.assertTrue(any(nonzero))
        for radius in (1.0, 2.0):
            audit = audit_chiral_annular_packet(radius, 1.0, 2.0, 0.75, 0.5, 1.0, 1.0)
            self.assertTrue(audit["axial_current_vanishes_at_this_outside_support_point"])

    def test_strict_invalid_input_and_complete_nonclaims(self) -> None:
        for bad in (True, "0.5", complex(math.inf, 0.0), math.inf, math.nan):
            with self.assertRaises(ValueError):
                ventrella_chiral_axial_bilinears(bad, 0.5, 1.0, 1.0)  # type: ignore[arg-type]
        for args in ((1.5, 0.0, 2.0, 0.75, 0.5, 1.0, 1.0), (1.5, 2.0, 1.0, 0.75, 0.5, 1.0, 1.0)):
            with self.assertRaises(ValueError):
                audit_chiral_annular_packet(*args)
        with self.assertRaises(ValueError):
            minimal_ecd_contact_scalar(0.0, -1.0)
        nonclaims = required_nonclaims()
        self.assertTrue(all(value is False for value in nonclaims.values()))
        json.dumps(nonclaims, allow_nan=False)

    def test_reproducer_record_is_json_safe_and_bounded(self) -> None:
        from scripts.reproduce_gmf1b_ecd_symmetry import record

        payload = record()
        self.assertEqual(payload["schema_version"], 1)
        self.assertEqual(payload["project_version"], "0.11.0")
        self.assertEqual(payload["model_id"], "GMF-1B-ECD-SYM1")
        self.assertFalse(payload["global_solution_constructed"])
        finster = payload["finster_full_dirac_singlet_control"]
        self.assertTrue(finster["pauli_doublet_cancellation_certified_to_binary64_roundoff"])
        self.assertTrue(finster["single_spinor_1_axial_current_nonzero"])
        self.assertLess(finster["maximum_direct_axial_residual"], 1.0e-12)
        controls = payload["ventrella_left_chiral_pair"]["one_amplitude_zero_controls"]
        for control in controls.values():
            self.assertTrue(control["algebraic_torsion_source_nonzero"])
            self.assertFalse(control["axial_contact_invariant_nonzero"])
        annular = payload["annular_packet_gate"]
        self.assertTrue(annular["centre_spinor_zero_buffer"])
        self.assertTrue(annular["outer_spinor_zero_buffer"])
        self.assertTrue(annular["sampled_axial_bilinears_and_contact_scalar_finite"])
        self.assertFalse(payload["constraint_compatible_vacuum_region_verified"])
        json.dumps(payload, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
