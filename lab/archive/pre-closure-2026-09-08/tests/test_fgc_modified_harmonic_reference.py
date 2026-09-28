from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.reproduce_fgc_hyp1_reduction import load_config  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_reference import (  # noqa: E402
    modified_harmonic_gauge_constraint,
    modified_harmonic_reference_certificate,
    modified_harmonic_reference_fixture_certificate,
    physical_connection_data,
    required_ref1_nonclaims,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_connection,
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    state_from_generalized_adm_pg_fixture,
)


class FGCModifiedHarmonicReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fixtures = load_config().configuration["fixtures"]
        cls.reference = flat_spherical_annulus_reference(
            radial_domain_minimum=Q(1, 2)
        )
        cls.radii = {
            "GR0_regular": Q(4),
            "SGBL_constant_f": Q(5),
            "FGCQR_schwarzschild_exterior": Q(8),
            "FGCQR_activated_generic": Q(9, 2),
        }
        cls.certificate = modified_harmonic_reference_certificate(
            {
                "fixtures": cls.fixtures,
                "coordinate_radii": cls.radii,
                "reference": cls.reference,
                "tilde_normal_factor": Q(4),
                "hat_normal_factor": Q(9),
                "gauge_surface_fixture_ids": (
                    "GR0_regular",
                    "SGBL_constant_f",
                ),
                "off_gauge_fixture_ids": (
                    "FGCQR_schwarzschild_exterior",
                    "FGCQR_activated_generic",
                ),
            }
        )

    def test_flat_reference_matches_independent_physical_connection(self) -> None:
        fixture = self.fixtures[0]
        state = state_from_generalized_adm_pg_fixture(fixture)
        physical = physical_connection_data(state)
        reference = flat_spherical_annulus_connection(
            self.reference, coordinate_radius=Q(4)
        )
        self.assertEqual(physical.christoffel, reference.christoffel)
        self.assertEqual(
            physical.christoffel_derivative,
            reference.christoffel_derivative,
        )

    def test_gauge_surface_and_off_gauge_controls(self) -> None:
        records = {
            record["fixture_id"]: record
            for record in self.certificate["fixture_records"]
        }
        for fixture_id in ("GR0_regular", "SGBL_constant_f"):
            record = records[fixture_id]
            self.assertTrue(record["gauge_surface_at_evaluation_jet"])
            self.assertFalse(record["extension_nonzero"])
            self.assertEqual(
                record["full_residual_vector"],
                record["unredefined_residual_vector"],
            )

        schwarzschild = records["FGCQR_schwarzschild_exterior"]
        activated = records["FGCQR_activated_generic"]
        self.assertEqual(schwarzschild["constraint_up"][:2], (Q(3, 32), Q(1, 64)))
        self.assertTrue(all(value == 0 for value in schwarzschild["unredefined_residual_vector"]))
        for record in (schwarzschild, activated):
            self.assertFalse(record["gauge_surface_at_evaluation_jet"])
            self.assertTrue(record["extension_nonzero"])
            self.assertEqual(record["constraint_up"][2:], (Q(0), Q(0)))

    def test_full_extension_regresses_to_every_mhg1_principal_gauge_entry(self) -> None:
        self.assertTrue(self.certificate["all_declared_exact_checks_pass"])
        self.assertTrue(all(self.certificate["verified_exact_checks"].values()))
        for record in self.certificate["fixture_records"]:
            with self.subTest(fixture=record["fixture_id"]):
                self.assertTrue(
                    record["reference_extension_principal_regression_exact"]
                )
                self.assertEqual(
                    record["reference_extension_principal_symbol"],
                    record["mhg1_principal_gauge_symbol"],
                )
                self.assertTrue(record["extension_tensor_symmetric"])
                self.assertTrue(record["equatorial_angular_isotropy"])
                self.assertTrue(record["scalar_equations_unmodified"])

    def test_nonzero_shift_keeps_contravariant_constraint_order_visible(self) -> None:
        fixture = self.fixtures[2]
        state = state_from_generalized_adm_pg_fixture(fixture)
        gauge = modified_harmonic_gauge_constraint(
            state,
            reference=self.reference,
            coordinate_radius=Q(8),
            tilde_normal_factor=Q(4),
        )
        reconstructed_down = tuple(
            sum(
                gauge.physical.metric[a][b] * gauge.constraint_up[b]
                for b in range(4)
            )
            for a in range(4)
        )
        self.assertEqual(gauge.constraint_down, reconstructed_down)
        self.assertNotEqual(gauge.constraint_up[:2], gauge.constraint_down[:2])

    def test_auxiliary_factor_mutation_changes_the_full_gauge_data(self) -> None:
        fixture = self.fixtures[3]
        state = state_from_generalized_adm_pg_fixture(fixture)
        factor_4 = modified_harmonic_gauge_constraint(
            state,
            reference=self.reference,
            coordinate_radius=Q(9, 2),
            tilde_normal_factor=Q(4),
        )
        factor_5 = modified_harmonic_gauge_constraint(
            state,
            reference=self.reference,
            coordinate_radius=Q(9, 2),
            tilde_normal_factor=Q(5),
        )
        self.assertNotEqual(factor_4.constraint_up, factor_5.constraint_up)
        self.assertNotEqual(
            factor_4.covariant_constraint_derivative,
            factor_5.covariant_constraint_derivative,
        )

    def test_fail_closed_annulus_configuration_and_nonclaims(self) -> None:
        with self.assertRaisesRegex(ValueError, "minimum must be positive"):
            flat_spherical_annulus_reference(radial_domain_minimum=Q(0))
        with self.assertRaisesRegex(ValueError, "strictly inside"):
            flat_spherical_annulus_connection(
                self.reference, coordinate_radius=Q(1, 2)
            )
        with self.assertRaisesRegex(ValueError, "unexpected keys"):
            modified_harmonic_reference_certificate(
                {
                    "fixtures": self.fixtures,
                    "coordinate_radii": self.radii,
                    "reference": self.reference,
                    "tilde_normal_factor": Q(4),
                    "hat_normal_factor": Q(9),
                    "gauge_surface_fixture_ids": (),
                    "off_gauge_fixture_ids": (),
                    "extra": True,
                }
            )
        self.assertTrue(all(value is False for value in required_ref1_nonclaims().values()))

    def test_single_fixture_certificate_is_exact(self) -> None:
        record = modified_harmonic_reference_fixture_certificate(
            self.fixtures[0],
            reference=self.reference,
            coordinate_radius=Q(4),
            tilde_normal_factor=Q(4),
            hat_normal_factor=Q(9),
        )
        self.assertTrue(record["gauge_surface_at_evaluation_jet"])
        self.assertTrue(record["reference_extension_principal_regression_exact"])


if __name__ == "__main__":
    unittest.main()
