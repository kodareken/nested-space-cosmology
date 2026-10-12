from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_reduction import load_config  # noqa: E402
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    state_from_generalized_adm_pg_fixture,
)
from recursive_horizons.fgc.exact_linear_algebra import (  # noqa: E402
    matrix_mul,
    poly_mul,
    poly_scale,
)
from recursive_horizons.fgc.modified_harmonic import (  # noqa: E402
    auxiliary_inverse_metric,
    gauge_constraint_propagation_principal_symbol,
    modified_harmonic_certificate,
    required_mhg1_nonclaims,
)


class FGCModifiedHarmonicTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fixtures = load_config().configuration["fixtures"]
        cls.certificate = modified_harmonic_certificate(
            {
                "fixtures": cls.fixtures,
                "tilde_normal_factor": Q(4),
                "hat_normal_factor": Q(9),
            }
        )

    def test_complete_factorization_cones_and_pointwise_eigenbases(self) -> None:
        self.assertTrue(self.certificate["all_declared_exact_checks_pass"])
        self.assertTrue(all(self.certificate["verified_exact_checks"].values()))
        for record in self.certificate["fixture_records"]:
            with self.subTest(fixture=record["fixture_id"]):
                factorization = record["determinant_factorization"]
                expected = poly_mul(
                    poly_mul(
                        poly_mul(
                            factorization["tilde_null_monic"],
                            factorization["tilde_null_monic"],
                        ),
                        poly_mul(
                            factorization["hat_null_monic"],
                            factorization["hat_null_monic"],
                        ),
                    ),
                    factorization["physical_scalar_determinant_monic"],
                )
                self.assertEqual(
                    record["complete_mhg_determinant"],
                    poly_scale(expected, factorization["nonzero_constant_quotient"]),
                )
                self.assertTrue(record["auxiliary_cones_separated_from_all_physical_modes"])
                self.assertTrue(record["pointwise_complete_real_first_order_eigenbasis"])
                self.assertTrue(
                    record["coordinate_time_kinetic_matches_characteristic_leading_coefficient"]
                )
                self.assertTrue(
                    record[
                        "first_order_characteristic_determinant_equals_second_order"
                    ]
                )
                self.assertEqual(
                    record["standard_first_order_principal_system"][
                        "characteristic_pencil_determinant"
                    ],
                    record["complete_mhg_determinant"],
                )

    def test_repeated_rational_modes_are_semisimple_and_activated_roots_are_exact(self) -> None:
        records = {
            record["fixture_id"]: record
            for record in self.certificate["fixture_records"]
        }
        for record in records.values():
            for mode in record["rational_characteristic_modes"].values():
                self.assertEqual(
                    mode["kernel_dimensions"],
                    (mode["expected_kernel_dimension"],) * 2,
                )
                self.assertEqual(
                    mode["second_order_kernel_dimensions"],
                    mode["first_order_kernel_dimensions"],
                )
                self.assertEqual(
                    mode["first_order_kernel_dimensions"],
                    mode["kernel_dimensions"],
                )
                self.assertTrue(mode["second_to_first_order_kernel_map_exact"])
                self.assertTrue(mode["semisimple_at_each_rational_root"])

        activated = records["FGCQR_activated_generic"]
        self.assertEqual(
            activated["rational_characteristic_modes"]["tilde_gauge"]["roots"],
            (Q(-17, 84), Q(73, 84)),
        )
        self.assertEqual(
            activated["rational_characteristic_modes"]["hat_constraint"]["roots"],
            (Q(-1, 42), Q(29, 42)),
        )
        self.assertEqual(
            activated["rational_characteristic_modes"]["physical_metric"]["roots"],
            (Q(-31, 42), Q(59, 42)),
        )
        self.assertNotEqual(
            activated["cone_resultants"]["physical_vs_regulator"], Q(0)
        )

    def test_standard_first_order_block_reconstructs_exactly(self) -> None:
        for record in self.certificate["fixture_records"]:
            first_order = record["standard_first_order_principal_system"]
            self.assertNotEqual(first_order["A_time_determinant"], Q(0))
            self.assertEqual(
                matrix_mul(
                    first_order["A_time"],
                    first_order["normal_principal_matrix"],
                ),
                first_order["A_radial"],
            )

    def test_hat_cone_propagation_x_zero_controls_and_nonclaims(self) -> None:
        records = self.certificate["fixture_records"]
        for record in records:
            with self.subTest(fixture=record["fixture_id"]):
                propagation = record[
                    "gauge_constraint_propagation_principal_symbol"
                ]
                self.assertTrue(record["gauge_constraint_propagation_hat_null_exact"])
                self.assertTrue(propagation["full_projector_contraction_exact"])
                self.assertEqual(
                    propagation["full_projector_contraction"],
                    propagation["expected_full_hat_wave_operator"],
                )
                self.assertEqual(
                    propagation["expected_hat_wave_factor_monic"],
                    record["determinant_factorization"]["hat_null_monic"],
                )
        self.assertEqual(
            [
                record["fixture_id"]
                for record in records
                if record["x_zero_direct_var1_coefficients_regular"]
            ],
            ["GR0_regular", "SGBL_constant_f", "FGCQR_schwarzschild_exterior"],
        )
        self.assertTrue(all(value is False for value in required_mhg1_nonclaims().values()))

    def test_gauge_constraint_projector_contraction_tracks_hat_metric(self) -> None:
        state = state_from_generalized_adm_pg_fixture(self.fixtures[0])
        factor_9 = gauge_constraint_propagation_principal_symbol(
            state, hat_normal_factor=Q(9)
        )
        factor_16 = gauge_constraint_propagation_principal_symbol(
            state, hat_normal_factor=Q(16)
        )
        self.assertTrue(factor_9["full_projector_contraction_exact"])
        self.assertTrue(factor_16["full_projector_contraction_exact"])
        self.assertNotEqual(
            factor_9["expected_hat_wave_factor_monic"],
            factor_16["expected_hat_wave_factor_monic"],
        )

    def test_fail_closed_auxiliary_cones_and_configuration(self) -> None:
        physical_inverse = (
            (Q(-1), Q(0), Q(0), Q(0)),
            (Q(0), Q(1), Q(0), Q(0)),
            (Q(0), Q(0), Q(1), Q(0)),
            (Q(0), Q(0), Q(0), Q(1)),
        )
        with self.assertRaisesRegex(ValueError, "greater than one"):
            auxiliary_inverse_metric(physical_inverse, Q(1))
        with self.assertRaisesRegex(ValueError, "1 < tilde"):
            modified_harmonic_certificate(
                {
                    "fixtures": self.fixtures,
                    "tilde_normal_factor": Q(9),
                    "hat_normal_factor": Q(4),
                }
            )
        with self.assertRaisesRegex(ValueError, "unexpected keys"):
            modified_harmonic_certificate(
                {
                    "fixtures": self.fixtures,
                    "tilde_normal_factor": Q(4),
                    "hat_normal_factor": Q(9),
                    "extra": True,
                }
            )

        self.assertTrue(self.certificate["fixture_records"][0]["gauge_fixing_nonzero"])


if __name__ == "__main__":
    unittest.main()
