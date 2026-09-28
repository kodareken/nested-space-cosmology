from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_reduction import load_config  # noqa: E402
from recursive_horizons.fgc.exact_linear_algebra import (  # noqa: E402
    matrix_inverse,
    matrix_mul,
    poly_eval,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    principal_matrix,
    state_from_generalized_adm_pg_fixture,
)
from recursive_horizons.fgc.spherical_symbol import (  # noqa: E402
    ZERO,
    _polynomial_matrix_multiply,
    adm_normal_principal_matrix,
    naive_fixed_gauge_comparator,
    quadratic_companion_certificate,
    quotient_symbol,
    radial_symbol_from_principal_matrix,
    required_sym1_nonclaims,
    right_gauge_generators,
    spherical_symbol_certificate,
    symbol_fixture_certificate,
)


class FGCSphericalSymbolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fixtures = load_config().configuration["fixtures"]
        cls.certificate = spherical_symbol_certificate({"fixtures": cls.fixtures})

    def test_exact_gauge_bianchi_atlas_and_scalar_factorization(self) -> None:
        self.assertTrue(self.certificate["all_declared_exact_checks_pass"])
        self.assertTrue(all(self.certificate["verified_exact_checks"].values()))
        for record in self.certificate["fixture_records"]:
            with self.subTest(fixture=record["fixture_id"]):
                self.assertTrue(record["right_gauge_null_identity_exact"])
                self.assertTrue(record["left_bianchi_null_identity_exact"])
                self.assertTrue(
                    record["quotient_atlas_scalar_determinants_agree_exactly"]
                )
                self.assertTrue(record["quotient_atlas_covers_patch_boundary"])
                self.assertTrue(record["scalar_determinant_factorization_exact"])
                self.assertTrue(record["metric_block_invertible_at_physical_roots"])
                for patch in ("rr_patch", "tt_patch"):
                    schur = record["quotient_atlas"][patch]["metric_schur"]
                    self.assertTrue(schur["schur_determinant_identity"])
                    self.assertTrue(schur["scalar_offdiagonal_numerators_zero"])

    def test_controls_share_metric_cone_and_activated_regulator_is_distinct(self) -> None:
        records = {
            record["fixture_id"]: record
            for record in self.certificate["fixture_records"]
        }
        for fixture_id in (
            "GR0_regular",
            "SGBL_constant_f",
            "FGCQR_schwarzschild_exterior",
        ):
            self.assertTrue(records[fixture_id]["regulator_cone_proportional_to_metric"])
            self.assertEqual(records[fixture_id]["cone_resultant"], 0)

        activated = records["FGCQR_activated_generic"]
        self.assertFalse(activated["regulator_cone_proportional_to_metric"])
        self.assertNotEqual(activated["cone_resultant"], 0)
        self.assertEqual(
            activated["metric_null_certificate"]["roots"],
            ({"exact": Q(-31, 42)}, {"exact": Q(59, 42)}),
        )
        for interval in activated["regulator_certificate"]["roots"]:
            self.assertIn("lower", interval)
            self.assertIn("upper", interval)
            self.assertLess(interval["lower"], interval["upper"])
            polynomial = activated["regulator_factor"]
            self.assertLess(
                poly_eval(polynomial, interval["lower"])
                * poly_eval(polynomial, interval["upper"]),
                0,
            )
        self.assertGreater(activated["regulator_certificate"]["discriminant"], 0)
        self.assertTrue(
            activated["regulator_certificate"]["symmetrizer_positive_definite"]
        )
        self.assertTrue(
            activated["regulator_certificate"]["product_is_symmetric"]
        )

    def test_adm_candidate_constraints_kinetic_block_and_rejected_shortcut(self) -> None:
        for fixture in self.fixtures:
            with self.subTest(fixture=fixture["fixture_id"]):
                data = adm_normal_principal_matrix(fixture)
                self.assertTrue(
                    data["candidate_constraints_free_of_normal_accelerations"]
                )
                self.assertTrue(
                    data["candidate_constraints_free_of_gauge_source_second_jets"]
                )
                self.assertNotEqual(data["kinetic_determinant"], 0)

        gr = adm_normal_principal_matrix(self.fixtures[0])
        self.assertEqual(
            [
                (index, value)
                for index, value in enumerate(gr["candidate_constraint_rows"][0])
                if value
            ],
            [(5, Q(-8))],
        )
        self.assertEqual(
            [
                (index, value)
                for index, value in enumerate(gr["candidate_constraint_rows"][1])
                if value
            ],
            [(4, Q(8))],
        )
        naive = naive_fixed_gauge_comparator(self.fixtures[0])
        self.assertEqual(naive["zero_speed_geometric_multiplicity"], 1)
        self.assertEqual(naive["zero_speed_generalized_multiplicity"], 4)
        self.assertTrue(naive["defective_zero_speed_sector"])
        self.assertTrue(naive["strong_hyperbolicity_rejected"])

    def test_naive_comparator_converts_normal_to_coordinate_time_at_nonzero_shift(
        self,
    ) -> None:
        fixture = self.fixtures[-1]
        shift = fixture["state"]["shift"]
        shift = Q(shift)
        data = adm_normal_principal_matrix(fixture)
        matrix = data["matrix"][2:]
        velocity_perp = (0, 3, 6, 9)
        velocity_r = (1, 4, 7, 10)
        gradient_r = (2, 5, 8, 11)
        kinetic = tuple(
            tuple(row[column] for column in velocity_perp) for row in matrix
        )
        velocity = tuple(
            tuple(row[column] for column in velocity_r) for row in matrix
        )
        gradient = tuple(
            tuple(row[column] for column in gradient_r) for row in matrix
        )
        kinetic_inverse = matrix_inverse(kinetic)
        expected_velocity = tuple(
            tuple(
                entry - (shift if row == column else Q(0))
                for column, entry in enumerate(source_row)
            )
            for row, source_row in enumerate(matrix_mul(kinetic_inverse, velocity))
        )
        expected_gradient = matrix_mul(kinetic_inverse, gradient)

        principal = naive_fixed_gauge_comparator(fixture)["normal_principal_matrix"]
        self.assertEqual(
            tuple(tuple(row[:4]) for row in principal[:4]),
            expected_velocity,
        )
        self.assertEqual(
            tuple(tuple(row[4:]) for row in principal[:4]),
            expected_gradient,
        )
        self.assertEqual(
            tuple(tuple(row[4:]) for row in principal[4:]),
            tuple(
                tuple(-shift if row == column else Q(0) for column in range(4))
                for row in range(4)
            ),
        )

    def test_fail_closed_mutation_patch_and_quadratic_checks(self) -> None:
        fixture = self.fixtures[-1]
        state = state_from_generalized_adm_pg_fixture(fixture)
        data = principal_matrix(state)
        rows = [list(row) for row in data["matrix"]]
        rows[0][0] += 1
        mutated = dict(data)
        mutated["matrix"] = tuple(tuple(row) for row in rows)
        symbol = radial_symbol_from_principal_matrix(mutated)
        residual = _polynomial_matrix_multiply(symbol, right_gauge_generators())
        self.assertTrue(any(entry != ZERO for row in residual for entry in row))

        with self.assertRaisesRegex(ValueError, "row_patch"):
            quotient_symbol(symbol, row_patch="unknown")
        with self.assertRaisesRegex(ValueError, "positive discriminant"):
            quadratic_companion_certificate((Q(1), Q(0), Q(1)))
        self.assertTrue(all(value is False for value in required_sym1_nonclaims().values()))
        with self.assertRaisesRegex(ValueError, "accepts only fixtures"):
            spherical_symbol_certificate({"fixtures": self.fixtures, "extra": 1})


if __name__ == "__main__":
    unittest.main()
