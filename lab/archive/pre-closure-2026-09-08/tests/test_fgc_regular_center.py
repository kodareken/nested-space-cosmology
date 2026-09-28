from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.regular_center import (  # noqa: E402
    LaurentSeries,
    analytic_defect_ledger,
    annular_regularized_evaluation,
    first_grid_point_convergence_certificate,
    regular_center_series_certificate,
    validate_regular_profile,
)
from scripts.reproduce_fgc_ctr1_reg1 import load_config  # noqa: E402


class FGCRegularCenterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        configuration = load_config()
        cls.minkowski_profile = configuration["profiles"]["minkowski"]
        cls.nontrivial_profile = configuration["profiles"]["nontrivial"]
        cls.minkowski = regular_center_series_certificate(cls.minkowski_profile)
        cls.nontrivial = regular_center_series_certificate(cls.nontrivial_profile)

    def test_exact_laurent_arithmetic_and_reference_quotient(self) -> None:
        radius = LaurentSeries.monomial(1, 1)
        even = LaurentSeries.from_mapping({0: Q(2), 2: Q(3), 4: -Q(1, 5)})
        product = (radius * even) / radius
        for power in range(-4, 5):
            self.assertEqual(product.coefficient(power), even.coefficient(power))
        reciprocal = even.reciprocal()
        identity = even * reciprocal
        self.assertEqual(identity.coefficient(0), 1)
        self.assertTrue(all(identity.coefficient(power) == 0 for power in range(1, 5)))

    def test_negative_guard_edge_fails_closed_instead_of_truncating(self) -> None:
        edge = LaurentSeries.monomial(1, -12)
        with self.assertRaisesRegex(ArithmeticError, "negative-power guard edge"):
            edge.derivative()
        with self.assertRaisesRegex(ArithmeticError, "negative-power guard edge"):
            edge * LaurentSeries.monomial(1, -1)

    def test_minkowski_and_nontrivial_exact_center_controls(self) -> None:
        self.assertEqual(self.minkowski["center_limits"], (Q(0),) * 6)
        self.assertEqual(self.minkowski["gauge_center_limits"], (Q(0),) * 2)
        self.assertTrue(all(not value for value in self.minkowski["curvature_series"].values()))
        self.assertTrue(all(self.nontrivial["checks"].values()))
        self.assertEqual(self.nontrivial["center_limits"][2], self.nontrivial["center_limits"][3])
        self.assertEqual(
            self.nontrivial["curvature_series"]["Kretschmann"][0],
            Q(18910513, 17280000),
        )
        reference = self.nontrivial["reference_regularization"]
        self.assertTrue(reference["R_r_over_R_minus_one_over_r_equals_A_r_over_A"])
        self.assertTrue(reference["radial_derivative_identity_exact"])
        self.assertTrue(reference["all_connection_differences_have_no_certified_negative_powers"])

    def test_first_grid_points_converge_to_independent_formal_limits(self) -> None:
        certificate = first_grid_point_convergence_certificate(
            self.nontrivial_profile,
            formal_certificate=self.nontrivial,
        )
        self.assertTrue(certificate["all_regular_equation_and_gauge_controls_passed"])
        self.assertTrue(certificate["not_a_PDE_discretization_convergence_claim"])
        for record in certificate["records"]:
            self.assertTrue(record["at_least_second_order_control_passed"])
            for ratio in record["coarse_to_fine_error_ratios"]:
                if ratio is not None:
                    self.assertGreaterEqual(ratio, Q(15, 4))

    def test_conical_parity_and_time_flatness_defects_fail_closed(self) -> None:
        conical = deepcopy(self.nontrivial_profile)
        conical["fields"]["A"]["value"][0] = Q(11, 10)
        self.assertEqual(
            analytic_defect_ledger(conical)["Ricci_scalar_r_minus_2_from_conical_defect"],
            -Q(42, 121),
        )
        with self.assertRaisesRegex(ValueError, "elementary flatness"):
            validate_regular_profile(conical)

        odd = deepcopy(self.nontrivial_profile)
        odd["fields"]["phi"]["value"][1] = Q(1, 17)
        self.assertEqual(
            analytic_defect_ledger(odd)["box_phi_r_minus_1_from_odd_scalar_term"],
            Q(2, 17),
        )
        with self.assertRaisesRegex(ValueError, "even centre parity"):
            validate_regular_profile(odd)

        time_defect = deepcopy(self.nontrivial_profile)
        time_defect["fields"]["A"]["dt"][0] = Q(1, 39)
        with self.assertRaisesRegex(ValueError, "elementary flatness"):
            validate_regular_profile(time_defect)

    def test_annular_control_never_evaluates_reference_at_center(self) -> None:
        with self.assertRaisesRegex(ValueError, "0 < reference minimum < radius"):
            annular_regularized_evaluation(
                self.minkowski_profile,
                coordinate_radius=Q(1, 16),
                reference_radial_minimum=Q(1, 16),
            )


if __name__ == "__main__":
    unittest.main()
