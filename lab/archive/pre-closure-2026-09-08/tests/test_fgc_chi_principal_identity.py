from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.chi_principal_identity import (
    FormalPolynomial,
    ZERO,
    chi_principal_identity,
)


class ChiPrincipalIdentityTests(unittest.TestCase):
    def test_formal_polynomial_canonical_arithmetic(self) -> None:
        a = FormalPolynomial.atom("a")
        b = FormalPolynomial.atom("b")
        self.assertEqual(a * b - b * a, ZERO)
        self.assertEqual((a + b) * (a - b), a * a - b * b)
        self.assertEqual((Q(2) * a - a).canonical(), "1/1:a")
        with self.assertRaises(TypeError):
            FormalPolynomial.constant(True)
        with self.assertRaises(ValueError):
            FormalPolynomial.atom("")

    def test_arbitrary_jet_cross_block_and_factor_identity(self) -> None:
        certificate = chi_principal_identity()
        ledger = certificate["principal_dependency_ledger"]
        self.assertEqual(
            certificate["classification"],
            "universal_formal_uneliminated_ACT1_chi_principal_identity",
        )
        self.assertTrue(ledger["metric_and_phi_rows_have_no_chi_second_jet_dependence"])
        self.assertTrue(ledger["chi_row_has_no_metric_or_phi_second_jet_dependence"])
        self.assertTrue(ledger["chi_row_non_chi_principal_coefficients_are_formally_zero"])
        self.assertTrue(ledger["non_chi_rows_chi_principal_coefficients_are_formally_zero"])
        self.assertEqual(
            certificate["chi_row_coefficients"],
            {"chi.dtt": "1/1:g_inv_tt", "chi.dtr": "2/1:g_inv_tr", "chi.drr": "1/1:g_inv_rr"},
        )
        self.assertEqual(certificate["radial_factor_difference"], "")
        self.assertTrue(certificate["radial_chi_factor_equals_physical_metric_null_formally"])

    def test_scope_is_uneliminated_and_nonclaims_remain_closed(self) -> None:
        certificate = chi_principal_identity()
        self.assertFalse(certificate["assumptions"]["constraint_or_gauge_elimination_performed"])
        self.assertTrue(certificate["assumptions"]["arbitrary_regular_two_jets"])
        self.assertEqual(certificate["assumptions"]["arbitrary_radial_covector"], "xi_A=(-c,1); covariant orbit polynomial is also retained")
        self.assertTrue(all(value is False for value in certificate["nonclaims"].values()))
        self.assertEqual(
            certificate["covariant_orbit_principal_polynomial"],
            "1/1:g_inv_rr*xi_r*xi_r;2/1:g_inv_tr*xi_r*xi_t;1/1:g_inv_tt*xi_t*xi_t",
        )


if __name__ == "__main__":
    unittest.main()
