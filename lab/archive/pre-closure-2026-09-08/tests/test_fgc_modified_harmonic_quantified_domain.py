from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config
from recursive_horizons.fgc.exact_interval import interval
from recursive_horizons.fgc.modified_harmonic_quantified_domain import (
    QIFT1_ACCELERATION_ORDER,
    QIFT1_PARAMETER_ORDER,
    krawczyk_acceleration_box,
    quantified_implicit_branch_certificate,
    required_qift1_nonclaims,
)
from recursive_horizons.fgc.reference_connection import (
    flat_spherical_annulus_reference,
)


class FGCModifiedHarmonicQuantifiedDomainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        predecessor = load_config()
        cls.configuration = {
            "flat_fixture": predecessor.flat_fixture,
            "reference": flat_spherical_annulus_reference(
                radial_domain_minimum=predecessor.radial_domain_minimum
            ),
            "coordinate_radius": predecessor.coordinate_radii[
                "FGCQR_flat_reference_vacuum"
            ],
            "tilde_normal_factor": predecessor.tilde_normal_factor,
            "hat_normal_factor": predecessor.hat_normal_factor,
            "parameter_half_width": Q(1, 2**16),
            "acceleration_half_width": Q(1, 2**7),
        }
        cls.certificate = quantified_implicit_branch_certificate(cls.configuration)

    def test_full_dimensional_uniform_contraction_certificate_passes(self) -> None:
        certificate = self.certificate
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(
            certificate[
                "quantified_full_dimensional_local_implicit_branch_box_derived"
            ]
        )
        self.assertEqual(certificate["formulation"]["parameter_dimension"], 30)
        self.assertEqual(certificate["formulation"]["acceleration_dimension"], 6)
        self.assertEqual(
            certificate["formulation"]["parameter_order"], QIFT1_PARAMETER_ORDER
        )
        self.assertEqual(
            certificate["formulation"]["acceleration_order"],
            QIFT1_ACCELERATION_ORDER,
        )
        contraction = certificate["krawczyk_contraction_certificate"]
        self.assertLess(contraction["contraction_infinity_norm_upper_bound"], 1)
        self.assertGreater(
            contraction["minimum_strict_interior_inclusion_margin"], 0
        )
        self.assertTrue(
            contraction["unique_acceleration_root_for_every_declared_parameter_point"]
        )

    def test_regular_margins_and_center_jacobian_enclosure_are_strict(self) -> None:
        certificate = self.certificate
        self.assertTrue(
            certificate["verified_exact_checks"][
                "interval_acceleration_jacobian_contains_exact_center"
            ]
        )
        margins = certificate["regular_domain_ledger"]["strict_positive_margins"]
        self.assertEqual(len(margins), 7)
        self.assertTrue(all(value > 0 for value in margins.values()))
        residual = certificate["interval_residual_and_jacobian"][
            "residual_at_center_acceleration_parameter_box"
        ]
        self.assertTrue(all(value.contains_zero() for value in residual))

    def test_krawczyk_kernel_accepts_identity_toy_and_rejects_unit_bound(self) -> None:
        identity = tuple(
            tuple(Q(int(row == column)) for column in range(6))
            for row in range(6)
        )
        residual = tuple(interval(Q(-1, 100), Q(1, 100)) for _ in range(6))
        jacobian = tuple(
            tuple(interval(int(row == column)) for column in range(6))
            for row in range(6)
        )
        accepted = krawczyk_acceleration_box(
            residual_at_center_acceleration=residual,
            acceleration_jacobian_box=jacobian,
            inverse_center_jacobian=identity,
            acceleration_half_width=Q(1, 10),
        )
        self.assertEqual(accepted["contraction_infinity_norm_upper_bound"], 0)
        self.assertTrue(
            accepted["krawczyk_image_strictly_inside_acceleration_box"]
        )

        ambiguous = tuple(
            tuple(
                interval(0, 2) if row == column else interval(0)
                for column in range(6)
            )
            for row in range(6)
        )
        with self.assertRaisesRegex(ValueError, "not below one"):
            krawczyk_acceleration_box(
                residual_at_center_acceleration=residual,
                acceleration_jacobian_box=ambiguous,
                inverse_center_jacobian=identity,
                acceleration_half_width=Q(1, 10),
            )

    def test_wide_or_malformed_domains_fail_closed(self) -> None:
        malformed = dict(self.configuration)
        malformed["unknown"] = True
        with self.assertRaisesRegex(ValueError, "unexpected keys"):
            quantified_implicit_branch_certificate(malformed)

        wide = dict(self.configuration)
        wide["parameter_half_width"] = Q(2)
        with self.assertRaises((ValueError, ZeroDivisionError, TypeError)):
            quantified_implicit_branch_certificate(wide)

    def test_nonclaims_remain_false(self) -> None:
        expected = required_qift1_nonclaims()
        self.assertEqual(self.certificate["nonclaims"], expected)
        self.assertTrue(expected)
        self.assertTrue(all(value is False for value in expected.values()))


if __name__ == "__main__":
    unittest.main()
