from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.metric_null_observable import (  # noqa: E402
    SphericalADMGeometryJet,
    metric_null_raychaudhuri_point_certificate,
)
from recursive_horizons.fgc.spherical_reduction import Jet2  # noqa: E402


class MetricNullObservableTests(unittest.TestCase):
    def test_minkowski_spherical_affine_focusing_identity(self) -> None:
        geometry = SphericalADMGeometryJet(
            lapse=Jet2.constant(1),
            shift=Jet2.constant(0),
            radial_scale=Jet2.constant(1),
            areal_radius=Jet2(4, dr=1),
        )
        output = metric_null_raychaudhuri_point_certificate(
            geometry, branch="outgoing", affine_scale=Jet2.constant(1)
        )
        self.assertEqual(output["null_frame"]["round_sphere_classification"], "normal")
        self.assertEqual(output["raychaudhuri"]["theta"], Q(1, 2))
        self.assertEqual(output["raychaudhuri"]["complete_rhs"], Q(-1, 8))
        self.assertFalse(output["raychaudhuri"]["locally_defocusing_at_point"])
        self.assertTrue(all(value is False for value in output["nonclaims"].values()))

    def test_contracting_flat_de_sitter_control_is_trapped_not_defocusing(self) -> None:
        geometry = SphericalADMGeometryJet(
            lapse=Jet2.constant(1),
            shift=Jet2.constant(0),
            radial_scale=Jet2(1, dt=-1, dtt=1),
            areal_radius=Jet2(2, dt=-2, dr=1, dtt=2, dtr=-1),
        )
        output = metric_null_raychaudhuri_point_certificate(
            geometry,
            branch="outgoing",
            affine_scale=Jet2(1, dt=1, dtt=1),
        )
        self.assertEqual(output["null_frame"]["round_sphere_classification"], "trapped")
        self.assertLess(output["raychaudhuri"]["complete_rhs"], 0)
        self.assertEqual(output["raychaudhuri"]["minus_R_ab_k_a_k_b"], 0)

    def test_synthetic_negative_Rkk_point_exercises_positive_full_rhs(self) -> None:
        geometry = SphericalADMGeometryJet(
            lapse=Jet2.constant(1),
            shift=Jet2.constant(0),
            radial_scale=Jet2.constant(1),
            areal_radius=Jet2(4, dt=-1, dr=1, dtt=1),
        )
        output = metric_null_raychaudhuri_point_certificate(
            geometry, branch="outgoing", affine_scale=Jet2.constant(1)
        )
        self.assertEqual(output["raychaudhuri"]["theta"], 0)
        self.assertEqual(output["raychaudhuri"]["minus_R_ab_k_a_k_b"], Q(1, 2))
        self.assertEqual(output["raychaudhuri"]["complete_rhs"], Q(1, 2))
        self.assertTrue(output["raychaudhuri"]["locally_defocusing_at_point"])

    def test_nonaffine_scale_and_auxiliary_claims_fail_closed(self) -> None:
        geometry = SphericalADMGeometryJet(
            lapse=Jet2.constant(1),
            shift=Jet2.constant(0),
            radial_scale=Jet2.constant(1),
            areal_radius=Jet2(4, dr=1),
        )
        with self.assertRaisesRegex(ValueError, "not affine"):
            metric_null_raychaudhuri_point_certificate(
                geometry,
                branch="outgoing",
                affine_scale=Jet2(1, dt=1),
            )
        with self.assertRaisesRegex(ValueError, "branch"):
            metric_null_raychaudhuri_point_certificate(
                geometry,
                branch="hat",
                affine_scale=Jet2.constant(1),
            )


if __name__ == "__main__":
    unittest.main()
