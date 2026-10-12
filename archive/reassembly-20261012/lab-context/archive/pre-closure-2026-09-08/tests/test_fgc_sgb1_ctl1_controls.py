from __future__ import annotations

import ast
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.sgb1_ctl1_controls import (  # noqa: E402
    SGBLControlRecord,
    decoupling_limit_linear_gb_scalar_jets,
    schwarzschild_gauss_bonnet_exact,
    sgbl_branch_controls,
    sgbl_constant_coupling_topological_control,
    sgbl_established_decoupling_control,
    sgbl_schwarzschild_vacuum_control,
    sgbl_weak_field_control,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    SGBLSourceInputs,
    sgbl_source_residual,
    sgbl_source_state,
)
from recursive_horizons.fgc.spherical_reduction import residuals  # noqa: E402


ZERO6 = (Q(0), Q(0), Q(0), Q(0), Q(0), Q(0))


class SGBLControlIdentityTests(unittest.TestCase):
    def test_weak_field_gb_vanishes_on_minkowski(self) -> None:
        record = sgbl_weak_field_control()
        self.assertEqual(record.name, "weak_field")
        self.assertTrue(record.local_identity_holds)
        self.assertEqual(record.payload["GB_constant"], 0)
        self.assertEqual(record.payload["GB_gradient"], 0)
        self.assertEqual(
            record.payload["phi_residual_constant"],
            -record.payload["potential_prime"],
        )
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)

    def test_constant_coupling_kills_metric_gb_stress(self) -> None:
        record = sgbl_constant_coupling_topological_control()
        self.assertTrue(record.payload["gb_residual_term_zero"])
        self.assertEqual(record.payload["Ricci_scalar"], 0)
        self.assertEqual(
            record.payload["GB"],
            schwarzschild_gauss_bonnet_exact(3, 4),
        )
        self.assertTrue(record.not_a_full_backreacted_solution)

    def test_schwarzschild_phi_zero_is_not_a_full_sgbl_solution(self) -> None:
        record = sgbl_schwarzschild_vacuum_control()
        gauss_bonnet = schwarzschild_gauss_bonnet_exact(3, 4)
        self.assertEqual(gauss_bonnet, Q(27, 1024))
        self.assertEqual(record.payload["phi_residual"], -Q(1, 4) * gauss_bonnet)
        self.assertTrue(record.payload["metric_residual_zero"])
        self.assertTrue(record.not_a_full_backreacted_solution)
        point = SGBLSourceInputs(
            coordinate_radius=4,
            alpha=(Q(1, 2), 0, Q(3, 16), 0, Q(-21, 128)),
            shift=(0, 0, 0, 0, 0),
            radial_metric=(2, 0, Q(-3, 4), 0, Q(39, 32)),
            areal_radius=(4, 0, 1, 0, 0),
            phi=(0, 0, 0, 0, 0),
            chi=(0, 0, 0, 0, 0),
            planck_mass=2,
            scalar_mass=3,
            quartic_coupling=Q(1, 2),
            alpha_gb=Q(-1, 4),
        )
        unredefined = residuals(sgbl_source_state(point), use_warped=True)
        self.assertEqual(unredefined["phi"], record.payload["phi_residual"])
        self.assertNotEqual(sgbl_source_residual(point, ZERO6), ZERO6)

    def test_established_decoupling_scalar_is_exact_and_not_backreacted(self) -> None:
        jets = decoupling_limit_linear_gb_scalar_jets()
        self.assertEqual(jets["value"], -Q(25, 192))
        self.assertEqual(jets["dr"], Q(37, 768))
        self.assertEqual(jets["drr"], -Q(13, 384))
        self.assertEqual(jets["box_phi_plus_alpha_gb_GB"], 0)
        self.assertEqual(jets["gauss_bonnet"], Q(27, 1024))
        record = sgbl_established_decoupling_control()
        self.assertEqual(
            record.classification,
            "decoupling_limit_linear_gb_scalar_on_schwarzschild",
        )
        self.assertIn("decoupling", record.label)
        self.assertIn("not a full backreacted solution", record.label)
        self.assertTrue(record.not_a_full_backreacted_solution)
        self.assertEqual(record.payload["box_phi_plus_alpha_gb_GB"], 0)
        self.assertEqual(record.payload["unredefined_phi_plus_potential_prime"], 0)
        self.assertFalse(record.payload["metric_residual_zero"])
        self.assertFalse(record.payload["backreacted"])
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(record.copied_fgcqr_health_evidence, False)


class SGBLControlBundleTests(unittest.TestCase):
    def test_bundle_keeps_aggregate_health_false(self) -> None:
        bundle = sgbl_branch_controls()
        self.assertEqual(
            bundle["order"],
            (
                "weak_field",
                "constant_coupling_topological",
                "schwarzschild",
                "established_sgb_decoupling",
            ),
        )
        self.assertTrue(bundle["all_local_identities_hold"])
        self.assertIs(bundle["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(bundle["holdout_authorized"], False)
        self.assertIs(bundle["FRZ1"], False)
        self.assertIs(bundle["PREF1"], False)
        self.assertTrue(bundle["established_control_is_decoupling_not_backreacted"])
        self.assertIs(SGBLControlRecord.__dataclass_params__.frozen, True)

    def test_module_does_not_import_qr_health_or_proto4(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_controls as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        forbidden = {
            "FGCQRActionParameters",
            "solve_initial_data",
            "newton_residual_limit",
            "proto1_compactness_upper_bound",
            "weak_coupling_health_certificate",
            "run_fgc_gr0_calibration_v1",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("1e-12", source)
        self.assertNotIn("HYP2", source)


if __name__ == "__main__":
    unittest.main()
