from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution.covariant_principal_health import (  # noqa: E402
    second_order_coefficient_matrices,
)
from recursive_horizons.fgc.modified_harmonic import modified_harmonic_symbol  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_principal import (  # noqa: E402
    ANGULAR_DIRECTION,
    FIELD_ORDER_12,
    POINT_SPECTRUM_CLASSIFICATION,
    RADIAL_DIRECTION,
    SGBLPrincipalPointFacts,
    sgbl_covariant_principal_background,
    sgbl_principal_point_facts,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    HAT_NORMAL_FACTOR,
    SGBLSourceInputs,
    TILDE_NORMAL_FACTOR,
    sgbl_source_solve,
    sgbl_source_state,
)


def _source_point() -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(2, Q(1, 3), Q(-1, 2), Q(1, 5), Q(-1, 7)),
        shift=(Q(1, 3), Q(-2, 5), Q(1, 4), Q(-1, 6), Q(2, 9)),
        radial_metric=(Q(3, 2), Q(1, 8), Q(-1, 3), Q(1, 7), Q(-1, 4)),
        areal_radius=(4, Q(-1, 5), Q(3, 2), Q(1, 9), Q(-2, 5)),
        phi=(Q(1, 2), Q(2, 3), Q(-4, 5), Q(1, 4), Q(-1, 8)),
        chi=(Q(-3, 4), Q(1, 2), Q(3, 5), Q(-2, 7), Q(1, 6)),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=Q(-1, 4),
    )


def _solved_state():
    point = _source_point()
    return sgbl_source_state(point, sgbl_source_solve(point))


def _spherical_embedding() -> np.ndarray:
    embedding = np.zeros((12, 6))
    for column, row in enumerate((0, 1, 4)):
        embedding[row, column] = 1.0
    embedding[7, 3] = embedding[9, 3] = 1.0
    embedding[10, 4] = embedding[11, 5] = 1.0
    return embedding


def _companion(a_matrix, b_matrix, c_matrix) -> np.ndarray:
    size = a_matrix.shape[0]
    return np.block([
        [np.zeros((size, size)), np.eye(size)],
        [-np.linalg.solve(a_matrix, c_matrix), -np.linalg.solve(a_matrix, b_matrix)],
    ])


class SGBLPrincipalAdapterTests(unittest.TestCase):
    def test_linear_branch_background_has_the_declared_action_derivatives(self) -> None:
        state = _solved_state()
        background = sgbl_covariant_principal_background(state)
        self.assertEqual(background.effective_planck_squared, 4.0)
        self.assertEqual(background.effective_planck_prime, 0.0)
        self.assertEqual(background.gb_coupling_prime, -0.25)
        self.assertGreater(float(np.max(np.abs(background.riemann_lower))), 0.0)
        self.assertGreater(float(np.max(np.abs(background.hessian_gb_lower))), 0.0)

    def test_radial_restriction_matches_independent_exact_sgbl_symbol(self) -> None:
        state = _solved_state()
        background = sgbl_covariant_principal_background(state)
        a_matrix, b_matrix, c_matrix = second_order_coefficient_matrices(
            background, (1.0, 0.0, 0.0)
        )
        embedding = _spherical_embedding()
        covariant_roots = np.linalg.eigvals(_companion(
            embedding.T @ a_matrix @ embedding,
            embedding.T @ b_matrix @ embedding,
            embedding.T @ c_matrix @ embedding,
        ))

        exact_symbol = modified_harmonic_symbol(
            state,
            tilde_normal_factor=TILDE_NORMAL_FACTOR,
            hat_normal_factor=HAT_NORMAL_FACTOR,
        )
        exact_a = np.asarray([
            [float(entry[2] if len(entry) > 2 else 0) for entry in row]
            for row in exact_symbol
        ])
        exact_b = np.asarray([
            [float(entry[1] if len(entry) > 1 else 0) for entry in row]
            for row in exact_symbol
        ])
        exact_c = np.asarray([
            [float(entry[0] if entry else 0) for entry in row]
            for row in exact_symbol
        ])
        coordinate_speeds = np.linalg.eigvals(_companion(exact_a, exact_b, exact_c))

        orbit = np.asarray([
            [float(state.h_tt.value), float(state.h_tr.value)],
            [float(state.h_tr.value), float(state.h_rr.value)],
        ])
        inverse = np.linalg.inv(orbit)
        lapse = 1.0 / (-inverse[0, 0]) ** 0.5
        normal = -lapse * inverse[:, 0]
        radial_metric = float(state.h_rr.value) ** 0.5
        expected = radial_metric * (normal[1] - normal[0] * coordinate_speeds)
        np.testing.assert_allclose(
            np.sort_complex(covariant_roots),
            np.sort_complex(expected),
            rtol=0.0,
            atol=2.0e-12,
        )

    def test_non_sgbl_and_action_drift_refuse_before_symbol_construction(self) -> None:
        state = _solved_state()
        gr0 = replace(state, branch="GR-0", alpha=Q(0))
        with self.assertRaisesRegex(ValueError, "linear branch"):
            sgbl_covariant_principal_background(gr0)
        with self.assertRaises(TypeError):
            sgbl_covariant_principal_background(object())  # type: ignore[arg-type]

    def test_point_background_does_not_claim_interval_or_health(self) -> None:
        background = sgbl_covariant_principal_background(_solved_state())
        self.assertFalse(hasattr(background, "SGBL_branch_owned_and_healthy"))
        self.assertFalse(hasattr(background, "interval_invertibility"))
        self.assertFalse(hasattr(background, "continuous_no_initial_trapped_sphere"))

    def test_point_facts_retain_tensors_and_keep_health_closed(self) -> None:
        state = _solved_state()
        facts = sgbl_principal_point_facts(state)
        background = facts.background
        self.assertEqual(background.effective_planck_squared, 4.0)
        self.assertEqual(background.effective_planck_prime, 0.0)
        self.assertEqual(background.gb_coupling_prime, -0.25)
        self.assertTrue(facts.linear_branch_background_has_declared_derivatives)
        self.assertEqual(FIELD_ORDER_12[-2:], ("phi", "chi"))
        self.assertTrue(facts.ten_metric_polarizations_plus_phi_chi_retained)
        self.assertTrue(facts.radial_restriction_agrees_with_spherical_symbol)
        self.assertTrue(facts.angular_spatial_symbol_differs_from_radial)
        self.assertTrue(facts.reconstructed_symbol_agrees_at_mixed_covector)
        self.assertTrue(facts.finite_coefficient_tensors)
        self.assertTrue(facts.binary64_uncertainty_visible)
        self.assertEqual(
            facts.radial_spectrum.spatial_covector, RADIAL_DIRECTION
        )
        self.assertEqual(
            facts.angular_spectrum.spatial_covector, ANGULAR_DIRECTION
        )
        self.assertEqual(
            facts.radial_spectrum.classification, POINT_SPECTRUM_CLASSIFICATION
        )
        self.assertGreater(
            float(np.linalg.norm(
                facts.angular_spectrum.spatial_symbol
                - facts.radial_spectrum.spatial_symbol
            )),
            1.0,
        )
        self.assertLess(
            facts.radial_spectrum.maximum_abs_imaginary_part, 1.0e-12
        )
        self.assertIs(facts.strongly_hyperbolic, False)
        self.assertIs(facts.interval_invertible, False)
        self.assertIs(facts.SGBL_branch_owned_and_healthy, False)
        self.assertIs(
            facts.quantitative_all_covector_weak_coupling_health_envelope_passed,
            False,
        )
        self.assertTrue(facts.kr_sufficient_envelope_is_not_this_certificate)
        self.assertEqual(facts.cone_certificate, "unqualified")
        self.assertEqual(
            facts.coefficient_tensor_sha256,
            SGBLPrincipalPointFacts(
                background=facts.background,
                tilde_normal_factor=facts.tilde_normal_factor,
                hat_normal_factor=facts.hat_normal_factor,
                coefficient_tensors=facts.coefficient_tensors,
                coefficient_tensor_sha256=facts.coefficient_tensor_sha256,
                deformation=facts.deformation,
                radial_spectrum=facts.radial_spectrum,
                angular_spectrum=facts.angular_spectrum,
                mixed_spectrum=facts.mixed_spectrum,
                radial_restriction_agrees_with_spherical_symbol=True,
                angular_spatial_symbol_differs_from_radial=True,
                reconstructed_symbol_agrees_at_mixed_covector=True,
                action_hessian_symmetry_defect=facts.action_hessian_symmetry_defect,
                gauge_identity_defect=facts.gauge_identity_defect,
                finite_coefficient_tensors=True,
            ).coefficient_tensor_sha256,
        )
        with self.assertRaises(TypeError):
            replace(facts, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(TypeError):
            replace(facts, strongly_hyperbolic=True)
        with self.assertRaises(ValueError):
            facts.radial_spectrum.spatial_symbol[0, 0] = 1.0
        with self.assertRaises(ValueError):
            facts.coefficient_tensors.time_time[0, 0] = 1.0

    def test_module_does_not_import_qr_health_or_src1_solve(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_principal as owner
        import ast

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
        forbidden = {
            "background_from_spherical_state",
            "weak_coupling_health_certificate",
            "canonical_health_monitor_values",
            "solve_accelerations",
            "quantified_implicit_branch_certificate",
            "FGCQRActionParameters",
            "solve_initial_data",
        }
        self.assertTrue(forbidden.isdisjoint(imported))


if __name__ == "__main__":
    unittest.main()
