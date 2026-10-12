from __future__ import annotations

from math import cos, sin
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.covariant_principal_health import (  # noqa: E402
    FIELD_COUNT,
    background_from_spherical_state,
    coefficient_tensor_symbol,
    complete_modified_harmonic_principal_symbol,
    esf_background,
    first_order_principal_matrix,
    frobenius_normalized_first_order_matrix,
    frobenius_normalized_principal_symbol,
    normalized_first_order_matrix_from_coefficients,
    point_spectrum_diagnostics,
    principal_coefficient_tensors,
    second_order_coefficient_matrices,
    ungauged_action_principal_symbol,
)
from recursive_horizons.fgc.modified_harmonic import (  # noqa: E402
    modified_harmonic_symbol,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from scripts.reproduce_fgc_hyp1_dom1_qift1 import (  # noqa: E402
    load_config as load_qift1_config,
)


class FGCCovariantPrincipalHealthTests(unittest.TestCase):
    def test_flat_esf_spectrum_has_all_expected_gauge_and_physical_multiplicities(self) -> None:
        background = esf_background()
        expected = np.asarray(
            [-1.0] * 4
            + [-0.5] * 4
            + [-1.0 / 3.0] * 4
            + [1.0 / 3.0] * 4
            + [0.5] * 4
            + [1.0] * 4
        )
        for angle in (0.0, 0.173, 0.731, 1.217):
            with self.subTest(angle=angle):
                diagnostics = point_spectrum_diagnostics(
                    background, (cos(angle), sin(angle), 0.0)
                )
                eigenvalues = np.asarray(diagnostics["eigenvalues"])
                self.assertLess(diagnostics["maximum_abs_imaginary_part"], 1.0e-12)
                np.testing.assert_allclose(
                    np.sort(eigenvalues.real), expected, rtol=0.0, atol=2.0e-12
                )
                self.assertGreater(
                    diagnostics["kinetic"]["kinetic_smallest_singular_value"],
                    0.0,
                )

    def test_ungauged_action_hessian_and_independent_chi_block_are_symmetric(self) -> None:
        background = esf_background()
        covector = (0.37, 0.8, -0.3, 0.4)
        star = ungauged_action_principal_symbol(background, covector)
        self.assertEqual(star.shape, (FIELD_COUNT, FIELD_COUNT))
        np.testing.assert_allclose(star, star.T, rtol=0.0, atol=0.0)
        self.assertTrue(np.all(star[11, :11] == 0.0))
        self.assertTrue(np.all(star[:11, 11] == 0.0))
        xi_squared = -covector[0] ** 2 + sum(value * value for value in covector[1:])
        self.assertEqual(star[11, 11], xi_squared)

    def test_gauge_extension_changes_only_metric_principal_block(self) -> None:
        background = esf_background()
        covector = (0.21, 1.0, 0.0, 0.0)
        star = ungauged_action_principal_symbol(background, covector)
        complete = complete_modified_harmonic_principal_symbol(background, covector)
        difference = complete - star
        self.assertGreater(float(np.max(np.abs(difference[:10, :10]))), 0.0)
        self.assertTrue(np.all(difference[10:, :] == 0.0))
        self.assertTrue(np.all(difference[:, 10:] == 0.0))

    def test_fail_closed_input_contract(self) -> None:
        background = esf_background()
        with self.assertRaisesRegex(ValueError, "unit norm"):
            point_spectrum_diagnostics(background, (2.0, 0.0, 0.0))
        with self.assertRaisesRegex(ValueError, "1<tilde<hat"):
            point_spectrum_diagnostics(
                background,
                (1.0, 0.0, 0.0),
                tilde_normal_factor=9.0,
                hat_normal_factor=4.0,
            )

    def test_extracted_covector_tensor_reassembles_symbol_and_companion(self) -> None:
        background = esf_background()
        coefficients = principal_coefficient_tensors(background)
        covector = (0.27, 0.3, -0.4, (0.75) ** 0.5)
        direct = frobenius_normalized_principal_symbol(
            complete_modified_harmonic_principal_symbol(background, covector)
        )
        reconstructed = coefficient_tensor_symbol(coefficients, covector)
        np.testing.assert_allclose(
            reconstructed, direct, rtol=0.0, atol=2.0e-14
        )

        direction = np.asarray(covector[1:])
        direction /= np.linalg.norm(direction)
        raw_companion, _raw_kinetic = first_order_principal_matrix(
            background, direction
        )
        normalized_direct = frobenius_normalized_first_order_matrix(
            raw_companion
        )
        normalized_reconstructed, _normalized_kinetic = (
            normalized_first_order_matrix_from_coefficients(
                coefficients, direction
            )
        )
        np.testing.assert_allclose(
            normalized_reconstructed,
            normalized_direct,
            rtol=0.0,
            atol=2.0e-13,
        )

    def test_covariant_radial_restriction_matches_independent_exact_symbol(self) -> None:
        exact_state = activated_compatible_state(
            load_qift1_config().flat_fixture
        )["state"]
        background = background_from_spherical_state(exact_state)
        a_matrix, b_matrix, c_matrix = second_order_coefficient_matrices(
            background, (1.0, 0.0, 0.0)
        )

        # Spherical perturbations retain h_00, h_01, h_11, the equal angular
        # pair h_22=h_33, phi, and chi.  T^T P T is the action-conjugate
        # restriction; invertible row/field scalings do not change its roots.
        embedding = np.zeros((12, 6))
        for column, row in enumerate((0, 1, 4)):
            embedding[row, column] = 1.0
        embedding[7, 3] = embedding[9, 3] = 1.0
        embedding[10, 4] = embedding[11, 5] = 1.0
        reduced_a = embedding.T @ a_matrix @ embedding
        reduced_b = embedding.T @ b_matrix @ embedding
        reduced_c = embedding.T @ c_matrix @ embedding
        reduced = np.block(
            [
                [np.zeros((6, 6)), np.eye(6)],
                [
                    -np.linalg.solve(reduced_a, reduced_c),
                    -np.linalg.solve(reduced_a, reduced_b),
                ],
            ]
        )
        covariant_roots = np.linalg.eigvals(reduced)

        exact_symbol = modified_harmonic_symbol(
            exact_state,
            tilde_normal_factor=4,
            hat_normal_factor=9,
        )
        exact_a = np.asarray(
            [
                [float(entry[2] if len(entry) > 2 else 0) for entry in row]
                for row in exact_symbol
            ]
        )
        exact_b = np.asarray(
            [
                [float(entry[1] if len(entry) > 1 else 0) for entry in row]
                for row in exact_symbol
            ]
        )
        exact_c = np.asarray(
            [
                [float(entry[0] if entry else 0) for entry in row]
                for row in exact_symbol
            ]
        )
        exact_companion = np.block(
            [
                [np.zeros((6, 6)), np.eye(6)],
                [
                    -np.linalg.solve(exact_a, exact_c),
                    -np.linalg.solve(exact_a, exact_b),
                ],
            ]
        )
        coordinate_speeds = np.linalg.eigvals(exact_companion)
        radial_metric = float(exact_state.h_rr.value) ** 0.5
        # The covariant symbol uses xi=(xi_0,n_hat); the exact radial symbol
        # uses xi=(-c,1).  Hence xi_0=-lambda*c on this zero-shift unit-lapse
        # fixture.
        expected_covariant_roots = -radial_metric * coordinate_speeds
        np.testing.assert_allclose(
            np.sort_complex(covariant_roots),
            np.sort_complex(expected_covariant_roots),
            rtol=0.0,
            atol=2.0e-13,
        )


if __name__ == "__main__":
    unittest.main()
