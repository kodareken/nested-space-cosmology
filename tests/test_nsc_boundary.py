"""Focused tests for the first-order radial Dirac boundary maps.

These tests do not rebuild the full published record. All-field reproduction
is `python3 scripts/check_nsc_boundary_response.py --check`.
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.linalg import eigh


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from recursive_horizons.nsc_boundary import (  # noqa: E402
    DEFAULT_RESULT,
    PRIMARY_KAPPA,
    PRIMARY_RADIUS,
    PROBE_ENERGIES,
    SCHEMA,
    SIDE_INTERVAL,
    _canonical_spectrum,
    compare_discrete_continuum,
    dense_svd_half_map_control,
    discrete_half_map,
    doubling_control,
    exact_relations,
    inverse_radius_scaling,
    joined_operator,
    naive_centered_operator,
    oriented_dtn_jump,
    real_energy_probability_current,
    schur_resolvent_test,
    staggered_forward,
    superpotential,
)


class ExactDictionaryTests(unittest.TestCase):
    def test_partner_and_current_residuals_are_zero(self) -> None:
        algebra = exact_relations()
        self.assertEqual(algebra["partners"]["plus_residual"], "0")
        self.assertEqual(algebra["partners"]["minus_residual"], "0")
        self.assertEqual(algebra["operator"]["factorization_residual"], "0")
        self.assertEqual(algebra["current"]["residual"], "0")
        self.assertTrue(algebra["current"]["identity_matching_preserves_current"])
        self.assertFalse(
            algebra["current"]["historical_sigma3_unitary_preserves_this_current"]
        )
        self.assertTrue(algebra["current"]["SO2_family_preserves_current"])
        self.assertFalse(algebra["current"]["abel_wronskian_is_probability_current"])

    def test_orientation_and_oriented_jump_identity(self) -> None:
        orientation = exact_relations()["orientation"]
        self.assertIn("rho_in_[0, R]", orientation["parent_interval"])
        self.assertIn("rho_in_[-R, 0]", orientation["child_interval"])
        self.assertEqual(orientation["oriented_jump_residual"], "0")
        self.assertTrue(orientation["w0_cancels_in_the_jump"])
        self.assertIn("minus_d_rho", orientation["parent_outward_normal"])
        self.assertIn("plus_d_rho", orientation["child_outward_normal"])

    def test_dtn_is_energy_times_weyl_and_is_not_phi(self) -> None:
        algebra = exact_relations()
        relation = algebra["weyl_calderon_versus_DtN"]
        self.assertEqual(relation["identity_residual"], "0")
        self.assertFalse(relation["spectral_DtN_is_Phi"])
        self.assertIn("not_the_continuum_Calderon_map", relation["staggered_edge_trace"])

    def test_schur_order_and_one_over_omega_convention(self) -> None:
        conventions = exact_relations()["resolvent_conventions"]
        self.assertEqual(
            conventions["noncommuting_block_Schur_residual"], "Matrix([[0, 0], [0, 0]])"
        )
        self.assertTrue(conventions["reversed_multiplication_differs"])
        self.assertFalse(conventions["conventions_may_be_mixed"])


class OperatorTests(unittest.TestCase):
    def test_staggered_hamiltonian_is_hermitian(self) -> None:
        assembled = joined_operator(PRIMARY_RADIUS, 40, PRIMARY_KAPPA)
        self.assertEqual(assembled["hermitian_residual"], 0.0)
        values = eigh(assembled["H"], eigvals_only=True)
        self.assertTrue(np.allclose(values.imag, 0.0))

    def test_superpotential_decays(self) -> None:
        near = float(superpotential(0.0, PRIMARY_KAPPA))
        far = float(superpotential(50.0, PRIMARY_KAPPA))
        self.assertAlmostEqual(near, 1.0)
        self.assertLess(far, 0.03)

    def test_incidence_is_node_to_edge(self) -> None:
        nodes = np.linspace(0.0, 1.0, 5)
        forward = staggered_forward(nodes, PRIMARY_KAPPA)
        self.assertEqual(forward.shape, (4, 5))

    def test_canonical_spectrum_uses_signed_bands(self) -> None:
        values = np.array([-0.5, -0.5, -1e-16, 0.5, 0.5, 0.8])
        record = _canonical_spectrum(values, 4.0, 0.0)
        self.assertEqual(record["index_kernel_count"], 1)
        self.assertEqual(record["first_positive"], 0.5)
        self.assertEqual(record["first_negative"], -0.5)
        self.assertEqual(record["positive_ascending"][0], 0.5)
        self.assertNotIn("lowest_values", record)


class MapAndResolventTests(unittest.TestCase):
    def test_parent_is_positive_rho_and_child_is_negative(self) -> None:
        parent = discrete_half_map(
            PRIMARY_RADIUS, 40, PROBE_ENERGIES[0], PRIMARY_KAPPA, "parent"
        )
        child = discrete_half_map(
            PRIMARY_RADIUS, 40, PROBE_ENERGIES[0], PRIMARY_KAPPA, "child"
        )
        self.assertGreater(parent["edge_rho"], 0.0)
        self.assertLess(child["edge_rho"], 0.0)
        self.assertEqual(parent["u_throat"], 1.0 + 0.0j)
        self.assertEqual(child["u_throat"], 1.0 + 0.0j)
        self.assertIn("[0, R]", SIDE_INTERVAL["parent"])
        self.assertIn("[-R, 0]", SIDE_INTERVAL["child"])

    def test_discrete_maps_track_solve_ivp_at_edge_and_throat(self) -> None:
        comparison = compare_discrete_continuum(
            PRIMARY_RADIUS, 80, PROBE_ENERGIES[0], PRIMARY_KAPPA, "parent"
        )
        self.assertLess(comparison["edge_trace_relative_residual"], 5.0e-3)
        self.assertLess(comparison["throat_map_relative_residual"], 5.0e-3)
        self.assertLess(comparison["linear_solve_residual"], 1.0e-12)
        self.assertFalse(comparison["staggered_edge_projector_is_continuum_Calderon"])
        self.assertFalse(comparison["abel_wronskian_is_probability_current"])

    def test_oriented_dtn_jump_cancels_w0(self) -> None:
        jump = oriented_dtn_jump(
            PRIMARY_RADIUS, 40, PROBE_ENERGIES[0], PRIMARY_KAPPA
        )
        self.assertLess(jump["algebraic_jump_residual"], 1.0e-12)
        self.assertLess(jump["continuum_jump_identity_residual"], 1.0e-12)
        self.assertTrue(jump["w0_cancels_in_the_jump"])
        self.assertGreater(jump["parent"]["edge_rho"], 0.0)
        self.assertLess(jump["child"]["edge_rho"], 0.0)

    def test_sparse_lu_matches_dense_svd_control(self) -> None:
        energy = PROBE_ENERGIES[0]
        for side in ("parent", "child"):
            sparse = discrete_half_map(PRIMARY_RADIUS, 20, energy, PRIMARY_KAPPA, side)
            dense = dense_svd_half_map_control(
                PRIMARY_RADIUS, 20, energy, PRIMARY_KAPPA, side
            )
            residual = abs(sparse["weyl_m_edge"] - dense) / abs(dense)
            self.assertLess(residual, 1.0e-10)

    def test_joined_schur_sparse_and_imaginary_resolvent_sign(self) -> None:
        result = schur_resolvent_test(
            PRIMARY_RADIUS, 20, PROBE_ENERGIES[0], PRIMARY_KAPPA, dense_control=True
        )
        self.assertEqual(result["solver"], "sparse_lu")
        self.assertFalse(result["dense_inversion_is_primary"])
        self.assertLess(result["direct_versus_eliminated_residual"], 1.0e-10)
        self.assertLess(result["dense_versus_sparse_residual"], 1.0e-10)
        self.assertTrue(result["imag_resolvent_sign_matches_energy"])

    def test_naive_centered_is_not_the_continuum_proof(self) -> None:
        control = doubling_control(PRIMARY_RADIUS, 40, PRIMARY_KAPPA)
        self.assertFalse(control["continuum_proof_uses_naive_centered"])
        naive = naive_centered_operator(PRIMARY_RADIUS, 40, PRIMARY_KAPPA)
        self.assertEqual(naive["method"], "naive_collocated_centered_derivative")
        self.assertNotEqual(control["staggered_max_absolute"], control["naive_max_absolute"])

    def test_probability_current_is_not_the_abel_wronskian(self) -> None:
        current = real_energy_probability_current(
            PRIMARY_RADIUS, 0.4, PRIMARY_KAPPA, "parent"
        )
        self.assertFalse(current["identified_with_abel_wronskian"])
        self.assertLess(current["relative_spread"], 1.0e-8)
        self.assertGreater(min(current["current_values"]), .99)

    def test_inverse_radius_scaling_uses_band_products(self) -> None:
        rows = [
            {"radius": 4.0, "first_positive": 0.705},
            {"radius": 8.0, "first_positive": 0.382},
            {"radius": 12.0, "first_positive": 0.259},
        ]
        scaling = inverse_radius_scaling(rows)
        self.assertTrue(scaling["scales_as_inverse_radius"])
        self.assertTrue(scaling["e1_decreasing"])
        self.assertGreater(scaling["first_positive_times_radius"][-1], 2.0)


class PublishedRecordTests(unittest.TestCase):
    def test_output_refuses_to_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nsc-3-boundary-response.json"
            path.write_bytes(b"{}")
            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(ROOT / "scripts/check_nsc_boundary_response.py"),
                    "--output",
                    str(path),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("refusing to overwrite", completed.stderr)
            self.assertEqual(path.read_bytes(), b"{}")

    def test_published_record_flags_if_present(self) -> None:
        if not DEFAULT_RESULT.is_file():
            self.skipTest("machine result has not been published yet")
        record = json.loads(DEFAULT_RESULT.read_bytes())
        self.assertEqual(record["schema"], SCHEMA)
        self.assertTrue(record["gate"]["parent_is_positive_rho"])
        self.assertTrue(record["gate"]["child_is_negative_rho"])
        self.assertTrue(record["gate"]["sparse_lu_is_primary"])
        self.assertTrue(record["gate"]["oriented_jump_identity_exact"])
        self.assertFalse(record["nonclaims"]["staggered_edge_trace_is_continuum_Calderon"])
        self.assertFalse(record["nonclaims"]["compact_zero_mode_is_exact_invariant_subspace"])
        self.assertFalse(record["gate"]["physical_zeta_selected"])
        self.assertFalse(record["nonclaims"]["stationary_scale_selected"])
        self.assertFalse(record["warped_spatial_coupling"]["compact_gap_theorem_from_this_finite_matrix"])
        self.assertGreater(record["half_domain_maps"][0]["parent"]["edge_rho"], 0.0)
        self.assertLess(record["half_domain_maps"][0]["child"]["edge_rho"], 0.0)


if __name__ == "__main__":
    unittest.main()
