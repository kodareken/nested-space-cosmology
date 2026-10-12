"""Focused controls for the compact local-observer correspondence."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from recursive_horizons.nsc_energy_transfer import regional_operators, schur_state_data
from recursive_horizons.nsc_local_observer_correspondence import (
    CLAIM_STATUS,
    EQUATIONS,
    LAYERS,
    MISSING_COMMON_ACTION,
    NONCLAIMS,
    NUMERIC,
    RESOLVENT,
    block_schur_resolvent_reduction,
    common_action_variation_conservation,
    correspondence_record,
    normalized_inherited_recursion,
    retarded_memory_initial_state_noise,
    toy_block_hamiltonian,
    two_sheet_dirac_charge_conjugation,
)
from recursive_horizons.nsc_regulated import direct_chain, recursive_response
from recursive_horizons.nsc_spinor_bridge import exact_checks, pauli, weyl_matrices
import derive_nsc_local_observer_correspondence as OWNER


class TestBlockSchurNegativeControls(unittest.TestCase):
    def test_exact_schur_matches_joined_resolvent_and_wrong_products_fail(self):
        result = block_schur_resolvent_reduction()
        assert result["exact_joined_versus_schur_residual"] == "Matrix([[0, 0], [0, 0]])"
        assert result["reversed_multiplication_differs"]
        assert result["dropped_adjoint_differs"]
        assert result["naive_bbdagger_over_z_differs"]
        assert result["numeric_joined_versus_schur_residual"] < NUMERIC
        assert result["wrong_block_order_residual"] > 0.05
        assert result["novelty_claimed"] is False

    def test_independent_numpy_split_rejects_transposed_coupling(self):
        hamiltonian, mask = toy_block_hamiltonian()
        aa = hamiltonian[np.ix_(mask, mask)]
        bb = hamiltonian[np.ix_(~mask, ~mask)]
        ab = hamiltonian[np.ix_(mask, ~mask)]
        z = RESOLVENT
        sigma = ab @ np.linalg.solve(z * np.eye(len(bb)) - bb, ab.T.conj())
        reduced = np.linalg.inv(z * np.eye(len(aa)) - aa - sigma)
        direct = np.linalg.inv(z * np.eye(len(hamiltonian)) - hamiltonian)[np.ix_(mask, mask)]
        np.testing.assert_allclose(reduced, direct, atol=1e-12)
        transposed = ab.T @ np.linalg.solve(z * np.eye(len(bb)) - bb, ab.conj())
        wrong = np.linalg.inv(z * np.eye(len(aa)) - aa - transposed)
        assert np.max(np.abs(wrong - direct)) > 0.05


class TestRetardedMemoryAndNoise(unittest.TestCase):
    def test_initial_child_source_is_required_in_the_resolvent(self):
        result = retarded_memory_initial_state_noise()
        assert result["resolvent_with_initial_child_residual"] < NUMERIC
        assert result["omitted_initial_child_residual"] > 0.05
        assert result["self_energy_independent_of_occupation"]
        assert result["retarded_map_selects_occupation"] is False
        assert abs(result["excited_self_energy_norm"] - result["vacuum_self_energy_norm"]) < NUMERIC
        assert abs(result["excited_initial_noise_norm"] - result["vacuum_initial_noise_norm"]) > 0.05

    def test_owned_schur_state_data_keeps_noise_off_the_retarded_map(self):
        hamiltonian, mask = toy_block_hamiltonian()
        empty = np.zeros_like(hamiltonian)
        filled = np.eye(len(hamiltonian), dtype=complex)
        empty_data = schur_state_data(hamiltonian, mask, RESOLVENT, empty)
        filled_data = schur_state_data(hamiltonian, mask, RESOLVENT, filled)
        assert empty_data["direct_reduced_error"] < NUMERIC
        assert empty_data["self_energy_norm"] == filled_data["self_energy_norm"]
        assert filled_data["initial_noise_norm"] > empty_data["initial_noise_norm"] + 0.05
        ops = regional_operators(hamiltonian, mask)
        packet = expm(-0.6j * hamiltonian) @ np.array([1.0, 0.0, 0.0, 0.0], complex)
        assert abs(np.trace(empty @ ops["current_a"]).real) < 1e-12
        assert abs(np.vdot(packet, ops["current_a"] @ packet).real) > 0.05


class TestCommonActionOwnership(unittest.TestCase):
    def test_finite_schur_variation_counts_the_child_block_once(self):
        result = common_action_variation_conservation()
        assert result["finite_logdet_schur_residual"] == "0"
        assert result["finite_schur_variation_residual"] == "0"
        assert result["omitted_child_variation_residual"] != "0"
        assert result["owned_ctp_max_abs_ward_residual"] < 2e-14
        assert result["owned_ctp_unresolved_terms"] == ["full continuum coefficient match"]
        assert result["induced_source_counted_once"]
        assert result["stronger_result_blocked_by_missing_common_action_evidence"]
        assert "complete_LambdaCDM_background_H_of_z_and_perturbation_growth_lensing_match" in result[
            "missing_common_action_evidence"
        ]

    def test_two_level_work_identity_needs_the_explicit_operator_derivative(self):
        h0 = np.array([[1.0, 0.2], [0.2, 1.3]])
        vertex = np.diag([0.4, -0.1])
        ops = regional_operators(h0, [1, 0])
        work = regional_operators(vertex, [1, 0])["energy_a"]
        psi = np.array([1.0, 0.0])
        dt = 1e-6
        values = []
        for time in (-dt, dt):
            state = expm(-1j * h0 * time) @ psi
            operator = ops["energy_a"] + time * work
            values.append(np.vdot(state, operator @ state).real)
        exact = np.vdot(psi, (ops["current_a"] + work) @ psi).real
        assert abs((values[1] - values[0]) / (2 * dt) - exact) < 1e-8


class TestNormalizedRecursion(unittest.TestCase):
    def test_one_over_omega_chain_matches_direct_inverse(self):
        result = normalized_inherited_recursion()
        assert result["direct_reduced_error"] < NUMERIC
        assert result["unweighted_versus_normalized_error"] > 0.05
        assert result["doubled_link_error"] > 0.05
        assert result["selects_physical_Omega"] is False

    def test_independent_depth_three_matrix_rejects_missing_omega_weight(self):
        local = np.array([[0.6, 0.1j], [-0.1j, 1.3]], complex)
        link = np.array([[0.2, 0.04], [0.07j, 0.11]], complex)
        energy = 0.8 + 0.25j
        omega = 1.4
        np.testing.assert_allclose(
            recursive_response(local, link, energy, omega, 3),
            direct_chain(local, link, energy, omega, 3),
            atol=1e-12,
        )
        unweighted = recursive_response(local, link, energy, 1.0, 3)
        assert np.linalg.norm(unweighted - direct_chain(local, link, energy, omega, 3)) > 0.01


class TestTwoSheetEmbedding(unittest.TestCase):
    def test_sheet_exchange_is_not_charge_conjugation(self):
        result = two_sheet_dirac_charge_conjugation()
        assert result["all_exact_identities_zero"]
        assert result["sheet_exchange_versus_C_state_residual"] > 0.05
        assert result["sheet_exchange_is_antimatter"] is False
        assert result["sheet_exchange_is_charge_conjugation"] is False
        assert result["selected_by_actual_throat_or_action"] is False
        assert result["identity_link_gap_at_unit_momentum"] < 1e-12
        algebra = exact_checks()
        assert algebra["candidate_restriction"]["selected_by_actual_throat_or_action"] is False

    def test_independent_spinor_control_distinguishes_swap_from_C(self):
        w = weyl_matrices()
        charge = np.array(w["charge_conjugation_matrix"], complex)
        swap = np.kron(np.array(pauli()[0], complex), np.eye(4))
        paired = np.kron(np.array(pauli()[0], complex), charge)
        psi = np.ones(8, complex)
        psi[1] = 1j
        assert np.linalg.norm(swap @ psi - paired @ psi.conj()) > 0.5
        assert np.linalg.norm(swap - paired) > 0.5


class TestRecordLayersAndReplay(unittest.TestCase):
    def test_layers_and_nonclaims_are_explicit(self):
        record = correspondence_record(ROOT)
        assert record["schema"] == "NSC-LOCAL-OBSERVER-CORRESPONDENCE-v1"
        assert record["physical_claim_verified"] is False
        assert record["claim_status"] == CLAIM_STATUS
        assert record["nonclaims"] == NONCLAIMS
        assert record["claim_status"]["lambdacdm_background_perturbation_match"] == "NONCLAIM"
        assert record["claim_status"]["sheet_exchange_equals_antimatter"] == "NONCLAIM"
        assert record["nonclaims"]["complete_LambdaCDM_background_perturbation_match"] is False
        assert record["nonclaims"]["sheet_exchange_is_antimatter"] is False
        assert set(record["layers"]) == set(LAYERS)
        assert set(EQUATIONS) <= set(record["equations"])
        assert list(MISSING_COMMON_ACTION) == record["missing_common_action_evidence"]
        assert record["scope"]["arxiv_publication"] is False
        assert record["scope"]["scientific_campaign_or_expensive_generator"] is False

    def test_recorded_artifact_replays_when_present(self):
        if not OWNER.OUTPUT.is_file():
            return
        from hashlib import sha256

        live = OWNER.compute()
        saved = json.loads(OWNER.OUTPUT.read_text())
        assert live == saved
        for collection in ("source_hashes", "input_hashes"):
            for path, expected in saved[collection].items():
                assert sha256((ROOT / path).read_bytes()).hexdigest() == expected


if __name__ == "__main__":
    unittest.main()
