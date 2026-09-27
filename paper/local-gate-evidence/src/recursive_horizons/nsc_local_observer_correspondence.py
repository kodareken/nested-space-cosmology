"""Compact local-observer correspondence of bound finite identities.

This owner independently checks small-matrix Schur/resolvent, retarded
memory with initial-state/noise data, common-action variation and
conservation, normalized inherited recursion, and the two-sheet
Dirac/charge-conjugation embedding. It binds existing records by hash.
It does not claim Schur/Feshbach novelty, LambdaCDM equivalence, or
identity of sheet exchange with antimatter or charge conjugation.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp
from scipy.linalg import expm

from .nsc_causal_common import (
    BoundaryDomain,
    CausalCommonFunctional,
    GaugeHistory,
    LinkHistory,
    MetricHistory,
    SpectralInducedSource,
)
from .nsc_energy_transfer import expectation, regional_operators, schur_state_data
from .nsc_influence import ground_covariance
from .nsc_regulated import direct_chain, recursive_response
from .nsc_spinor_bridge import exact_checks, pauli, weyl_matrices


SCHEMA = "NSC-LOCAL-OBSERVER-CORRESPONDENCE-v1"
ARTIFACT_ID = "NSC-LOCAL-OBSERVER-CORRESPONDENCE"
EXACT = 0
NUMERIC = 1e-12
NEGATIVE = 0.05
RESOLVENT = 0.4 + 0.5j

OWNERS = (
    "src/recursive_horizons/nsc_local_observer_correspondence.py",
    "scripts/derive_nsc_local_observer_correspondence.py",
    "tests/test_nsc_local_observer_correspondence.py",
    "docs/nsc-local-observer-correspondence.md",
    "src/recursive_horizons/evidence_io.py",
)
INPUTS = (
    "results/nsc-3-boundary-response.json",
    "results/nsc-3-regulated-recursion.json",
    "results/nsc-6-energy-transfer.json",
    "results/nsc-6-vacuum-work.json",
    "results/nsc-7-observable-bridge.json",
    "docs/nsc-common-source-derivation.md",
    "docs/nsc-declared-action-scope.md",
    "docs/nsc-energy-transfer.md",
    "docs/nsc-regulated-recursion.md",
    "docs/nsc-observable-bridge.md",
    "src/recursive_horizons/nsc_energy_transfer.py",
    "src/recursive_horizons/nsc_regulated.py",
    "src/recursive_horizons/nsc_spinor_bridge.py",
    "src/recursive_horizons/nsc_causal_common.py",
    "src/recursive_horizons/nsc_influence.py",
)

EQUATIONS = {
    "block_schur": "S = K_A - B K_B^{-1} B^\\dagger,  G_AA = S^{-1} = ((zI-H)^{-1})_{AA}",
    "feshbach_self_energy": "Sigma_A(z) = B (zI-H_BB)^{-1} B^\\dagger",
    "retarded_memory": (
        "i dpsi_A/dt = H_AA psi_A - i \\int_0^t B exp[-i H_BB(t-s)] B^\\dagger psi_A(s) ds "
        "+ B exp[-i H_BB t] psi_B(0)"
    ),
    "initial_noise": "B exp[-i H_BB t] C_BB(0) exp[i H_BB s] B^\\dagger, plus initial C_AB",
    "energy_balance": "dot E_A = Tr(C J_A) + Tr(C dot h_A),  J_A = i[H,h_A]",
    "total_conservation": "dot E_tot = Tr(C dot H);  Tr([H,C] H) = 0",
    "gaussian_metric_variation": (
        "delta Gamma_F / delta X^A = -Tr(C delta H / delta X^A) "
        "on the physical history at fixed preparation"
    ),
    "gaussian_schur_variation": (
        "delta log det K = Tr(G_c delta K_c) + Tr(S^{-1} delta S),  "
        "S = K_p - V G_c W"
    ),
    "normalized_recursion": (
        "Gamma_n(x) = x I - H - (1/Omega) b Gamma_{n+1}(x/Omega)^{-1} b^\\dagger,  "
        "H_n = Omega^n H, B_n = Omega^n b"
    ),
    "scalar_sheet_embedding": (
        "Pi_- = (I_8 - tau_3 \\otimes gamma^5)/2,  "
        "W^\\dagger H_8 W = alpha·p + Phi beta"
    ),
    "charge_conjugation": "psi^c = i gamma^2 psi^* with the conjugate gauge representation",
}

LAYERS = {
    "imported_algebra": [
        "Block Schur/Feshbach reduction of a finite resolvent is standard linear algebra.",
        "Retarded memory plus an eliminated initial state is standard nonequilibrium Green-function structure (Jauho–Wingreen–Meir; Martín–Verdaguer as methodological prior art).",
        "Klich's Fock determinant identity is imported; this record does not re-prove it.",
        "Dirac charge conjugation and Clifford algebra are imported; the scalar two-sheet bilinear is a kinematic benchmark, not a new Clifford theorem.",
    ],
    "repository_derived_identities": [
        "Noncommuting block Schur agrees with the joined resolvent restriction; reversed multiplication differs.",
        "The Laplace/resolvent identity retains B (zI-H_BB)^{-1} psi_B(0); omitting it changes the solution.",
        "The retarded self-energy is independent of the occupation, while the initial-noise term B C_BB B^\\dagger is not.",
        "Finite regional energy obeys dot E_A = Tr(C J_A)+Tr(C dot h_A) for h_A={P,H}/2.",
        "Finite Gaussian block elimination counts log det K_c once; all link variations sit in delta S.",
        "Normalized recursion with 1/Omega agrees with the direct multi-room inverse on the declared finite chain.",
        "A local scalar sheet coupling admits an invariant rank-four Dirac sector; the complementary sector is equally invariant.",
        "Combined sheet parity and (tau_1 ⊗ i gamma_2)K restrict to ordinary Dirac parity and C; sheet exchange alone does neither.",
    ],
    "nsc_interpretation": [
        "A local observer in one nested region is described by the reduced block, its retarded memory, the retained initial-state/noise data, and the same-action variation of the remaining Hamiltonian.",
        "Inheritance T_Theta^* D_Theta = D_Theta is the organizing invariance; the finite normalized chain is a control of that law, not an eternal cosmology.",
        "Inside/outside are relational roles of the region map, not gauge-charge assignments and not antimatter labels.",
        "Metric, link and relative-scale equations are different variations of one declared action; no independent dark-fluid function is added here.",
    ],
    "nonclaims": [
        "Schur/Feshbach algebra is not NSC novelty.",
        "No complete LambdaCDM background or perturbation match is established.",
        "Sheet exchange is not charge conjugation and is not antimatter.",
        "Antimatter is not mechanical pressure.",
        "The retarded map alone does not select an occupation or an energy current.",
        "The local incoming gate under evolved C_Sigma[g] remains OPEN.",
        "Prescribed geometry work is not same-action backreaction.",
    ],
}

CLAIM_STATUS = {
    "exact_block_schur_resolvent_reduction": "PASS",
    "retained_retarded_memory_initial_state_noise": "PASS",
    "common_action_metric_variation_conservation": "PASS_FINITE_OWNERSHIP",
    "normalized_inherited_recursion": "PASS",
    "two_sheet_dirac_charge_conjugation_embedding": "PASS",
    "lambdacdm_background_perturbation_match": "NONCLAIM",
    "sheet_exchange_equals_antimatter": "NONCLAIM",
    "sheet_exchange_equals_charge_conjugation": "NONCLAIM",
    "schur_feshbach_algebra_is_nsc_novelty": "NONCLAIM",
    "antimatter_is_mechanical_pressure": "NONCLAIM",
    "retarded_map_selects_occupation_or_current": "NONCLAIM",
    "local_incoming_gate": "OPEN",
    "same_action_cosmological_backreaction": "OPEN",
    "absolute_continuum_stress_match": "OPEN",
}

NONCLAIMS = {
    "complete_LambdaCDM_background_perturbation_match": False,
    "sheet_exchange_is_antimatter": False,
    "sheet_exchange_is_charge_conjugation": False,
    "schur_feshbach_algebra_is_nsc_novelty": False,
    "antimatter_is_mechanical_pressure": False,
    "retarded_map_selects_occupation_or_energy_current": False,
    "local_incoming_EXISTENCE": False,
    "local_incoming_NONEXISTENCE": False,
    "same_action_cosmological_backreaction_solved": False,
    "absolute_continuum_stress_matched": False,
    "physical_sector_selected_by_actual_throat": False,
}

MISSING_COMMON_ACTION = (
    "absolute_covariantly_renormalized_continuum_stress_on_the_varying_throat",
    "same_action_metric_backreaction_of_an_evolving_geometry",
    "full_continuum_CTP_coefficient_match",
    "complete_LambdaCDM_background_H_of_z_and_perturbation_growth_lensing_match",
    "closed_local_incoming_gate_under_evolved_C_Sigma",
)


def digest(root: Path, relative: str) -> str:
    return sha256((root / relative).read_bytes()).hexdigest()


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise AssertionError(message)


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.complexfloating, complex)):
        return [float(value.real), float(value.imag)]
    if value is None or isinstance(value, (str, int)):
        return value
    raise TypeError(f"unsupported record value {type(value)!r}")


def canonical_record(value: dict[str, Any]) -> dict[str, Any]:
    import json

    return json.loads(json.dumps(_jsonable(value), sort_keys=True, allow_nan=False))


def toy_block_hamiltonian() -> tuple[np.ndarray, np.ndarray]:
    """Independent 2+2 Hermitian control, not the published Dirac lattice."""
    aa = np.array([[1.2, 0.1 - 0.2j], [0.1 + 0.2j, 0.8]], complex)
    bb = np.array([[0.7, 0.3j], [-0.3j, 1.4]], complex)
    coupling = np.array([[0.25, 0.1j], [0.05, 0.2 - 0.1j]], complex)
    hamiltonian = np.block([[aa, coupling], [coupling.conj().T, bb]])
    _require(np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-15), "toy H must be Hermitian")
    mask = np.array([True, True, False, False])
    return hamiltonian, mask


def block_schur_resolvent_reduction() -> dict[str, Any]:
    """Exact finite Schur versus the joined resolvent, with block negative controls."""
    hp = sp.Matrix([[2, sp.I / 5], [-sp.I / 5, 3]])
    hc = sp.Matrix([[1, sp.Rational(1, 3)], [sp.Rational(1, 3), sp.Rational(5, 2)]])
    coupling = sp.Matrix([[sp.Rational(1, 4), sp.I / 8], [sp.Rational(2, 7), -sp.I / 6]])
    z = 3 + 2 * sp.I
    rp = z * sp.eye(2) - hp
    rc = z * sp.eye(2) - hc
    full = rp.row_join(-coupling).col_join((-coupling.H).row_join(rc))
    schur = rp - coupling * rc.inv() * coupling.H
    residual = (full.inv()[:2, :2] - schur.inv()).applyfunc(sp.simplify)
    reversed_order = rp - coupling.H * rc.inv() * coupling
    order_difference = (schur - reversed_order).applyfunc(sp.simplify)
    dropped_adjoint = rp - coupling * rc.inv() * coupling
    adjoint_difference = (schur - dropped_adjoint).applyfunc(sp.simplify)
    naive = rp - coupling * coupling.H / z
    naive_difference = (schur - naive).applyfunc(sp.simplify)
    _require(residual == sp.zeros(2), "exact Schur disagrees with the joined resolvent")
    _require(order_difference != sp.zeros(2), "reversed multiplication collapsed")
    _require(adjoint_difference != sp.zeros(2), "dropping B^\\dagger collapsed")
    _require(naive_difference != sp.zeros(2), "naive B B^\\dagger/z collapsed")

    hamiltonian, mask = toy_block_hamiltonian()
    aa, bb, ab = hamiltonian[np.ix_(mask, mask)], hamiltonian[np.ix_(~mask, ~mask)], hamiltonian[np.ix_(mask, ~mask)]
    energy = RESOLVENT
    sigma = ab @ np.linalg.solve(energy * np.eye(len(bb)) - bb, ab.T.conj())
    reduced = np.linalg.inv(energy * np.eye(len(aa)) - aa - sigma)
    direct = np.linalg.inv(energy * np.eye(len(hamiltonian)) - hamiltonian)[np.ix_(mask, mask)]
    numeric = float(np.max(np.abs(reduced - direct)))
    wrong_order = ab.T.conj() @ np.linalg.solve(energy * np.eye(len(bb)) - bb, ab)
    swapped = np.linalg.inv(energy * np.eye(len(aa)) - aa - wrong_order)
    swapped_error = float(np.max(np.abs(swapped - direct)))
    _require(numeric < NUMERIC, "numeric Schur residual exceeds tolerance")
    _require(swapped_error > NEGATIVE, "wrong-order block control was not negative")
    return {
        "equation": EQUATIONS["block_schur"],
        "layer": "imported_algebra_independently_verified",
        "exact_joined_versus_schur_residual": str(residual),
        "reversed_multiplication_differs": True,
        "reversed_multiplication_residual": str(order_difference),
        "dropped_adjoint_differs": True,
        "naive_bbdagger_over_z_differs": True,
        "numeric_joined_versus_schur_residual": numeric,
        "wrong_block_order_residual": swapped_error,
        "negative_controls": [
            "reversed_multiplication",
            "dropped_B_adjoint",
            "naive_B_Bdagger_over_z",
            "wrong_block_order_on_numeric_split",
        ],
        "novelty_claimed": False,
    }


def retarded_memory_initial_state_noise() -> dict[str, Any]:
    """Retarded self-energy plus the eliminated initial state and its noise."""
    hamiltonian, mask = toy_block_hamiltonian()
    aa, bb, ab = hamiltonian[np.ix_(mask, mask)], hamiltonian[np.ix_(~mask, ~mask)], hamiltonian[np.ix_(mask, ~mask)]
    z = RESOLVENT
    rng = np.random.default_rng(73)
    initial = rng.normal(size=len(hamiltonian)) + 1j * rng.normal(size=len(hamiltonian))
    child_resolvent = np.linalg.inv(z * np.eye(len(bb)) - bb)
    effective = z * np.eye(len(aa)) - aa - ab @ child_resolvent @ ab.T.conj()
    retained = np.linalg.solve(effective, initial[mask] + ab @ child_resolvent @ initial[~mask])
    direct = np.linalg.solve(z * np.eye(len(hamiltonian)) - hamiltonian, initial)[mask]
    omitted = np.linalg.solve(effective, initial[mask])
    retained_error = float(np.max(np.abs(retained - direct)))
    omitted_error = float(np.linalg.norm(omitted - direct))
    _require(retained_error < NUMERIC, "initial-state resolvent identity failed")
    _require(omitted_error > NEGATIVE, "omitted initial child source was not a negative control")

    energies, vectors = np.linalg.eigh(hamiltonian)
    _require(float(np.min(np.abs(energies))) > 1e-10, "toy spectrum must avoid zero modes")
    empty = np.zeros_like(hamiltonian)
    filled = np.eye(len(hamiltonian), dtype=complex)
    empty_data = schur_state_data(hamiltonian, mask, z, empty)
    filled_data = schur_state_data(hamiltonian, mask, z, filled)
    _require(empty_data["direct_reduced_error"] < NUMERIC, "empty-state Schur identity failed")
    _require(
        abs(empty_data["self_energy_norm"] - filled_data["self_energy_norm"]) < NUMERIC,
        "self-energy must be independent of occupation",
    )
    _require(
        abs(filled_data["initial_noise_norm"] - empty_data["initial_noise_norm"]) > NEGATIVE,
        "initial noise must change with occupation",
    )
    ops = regional_operators(hamiltonian, mask)
    vacuum_current = abs(expectation(empty, ops["current_a"]))
    packet = expm(-0.6j * hamiltonian) @ np.array([1.0, 0.0, 0.0, 0.0], complex)
    excited = np.outer(packet, packet.conj())
    excited_current = abs(expectation(excited, ops["current_a"]))
    _require(vacuum_current < 1e-12, "empty-state current must vanish")
    _require(excited_current > NEGATIVE, "a transported packet must source a current")
    return {
        "equation": EQUATIONS["retarded_memory"],
        "noise_equation": EQUATIONS["initial_noise"],
        "layer": "repository_derived_identities",
        "resolvent_with_initial_child_residual": retained_error,
        "omitted_initial_child_residual": omitted_error,
        "self_energy_independent_of_occupation": True,
        "vacuum_self_energy_norm": empty_data["self_energy_norm"],
        "excited_self_energy_norm": filled_data["self_energy_norm"],
        "vacuum_initial_noise_norm": empty_data["initial_noise_norm"],
        "excited_initial_noise_norm": filled_data["initial_noise_norm"],
        "vacuum_current": vacuum_current,
        "excited_current": excited_current,
        "negative_controls": [
            "omit_initial_child_source",
            "occupation_changes_noise_not_self_energy",
            "same_schur_map_zero_vacuum_current_nonzero_excited_current",
        ],
        "retarded_map_selects_occupation": False,
    }


def common_action_variation_conservation() -> dict[str, Any]:
    """Same-action finite variation and conservation; continuum match remains open."""
    hamiltonian_sym = sp.Matrix([[2, sp.I], [-sp.I, 3]])
    covariance_sym = sp.Matrix([[sp.Rational(1, 2), sp.I / 8], [-sp.I / 8, sp.Rational(1, 3)]])
    commutator_work = sp.simplify(((hamiltonian_sym * covariance_sym - covariance_sym * hamiltonian_sym) * hamiltonian_sym).trace())
    _require(commutator_work == 0, "Tr([H,C]H) failed as a cyclic identity")

    kp = sp.Matrix([[2, sp.I / 7], [-sp.I / 7, 3]])
    kc = sp.Matrix([[4, sp.Rational(1, 5)], [sp.Rational(1, 5), sp.Rational(5, 2)]])
    v = sp.Matrix([[sp.Rational(1, 6), sp.I / 9], [sp.Rational(1, 8), 0]])
    w = v.H
    dkp = sp.Matrix([[sp.Rational(1, 10), -sp.I / 11], [sp.I / 11, -sp.Rational(1, 12)]])
    dkc = sp.Matrix([[0, sp.Rational(1, 13)], [sp.Rational(1, 13), sp.Rational(1, 14)]])
    dv = sp.Matrix([[sp.I / 15, sp.Rational(1, 16)], [0, -sp.I / 17]])
    dw = dv.H
    kernel = kp.row_join(-v).col_join((-w).row_join(kc))
    dkernel = dkp.row_join(-dv).col_join((-dw).row_join(dkc))
    child_green = kc.inv()
    schur = kp - v * child_green * w
    dschur = dkp - dv * child_green * w - v * child_green * dw + v * child_green * dkc * child_green * w
    det_residual = sp.simplify(kernel.det() - kc.det() * schur.det())
    variation = sp.simplify(
        (kernel.inv() * dkernel).trace() - (child_green * dkc).trace() - (schur.inv() * dschur).trace()
    )
    omitted_child = sp.simplify((kernel.inv() * dkernel).trace() - (schur.inv() * dschur).trace())
    _require(det_residual == 0, "finite log-det Schur factorization failed")
    _require(variation == 0, "finite Gaussian Schur variation failed")
    _require(omitted_child != 0, "omitting the child variation was not a negative control")

    hamiltonian, mask = toy_block_hamiltonian()
    ops = regional_operators(hamiltonian, mask)
    energies, vectors = np.linalg.eigh(hamiltonian)
    vacuum = vectors[:, energies < 0] @ vectors[:, energies < 0].T.conj()
    _require(abs(expectation(vacuum, ops["current_a"])) < 1e-12, "vacuum regional current failed")
    two = np.array([[1.1, 0.4], [0.4, 1.1]], float)
    split = regional_operators(two, [1, 0])
    psi = np.array([1.0, 0.0])
    time = 0.37
    state = expm(-1j * time * two) @ psi
    energy_a = float(np.vdot(state, split["energy_a"] @ state).real)
    expected = float(two[0, 0] * np.cos(two[0, 1] * time) ** 2)
    _require(abs(energy_a - expected) < 1e-12, "two-level regional energy failed")

    times = np.array([0.0, 0.13, 0.31])
    parent = np.repeat(np.array([[[-0.8]]], complex), len(times), axis=0)
    child = np.repeat(np.array([[[0.9]]], complex), len(times), axis=0)
    links = np.repeat(np.array([[[0.23 + 0.04j]]], complex), len(times), axis=0)
    full = np.array([[-0.8, 0.23 + 0.04j], [0.23 - 0.04j, 0.9]])
    _, _, covariance = ground_covariance(full)
    vertices = {
        "N": np.array([[0.4, 0.1], [0.1, -0.2]]),
        "beta": np.array([[0.0, 0.2j], [-0.2j, 0.0]]),
        "q": np.diag([0.3, -0.1]),
        "r": np.array([[0.0, 0.07], [0.07, 0.0]]),
    }
    induced = SpectralInducedSource(
        forces={"N": 0.01, "beta": -0.02, "q": 0.03, "r": -0.04},
        coefficient_owner="finite common spectral control",
        unresolved_terms=("full continuum coefficient match",),
    )
    domain = BoundaryDomain(
        name="two-mode common-action ownership control",
        parent_dimension=1,
        child_dimension=1,
        spin_frame="common canonical frame",
        units="hbar=1",
        stress_factors={
            "rho": ("N", 2.0),
            "T_01": ("beta", 3.0),
            "p_parallel": ("q", -4.0),
            "p_perp": ("r", -5.0),
        },
        power_weights={"N": -0.5, "beta": 0.25, "q": 0.75},
    )
    owned = CausalCommonFunctional().evaluate(
        MetricHistory(times, parent, child, vertices),
        GaugeHistory(np.zeros_like(parent), np.zeros_like(child)),
        LinkHistory(links),
        covariance,
        induced,
        domain,
    )
    ward = {name: float(abs(value)) for name, value in owned.ward_residuals.items()}
    _require(ward["equal_history_action"] < 2e-15, "equal-history CTP action failed")
    _require(ward["energy_work_balance"] < 2e-14, "finite energy/work identity failed")
    force_split = max(
        abs(owned.forces[name] - owned.matter_forces[name] - owned.induced_forces[name])
        for name in owned.forces
    )
    _require(force_split < 1e-15, "induced source was not counted once")
    _require(owned.unresolved_terms == ("full continuum coefficient match",), "missing continuum term was dropped")
    return {
        "equations": {
            "energy_balance": EQUATIONS["energy_balance"],
            "total_conservation": EQUATIONS["total_conservation"],
            "gaussian_metric_variation": EQUATIONS["gaussian_metric_variation"],
            "gaussian_schur_variation": EQUATIONS["gaussian_schur_variation"],
        },
        "layer": "repository_derived_identities_with_imported_determinant_algebra",
        "cyclic_trace_identity": "Tr([H,C]H)=0",
        "finite_logdet_schur_residual": str(det_residual),
        "finite_schur_variation_residual": str(variation),
        "omitted_child_variation_residual": str(omitted_child),
        "two_level_regional_energy_residual": abs(energy_a - expected),
        "owned_ctp_max_abs_ward_residual": max(ward.values()),
        "owned_ctp_unresolved_terms": list(owned.unresolved_terms),
        "induced_source_counted_once": True,
        "negative_controls": [
            "omit_child_block_variation_from_logdet",
            "drop_induced_or_matter_half_of_the_common_force",
        ],
        "missing_common_action_evidence": list(MISSING_COMMON_ACTION),
        "stronger_result_blocked_by_missing_common_action_evidence": True,
    }


def normalized_inherited_recursion() -> dict[str, Any]:
    """Parent-normalized 1/Omega chain versus the direct multi-room inverse."""
    local = np.array([[1.1, 0.2j], [-0.2j, 0.9]], complex)
    link = np.array([[0.15, 0.05j], [0.08, 0.22]], complex)
    energy = 0.5 + 0.3j
    omega = 1.6
    depth = 5
    reduced = recursive_response(local, link, energy, omega, depth)
    direct = direct_chain(local, link, energy, omega, depth)
    error = float(np.linalg.norm(reduced - direct))
    _require(error < NUMERIC, "normalized recursion disagrees with the direct chain")
    unweighted = recursive_response(local, link, energy, 1.0, depth)
    unweighted_error = float(np.linalg.norm(unweighted - direct))
    _require(unweighted_error > NEGATIVE, "dropping 1/Omega was not a negative control")
    doubled_error = float(
        np.linalg.norm(recursive_response(local, 2 * link, energy, omega, depth) - direct)
    )
    _require(doubled_error > NEGATIVE, "rescaled link was not a negative control")
    unit_error = float(
        np.linalg.norm(
            recursive_response(local, link, energy, 1.0, depth)
            - direct_chain(local, link, energy, 1.0, depth)
        )
    )
    _require(unit_error < NUMERIC, "Omega=1 recursion failed")
    imag = (reduced - reduced.conj().T) / (2j)
    imag_min = float(np.linalg.eigvalsh(imag).min())
    _require(imag_min > 0, "reduced response lost positive imaginary part")
    return {
        "equation": EQUATIONS["normalized_recursion"],
        "layer": "repository_derived_identities",
        "omega": omega,
        "depth": depth,
        "direct_reduced_error": error,
        "unweighted_versus_normalized_error": unweighted_error,
        "doubled_link_error": doubled_error,
        "omega_one_direct_reduced_error": unit_error,
        "minimum_imaginary_eigenvalue": imag_min,
        "negative_controls": ["drop_one_over_omega_weight", "double_link_amplitude"],
        "selects_physical_Omega": False,
        "eternal_time_dependent_solution": False,
    }


def two_sheet_dirac_charge_conjugation() -> dict[str, Any]:
    """Scalar-sheet Dirac embedding, with sheet exchange distinct from C."""
    algebra = exact_checks()
    _require(set(algebra["identities"].values()) == {"0"}, "sheet/Dirac identities failed")
    _require(not algebra["candidate_restriction"]["selected_by_actual_throat_or_action"], "sector selection was overclaimed")
    w = weyl_matrices()
    charge = np.array(w["charge_conjugation_matrix"], complex)
    tau1 = np.array(pauli()[0], complex)
    swap = np.kron(tau1, np.eye(4, dtype=complex))
    paired_c = np.kron(tau1, charge)
    psi = np.array([1, 0.3j, -0.2, 0.5j, 0.1, -0.4j, 0.7, 0.2j], complex)
    exchange_versus_c = float(np.linalg.norm(swap @ psi - paired_c @ psi.conj()))
    exchange_versus_linear_c = float(np.linalg.norm(swap @ psi - paired_c @ psi))
    _require(exchange_versus_c > NEGATIVE, "sheet exchange coincided with C")
    _require(exchange_versus_linear_c > NEGATIVE, "sheet exchange coincided with the C matrix")

    kinetic = sum(p * np.array(matrix, complex) for p, matrix in zip((0.3, -0.7, 1.1), w["alpha"]))
    beta = np.array(w["gamma"][0], complex)
    phi = 0.6
    full = np.kron(np.eye(2), kinetic) + phi * np.kron(tau1, beta)
    embedding = np.eye(8)[:, [0, 1, 6, 7]]
    initial = np.array([1, 2j, -0.3, 0.7j], complex)
    initial /= np.linalg.norm(initial)
    evolved_full = expm(-0.8j * full) @ embedding @ initial
    evolved_reduced = embedding @ expm(-0.8j * (kinetic + phi * beta)) @ initial
    embedding_error = float(np.max(np.abs(evolved_full - evolved_reduced)))
    _require(embedding_error < NUMERIC, "restricted-sector evolution failed")
    identity_link = np.kron(np.eye(2), np.array(w["alpha"][2], complex)) + np.kron(tau1, np.eye(4))
    identity_gap = float(np.min(np.abs(np.linalg.eigvalsh(identity_link))))
    _require(identity_gap < 1e-12, "identity sheet link must be gapless at |p|=Phi=1")
    return {
        "equations": {
            "scalar_sheet_embedding": EQUATIONS["scalar_sheet_embedding"],
            "charge_conjugation": EQUATIONS["charge_conjugation"],
        },
        "layer": "repository_derived_identities_on_imported_clifford_algebra",
        "all_exact_identities_zero": True,
        "joint_projector_ranks": algebra["joint_projector_ranks"],
        "restricted_sector": algebra["candidate_restriction"],
        "sheet_exchange_operator": "tau_1 ⊗ I_4",
        "charge_conjugation_operator": "(tau_1 ⊗ i gamma_2) K",
        "sheet_exchange_versus_C_state_residual": exchange_versus_c,
        "sheet_exchange_versus_C_matrix_residual": exchange_versus_linear_c,
        "restricted_sector_evolution_residual": embedding_error,
        "identity_link_gap_at_unit_momentum": identity_gap,
        "negative_controls": [
            "sheet_exchange_is_not_charge_conjugation",
            "identity_sheet_link_is_not_a_mass_gap",
            "complementary_rank_four_sector_is_equally_invariant",
        ],
        "sheet_exchange_is_antimatter": False,
        "sheet_exchange_is_charge_conjugation": False,
        "selected_by_actual_throat_or_action": False,
    }


def correspondence_record(root: Path | None = None) -> dict[str, Any]:
    """Assemble the compact v1 correspondence record."""
    root = Path(root) if root is not None else Path(__file__).resolve().parents[2]
    schur = block_schur_resolvent_reduction()
    memory = retarded_memory_initial_state_noise()
    action = common_action_variation_conservation()
    recursion = normalized_inherited_recursion()
    embedding = two_sheet_dirac_charge_conjugation()
    record = {
        "schema": SCHEMA,
        "artifact_id": ARTIFACT_ID,
        "accountable_author": "Douglas Ek",
        "classification": (
            "compact_local_observer_correspondence_of_bound_finite_identities_"
            "not_lambdacdm_or_antimatter_closure"
        ),
        "status": (
            "PASS: compact local-observer correspondence of bound finite identities; "
            "no LambdaCDM match and no sheet-exchange/antimatter identity; "
            "local incoming gate OPEN"
        ),
        "physical_claim_verified": False,
        "claim_status": dict(CLAIM_STATUS),
        "layers": {name: list(items) for name, items in LAYERS.items()},
        "equations": dict(EQUATIONS),
        "independent_controls": {
            "block_schur_resolvent_reduction": schur,
            "retarded_memory_initial_state_noise": memory,
            "common_action_metric_variation_conservation": action,
            "normalized_inherited_recursion": recursion,
            "two_sheet_dirac_charge_conjugation_embedding": embedding,
        },
        "bound_source_records": {
            "block_schur": "results/nsc-3-boundary-response.json",
            "retarded_memory_and_conservation": "results/nsc-6-energy-transfer.json",
            "prescribed_geometry_work": "results/nsc-6-vacuum-work.json",
            "normalized_recursion": "results/nsc-3-regulated-recursion.json",
            "two_sheet_embedding": "results/nsc-7-observable-bridge.json",
            "common_action_derivation": "docs/nsc-common-source-derivation.md",
            "declared_action_scope": "docs/nsc-declared-action-scope.md",
        },
        "missing_common_action_evidence": list(MISSING_COMMON_ACTION),
        "nonclaims": dict(NONCLAIMS),
        "scope": {
            "local_observer_correspondence": True,
            "finite_matrix_identities": True,
            "imported_algebra_not_claimed_as_novelty": True,
            "nsc_interpretation_separated_from_identities": True,
            "lambdacdm_match": False,
            "sheet_exchange_antimatter_identity": False,
            "scientific_campaign_or_expensive_generator": False,
            "local_incoming_gate": "OPEN",
            "arxiv_publication": False,
            "metric_timestep": False,
            "global_parent_child_matching": False,
        },
        "source_hashes": {path: digest(root, path) for path in OWNERS},
        "input_hashes": {path: digest(root, path) for path in INPUTS},
        "reproducer": "python3 scripts/derive_nsc_local_observer_correspondence.py --check",
    }
    return canonical_record(record)
