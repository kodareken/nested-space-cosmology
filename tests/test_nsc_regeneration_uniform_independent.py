"""Uniform six occupations, the projected Galerkin generator, and the G bracket.

A rank-6 preparation Phi, with Phi† Phi = I_6, defines

    C(c) = Phi diag(c) Phi†.

Equal weights c = (1/2, ..., 1/2) are C = (1/2) P_6. The physical column
evolution is the Galerkin generator on the fermion band,

    H_G = U_f† H_f U_f,    dimension 2 n_f.

C = (1/2) I_{2 n_f} commutes with H_G. The embedded quadrature matrix
C_f = (1/2) U_f U_f† does not commute with the unprojected H_f, because
H_f reaches modes the prolongation discards. That quadrature commutator
is not a flow of the finite Galerkin half-identity. C = (1/2) I_{2 n_q}
remains an algebraic fine-grid source control. It is not a demand to
occupy the quadrature complement as the physical model.

The owned one-block image is

    H = sigma_2 ⊗ {L/Q, P}/2 + sigma_1 ⊗ diag(kappa L) - I ⊗ {beta, P}/2,

and the mean-field source is M times its nodal derivatives, M = 4 kappa once.
On an even antiperiodic grid the half-integer momenta sum to zero and the
Fourier modes have flat modulus, so Tr({beta, P}) = 0 for every real shift.
Both Pauli matrices are traceless, so Tr H_f(L, Q, beta, kappa) = 0 on the
positive chart with or without the angular term. Therefore every nodal
derivative of M Tr((1/2) I_{2 n_q} H_f) vanishes. That is the fine-grid
algebra. The same column formula on (1/2) P_6 keeps the kinetic density K
of those six columns.

This is a raw finite Gamma_one statement. Gamma_one is the induced chart
counted once plus the canonical Gaussian of one finite covariance. It does
not add a sea, a heat kernel, or a cosmological counterterm. A zero raw
source is not a renormalized stress.

nf = 16 is inside the requested cap but the owned odd lobe is sampled only
at its midpoint and prepare_rank6 rejects it. nf = 14 is the finest even
count at or below 16 on which that packet resolves. Quadrature is 64.
A passing test is this control identity, not a regeneration result.
"""
from __future__ import annotations

import json
from dataclasses import replace
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from recursive_horizons import nsc_regeneration_controls as controls
from recursive_horizons import nsc_spherical_cauchy_weak as weak
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_conformal_adm_source import (
    conformal_source,
    direct_hamiltonian,
    sector_multiplicity,
    spin_identity_covariance,
)
from recursive_horizons.nsc_covariant_operator import CovariantStaticMetric
from recursive_horizons.nsc_influence import _covariance, influence
from recursive_horizons.nsc_spherical_coupling import (
    apply_dirac,
    chart_failure,
    hamilton_constraint,
    nodal_forces,
    source_from_columns,
)
from recursive_horizons.nsc_spherical_feedback_action import (
    feedback_F,
    feedback_V,
    feedback_Z,
    partial_F,
)

RECORD = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-regeneration-controls-v1.json"

NF = 14
NQ = 64
DT = 1.0e-3
MOMENTUM_SHIFT = 1.0e-3
ROUND = 1e-9


def _frobenius(matrix):
    return float(np.linalg.norm(np.asarray(matrix), ord="fro"))


def _covariance_from_columns(phi0, phi1, occupations):
    vectors = np.vstack((phi0, phi1))
    weights = np.asarray(occupations, dtype=float)
    return (vectors * weights) @ vectors.conj().T


def _metric(system, state):
    return CovariantStaticMetric(
        system.length,
        system.length_density * state.r,
        state.Q * state.r,
        state.r,
        eta=0.5,
    )


def _commutator(hamiltonian, covariance):
    return hamiltonian @ covariance - covariance @ hamiltonian


def _forces(source):
    return source["force_L"], source["force_Q"], source["force_beta"]


def _max_force(source):
    return float(max(np.max(np.abs(force)) for force in _forces(source)))


@lru_cache(maxsize=1)
def _measure():
    grid = galerkin.build_grid(NF, quadrature=NQ)
    phi0, phi1, _preparation, problems = galerkin.load_physical_columns(NF)
    coarse = galerkin.blank_state(grid, phi0, phi1)
    fine = galerkin.prolong_state(grid, coarse)
    system = grid.fine
    points = system.points
    metric = _metric(system, fine)
    hamiltonian = direct_hamiltonian(metric, system.shift, system.kappa)
    massless = direct_hamiltonian(metric, system.shift, 0)
    image0, image1 = apply_dirac(
        fine.phi0, fine.phi1, system.length_density, fine.Q, system.shift,
        system.kappa, system.momentum,
    )
    columns = np.vstack((fine.phi0, fine.phi1))
    uniform = np.full(6, 0.5)
    recorded = np.asarray(coupling.OCCUPATIONS, dtype=float)
    covariance_packet = _covariance_from_columns(fine.phi0, fine.phi1, uniform)
    covariance_recorded = _covariance_from_columns(fine.phi0, fine.phi1, recorded)
    covariance_identity = 0.5 * np.eye(2 * points, dtype=complex)
    eye = np.eye(points, dtype=complex)
    zero = np.zeros((points, points), dtype=complex)
    phi0_identity = np.concatenate((eye, zero), axis=1)
    phi1_identity = np.concatenate((zero, eye), axis=1)
    occupations_identity = np.full(2 * points, 0.5)
    band = grid.U_f
    phi0_band = np.concatenate((band, np.zeros_like(band)), axis=1)
    phi1_band = np.concatenate((np.zeros_like(band), band), axis=1)
    occupations_band = np.full(phi0_band.shape[1], 0.5)
    covariance_band = _covariance_from_columns(phi0_band, phi1_band, occupations_band)
    system_identity = replace(system, occupations=occupations_identity)
    system_uniform = replace(system, occupations=uniform)
    system_band = replace(system, occupations=occupations_band)
    state_identity = replace(fine, phi0=phi0_identity, phi1=phi1_identity)
    state_band = replace(fine, phi0=phi0_band, phi1=phi1_band)
    source_identity = source_from_columns(system_identity, state_identity)
    source_uniform = source_from_columns(system_uniform, fine)
    source_recorded = source_from_columns(system, fine)
    source_band = source_from_columns(system_band, state_band)
    bumped_q = fine.Q * (1.0 + 0.01 * np.cos(2.0 * np.pi * system.xi / system.length))
    bumped_state = replace(fine, Q=bumped_q)
    bumped_metric = _metric(system, bumped_state)
    bumped_hamiltonian = direct_hamiltonian(bumped_metric, system.shift, system.kappa)
    source_uniform_bumped = source_from_columns(system_uniform, bumped_state)
    source_band_bumped = source_from_columns(
        system_band, replace(state_band, Q=bumped_q)
    )
    source_identity_bumped = source_from_columns(
        system_identity, replace(state_identity, Q=bumped_q)
    )
    dense_identity = conformal_source(
        metric, system.shift, spin_identity_covariance(points), kappa=system.kappa,
    )
    dense_uniform = conformal_source(
        metric, system.shift, covariance_packet, kappa=system.kappa,
    )
    evolution = expm(-1j * hamiltonian * DT)
    shifted_momentum = system.momentum + MOMENTUM_SHIFT * np.eye(points)
    _kinetic, _mass, shifted_momentum_density = coupling.column_moments(
        phi0_identity, phi1_identity, occupations_identity, shifted_momentum,
    )
    shifted_forces = nodal_forces(
        _kinetic, _mass, shifted_momentum_density, system.length_density, fine.Q,
        system.kappa, system.multiplicity,
    )
    projector = columns @ columns.conj().T
    outside = (np.eye(2 * points) - projector) @ (hamiltonian @ columns)
    image = hamiltonian @ columns
    return {
        "grid": grid,
        "fine": fine,
        "system": system,
        "problems": problems,
        "hamiltonian": hamiltonian,
        "massless": massless,
        "bumped_hamiltonian": bumped_hamiltonian,
        "image_gap": float(np.max(np.abs(np.vstack((image0, image1)) - image))),
        "covariance_packet": covariance_packet,
        "covariance_recorded": covariance_recorded,
        "covariance_identity": covariance_identity,
        "covariance_band": covariance_band,
        "source_identity": source_identity,
        "source_uniform": source_uniform,
        "source_recorded": source_recorded,
        "source_band": source_band,
        "source_uniform_bumped": source_uniform_bumped,
        "source_band_bumped": source_band_bumped,
        "source_identity_bumped": source_identity_bumped,
        "dense_identity": dense_identity,
        "dense_uniform": dense_uniform,
        "evolution": evolution,
        "shifted_forces": shifted_forces,
        "outside": _frobenius(outside),
        "image": _frobenius(image),
        "points": points,
        "uniform": uniform,
        "recorded": recorded,
    }


def test_nf16_odd_lobe_is_unresolved_and_nf14_packet_is_orthonormal():
    with pytest.raises(ValueError, match="unresolved"):
        galerkin.load_physical_columns(16)
    measured = _measure()
    gram = measured["fine"].phi0.conj().T @ measured["fine"].phi0
    gram = gram + measured["fine"].phi1.conj().T @ measured["fine"].phi1
    assert measured["problems"] == []
    assert measured["grid"].nf == NF and measured["grid"].nq == NQ
    assert measured["grid"].nf <= 16 and measured["grid"].nq <= 64
    assert float(np.max(np.abs(gram - np.eye(6)))) < ROUND
    assert measured["image_gap"] < ROUND


def test_equal_six_occupations_differ_from_the_full_half_identity():
    measured = _measure()
    packet = measured["covariance_packet"]
    identity = measured["covariance_identity"]
    spatial_packet = np.real(np.diag(packet[: measured["points"], : measured["points"]]))
    spatial_packet = spatial_packet + np.real(
        np.diag(packet[measured["points"] :, measured["points"] :])
    )
    spatial_identity = np.real(np.diag(identity[: measured["points"], : measured["points"]]))
    spatial_identity = spatial_identity + np.real(
        np.diag(identity[measured["points"] :, measured["points"] :])
    )
    assert _frobenius(packet - identity) > 1.0
    assert abs(float(np.trace(packet).real) - 3.0) < ROUND
    assert abs(float(np.trace(identity).real) - float(measured["points"])) < ROUND
    assert float(np.max(spatial_packet) - np.min(spatial_packet)) > 1.0e-2
    assert float(np.max(spatial_identity) - np.min(spatial_identity)) < ROUND
    assert np.array_equal(measured["uniform"], np.full(6, 0.5))
    assert not np.array_equal(measured["recorded"], measured["uniform"])


def test_balanced_ap_trace_and_half_identity_forces_vanish():
    measured = _measure()
    system = measured["system"]
    hamiltonian = measured["hamiltonian"]
    identity = measured["covariance_identity"]
    multiplicity = int(system.multiplicity)
    assert multiplicity == 4 * int(system.kappa) == int(sector_multiplicity(system.kappa))
    metric_momenta = _metric(system, measured["fine"]).momenta
    assert abs(float(np.sum(metric_momenta))) < ROUND
    assert float(np.max(np.abs(metric_momenta + metric_momenta[::-1]))) < ROUND
    assert abs(float(np.sum(measured["grid"].modes_f))) < ROUND
    assert float(np.max(np.abs(np.diag(system.momentum)))) < ROUND
    assert float(np.max(np.abs(hamiltonian - hamiltonian.conj().T))) < ROUND
    assert abs(complex(np.trace(hamiltonian))) < ROUND
    assert abs(complex(np.trace(measured["massless"]))) < ROUND
    assert abs(complex(np.trace(measured["bumped_hamiltonian"]))) < ROUND
    assert _frobenius(_commutator(hamiltonian, identity)) < ROUND
    step = measured["evolution"] @ identity @ measured["evolution"].conj().T - identity
    assert _frobenius(step) < ROUND
    assert _max_force(measured["source_identity"]) < ROUND
    assert _max_force(measured["source_identity_bumped"]) < ROUND
    assert measured["source_identity"]["multiplicity_applied_once"] is True
    for name, column_name in (("L", "force_L"), ("Q", "force_Q"), ("beta", "force_beta")):
        gap = measured["dense_identity"]["nodal"][name] - measured["source_identity"][column_name]
        assert float(np.max(np.abs(gap))) < ROUND
    assert int(measured["dense_identity"]["multiplicity"]) == multiplicity
    assert measured["dense_identity"]["vacuum_subtracted"] is False
    assert abs(float(measured["dense_identity"]["energy"])) < ROUND
    shifted = measured["shifted_forces"]
    expected_current = -multiplicity * MOMENTUM_SHIFT
    assert float(np.max(np.abs(shifted[0]))) < ROUND
    assert float(np.max(np.abs(shifted[1]))) < ROUND
    assert abs(float(np.mean(shifted[2])) - expected_current) < ROUND
    assert abs(float(np.max(np.abs(shifted[2]))) - abs(expected_current)) < ROUND


def test_rank6_uniform_still_flows_and_imbalance_remains_a_separate_increment():
    measured = _measure()
    hamiltonian = measured["hamiltonian"]
    packet = measured["covariance_packet"]
    recorded = measured["covariance_recorded"]
    packet_commutator = _frobenius(_commutator(hamiltonian, packet))
    recorded_commutator = _frobenius(_commutator(hamiltonian, recorded))
    packet_step = measured["evolution"] @ packet @ measured["evolution"].conj().T - packet
    recorded_step = measured["evolution"] @ recorded @ measured["evolution"].conj().T - recorded
    uniform = measured["source_uniform"]
    recorded_source = measured["source_recorded"]
    bumped = measured["source_uniform_bumped"]
    assert packet_commutator > 1.0
    assert recorded_commutator > 1.0
    assert _frobenius(packet_step) > 1.0e-3
    assert _frobenius(recorded_step) > 1.0e-3
    assert measured["outside"] > 0.5 * measured["image"]
    assert float(np.max(np.abs(uniform["force_L"]))) > 1.0
    assert float(np.max(np.abs(uniform["force_Q"]))) > 1.0
    assert float(np.max(np.abs(uniform["force_beta"]))) < ROUND
    assert float(np.max(np.abs(recorded_source["force_L"] - uniform["force_L"]))) > 0.2
    assert float(np.max(np.abs(recorded_source["force_Q"] - uniform["force_Q"]))) > 0.2
    assert float(np.max(np.abs(bumped["force_Q"] - uniform["force_Q"]))) > 1.0e-3
    assert float(np.max(np.abs(uniform["K"] - bumped["K"]))) < ROUND
    for name, column_name in (("L", "force_L"), ("Q", "force_Q"), ("beta", "force_beta")):
        gap = measured["dense_uniform"]["nodal"][name] - uniform[column_name]
        assert float(np.max(np.abs(gap))) < ROUND
    massless_images = apply_dirac(
        measured["fine"].phi0, measured["fine"].phi1, measured["system"].length_density,
        measured["fine"].Q, measured["system"].shift, 0, measured["system"].momentum,
    )
    massive_images = apply_dirac(
        measured["fine"].phi0, measured["fine"].phi1, measured["system"].length_density,
        measured["fine"].Q, measured["system"].shift, measured["system"].kappa,
        measured["system"].momentum,
    )
    assert float(np.max(np.abs(massless_images[0] - massive_images[0]))) > 1.0e-6
    assert float(np.max(np.abs(uniform["S1"]))) < ROUND
    assert float(np.max(np.abs(uniform["Pmom"]))) < ROUND


def test_projected_half_identity_commutes_while_the_quadrature_embedding_does_not():
    """H_G acts on 2 n_f. The embedded band projector is a different operator."""
    measured = _measure()
    system = measured["system"]
    hamiltonian = measured["hamiltonian"]
    band = measured["covariance_band"]
    points = measured["points"]
    isometry = measured["grid"].U_f
    embedding = np.zeros((2 * points, 2 * NF), dtype=complex)
    embedding[:points, :NF] = isometry
    embedding[points:, NF:] = isometry
    projected = embedding.conj().T @ hamiltonian @ embedding
    half_identity = 0.5 * np.eye(2 * NF, dtype=complex)
    embedded = embedding @ half_identity @ embedding.conj().T
    eye = np.eye(NF, dtype=complex)
    zero = np.zeros((NF, NF), dtype=complex)
    phi0 = np.concatenate((eye, zero), axis=1)
    phi1 = np.concatenate((zero, eye), axis=1)
    image0, image1 = apply_dirac(
        isometry @ phi0, isometry @ phi1, system.length_density, measured["fine"].Q,
        system.shift, system.kappa, system.momentum,
    )
    pulled = np.vstack((isometry.conj().T @ image0, isometry.conj().T @ image1))
    projected_step = expm(-1j * projected * DT)
    projected_change = projected_step @ half_identity @ projected_step.conj().T - half_identity
    quadrature_change = measured["evolution"] @ band @ measured["evolution"].conj().T - band
    assert projected.shape == (2 * NF, 2 * NF)
    assert hamiltonian.shape == (2 * points, 2 * points)
    assert float(np.max(np.abs(embedding.conj().T @ embedding - np.eye(2 * NF)))) < ROUND
    assert float(np.max(np.abs(pulled - projected))) < ROUND
    assert _frobenius(_commutator(projected, half_identity)) < ROUND
    assert _frobenius(projected_change) < ROUND
    assert _frobenius(embedded - band) < ROUND
    assert _frobenius(_commutator(hamiltonian, embedded)) > 1.0
    assert _frobenius(quadrature_change) > 1.0e-3
    assert _max_force(measured["source_band"]) < ROUND
    assert _max_force(measured["source_band_bumped"]) < ROUND
    assert abs(float(np.trace(band).real) - float(NF)) < ROUND


def _hamilton_increment(system, evolved, ansatz):
    """Terms present once chi or the momenta leave the initial radius ansatz."""
    force_r, f_chi = partial_F(evolved.r, system.A, system.C_W)
    f_chi = float(f_chi)
    pi = evolved.p_r - (force_r / f_chi) * evolved.p_chi
    stiffness = float(feedback_Z(system.A))
    potential = feedback_V(
        evolved.r, evolved.chi, system.A, system.C_W, system.C_F, system.flux,
    )
    potential_0 = feedback_V(
        ansatz.r, ansatz.chi, system.A, system.C_W, system.C_F, system.flux,
    )
    shell = feedback_F(evolved.r, evolved.chi, system.A, system.C_W)
    shell_0 = feedback_F(ansatz.r, ansatz.chi, system.A, system.C_W)
    derivative = system.derivative
    return (
        evolved.p_Q * evolved.p_chi / (2.0 * f_chi)
        + pi ** 2 / (4.0 * stiffness * evolved.Q)
        - evolved.Q * (potential - potential_0)
        - 2.0 * (derivative @ ((derivative @ (shell - shell_0)) / evolved.Q))
    )


def test_negative_bracket_is_not_a_lorentzian_chart_exit():
    """G>0 is the initial convexity hypothesis. The evolved chart is r, Q, L > 0."""
    text = weak.__doc__
    assert "constant ``Q``" in text
    assert "chi = p_r = p_chi = 0" in text
    assert "G > 0" in text
    assert "strictly convex" in text
    measured = _measure()
    grid = measured["grid"]
    system = measured["system"]
    q0 = float(system.calibration["b0"] / system.calibration["a0"])
    angle = 2.0 * np.pi * grid.xi_g / grid.length
    evolved = galerkin.blank_state(grid, np.zeros((NF, 1)), np.zeros((NF, 1)))
    evolved.r = np.full(grid.ng, 2.0)
    evolved.Q = np.full(grid.ng, q0)
    evolved.chi = 0.15 * np.cos(angle)
    evolved.p_chi = np.full(grid.ng, 0.04)
    evolved.p_r = 0.02 * np.sin(angle)
    evolved.p_Q = 0.01 * np.cos(angle)
    fine = galerkin.prolong_state(grid, evolved)
    ansatz = replace(
        fine,
        chi=np.zeros_like(fine.chi),
        p_Q=np.zeros_like(fine.p_Q),
        p_r=np.zeros_like(fine.p_r),
        p_chi=np.zeros_like(fine.p_chi),
    )
    blank = galerkin.blank_state(grid, np.zeros((NF, 1)), np.zeros((NF, 1)))
    blank.r = np.full(grid.ng, 2.0)
    blank.Q = np.full(grid.ng, q0)
    fine_blank = galerkin.prolong_state(grid, blank)
    increment = hamilton_constraint(system, fine) - hamilton_constraint(system, ansatz)
    ansatz_increment = hamilton_constraint(system, fine_blank) - hamilton_constraint(system, fine_blank)
    rho = np.full(system.points, -1.0)
    chart = controls.chart_from_fields(grid, {"Q": fine.Q, "rho": rho, "r": fine.r})
    lapse = fine.r * system.length_density
    assert chart_failure(system, fine) is None
    assert float(np.min(fine.r)) > 0.0
    assert float(np.min(fine.Q)) > 0.0
    assert float(np.min(lapse)) > 0.0
    assert chart["positive_r"] is True and chart["positive_Q"] is True
    assert chart["positive_G"] is False
    assert float(chart["G_min"]) < 0.0
    assert float(np.max(np.abs(
        controls.bracket_g(fine.Q, rho, grid) - controls.bracket_g(ansatz.Q, rho, grid)
    ))) < ROUND
    assert float(np.max(np.abs(ansatz_increment))) < ROUND
    assert float(np.max(np.abs(increment))) > 1.0
    assert float(np.max(np.abs(increment - _hamilton_increment(system, fine, ansatz)))) < ROUND
    assert "min G" in controls.CRITERION["chart_stop"]
    assert "when G, r, and Q are positive" in controls.CRITERION["evolution_gate"]


def test_window_proxies_do_not_erase_a_balanced_exchange_or_a_maintained_drift():
    content = np.array([[1.0, 1.0, 0.0, 0.0], [1.0, 1.0, 0.0, 0.0]])
    flux = np.array([[4.0, -4.0, 0.0, 0.0], [4.0, -4.0, 0.0, 0.0]])
    split = controls.assess_windows(content, flux)
    assert float(np.max(np.abs(np.sum(flux, axis=1)))) < ROUND
    assert float(np.max(np.abs(flux[:, 0]))) > 1.0
    assert float(np.min(content[:, 2:])) == 0.0
    assert split["end_share"] == 0.5
    assert split["renewed"] is False
    assert split["stable_throughflow"] is False
    assert split["candidate_regime"] is False
    assert split["packet_flux_end"] == 0.0
    drift = controls.assess_windows(
        np.array([[6.096, 2.0, 1.0, 0.9], [6.213, 2.0, 1.0, 0.8]]),
        np.array([[-5.5, 3.0, 1.5, 1.0], [-5.534, 3.2, 1.4, 0.934]]),
    )
    assert drift["reversal"] == 0.0
    assert drift["leader_content_change"] > 0.1
    assert drift["maintained"] is True
    assert drift["renewed"] is False
    assert drift["stable_throughflow"] is False
    assert drift["candidate_regime"] is False
    assert drift["packet_flux_end"] < -1.0


def _omega_upper(fermions, q_min, calibration, length):
    """max(L)/min(Q) bound. The two extrema need not sit on the same node."""
    wavenumber = (fermions / 2.0 - 0.5) * 2.0 * np.pi / length
    scale = float(calibration["Omega"]) ** 4.0
    lapse = float(calibration["b0"]) * scale
    shift = abs(float(calibration["beta0"])) * scale
    return (lapse / float(q_min) + shift) * wavenumber + lapse


def test_recorded_continuation_keeps_its_signal_while_G_stays_positive():
    record = json.loads(RECORD.read_text())
    assert record["conclusion"]["maintained_structure"] is True
    assert record["conclusion"]["renewed_structure"] is False
    assert record["conclusion"]["one_drift"] is True
    assert record["conclusion"]["candidate_regime"] is False
    limit = 2.0 * np.sqrt(2.0)
    for name in ("continuation_nf512", "continuation_nf256"):
        stage = record["stages"][name]
        start = stage["checkpoints"][0]
        end = stage["checkpoints"][-1]
        structure = stage["structure"]
        ratio = stage["balance"]["versus_field_change"]["balance_over_exchange"]
        assert stage["attained_time"] == 0.1
        assert stage["stop_reason"] is None
        assert stage["chart_stop"] is False
        assert stage["positive_G"] is True
        assert stage["G_min_end"] > 0.0
        assert end["G_min"] < start["G_min"]
        assert end["Q_min"] < start["Q_min"]
        assert end["positive_r"] is True and end["positive_Q"] is True
        assert start["unitarity"] < end["unitarity"] < 1.0e-8
        assert abs(float(ratio)) < 1.0e-6
        assert structure["maintained"] is True
        assert structure["renewed"] is False
        assert structure["reversal"] == 0.0
        assert structure["leader_content_change"] > 0.1
        assert structure["stable_throughflow"] is False
        assert structure["packet_flux_end"] < -1.0
        assert structure["reservoir_flux_end"] > 1.0
        assert abs(structure["flux_sum_end"]) < 1.0e-12
        assert end["chi_max"] > 80.0
    fine = record["stages"]["continuation_nf512"]
    coarse = record["stages"]["continuation_nf256"]
    assert fine["full_h_end"] < 1.0e-3
    assert coarse["full_h_end"] > 1.0e-2
    calibration = measured_calibration()
    for fermions, stage in ((512, fine), (256, coarse)):
        frequency = _omega_upper(fermions, stage["Q_min_end"], calibration, 8.0)
        assert stage["dt"] * frequency < limit
        assert stage["dt"] < 1.4 / frequency


def measured_calibration():
    return _measure()["system"].calibration


def test_half_identity_fits_the_declared_finite_gaussian_without_a_new_action():
    measured = _measure()
    identity = measured["covariance_identity"]
    accepted = _covariance(identity)
    assert accepted.shape == identity.shape
    report = influence(identity, measured["evolution"], np.eye(identity.shape[0]))
    assert np.isfinite(report["log_modulus"])
    text = coupling.__doc__
    assert "Gamma_one" in text
    assert "canonical" in text and "finite covariance" in text
    assert coupling.SEA_SUBTRACTED is False
    assert coupling.GAMMA_REST_INCLUDED is False
    assert coupling.VACUUM_MATCHED_PREREQUISITE is False
    limits = measured["dense_identity"]["limits"]
    assert limits["vacuum_branch_accounted_in_this_module"] is False
    assert limits["naive_spin_identity_means_zero_renormalized_source"] is False
    assert measured["dense_identity"]["vacuum_subtracted"] is False


def test_population_control_accepts_six_weights_and_refuses_the_full_column_count():
    measured = _measure()
    grid = measured["grid"]
    full_count = 2 * measured["points"]
    band_count = 2 * NF
    owned = controls.with_occupations(grid, np.full(6, 0.5))
    assert owned.fine.occupations.shape == (6,)
    assert np.array_equal(owned.fine.occupations, np.full(6, 0.5))
    assert grid.fine.occupations.shape == (6,)
    assert not np.array_equal(grid.fine.occupations, owned.fine.occupations)
    with pytest.raises(ValueError, match="six"):
        controls.with_occupations(grid, np.full(full_count, 0.5))
    with pytest.raises(ValueError, match="six"):
        controls.with_occupations(grid, np.full(band_count, 0.5))
    assert np.array_equal(np.asarray(controls.UNIFORM_OCCUPATIONS), np.full(6, 0.5))
    assert measured["source_identity"]["force_L"].shape == (measured["points"],)
    assert measured["covariance_identity"].shape == (full_count, full_count)
    assert measured["covariance_band"].shape[0] == full_count
    assert int(measured["source_band"]["K"].shape[0]) == measured["points"]
