"""State-defined Killing-energy transfer on the existing static Dirac lattice.

Units hbar=c=a_throat=1. One angular channel, compact axial geometry.
The regional energy assigns half of each crossing link to each region.
This is neither a baryonic sector definition nor an expanding metric solver.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh
from .nsc_shape_response import staggered_operator


def dirac_hamiltonian(metric, kappa=1, eta=.5):
    a, _ = staggered_operator(metric, kappa, eta)
    zero = np.zeros_like(a)
    return np.block([[zero, a.T.conj()], [a, zero]])


def positions(metric):
    """Upper components on nodes; lower components on edge centers."""
    return np.concatenate((metric.x, metric.x + metric.spacing / 2))


def regional_operators(h, mask):
    """Return h_A,h_B,h_link and J_A, positive for energy entering A.

    h_A={P,H}/2, J_A=i[H,h_A]=i[H²,P]/2. Bare energies plus h_link
    also sum to H. Their currents differ by interaction-energy storage.
    """
    h = np.asarray(h)
    p = np.asarray(mask, dtype=float)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or p.shape != (len(h),):
        raise ValueError("square Hamiltonian and matching region required")
    if not np.allclose(h, h.T.conj(), rtol=0, atol=1e-12):
        raise ValueError("Hermitian Hamiltonian required")
    if not np.isin(p, [0, 1]).all() or not 0 < p.sum() < len(p):
        raise ValueError("nontrivial orthogonal spatial partition required")
    bare_a = p[:, None] * h * p[None, :]
    bare_b = (1-p[:, None]) * h * (1-p[None, :])
    link = h - bare_a - bare_b
    half = (p[:, None] * h + h * p[None, :]) / 2
    current = 1j * (h @ half - half @ h)
    return {"energy_a": half, "energy_b": h-half,
            "bare_a": bare_a, "bare_b": bare_b, "link": link,
            "current_a": current}


def split_cut_currents(current, coordinates, period):
    """Split the two short-range boundary currents on this compact circle.

    A is x>=0; its boundaries are the throat x=0 and the compact seam.
    Every nonzero current entry crossing the seam has raw distance >L/2.
    This routine is for the local staggered stencil, not arbitrary dense H.
    """
    seam = np.abs(coordinates[:, None] - coordinates[None, :]) > period / 2
    return {"throat": np.where(seam, 0, current),
            "seam": np.where(seam, current, 0)}


def vacuum(h):
    energies, vectors = eigh(h, check_finite=False, driver="evd")
    if np.min(np.abs(energies)) < 1e-10:
        raise ValueError("zero-mode occupation must be specified explicitly")
    negative = vectors[:, energies < 0]
    return energies, vectors, negative @ negative.T.conj()


def expectation(covariance, operator):
    """C_ij=<c_j† c_i>; quadratic expectation Tr(C operator)."""
    return float(np.einsum("ij,ji->", covariance, operator).real)


def positive_packet(metric, energies, vectors, center=-.6, width=.35, momentum=4.):
    """Smooth finite-energy excitation of ONE mode above the static vacuum.

    A compact C-infinity bump is projected to the positive spectral subspace.
    Projection is nonlocal; no claim of compact support or vacuum generation.
    The +sigma2 seed travels toward increasing x in the continuum principal part.
    """
    x = positions(metric)
    y = (x-center)/width
    envelope = np.zeros_like(x)
    inside = np.abs(y) < 1
    envelope[inside] = np.exp(-1/(1-y[inside]**2))
    spin = np.concatenate((np.ones(metric.points), 1j*np.ones(metric.points)))
    seed = envelope * np.exp(1j*momentum*x) * spin * np.sqrt(metric.spacing)
    positive = vectors[:, energies > 0]
    packet = positive @ (positive.T.conj() @ seed)
    norm = np.linalg.norm(packet)
    if norm < 1e-12:
        raise ValueError("positive-energy seed vanished")
    return packet/norm


def evolve_packet(energies, vectors, initial, times):
    coefficients = vectors.T.conj() @ initial
    return vectors @ (coefficients[:, None] * np.exp(-1j*energies[:, None]*np.asarray(times)))


def packet_expectations(states, operator):
    return np.einsum("it,it->t", states.conj(), operator @ states).real


def integrated_current(energies, vectors, initial, current, time):
    """Exact finite spectral time integral; used against a quadrature control."""
    coefficient = vectors.T.conj() @ initial
    frequencies = energies[:, None]-energies[None, :]
    integral = time * np.exp(.5j*frequencies*time)*np.sinc(frequencies*time/(2*np.pi))
    rotated = vectors.T.conj() @ current @ vectors
    return float(np.sum(coefficient.conj()[:, None]*coefficient[None, :]*rotated*integral).real)


def schur_state_data(h, mask, energy, covariance):
    """Retarded elimination plus initial-state data; neither replaces the other.

    i dpsi_A/dt=H_AA psi_A -i int B exp[-iH_BB(t-s)] B† psi_A(s) ds
                  + B exp[-iH_BB t] psi_B(0).
    The last term and initial cross-correlations are needed for state averages.
    """
    p = np.asarray(mask, dtype=bool)
    aa, bb, ab = h[np.ix_(p,p)], h[np.ix_(~p,~p)], h[np.ix_(p,~p)]
    sigma = ab @ np.linalg.solve(energy*np.eye(len(bb))-bb, ab.T.conj())
    reduced = np.linalg.inv(energy*np.eye(len(aa))-aa-sigma)
    direct = np.linalg.inv(energy*np.eye(len(h))-h)[np.ix_(p,p)]
    cbb = covariance[np.ix_(~p,~p)]
    noise_equal_time = ab @ cbb @ ab.T.conj()
    return {"self_energy_norm": float(np.linalg.norm(sigma)),
            "direct_reduced_error": float(np.max(np.abs(reduced-direct))),
            "initial_noise_norm": float(np.linalg.norm(noise_equal_time)),
            "initial_cross_covariance_norm": float(np.linalg.norm(covariance[np.ix_(p,~p)]))}
