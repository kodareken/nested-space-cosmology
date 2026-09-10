"""Coupled spherical lapse/shift/radial Dirac source vertices.

The existing factorized spatial discretization is retained. Continuum
diffeomorphism identities are separate from a finite Fourier product rule.
All time blocks act in the canonical half-density frame u=r sqrt(q) psi.
"""
import numpy as np
from scipy.linalg import block_diag

from .nsc_covariant_operator import SIGMA1, SIGMA2, SIGMA3
from .nsc_influence import canonical_hamiltonian, _covariance
from .nsc_regulated import hermitian


def shift_matrix(metric, shift):
    b = np.asarray(shift, dtype=complex)
    if b.shape != (metric.points,) or not np.isfinite(b).all():
        raise ValueError('finite shift on the existing spatial grid required')
    p = metric.momentum_matrix
    return np.kron(np.eye(2), .5*(b[:, None]*p+p*b[None, :]))


def adm_hamiltonian(metric, shift, kappa=1):
    if np.iscomplexobj(shift) and np.max(abs(np.imag(shift))) > 0:
        raise ValueError('physical Lorentzian shift must be real')
    return canonical_hamiltonian(metric, kappa)-shift_matrix(metric, shift)


def perpendicular_tangents(metric, direction, kappa=1):
    """Exact derivatives of the retained finite H_perp on a log-metric path.

    direction=(delta log N, delta log q, delta log r); second derivatives
    include the contacts from this coordinate choice.
    """
    direction = np.asarray(direction, dtype=float)
    if direction.shape != (3, metric.points) or not np.isfinite(direction).all():
        raise ValueError('three finite log-metric directions required')
    a, b, c = direction
    p, n, q, r = metric.momentum_matrix, metric.lapse, metric.radial_scale, metric.sphere_radius
    h = canonical_hamiltonian(metric, kappa)
    aa = np.diag(np.r_[a, a])
    root = np.sqrt(np.r_[n, n])
    def core(order):
        v = (-b)**order/q
        mass = kappa*(-c)**order/r
        return (np.kron(SIGMA2, .5*(v[:, None]*p+p*v[None, :]))
                + np.kron(SIGMA1, np.diag(mass)))
    v1, v2 = (root[:, None]*core(order)*root[None, :] for order in (1, 2))
    h1 = .5*(aa@h+h@aa)+v1
    h2 = .25*(aa@aa@h+h@aa@aa)+.5*aa@h@aa+aa@v1+v1@aa+v2
    return h, h1, h2


def covariant_history_operator(metrics, euclidean_shifts, time_momentum,
                               directions=None, shift_directions=None, kappa=1):
    """D_E=W beta[P_tau-B_E-i H_perp]W, W=N^-1/2.

    Real Euclidean shifts give a Hermitian matrix. Physical Lorentzian
    shifts enter the analytic prescription as b_E=-i beta_L; this routine
    does not apply a Hermitian eigensolver after such a continuation.
    The sign of i H differs from the old static fiber by constant beta
    conjugation, preserving that fiber's modulus.
    """
    if not metrics or len(metrics) != len(euclidean_shifts):
        raise ValueError('one shift array per time slice required')
    first = metrics[0]
    for m in metrics:
        if (m.points, m.length, m.eta) != (first.points, first.length, first.eta):
            raise ValueError('a common spatial domain and basis are required')
    pt = hermitian(time_momentum)
    if pt.shape != (len(metrics), len(metrics)):
        raise ValueError('time momentum and slice count differ')
    size = 2*first.points
    beta = np.kron(np.eye(len(metrics)), np.kron(SIGMA3, np.eye(first.points)))
    if directions is None:
        directions = np.zeros((len(metrics), 3, first.points))
    directions = np.asarray(directions, dtype=float)
    if directions.shape != (len(metrics), 3, first.points):
        raise ValueError('three log-metric directions per time slice required')
    if shift_directions is None:
        shift_directions = np.zeros((len(metrics), first.points))
    if np.shape(shift_directions) != (len(metrics), first.points):
        raise ValueError('one additive shift direction per time slice required')
    triples = [perpendicular_tangents(m, d, kappa) for m, d in zip(metrics, directions)]
    h, h1, h2 = (block_diag(*(row[j] for row in triples)) for j in range(3))
    b = block_diag(*(shift_matrix(m, s) for m, s in zip(metrics, euclidean_shifts)))
    b1 = block_diag(*(shift_matrix(m, s) for m, s in zip(metrics, shift_directions)))
    w = np.concatenate([1/np.sqrt(np.r_[m.lapse, m.lapse]) for m in metrics])
    aa = np.diag(np.concatenate([np.r_[d[0], d[0]] for d in directions]))
    sandwich = lambda x: w[:, None]*(beta@x)*w[None, :]
    d = sandwich(np.kron(pt, np.eye(size))-b-1j*h)
    v1, v2 = sandwich(-b1-1j*h1), sandwich(-1j*h2)
    d1 = -.5*(aa@d+d@aa)+v1
    d2 = .25*(aa@aa@d+d@aa@aa)+.5*aa@d@aa-aa@v1-v1@aa+v2
    return {'operator': d, 'first': d1, 'second': d2}


def canonical_force_gradients(metric, shift, covariance, kappa=1):
    """Nodal derivatives of Tr(C H_ADM), holding the canonical C fixed.

    Divide by the spatial spacing before treating these as functional
    derivative densities. No vacuum subtraction or angular sum is implied.
    """
    c = _covariance(covariance)
    h = canonical_hamiltonian(metric, kappa)
    if c.shape != h.shape:
        raise ValueError('covariance and spatial operator dimensions differ')
    n, q, r, p = metric.lapse, metric.radial_scale, metric.sphere_radius, metric.momentum_matrix
    root = np.sqrt(np.r_[n, n])
    weighted = root[:, None]*c*root[None, :]
    def spin_trace(matrix, gamma):
        blocks = matrix.reshape(2, metric.points, 2, metric.points)
        return np.einsum('ba,aibj->ij', gamma, blocks)
    trace0 = spin_trace(c, np.eye(2))
    trace1, trace2 = spin_trace(weighted, SIGMA1), spin_trace(weighted, SIGMA2)
    gn = .5*np.diag(c@h+h@c).reshape(2, metric.points).sum(axis=0).real/n
    gb = -.5*np.diag(p@trace0+trace0@p).real
    gq = -.5*np.diag(p@trace2+trace2@p).real/q**2
    gr = -kappa*np.diag(trace1).real/r**2
    return {'N': gn, 'beta': gb, 'q': gq, 'r': gr,
            'energy': float(np.trace(c@adm_hamiltonian(metric, shift, kappa)).real)}
