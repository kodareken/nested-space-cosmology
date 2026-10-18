"""Stress and metric vertices of the same minimally coupled radial Dirac action.

Modes xi=sqrt(a)*r*psi have canonical norm integral xi^dagger*xi dy.
The full angular multiplet (2*abs(k) magnetic states) is filled with the
same radial covariance. Arrays end in (spinor=2, mode_count); preceding
axes are physical sample points. No covariance or renormalization is
silently chosen. Units hbar=c=1; stress matrices scale as length^-4.
"""
from __future__ import annotations

import math

import numpy as np

from .nsc_atlas_dirac_evolution import ALPHA

SIGMA1 = np.array([[0., 1.], [1., 0.]])
IDENTITY = np.eye(2)


def mode_bilinears(values, derivatives):
    """Hermitian matrices for arbitrary pairs of radial modes, not just densities."""
    u, du = (np.asarray(v, dtype=complex) for v in (values, derivatives))
    if (u.ndim < 2 or u.shape[-2] != 2 or u.shape[-1] < 1 or du.shape != u.shape
            or not np.all(np.isfinite(u)) or not np.all(np.isfinite(du))):
        raise ValueError('finite mode arrays with identical (...,2,modes) shape required')

    def product(left, matrix, right):
        return np.einsum('...ia,ij,...jb->...ab', left.conj(), matrix, right)

    def derivative(matrix):
        return -.5j*(product(u, matrix, du)-product(du, matrix, u))

    return dict(energy_derivative=derivative(ALPHA), momentum_derivative=derivative(IDENTITY),
        angular=product(u, SIGMA1, u), norm=product(u, IDENTITY, u),
        charge_current=product(u, ALPHA, u))


def _metric(value, shape, name, positive=True):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf':
        raise ValueError(name+' must be finite real data')
    try:
        field = np.broadcast_to(np.asarray(raw, dtype=float), shape)
    except ValueError as error:
        raise ValueError(name+' does not broadcast to the mode sample axes') from error
    if not np.all(np.isfinite(field)) or (positive and np.any(field <= 0)):
        raise ValueError(name+' must be finite'+(' and positive' if positive else ''))
    return field[..., None, None]


def _angular(k):
    if isinstance(k, bool) or not isinstance(k, int) or not k:
        raise ValueError('nonzero integer angular eigenvalue required')
    return 2*abs(k)


def stress_matrices(values, derivatives, *, radial_length, radius, angular_k=1):
    """Orthonormal rho, pressures, outward T^{hat0 hat1} and charge-current matrices.

    Trace against a finite radial covariance for expectation values, or
    against its tangent/difference for a response. Global normalization,
    CAR occupation bounds, state preparation and UV subtraction are caller
    responsibilities. These are finite-mode insertions, not the anomaly.
    """
    b = mode_bilinears(values, derivatives)
    shape = b['norm'].shape[:-2]
    a = _metric(radial_length, shape, 'radial_length')
    r = _metric(radius, shape, 'radius')
    weight = _angular(angular_k)/(4*math.pi)
    qa, qi, s1 = b['energy_derivative'], b['momentum_derivative'], b['angular']
    return dict(energy_density=weight*(qa/a+angular_k*s1/r)/(a*r*r),
        radial_pressure=weight*qa/(a*a*r*r),
        angular_pressure=weight*angular_k*s1/(2*a*r**3),
        outward_energy_flux=weight*qi/(a*a*r*r),
        charge_density=weight*b['norm']/(a*r*r),
        outward_charge_current=weight*b['charge_current']/(a*r*r))


def canonical_metric_vertices(values, derivatives, *, lapse, radial_length,
                              radius, shift=0., angular_k=1):
    """First/second variations of the real canonical Hamiltonian density.

    h=d*[N*(q_alpha/a+k*sigma1/r)-beta*q_I], per coordinate dy.
    The omitted integration endpoint is -i*d*[(N/a)*xi^dagger*alpha*xi
    -beta*xi^dagger*xi]/2. Keep it for an open port. The bulk action is
    varied at fixed canonical half-density, with initial preparation held
    fixed. These bare contacts do not replace physical stress/measure
    variation or covariant continuum/interface renormalization.
    """
    b = mode_bilinears(values, derivatives)
    shape = b['norm'].shape[:-2]
    N = _metric(lapse, shape, 'lapse')
    a = _metric(radial_length, shape, 'radial_length')
    r = _metric(radius, shape, 'radius')
    beta = _metric(shift, shape, 'shift', positive=False)
    d = _angular(angular_k)
    qa, qi, s1 = b['energy_derivative'], b['momentum_derivative'], b['angular']
    zero = np.zeros_like(qa)
    first = dict(lapse=d*(qa/a+angular_k*s1/r), radial_length=-d*N*qa/a**2,
                 radius=-d*N*angular_k*s1/r**2, shift=-d*qi)
    # Both orders are supplied; absent pairs are zero rather than hidden.
    keys = ('lapse', 'radial_length', 'radius', 'shift')
    second = {left: {right: zero for right in keys} for left in keys}
    second['lapse']['radial_length'] = second['radial_length']['lapse'] = -d*qa/a**2
    second['lapse']['radius'] = second['radius']['lapse'] = -d*angular_k*s1/r**2
    second['radial_length']['radial_length'] = 2*d*N*qa/a**3
    second['radius']['radius'] = 2*d*N*angular_k*s1/r**3
    return dict(hamiltonian_density=N*first['lapse']-d*beta*qi,
                first=first, second=second)
