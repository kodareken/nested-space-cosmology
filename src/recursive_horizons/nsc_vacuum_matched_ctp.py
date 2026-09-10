"""Vacuum-matched canonical-state prescription in one finite Dirac channel.

The real geometric branch is calculated, not independently weighted. Its
quadratic vertex retains full frequency dependence. It is not itself a
retarded matter correlator, and no coupled metric-causality claim is made.
"""
from math import pi, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.special import erf

from .nsc_regulated import hermitian
from .nsc_influence import influence, with_real_branch_action


class VacuumMatchedDiracAction:
    """H(J)=H+J V; unit-lapse Dirac family with a common anticommuting beta."""

    def __init__(self, h, vertex, beta, cutoff=4., quadrature=40):
        self.h, self.vertex, beta = map(hermitian, (h, vertex, beta))
        if self.h.shape != self.vertex.shape or self.h.shape != beta.shape:
            raise ValueError('matching finite Dirac matrices required')
        if not np.isfinite(cutoff) or cutoff <= 0:
            raise ValueError('positive finite cutoff required')
        if not isinstance(quadrature, int) or quadrature < 16:
            raise ValueError('at least 16 proper-time/Feynman nodes required')
        identity = np.eye(len(beta))
        if not np.allclose(beta@beta, identity, atol=1e-11, rtol=0):
            raise ValueError('beta must be an involution')
        for a in (self.h, self.vertex):
            if not np.allclose(beta@a+a@beta, 0, atol=1e-11, rtol=0):
                raise ValueError('this owner requires the paired unit-lapse Dirac channel')
        self.cutoff = cutoff
        e, u = eigh(self.h, driver='evd')
        if min(abs(e)) < 1e-9:
            raise ValueError('zero-mode preparation requires separate treatment')
        rotated = u.conj().T@self.vertex@u
        self.e, self.v2 = e, abs(rotated)**2
        self.reference_covariance = u[:, e < 0]@u[:, e < 0].conj().T
        self.data = {'energies': e, 'vertex_squared': self.v2}
        self.linear = float(.5*np.dot(np.sign(e)*erf(abs(e)/cutoff), rotated.diagonal().real))
        self.constant = float(.5*np.sum(cutoff/sqrt(pi)*np.exp(-(e/cutoff)**2)
                                           + abs(e)*erf(abs(e)/cutoff)))
        self.contact = float(np.dot(self.v2.sum(axis=1), erf(abs(e)/cutoff)/(2*abs(e))))
        x, w = leggauss(quadrature)
        self.nodes, self.weights = (x+1)/2, w/2

    def euclidean_counter_kernel(self, frequency_squared):
        """B2(z)=K_heat(z)-K_canonical(z), entire in external z=nu².

        Calculated directly over 0<s<Lambda^-2, without subtracting two
        continued propagators or moving a loop contour across its poles.
        """
        z = complex(frequency_squared)
        if not np.isfinite(z):
            raise ValueError('finite squared external frequency required')
        e, v2 = self.e, self.v2
        alpha = self.nodes[:, None, None]
        mass = ((1-alpha)*e[None, :, None]**2
                + alpha*e[None, None, :]**2 + alpha*(1-alpha)*z)
        numerator = ((e[:, None]+e[None, :])**2+z)*v2
        bubble = 0j
        for u, w in zip(self.nodes, self.weights):
            integrand = np.sum(np.exp(-u*u*mass/self.cutoff**2)*numerator, axis=(1, 2))
            bubble += w*u*u*np.dot(self.weights, integrand)
        result = self.contact-bubble/(2*sqrt(pi)*self.cutoff**3)
        if not np.isfinite(result):
            raise ArithmeticError('frequency/quadrature exceeds numerical range')
        return complex(result)

    def continued_force_kernel(self, frequency):
        """Real bare geometric vertex B2(-omega²), not a retarded propagator."""
        if not np.isreal(frequency) or not np.isfinite(frequency):
            raise ValueError('real finite Lorentzian frequency required')
        value = self.euclidean_counter_kernel(-float(frequency)**2)
        if abs(value.imag) > 1e-10:
            raise ArithmeticError('continued geometric vertex lost reality')
        return float(value.real)

    def relative_branch_action(self, period, positive_fourier):
        """S_B[J]-S_B[0] through O(J²), for a real finite Fourier history.

        J(t)=J0+sum_(n>0)[Jn exp(-i omega_n t)+Jn* exp(i omega_n t)].
        No expansion in omega/Lambda is made. The constant vacuum energy
        is separately available; its common reference phase cancels here.
        """
        if not np.isfinite(period) or period <= 0:
            raise ValueError('positive common history period required')
        coefficients = dict(positive_fourier)
        for n, value in coefficients.items():
            if not isinstance(n, int) or n < 0 or not np.isfinite(value):
                raise ValueError('finite nonnegative Fourier modes required')
        zero = complex(coefficients.get(0, 0.))
        if abs(zero.imag) > 1e-12:
            raise ValueError('the zero Fourier coefficient must be real')
        energy = self.linear*zero.real + .5*self.continued_force_kernel(0.)*zero.real**2
        for n, value in coefficients.items():
            if n:
                energy += self.continued_force_kernel(2*pi*n/period)*abs(value)**2
        return float(-period*energy)


def vacuum_matched_influence(covariance, plus, minus, branch_plus, branch_minus):
    """Canonical state factor times the calculated real geometric phase.

    Positivity of the history kernel is preserved by diagonal phase
    multiplication. This says nothing by itself about metric causality.
    """
    canonical = influence(covariance, plus, minus)
    if not np.isreal(branch_plus) or not np.isreal(branch_minus):
        raise ValueError('real geometric branch actions required')
    delta = float(branch_plus-branch_minus)
    phase = canonical['principal_phase']
    phase = None if phase is None else float(np.angle(np.exp(1j*(phase+delta))))
    return {**canonical,
            'amplitude': with_real_branch_action(canonical['amplitude'],
                                                 float(branch_plus), float(branch_minus)),
            'canonical_amplitude': canonical['amplitude'],
            'principal_phase': phase,
            'principal_action': None if phase is None else complex(phase, -canonical['log_modulus']),
            'geometric_phase_difference': delta,
            'action_convention': 'canonical action plus a real branch difference; branch winding retained separately'}
