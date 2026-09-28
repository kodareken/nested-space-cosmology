"""Exact finite Gaussian representation of the retained proper-time profile.

This constructs a Euclidean kinetic function and its ordered Jacobian.
It does not replace the canonical real-time Hamiltonian or choose a new
continuum quantum measure. M and Lambda remain distinct normalizations.
"""
import numpy as np
from scipy.linalg import eigh, expm
from scipy.special import exp1

from .nsc_regulated import hermitian


def ein_positive(argument):
    """Entire Ein on the nonnegative real axis, with a regular value at zero."""
    x = np.asarray(argument, dtype=float)
    if not np.isfinite(x).all() or (x < 0).any():
        raise ValueError('finite nonnegative spectral arguments required')
    result = np.empty_like(x)
    small = x < .5
    # Ein(x)=sum_(k>=1) (-1)^(k+1) x^k/(k*k!).
    z = x[small]
    term, total = z.copy(), z.copy()
    for k in range(2, 32):
        term *= -z*(k-1)/(k*k)
        total += term
    result[small] = total
    z = x[~small]
    result[~small] = np.euler_gamma + np.log(z) + exp1(z)
    return result


def kinetic_representation(matrix, cutoff=1.):
    """K=D F with F=exp[-Ein(D²/Lambda²)/2], evaluated by spectral calculus."""
    d = hermitian(matrix)
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('positive finite cutoff required')
    values, vectors = eigh(d)
    argument = (values/cutoff)**2
    ein = ein_positive(argument)
    factor = np.exp(-.5*ein)
    build = lambda diagonal: (vectors*diagonal) @ vectors.conj().T
    return {'kinetic': build(values*factor), 'factor': build(factor),
            'factor_root': build(np.sqrt(factor)),
            'ein_trace': float(np.sum(ein)),
            'euclidean_kinetic_bound': float(cutoff*np.exp(-np.euler_gamma/2))}


def warped_gaussian_measure(reference, sigma, cutoff=1., normalization=1.):
    """Finite congruence D_s=W D0 W, W=exp(-sigma/2).

    With A=W F_s^(1/2), bar(chi)=bar(psi) A† and chi=A psi,
    A† D0 A=K_s. The Grassmann Jacobian is |det A|² and
    -log J=Tr(sigma)+Tr Ein(D_s²/Lambda²)/2.
    No commutation between W and F_s is assumed.
    """
    d0, generator = hermitian(reference), hermitian(sigma)
    if generator.shape != d0.shape:
        raise ValueError('warp and reference dimensions differ')
    if not np.isfinite(normalization) or normalization <= 0:
        raise ValueError('positive independent normalization scale required')
    w = expm(-generator/2)
    ds = w @ d0 @ w
    data = kinetic_representation(ds, cutoff)
    a = w @ data['factor_root']
    n = len(d0)
    sign0, log0 = np.linalg.slogdet(d0/normalization)
    signk, logk = np.linalg.slogdet(data['kinetic']/normalization)
    if sign0 == 0 or signk == 0:
        raise ValueError('zero modes require their own determinant/state prescription')
    measure_action = float(np.trace(generator).real + data['ein_trace']/2)
    shift = float(n*(np.log(normalization/cutoff)+np.euler_gamma/2))
    return {**data, 'covariant_operator': ds, 'warp': w,
            'canonical_transform': a, 'canonical_action': float(-log0),
            'kinetic_action': float(-logk), 'measure_action': measure_action,
            'normalization_shift': shift,
            'raw_heat_action_from_measure': float(-log0+measure_action-shift)}
