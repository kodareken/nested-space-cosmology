"""Exact arithmetic, derivative growth and coherent observable-error checks."""
from fractions import Fraction as Q

import numpy as np
import pytest

from recursive_horizons.nsc_ks_residual_error import (
    sqrt_upper, propagate_h2, pointwise_field_bounds,
    pure_radius_commutator_integrals, matter_error_bounds,
)
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter


def test_directed_square_roots_and_missing_residual_bounds():
    for value in (Q(0), Q(1, 10**30), Q(2, 5), Q(7), Q(10**20)):
        result = sqrt_upper(value)
        assert result**2 >= value
        if result:
            assert (result - Q(1, 2**100))**2 < value
    with pytest.raises(ValueError):
        propagate_h2((0, 0, 0), (0, None, 0), 0, 0)


def test_triangular_commutator_growth_matches_exact_comparison_system():
    # e0'=0, e1'=b1 e0, e2'=2 b1 e1+b2 e0, with constant b1,b2.
    initial, t, b1, b2 = (Q(1), Q(2), Q(3)), Q(1, 3), Q(2), Q(5)
    bound = propagate_h2(initial, (0, 0, 0), b1 * t, b2 * t)
    expected = (initial[0], initial[1] + b1*t*initial[0],
                initial[2] + 2*b1*t*initial[1] + (b1*b1*t*t + b2*t)*initial[0])
    assert bound == expected
    assert propagate_h2(bound, (0, 0, 0), b1*t, b2*t) == propagate_h2(initial, (0, 0, 0), 2*b1*t, 2*b2*t)


def test_periodic_embedding_controls_field_and_carrier_derivative():
    L, energy = Q(2, 5), Q(13, 4)
    wave = 2*np.pi/float(L)
    # Envelope error 1+0.2 exp(i k z), with explicit L2 derivative norms.
    norms = (np.sqrt(float(L)*1.04), .2*wave*np.sqrt(float(L)), .2*wave**2*np.sqrt(float(L)))
    bounds = pointwise_field_bounds(tuple(np.nextafter(v, np.inf).item() for v in norms), L, energy)
    z = np.linspace(0, float(L), 1001)
    err = 1 + .2*np.exp(1j*wave*z)
    errz = .2j*wave*np.exp(1j*wave*z) - 1j*float(energy)*err
    assert float(bounds['F']) >= np.max(abs(err))
    assert float(bounds['F_z']) >= np.max(abs(errz))


def test_coherent_matter_perturbations_are_enclosed():
    rng = np.random.default_rng(18)
    F = rng.normal(size=(4, 2, 3)) + 1j*rng.normal(size=(4, 2, 3))
    Fz = rng.normal(size=F.shape) + 1j*rng.normal(size=F.shape)
    dF, dFz = [1e-4*(rng.normal(size=F.shape) + 1j*rng.normal(size=F.shape)) for _ in range(2)]
    Qmat, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)))
    C = Qmat @ np.diag([.1, .4, .8]) @ Qmat.conj().T
    zeros = np.zeros((0, *F.shape), complex)
    params = dict(mass=1.5, angular=-2., axial_scale=.8, radius=1.4, multiplicity=5.)
    approximate = source_column_matter(F, Fz, C, zeros, zeros, **params)['action_gradient']
    exact = source_column_matter(F+dF, Fz+dFz, C, zeros, zeros, **params)['action_gradient']
    upper_norm = lambda x: np.nextafter(np.linalg.norm(x, axis=(1, 2)).max(), np.inf).item()
    bound = matter_error_bounds(upper_norm(F), upper_norm(Fz), upper_norm(dF), upper_norm(dFz),
        mass=1.5, absolute_angular=2., axial_lower=.8, radius_lower=1.4, multiplicity=5.)
    assert np.all(np.max(abs(exact-approximate), axis=0) <= [float(bound['N']), float(bound['beta'])])


def test_pure_radius_integrals_keep_cross_term_and_zero_angular_identity():
    args = ((1, 2, 3), (4, 5, 6), Q(3, 100), Q(9, 10))
    assert pure_radius_commutator_integrals(*args, 0) == (0, 0)
    first, second = pure_radius_commutator_integrals(*args, 2)
    assert first > 0 and second > first
