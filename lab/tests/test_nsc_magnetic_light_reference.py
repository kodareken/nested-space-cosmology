"""Charged spectral continuation and same-action trace/conservation checks."""
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

from recursive_horizons.nsc_magnetic_light_reference import (
    magnetic_light_spectrum, derivative_series_tail_bound, homogeneous_magnetic_restoration,
)
from recursive_horizons.nsc_angular_stress import (
    static_cylinder_density, profile_jets, curvature_squared_tensor, physical_source,
)
from recursive_horizons.nsc_horizon_source import conformal_stress


@pytest.fixture(scope='module')
def charged():
    return magnetic_light_spectrum(4, series_tolerance='1e-35', decimal_precision=80)


def continued_zeta(s, charge, terms):
    """Independent off-pole evaluation before taking the special-point limit.

    Differentiating the rising factorial numerically, rather than inserting
    its derivative formula, independently tests the pole-times-zero term.
    """
    b = mp.mpf(charge)/2; a = b+1
    return 4*mp.fsum(mp.rf(s, j)/mp.factorial(j)*b**(2*j)*mp.zeta(2*s+2*j-1, a)
                    for j in range(terms+1))


def test_charged_special_values_and_derivative_continuation(charged):
    assert charged.z_minus_one == Fraction(1, 30)-Fraction(16, 6)
    assert charged.z_nonzero_zero == -Fraction(13, 3)
    assert charged.harmonic_number == Fraction(25, 12)
    assert charged.derivative_tail_bound < Fraction('1e-35')
    with mp.workdps(80):
        h = mp.mpf('1e-15'); terms = charged.last_series_index+12
        minus_one = continued_zeta(-1+1j*h, 4, terms)
        zero = continued_zeta(1j*h, 4, terms)
        one_finite = continued_zeta(1+1j*h, 4, terms)-2/(1j*h)
        assert abs(minus_one.real-mp.mpf(-79)/30) < mp.mpf('1e-26')
        assert abs(minus_one.imag/h-charged.zprime_minus_one) < mp.mpf('1e-26')
        assert abs(zero.real+mp.mpf(13)/3) < mp.mpf('1e-26')
        assert abs(one_finite.real-charged.finite_part_one) < mp.mpf('1e-26')


def test_series_remainder_and_exact_neutral_normalization(charged):
    coarse = magnetic_light_spectrum(4, series_tolerance='1e-22', decimal_precision=60)
    with mp.workdps(80):
        bound = mp.mpf(coarse.derivative_tail_bound.numerator)/coarse.derivative_tail_bound.denominator
        # Omitted terms are positive and SUBTRACTED in Zprime.
        difference = coarse.zprime_minus_one-charged.zprime_minus_one
        assert 0 < difference < bound
        assert abs(coarse.cylinder_density()-charged.cylinder_density()) < coarse.cylinder_series_error_bound()
    neutral = magnetic_light_spectrum(0)
    assert neutral.z_minus_one == Fraction(1, 30)
    assert neutral.z_nonzero_zero == -Fraction(1, 3)
    assert neutral.derivative_tail_bound == 0 and neutral.series_terms == 0
    assert abs(float(neutral.cylinder_density())-static_cylinder_density()) < 3e-15
    assert abs(float(neutral.fourth_order_harmonic)-float(mp.euler)/2) < 3e-15
    assert derivative_series_tail_bound(4, charged.last_series_index+2) < charged.derivative_tail_bound
    with pytest.raises(ValueError, match='retained'):
        magnetic_light_spectrum(2)
    with pytest.raises(ArithmeticError, match='tail tolerance'):
        magnetic_light_spectrum(4, max_terms=1)


def physical_lll_geometry(theta, charge):
    """Existing conformal-stress owner, transformed into its physical KS frame."""
    W, W1, W2, _, _ = profile_jets(theta)
    r = 1/np.sin(theta); s1 = -np.cos(theta)/np.sin(theta); s2 = r*r
    A = -r*r*W; Ap = -W1-2*s1*W
    App = -(W2+2*W1*s1+2*W*s2)/(r*r)
    uu, uv, vv = conformal_stress(A, Ap, App, 0., 0., central_charge=charge)
    return np.array([(uu-2*uv+vv)/(-A), (uu+2*uv+vv)/(-A), 0.])/(4*np.pi*r*r)


def test_homogeneous_trace_and_lll_count_once(charged):
    for theta in (1.4, np.pi/2, 2., 3*np.pi/4):
        row = homogeneous_magnetic_restoration(theta, charged)
        lll = physical_lll_geometry(theta, 4)
        full_trace = row['rho']+lll[0]-row['p_parallel']-lll[1]-2*row['p_sphere']
        assert abs(full_trace-row['neutral_gravitational_trace']-row['gauge_trace']) < 3e-14
        assert abs(lll[0]-lll[1]-row['excluded_physical_LLL_trace']) < 3e-14
        assert row['scope']['LLL_state_or_geometric_source_included'] is False
        assert row['scope']['compact_complement_included'] is False
    neutral = magnetic_light_spectrum(0)
    theta = 3*np.pi/4; jets = profile_jets(theta); e4, p4 = curvature_squared_tensor(jets)
    expected = physical_source({'q': theta, 'W_jets': jets, 'rho_coordinate': 1., 'sphere_radius': np.sqrt(2),
        'bar_rho': static_cylinder_density()+float(mp.euler)/2*e4/(480*np.pi**2),
        'bar_p_parallel': -static_cylinder_density()-float(mp.euler)/2*p4/(480*np.pi**2)})
    got = homogeneous_magnetic_restoration(theta, neutral)
    assert max(abs(got[k]-expected[k]) for k in ('rho', 'p_parallel', 'p_sphere', 'trace')) < 3e-14


def test_homogeneous_restoration_obeys_source_conservation(charged):
    # The existing homogeneous Ward equation, expressed in theta. Both the
    # nonzero-angular allocation and its sum with the physical LLL are tested.
    for theta in (1.4, np.pi/2, 2., 3*np.pi/4):
        for include_lll in (False, True):
            def tensor(x):
                row = homogeneous_magnetic_restoration(x, charged)
                value = np.array([row['rho'], row['p_parallel'], row['p_sphere']])
                return value+(physical_lll_geometry(x, 4) if include_lll else 0)
            h = 8e-5
            derivative = (tensor(theta-2*h)[0]-8*tensor(theta-h)[0]
                          +8*tensor(theta+h)[0]-tensor(theta+2*h)[0])/(12*h)
            rho, parallel, sphere = tensor(theta)
            W, W1, _, _, _ = profile_jets(theta); s1 = -np.cos(theta)/np.sin(theta)
            connection = (s1+W1/(2*W))*(rho+parallel)+2*s1*(rho+sphere)
            assert abs(derivative+connection) < 3e-9*(1+abs(derivative)+abs(connection))
