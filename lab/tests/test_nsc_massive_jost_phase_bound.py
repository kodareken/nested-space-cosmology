"""Analytic and independent controls for an exterior phase certificate."""
import numpy as np
import pytest
from flint import arb, ctx

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_massive_jost_modes import outgoing_ratio
from recursive_horizons.nsc_massive_jost_phase_bound import (
    _metric, _series_context, _variable, subgap_phase_bound,
)


def bound(energy=.0013248831260437577, radius=60., sign=1):
    mass, angular = np.pi/2, sign*np.sqrt(5.)
    value, _, coefficients = outgoing_ratio(energy, mass, angular, radius, 8)
    return subgap_phase_bound(energy, mass, angular, radius, coefficients,
                              np.angle(value))


def test_original_low_energy_pair_has_small_exterior_bound():
    for sign in (-1, 1):
        result = bound(sign=sign)
        assert restored_upper(result['combined_phase_initial_error_upper']) < arb('1e-13')
        assert restored_upper(result['phase_derivative_lower']) > 3
        assert result['all_exterior_radii'] is True
        assert result['physical_rho1_source_error'] is None


def test_exact_zero_energy_zero_angular_phase_control():
    # theta=0 is an exact solution for any positive A(r) if E=ell=0.
    result = subgap_phase_bound(0., 1., 0., 60., [1.+0j], 0.)
    assert restored_upper(result['combined_phase_initial_error_upper']).is_zero()


def test_polynomial_error_and_initializer_rounding_are_both_enclosed():
    # Exact solution remains theta=0. An intentionally wrong ratio polynomial
    # has phi(R)=atan(1e-4/R); no numerical reference solver supplies the answer.
    radius = 60.
    coefficients = np.array([1.+0j, 1e-4j])
    phase = float(np.arctan(1e-4/radius))
    result = subgap_phase_bound(0., 1., 0., radius, coefficients, phase)
    total = restored_upper(result['combined_phase_initial_error_upper'])
    assert total >= abs(arb(phase))
    assert total < arb('1e-5')
    shifted = subgap_phase_bound(0., 1., 0., radius, coefficients, phase+1e-7)
    assert restored_upper(shifted['initializer_rounding_upper']) > arb('9e-8')
    assert restored_upper(shifted['combined_phase_initial_error_upper']) >= abs(arb(phase+1e-7))


def test_phase_rhs_is_original_complex_riccati_rhs():
    # Independent direct check of the reduction, with both angular signs.
    for sign in (-1, 1):
        radius, energy, mass, theta = 10., .4, 1.2, -.31
        angular = sign*2.1
        a = 1-3*((1+radius**2)*np.arctan(1/radius)-radius)
        v1, v2 = angular*np.sqrt(a)/np.sqrt(1+radius**2), mass*np.sqrt(a)
        z = np.exp(1j*theta)
        dr = (1j*(v1+1j*v2)-2j*energy*z+1j*(v1-1j*v2)*z*z)/a
        phase_rhs = 2/a*(v1*np.cos(theta)+v2*np.sin(theta)-energy)
        assert abs(dr/(1j*z)-phase_rhs) < 1e-14


def test_original_series_refinement_uses_two_independent_enclosures():
    coarse, fine = bound(radius=60.), bound(radius=120.)
    assert restored_upper(fine['combined_phase_initial_error_upper']) < restored_upper(coarse['combined_phase_initial_error_upper'])


def test_arb_parameter_intervals_cover_constant_exact_solution():
    result = subgap_phase_bound(0., arb('1 +/- 0.01'), 0., 60., [1.+0j], 0.)
    assert restored_upper(result['combined_phase_initial_error_upper']).is_zero()


def test_metric_tail_encloses_closed_form_derivatives():
    # The proof uses a power-series majorant. This control computes derivatives
    # independently from the exact arctan expression at several positive centers.
    with _series_context(256, 12):
        upper = arb(1)/60
        for denominator in (60, 90, 200):
            x = _variable(arb(1)/denominator, 12)
            exact = 1-3*((1+1/(x*x))*x.atan()-1/x)
            enclosed = _metric(x, upper, 12, 12, enclose_tail=True)
            for k in range(13):
                assert enclosed[k].contains(exact[k])


@pytest.mark.parametrize('energy,mass,radius,coefficients', [
    (1., 1., 60., [1.+0j]), (1.1, 1., 60., [1.+0j]),
    (-.1, 1., 60., [1.+0j]), (0., 1., 1., [1.+0j]),
    (0., 1., 60., [complex('nan')]),
])
def test_inadmissible_input_cannot_emit_bound(energy, mass, radius, coefficients):
    with pytest.raises(ValueError):
        subgap_phase_bound(energy, mass, 0., radius, coefficients, 0.)


def test_wrong_phase_chart_and_large_defect_are_rejected():
    with pytest.raises(ArithmeticError, match='phase chart'):
        subgap_phase_bound(0., 1., 0., 60., [-1.+0j], np.pi)
    with pytest.raises(ArithmeticError, match='tube'):
        subgap_phase_bound(0., 1., 0., 60., [1.+0j, .1j], 0.)


def test_flint_precision_and_series_cap_are_restored():
    before = (ctx.prec, ctx.cap)
    bound()
    assert (ctx.prec, ctx.cap) == before
    with pytest.raises(ArithmeticError):
        subgap_phase_bound(0., 1., 0., 60., [-1.+0j], 0.)
    assert (ctx.prec, ctx.cap) == before
