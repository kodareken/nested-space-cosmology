"""Manufactured radial controls for optional DOP853 panel startup steps."""
import re

import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_variational_reference as radial

CONTROLS = ControlBudget(maximum_evaluations=10_000_000, cpu_seconds=15.)
BREAKS = np.array([0., .1, .25, .5, .8, 1.])
PANELS = np.diff(BREAKS)
ADJOINT_SAMPLES = 32*len(PANELS)
RTOL, ATOL = 2e-13, 2e-15
# Affine sources on a constant positive metric. The sensitivity equation is
# linear, and mass, lapse, and both sensitivities have closed forms.
A0, RHO0, RHO1 = .25, .25, .125
P0 = np.array([.5, -.25])
P1 = np.array([.125, .0625])
CENTRE_LAPSE = -.0625


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module', autouse=True)
def report():
    yield
    print('radial startup controls:', CONTROLS.budget.evaluations,
          'charged values;', CONTROLS.active_cpu_seconds, 'active CPU seconds')


def _exact(x):
    mass = RHO0*x+RHO1*x*x/2
    lapse = CENTRE_LAPSE+A0*x
    gradient = (P0/A0+P1*(x/A0-1/A0**2))-np.exp(-A0*x)*(P0/A0-P1/A0**2)
    return mass, lapse, gradient


def _coefficients(cell, point, mass):
    partials = P0+P1*point
    return point, 1., 1., A0, RHO0+RHO1*point, partials


def _zero_coefficients(cell, point, mass):
    return 0., 1., 1., 0., 0., np.zeros(2)


def _solve(budget, fraction, coefficients, nonzero, **extra):
    return radial._radial_variations(
        BREAKS, coefficients, centre_lapse=CENTRE_LAPSE if nonzero else 0.,
        derivative_count=2, nonzero_field=nonzero, budget=budget,
        rtol=RTOL, atol=ATOL, first_step_fraction=fraction, **extra)


def _rhs_calls(result):
    assert result['field_sites'] >= ADJOINT_SAMPLES
    return result['field_sites']-ADJOINT_SAMPLES


def _watch_steps(monkeypatch):
    seen = []
    original = radial.solve_ivp

    def wrapped(fun, t_span, y0, **kwargs):
        seen.append(kwargs.get('first_step', None))
        return original(fun, t_span, y0, **kwargs)

    monkeypatch.setattr(radial, 'solve_ivp', wrapped)
    return seen


def test_linear_default_and_full_panel_match_exact_solution(budget, monkeypatch):
    steps = _watch_steps(monkeypatch)
    default = _solve(budget, None, _coefficients, True)
    default_steps = list(steps)
    steps.clear()
    full = _solve(budget, 1.0, _coefficients, True)
    mass, lapse, gradient = _exact(BREAKS[-1])
    for result in (default, full):
        assert result['mass'] == pytest.approx(mass, abs=5e-12)
        assert result['lapse'] == pytest.approx(lapse, abs=5e-12)
        np.testing.assert_allclose(result['gradient'], gradient, rtol=0., atol=5e-12)
        np.testing.assert_allclose(result['adjoint_gradient'], gradient, rtol=0., atol=5e-12)
        assert result['gradient_discrepancy'] < 1e-10
        assert result['sampled_minimum_A'] == 1.
    assert abs(default['mass']-full['mass']) < 5e-12
    assert abs(default['lapse']-full['lapse']) < 5e-12
    np.testing.assert_allclose(default['gradient'], full['gradient'], rtol=0., atol=5e-12)
    assert default_steps == [None]*len(PANELS)
    np.testing.assert_array_equal(steps, PANELS)
    default_rhs, full_rhs = _rhs_calls(default), _rhs_calls(full)
    assert full['field_sites'] < default['field_sites']
    assert full_rhs*2 <= default_rhs
    print('linear startup mass/lapse/gradient gaps:',
          abs(full['mass']-mass), abs(full['lapse']-lapse),
          float(np.max(np.abs(full['gradient']-gradient))),
          'default-full gaps:', abs(default['mass']-full['mass']),
          abs(default['lapse']-full['lapse']),
          float(np.max(np.abs(default['gradient']-full['gradient']))),
          'samples:', default['field_sites'], full['field_sites'],
          'rhs:', default_rhs, full_rhs)


def test_exact_zero_rhs_keeps_zero_state_with_fewer_samples(budget, monkeypatch):
    steps = _watch_steps(monkeypatch)
    default = _solve(budget, None, _zero_coefficients, False)
    steps.clear()
    full = _solve(budget, 1.0, _zero_coefficients, False)
    for result in (default, full):
        assert result['mass'] == result['lapse'] == 0.
        assert np.all(result['gradient'] == 0.)
        assert np.all(result['adjoint_gradient'] == 0.)
        assert result['gradient_discrepancy'] == 0.
    np.testing.assert_array_equal(steps, PANELS)
    default_rhs, full_rhs = _rhs_calls(default), _rhs_calls(full)
    assert full['field_sites'] < default['field_sites']
    assert full_rhs*3 <= default_rhs
    print('zero startup samples:', default['field_sites'], full['field_sites'],
          'rhs:', default_rhs, full_rhs)


def test_invalid_first_step_fraction_rejected_before_sampling(budget):
    def forbidden(*_args):
        raise AssertionError('radial solve started')

    for fraction in (0., -1e-12, 1.+1e-12, np.nan, np.inf, -np.inf, 'panel'):
        with pytest.raises(ValueError, match='first-step fraction'):
            radial._radial_variations(
                BREAKS, forbidden, centre_lapse=0., derivative_count=2,
                nonzero_field=False, budget=budget, rtol=RTOL, atol=ATOL,
                first_step_fraction=fraction)


def test_adjoint_orders_32_and_64_match_exact_linear_solution(budget):
    assert radial.GRADIENT_ABSOLUTE_TOLERANCE == 2e-8
    mass, lapse, gradient = _exact(BREAKS[-1])
    default = _solve(budget, 1., _coefficients, True)
    solved = {}
    for order in (32, 64):
        result = _solve(budget, 1., _coefficients, True, adjoint_order=order)
        solved[order] = result
        assert result['adjoint_order'] == order
        assert result['mass'] == pytest.approx(mass, abs=5e-12)
        assert result['lapse'] == pytest.approx(lapse, abs=5e-12)
        np.testing.assert_allclose(result['gradient'], gradient, rtol=0., atol=5e-12)
        np.testing.assert_allclose(result['adjoint_gradient'], gradient, rtol=0., atol=5e-12)
        assert result['gradient_discrepancy'] < 1e-10
    assert default['adjoint_order'] == 32
    assert default['field_sites'] == solved[32]['field_sites']
    np.testing.assert_array_equal(default['gradient'], solved[32]['gradient'])
    assert solved[64]['field_sites'] == solved[32]['field_sites']+(64-32)*len(PANELS)
    print('adjoint order 32/64 discrepancies:',
          solved[32]['gradient_discrepancy'], solved[64]['gradient_discrepancy'],
          'samples:', solved[32]['field_sites'], solved[64]['field_sites'])


def test_invalid_adjoint_order_rejected_before_sampling(budget):
    def forbidden(*_args):
        raise AssertionError('radial solve started')

    for order in (7, 257, 0, -4, True, False, 32.0, np.nan, np.inf, '32'):
        with pytest.raises(ValueError, match='adjoint quadrature order'):
            radial._radial_variations(
                BREAKS, forbidden, centre_lapse=0., derivative_count=2,
                nonzero_field=False, budget=budget, rtol=RTOL, atol=ATOL,
                adjoint_order=order)


def test_underresolved_adjoint_failure_names_the_worst_component(budget):
    rtol, atol = 1.25e-11, 3.5e-14

    def coefficients(cell, point, mass):
        return 0., 1., 1., 0., .25, np.array([np.cos(4*np.pi*point), .01])

    with pytest.raises(ValueError, match='declared tolerance') as caught:
        radial._radial_variations(
            np.array([0., 1.]), coefficients, centre_lapse=0., derivative_count=2,
            nonzero_field=True, budget=budget, rtol=rtol, atol=atol,
            first_step_fraction=1., adjoint_order=8)
    fields = dict(re.findall(r'([A-Za-z_]+)=([^,\s]+)', str(caught.value)))
    discrepancy = float(fields['discrepancy'])
    forward, adjoint = float(fields['forward']), float(fields['adjoint'])
    assert discrepancy > radial.GRADIENT_ABSOLUTE_TOLERANCE
    assert int(fields['worst_component']) == 0
    assert int(fields['order']) == 8
    assert abs(forward-adjoint) == pytest.approx(discrepancy, abs=0., rel=0.)
    assert float(fields['largest_absolute_gradient']) == pytest.approx(.01, abs=1e-12)
    assert float(fields['rtol']) == rtol
    assert float(fields['atol']) == atol
    print('underresolved adjoint discrepancy:', discrepancy,
          'forward/adjoint:', forward, adjoint,
          'largest gradient:', float(fields['largest_absolute_gradient']))
