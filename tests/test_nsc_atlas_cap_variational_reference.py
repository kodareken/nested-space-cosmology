"""Manufactured finite-backreaction gradients, not physical cap histories."""
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_variational_reference as model
from recursive_horizons.nsc_atlas_cap_geometric_reference import geometric_cap_reference

CONTROLS = ControlBudget(maximum_evaluations=2000000, cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module', autouse=True)
def report():
    yield
    print('variational cap controls:', CONTROLS.budget.evaluations,
          'coefficient/derivative values;', CONTROLS.active_cpu_seconds,
          'active CPU seconds')


def inputs(**changed):
    result = dict(z_nodes=np.array([0., .07, .21, .5]),
                  chi=np.array([0., .03, .12, .13]),
                  canonical_density=np.array([0., -.02, .07, .08]),
                  alpha=.002,T_x=.7,incoming_w_plus=.4)
    result.update(changed)
    return result


def evaluate(budget, **changed):
    return model.variational_cap_reference(**inputs(**changed),budget=budget)


def test_zero_alpha_recovers_independent_exact_element_mass_and_gradients(budget):
    data = inputs(alpha=0.)
    row = model.variational_cap_reference(**data,budget=budget)
    old = geometric_cap_reference(data['z_nodes'],data['chi'],data['canonical_density'],
                                  T_x=data['T_x'],incoming_w_plus=data['incoming_w_plus'],budget=budget)
    assert row['mass'] == pytest.approx(old['mass'],abs=2e-9)
    for newkey,oldkey in [('chi','chi'),('canonical_density','velocity')]:
        np.testing.assert_allclose(row['gradient'][newkey],old['gradient'][oldkey],rtol=0,atol=2e-9)
        np.testing.assert_allclose(row['rates'][newkey],old['rates'][oldkey],rtol=0,atol=2e-9)
    assert row['lapse'] == 0
    assert row['gradient_discrepancy'] < 2e-9
    assert abs(row['wall']['momentum_discrepancy']) < 2e-9


def test_nonzero_alpha_forward_adjoint_and_tolerance_checks(budget):
    coarse = evaluate(budget,rtol=2e-9,atol=2e-11)
    fine = evaluate(budget,rtol=2e-11,atol=2e-13)
    assert coarse['gradient_discrepancy'] < 2e-8
    assert fine['gradient_discrepancy'] < 2e-8
    assert abs(fine['mass']-coarse['mass']) < 2e-8
    for name in ['chi','canonical_density']:
        np.testing.assert_allclose(fine['gradient'][name],coarse['gradient'][name],rtol=0,atol=2e-8)
        np.testing.assert_allclose(fine['rates'][name],coarse['rates'][name],rtol=0,atol=2e-8)
    assert fine['lapse'] > fine['centre']['delta']
    assert abs(fine['balance_defect']) < 1e-10
    assert not fine['numerical_error_enclosure']
    assert abs(fine['wall']['momentum_discrepancy']) > 1e-8
    assert abs(fine['wall']['normal_discrepancy']) > .01
    print('finite-alpha manufactured gradient discrepancy:',fine['gradient_discrepancy'],
          'balance defect:',fine['balance_defect'],
          'metric lower bound:',fine['metric_certificate']['A_lower_float_down'],
          'momentum discrepancy:',fine['wall']['momentum_discrepancy'],
          'normal discrepancy:',fine['wall']['normal_discrepancy'])


def test_inherited_centre_clock_and_lapse_integral_fix_the_rate_normalization(budget):
    data = inputs()
    row = model.variational_cap_reference(**data,budget=budget)
    z, chi, k = data['z_nodes'],data['chi'],data['canonical_density']
    with mp.workdps(65):
        expected_centre = float(-mp.log(mp.cosh(mp.sqrt(mp.mpf(data['alpha'])/6)
                                               *mp.mpf(chi[1])/mp.mpf(z[1]))))
    assert expected_centre < 0
    assert row['centre']['delta'] == pytest.approx(expected_centre,abs=1e-18)
    nodes,weights = np.polynomial.legendre.leggauss(32)
    increment = 0.
    M = np.zeros((len(z),len(z)))
    for j,(lo,hi) in enumerate(zip(z[:-1],z[1:])):
        budget.check(9*len(nodes)+4)
        points = (lo+hi)/2+(hi-lo)*nodes/2
        r, A0 = np.tanh(points),1/np.cosh(points)**2
        basis = (points-lo)/(hi-lo)
        field = (1-basis)*chi[j]+basis*chi[j+1]
        density = (1-basis)*k[j]+basis*k[j+1]
        p = density/A0
        h = (chi[j+1]-chi[j])/((hi-lo)*A0)-field/r
        integrand = data['alpha']*A0*(p*p+h*h)/(2*r)
        increment += (hi-lo)*(weights@integrand)/2
        M[j:j+2,j:j+2] += (hi-lo)*np.array([[2.,1.],[1.,2.]])/6
    expected_wall = expected_centre+increment
    assert row['lapse'] == pytest.approx(expected_wall,abs=1e-13)
    expected_rate = 2*data['T_x']*np.exp(expected_wall)*np.linalg.solve(
        M[1:,1:],row['gradient']['canonical_density'])
    np.testing.assert_allclose(row['rates']['chi'][1:],expected_rate,rtol=0,atol=1e-12)
    assert row['wall']['B'] == pytest.approx(row['wall']['A']*np.exp(expected_wall),abs=1e-13)


def test_directional_mass_change_is_independent_of_the_gradient_solve(budget):
    base = inputs()
    dc, dk = np.array([0., .1, -.15, .2]), np.array([0., -.2, .05, .1])
    epsilon = 2e-4
    row = model.variational_cap_reference(**base,budget=budget)
    plus = evaluate(budget,chi=base['chi']+epsilon*dc,canonical_density=base['canonical_density']+epsilon*dk)
    minus = evaluate(budget,chi=base['chi']-epsilon*dc,canonical_density=base['canonical_density']-epsilon*dk)
    difference = (plus['mass']-minus['mass'])/(2*epsilon)
    predicted = row['gradient']['chi']@dc[1:]+row['gradient']['canonical_density']@dk[1:]
    assert difference == pytest.approx(predicted,abs=2e-8)


def test_exact_binary_metric_certificate_and_incoming_boundary_flux(budget):
    precision = model.ctx.prec
    row = evaluate(budget,incoming_w_plus=0.)
    assert model.ctx.prec == precision
    certificate = row['metric_certificate']
    exact = certificate['A_lower']
    bound = Fraction(int(exact['mantissa']))*Fraction(2)**exact['exponent']
    assert Fraction(certificate['A_lower_float_down']) <= bound
    assert bound > 0
    assert row['sampled_minimum_A'] > float(bound)
    assert row['wall']['A'] > float(bound)
    wall = row['wall']
    direct = wall['A']*wall['B']*wall['p_weak']*(wall['q_weak']-wall['f'])/2
    assert direct == pytest.approx(row['physical_weak_flux'],abs=1e-14)
    assert row['mass_rate'] == pytest.approx(row['T_x']*direct,abs=1e-10)
    assert row['mass_rate'] <= 0
    assert row['incoming_flux_bound'] == 0
    assert wall['incoming_discrepancy'] == pytest.approx(
        (wall['momentum_discrepancy']+wall['normal_discrepancy'])/row['radius'])


def test_independent_adjoint_rejects_corrupted_forward_gradient_before_weak_rhs(budget,monkeypatch):
    original = model.solve_ivp
    def corrupted(*args,**kwargs):
        result = original(*args,**kwargs)
        if args[1][1] == .5:
            result.y[2:,-1] += .25
        return result
    def forbidden(*args,**kwargs):
        raise AssertionError('known bad gradient reached the weak RHS')
    monkeypatch.setattr(model,'solve_ivp',corrupted)
    monkeypatch.setattr(model,'solve_banded',forbidden)
    with pytest.raises(ValueError,match='gradients exceed'):
        evaluate(budget)


def test_exact_skew_pairing_and_geometric_field_identity(budget):
    budget.check(16)
    c, v, gc, gk, d = s.symbols('c v gc gk d',nonzero=True)
    mass_matrix = s.symbols('M',positive=True)
    V = 2*c*gk/mass_matrix
    k_rate = (-2*c*gc+d-V)/mass_matrix
    assert s.expand(gc*V+gk*k_rate-V*(d-V)/(2*c)) == 0
    z, r, A0, a, b = s.symbols('z r A0 a b',nonzero=True)
    f, q = (a*z+b)/r, a/A0
    gap = 1/z-A0/r
    assert s.expand(q-f-(a*z*gap/A0-b/r)) == 0


def test_quiet_nonzero_alpha_and_input_immutability(budget):
    data = inputs(chi=np.zeros(4),canonical_density=np.zeros(4),incoming_w_plus=0.)
    before = {key:value.copy() for key,value in data.items() if isinstance(value,np.ndarray)}
    row = model.variational_cap_reference(**data,budget=budget)
    assert row['mass'] == row['lapse'] == row['mass_rate'] == 0
    for key in ['chi','canonical_density']:
        np.testing.assert_array_equal(row['gradient'][key],np.zeros(3))
        np.testing.assert_array_equal(row['rates'][key],np.zeros(4))
    for key,old in before.items():
        np.testing.assert_array_equal(data[key],old)
        assert not np.shares_memory(row[key],data[key])


@pytest.mark.parametrize('changed',[
    dict(alpha=-.1),dict(alpha=100.),dict(alpha=True),dict(T_x=0.),
    dict(z_nodes=[0.,.07,.07,.5]),dict(z_nodes=[0.,.07,.21,100.]),
    dict(chi=[.01,.03,.12,.13]),dict(canonical_density=[0.,np.nan,.07,.08]),
    dict(chi=[0.,1.]),dict(rtol=1e-20),
])
def test_rejected_input_or_sufficient_bound_does_not_start_radial_integration(changed,budget,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('rejected data started an ODE')
    monkeypatch.setattr(model,'solve_ivp',forbidden)
    with pytest.raises(ValueError):
        evaluate(budget,**changed)


def test_shared_budget_prevents_an_unfunded_radial_solve(budget,monkeypatch):
    original = budget.maximum_evaluations
    budget.maximum_evaluations = budget.evaluations+1
    def forbidden(*args,**kwargs):
        raise AssertionError('exhausted budget started an ODE')
    monkeypatch.setattr(model,'solve_ivp',forbidden)
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            evaluate(budget)
    finally:
        budget.maximum_evaluations = original
