"""Independent manufactured checks for the polynomial moving-cap reference."""
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_polynomial_reference as model
from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary


CONTROLS = ControlBudget(maximum_evaluations=2000000,cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('polynomial cap controls:',CONTROLS.budget.evaluations,
          'coefficient/derivative values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def inputs(**changes):
    data = dict(Z=.5,chi_coefficients=np.array([.03,.02,-.005]),
                density_coefficients=np.array([.01,.015,-.002]),
                alpha=.002,T_x=.7,radius_rate=.04,incoming_w_plus=.4)
    data.update(changes)
    return data


def evaluate(budget,**changes):
    return model.polynomial_cap_reference(**inputs(**changes),budget=budget)


def test_zero_alpha_matches_independently_integrated_quadratic_functional(budget):
    data = inputs(alpha=0.)
    row = model.polynomial_cap_reference(**data,budget=budget)
    Z,c,k = data['Z'],data['chi_coefficients'],data['density_coefficients']
    xi = s.symbols('xi',real=True)
    powers = [1,3,5]
    basis = s.Matrix([xi**p for p in powers])
    budget.check(200)
    M = Z*np.array(s.integrate(basis*basis.T,(xi,0,1)),dtype=float)
    K = np.array(s.integrate(basis.diff(xi)*basis.diff(xi).T,(xi,0,1)),dtype=float)/Z
    R,A0 = np.tanh(Z),1/np.cosh(Z)**2
    e = np.ones(3)
    S = K-A0/R*np.outer(e,e)
    mass = (c@S@c+k@M@k)/4
    gc,gk = S@c/2,M@k/2
    shape = (k@M@k/Z-c@K@c/Z+(2*A0+(A0/R)**2)*sum(c)**2)/4
    assert row['mass'] == pytest.approx(mass,abs=2e-9)
    np.testing.assert_allclose(row['gradient']['chi'],gc,atol=2e-9,rtol=0)
    np.testing.assert_allclose(row['gradient']['canonical_density'],gk,atol=2e-9,rtol=0)
    assert row['gradient']['Z'] == pytest.approx(shape,abs=2e-9)
    expected = data['T_x']*k+data['radius_rate']/A0*np.array(powers)*c/Z
    np.testing.assert_allclose(row['rates']['chi'],expected,atol=2e-9,rtol=0)
    assert abs(row['wall']['momentum_discrepancy']) < 2e-9
    assert abs(row['balance_defect']) < 1e-10


def test_finite_alpha_shape_and_two_radial_tolerances(budget):
    coarse = evaluate(budget,rtol=2e-9,atol=2e-11)
    fine = evaluate(budget)
    assert abs(fine['shape_discrepancy']) < 2e-8
    assert fine['gradient_discrepancy'] < 2e-8
    assert abs(fine['balance_defect']) < 1e-10
    assert fine['mass'] == pytest.approx(coarse['mass'],abs=2e-8)
    for name in ['chi','canonical_density','Z']:
        np.testing.assert_allclose(fine['gradient'][name],coarse['gradient'][name],atol=2e-8,rtol=0)
    for name in ['chi','canonical_density']:
        np.testing.assert_allclose(fine['rates'][name],coarse['rates'][name],atol=2e-8,rtol=0)
    assert fine['wall']['normal_discrepancy'] != 0
    assert not fine['numerical_error_enclosure']
    print('manufactured polynomial:',{name:fine[name] for name in
          ['mass','shape_discrepancy','gradient_discrepancy','balance_defect','gram_condition','interface_balance_defect']},
          'metric floor:',fine['metric_certificate']['A_lower_float_down'],
          'wall discrepancies:',{name:fine['wall'][name] for name in
          ['momentum_discrepancy','normal_discrepancy','scalar_tangent_discrepancy']})


def test_joint_radius_and_coefficient_direction_is_independent(budget):
    data = inputs()
    dc,dk,dZ = np.array([.1,-.05,.02]),np.array([-.03,.07,.01]),.04
    epsilon = 1e-4
    row = model.polynomial_cap_reference(**data,budget=budget)
    plus = evaluate(budget,Z=data['Z']+epsilon*dZ,chi_coefficients=data['chi_coefficients']+epsilon*dc,
                    density_coefficients=data['density_coefficients']+epsilon*dk)
    minus = evaluate(budget,Z=data['Z']-epsilon*dZ,chi_coefficients=data['chi_coefficients']-epsilon*dc,
                     density_coefficients=data['density_coefficients']-epsilon*dk)
    difference = (plus['mass']-minus['mass'])/(2*epsilon)
    expected = row['gradient']['chi']@dc+row['gradient']['canonical_density']@dk+row['gradient']['Z']*dZ
    assert difference == pytest.approx(expected,abs=2e-8)


def test_inherited_clock_and_gram_solution_have_independent_checks(budget):
    data = inputs()
    row = model.polynomial_cap_reference(**data,budget=budget)
    budget.check(800)
    with mp.workdps(75):
        Z,alpha = mp.mpf(data['Z']),mp.mpf(data['alpha'])
        coefficients = [mp.mpf(v) for v in data['chi_coefficients']]
        expected_centre = -mp.log(mp.cosh(mp.sqrt(alpha/6)*coefficients[0]/Z))
        assert row['centre']['delta'] == pytest.approx(float(expected_centre),abs=1e-18)
        M = mp.matrix([[Z/(i+j+1) for j in [1,3,5]] for i in [1,3,5]])
        velocity = 2*mp.exp(mp.mpf(row['lapse']))*mp.lu_solve(M,mp.matrix(row['gradient']['canonical_density'].tolist()))
        transported = data['radius_rate']/(1/np.cosh(data['Z'])**2)*np.array([1,3,5])*data['chi_coefficients']/data['Z']
        np.testing.assert_allclose(row['rates']['chi'],
                                   data['T_x']*np.array(list(velocity),dtype=float)+transported,atol=1e-11,rtol=0)
    # Separate lapse integral evaluates the original polynomial expressions,
    # with no call into the field evaluator or radial solution.
    nodes,weights = np.polynomial.legendre.leggauss(40)
    xi = (nodes+1)/2
    z = data['Z']*xi
    r,A0 = np.tanh(z),1/np.cosh(z)**2
    c,k = data['chi_coefficients'],data['density_coefficients']
    chi = sum(cj*xi**p for cj,p in zip(c,[1,3,5]))
    q = sum(p*cj*xi**(p-1) for cj,p in zip(c,[1,3,5]))/(data['Z']*A0)
    p = sum(kj*xi**p for kj,p in zip(k,[1,3,5]))/A0
    a = data['alpha']*(p*p+(q-chi/r)**2)/(2*r)
    expected_wall = float(expected_centre)+data['Z']*(weights@(A0*a))/2
    assert row['lapse'] == pytest.approx(expected_wall,abs=1e-13)
    assert row['gram_normalized_residual'] < 64*np.finfo(float).eps


def test_global_polynomial_slope_certificate_is_outward(budget):
    data = inputs()
    precision = model.ctx.prec
    row = model.polynomial_cap_reference(**data,budget=budget)
    assert model.ctx.prec == precision
    budget.check(256)
    certificate = row['metric_certificate']
    bound = certificate['A_lower']
    dyadic = Fraction(int(bound['mantissa']))*Fraction(2)**bound['exponent']
    assert Fraction(certificate['A_lower_float_down']) <= dyadic
    with mp.workdps(80):
        Z,alpha = mp.mpf(data['Z']),mp.mpf(data['alpha'])
        R,A0 = mp.tanh(Z),1/mp.cosh(Z)**2
        lc = sum(p*abs(mp.mpf(v)) for p,v in zip([1,3,5],data['chi_coefficients']))/Z
        lk = sum(p*abs(mp.mpf(v)) for p,v in zip([1,3,5],data['density_coefficients']))/Z
        F = lc*Z/R
        argument = mp.sqrt(alpha/6)*F
        P = F*F*(mp.sinh(argument)/argument)**2*(mp.cosh(argument)**2+1)/2
        exact = A0-2*alpha*(R*R*P/2+((lk*Z/A0)**2+(lc/A0+F)**2)/4)
        endpoint = mp.mpf(dyadic.numerator)/dyadic.denominator
        assert 0 < endpoint <= exact
        assert exact-endpoint < mp.mpf('1e-40')
    assert row['sampled_minimum_A'] > float(dyadic)


@pytest.mark.parametrize('rx',[-.04,0.,.04])
def test_boundary_mass_and_scalar_tangent_keep_their_distinct_traces(rx,budget):
    row = evaluate(budget,radius_rate=rx)
    w = row['wall']
    flux = w['A']*w['B']*w['p_weak']*(w['q_weak']-w['f'])/2
    assert row['physical_weak_flux'] == pytest.approx(flux,abs=1e-14)
    assert row['mass_rate'] == pytest.approx(row['T_x']*flux+rx*w['rho'],abs=1e-10)
    assert w['weak_scalar_tangent_discrepancy'] == pytest.approx(rx*w['normal_discrepancy']/row['radius'],abs=1e-14)
    expected = ((rx-row['T_x']*w['B'])*w['momentum_discrepancy']
                +(rx+row['T_x']*w['B'])*w['normal_discrepancy'])/(2*row['radius'])
    assert w['scalar_tangent_discrepancy'] == pytest.approx(expected,abs=1e-14)
    direct_interface = (row['T_x']*w['B']*w['p_interface']+rx*(w['q_interface']-w['f']))/row['radius']
    assert w['transmitted_scalar_tangent'] == pytest.approx(direct_interface,abs=1e-14)
    assert row['interface_balance_defect'] == pytest.approx(
        row['balance_defect']+row['interface_flux_trace_gap']+row['interface_moving_density_gap'],abs=1e-14)
    assert w['incoming_discrepancy'] == pytest.approx((w['momentum_discrepancy']+w['normal_discrepancy'])/row['radius'])


def test_shape_guard_rejects_jointly_corrupted_forward_and_adjoint_values(budget,monkeypatch):
    original = model.radial._radial_variations
    def corrupted(*args,**kwargs):
        row = original(*args,**kwargs)
        row['gradient'][-1] += .1
        row['adjoint_gradient'][-1] += .1
        return row
    monkeypatch.setattr(model.radial,'_radial_variations',corrupted)
    def forbidden(*args,**kwargs):
        raise AssertionError('known bad shape reached a weak Gram solve')
    monkeypatch.setattr(model.np.linalg,'solve',forbidden)
    with pytest.raises(ValueError,match='shape and dilation'):
        evaluate(budget)


def test_actual_hybrid_boundary_recovers_reported_transmission(budget):
    data = inputs()
    row = model.polynomial_cap_reference(**data,budget=budget)
    w,R,alpha = row['wall'],row['radius'],data['alpha']
    tx,rx = data['T_x'],data['radius_rate']
    budget.check(64)
    proper_rate = np.sqrt((tx*w['B'])**2-rx*rx)/np.sqrt(w['A'])
    H = (tx*w['B']+rx)/proper_rate
    qf = np.sqrt(alpha/6)*w['f']
    u = w['f']*np.sinh(qf)/qf
    D = 1+alpha*u*u/6
    tension = 2*D*np.sqrt(D)
    # Construct one admissible local exterior-jet input, not a source history.
    boundary = hybrid_wall_boundary(
        dict(r=R,c=R*u,eta_plus=0.),
        dict(r_V=(H-tension*R)/2,c_V=H*(u+R*np.sqrt(D)*data['incoming_w_plus'])/2,log_C_V=0.),
        dict(A=w['A'],delta=row['lapse'],w_minus=w['outgoing_w_minus']),
        dict(C_wall=1.,log_C_s=0.),alpha=alpha)
    assert boundary['T_s']*proper_rate == pytest.approx(tx,abs=1e-13)
    assert boundary['r_s']*proper_rate == pytest.approx(rx,abs=1e-13)
    assert boundary['w_plus'] == pytest.approx(data['incoming_w_plus'],abs=1e-13)
    assert boundary['Pi']*R == pytest.approx(w['p_interface'],abs=1e-13)
    assert boundary['Phi']*R+w['f'] == pytest.approx(w['q_interface'],abs=1e-13)
    assert boundary['f_s']*proper_rate == pytest.approx(w['transmitted_scalar_tangent'],abs=1e-13)


def test_single_odd_mode_agrees_with_the_existing_linear_cap(budget):
    data = inputs(chi_coefficients=np.array([.03]),density_coefficients=np.array([.01]),radius_rate=0.)
    new = model.polynomial_cap_reference(**data,budget=budget)
    old = model.radial.variational_cap_reference(
        [0.,data['Z']],[0.,.03],[0.,.01],alpha=data['alpha'],T_x=data['T_x'],
        incoming_w_plus=data['incoming_w_plus'],budget=budget)
    for name in ['mass','lapse','mass_rate','physical_weak_flux']:
        assert new[name] == pytest.approx(old[name],abs=2e-9)
    for name in ['chi','canonical_density']:
        np.testing.assert_allclose(new['gradient'][name],old['gradient'][name],atol=2e-9,rtol=0)
        np.testing.assert_allclose(new['rates'][name],old['rates'][name][1:],atol=2e-9,rtol=0)


def test_quiet_data_and_input_immutability(budget):
    data = inputs(chi_coefficients=np.zeros(3),density_coefficients=np.zeros(3),incoming_w_plus=0.)
    row = model.polynomial_cap_reference(**data,budget=budget)
    assert row['mass'] == row['lapse'] == row['mass_rate'] == row['gradient']['Z'] == 0
    for name in ['chi','canonical_density']:
        np.testing.assert_array_equal(row['rates'][name],np.zeros(3))
    for name in ['chi_coefficients','density_coefficients']:
        np.testing.assert_array_equal(data[name],np.zeros(3))
        assert not np.shares_memory(row[name],data[name])


def test_shape_normal_derivative_gap_identity(budget):
    budget.check(16)
    Z,xi,r,A0,f,h = s.symbols('Z xi r A0 f h',nonzero=True)
    gap = 1/(Z*xi)-A0/r
    direct = (-1/Z+2*r*xi)*(h+f)+f*xi*A0/r
    stable = (-1/Z+2*r*xi)*h+xi*(2*r-gap)*f
    assert s.expand(direct-stable) == 0


def test_actual_interface_scalar_pair_has_both_trace_errors(budget):
    budget.check(16)
    p,h,dp,dh,R,B,tx,rx = s.symbols('p h dp dh R B tx rx',nonzero=True)
    # p,h denote the weak pair; trace values are p+dp,h+dh.
    incoming,outgoing = (p+h)/R,(p+dp-h-dh)/R
    pi,hi = R*(incoming+outgoing)/2,R*(incoming-outgoing)/2
    computed = (tx*B*p+rx*(h+dh))/R
    transmitted = (tx*B*pi+rx*hi)/R
    expected = ((rx-tx*B)*dp+(rx+tx*B)*dh)/(2*R)
    assert s.expand(computed-transmitted-expected) == 0


@pytest.mark.parametrize('changes',[
    dict(Z=0.),dict(alpha=-.1),dict(alpha=10000.),dict(T_x=0.),dict(rtol=1e-20),
    dict(chi_coefficients=[.1,.2,.3,.4]),dict(density_coefficients=[.1]),
    dict(chi_coefficients=[np.nan,.02,-.005]),
])
def test_invalid_inputs_reject_before_radial_solve(changes,budget,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('invalid input started a radial solve')
    monkeypatch.setattr(model.radial,'_radial_variations',forbidden)
    with pytest.raises(ValueError):
        evaluate(budget,**changes)


def test_nontimelike_wall_is_rejected(budget):
    with pytest.raises(ValueError,match='not timelike'):
        evaluate(budget,radius_rate=1.)


def test_shared_budget_stops_before_radial_solve(budget,monkeypatch):
    maximum = budget.maximum_evaluations
    budget.maximum_evaluations = budget.evaluations+1
    def forbidden(*args,**kwargs):
        raise AssertionError('unfunded reference started radial integration')
    monkeypatch.setattr(model.radial,'_radial_variations',forbidden)
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            evaluate(budget)
    finally:
        budget.maximum_evaluations = maximum
