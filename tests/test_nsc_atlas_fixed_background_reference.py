"""Actual constrained time derivative for retained history and fixed departures."""
from fractions import Fraction as Q

import numpy as np
from numpy.polynomial import Polynomial
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_transfer import _transfer_coefficients
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
from recursive_horizons.nsc_atlas_cap_fixed_hybrid import fixed_cap_hybrid_reference
from recursive_horizons.nsc_atlas_cap_fixed_bound import fixed_background_metric_bound

CONTROLS=ControlBudget(maximum_evaluations=6000000,cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('fixed background controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def history(poly,knots=(-1.,0.,.5,1.)):
    knots=np.array(knots)
    return RetainedCapBackground(knots,poly(knots),poly.deriv()(knots),poly.deriv(2)(knots))


def test_same_total_polynomial_and_independent_time_mass_difference(budget):
    basis=FixedCapBasis([0,1],budget=budget);Z=Q(1,2);T=Q(1,8);a=Q(1,64)
    bg=history(Polynomial([0,0,0,float(a)]))
    def coefficients(power):
        return np.array(_transfer_coefficients(basis,[tuple(power)],budget),dtype=float)
    cb=coefficients([0,6*a*T*T,0,2*a,0,0])
    kb=coefficients([0,12*a*T,0,0,0,0])
    bb=coefficients([0,12*a,0,0,0,0])
    c=np.array([.001,-.002,.003]);d=np.array([.003,.001,-.002])
    opts=dict(alpha=.002,T_x=.7,radius_rate=.001,incoming_w_plus=.2,budget=budget)
    affine=fixed_cap_reference(basis,c,d,Z,background=bg,background_time=float(T),**opts)
    total=fixed_cap_reference(basis,c+cb,d+kb,Z,**opts)
    for key in ('mass','lapse','mass_rate','physical_weak_flux','shape_derivative_Z'):
        assert affine[key]==pytest.approx(total[key],abs=2e-10)
    for key in ('chi','canonical_density'):
        np.testing.assert_allclose(affine['gradient'][key],total['gradient'][key],rtol=0.,atol=2e-10)
    predicted=total['gradient']['chi']@kb+total['gradient']['canonical_density']@bb
    assert affine['gradient']['T']==pytest.approx(predicted,abs=2e-10)
    assert affine['adjoint_gradient']['T']==pytest.approx(predicted,abs=2e-10)
    for key,rate in (('chi',kb),('canonical_density',bb)):
        np.testing.assert_allclose(affine['rates'][key]+opts['T_x']*rate,total['rates'][key],rtol=0.,atol=2e-9)
    assert abs(affine['explicit_time_remainder'])<2e-10
    assert abs(affine['remainder_identity_defect'])<2e-12
    h=1e-4
    plus=fixed_cap_reference(basis,c,d,Z,background=bg,background_time=float(T)+h,**opts)
    minus=fixed_cap_reference(basis,c,d,Z,background=bg,background_time=float(T)-h,**opts)
    assert (plus['mass']-minus['mass'])/(2*h)==pytest.approx(affine['gradient']['T'],abs=2e-8)


def test_incoming_quartic_front_carries_exact_mass_and_hybrid_trace(budget):
    basis=FixedCapBasis([0,1],budget=budget);zero=np.zeros(basis.dimension)
    T=.125;Z=Q(1,2);a=1/64;front=.5;d=T+float(Z)-front
    bg=history(a*Polynomial([-front,1])**4,knots=(front,.75,1.))
    R=np.tanh(float(Z));A0=1-R*R;b=A0/R
    F=a*d**4;k=4*a*d**3;H=-.22+2*R
    def exterior(common):
        assert common['c']==pytest.approx(F,abs=2e-16)
        return dict(r_V=-.11,c_V=k*H/A0,log_C_V=0.)
    out=fixed_cap_hybrid_reference(basis,zero,zero,Z,alpha=0.,exterior_wall=exterior,
        C_wall=1.,eta_plus=0.,eta_plus_x=0.,background=bg,background_time=T,budget=budget)
    cap=out['cap'];N=8*a*a*d**7/7-b*F*F/4;NT=k*(k-b*F)/2
    assert cap['mass']==pytest.approx(N,abs=2e-14)
    assert cap['gradient']['T']==pytest.approx(NT,abs=2e-13)
    assert cap['metric_certificate']['interval_count']==2
    assert cap['background']['radial_intervals']==3
    assert abs(cap['explicit_time_remainder'])<2e-13
    assert abs(cap['remainder_identity_defect'])<2e-13
    assert abs(out['residual']['conformal_wall_flux'])<2e-12
    assert abs(out['residual']['temporal_mass_actual'])<2e-12
    assert cap['mass']<=cap['metric_certificate']['mass_upper']
    gram=basis.gram(Z,budget=budget)
    for rate in cap['rates'].values():
        assert np.sqrt(gram.squared_norm(rate,budget=budget))<2e-10
        assert abs(gram.wall['value']@rate[gram.wall['indices']])<2e-10
    print('incoming mass/time/flux:',cap['mass'],cap['gradient']['T'],
          out['residual']['conformal_wall_flux'],out['residual']['temporal_mass_actual'])


def test_nonrepresentable_history_retains_actual_finite_remainder(budget):
    basis=FixedCapBasis([0,1],budget=budget);zero=np.zeros(basis.dimension)
    bg=history(Polynomial([-.5,1])**4,knots=(.5,.75,1.))
    row=fixed_cap_reference(basis,zero,zero,Q(1,2),alpha=.002,T_x=1.,radius_rate=.001,
        incoming_w_plus=.03,background=bg,background_time=.125,budget=budget)
    assert abs(row['remainder_identity_defect'])<2e-12
    assert abs(row['explicit_time_remainder'])>1e-10
    assert row['gradient']['T']==pytest.approx(row['adjoint_gradient']['T'],abs=2e-10)
    assert row['sampled_minimum_A']>=row['metric_certificate']['A_lower_float_down']
    assert row['mass']<=row['metric_certificate']['mass_upper']


def test_quiet_history_and_contract_failures(budget):
    basis=FixedCapBasis([0,1],budget=budget);zero=np.zeros(basis.dimension)
    bg=history(Polynomial([-.5,1])**4,knots=(.5,.75,1.))
    opts=dict(alpha=.002,T_x=1.,radius_rate=0.,incoming_w_plus=0.,budget=budget)
    quiet=fixed_cap_reference(basis,zero,zero,Q(1,2),background=bg,background_time=-1.,**opts)
    assert quiet['mass']==quiet['lapse']==quiet['gradient']['T']==quiet['explicit_time_remainder']==0.
    linear=history(Polynomial([0.,1/64]))
    cancelling=np.array(_transfer_coefficients(basis,[(0,Q(-1,32),0,0,0,0)],budget),dtype=float)
    # Separate rounded representations do not exactly cancel. A mass below
    # the radial solver's numerical resolution must not be declared zero.
    try:
        cancelled=fixed_cap_reference(basis,cancelling,zero,Q(1,2),background=linear,background_time=0.,**opts)
    except ValueError as error:
        assert 'mass is unresolved' in str(error)
    else:
        assert 0.<=cancelled['mass']<1e-28
        assert abs(cancelled['lapse'])<1e-28
    constant=history(Polynomial([1.]))
    cancelled=fixed_cap_reference(basis,zero,zero,Q(1,2),background=constant,background_time=0.,**opts)
    assert cancelled['mass']==cancelled['lapse']==0.
    for update in (dict(background_time=.1),dict(background=object(),background_time=.1),
                   dict(background=bg,background_time=1.)):
        with pytest.raises(ValueError): fixed_cap_reference(basis,zero,zero,Q(1,2),**opts,**update)
    with pytest.raises(ValueError,match='not positive|potential bound'):
        fixed_background_metric_bound(basis,np.ones(3)*1e6,zero,Q(1,2),background=bg,T=.125,alpha=.002,budget=budget)
    tiny=history(Polynomial([0.,1e-200]))
    with pytest.raises(ValueError,match='mass is unresolved'):
        fixed_cap_reference(basis,zero,zero,Q(1,2),background=tiny,background_time=0.,**opts)
