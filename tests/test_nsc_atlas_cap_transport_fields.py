"""Manufactured fixed-radius/material jets; no retained cap or radial solve."""
import copy

import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_background_cap import _total as native_fields
from recursive_horizons.nsc_atlas_reduced_cap import _reconstruct
from recursive_horizons.nsc_atlas_cap_transport_fields import CapTransportFields

CONTROLS = ControlBudget(maximum_evaluations=200000, cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("transport field controls:", CONTROLS.budget.evaluations,
          "charged values;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def profile(zero=True):
    z = np.array([-1.,0.,1.]); c=0. if zero else 1.
    return CharacteristicProfile(z,c*z**3,3*c*z*z,6*c*z)


def context(budget, *, R=.4, T=.1, alpha=0., state=None, source=None):
    """Only the frozen field inputs that the adapter reads; no metric solve."""
    y=np.linspace(0,1,5)
    return dict(y_faces=y,radius=R,r_faces=R*y,T=T,alpha=alpha,
        state=state or dict(Ip=np.zeros(4),Iq=np.zeros(4)),
        profile=source or profile(),budget=budget)


def direction(*, R=0., T=0., Ip=None, Iq=None):
    return dict(R=R,T=T,Ip=np.zeros(4) if Ip is None else np.asarray(Ip),
                Iq=np.zeros(4) if Iq is None else np.asarray(Iq))


def backward_values(base, rates, cell, radii, h, budget, *, material=False):
    result=[]
    for step in (0.,-h,-2*h):
        R=base['radius']+step*rates['R']; T=base['T']+step*rates['T']
        state={k:base['state'][k]+step*rates[k] for k in ('Ip','Iq')}
        rec=_reconstruct(R*base['y_faces'],state)
        at=radii*R/base['radius'] if material else radii
        result.append(native_fields(base['profile'],T,rec,cell,at,budget))
    return {k:(3*result[0][k]-4*result[1][k]+result[2][k])/(2*h) for k in ('f','Pi','Phi')}


@pytest.mark.parametrize('alpha,q',[(0.,.2),(.02,.2),(.02,2e-4)])
def test_regular_polynomial_oracle_with_moving_cells_and_both_potential_branches(alpha,q,budget):
    R,p,rx,pdot,qdot=.4,.3,.07,.11,-.04
    faces=R*np.linspace(0,1,5)
    state=dict(Ip=p*np.diff(faces*faces)/2,Iq=q*np.diff(faces))
    rates=direction(R=rx,T=.2,Ip=(pdot/p+2*rx/R)*state['Ip'],Iq=(qdot/q+rx/R)*state['Iq'])
    view=CapTransportFields(context(budget,state=state,alpha=alpha),rates,budget=budget)
    z=np.sqrt(alpha/6)*q
    if abs(z)<1e-4:
        sinhc=1+z*z/6+z**4/120
        U=q*sinhc*((np.cosh(z)**2+1)*(1+z*z/2+z**4/24)+z*sinhc*np.cosh(z)*np.sinh(z))
    else:
        sinhc=np.sinh(z)/z
        U=2*q*sinhc*np.cosh(z)**3
    potential=q*q*sinhc*sinhc*(np.cosh(z)**2+1)/2
    for cell,r in ((0,.003),(2,.25),(3,.35)):
        row=view.sample(cell,[r],budget=budget)
        for k,v in dict(f=q,Pi=p,Phi=0.,f_x=qdot,Pi_x=pdot,Phi_x=0.,f_r=0.,Pi_r=0.,Phi_r=0.).items():
            np.testing.assert_allclose(row[k],v,rtol=2e-12,atol=4e-13)
        expected=dict(a=alpha*r*p*p/2,a_x=alpha*r*p*pdot,a_r=alpha*p*p/2,
            b=r*r*potential/2+r*r*(1-r*r)*p*p/4,
            b_x=r*r*U*qdot/2+r*r*(1-r*r)*p*pdot/2,
            b_r=r*potential+(r/2-r**3)*p*p)
        for key,value in expected.items():
            np.testing.assert_allclose(row[key],value,rtol=2e-12,atol=4e-13)


@pytest.mark.parametrize('cell,r',[(0,.003),(2,.27)])
def test_full_direction_matches_backward_native_fields(cell,r,budget):
    state=dict(Ip=np.array([.002,.005,.013,.025]),Iq=np.array([.003,.001,-.0005,-.002]))
    base=context(budget,state=state,source=profile(False),alpha=.02)
    rates=direction(R=.02,T=.3,Ip=[.01,.03,.02,.04],Iq=[-.01,.02,.04,.01])
    row=CapTransportFields(base,rates,budget=budget).sample(cell,[r],budget=budget)
    oracle=backward_values(base,rates,cell,np.array([r]),1e-5,budget)
    for key in ('f','Pi','Phi'):
        np.testing.assert_allclose(row[key+'_x'],oracle[key],rtol=3e-6,atol=3e-7)


@pytest.mark.parametrize('cell,r',[(0,.003),(2,.27)])
def test_radial_column_matches_nonconstant_native_fields(cell,r,budget):
    state=dict(Ip=np.array([.002,.005,.013,.025]),Iq=np.array([.003,.001,-.0005,-.002]))
    base=context(budget,state=state,source=profile(False),alpha=.02)
    row=CapTransportFields(base,direction(),budget=budget).sample(cell,[r],budget=budget)
    rec=_reconstruct(base['r_faces'],state); h=1e-6
    minus=native_fields(base['profile'],base['T'],rec,cell,np.array([r-h]),budget)
    plus=native_fields(base['profile'],base['T'],rec,cell,np.array([r+h]),budget)
    for key in ('f','Pi','Phi'):
        np.testing.assert_allclose(row[key+'_r'],(plus[key]-minus[key])/(2*h),rtol=2e-7,atol=2e-8)


def test_upstream_primitive_and_first_cell_moment_are_differentiated(budget):
    base=context(budget)
    row=CapTransportFields(base,direction(Iq=[1.,0.,0.,0.]),budget=budget).sample(3,[.35],budget=budget)
    np.testing.assert_allclose(row['f_x'],1/.35,rtol=1e-13)
    np.testing.assert_allclose(row['Phi_x'],-1/.35**2,rtol=1e-13)
    np.testing.assert_array_equal(row['Pi_x'],0.)
    row=CapTransportFields(base,direction(Iq=[0.,1.,0.,0.]),budget=budget).sample(0,[.05],budget=budget)
    np.testing.assert_allclose(row['f_x'],-1.25,rtol=1e-13)
    np.testing.assert_allclose(row['Phi_x'],50/3,rtol=1e-13)


def test_zero_limiter_opens_on_the_declared_backward_direction(budget):
    base=context(budget)
    rates=direction(Ip=np.diff(base['r_faces']**3)/3)
    view=CapTransportFields(base,rates,budget=budget)
    assert any(v!='zero' for v in view.limiter_branches['plus'][1:-1])
    row=view.sample(2,[.27],budget=budget)
    oracle=backward_values(base,rates,2,np.array([.27]),1e-5,budget)
    np.testing.assert_allclose(row['Pi_x'],oracle['Pi'],rtol=1e-12,atol=1e-13)
    quiet=CapTransportFields(base,direction(),budget=budget).sample(2,[.27],budget=budget)
    for key in ('f_x','Pi_x','Phi_x','a_x','b_x'):
        np.testing.assert_array_equal(quiet[key],0.)


def test_both_face_traces_keep_the_material_chain_rule(budget):
    state=dict(Ip=np.array([.002,.005,.013,.025]),Iq=np.array([.003,.001,-.0005,-.002]))
    base=context(budget,state=state,source=profile(False))
    rates=direction(R=.02,T=.3,Ip=[.01,.03,.02,.04],Iq=[-.01,.02,.04,.01])
    view=CapTransportFields(base,rates,budget=budget)
    radius=base['r_faces'][2]; speed=rates['R']*base['y_faces'][2]
    sides=[]
    for cell in (1,2):
        row=view.sample(cell,[radius],budget=budget); sides.append(row)
        oracle=backward_values(base,rates,cell,np.array([radius]),1e-5,budget,material=True)
        for key in ('f','Pi','Phi'):
            material=row[key+'_x']+speed*row[key+'_r']
            np.testing.assert_allclose(material,oracle[key],rtol=3e-6,atol=3e-7)
        assert abs(row['Pi_x'][0]-oracle['Pi'][0])>1e-5
    assert sides[0]['Pi'][0] != sides[1]['Pi'][0]


def test_frozen_inputs_and_returned_arrays_do_not_modify_the_adapter(budget):
    base=context(budget,source=profile(False)); rates=direction(T=.2,Ip=[1.,2.,3.,4.])
    before=copy.deepcopy(base['state'])
    view=CapTransportFields(base,rates,budget=budget)
    original=view.sample(1,[.15],budget=budget)
    for k in before: np.testing.assert_array_equal(base['state'][k],before[k])
    base['state']['Iq'][:]=17.; rates['Ip'][:]=9.
    base['profile']._spline.c[:]=0.
    changed=view.sample(1,[.15],budget=budget)
    for key in original: np.testing.assert_array_equal(changed[key],original[key])
    changed['f'][:]=42.
    np.testing.assert_array_equal(view.sample(1,[.15],budget=budget)['f'],original['f'])


@pytest.mark.parametrize('cell,r',[(True,[.05]),(-1,[.05]),(4,[.05]),(0,[0.]),
    (0,[.11]),(1,[.05]),(0,[np.nan]),(0,[]),(0,[.05]*514),(0,[.01])])
def test_invalid_queries_fail_before_any_unadmitted_field_use(cell,r,budget):
    view=CapTransportFields(context(budget),direction(),budget=budget)
    with pytest.raises(ValueError): view.sample(cell,r,budget=budget)


def test_potential_switch_and_leaving_profile_directions_are_rejected(budget):
    faces=.5*np.linspace(0,1,5)
    base=context(budget,R=.5,alpha=6.,state=dict(Ip=np.zeros(4),Iq=1e-4*np.diff(faces)))
    view=CapTransportFields(base,direction(),budget=budget)
    with pytest.raises(ValueError,match='potential formula switch'):
        view.sample(0,[.05],budget=budget)
    base=context(budget,T=1-np.arctanh(.4))
    with pytest.raises(ValueError,match='upper profile'):
        CapTransportFields(base,direction(T=-1.),budget=budget)
    # Material upper-label motion is inward, but the fixed-r partial would
    # still need an unavailable extension. Do not silently extrapolate it.
    view=CapTransportFields(base,direction(R=1.,T=-1.),budget=budget)
    with pytest.raises(ValueError,match='fixed-radius backward'):
        view.sample(3,[.4],budget=budget)
    allowed=CapTransportFields(base,direction(T=1.),budget=budget)
    assert np.isfinite(allowed.sample(3,[.4],budget=budget)['f_x']).all()


def test_same_budget_and_exhaustion_are_enforced(budget):
    base=context(budget); view=CapTransportFields(base,direction(),budget=budget)
    with pytest.raises(ValueError,match='same caller budget'):
        view.sample(0,[.05],budget=Budget())
    previous=budget.maximum_evaluations
    budget.maximum_evaluations=budget.evaluations+1
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            view.sample(0,[.05],budget=budget)
    finally:
        budget.maximum_evaluations=previous
