"""Numerical geometry controls only: no physical profile, duration or C0 claim."""
from math import factorial

import numpy as np
import pytest

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleRadiusDirection,CompatibleIncomingMetric,INCOMING_SLOT_DERIVATIVES,
)
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_pg_ks_metric_pullback import reference_chart,ks_to_pg,pg_to_ks


def polynomial(coefficients):
    # Coherent analytic control; coefficients are ordinary monomial values.
    def derivative(z,order):
        return sum(c*factorial(n)/factorial(n-order)*z**(n-order)
                   for n,c in enumerate(coefficients) if n>=order)
    return derivative


def wave(k,phase):
    return lambda z,order: float(np.real((1j*k)**order*np.exp(1j*(k*z+phase))))


def direction():
    return CompatibleRadiusDirection(polynomial((.3,-.2,.07,.03)),polynomial((-.4,.09)),.15,.3)


def derivative_weights(order):
    nodes=np.arange(-2,3,dtype=float)
    moments=np.array([factorial(order) if n==order else 0 for n in range(5)],dtype=float)
    return nodes,np.linalg.solve(np.vstack([nodes**n for n in range(5)]),moments)


def test_six_incoming_derivatives_and_factorials_in_exact_plateau():
    control=direction();z=.27;step=1/32;slots=control.incoming_slots(z)
    assert set(slots)==set(INCOMING_SLOT_DERIVATIVES)
    for key,(_,dt,dz) in INCOMING_SLOT_DERIVATIVES.items():
        tn,tw=derivative_weights(dt);zn,zw=derivative_weights(dz)
        sampled=np.array([[control.value(t*step,z+u*step) for u in zn] for t in tn])
        observed=tw@sampled@zw/step**(dt+dz)
        assert abs(observed-slots[key])<3e-11
    assert control.value(0.,z)==0.
    for s in (-.1,.04,.12):
        assert control.value(s,z)==s*control.w(z,0)+s**3*control.U(z,0)/6


def test_owned_chart_centers_at_rho1_and_keeps_raw_KS_fields():
    control=direction();provider=CompatibleIncomingMetric((.2,),(control,))
    rho=np.array([.8,1.,1.15]);tau=.17;metric=provider.values(tau,rho)
    T1=chart_coordinates(1.)[0]
    assert T1!=0.  # Failure control: PG time and the global chart origin differ.
    for j,x in enumerate(rho):
        T,S=chart_coordinates(x);a=reference_chart(x)[2]
        radius=np.sqrt(1+x*x)+.2*control.value(T-T1,tau+S)
        expected=ks_to_pg(x,1.,0.,a,radius)
        np.testing.assert_array_equal(metric[:,j],expected)
        np.testing.assert_allclose(pg_to_ks(x,*metric[:,j]),[1.,0.,a,radius],rtol=0,atol=3e-14)
    assert metric[3,1]==np.sqrt(2.)
    wrong_origin=np.sqrt(2.)+.2*control.value(0.-T1,tau+chart_coordinates(1.)[1])
    # A wide-window control makes a missing T(1) subtraction visibly wrong.
    wide=CompatibleRadiusDirection(control.w,control.U,2.,3.)
    assert abs(wide.value(T1,tau+chart_coordinates(1.)[1]))>1e-3
    assert provider.log_directions(tau,rho,metric)[0,3,1]==0.


def test_actual_radius_log_tangents_match_independent_amplitude_differences():
    controls=(CompatibleRadiusDirection(wave(.4,.1),wave(.3,-.2),.08,.4),direction())
    amps=np.array([.18,-.13]);rho=np.array([.75,.9,1.,1.1,1.22]);tau=.11
    provider=CompatibleIncomingMetric(amps,controls);metric=provider.values(tau,rho)
    analytic=provider.log_directions(tau,rho,metric)
    assert analytic.shape==(2,4,len(rho));np.testing.assert_array_equal(analytic[:,:3],0.)
    step=2e-5
    def logarithmic(fields):
        out=fields.copy();out[[0,2,3]]=np.log(out[[0,2,3]]);return out
    for i in range(2):
        change=np.zeros(2);change[i]=step
        finite=(logarithmic(CompatibleIncomingMetric(amps+change,controls).values(tau,rho))
                -logarithmic(CompatibleIncomingMetric(amps-change,controls).values(tau,rho)))/(2*step)
        np.testing.assert_allclose(analytic[i],finite,rtol=0,atol=3e-10)
    basis=np.array([[d.value(chart_coordinates(x)[0]-chart_coordinates(1.)[0],tau+chart_coordinates(x)[1]) for x in rho] for d in controls])
    np.testing.assert_array_equal(analytic[:,3],basis/metric[3])
    wrong=basis/np.sqrt(1+rho*rho)
    assert np.max(abs(analytic[:,3]-wrong))>1e-5


def test_smooth_even_window_is_flat_and_exactly_compact():
    control=CompatibleRadiusDirection(polynomial((1.,)),polynomial((0.,)),.2,.4)
    for s in (0.,.1,.2,-.1,-.2):assert control.value(s,.3)==s
    for s in (.4,.5,-.4,-.5):assert control.value(s,.3)==0.
    assert abs(control.value(.3,.3)/.3-.5)<1e-14
    grid=np.linspace(.2,.4,101);window=control.value(grid,.3)/grid
    assert np.all(np.diff(window)<=0)
    np.testing.assert_array_equal(control.value(-grid,.3),-control.value(grid,.3))
    # Flat endpoint derivatives of the numerical control; C-infinity itself
    # follows from its exp(-1/t) construction, not from finitely many samples.
    h=1e-3
    for edge in (.2,.4):
        for order in (1,2,3):
            nodes,weights=derivative_weights(order)
            values=np.array([control.value(edge+h*n,0.)/(edge+h*n) for n in nodes])
            assert abs(weights@values/h**order)<3e-7


@pytest.mark.parametrize('inner,outer',[(0.,1.),(-1.,1.),(1.,1.),(2.,1.),(np.nan,1.),(.1,np.inf),(True,2.)])
def test_invalid_window_radii_are_rejected(inner,outer):
    with pytest.raises(ValueError):CompatibleRadiusDirection(polynomial((1.,)),polynomial((0.,)),inner,outer)


def test_invalid_callbacks_dimensions_values_and_actual_metric_are_rejected():
    with pytest.raises(ValueError,match='callbacks'):CompatibleRadiusDirection(1,polynomial((0.,)),.1,.2)
    bad=CompatibleRadiusDirection(lambda z,n:[1.,2.],polynomial((0.,)),.1,.2)
    with pytest.raises(ValueError,match='scalar'):bad.incoming_slots(0.)
    bad=CompatibleRadiusDirection(lambda z,n:np.nan,polynomial((0.,)),.1,.2)
    with pytest.raises(ValueError,match='finite'):bad.value(.05,0.)
    with pytest.raises(ValueError,match='real'):direction().value(.1+0j,0.)
    with pytest.raises(ValueError,match='finite'):direction().incoming_slots(np.inf)
    with pytest.raises(ValueError,match='one-dimensional'):CompatibleIncomingMetric([[.1]],(direction(),))
    with pytest.raises(ValueError,match='per amplitude'):CompatibleIncomingMetric((.1,.2),(direction(),))
    provider=CompatibleIncomingMetric((.2,),(direction(),))
    with pytest.raises(ValueError,match='one-dimensional'):provider.values(0.,[[.9,1.]])
    with pytest.raises(ValueError,match='scalar'):provider.values([0.],[.9,1.])
    with pytest.raises(ValueError,match='trapped'):provider.values(0.,[3.])
    metric=provider.values(.1,[.9,1.1])
    with pytest.raises(ValueError,match='shape'):provider.log_directions(.1,[.9,1.1],metric[:3])
    changed=metric.copy();changed[3]+=.01
    with pytest.raises(ValueError,match='actual amplitude'):provider.log_directions(.1,[.9,1.1],changed)
    changed=metric.copy();changed[3,0]=0.
    with pytest.raises(ValueError,match='positive'):provider.log_directions(.1,[.9,1.1],changed)
    huge=CompatibleIncomingMetric((-1000.,),(CompatibleRadiusDirection(polynomial((1.,)),polynomial((0.,)),.4,.5),))
    with pytest.raises(ValueError,match='positive'):huge.values(0.,[.9])


def test_explicit_empty_direction_family_and_broadcast_direction_values():
    base=CompatibleIncomingMetric((),());metric=base.values(.1,[.9,1.])
    assert metric.shape==(4,2) and base.log_directions(.1,[.9,1.],metric).shape==(0,4,2)
    control=direction();s=np.array([[-.1],[.1]]);z=np.array([0.,.2,.4])
    values=control.value(s,z)
    assert values.shape==(2,3)
    for i in range(2):
        for j in range(3):assert values[i,j]==control.value(s[i,0],z[j])
