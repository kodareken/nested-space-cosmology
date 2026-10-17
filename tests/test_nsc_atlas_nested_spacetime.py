"""Independent curvature, orientation and detector-clock controls."""
import math

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_nested_spacetime import (
    central_inward_flash, regular_kottler_chart, static_slice_frame,
)


def test_regular_metric_einstein_equations_from_warped_curvature():
    # Generic exact symbols independently differentiate the metric and r;
    # no target host parameters or numerical wall application is evaluated.
    u,w,C,L=s.symbols('u w C L', nonzero=True, real=True)
    k=(1-3*L*C*C)/(2*C);r=C+k*u*w
    F=w*w/r*(1/C+L*(r-C))
    Fu,Fw,Fww=s.diff(F,u),s.diff(F,w),s.diff(F,w,2)
    hu=Fw*k*w/2-(Fu+F*Fw)*k*u/2
    huw=k-Fw*k*u/2
    radial_norm=2*k*k*u*w-F*k*k*u*u
    f=1-(C-L*C**3)/r-L*r*r
    for residual in (radial_norm-f, Fww*F/2-2*hu/r-3*L*F,
                     Fww/2-2*huw/r-3*L, 1-f-2*r*huw-3*L*r*r):
        assert s.factor(residual)==0


@pytest.mark.parametrize('C,H', [(2.,0.),(1.,.25),(3.,.25)])
def test_chart_jets_and_mass_invariant(C,H):
    # Includes a negative-surface-gravity cosmological horizon.
    mass=C*(1-H*H*C*C)/2
    u,w=s.symbols('u w');k=(1-3*H*H*C*C)/(2*C)
    r=C+k*u*w;F=w*w/r*(1/C+H*H*(r-C))
    row=regular_kottler_chart(.2,-.3,mass=mass,H=H,horizon=C)
    for name,expression in [('F',F),('F_u',s.diff(F,u)),
                            ('F_w',s.diff(F,w)),('F_ww',s.diff(F,w,2))]:
        expected=float(expression.subs({u:.2,w:-.3}))
        assert float(row[name])==pytest.approx(expected,rel=2e-13,abs=1e-14)
    norm=2*row['R_u']*row['R_w']-row['F']*row['R_w']**2
    assert norm==pytest.approx(1-2*mass/row['R']-(H*row['R'])**2,abs=3e-15)
    axis=regular_kottler_chart(np.array([0.,.1,0.]),np.array([.2,0.,0.]),
                               mass=mass,H=H,horizon=C)
    assert np.all(axis['R']==C)
    assert np.all(np.isfinite(axis['F']))
    assert axis['radial_metric_determinant']==-1


@pytest.mark.parametrize('orientation',[-1,1])
def test_oriented_frame_preserves_one_killing_time(orientation):
    q=static_slice_frame(3.,.2,mass=1.,H=0.,orientation=orientation)
    f=q['f'];metric=np.diag([-f,1/f])
    normal=np.array([q['N_T'],q['N_R']]);tangent=np.array([q['s_T'],q['s_R']])
    assert normal@metric@normal==pytest.approx(-1,abs=1e-15)
    assert tangent@metric@tangent==pytest.approx(1,abs=1e-15)
    assert normal@metric@tangent==pytest.approx(0,abs=1e-15)
    assert np.linalg.det(np.column_stack([normal,tangent]))==pytest.approx(1)
    assert math.copysign(1,q['N_T'])==orientation
    assert q['outgoing_l']==pytest.approx(q['T_l']-q['Rprime']/f)


def test_flash_flat_limit_and_independent_doppler_factor():
    v=.3;g=1/math.sqrt(1-v*v)
    flat=central_inward_flash(radius=2.,H=0.,beta_minus=g,Rdot=g*v)
    assert flat['central_proper_time']==2.
    assert flat['one_plus_redshift']==pytest.approx(math.sqrt((1+v)/(1-v)))
    cap=central_inward_flash(radius=.5,H=.5,beta_minus=1.,Rdot=.25)
    assert cap['central_proper_time']==pytest.approx(2*math.log(4/3))
    assert cap['one_plus_redshift']==pytest.approx(4/3)


def test_flash_and_chart_length_scaling():
    flash=central_inward_flash(radius=.5,H=.5,beta_minus=1.,Rdot=.25)
    scaled=central_inward_flash(radius=3.5,H=.5/7,beta_minus=1.,Rdot=.25)
    assert scaled['central_proper_time']==pytest.approx(7*flash['central_proper_time'])
    assert scaled['one_plus_redshift']==flash['one_plus_redshift']
    base=regular_kottler_chart(.2,-.3,mass=1.,H=0.,horizon=2.)
    big=regular_kottler_chart(1.4,-2.1,mass=7.,H=0.,horizon=14.)
    assert big['R']==pytest.approx(7*base['R'])
    assert big['F']==pytest.approx(base['F'])
    assert big['F_ww']==pytest.approx(base['F_ww']/49)


def test_invalid_chart_and_flash_domains_stay_explicit():
    with pytest.raises(ValueError,match='root'):
        regular_kottler_chart(0.,0.,mass=1.,H=0.,horizon=3.)
    with pytest.raises(ValueError,match='simple'):
        regular_kottler_chart(0.,0.,mass=1/3,H=1/math.sqrt(3),horizon=1.)
    with pytest.raises(ValueError,match='R>0'):
        regular_kottler_chart(10.,-10.,mass=1.,H=0.,horizon=2.)
    with pytest.raises(ValueError,match='f>0'):
        static_slice_frame(2.,0.,mass=1.,H=0.,orientation=1)
    with pytest.raises(ValueError,match=r'H\*R<1'):
        central_inward_flash(radius=1.,H=1.,beta_minus=1.,Rdot=1.)
    with pytest.raises(ValueError,match='factors'):
        central_inward_flash(radius=.5,H=.5,beta_minus=2.,Rdot=.25)
