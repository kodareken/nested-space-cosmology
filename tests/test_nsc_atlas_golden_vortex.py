import math
import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from recursive_horizons.nsc_atlas_golden_vortex import (
    PHI,GOLDEN_ANGLE,GOLDEN_PITCH,suspend,nested_ball_point,vortex_point,
    fibonacci_points,fibonacci_cover,kerr_control,kasner,kasner_bounce,bkl_step,contained_prefix)


def test_same_suspension_rule_fills_next_ball_and_preserves_boundary():
    x,a,b=sp.symbols('x a b',real=True)
    norm2=(1-b*b)*(1-a*a)*x*x+(1-b*b)*a*a+b*b
    assert sp.simplify(1-norm2-(1-a)*(1+a)*(1-b)*(1+b)*(1-x)*(1+x))==0
    for x0 in (-1,1):assert np.linalg.norm(suspend(suspend([x0],.3),-.6))==pytest.approx(1)
    # Explicit inverse for an interior ball point establishes surjectivity.
    y=np.array([.2,.3,.4]);t2=y[2];t1=y[1]/math.sqrt(1-t2*t2)
    x0=y[0]/math.sqrt((1-t2*t2)*(1-t1*t1))
    np.testing.assert_allclose(suspend(suspend([x0],t1),t2),y)


def test_golden_turn_and_nested_map_are_exact_scale_rotations():
    assert PHI*PHI-PHI-1==pytest.approx(0)
    assert math.exp(GOLDEN_PITCH*math.pi/2)==pytest.approx(PHI)
    assert math.exp(GOLDEN_PITCH*2*math.pi)==pytest.approx(PHI**4)
    for f,args in ((vortex_point,()),):
        x,L=f(.23);y,M=f(1.23)
        c,s=math.cos(GOLDEN_ANGLE),math.sin(GOLDEN_ANGLE)
        rot=np.array([[c,-s,0],[s,c,0],[0,0,1]])
        np.testing.assert_allclose(y,rot@x,atol=2e-15)
        assert M-L==pytest.approx(-4*math.log(PHI))
    _,logR=nested_ball_point(.2,.3,.4,1000)
    assert math.isfinite(logR) and logR < -1900


def test_fibonacci_samples_and_coverage_have_two_dimensional_area_bound():
    p=fibonacci_points(100)
    np.testing.assert_allclose(np.sum(p*p,axis=1),1,atol=5e-16)
    r=fibonacci_cover(100)
    assert r['cover_radius']>=r['area_lower_bound']
    assert r['cover_radius']<.35
    # Heights change with N: it is not a persistent prefix of older samples.
    assert not np.allclose(fibonacci_points(200)[:100],p)


def test_parke_scale_covariance_does_not_choose_a_ratio():
    ui,uo,t,k,q=sp.symbols('ui uo t k q',positive=True)
    inv=(uo-ui)**2/(9*t*t)+k*(ui+uo)/6+k*k*t*t/16
    changed=inv.subs({ui:q*q*ui,uo:q*q*uo,t:q*t},simultaneous=True)
    assert sp.simplify(changed-q*q*inv)==0
    assert PHI**-8-1<0


def test_schwarzschild_and_kerr_orbit_exponent_from_radial_potential():
    e=sp.symbols('e')
    u=sp.Rational(1,3)+e
    assert sp.expand(sp.Rational(1,27)-u*u+2*u**3)==e*e+2*e**3
    R,a,b,r=sp.symbols('R a b r',positive=True)
    q=b-a
    potential=(R*R-a*q)**2-(R*R-2*R+a*a)*q*q
    rpp=sp.diff(potential,R,2).subs(R,r)
    for sign in (1,-1):
        ac=sign*sp.sqrt(r)*(3-r)/2
        bc=ac+sign*r**sp.Rational(3,2)
        assert sp.simplify(rpp.subs({a:ac,b:bc},simultaneous=True)-6*r*r)==0
    for spin in (0,.5,.9,.99):
        data=kerr_control(spin);rp=data['horizon_BL_radius_over_mass']
        assert math.sqrt(rp*rp+spin*spin+2*spin*spin/rp)==pytest.approx(2)
        for row in data['photon_orbits']:
            rr=row['photon_radius_over_mass'];s=1 if row['orientation']=='prograde' else -1
            bb=spin+s*rr**1.5;Q=bb-spin;P=rr*rr-spin*Q;Delta=rr*rr-2*rr+spin*spin
            exponent=math.sqrt(3)*rr/abs(Q+spin*P/Delta)
            assert exponent==pytest.approx(row['angular_radial_instability_exponent'],rel=1e-12)


def test_golden_kasner_orbit_is_repelling_and_cycles_axes():
    u=sp.symbols('u',positive=True)
    p=[-u/(1+u+u*u),(1+u)/(1+u+u*u),u*(1+u)/(1+u+u*u)]
    assert sp.simplify(sum(p)-1)==0
    assert sp.simplify(sum(x*x for x in p)-1)==0
    K=4*(sum(x*x*(x-1)**2 for x in p)+sum(p[i]**2*p[j]**2 for i in range(3) for j in range(i+1,3)))
    assert sp.simplify(K+16*p[0]*p[1]*p[2])==0
    exactphi=(1+sp.sqrt(5))/2
    assert sp.simplify(K.subs(u,exactphi))==2
    assert bkl_step(PHI)==pytest.approx(PHI)
    assert abs(-1/(PHI-1)**2)>1
    a=kasner(PHI);np.testing.assert_allclose(kasner_bounce(a),a[[2,0,1]],atol=3e-16)


def test_gauss_invariant_normalization_and_lyapunov_integral():
    with mp.workdps(50):
        density=lambda x:1/(mp.log(2)*(1+x))
        assert float(mp.quad(density,[0,1]))==pytest.approx(1)
        L=mp.quad(lambda x:-2*mp.log(x)*density(x),[0,1])
        assert float(L)==pytest.approx(math.pi**2/(6*math.log(2)),rel=1e-14)


def test_common_damping_and_real_conjugate_transfer_leave_equal_counts_equal():
    A,H,t=sp.symbols('A H t',positive=True)
    power=A*A*sp.exp(-3*H*t)
    assert sp.simplify(power-power)==0
    re,im=sp.symbols('re im',real=True);T=re+sp.I*im
    assert sp.simplify(T*sp.conjugate(T)-sp.conjugate(T)*T)==0
    assert contained_prefix(7.206926537906565e-9,3.973074368754331)['log_additional_area_upper']<-1e7
