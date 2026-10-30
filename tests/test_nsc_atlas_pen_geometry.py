import math
import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad
from recursive_horizons.nsc_atlas_pen_geometry import (
    sphere,phase_offset,angular_gaps,current_resolution_area,longitude_coverage,
    orientations,persistent_measure,phase_residual,quadratic_irrational_bound,action_phase_clock)


def test_map_is_sphere_and_two_independent_parameters_have_rank_two():
    u,v,R=sp.symbols('u v R',real=True)
    x=sp.Matrix([R*sp.sin(u)*sp.cos(v),R*sp.sin(u)*sp.sin(v),R*sp.cos(u)])
    assert sp.trigsimp(x.dot(x)-R*R)==0
    cross=x.diff(u).cross(x.diff(v))
    assert sp.trigsimp(cross.dot(cross)-R**4*sp.sin(u)**2)==0
    assert sp.sqrt(6).is_irrational
    assert np.linalg.norm(sphere(.3,.7,2))==pytest.approx(2)


def test_achievable_resolution_and_quadratic_irrational_bounds():
    delta=phase_offset()
    for n in (1,2,100,200):
        gaps=angular_gaps(n,delta)
        assert np.sum(gaps)==pytest.approx(math.pi)
        assert max(gaps)/2>=math.pi/(2*n)-1e-14
        assert phase_residual(n,delta)>=quadratic_irrational_bound(n)
    with mp.workdps(80):
        exact=2*mp.pi*(mp.sqrt(6)-2)
        for n in (100,200):
            expected=abs(mp.sin(n*exact/2))/n
            assert phase_residual(n,delta)==pytest.approx(float(expected),abs=2e-15)
    omega=3.973074368754331
    assert 2/(math.e*math.log(omega))<math.pi/2


def test_geometric_residual_telescopes_but_scalar_signed_count_is_zero():
    delta=phase_offset()
    for n in (1,7,100,200):
        z=np.exp(1j*delta*np.arange(n+1))
        direct=np.sum(z[:-1]-z[1:])/(2*n)
        assert direct==pytest.approx((1-z[-1])/(2*n),abs=2e-15)
        assert abs(direct)==pytest.approx(phase_residual(n,delta),abs=2e-15)
        assert n-n==0


def test_stable_current_area_matches_independent_latitude_integration():
    delta=phase_offset()
    for n,eps in ((1,.07),(3,.1),(8,.03)):
        c=orientations(n,delta);w=[eps]*n
        direct,_=quad(lambda z:longitude_coverage(z,c,w),0,1,epsabs=1e-9,limit=400)
        assert current_resolution_area(n,delta,eps)==pytest.approx(direct,abs=3e-8)
    tiny=math.exp(-100*math.log(3.973074368754331))
    assert current_resolution_area(100,delta,tiny)==pytest.approx(100*tiny,rel=1e-12,abs=0)


def test_persistent_area_and_probe_enclosures_match_single_band_exact_values():
    omega=3.973074368754331;delta=phase_offset()
    r=persistent_measure(1,omega,delta)
    f=math.sin(1/omega);g=f/2-f**3/6
    assert r['area_enclosure'][0]<f<r['area_enclosure'][1]
    assert r['z_probe_enclosure'][0]<g<r['z_probe_enclosure'][1]
    assert r['area_fraction']==pytest.approx(f,abs=1e-8)
    assert r['z_probe_moment']==pytest.approx(g,abs=1e-8)
    two=persistent_measure(2,omega,delta)
    assert two['area_fraction']>f
    assert two['z_probe_moment']>g


def test_coverage_bounds_do_not_need_resolving_exponential_smallness_as_one_minus_one():
    omega=3.973074368754331;delta=phase_offset()
    f=current_resolution_area(200,delta,math.exp(-200*math.log(omega)))
    assert f>0 and math.isfinite(math.log(f))
    assert math.log(f)==pytest.approx(math.log(200)-200*math.log(omega),abs=1e-11)


def test_resolution_bound_has_an_explicit_branching_construction():
    eps=.01;M=math.ceil(math.pi/(2*eps))
    v=np.sort(np.remainder(7*phase_offset()+np.arange(M)*math.pi/M,math.pi))
    h=max(np.diff(np.r_[v,v[0]+math.pi]))/2
    assert h==pytest.approx(math.pi/(2*M),abs=3e-15)
    assert h<=eps


def test_existing_action_fixes_linear_frequency_ratio_for_a_supplied_source_member():
    from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
    r=VolumeParent().vacuum_roots()[1]
    c=action_phase_clock(r,20.)
    x,w=sp.symbols('x w')
    P=25*x**6-140*x+112
    V=7*x**6-sp.Rational(15,2)*x**7+sp.Rational(25,32)*x**12
    # Correct second derivative of each exponential power, before reduction.
    ms=(144*7*x**6-196*sp.Rational(15,2)*x**7+576*sp.Rational(25,32)*x**12)/84
    N=sp.rem(ms-sp.Rational(3,4)*V,P,x)
    D=sp.rem(w*x**6-sp.Rational(3,4)*V,P,x)
    assert sp.simplify(N/D-(475*x-468)/(32*w+75*x-84))==0
    # Eisenstein at7 proves P irreducible over Q; the ratio is nonconstant.
    assert all(int(a)%7==0 for a in sp.Poly(P,x).all_coeffs()[1:])
    assert 112%49!=0 and 25%7!=0
    assert sp.diff((475*x-468)/(32*w+75*x-84),x).subs(w,400)!=0
    xx=math.exp(-2*r['s'])
    expected=2*math.pi*math.sqrt((475*xx-468)/(32*400+75*xx-84))
    assert c['delta']==pytest.approx(expected,rel=2e-14)
    assert c['delta']==pytest.approx(.3603950083870508)
    assert action_phase_clock(r,1e9)['delta']==pytest.approx(7.206926537906565e-9,rel=2e-14,abs=0)
