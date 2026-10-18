"""Independent curvature, profile and conserved-source controls."""
import math
import time

import mpmath as mp
import pytest
import sympy as sp

from recursive_horizons.nsc_covariant_spectral_response import (
    curvature_profile_kernels, static_gravitational_response,
)
from recursive_horizons.nsc_covariant_identities import dirac_heat_coefficients
from recursive_horizons.nsc_spectral_gravitational_response import flat_spin2_profile_response


ACTIVE_CPU=[0.0]


@pytest.fixture(autouse=True)
def bounded():
    started=time.process_time()
    yield
    ACTIVE_CPU[0]+=time.process_time()-started
    assert ACTIVE_CPU[0]<60,'covariant response controls exhausted their allocation'


def test_dirac_projection_and_local_scalar_cancellation():
    z,f=sp.symbols('z f',positive=True)
    ric=1/(6*z)+(f-1)/z**2
    rr=f/32+f/(8*z)-sp.Rational(7,48)/z-(f-1)/(8*z**2)
    ru=-f/4-(f-1)/(2*z);uu=f/2;om=-(f-1)/(2*z)
    cc=2*ric-om
    scalar=sp.Rational(4,3)*ric+4*rr+ru+uu/4-om/6
    assert sp.cancel(cc-((3*f-1)/(6*z)+2*(f-1)/z**2))==0
    assert sp.cancel(scalar-((3*f+2)/(36*z)+5*(f-1)/(6*z**2)))==0
    series=1-z/6+z**2/60-z**3/840+z**4/15120
    assert sp.limit(cc.subs(f,series),z,0)==-sp.Rational(1,20)
    assert sp.limit(scalar.subs(f,series),z,0)==0
    assert sp.limit(scalar.subs(f,series)/z,z,0)==sp.Rational(1,2520)
    owned=dirac_heat_coefficients()
    assert sp.sympify(owned['spin_curvature_riemann_ratio'])==-sp.Rational(1,2)
    assert owned['spin_curvature_contraction_residual']=='0'


def test_profile_integral_exchange_and_scalar_zero_mode():
    x,u,t,v=sp.symbols('x u t v',positive=True)
    reach=sp.sqrt(1-4*u/x)
    # Exchange the two antiderivative integrals of the imported scalar heat function.
    domain=sp.integrate(sp.Rational(1,2),(t,-reach,reach))
    linear=sp.integrate((x*(1-t*t)/4-u)/2,(t,-reach,reach))
    weight=-domain/(12*x)+5*linear/(6*x*x)
    assert sp.simplify(weight-reach*(1-10*u/x)/(18*x))==0
    assert sp.integrate((1-sp.Rational(5,2)*v)*sp.sqrt(1-v),(v,0,1))==0
    assert sp.integrate(v*(1-sp.Rational(5,2)*v)*sp.sqrt(1-v),(v,0,1))/288==-sp.Rational(1,2520)


def test_profile_kernels_against_independent_heat_integral():
    with mp.workdps(45):
        x=mp.mpf('0.6');lo=mp.mpf('0.2');hi=mp.mpf('0.7')
        # Closed value of the independent heat parameter integral.
        def basic(z):return mp.sqrt(mp.pi/z)*mp.exp(-z/4)*mp.erfi(mp.sqrt(z)/2)
        def heat(t,scalar):
            z=t*x;f=basic(z)
            return ((3*f+2)/(36*z)+5*(f-1)/(6*z*z) if scalar else
                    (3*f-1)/(6*z)+2*(f-1)/(z*z))
        # A finite proper-time window has a signed-subtracted E1 profile.
        expected=[float(mp.quad(lambda t:heat(t,s)/t,[lo,hi])) for s in (False,True)]
        def profile(y):return float(mp.e1(lo*y)-mp.e1(hi*y))
        h0=float(mp.log(hi/lo))
        result=curvature_profile_kernels(float(x),profile,h0,quadrature_order=24)
        assert result['weyl_kernel']==pytest.approx(expected[0],abs=2e-13)
        assert result['scalar_kernel']==pytest.approx(expected[1],abs=2e-13)


def _linear_curvatures(h,p):
    d=4
    christ=lambda a,b,c:sp.I*(p[b]*h[a,c]+p[c]*h[a,b]-p[a]*h[b,c])/2
    riem={(a,b,c,e):sp.expand(sp.I*p[c]*christ(a,e,b)-sp.I*p[e]*christ(a,c,b))
          for a in range(d) for b in range(d) for c in range(d) for e in range(d)}
    ric=sp.Matrix(d,d,lambda b,e:sum(riem[a,b,a,e] for a in range(d)))
    R=sp.trace(ric)
    C2=sum(v*v for v in riem.values())-2*sp.trace(ric*ric)+R*R/3
    return sp.simplify(R),sp.factor(C2)


def test_static_metric_invariants_and_pure_gauge_control():
    phi,psi,k=sp.symbols('Phi Psi k',real=True)
    h=sp.diag(2*phi,-2*psi,-2*psi,-2*psi)
    R,C2=_linear_curvatures(h,[0,0,0,k])
    assert sp.expand(R-2*k*k*(phi-2*psi))==0
    assert sp.expand(C2-sp.Rational(4,3)*k**4*(phi+psi)**2)==0
    xi=sp.Matrix(sp.symbols('xi0:4',real=True));p=sp.Matrix([0,0,0,k])
    gauge=sp.I*(p*xi.T+xi*p.T)
    assert _linear_curvatures(gauge,p)==(0,0)


def test_conserved_static_source_response_by_action_variation():
    phi,psi,rho,A,x,fc,fr=sp.symbols('Phi Psi rho A x fC fR',nonzero=True,real=True)
    norm=sp.symbols('norm',positive=True)
    action=A*x*(4*phi*psi-2*psi**2)+(sp.Rational(4,3)*x*x*fc*(phi+psi)**2
            +4*x*x*fr*(phi-2*psi)**2)/(2*norm)+rho*phi
    a2,a0=1+x*fc/(norm*A),1-3*x*fr/(norm*A)
    base=-rho/(4*A*x)
    values={phi:base*(4/(3*a2)-1/(3*a0)),psi:base*(2/(3*a2)+1/(3*a0))}
    assert sp.cancel(sp.diff(action,phi).subs(values,simultaneous=True))==0
    assert sp.cancel(sp.diff(action,psi).subs(values,simultaneous=True))==0
    chosen={x:sp.Rational(1,5),fc:-sp.Rational(1,4),fr:sp.Rational(1,20),A:2,norm:3}
    result=static_gravitational_response(float(chosen[x]),float(chosen[fc]),float(chosen[fr]),36)
    assert result['Newton_potential_ratio']==pytest.approx(float((values[phi]/base).subs(chosen)))
    assert result['spatial_potential_ratio']==pytest.approx(float((values[psi]/base).subs(chosen)))
    assert result['lensing_potential_ratio']==pytest.approx(1/result['spin2_factor'])


def test_Einstein_density_from_connection_with_explicit_boundary_flux():
    from functools import lru_cache
    z,eps,A=sp.symbols('z epsilon A',real=True)
    N,S=sp.Function('N')(z),sp.Function('S')(z)
    phi,psi=sp.Function('Phi')(z),sp.Function('Psi')(z)
    metric=sp.diag(N*N,S*S,S*S,S*S)
    inverse=[1/metric[i,i] for i in range(4)]
    derivative=lambda i,f:sp.diff(f,z) if i==3 else 0
    @lru_cache(None)
    def connection(a,b,c):
        return inverse[a]*(derivative(b,metric[a,c])+derivative(c,metric[a,b])
                           -derivative(a,metric[b,c]))/2
    scalar=0
    for b in range(4):
        ricci=sum(derivative(a,connection(a,b,b))-derivative(b,connection(a,a,b))
            +sum(connection(a,a,d)*connection(d,b,b)-connection(a,b,d)*connection(d,a,b)
                 for d in range(4)) for a in range(4))
        scalar+=inverse[b]*ricci
    density=(-A*N*S**3*scalar).subs({N:sp.exp(eps*phi),S:sp.exp(-eps*psi)},simultaneous=True).doit()
    quadratic=sp.diff(density,eps,2).subs(eps,0)/2
    bulk=A*(4*sp.diff(phi,z)*sp.diff(psi,z)-2*sp.diff(psi,z)**2)
    flux=A*sp.diff((phi-psi)*(2*sp.diff(phi,z)-4*sp.diff(psi,z)),z)
    assert sp.expand(quadratic-bulk-flux)==0


def test_TT_recovery_units_and_known_limits():
    profile=lambda y:math.exp(-0.8*y);q1=1/0.8;x=0.3
    kernels=curvature_profile_kernels(x,profile,1)
    response=static_gravitational_response(x,kernels['weyl_kernel'],kernels['scalar_kernel'],q1)
    original=flat_spin2_profile_response(x,profile,q1)
    assert response['spin2_factor']==pytest.approx(original['response_over_Einstein'],abs=1e-14)
    scale=2.3
    changed=curvature_profile_kernels(x/scale**2,lambda y:profile(scale**2*y),1)
    other=static_gravitational_response(x/scale**2,changed['weyl_kernel'],changed['scalar_kernel'],q1/scale**2)
    assert other==pytest.approx(response,abs=1e-14)
    zero=curvature_profile_kernels(0,lambda _:pytest.fail('no zero sample needed'),1)
    assert zero['weyl_kernel']==-0.05 and zero['scalar_kernel']==0
    GR=static_gravitational_response(x,0,0,q1)
    assert all(v==1 for v in GR.values())


def test_invalid_domains_and_singular_response():
    for x,h0,n in [(-1,1,8),(1,float('nan'),8),(1,1,True),(1,1,1)]:
        with pytest.raises(ValueError):curvature_profile_kernels(x,lambda _:1,h0,quadrature_order=n)
    with pytest.raises(ValueError):curvature_profile_kernels(1,lambda _:float('nan'),1)
    with pytest.raises(ValueError):static_gravitational_response(1,-1/6,0,1)
    with pytest.raises(ValueError):static_gravitational_response(1,0,1/18,1)
