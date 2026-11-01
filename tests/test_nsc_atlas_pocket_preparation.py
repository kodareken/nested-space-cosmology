"""Independent local-force and dense initial-constraint controls; no evolution."""
import math

import numpy as np
import pytest
from scipy.special import eval_gegenbauer,roots_jacobi
import sympy as sp
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_sphere_flux_completion import conditional_volume_barrier
from recursive_horizons.nsc_atlas_valve_neck import solve_pocket_data
from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve
from recursive_horizons.nsc_atlas_shape_transition import flat_horizon_pp_record


def test_conditional_top_and_analytic_derivative_satisfy_implicit_force_equation():
    chi=np.array([-.02,-.007,0.,.015]);top,derivative=conditional_volume_barrier(chi)
    A=7+200*chi**2
    force=-12*A*np.exp(-12*top)+105*np.exp(-14*top)-75/4*np.exp(-24*top)
    stiffness=144*A*np.exp(-12*top)-1470*np.exp(-14*top)+450*np.exp(-24*top)
    mixed=-4800*chi*np.exp(-12*top)
    np.testing.assert_allclose(force,0.,atol=2e-12,rtol=0)
    assert np.all(stiffness<0)
    np.testing.assert_allclose(derivative,-mixed/stiffness,rtol=3e-13,atol=3e-15)
    h=1e-6;point=.012
    numeric=(conditional_volume_barrier(point+h)[0]-conditional_volume_barrier(point-h)[0])/(2*h)
    assert conditional_volume_barrier(point)[1]==pytest.approx(numeric,rel=3e-8,abs=2e-9)


def test_conditional_fold_requires_dynamics_instead_of_a_static_top():
    B=15/2;T=25/32;turn=(7*B/(72*T))**.2;critical_A=35*B*turn/36
    fold=math.sqrt((critical_A-7)/200)
    with pytest.raises(ValueError,match='dynamical preparation'):
        conditional_volume_barrier(1.001*fold)


@pytest.fixture(scope='module')
def preparations():
    with threadpool_limits(limits=1):
        return dict(pure=solve_pocket_data(40,epsilon=.02),
            carrier=solve_pocket_data(40,epsilon=.02,carrier_scale=.25),
            compressed=solve_pocket_data(40,epsilon=0.,carrier_scale=.25,compression=.01))


def independent_fields(data,order=401):
    # Independent higher Gauss-Jacobi rule; the producer used241 nodes.
    u,w=roots_jacobi(order,.5,.5);w*=4*math.pi;v=1-u*u;H=data.parent.H
    psi=data.values(u);du=data.values(u,1);duu=data.values(u,2)
    state=data.state_values(u);S=data.S_values(u);Su=data.S_values(u,1)
    labels=data.S_labels;normal=math.sqrt(2*math.pi**2)
    # Analytic third derivatives of Gegenbauer polynomials independently
    # differentiate the longitudinal potential; no waveform samples are used.
    Fuu=sum(c/(l*(l+2)-3)*8*eval_gegenbauer(int(l)-2,3,u)/normal
        for c,l in zip(data.S_coefficients,labels) if l>=2)
    Fuuu=sum(c/(l*(l+2)-3)*48*eval_gegenbauer(int(l)-3,4,u)/normal
        for c,l in zip(data.S_coefficients,labels) if l>=3)
    coverS=data.parent.kappa*data.PiS_conformal*Fuu
    coverSu=data.parent.kappa*data.PiS_conformal*Fuuu
    C=-data.parent.kappa*data.carrier_A*data.carrier_B*math.sqrt(3/5)
    ArS=v*coverS;ArChi=C*u*v/6;Ar=ArS+ArChi
    Aru=-2*u*coverS+v*coverSu+C*(1-3*u*u)/6
    Kr=state[2]+v*state[3];Kt=state[2];PiS,PiChi=state[6:8]
    Ktu=-.5*Aru*psi**-6+3*Ar*psi**-7*du
    radius_scale=psi**2/H
    s=data.sstar+S/math.sqrt(84);chi=data.carrier_A*u
    V=7*np.exp(-12*s)-7.5*np.exp(-14*s)+25/32*np.exp(-24*s)+200*np.exp(-12*s)*chi**2
    gradient=H*H*v*psi**-4*(Su*Su+data.carrier_A**2)
    rho=.5*(PiS**2+PiChi**2+gradient)+V
    scalarR=H*H*psi**-5*(-8*(v*duu-3*u*du)+6*psi)
    Ham=scalarR+(Kr+2*Kt)**2-(Kr*Kr+2*Kt*Kt)-2*data.parent.kappa*rho
    j=np.sqrt(v)/radius_scale*(PiS*Su+PiChi*data.carrier_A)
    momentum=2*((u/np.sqrt(v)-2*np.sqrt(v)*du/psi)*(Kr-Kt)+np.sqrt(v)*Ktu)/radius_scale-data.parent.kappa*j
    measure=w*psi**6/H**3
    return dict(u=u,w=w,v=v,H=H,psi=psi,state=state,S=S,Su=Su,ArS=ArS,ArChi=ArChi,Ar=Ar,Aru=Aru,
        Kr=Kr,Kt=Kt,PiS=PiS,PiChi=PiChi,Ham=Ham,momentum=momentum,measure=measure)


@pytest.mark.parametrize('name',('pure','carrier','compressed'))
def test_fresh_preparations_satisfy_independent_dense_Hamiltonian_and_momentum(preparations,name):
    data=preparations[name];f=independent_fields(data)
    rms=lambda a:math.sqrt(float(f['measure']@(a*a)/sum(f['measure'])))/f['H']**2
    assert rms(f['Ham'])<2e-8
    assert rms(f['momentum'])<5e-10
    assert min(f['psi'])>0 and data.positivity_lower_bound>0
    np.testing.assert_allclose(f['Kr'],f['Ar']*f['psi']**-6,rtol=3e-13,atol=3e-15)
    np.testing.assert_allclose(f['Kt'],-.5*f['Ar']*f['psi']**-6,rtol=3e-13,atol=3e-15)


def test_joint_York_sources_and_cross_shear_are_retained(preparations):
    data=preparations['compressed'];f=independent_fields(data);u=f['u'];v=f['v'];psi=f['psi'];H=f['H']
    div=H*(-np.sqrt(v)*f['Aru']+3*u/np.sqrt(v)*f['Ar'])
    PiSbar=f['PiS']*psi**6;PiChibar=f['PiChi']*psi**6
    source=data.parent.kappa*H*np.sqrt(v)*(PiSbar*f['Su']+PiChibar*data.carrier_A)
    np.testing.assert_allclose(div,source,rtol=2e-11,atol=3e-13)
    np.testing.assert_allclose(data.radial_shear_over_v(u),(f['ArS']+f['ArChi'])/v,rtol=3e-13,atol=3e-15)
    combined=f['Kr']**2+2*f['Kt']**2
    separate=1.5*(f['ArS']**2+f['ArChi']**2)*psi**-12
    cross=3*f['ArS']*f['ArChi']*psi**-12
    assert np.max(abs(cross))>1e-10
    np.testing.assert_allclose(combined-separate,cross,rtol=2e-10,atol=3e-17)


def test_compressed_carrier_keeps_odd_geometry_and_nonzero_radion_momentum(preparations):
    data=preparations['compressed'];u=np.linspace(-1.,1.,101)
    assert np.array_equal(data.labels,np.arange(41))
    odd=data.coefficients[data.labels%2==1]
    assert np.linalg.norm(odd)>1000*np.finfo(float).eps*np.linalg.norm(data.coefficients)
    assert np.max(abs(data.values(u)-data.values(-u)))>1e-10
    PiS=data.state_values(u)[6]
    assert np.max(PiS)<0
    np.testing.assert_allclose(PiS*data.values(u)**6,-math.sqrt(84)*.01,rtol=3e-15,atol=3e-17)


def test_directed_carrier_refuses_old_reflection_mask(preparations):
    data=preparations['carrier']
    Pi=data.state_values(np.linspace(-.9,.9,21))[7]
    np.testing.assert_allclose(Pi,Pi[::-1],rtol=2e-13,atol=3e-15)
    assert np.max(abs(Pi))>0
    with pytest.raises(ValueError,match='both reflection parities'):SphericalValve(64,initial_data=data)
    system=SphericalValve(64,initial_data=data,reflection_symmetry=False)
    np.testing.assert_array_equal(system.initial,data.state_values(system.u))
    assert system.calls==0


def test_flat_departure_has_divergent_parallel_tides_despite_finite_scalar_limits():
    row=flat_horizon_pp_record();stop=conditional_volume_barrier(0.)[0]
    V=7*math.exp(-12*stop)-7.5*math.exp(-14*stop)+25/32*math.exp(-24*stop)
    Vss=1008*math.exp(-12*stop)-1470*math.exp(-14*stop)+450*math.exp(-24*stop)
    alpha=-Vss/(84*(V/3));p=(-3+math.sqrt(9+4*alpha))/2
    assert row['barrier_s']==pytest.approx(stop,abs=2e-14)
    assert row['minus_mass_squared_over_H_squared']==pytest.approx(alpha,rel=3e-14)
    assert row['departure_exponent']==pytest.approx(p,rel=3e-14)
    # At a stationary point, imposing alpha=4 independently gives x^5
    # =7B/(132T), and A=x(7B/6-2T*x^5). No radius is hardcoded here.
    B=sp.Rational(15,2);T=sp.Rational(25,32)
    x5=7*B/(132*T);prefactor=7*B/6-2*T*x5
    threshold5=prefactor**5*x5
    merger5=(sp.Rational(35,36)*B)**5*(7*B/(72*T))
    assert threshold5<7**5<merger5
    assert sp.Rational(row['exact_p1_threshold_A_fifth'])==threshold5
    assert sp.Rational(row['exact_merger_A_fifth'])==merger5
    assert 0<p<1 and row['pp_power']==pytest.approx(p-2,abs=3e-15)
    assert row['pp_power']<-1 and row['internal_pp_coefficient_per_C']>0

    lam,power,C,s0=sp.symbols('lambda power C s0',positive=True)
    s=s0+C*lam**power;b=sp.exp(s)
    A=lam*(1+6*(1-power)/(power+1)*C*lam**power)
    normalization=C*lam**(power-2)
    internal=sp.simplify(-sp.diff(b,lam,2)/b/normalization)
    external=sp.simplify(-sp.diff(A,lam,2)/A/normalization)
    internal_limit=sp.limit(internal,lam,0,dir='+')
    external_limit=sp.limit(external,lam,0,dir='+')
    assert sp.simplify(12*internal_limit+2*external_limit)==0
    ci=float(internal_limit.subs(power,p));ce=float(external_limit.subs(power,p))
    assert row['internal_pp_coefficient_per_C']==pytest.approx(ci,rel=3e-14)
    assert row['external_pp_coefficient_per_C']==pytest.approx(ce,rel=3e-14)
    assert row['leading_null_Ricci_trace']==0.
    normalized_i=sp.lambdify((lam,power,C),internal,'numpy')
    normalized_e=sp.lambdify((lam,power,C),external,'numpy')
    errors_i=[];errors_e=[];tidal=[]
    for affine in (.01,.001):
        actual_i=float(normalized_i(affine,p,.2));actual_e=float(normalized_e(affine,p,.2))
        errors_i.append(abs(actual_i-ci));errors_e.append(abs(actual_e-ce))
        tidal.append(abs(.2*affine**(p-2)*actual_i))
    assert errors_i[1]<errors_i[0]/5 and errors_i[1]<.01*abs(ci)
    assert errors_e[1]<errors_e[0]/5 and errors_e[1]<.001*abs(ce)
    assert tidal[1]>10*tidal[0]
    # Independent product-space scalar limits: dS4 plus four static S3s.
    Hext2=math.exp(12*stop)*V/3;section=(5/8)*math.exp(-2*stop)
    scalar=4*3*Hext2+4*(3*2*section)
    kretschmann=2*4*3*Hext2**2+4*(2*3*2*section**2)
    assert row['full_Ricci_limit']==pytest.approx(scalar,rel=3e-14)
    assert row['full_Kretschmann_limit']==pytest.approx(kretschmann,rel=3e-14)
    assert math.isfinite(scalar) and math.isfinite(kretschmann)
    # Finite contractions coexist with unbounded parallel-frame components;
    # they cannot certify a C2 horizon extension of the departing solution.
