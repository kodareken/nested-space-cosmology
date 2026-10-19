"""Causal cap, conformal-action and detector checks of the inherited local limit."""
import math

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad

from recursive_horizons.nsc_atlas_causal_reference import signed_state
from recursive_horizons.nsc_atlas_nested_spacetime import expanding_cap_wall_clock
from recursive_horizons.nsc_atlas_dirac_trace import inherited_late_cap_vacuum
from recursive_horizons.nsc_atlas_dirac_preparation import (
    cap_dirac_detector_shape,cap_dirac_limit_factor,cap_inheritance_scales,
    inherit_cap_stress,inherit_cap_two_point,
)


def test_actual_signed_wall_reaches_the_cap_horizon_with_regular_future_clock():
    r=sp.symbols('r',positive=True)
    bm=3*r/4+1/(5*r**2);bp=bm-2*r
    speed2=sp.expand(bm**2-(1-r*r))
    assert sp.simplify(speed2-bp**2-(sp.Rational(4,5)/r-1))==0
    # Positivity on [.45,.8] follows from this exterior identity; on [.8,1]
    # bm>=.8 and fcap<=.36 give v²>=.28; for R>=1, fcap<=0.
    assert bp.subs(r,sp.Rational(4,5))!=0
    assert sp.simplify(speed2.subs(r,1))==sp.Rational(361,400)
    for R in (.45,.8,1.,1.1,1.5):
        state=signed_state(R);v=float(state['Rdot']);beta=float(state['beta_minus'])
        frame=expanding_cap_wall_clock(radius=R,H=1.,beta_minus=beta,Rdot=v)
        t=frame['d_cap_time_d_wall_time']
        # Independently check normalization directly in the areal FLRW metric.
        assert -(1-R*R)*t*t-2*R*t*v+v*v==pytest.approx(-1.,abs=3e-13)
        assert frame['wall_speed_in_cap_time']>R-1
        assert frame['wall_speed_in_cap_time']<R+1
    horizon=expanding_cap_wall_clock(radius=1.,H=1.,beta_minus=.95,Rdot=.95)
    assert horizon['d_cap_time_d_wall_time']==pytest.approx(761/760,abs=3e-14)


def test_massless_full_spinor_conformal_action_and_smooth_inherited_difference():
    # Two co-propagating exact four-component flat Dirac plane solutions.
    # They check propagation/weights, not a selected finite parent vacuum.
    alpha3=np.array([[0,0,1,0],[0,0,0,-1],[1,0,0,0],[0,-1,0,0]],dtype=complex)
    u=np.array([1,0,1,0],dtype=complex)/math.sqrt(2)
    assert np.allclose(alpha3@u,u,atol=3e-14,rtol=0)
    p,q,phase=1.7,.6,.4
    def mode(eta):
        return ((np.exp(-1j*p*eta)+np.exp(1j*phase-1j*q*eta))/math.sqrt(2))[...,None]*u
    def kernel(a,b):
        return mode(a)[..., :,None]*mode(b).conj()[...,None,:]
    n=np.array([0.,.4,1.,2.]);m=n+.2;L=1.3
    transported=inherit_cap_two_point(n,m,kernel,length=L)
    physical_a=np.exp(-1.5*n)[:,None]*mode(-L*np.exp(-n))
    physical_b=np.exp(-1.5*m)[:,None]*mode(-L*np.exp(-m))
    assert transported==pytest.approx(physical_a[:,:,None]*physical_b.conj()[:,None,:],abs=3e-14)
    assert np.all(np.max(np.abs(transported),axis=(-1,-2))<=2*np.exp(-1.5*(n+m)))


def test_stress_weight_from_direct_dirac_energy_and_finite_memory_bound():
    # For chi=a^(3/2)psi, the physical Dirac energy density is
    # a^-3 Re(chi† i*d_t chi)=a^-4 Re(chi† H_flat chi).
    p,q,phase=1.7,.6,.4;L=1.3
    null=np.array([1.,0.,0.,-1.]);tensor=np.outer(null,null)
    def flat_stress(eta):
        density=(p+q)/2*(1+np.cos((p-q)*eta+phase))
        return density[...,None,None]*tensor
    n=np.array([0.,.3,1.,2.,3.]);a=np.exp(n);eta=-L/a
    bound=(p+q)*abs(p-q)/2
    result=inherit_cap_stress(n,flat_stress,length=L,derivative_bound=bound)
    chi=(np.exp(-1j*p*eta)+np.exp(1j*phase-1j*q*eta))/math.sqrt(2)
    i_dt_chi=(p*np.exp(-1j*p*eta)+q*np.exp(1j*phase-1j*q*eta))/(math.sqrt(2)*a)
    direct_energy=(chi.conj()*i_dt_chi).real/a**3
    assert result['stress_difference'][:,0,0]==pytest.approx(direct_energy,abs=3e-13)
    error=np.max(np.abs(result['stress_difference']-result['leading_stress_difference']),axis=(-1,-2))
    assert np.all(error<=result['component_error_bound'])
    # The stated error concerns a supplied derivative norm, not the value of
    # an unspecified parent state. No norm means no numerical error bound.
    assert inherit_cap_stress(0.,flat_stress,length=L)['component_error_bound'] is None


def test_full_four_dimensional_kernel_from_conformal_times_and_KMS_antiperiod():
    L=1.4;s=.37-.19j
    for mean in (0.,.7,2.):
        left,right=mean+s/2,mean-s/2
        eta_left,eta_right=-L*np.exp(-left/L),-L*np.exp(-right/L)
        direct=np.exp(-1.5*(left+right)/L)/(eta_left-eta_right)**3
        assert cap_dirac_limit_factor(s,length=L)==pytest.approx(direct,rel=3e-13)
    assert cap_dirac_limit_factor(s-2j*math.pi*L,length=L)==pytest.approx(
        -cap_dirac_limit_factor(s,length=L),rel=3e-13)
    assert cap_dirac_limit_factor(3*s,length=3*L)==pytest.approx(
        cap_dirac_limit_factor(s,length=L)/27,rel=3e-13)


def test_bilinear_detector_spectrum_from_independent_shifted_contour_integral():
    # Shift s->y-i*pi: [2i*sinh((s-i0)/2)]^-6 becomes
    # [2*cosh(y/2)]^-6 and the Fourier weight acquires exp(-pi*x).
    # This nonsingular quadrature is independent of the Planck polynomial.
    for x in (0.,.35,1.,-.7):
        integral,error=quad(lambda y:2*math.cos(x*y)/(2*math.cosh(y/2))**6,
                            0.,30.,epsabs=2e-13,epsrel=2e-13)
        independent=60/math.pi*math.exp(-math.pi*x)*integral
        assert error<3e-14
        assert cap_dirac_detector_shape(x)==pytest.approx(independent,rel=8e-13,abs=4e-14)
    x=np.array([.1,1.,3.])
    ratio=cap_dirac_detector_shape(x)/cap_dirac_detector_shape(-x)
    assert ratio==pytest.approx(np.exp(-2*math.pi*x),rel=3e-13,abs=0)
    scales=cap_inheritance_scales(10000.)
    assert scales['kms_time_s']==pytest.approx(2*math.pi*10000/299792458.,rel=3e-14,abs=0)
    assert scales['temperature_K']*scales['kms_time_s']==pytest.approx(
        6.62607015e-34/(2*math.pi*1.380649e-23),rel=3e-14,abs=0)
    scaled=cap_inheritance_scales(70000.)
    assert scaled['temperature_K']==pytest.approx(scales['temperature_K']/7,rel=3e-14,abs=0)
    assert scaled['stress_decay_time_s']==pytest.approx(7*scales['stress_decay_time_s'],rel=3e-14,abs=0)


def test_inheritance_calculator_keeps_its_domains():
    with pytest.raises(ValueError):
        inherit_cap_two_point(-1.,0.,lambda a,b:np.eye(4))
    with pytest.raises(ValueError):
        inherit_cap_two_point(0.,0.,lambda a,b:np.eye(2))
    with pytest.raises(ValueError):
        inherit_cap_stress(0.,lambda a:np.eye(4),derivative_bound=-1.)
    with pytest.raises(ValueError,match='pole'):
        cap_dirac_limit_factor(0.)
    with pytest.raises(ValueError):
        cap_dirac_detector_shape(1j)
    with pytest.raises(ValueError):
        cap_inheritance_scales(0.)


def test_late_vacuum_form_follows_trace_conservation_and_fixed_couplings():
    # Independent FLRW Ward identity: -rho+3p=A and
    # d_t rho+3H(rho+p)=0 imply d_t rho+4H rho=-H A.
    # A constant late limit is rho=-A/4, p=A/4. The free radiation
    # integration constant is precisely the decaying a^-4 state term.
    H,A,C,t=sp.symbols('H A C t',real=True)
    rho=-A/4+C*sp.exp(-4*H*t);p=(A+rho)/3
    assert sp.simplify(sp.diff(rho,t)+3*H*(rho+p))==0
    result=inherited_late_cap_vacuum(ell_m=10000.)
    trace=-11/(240*math.pi**2)
    assert result['trace_times_ell4_over_hbar_c']==pytest.approx(trace,rel=3e-14,abs=0)
    assert result['energy_density_times_ell4_over_hbar_c']==pytest.approx(-trace/4,rel=3e-14,abs=0)
    lP2=6.62607015e-34*6.67430e-11/(2*math.pi*299792458.**3)
    assert result['delta_Lambda_times_ell2']==pytest.approx(
        11*lP2/(120*math.pi*10000.**2),rel=3e-14,abs=0)
    assert result['nested_minus_smooth_late_vacuum_source']==0
    assert inherited_late_cap_vacuum(ell_m=30000.)['delta_Lambda_times_ell2']==pytest.approx(
        result['delta_Lambda_times_ell2']/9,rel=3e-14,abs=0)
