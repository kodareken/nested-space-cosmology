"""Independent action identities, exact enclosures and physical clock controls."""
from fractions import Fraction as F
import math

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_inherited_matter import (
    exact_initial_geometry,uniform_bound,dust_lightcone,scalar_lightcone_bound,
    outward_decimal_interval,solve_parent,
)
from recursive_horizons.nsc_atlas_shared_scalar_slice import SharedScalarSlice


def test_metric_source_corrected_energy_identity():
    x,p,A,H,w,K,kap=sp.symbols('x p A H w K kap',real=True)
    E=(p*p+w*w*x*x)/2;pressure=(p*p-w*w*x*x)/2
    D=A**3*E;S=A**3*x*p;I=D+sp.Rational(3,2)*H*S
    rates=(p,-3*H*p-w*w*x,A*H,K/A**2-kap*p*p/2)
    dt=lambda f:sum(sp.diff(f,z)*r for z,r in zip((x,p,A,H),rates))
    assert sp.expand(dt(S)-2*A**3*pressure)==0
    assert sp.expand(dt(D)+3*H*A**3*pressure)==0
    assert sp.simplify(dt(I)-sp.Rational(3,2)*rates[-1]*S)==0
    assert sp.expand(E*E-w*w*x*x*p*p-pressure*pressure)==0


def test_exact_junction_and_ancestral_energy_fix_the_normalization():
    d=exact_initial_geometry();atlas=SharedScalarSlice();p=atlas.parent_causal_domain()
    assert math.isclose(float(d['H02']),p['initial_H']**2,rel_tol=3e-15)
    assert math.isclose(float(d['curvature']),1/p['initial_scale_factor']**2,rel_tol=3e-14)
    eps,e=F(1,4),F(1,2)
    assert (1-eps)*(2*e)**2/2==d['excitation']==F(3,8)
    assert d['excitation']/d['vacuum']==F(4,3)


@pytest.mark.parametrize('omega',[20,2000,10**9])
def test_rational_comparison_is_positive_and_its_decimal_interval_is_outward(omega):
    b=uniform_bound(omega)
    assert 0<b['lower']<1<b['upper']
    q=scalar_lightcone_bound(omega)
    text=outward_decimal_interval(q['lower'],q['upper'])
    assert F(text[0])<=q['lower']<=q['upper']<=F(text[1])
    if omega==10**9:
        assert b['density_relative_error']<F(1097,10**12)
        assert q['relative_squared_distance_error']<F(7,10**10)


def test_independent_high_precision_metric_integral_is_inside_the_rational_photon_interval():
    d=exact_initial_geometry();bound=dust_lightcone()
    with mp.workdps(70):
        number=lambda x:mp.mpf(x.numerator)/x.denominator
        K,vac,E0=map(number,(d['curvature'],d['vacuum'],d['excitation']))
        eta=mp.sqrt(K)*mp.quad(lambda u:1/mp.sqrt(vac/3+E0*u**3/3-K*u*u),
            [mp.mpf('0.5'),mp.mpf('0.625')])
        Q=vac/K*(mp.mpf('2.5')*mp.sin(eta))**2
        assert number(bound['lower'])<Q<number(bound['upper'])
    assert dust_lightcone(redshift=F(0))['upper']==0
    with pytest.raises(ValueError):dust_lightcone(redshift=F(2))


def test_constant_H_analytic_solution_keeps_the_corrected_energy_exactly_constant():
    H,w,x0=.3,20.,.04;nu=math.sqrt(w*w-9*H*H/4)
    t=np.linspace(0,10,1001)
    x=x0*np.exp(-1.5*H*t)*(np.cos(nu*t)+1.5*H/nu*np.sin(nu*t))
    p=-x0*w*w/nu*np.exp(-1.5*H*t)*np.sin(nu*t)
    D=np.exp(3*H*t)*(p*p+w*w*x*x)/2
    I=D+1.5*H*np.exp(3*H*t)*x*p
    assert np.max(abs(I/(w*w*x0*x0/2)-1))<2e-15
    # Density oscillates at finite frequency even though its corrected action is constant.
    assert np.ptp(D)>0.01*(w*w*x0*x0/2)


def test_unaveraged_gravity_and_independent_clock_predict_the_same_lightcone():
    p=solve_parent(clock='proper',rtol=2e-11)
    e=solve_parent(clock='efolds',rtol=2e-12)
    for key in ('Lambda_DL_squared','energy_equality_scale','final_comoving_energy_ratio'):
        assert abs(p[key]-e[key])<2e-8
    assert p['independent_photon_integral_gap']<1e-11
    assert max(p['maximum_Friedmann_residual'],e['maximum_Friedmann_residual'])<3e-9
    assert max(p['maximum_energy_work_residual'],e['maximum_energy_work_residual'])<3e-9
