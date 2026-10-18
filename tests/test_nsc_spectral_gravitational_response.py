"""Independent heat/profile controls; no physical archive application."""
import math
import time

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad
from scipy.special import exp1

from recursive_horizons.nsc_spectral_gravitational_response import (
    flat_spin2_profile_response, profile_local_coefficients,
)


ACTIVE_CPU = [0.0]


@pytest.fixture(autouse=True)
def bounded():
    start = time.process_time()
    yield
    ACTIVE_CPU[0] += time.process_time()-start
    assert ACTIVE_CPU[0] < 60, "spectral response controls exhausted their allocation"


def heat_response(z):
    # Direct proper-time heat formula; cancellation-free small-z integrand.
    return quad(lambda a: math.expm1(-a*(1-a)*z)
                +z/4*math.exp(-a*(1-a)*z), 0, 1,
                epsabs=2e-13, epsrel=2e-13)[0]


def test_exchange_of_integrals_and_owned_local_coefficients():
    t, x, u = sp.symbols('t x u', positive=True)
    reach = sp.sqrt(1-4*u/x)
    weight = sp.integrate(x*t*t/8+u/2, (t, -reach, reach))
    assert sp.simplify(weight-(x/12)*(1+8*u/x)*reach) == 0
    v = sp.symbols('v', nonnegative=True)
    q0, q1, slope = sp.symbols('Q0 Q1 slope', positive=True)
    integral = sp.integrate((1+2*v)*sp.sqrt(1-v)*(q0+slope*x*v/4), (v, 0, 1))
    P = sp.expand(x*(1-x*integral/(4*q1)))
    assert sp.simplify(P-(x-3*q0*x*x/(10*q1)-slope*x**3/(28*q1))) == 0
    A, CW = q1/(6*(4*sp.pi)**2), -q0/(40*(4*sp.pi)**2)
    assert sp.simplify(2*CW/A+3*q0/(10*q1)) == 0
    # The existing half-trace coefficient is A/4, not a wave normalization.
    assert sp.simplify((q1/12)/(2*(4*sp.pi)**2)-A/4) == 0


@pytest.mark.parametrize('terms,x', [([(1.0, 0.7)], 1.2), ([(1.0, 0.7), (-0.2, 1.5)], 0.4)])
def test_exact_heat_superposition_including_signed_weights(terms, x):
    profile = lambda y: sum(a*math.exp(-t*y) for a,t in terms)
    q1 = sum(a/t for a,t in terms)
    result = flat_spin2_profile_response(x, profile, q1, quadrature_order=20)
    independent = 12/q1*sum(a/t**2*heat_response(t*x) for a,t in terms)
    assert result['response'] == pytest.approx(independent, rel=2e-12, abs=2e-13)


def test_finite_proper_time_profile_against_independent_double_integral():
    upper, lower = 3.0, 1.1
    profile = lambda y: float(exp1(y/upper**2)-exp1(y/lower**2))
    q1 = upper**2-lower**2
    x = 0.6
    direct = quad(lambda t: heat_response(t*x)/t**3, 1/upper**2, 1/lower**2,
                  epsabs=1e-12, epsrel=1e-12)[0]*12/q1
    result = flat_spin2_profile_response(x, profile, q1, quadrature_order=24)
    assert result['response'] == pytest.approx(direct, rel=3e-12, abs=3e-13)
    coeff = profile_local_coefficients(2*math.log(upper/lower), q1)
    assert 2*coeff['C_W']/coeff['A'] == pytest.approx(coeff['local_quartic_per_quadratic'])


def test_units_zero_momentum_and_local_error_order():
    profile = lambda y: math.exp(-0.8*y)
    q1 = 1/0.8
    scale = 2.3
    original = flat_spin2_profile_response(0.3, profile, q1)
    changed = flat_spin2_profile_response(0.3/scale**2,
                 lambda y: profile(scale**2*y), q1/scale**2)
    assert changed['response'] == pytest.approx(original['response']/scale**2, rel=1e-13)
    assert changed['response_over_Einstein'] == pytest.approx(original['response_over_Einstein'], abs=1e-14)
    assert flat_spin2_profile_response(0, lambda _: pytest.fail('zero evaluates no profile'), q1)['response']==0
    errors=[]
    for x in (0.02,0.01):
        exact=flat_spin2_profile_response(x, profile, q1)['response']
        local=x-3*x*x/(10*q1)
        errors.append(exact-local)
    assert errors[0]/errors[1] == pytest.approx(8, rel=0.003)


def test_invalid_domains_and_nonfinite_profile():
    for x,q1,order in [(-1,1,8),(float('nan'),1,8),(1,0,8),(1,1,True),(1,1,1)]:
        with pytest.raises(ValueError):flat_spin2_profile_response(x,lambda y:1,q1,quadrature_order=order)
    with pytest.raises(ValueError):flat_spin2_profile_response(1,lambda y:float('nan'),1)
    for q0,q1 in [(float('inf'),1),(1,0)]:
        with pytest.raises(ValueError):profile_local_coefficients(q0,q1)


def test_application_preserves_cutoff_couplings_and_independent_ruler():
    import importlib.util
    from pathlib import Path
    path=Path(__file__).resolve().parents[1]/'scripts/measure_nsc_spectral_gravitational_response.py'
    spec=importlib.util.spec_from_file_location('spectral_gravity_measure_control',path)
    owner=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(owner)
    binding={'development_branch':{'coefficients':{'A_Wilsonian_at_magnetic_scale':0.1,
        'C_Weyl':-0.002,'V_full_relational':0},'omega':4.0,'magnetic_flux':4,
        'mmp_map':{'r_e_squared':1.0}}}
    locked={'locked_inputs':{'A':0.1,'C_Weyl':-0.002,'Omega':4.0,'magnetic_flux':4}}
    reference={'domain':{'reference_length':3.0,'target_observation_used':False},
               'parameters':{'A':0.1,'C_W':-0.002}}
    assert owner.select_inputs(binding,locked,reference)[-1]==3.0
    reference['domain']['target_observation_used']=True
    with pytest.raises(ValueError):owner.select_inputs(binding,locked,reference)
    reference['domain']['target_observation_used']=False
    locked['locked_inputs']['C_Weyl']=-0.001
    with pytest.raises(ValueError):owner.select_inputs(binding,locked,reference)
