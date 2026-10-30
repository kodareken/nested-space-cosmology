"""Null geometry, exact source normalization and independent future controls."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_parent_circulation import Parent,horizon_flux,static_patch_energy
from recursive_horizons.nsc_atlas_horizon_response import (
    solve_source,solve_response,initial_area_coefficient,area_record,angular_record)


@pytest.fixture(scope='module')
def evolved():
    source=solve_source()
    return source,solve_response(source)


def test_rescaled_mode_equation_and_its_nonincreasing_energy():
    u,m,l=s.symbols('u m l',real=True);Y=s.Function('Y')(u)
    q=s.cosh(u)**(-s.Rational(3,2))*Y
    transformed=s.diff(q,u,2)+3*s.tanh(u)*s.diff(q,u)+(m*m+l/s.cosh(u)**2)*q
    expected=s.diff(Y,u,2)+(m*m-s.Rational(9,4)+(l+s.Rational(3,4))/s.cosh(u)**2)*Y
    assert s.simplify(s.expand_trig(transformed*s.cosh(u)**s.Rational(3,2)-expected))==0
    w2=m*m-s.Rational(9,4)+(l+s.Rational(3,4))/s.cosh(u)**2
    E=(s.diff(Y,u)**2+w2*Y*Y)/2
    derivative=s.diff(E,u).subs(s.diff(Y,u,2),-w2*Y)
    assert s.simplify(derivative+(l+s.Rational(3,4))*s.tanh(u)*Y*Y/s.cosh(u)**2)==0


def test_null_generator_nonaffinity_and_two_background_expansions():
    t,c,H=s.symbols('t c H',positive=True)
    a=s.cosh(H*t)/H;h=s.diff(a,t)/a
    k=s.Matrix([1,-1/a])
    # Nonzero connection coefficients of the radial -dt^2+a^2 dchi^2 metric.
    acceleration=s.Matrix([a*s.diff(a,t)*k[1]**2,s.diff(k[1],t)+2*h*k[1]])
    assert all(s.simplify(z)==0 for z in acceleration-h*k)
    ch=s.pi/2-s.atan(s.sinh(H*t))
    minus=2*(h-s.cot(ch)/a);plus=2*(h+s.cot(ch)/a)
    assert s.trigsimp(minus)==0
    assert s.trigsimp(plus-4*h)==0


def test_affine_parameter_gives_the_published_dynamical_area_operator():
    t,H=s.symbols('t H',positive=True)
    affine=s.sinh(H*t)/H
    coefficient=s.simplify(affine/s.diff(affine,t))
    assert coefficient==s.tanh(H*t)/H


def test_exact_initial_energy_and_area_coefficient_from_volume_polynomials():
    x,m2=s.symbols('x m2',positive=True)
    first=s.integrate(x*x*(3+(m2-1)*x*x),(x,0,1))
    second=s.integrate(x*x*(3-6*x*x+4*x**4+m2*(x*x-x**4)),(x,0,1))
    assert first==(m2+4)/5
    assert second==(2*m2+13)/35
    from_integrals=2*first/(m2+3)+12*second/(m2+8)
    expected=2*(m2+4)/(5*(m2+3))+12*(2*m2+13)/(35*(m2+8))
    assert s.simplify(from_integrals-expected)==0
    p=Parent.retained();z=p.initial()
    numeric=p.kappa*p.H*static_patch_energy(0,z,p)/(4*math.pi)
    assert numeric==pytest.approx(initial_area_coefficient(p),abs=1e-13)


def test_future_raychaudhuri_solver_against_an_exact_manufactured_source():
    # kappa=1; exact G=sech^2(u), D=1-tanh(u), D+tanh(u)G=1-tanh^3(u).
    def f(u):
        F=12*math.pi*math.tanh(u)/math.cosh(u)**2
        return np.array([math.tanh(u)*F,F/math.cosh(u),0.,0.,0.])
    sol=solve_response(SimpleNamespace(endpoint_u=14.,integrands=f))
    for u in (0.,.3,1.,2.):
        P,S=sol.sol(u)[:2]
        D=(P-math.sinh(u)*S)/(4*math.pi)
        G=math.cosh(u)*S/(4*math.pi)
        assert D==pytest.approx(1-math.tanh(u),abs=2e-11)
        assert G==pytest.approx(1/math.cosh(u)**2,abs=2e-11)
        assert D+math.tanh(u)*G==pytest.approx(1-math.tanh(u)**3,abs=2e-11)


def test_scalar_null_stress_matches_the_previous_parent_interface(evolved):
    source,_=evolved;p=source.parent
    for u in (0.,.25,1.,4.,8.):
        A=source.null_amplitudes(u)
        F=4*math.pi/(3*p.H**2)*np.dot(A,A)
        old,_=horizon_flux(u/p.H,source.state(u),p)
        assert F==pytest.approx(old,rel=2e-13,abs=1e-20)


def test_same_flux_fixes_area_and_the_finite_endpoint_entropy_correction(evolved):
    source,response=evolved;p=source.parent
    for u in (0.,.5,1.1846453242616595,3.):
        r=area_record(source,response,u)
        assert r['comparison_first_law_absolute_residual']<2e-10
        assert r['normalized_area_growth_rate_over_epsilon_squared']>0
        assert r['ordinary_area_deficit_over_epsilon_squared']>=0
    r=area_record(source,response,1.1846453242616595)
    assert r['finite_endpoint_correction_over_epsilon_squared']>0.09
    assert r['ordinary_area_growth_over_epsilon_squared']==pytest.approx(1.05838813260,abs=2e-10)
    assert r['dynamical_area_growth_over_epsilon_squared']==pytest.approx(.96050282298,abs=2e-10)


def test_angular_flux_matrix_and_cap_integrals(evolved):
    source,response=evolved;r=angular_record(source,response)
    Q=np.array(r['matrix']);axis=np.array(r['axis_in_n_coordinates'])
    assert min(np.linalg.eigvalsh(Q))>=-1e-14
    np.testing.assert_allclose(Q@axis,np.linalg.eigvalsh(Q)[-1]*axis,atol=1e-14)
    assert r['total_energy_from_trace']==pytest.approx(static_patch_energy(0,source.state(0),source.parent),rel=3e-11)
    # Independently integrate both caps after diagonalizing Q, with mu=cos(theta).
    nodes,w=np.polynomial.legendre.leggauss(32);c=math.cos(math.pi/6)
    mu=c+(nodes+1)*(1-c)/2
    values=np.linalg.eigvalsh(Q)[::-1]
    cap=2*2*math.pi*np.dot(w*(1-c)/2,values[0]*mu*mu+(values[1]+values[2])*(1-mu*mu)/2)
    total=4*math.pi*np.trace(Q)/3
    assert cap/total==pytest.approx(r['two_30_degree_caps_energy_fraction'],abs=1e-14)


def test_analytic_future_tails_control_exact_source_decay(evolved):
    source,_=evolved;p=source.parent;bound=source.tail_bounds()
    assert bound['normalized_area_deficit_tail_upper_bound']<1e-15
    for u in (2.,4.,6.,8.):
        F=source.integrands(u)[0]/math.tanh(u)
        assert F<=bound['F_sech5_coefficient']/math.cosh(u)**5
    assert source.nu_squared>0


def test_record_source_bindings_and_scope():
    root=Path(__file__).resolve().parents[1]
    path=root/'results/summaries/nsc-atlas-horizon-response-v1.json'
    if not path.exists():pytest.skip('record follows source freeze')
    r=json.loads(path.read_text(encoding='utf-8'))
    for name,digest in r['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
    assert r['no_observational_fit'] is True
    assert r['controls']['comparison_first_law_max_absolute_residual']<1e-9
