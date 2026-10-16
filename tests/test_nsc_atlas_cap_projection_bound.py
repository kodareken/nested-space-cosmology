"""Exact projection constants and independently enclosed manufactured norms."""
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_projection_bound as model
from recursive_horizons import nsc_atlas_cap_variational_reference as radial

CONTROLS = ControlBudget(maximum_evaluations=100000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('cap projection-bound controls:',CONTROLS.budget.evaluations,
          'coefficient values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def inputs(**changed):
    data = dict(z_nodes=np.array([0.,.07,.21,.5]),chi=np.array([0.,.03,.12,.13]),
                canonical_density=np.array([0.,-.02,.07,.08]),alpha=.002)
    data.update(changed)
    return data


def dyadic(row):
    return Fraction(int(row['mantissa']))*Fraction(2)**row['exponent']


def test_exact_interpolation_green_kernel_and_eighth_constant(budget):
    budget.check(16)
    h,x,t = s.symbols('h x t',positive=True)
    # 0<=x<=h: the absolute Green kernel is integrated on its two pieces.
    area = s.integrate(t*(h-x)/h,(t,0,x))+s.integrate(x*(h-t)/h,(t,x,h))
    assert s.expand(area-x*(h-x)/2) == 0
    assert s.expand(h*h/8-area-(x-h/2)**2/2) == 0


def test_centre_reduced_mass_matrix_has_the_required_dominance_margin(budget):
    budget.check(12)
    h0,h1 = s.symbols('h0 h1',positive=True)
    interior = (h0+h1)/3-(h0+h1)/6-(h0+h1)/6
    first = (h0+h1)/3-h1/6-(h0+h1)/6
    wall = h0/3-h0/6-h0/6
    single = h0/3-h0/6
    assert s.expand(interior) == 0
    assert s.expand(first) == h0/6
    assert s.expand(wall) == 0
    assert single == h0/6


def test_constrained_spatial_derivatives_of_metric_ratio(budget):
    budget.check(32)
    r,alpha = s.symbols('r alpha',positive=True)
    n,f,p,q,delta = (s.Function(k)(r) for k in ['n','f','p','q','delta'])
    P = s.Function('P')(f)
    U = s.diff(P,f)
    A0 = 1-r*r
    A = A0-2*alpha*n/r
    a = alpha*(p*p+(q-f)**2)/(2*r)
    rho = r*r*P/2+A*(p*p+(q-f)**2)/4
    beta = A*s.exp(delta)/A0
    C = 2*n/r**2-r*P-4*n/A0
    rules = {s.diff(n,r):rho,s.diff(delta,r):a,s.diff(f,r):(q-f)/r}
    first = (A0*s.diff(beta,r)).subs(rules)/s.exp(delta)-alpha*C
    assert s.cancel(s.expand(first)) == 0
    Cr = 2*rho/r**2-4*n/r**3-P-r*U*(q-f)/r-4*rho/A0-8*r*n/A0**2
    assert s.cancel(s.expand(s.diff(C,r).subs(rules)-Cr)) == 0
    second = (A0*s.diff(alpha*s.exp(delta)*C,r)).subs(rules)
    assert s.cancel(s.expand(second/s.exp(delta)-alpha*A0*(a*C+Cr))) == 0
    assert s.cancel(s.expand(C-(2*n/r**2-r*P-2*n/A0))) != 0


def test_manufactured_bound_encloses_independent_high_precision_formula(budget,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('an algebraic bound called a radial ODE')
    monkeypatch.setattr(radial,'solve_ivp',forbidden)
    data = inputs()
    precision = model.ctx.prec
    row = model.cap_momentum_projection_bound(**data,budget=budget)
    assert model.ctx.prec == precision
    budget.check(256)
    with mp.workdps(85):
        z,c,k = [[mp.mpf(float(v)) for v in data[name]]
                 for name in ['z_nodes','chi','canonical_density']]
        alpha,Z = mp.mpf(data['alpha']),z[-1]
        R,A0 = mp.tanh(Z),1/mp.cosh(Z)**2
        widths = [b-a for a,b in zip(z[:-1],z[1:])]
        lc = max(abs((b-a)/h) for a,b,h in zip(c[:-1],c[1:],widths))
        lk = max(abs((b-a)/h) for a,b,h in zip(k[:-1],k[1:],widths))
        F = lc*Z/R
        argument = mp.sqrt(alpha/6)*F
        sinc = mp.sinh(argument)/argument
        P = F*F*sinc*sinc*(mp.cosh(argument)**2+1)/2
        U = 2*F*sinc*mp.cosh(argument)**3
        hmax = lc/A0+F
        pi = lk*Z/(R*A0)
        phi = max(lc*(R/A0-Z)/(R*R),hmax/mp.tanh(z[1]))
        local_phi = []
        for j,width in enumerate(widths):
            slope = (c[j+1]-c[j])/width
            intercept = c[j]-slope*z[j]
            hi,ahi = mp.tanh(z[j+1]),1/mp.cosh(z[j+1])**2
            bound = abs(slope)*(hi/ahi-z[j+1])/(hi*hi)
            if j:
                bound += abs(intercept)/mp.tanh(z[j])**2
            local_phi.append(bound)
        phi = min(phi,max(local_phi))
        square = pi*pi+phi*phi
        cb,nu = P/2+square/4,(P/2+square/4)/3
        C = R*(2*nu+P+4*nu*R*R/A0)
        Cr = 2*cb+4*nu+P+R*U*phi+4*R*R*cb/A0+8*nu*R**4/A0**2
        delta0 = -mp.log(mp.cosh(mp.sqrt(alpha/6)*c[1]/z[1]))
        delta_upper = delta0+alpha*R*R*square/4
        beta1 = alpha*mp.exp(delta_upper)*C
        beta2 = alpha*mp.exp(delta_upper)*(alpha*R*square*C/2+Cr)
        g2 = beta2*lk*Z+2*beta1*lk
        velocity = 3*max(widths)**2*g2/8
        eta_fraction = dyadic(row['metric_certificate']['A_lower'])
        eta = mp.mpf(eta_fraction.numerator)/eta_fraction.denominator
        momentum = velocity/(eta*mp.exp(delta0))
        for key,exact in [('g_second_derivative_upper',g2),
                          ('wall_velocity_projection_upper',velocity),
                          ('wall_momentum_projection_upper',momentum)]:
            bound = dyadic(row[key])
            value = mp.mpf(bound.numerator)/bound.denominator
            assert value >= exact
            assert value-exact < abs(exact)*mp.mpf('1e-40')
            assert Fraction(row[key]['float_up']) >= bound
    assert not row['includes_gradient_error']
    assert not row['includes_radial_solver_error']
    assert not row['includes_normal_trace_error']
    assert not row['includes_mesh_motion']
    print('manufactured exact-model momentum projection bound:',row['wall_momentum_projection_upper']['float_up'])


def test_affine_element_normal_field_and_monotone_kernel(budget):
    budget.check(16)
    r,t,a,b = s.symbols('r t a b',positive=True)
    A0 = 1-r*r
    field = a*s.atanh(r)+b
    phi = (s.diff(field,r)-field/r)/r
    kernel = (r/A0-s.atanh(r))/(r*r)
    assert s.simplify(phi-a*kernel+b/r**2) == 0
    # K(r)=r*integral_0^1 2*t^2/(1-r^2*t^2)^2 dt. For 0<r<1
    # its derivative has the following strictly positive integrand.
    integrand = 2*r*t*t/(1-r*r*t*t)**2
    derivative = 2*t*t/(1-r*r*t*t)**2+8*r*r*t**4/(1-r*r*t*t)**3
    assert s.simplify(s.diff(integrand,r)-derivative) == 0
    # Fix the real primitive on 0<=r<1 by its derivative and zero value.
    # An unconstrained indefinite integral may choose log(r-1), differing
    # by a complex constant on this physical interval.
    primitive = r/A0-s.atanh(r)
    assert s.simplify(s.diff(primitive,r)-2*r**2/(1-r*r)**2) == 0
    assert s.limit(primitive,r,0) == 0


@pytest.mark.parametrize('changed',[dict(alpha=0.),dict(canonical_density=np.zeros(4))])
def test_projection_vanishes_exactly_in_its_known_limits(changed,budget):
    row = model.cap_momentum_projection_bound(**inputs(**changed),budget=budget)
    assert dyadic(row['wall_velocity_projection_upper']) == 0
    assert dyadic(row['wall_momentum_projection_upper']) == 0
    assert row['wall_momentum_projection_upper']['float_up'] == 0
    assert dyadic(row['B_wall_lower']) > 0


@pytest.mark.parametrize('changed',[
    dict(alpha=100.),dict(alpha=-.1),dict(chi=[1.,.03,.12,.13]),
    dict(z_nodes=[0.,.07,.07,.5]),dict(canonical_density=[0.,np.nan,.07,.08]),
])
def test_invalid_domain_is_rejected(changed,budget):
    with pytest.raises(ValueError):
        model.cap_momentum_projection_bound(**inputs(**changed),budget=budget)


def test_bound_uses_the_shared_cumulative_budget(budget):
    previous = budget.maximum_evaluations
    budget.maximum_evaluations = budget.evaluations+1
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            model.cap_momentum_projection_bound(**inputs(),budget=budget)
    finally:
        budget.maximum_evaluations = previous
