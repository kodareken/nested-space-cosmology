"""Exact basis identities and fixed high-precision arithmetic controls."""
import mpmath as mp
import numpy as np
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_polynomial_basis as model
from recursive_horizons.nsc_atlas_cap_geometric_reference import _positive_gap


CONTROLS = ControlBudget(maximum_evaluations=100000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def receipt():
    yield
    print('odd Legendre controls:',CONTROLS.budget.evaluations,'charged values;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_matrices_and_traces_against_exact_integrals(budget):
    row = model.odd_legendre_operators(3,.5,budget=budget)
    x = s.symbols('x',real=True)
    basis = s.Matrix([s.legendre(n,x) for n in [1,3,5]])
    budget.check(500)
    M = s.integrate(basis*basis.T,(x,0,1))/2
    K = 2*s.integrate(basis.diff(x)*basis.diff(x).T,(x,0,1))
    np.testing.assert_allclose(row['gram'],np.array(M,dtype=float),atol=0,rtol=1e-15)
    np.testing.assert_allclose(row['stiffness'],np.array(K,dtype=float),atol=0,rtol=1e-15)
    E = s.Matrix(row['dilation'].astype(int))
    assert s.simplify(x*basis.diff(x)-E.T*basis) == s.zeros(3,1)
    e = s.ones(3,1)
    t = basis.diff(x).subs(x,1)
    assert E.T*M+M*E-s.Rational(1,2)*e*e.T+M == s.zeros(3)
    assert E.T*K+K*E-K-2*t*t.T == s.zeros(3)
    centre = basis.diff(x).subs(x,0)
    assert e.T*E == t.T
    assert centre.T*E == centre.T
    np.testing.assert_array_equal(row['wall_value'],np.ones(3))
    np.testing.assert_allclose(row['wall_derivative_xi'],np.array(t,dtype=float).ravel(),rtol=0,atol=0)
    np.testing.assert_allclose(row['centre_derivative_xi'],np.array(basis.diff(x).subs(x,0),dtype=float).ravel(),rtol=0,atol=0)


def test_derivative_bound_follows_from_legendre_energy_and_expansion(budget):
    budget.check(500)
    x,n = s.symbols('x n',positive=True)
    P = s.Function('P')(x)
    energy = (1-x*x)*s.diff(P,x)**2+n*(n+1)*P**2
    ode_second = (2*x*s.diff(P,x)-n*(n+1)*P)/(1-x*x)
    assert s.simplify(s.diff(energy,x).subs(s.diff(P,x,2),ode_second)-2*x*s.diff(P,x)**2) == 0
    # At x=1 energy=n(n+1), so the nonnegative energy derivative bounds |P|.
    for degree in [1,3,5,31]:
        expansion = sum((2*j+1)*s.legendre(j,x) for j in range(0,degree,2))
        assert s.expand(s.diff(s.legendre(degree,x),x)-expansion) == 0
        assert sum(2*j+1 for j in range(0,degree,2)) == degree*(degree+1)//2


def test_sixteen_mode_fixed_points_match_high_precision_and_regular_normal_basis(budget):
    maximum = 0.
    with mp.workdps(85):
        for point in [0.,1e-12,.5,.91,1.]:
            row = model.odd_legendre_values(16,point,budget=budget)
            budget.check(16*16)
            x = mp.mpf(point)
            for j,degree in enumerate(range(1,32,2)):
                P = mp.legendre(degree,x)
                derivative = mp.diff(lambda value:mp.legendre(degree,value),x)
                Q = P/x if point else derivative
                Qt = (derivative-Q)/(2*x*x) if point else mp.diff(lambda value:mp.legendre(degree,value),x,3)/6
                for name,expected in [('value',P),('derivative_xi',derivative),('quotient',Q),('quotient_derivative_t',Qt)]:
                    error = abs(row[name][j]-float(expected))
                    maximum = max(maximum,error)
                    assert error <= 2e-11*max(1.,abs(float(expected)))
                if point:
                    Z = .5
                    z,r,A0 = Z*point,np.tanh(Z*point),1/np.cosh(Z*point)**2
                    numerical_h = (2*point*point*row['quotient_derivative_t'][j]
                                   +z*_positive_gap(z)*row['quotient'][j])/(Z*A0)
                    exact_h = derivative/(mp.mpf(Z)/mp.cosh(mp.mpf(Z)*x)**2)-P/mp.tanh(mp.mpf(Z)*x)
                    assert abs(numerical_h-float(exact_h)) <= 2e-11*max(abs(float(exact_h)),1e-30)
    print('sixteen-mode largest raw basis discrepancy:',maximum)


def test_diagonal_gram_has_the_declared_condition_and_dilation_matrix(budget):
    row = model.odd_legendre_operators(16,.5,budget=budget)
    M,E = row['gram'],row['dilation']
    assert row['exact_gram_condition'] == 21.
    assert np.linalg.cond(M) == pytest.approx(21.,abs=1e-13)
    np.testing.assert_allclose(E.T@M+M@E,.5*np.ones((16,16))-M,rtol=0,atol=2e-16)


@pytest.mark.parametrize('count',[0,17,True,1.5])
def test_invalid_mode_count(count,budget):
    with pytest.raises(ValueError):
        model.odd_legendre_operators(count,.5,budget=budget)


@pytest.mark.parametrize('point',[-.1,1.1,np.nan])
def test_invalid_point(point,budget):
    with pytest.raises(ValueError):
        model.odd_legendre_values(3,point,budget=budget)


def test_unrepresentable_radius_and_shared_budget(budget):
    with pytest.raises(ValueError):
        model.odd_legendre_operators(3,0.,budget=budget)
    maximum = budget.maximum_evaluations
    budget.maximum_evaluations = budget.evaluations+1
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            model.odd_legendre_values(3,.5,budget=budget)
    finally:
        budget.maximum_evaluations = maximum
