import pytest
import sympy as s
from recursive_horizons.nsc_atlas_boundary_work import work_coefficient


def test_covariant_boundary_work_has_the_required_clock_factors():
    t,eps,nu,V=s.symbols('t eps nu V',positive=True)
    E=s.Function('E')(t);b2=s.Function('b2')(t);D=s.symbols('D',positive=True)
    b=1+eps**2*b2
    density=2*eps**2*E/(V*nu**s.Rational(3,2))
    proper_derivative=lambda x:eps**2/s.sqrt(nu)*s.diff(x,t)
    H=proper_derivative(b)/b
    ward=proper_derivative(density)+3*H*(density+1/s.sqrt(nu))
    coefficient=s.expand(s.series(ward,eps,0,5).removeO()).coeff(eps,4)
    chosen={s.diff(E,t):-D,s.diff(b2,t):2*D/(3*V*nu)}
    assert s.simplify(coefficient.subs(chosen))==0
    assert s.simplify(coefficient.subs({s.diff(E,t):-D,s.diff(b2,t):0}))==-2*D/(V*nu**2)


def test_integrated_scale_law_and_work_sign():
    row=work_coefficient(5.4760189485245965,5.391737694582312,1.6573316678697152)
    assert row['b2']>0
    assert row['db2_dslowtime']>0
    assert row['flat_boundary_Ward_residual_over_epsilon4']+row['work_term_over_epsilon4']==pytest.approx(0)
    assert row['b2']<row['maximum_energy_budget_b2']
    assert work_coefficient(5,5,0)['b2']==0
