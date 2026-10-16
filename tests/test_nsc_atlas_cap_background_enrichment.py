"""Exact moving-knot obstruction to smooth tangent-enriched cap evolution.

The manufactured C2 profile is F(u)=(u-kappa)_+^3. Its second time
tangent is continuous, but its time derivative is not an H1 variation.
These controls do not change the physical source or evolve a cap.
"""
import pytest
import sympy as s

from control_budget import ControlBudget


CONTROLS = ControlBudget(maximum_evaluations=50_000, cpu_seconds=15.)


@pytest.fixture(scope='module')
def data():
    with CONTROLS.active() as budget:
        budget.check(12_000)
        x, L = s.symbols('x L', positive=True)
        q, beta = s.symbols('q beta', real=True)
        # x=z-ell(T), ell_T=-1, and L=Z-ell(T) at a fixed wall Z.
        chi, k = x**3+6*q*x, 3*x**2
        density = s.expand((s.diff(chi,x)**2+k**2)/4)
        chiW = chi.subs(x,L)
        boundary = -beta*chiW**2/4
        mass = s.integrate(density,(x,0,L))+boundary
        direct = s.diff(mass,L)
        classical = s.integrate(s.diff(density,x),(x,0,L))+s.diff(boundary,L)
        yield locals()
    print('moving tangent enrichment:',CONTROLS.budget.evaluations,
          'charged algebra values; zero field samples;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_full_cap_mass_has_a_moving_internal_interface_term(data):
    a=data
    L,q,beta=a['L'],a['q'],a['beta']
    expected=s.Rational(9,10)*L**5+3*q*L**3+9*q*q*L-beta*(L**3+6*q*L)**2/4
    assert s.expand(a['mass']-expected)==0
    # Reynolds: -ell_T [density] = +[density], including its sign.
    jump=a['density'].subs(a['x'],0)
    assert jump==9*q*q
    assert s.expand(a['direct']-a['classical']-jump)==0
    assert s.expand((a['direct']-a['classical']).subs(q,0))==0


def test_time_translation_of_the_added_tangent_is_not_h1_differentiable(data):
    q=data['q']
    epsilon,ell,z=s.symbols('epsilon ell z',positive=True)
    # For chi=6q(z-ell)_+, the spatial-gradient difference is 6q on
    # (ell-epsilon,ell), and zero elsewhere. Its quotient diverges in H1.
    quotient_norm_squared=s.integrate((6*q/epsilon)**2,(z,ell-epsilon,ell))
    assert s.simplify(quotient_norm_squared-36*q*q/epsilon)==0
    assert s.limit(quotient_norm_squared.subs(q,1),epsilon,0,dir='+')==s.oo


def test_third_time_tangent_leaves_the_frozen_three_function_space(data):
    z=s.symbols('z',real=True)
    x=z-s.Rational(1,2)
    # Z=1, ell=1/2; split exactly at the known kink, with no sampling.
    left=s.Matrix([z,0,0])
    right=s.Matrix([z,3*x*x,6*x])
    gram=s.integrate(left*left.T,(z,0,s.Rational(1,2)))
    gram+=s.integrate(right*right.T,(z,s.Rational(1,2),1))
    moment=s.integrate(6*right,(z,s.Rational(1,2),1))
    coefficients=gram.inv()*moment
    direct=s.integrate((left.dot(coefficients))**2,(z,0,s.Rational(1,2)))
    direct+=s.integrate((6-right.dot(coefficients))**2,(z,s.Rational(1,2),1))
    complement=18-(moment.T*coefficients)[0]
    assert gram.det()>0
    assert gram*coefficients==moment
    assert s.simplify(direct-complement)==0
    assert complement>0
    print('third tangent squared L2 projection error:',complement)
