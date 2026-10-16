"""Exact frozen-cap energy and time-operator controls, without trajectories."""
import sympy as s
import pytest

from control_budget import ControlBudget

CONTROLS=ControlBudget(maximum_evaluations=80000,cpu_seconds=20.)
Z,beta=s.symbols('Z beta',positive=True)
powers=(1,3,5)
# Independent monomials span the same odd-quintic one-cell space.
M=s.Matrix([[Z**(i+j+1)/s.Integer(i+j+1) for j in powers] for i in powers])
K=s.Matrix([[i*j*Z**(i+j-1)/s.Integer(i+j-1) for j in powers] for i in powers])
e=s.Matrix([Z**i for i in powers])


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('frozen cap time controls:',CONTROLS.budget.evaluations,
          'charged algebra values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def matrices():
    Minv=M.inv();S=K-beta*e*e.T
    J=s.zeros(3).row_join(s.eye(3)).col_join((-Minv*S).row_join(-Minv*e*e.T))
    return Minv,S,J,s.diag(S,M)


def test_endpoint_trace_and_energy_forms_are_positive_in_the_declared_domain(budget):
    budget.check(12000)
    for matrix in (M,K):
        for n in (1,2,3): assert s.factor(matrix[:n,:n].det()).is_positive
    trace=s.simplify(Z*K-e*e.T)
    assert trace[0,:]==s.zeros(1,3) and trace[:,0]==s.zeros(3,1)
    assert s.factor(trace[1,1]).is_positive
    assert s.factor(trace[1:,1:].det()).is_positive
    S=K-beta*e*e.T
    assert s.simplify(S-((1-beta*Z)*K+beta*trace))==s.zeros(3)
    # Thus beta*Z<1 makes S positive. For the physical beta this follows
    # from sinh(2Z)>2Z, since every non-linear Taylor coefficient is positive.
    remainder=s.series(s.sinh(2*Z)-2*Z,Z,0,9).removeO()
    assert all(c>0 for c in s.Poly(remainder,Z).coeffs())


def test_frozen_generator_dissipates_exact_mass_at_the_boundary(budget):
    budget.check(12000)
    _,_,J,H=matrices()
    expected=s.diag(s.zeros(3),-2*e*e.T)
    assert s.simplify(J.T*H+H*J-expected)==s.zeros(6)
    c=s.Matrix(s.symbols('c0:3'));d=s.Matrix(s.symbols('d0:3'));u=c.col_join(d)
    mass_rate=(u.T*(J.T*H+H*J)*u)[0]/4
    assert s.factor(mass_rate+(e.dot(d))**2/2)==0


def test_energy_weighted_operator_bound_is_not_the_gram_condition(budget):
    budget.check(12000)
    Minv,S,J,H=matrices()
    q=(e.T*Minv*e)[0]
    weighted=s.simplify(s.trace(H.inv()*J.T*H*J))
    exact=2*s.trace(Minv*S)+q*q
    bound=2*s.trace(Minv*K)+q*q
    assert s.simplify(weighted-exact)==0
    assert s.simplify(bound-weighted-2*beta*q)==0
    assert s.factor(q).is_positive
    assert s.simplify(bound.subs(Z,s.Rational(1,2)))==3444


def test_cubic_step_and_a_rational_energy_amplification_bound(budget):
    budget.check(24000)
    Minv,_,J,H=matrices()
    # One manufactured admissible Robin value, one step fixed by the
    # exact generator bound. This is not a physical pulse or step scan.
    at={Z:s.Rational(1,2),beta:s.Rational(1)}
    J,H=J.subs(at),H.subs(at)
    L=s.Integer(59)  # L^2=3481 exceeds the exact upper bound 3444.
    h=1/(16*L);identity=s.eye(6)
    first=identity+h*J
    second=3*identity/4+(identity+h*J)*first/4
    third=identity/3+2*(identity+h*J)*second/3
    polynomial=identity+h*J+(h*J)**2/2+(h*J)**3/6
    assert third==polynomial
    x=h*L
    # exp(x)<=1/(1-x), hence this rational delta bounds the exponential
    # Taylor remainder. Test the resulting energy inequality exactly.
    delta=x**4/(24*(1-x))
    remainder=(1+delta)**2*H-polynomial.T*H*polynomial
    for n in range(1,7): assert remainder[:n,:n].det()>0


def test_moving_cut_energy_keeps_one_domain_term_and_the_input_flux(budget):
    budget.check(16000)
    bZ,A0,Tx,Zx,C=s.symbols('bZ A0 Tx Zx C',real=True)
    Minv,S,J,H=matrices()
    ell=e.diff(Z);v=ell-beta*e
    SZ=S.diff(Z)+S.diff(beta)*bZ
    SZ=s.simplify(SZ.subs(bZ,-beta**2-2*A0))
    assert s.simplify(SZ-v*v.T-2*A0*e*e.T)==s.zeros(3)
    assert M.diff(Z)==e*e.T
    HZ=s.diag(SZ,e*e.T)
    gamma=s.trace(S.inv()*SZ)+s.trace(Minv*e*e.T)
    expected=(v.T*S.inv()*v)[0]+2*A0*(e.T*S.inv()*e)[0]+(e.T*Minv*e)[0]
    assert s.simplify(gamma-expected)==0
    c=s.Matrix(s.symbols('a0:3'));d=s.Matrix(s.symbols('d0:3'));u=c.col_join(d)
    forcing=s.zeros(3,1).col_join(Tx*Minv*e*C)
    rate=Tx*J*u+forcing
    mass_rate=(2*(u.T*H*rate)[0]+Zx*(u.T*HZ*u)[0])/4
    V=e.dot(d)
    target=Tx*(V*C-V**2)/2+Zx*(u.T*HZ*u)[0]/4
    assert s.simplify(mass_rate-target)==0
    assert s.expand((V*C-V**2)/2-(C**2/8-(V-C/2)**2/2))==0
    # The geometry term equals A0*rho_W, not an additional M-dot force.
    shape=(u.T*HZ*u)[0]/4
    assert s.simplify(shape-((v.dot(c))**2+V**2)/4-A0*(e.dot(c))**2/2)==0
