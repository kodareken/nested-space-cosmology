"""Exact finite-space boundary onset order; no physical trajectory replay."""
import sympy as s
import pytest
from control_budget import ControlBudget

CONTROLS=ControlBudget(maximum_evaluations=30000,cpu_seconds=15.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('onset lag controls:',CONTROLS.budget.evaluations,
          'charged algebra values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_cubic_port_has_one_order_of_strong_boundary_lag(budget):
    budget.check(18000)
    t,b=s.symbols('t b');Z=s.Rational(1,2);powers=(1,3,5);tx=s.Integer(2)
    M=s.Matrix([[Z**(i+j+1)/s.Integer(i+j+1) for j in powers] for i in powers])
    K=s.Matrix([[i*j*Z**(i+j-1)/s.Integer(i+j-1) for j in powers] for i in powers])
    e=s.Matrix([Z**i for i in powers]);ell=s.Matrix([i*Z**(i-1) for i in powers])
    inv=M.inv();S=K-b*e*e.T
    J=s.zeros(3).row_join(tx*s.eye(3)).col_join((-tx*inv*S).row_join(-tx*inv*e*e.T))
    forcing=s.zeros(3,1).col_join(tx*inv*e)
    # Independent Taylor recurrence from u'=J u+forcing*t^3, u(0)=0.
    u=forcing*t**4/4+J*forcing*t**5/20+J**2*forcing*t**6/120
    residual=u.diff(t)-J*u-forcing*t**3
    for entry in residual:
        assert all(s.expand(entry).coeff(t,j)==0 for j in range(6))
    chi=u[:3,0];density=u[3:,0]
    assert chi.applyfunc(lambda v:s.expand(v).coeff(t,5))==tx**2*inv*e/20
    assert density.applyfunc(lambda v:s.expand(v).coeff(t,4))==tx*inv*e/4
    incoming=s.expand(((ell-b*e).T*chi+e.T*density)[0])
    assert incoming.coeff(t,3)==0
    assert incoming.coeff(t,4)==tx*(e.T*inv*e)[0]/4
    boundary=incoming-t**3
    assert boundary.coeff(t,3)==-1


def test_strong_port_defect_maps_to_actual_moving_scalar_tangent(budget):
    budget.check(6000)
    R,A,tx,rx,chi,chiz,d,C=s.symbols('R A tx rx chi chiz d C',nonzero=True)
    p=d/A;h=chiz/A-chi/R
    incoming=C/A
    matched_p=(incoming+p-h)/2;matched_h=(incoming-p+h)/2
    actual=(tx*d+rx*h)/R
    target=(tx*A*matched_p+rx*matched_h)/R
    boundary=chiz-A*chi/R+d-C
    expected=(rx+tx*A)*boundary/(2*R*A)
    assert s.simplify(actual-target-expected)==0
    assert s.simplify(R*(actual-target)-(rx+tx*A)*boundary/(2*A))==0
