"""Actual source forcing controls, independent of saved trajectory samples."""
from fractions import Fraction as Q

import mpmath as mp
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_source_forcing_bound import _up,_down,source_forcing_bound

CONTROLS=ControlBudget(maximum_evaluations=30000,cpu_seconds=15.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('source forcing controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_source_frame_identities_and_canonical_gap(budget):
    budget.check(8000)
    sj,v,g,gx,al=s.symbols('sj v g gx al',positive=True)
    r=s.sqrt(sj*sj-al*g*g/6);F=1-al*g*g/(6*sj*sj)
    rx=s.diff(r,sj)*v+s.diff(r,g)*gx
    Lx=s.diff(s.log(F),sj)*v+s.diff(s.log(F),g)*gx
    assert s.simplify(r*rx-(sj*v-al*g*gx/6))==0
    assert s.simplify(Lx+al*(g*gx/sj**2-g*g*v/sj**3)/(3*F))==0
    a,u=s.symbols('a u',positive=True)
    j=s.asinh(a*u)/a-u/s.sqrt(1+a*a*u*u)
    assert s.simplify(s.diff(j,u)-a*a*u*u/(1+a*a*u*u)**s.Rational(3,2))==0
    assert s.limit(j,u,0)==0  # integrate 0<=j'<=a^2*u^2 to get |j|<=alpha|u|^3/18


def test_zero_source_and_zero_coupling_limits(budget):
    zero=source_forcing_bound(left=Q(2,5),right=Q(41,100),alpha_max=Q(1,1000),
                              g=0,g_x=0,g_xx=0,budget=budget)
    assert all(v==0 for v in zero['forcing'].values())
    probe=source_forcing_bound(left=Q(2,5),right=Q(41,100),alpha_max=0,
                               g=Q(1,100),g_x=Q(1,2),g_xx=3,budget=budget)
    assert all(probe['forcing'][k]==0 for k in ('r','r_x','log_C','log_C_x','r_times_r_x_minus_x'))
    assert probe['forcing']['X']>=Q(1,100) and probe['forcing']['Q']==Q(1,2)


def test_bound_contains_independent_high_precision_frame_values(budget):
    bound=source_forcing_bound(left=Q(2,5),right=Q(41,100),alpha_max=Q(1,1000),
                               g=Q(1,100),g_x=Q(1,2),g_xx=3,budget=budget)
    budget.check(4000)
    with mp.workdps(80):
        def m(q):q=Q(q);return mp.mpf(q.numerator)/q.denominator
        alpha=m(Q(1,1000));x=m(Q(81,200));g=m(Q(1,100));gx=-m(Q(1,2))
        sj=x-alpha*m(bound['radial']['h']);v=1-alpha*m(bound['radial']['h_x'])
        F=1-alpha*g*g/(6*sj*sj);r=mp.sqrt(sj*sj-alpha*g*g/6)
        rx=(sj*v-alpha*g*gx/6)/r
        Lx=-alpha*(g*gx/(sj*sj)-g*g*v/(sj**3))/(3*F)
        a=mp.sqrt(alpha/6);u=g/r;h=1/mp.sqrt(1+a*a*u*u)
        j=mp.asinh(a*u)/a-u*h
        values=dict(r=abs(r-x),r_x=abs(rx-1),log_C=abs(mp.log(F)),log_C_x=abs(Lx),
                    c=abs(g),c_x=abs(gx),X=abs(r*mp.asinh(a*u)/a),Q=abs(h*gx+j*rx),
                    r_times_r_x_minus_x=abs(r*rx-x))
        assert all(value<=m(bound['forcing'][key])+mp.mpf('1e-75') for key,value in values.items())


def test_fixed_dyadic_rounding_and_unclosed_source_tube(budget):
    budget.check(1000)
    for value in (Q(0),Q(1,3),-Q(1,3),Q(1,10**50)):
        assert _down(value)<=value<=_up(value)
        assert _up(value)-_down(value)<=Q(1,2**128)
    with pytest.raises(ValueError,match='Jordan tube'):
        source_forcing_bound(left=Q(2,5),right=Q(41,100),alpha_max=1,g=1,g_x=100,g_xx=100,budget=budget)
