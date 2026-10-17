"""Independent radial identities and potential controls, without source data."""
from fractions import Fraction as Q

from flint import ctx
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_departure_bound as m

CONTROLS=ControlBudget(maximum_evaluations=12000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('cap departure controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def quiet():
    return dict(r_min=Q(1,4),r_max=Q(1,2),A0_min=Q(3,4),A0_max=Q(15,16),
        length_max=Q(1,3),chi_max=0,E_plus_max=0,E_minus_max=0,N_max=0,D_max=0,
        mu_max=1,potential_max=0,potential_remainder=0,C_Z=0,C_chi=0,C_E_plus=0,C_E_minus=0)


def test_quiet_source_and_departures(budget):
    result=m.radial_comparison(quiet(),budget=budget)
    assert all(result[k]==0 for k in ('N','D','mu_second_order','A_vacuum_subtracted_second_order'))


def test_attenuation_survives_fixed_field_comparison(budget):
    v=quiet();v.update(E_plus_max=Q(3,2),N_max=2)
    result=m.radial_comparison(v,budget=budget)
    # p=h=1; Qc<=32; J<=5/16. Thus J*r*Qc*N/2=5.
    assert result['N']==result['attenuation_contribution']==5
    assert result['D']==0
    assert result['A_vacuum_subtracted_second_order']==40


@pytest.mark.parametrize('change,match',[
    ({'r_min':0},'positive common'),({'C_Z':-1},'nonnegative'),
    ({'length_max':.3},'exact rational'),({'chi_max':1},'zero-alpha')])
def test_domain_and_missing_potential_limit_stop(change,match,budget):
    v=quiet();v.update(change)
    with pytest.raises(ValueError,match=match):m.radial_comparison(v,budget=budget)


def test_normalized_radial_mass_equation(budget):
    budget.check(1024)
    r,alpha,V,Qc=s.symbols('r alpha V Qc');N=s.Function('N')(r)
    mass=r**3/2+alpha*N;A=1-2*mass/r;potential=3+alpha*V
    residual=s.diff(mass,r)-r*r*potential/2-alpha*r*r*A*Qc/4
    normalized=s.diff(N,r)+alpha*r*Qc*N/2-r*r*V/2-r*r*(1-r*r)*Qc/4
    assert s.expand(residual-alpha*normalized)==0


def test_potential_multiple_angle_and_quartic_coefficient(budget):
    budget.check(1024)
    q=s.symbols('q');expr=8*s.cosh(q)**4-s.cosh(4*q)-4*s.cosh(2*q)-3
    assert s.simplify(s.expand(expr.rewrite(s.exp)))==0
    assert s.expand(s.series(3*s.cosh(q)**4-3,q,0,6).removeO())==6*q*q+5*q**4
    with ctx.workprec(160):
        upper=m.potential_remainder_coefficient(m._ball(Q(1)),m._ball(Q(0)))
        assert Q(5,36)<=m.dyadic(upper)<Q(5,36)+Q(1,2**150)


def test_collar_jacobian_in_the_exact_radial_integral(budget):
    budget.check(1024)
    eta,L,z0=s.symbols('eta L z0');r=s.tanh(z0+eta*L)
    # A primitive with nonconstant density verifies the changing measure.
    N=r**3/3-s.tanh(z0)**3/3
    assert s.simplify(s.diff(N,eta)-r*r*L*(1-r*r))==0
    assert s.simplify(s.diff(s.diff(N,eta),L)-s.diff(r*r*L*(1-r*r),L))==0
    assert s.simplify(s.diff(N,eta)-r*r)!=0


def test_budget_stops_before_input_access(budget):
    class Spent:
        def check(self,count=0):raise RuntimeError('spent')
    with pytest.raises(RuntimeError,match='spent'):
        m.cap_departure_bound({}, {},budget=Spent())
