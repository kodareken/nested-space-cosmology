"""Exact frozen-forcing primitive and wall controls; no existence assertion."""
import pytest
import sympy as s

from control_budget import ControlBudget

CONTROLS=ControlBudget(maximum_evaluations=36000,cpu_seconds=30.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('null construction controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_fixed_source_and_wall_birth_preserve_one_primitive(budget):
    budget.check(10000)
    u,v,t,a=s.symbols('u v t a')
    force=u+v
    # Source u=0, quiet values at v=0, wall u=1+v.
    # The transverse derivative initially vanishes; q(u,0)=u^2/2.
    p_wall=v+v*v
    old=s.integrate(force.subs(v,t),(t,0,v))
    entry=u-1
    born=p_wall.subs(v,entry)+s.integrate(force.subs(v,t),(t,entry,v))
    q=s.integrate(force.subs(u,a),(a,0,u))
    y_old=s.integrate(old.subs(u,a),(a,0,u))
    y_born=s.integrate(old.subs(u,a),(a,0,1))+s.integrate(born.subs(u,a),(a,1,u))
    for p,y in ((old,y_old),(born,y_born)):
        assert s.expand(s.diff(p,v)-force)==0
        assert s.expand(s.diff(y,u)-p)==0
        assert s.expand(s.diff(y,v)-q)==0
    assert s.expand(s.diff(q,u)-force)==0
    assert s.expand((old-born).subs(u,1))==0
    assert s.expand(s.diff(old-born,u).subs(u,1))==0
    wall_value=y_born.subs(u,1+v)
    assert s.expand(s.diff(wall_value,v)-p_wall-q.subs(u,1+v))==0


def test_spacelike_initial_curve_and_wall_birth_preserve_one_primitive(budget):
    budget.check(10000)
    u,v,t,a=s.symbols('u v t a')
    force=u+v;lower=1-v;wall=3+v
    # Initial curve u+v=1 meets the wall at (2,-1).
    # y0'=lower'*p0+q0 with p0(u)=u^2 and q0(v)=v^2.
    y0=v*v-v;q0=v*v
    assert s.expand(s.diff(y0,v)-s.diff(lower,v)*lower**2-q0)==0
    p_wall=4+6*(v+1)+(v+1)**2
    old=u*u+s.integrate(force.subs(v,t),(t,1-u,v))
    entry=u-3
    born=p_wall.subs(v,entry)+s.integrate(force.subs(v,t),(t,entry,v))
    q=q0+s.integrate(force.subs(u,a),(a,lower,u))
    y_old=y0+s.integrate(old.subs(u,a),(a,lower,u))
    y_born=y0+s.integrate(old.subs(u,a),(a,lower,2))+s.integrate(born.subs(u,a),(a,2,u))
    for p,y in ((old,y_old),(born,y_born)):
        assert s.expand(s.diff(p,v)-force)==0
        assert s.expand(s.diff(y,u)-p)==0
        assert s.expand(s.diff(y,v)-q)==0
    assert s.expand((old-born).subs(u,2))==0
    assert s.expand(s.diff(old-born,u).subs(u,2))==0
    assert s.expand(s.diff(y_born.subs(u,wall),v)-p_wall-q.subs(u,wall))==0


def test_two_sided_incoming_wall_jets_transport_the_shared_values(budget):
    budget.check(10000)
    Up,Vp,Qm,Qp,r,Sw=s.symbols('Up Vp Qm Qp r Sw',nonzero=True)
    rp,rm,cp,cm,Lp,Lm,Gp,Gm=s.symbols('rp rm cp cm Lp Lm Gp Gm')
    Rp,Rm=2*Vp*rp,2*Qm*rm
    Cp,Cm=2*Vp*cp,2*Qm*cm
    velocity=(Rp+Rm+Sw*r)/2;scalar_rate=(Cp+Cm)/2
    incoming_plus=((Rm+Sw*r)/(2*Up),Cm/(2*Up),(Gp-Vp*Lp)/Up)
    incoming_minus=((Rp+Sw*r)/(2*Qp),Cp/(2*Qp),(Gm-Qm*Lm)/Qp)
    for incoming,outgoing,a,b,targets in (
        (incoming_plus,(rp,cp,Lp),Up,Vp,(velocity,scalar_rate,Gp)),
        (incoming_minus,(rm,cm,Lm),Qp,Qm,(velocity,scalar_rate,Gm))):
        for pin,qout,target in zip(incoming,outgoing,targets):
            assert s.expand(a*pin+b*qout-target)==0


def test_uniform_lipschitz_control_does_not_imply_derivative_convergence(budget):
    budget.check(4000)
    t=s.symbols('t',real=True);n=s.symbols('n',positive=True,integer=True)
    k=s.symbols('k',nonnegative=True,integer=True)
    f=(1-s.cos(n*t))/n
    assert f.subs(t,0)==0 and s.diff(f,t).subs(t,0)==0
    assert s.diff(f,t)==s.sin(n*t)  # Lipschitz bound 1 for every n.
    # |f|<=2/n -> 0, including its spatial L1 multiple. Its derivative stays 1
    # at a fixed time along n=4k+1. A small weak error cannot bound that trace.
    assert s.limit(2/(4*k+1),k,s.oo)==0
    assert s.simplify(s.diff(f,t).subs({n:4*k+1,t:s.pi/2}))==1
