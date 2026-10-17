"""Symbolic controls of the conditional derivative closure, not existence."""
import pytest
import sympy as s

from control_budget import ControlBudget

CONTROLS=ControlBudget(maximum_evaluations=30000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('derivative closure controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_exterior_quotients_and_triangular_transverse_system(budget):
    budget.check(8000)
    r,X,pr,rx,P,Q,C,pL,Lx,rxx,Qx,al=s.symbols('r X pr rx P Q C pL Lx rxx Qx alpha')
    a,b,c=s.symbols('a b c')
    K=(pr*rx+C/4)/r**2;f=X/r
    fU=(P-f*pr)/r;fx=(Q-f*rx)/r
    def derivative(expr,jets):return sum(s.diff(expr,key)*value for key,value in jets.items())
    dx={r:rx,X:Q,pr:-r*K,P:-K*X,C:C*Lx,rx:rxx,Q:Qx}
    du={r:pr,X:P,pr:a,P:c,C:C*pL,rx:-r*K,Q:-K*X,pL:b}
    cross=-(fU*rx+fx*pr)/r
    assert s.cancel(derivative(fU,dx)-cross)==0
    assert s.cancel(derivative(fx,du)-cross)==0
    assert s.cancel(derivative(fx,dx)-(Qx-2*rx*fx-f*rxx)/r)==0
    assert s.cancel(derivative(fU,du)-(c-2*pr*fU-f*a)/r)==0
    Kx=(-r*K*rx+pr*rxx+C*Lx/4)/r**2-2*K*rx/r
    FK=(-pr*r*K+C*pL/4)/r**2-2*K*pr/r
    KU=rx*a/r**2+FK
    assert s.cancel(derivative(K,dx)-Kx)==0
    assert s.cancel(derivative(K,du)-KU)==0
    targets=(-pr*K-r*KU,2*KU-al*((c-2*pr*fU-f*a)*fx/r+fU*cross),-X*KU-K*P)
    originals=(-r*K,2*K-al*fU*fx,-K*X)
    for original,target in zip(originals,targets):
        assert s.cancel(derivative(original,du)-target)==0
        assert all(s.diff(target,i,j)==0 for i in (a,b,c) for j in (a,b,c))
    assert not targets[0].has(b,c) and not targets[2].has(b,c) and not targets[1].has(b)


def test_cap_commutator_and_primitive_compatibility(budget):
    budget.check(8000)
    z,T=s.symbols('z T')
    beta=s.Function('beta')(z,T);Ep=s.Function('Ep')(z,T);Em=s.Function('Em')(z,T)
    A0=s.Function('A0')(z);S=s.Function('S')(z,T)
    epT=beta*s.diff(Ep,z)+s.diff(beta,z)*Ep-A0*S
    emT=-beta*s.diff(Em,z)-s.diff(beta,z)*Em-A0*S
    for e,rhs,sign in ((Ep,epT,1),(Em,emT,-1)):
        target=sign*beta*s.diff(e,z,2)+sign*2*s.diff(beta,z)*s.diff(e,z)+sign*s.diff(beta,z,2)*e-s.diff(A0*S,z)
        assert s.expand(s.diff(rhs,z)-target)==0
    # chi_z=(Ep-Em)/2 propagates from compatible data, including born points.
    chiT=beta*(Ep+Em)/2
    assert s.expand(s.diff(chiT,z)-(epT-emT)/2)==0


def test_radial_metric_time_derivative_and_beta_second_jet(budget):
    budget.check(8000)
    r,N,mu,f,al,ac,Nr,fz=s.symbols('r N mu f alpha ac Nr fz')
    A0=1-r*r;A=A0-2*al*N/r
    potential=s.Function('P')(f)
    D=2*N/r**2-r*potential-4*N/A0
    def dz(expr):
        return A0*s.diff(expr,r)+A0*Nr*s.diff(expr,N)+fz*s.diff(expr,f)+mu*A0*ac*s.diff(expr,mu)
    Dz=2*A0*Nr/r**2-4*A0*N/r**3-A0*potential-r*s.diff(potential,f)*fz-4*Nr-8*r*N/A0
    assert s.cancel(dz(al*mu*D)-al*mu*(A0*ac*D+Dz))==0
    energy,energyT,fT=s.symbols('energy energyT fT')
    a_c=al*energy/(2*r);b_c=r*r*potential/2+A0*energy/4
    def dt(expr):return s.diff(expr,energy)*energyT+s.diff(expr,f)*fT
    target=r*r*s.diff(potential,f)*fT/2+A*energyT/4
    assert s.expand(dt(b_c)-N*dt(a_c)-target)==0
