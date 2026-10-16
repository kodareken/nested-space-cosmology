"""Fixed physical-coordinate space: exact shape, weak mass and initial data.

These are algebraic controls for the next local-space construction, not a
finite-alpha trajectory or a field-approximation error certificate.
"""
import pytest
import sympy as s

from control_budget import ControlBudget


CONTROLS=ControlBudget(maximum_evaluations=60_000,cpu_seconds=20.)


@pytest.fixture(scope='module')
def data():
    with CONTROLS.active() as budget:
        budget.check(20_000)
        z,Z,lam,A0=s.symbols('z Z lam A0',positive=True)
        beta,C,Zx=s.symbols('beta C Zx',real=True)
        c=s.Matrix(s.symbols('c0 c1',real=True))
        d=s.Matrix(s.symbols('d0 d1',real=True))
        phi=s.Matrix([z,z**3+lam*z**5])
        e=phi.subs(z,Z)
        chi,k=phi.dot(c),phi.dot(d)
        M=s.integrate(phi*phi.T,(z,0,Z))
        wall=chi.subs(z,Z)
        slope=s.diff(chi,z).subs(z,Z)
        mass=s.integrate(s.expand(k*k+s.diff(chi,z)**2),(z,0,Z))/4-beta*wall**2/4
        gc=s.Matrix([s.diff(mass,v) for v in c])
        gd=s.Matrix([s.diff(mass,v) for v in d])
        yield locals()
    print('fixed physical-space identities:',CONTROLS.budget.evaluations,
          'charged algebra values; zero field samples;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_fixed_space_shape_is_boundary_density_without_dilation(data):
    q=data
    Z,A0,beta=q['Z'],q['A0'],q['beta']
    shape=s.diff(q['mass'],Z)-(2*A0+beta**2)*s.diff(q['mass'],beta)
    density=A0*q['wall']**2/2+(q['k'].subs(q['z'],Z)**2+(q['slope']-beta*q['wall'])**2)/4
    assert s.expand(shape-density)==0
    # This basis is intentionally not closed under z*d/dz: equal z^3
    # coefficients force a factor 3, which leaves a nonzero z^5 term.
    z,lam=q['z'],q['lam']
    assert s.expand(z*s.diff(q['phi'][1],z)-3*q['phi'][1]-2*lam*z**5)==0
    assert lam>0


def test_actual_quadratic_mass_and_weak_rates_give_the_same_flux(data):
    q=data
    M,d,gc,gd,e=q['M'],q['d'],q['gc'],q['gd'],q['e']
    assert all(s.expand(v)==0 for v in gd-M*d/2)
    # Positive definiteness is exact for every Z,lambda>0.
    assert M[0,0]>0
    assert s.expand(M.det()).is_positive
    Vw=e.dot(d)
    rhs=-2*gc+e*(q['C']-Vw)  # M*d_T; chi_T coefficients are d.
    direct=gc.dot(d)+d.dot(rhs)/2
    assert s.expand(direct-Vw*(q['C']-Vw)/2)==0


def test_field_coefficients_are_not_moving_domain_moments(data):
    q=data
    MZ=q['M'].diff(q['Z'])
    assert all(s.expand(v)==0 for v in MZ-q['e']*q['e'].T)
    # Incorrectly adding -M^-1*M_Z*d*Z_x to field coefficient rates
    # changes N_x by -Z_x*k_W^2/2. Moment rates instead include +M_Z*d.
    error=-q['Zx']*q['d'].dot(MZ*q['d'])/2
    assert s.expand(error+q['Zx']*q['e'].dot(q['d'])**2/2)==0
    assert error.subs({q['Zx']:1,q['d'][0]:1,q['d'][1]:0})==-q['Z']**2/2


def test_c1_quintic_elements_represent_both_local_initial_fields_exactly(data):
    x=s.symbols('x',real=True)
    coefficients=s.symbols('a0:6')
    g=sum(coefficients[i]*x**i for i in range(6))
    H=(1-3*x*x+2*x**3)*g.subs(x,0)+(x-2*x*x+x**3)*s.diff(g,x).subs(x,0)
    H+=(3*x*x-2*x**3)*g.subs(x,1)+(x**3-x*x)*s.diff(g,x).subs(x,1)
    bubble=x*x*(1-x)**2
    quotient,remainder=s.div(s.expand(g-H),s.expand(bubble),x)
    assert s.expand(remainder)==0
    assert s.degree(quotient,x)<=1
    assert s.expand(H+bubble*quotient-g)==0
    # A C2 quintic profile can give only C1 quartic momentum at a knot.
    h=s.symbols('h',positive=True)
    chi=(h*x)**3+(h*x)**5
    k=3*(h*x)**2+5*(h*x)**4
    assert s.degree(chi,x)==5 and s.degree(k,x)==4
    assert all(v.subs(x,0)==0 for v in (chi,k,s.diff(chi,x),s.diff(k,x)))
    assert s.diff(k,x,2).subs(x,0)==6*h*h
