"""Independent exact mass integral for a background plus regular departure."""
import pytest
import sympy as s

from control_budget import ControlBudget


CONTROLS = ControlBudget(maximum_evaluations=10000,cpu_seconds=10.)


@pytest.fixture(scope='module')
def data():
    with CONTROLS.active() as budget:
        budget.check(6000)
        z,Z,T,c,d,beta,C,A0 = s.symbols('z Z T c d beta C A0',real=True,nonzero=True)
        phi = z/Z
        chi0 = (T+z)**4-(T-z)**4
        k0 = s.diff(chi0,T)
        chi,k = chi0+c*phi,k0+d*phi
        chiW,kW,slope = [value.subs(z,Z) for value in (chi,k,s.diff(chi,z))]
        mass = s.integrate(s.expand(k*k+s.diff(chi,z)**2),(z,0,Z))/4-beta*chiW**2/4
        gc,gd,NT = [s.diff(mass,value) for value in (c,d,T)]
        M = s.integrate(phi*phi,(z,0,Z))
        a,b = s.diff(chi0,T),s.diff(k0,T)
        ac,bc = [s.integrate(phi*value,(z,0,Z))/M for value in (a,b)]
        aW = a.subs(z,Z)
        v = s.expand(2*gd/M-ac)
        V = aW+v
        w = s.expand((-2*gc+C-V)/M-bc)
        direct = s.expand(NT+gc*v+gd*w)
        flux = V*(C-V)/2
        yield locals()
    print('affine background identities:',CONTROLS.budget.evaluations,
          'charged algebra values; zero field samples;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_direct_mass_derivative_has_the_background_trace_remainder(data):
    q = data
    remainder = (q['aW']-q['ac'])*(q['kW']+q['slope']-q['beta']*q['chiW']-q['C'])/2
    assert s.expand(q['direct']-q['flux']-remainder) == 0


def test_moving_shape_keeps_background_in_physical_coordinates(data):
    q = data
    shape = s.diff(q['mass'],q['Z'])-(2*q['A0']+q['beta']**2)*s.diff(q['mass'],q['beta'])
    transport = (q['gc']*q['c']+q['gd']*q['d'])/q['Z']
    density = q['A0']*q['chiW']**2/2+(q['kW']**2+(q['slope']-q['beta']*q['chiW'])**2)/4
    assert s.expand(shape+transport-density) == 0


def test_incompatible_incoming_exposes_error_and_compatible_limit_removes_it(data):
    q = data
    at = {q['T']:0,q['c']:0,q['d']:0}
    Z = q['Z']
    assert s.expand(q['direct'].subs(at).subs(q['C'],0)) == -s.Rational(96,5)*Z**6
    assert s.expand(q['flux'].subs(at).subs(q['C'],0)) == -32*Z**6
    defect = s.expand((q['direct']-q['flux']).subs(at).subs(q['C'],0))
    assert defect == s.Rational(64,5)*Z**6
    compatible = {**at,q['C']:8*Z**3}
    assert s.expand(q['v'].subs(compatible)) == 0
    assert s.expand(q['w'].subs(compatible)) == 0
    assert s.expand((q['direct']-q['flux']).subs(compatible)) == 0
