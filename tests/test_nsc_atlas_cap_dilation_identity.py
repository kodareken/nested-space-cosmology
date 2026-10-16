"""Exact regular-polynomial moving-cap identities; no numerical fields."""
import time

import pytest
import sympy as s


CPU_LIMIT = 20.
started = None


@pytest.fixture(scope='module', autouse=True)
def receipt():
    global started
    started = time.process_time()
    yield
    elapsed = time.process_time()-started
    print('cap dilation identities: zero field samples;',elapsed,'active CPU seconds')
    assert elapsed < CPU_LIMIT, 'dilation algebra budget exhausted'


@pytest.fixture(autouse=True)
def algebra_budget():
    yield
    assert time.process_time()-started < CPU_LIMIT, 'dilation algebra budget exhausted'


@pytest.fixture(scope='module')
def data(receipt):
    xi,Z = s.symbols('xi Z', positive=True)
    powers = s.Matrix([1,3,5])
    basis = s.Matrix([xi**p for p in powers])
    M = Z*s.integrate(basis*basis.T,(xi,0,1))
    K = s.integrate(basis.diff(xi)*basis.diff(xi).T,(xi,0,1))/Z
    return xi,Z,powers,basis,M,K


def test_odd_fields_preserve_the_regular_centre(data):
    xi,Z,powers,basis,M,K = data
    z = s.symbols('z', real=True)
    c = s.Matrix(s.symbols('c1 c3 c5',real=True))
    d = s.Matrix(s.symbols('d1 d3 d5',real=True))
    chi,k = (v.dot(basis).subs(xi,z/Z) for v in (c,d))
    f_series = s.series(chi/s.tanh(z),z,0,4).removeO()
    pi_series = s.series(k*s.cosh(z)**2/s.tanh(z),z,0,4).removeO()
    assert s.simplify(f_series.subs(z,0)-c[0]/Z) == 0
    assert s.simplify(pi_series.subs(z,0)-d[0]/Z) == 0
    assert f_series.diff(z).subs(z,0) == 0
    assert pi_series.diff(z).subs(z,0) == 0
    E = s.diag(*powers)
    assert s.expand(z*s.diff(chi,z)-(E*c).dot(basis).subs(xi,z/Z)) == 0


def test_pure_dilation_leaves_physical_fields_and_centre_clock_fixed(data):
    xi,Z,powers,basis,M,K = data
    z,alpha,nu = s.symbols('z alpha nu',positive=True)
    c = s.Matrix(s.symbols('c1 c3 c5',real=True))
    d = s.Matrix(s.symbols('d1 d3 d5',real=True))
    E = s.diag(*powers)
    for coefficients in (c,d):
        field = coefficients.dot(basis).subs(xi,z/Z)
        transported = sum(s.diff(field,q)*rate for q,rate in zip(coefficients,nu*E*coefficients))
        assert s.simplify(transported+nu*Z*s.diff(field,Z)) == 0
    delta0 = -s.log(s.cosh(s.sqrt(alpha/6)*c[0]/Z))
    assert s.simplify(nu*c[0]*s.diff(delta0,c[0])+nu*Z*s.diff(delta0,Z)) == 0


def test_integrated_matrices_obey_dilation_identities(data):
    xi,Z,powers,basis,M,K = data
    e,E = s.ones(3,1),s.diag(*powers)
    # C/nu is integrated independently from the two mass pairings.
    transport = Z*s.integrate(basis*(xi*basis.diff(xi)).T,(xi,0,1))
    assert transport == M*E
    assert s.simplify(M*E+E*M-Z*e*e.T+M) == s.zeros(3)
    assert s.simplify(E*K+K*E-K-powers*powers.T/Z) == s.zeros(3)
    for n in range(1,4):
        assert s.simplify(M[:n,:n].det()/Z**n) > 0


def test_independent_quadratic_mass_shape_derivative(data):
    xi,Z,powers,basis,M,K = data
    c = s.Matrix(s.symbols('c1 c3 c5',real=True))
    d = s.Matrix(s.symbols('d1 d3 d5',real=True))
    b,A0 = s.symbols('b A0',positive=True)
    e,E = s.ones(3,1),s.diag(*powers)
    mass = (d.dot(M*d)+c.dot((K-b*e*e.T)*c))/4
    shape = s.diff(mass,Z)-(2*A0+b*b)*s.diff(mass,b)
    transport = sum(s.diff(mass,q)*rate for q,rate in zip(list(c)+list(d),list(E*c)+list(E*d)))/Z
    wall_chi,wall_k,slope = e.dot(c),e.dot(d),powers.dot(c)/Z
    density_times_A0 = A0*wall_chi**2/2+(wall_k**2+(slope-b*wall_chi)**2)/4
    assert s.expand(shape+transport-density_times_A0) == 0


def test_weak_pairing_plus_exact_shape_gives_full_moving_mass_balance(data):
    xi,Z,powers,basis,M,K = data
    c = s.Matrix(s.symbols('c1 c3 c5',real=True))
    d = s.Matrix(s.symbols('d1 d3 d5',real=True))
    gc = s.Matrix(s.symbols('gc1 gc3 gc5',real=True))
    gd = s.Matrix(s.symbols('gd1 gd3 gd5',real=True))
    clock,A0,rho,drive,zdot = s.symbols('clock A0 rho drive zdot',real=True,nonzero=True)
    e,E = s.ones(3,1),s.diag(*powers)
    V = 2*clock*M.inv()*gd
    rates_c = V+zdot*E*c/Z
    rates_d = M.inv()*(-2*clock*gc+e*(drive-e.dot(V)))+zdot*E*d/Z
    # The shape identity follows from fixed-physical-field transport, tested
    # above and derived for the constrained radial functional in the owner.
    shape = A0*rho-(gc.dot(E*c)+gd.dot(E*d))/Z
    rate = gc.dot(rates_c)+gd.dot(rates_d)+zdot*shape
    flux = e.dot(V)*(drive-e.dot(V))/(2*clock)
    assert s.expand(rate-flux-zdot*A0*rho) == 0
    assert s.expand(rate.subs(zdot,0)-flux) == 0


def test_weak_and_reconstructed_normal_traces_remain_distinct(data):
    xi,Z,powers,basis,M,K = data
    b,A0,drive = s.symbols('b A0 drive',real=True,nonzero=True)
    # alpha=0, chi=xi^3, k=0. The exact weak velocity is zero.
    chi = s.Matrix([0,1,0])
    slope = powers.dot(chi)/Z
    reconstructed_h = (slope-b)/A0
    weak_h = drive/A0
    assert s.simplify(reconstructed_h-weak_h-(3/Z-b-drive)/A0) == 0
    assert s.simplify((reconstructed_h-weak_h).subs(drive,3/Z-b)) == 0
