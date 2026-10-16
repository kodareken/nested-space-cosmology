"""Exact moving-domain controls; no field integration or production update."""
import time

import pytest
import sympy as s


CPU_LIMIT = 20.
active_cpu = 0.
started = None


@pytest.fixture(autouse=True)
def algebra_budget():
    yield
    assert time.process_time()-started < CPU_LIMIT, 'moving-projection algebra budget exhausted'


@pytest.fixture(scope='module', autouse=True)
def receipt():
    global started, active_cpu
    started = time.process_time()
    yield
    active_cpu = time.process_time()-started
    print('moving-projection identities: zero field samples;', active_cpu,
          'active CPU seconds')


def exact_elements(nodes):
    """Integrate basis products directly, then eliminate the centre node."""
    z = s.symbols('z', real=True)
    count = len(nodes)
    M,K,C = [s.zeros(count) for _ in range(3)]
    for i,(left,right) in enumerate(zip(nodes[:-1],nodes[1:])):
        phi = [(right-z)/(right-left),(z-left)/(right-left)]
        for a in range(2):
            for b in range(2):
                M[i+a,i+b] += s.integrate(phi[a]*phi[b],(z,left,right))
                K[i+a,i+b] += s.integrate(s.diff(phi[a],z)*s.diff(phi[b],z),(z,left,right))
                C[i+a,i+b] += s.integrate(phi[a]*z*s.diff(phi[b],z),(z,left,right))
    return M[1:,1:],K[1:,1:],C[1:,1:]


@pytest.fixture(scope='module')
def matrices(receipt):
    return exact_elements([s.S.Zero,s.Rational(1,2),s.S.One])


def test_general_element_transport_matrix_identity():
    z,l,h,nu = s.symbols('z l h nu', real=True, nonzero=True)
    phi = s.Matrix([(l+h-z)/h,(z-l)/h])
    M = s.integrate(phi*phi.T,(z,l,l+h))
    C = s.integrate(phi*(nu*z*phi.diff(z)).T,(z,l,l+h))
    endpoint = s.diag(-nu*l,nu*(l+h))
    assert s.simplify(C+C.T+nu*M-endpoint) == s.zeros(2)


def test_wall_potential_shape_derivative():
    Z = s.symbols('Z', positive=True)
    A0,b = 1/s.cosh(Z)**2,2/s.sinh(2*Z)
    assert s.simplify(s.diff(b,Z)+2*A0+b*b) == 0


def test_full_moving_mass_derivative_equals_flux_plus_explicit_remainder(matrices):
    M,K,C = matrices
    c1,c2,k1,k2,b,A0,d,v = s.symbols('c1 c2 k1 k2 b A0 d v', real=True)
    chi,k,e = s.Matrix([c1,c2]),s.Matrix([k1,k2]),s.Matrix([0,1])
    S = K-b*e*e.T
    advection = v*M.inv()*C*chi
    chi_rate = k+advection
    k_rate = M.inv()*(-S*chi+e*(d-k2)+v*C*k)
    shape = v*(k.dot(M*k)-chi.dot(K*chi)+(2*A0+b*b)*c2*c2)/4
    direct = (chi.dot(S*chi_rate)+k.dot(M*k_rate))/2+shape
    slope = 2*(c2-c1)
    flux = k2*(d-k2)/2
    moving_wall = v*(A0*c2*c2/2+(k2*k2+(slope-b*c2)**2)/4)
    remainder = (chi.dot(K*advection)/2-v*chi.dot(K*chi)/4
                 -v*slope*slope/4+b*c2*(v*slope-advection[-1])/2)
    assert s.expand(direct-flux-moving_wall-remainder) == 0
    assert s.expand(remainder).has(k1,k2,d) is False


def test_globally_affine_primitive_has_zero_moving_remainder(matrices):
    M,K,C = matrices
    slope,b,v = s.symbols('slope b v', real=True)
    chi = s.Matrix([slope/2,slope])
    advection = v*M.inv()*C*chi
    assert s.simplify(advection-v*chi) == s.zeros(2,1)
    remainder = (chi.dot(K*advection)/2-v*chi.dot(K*chi)/4
                 -v*slope*slope/4+b*chi[-1]*(v*slope-advection[-1])/2)
    assert s.expand(remainder) == 0


@pytest.mark.parametrize('speed',[s.Rational(1,2),-s.Rational(1,2),s.S.Zero])
def test_kinked_primitive_exposes_signed_motion_error(matrices,speed):
    M,K,C = matrices
    # chi_W=k_W=0; choose incoming d=-2 for the actual normal trace.
    # Thus the state can satisfy the incoming condition at this instant.
    # |Z_T|<1 makes the wall timelike in this zero-alpha reference.
    chi = s.Matrix([1,0])
    advection = speed*M.inv()*C*chi
    assert advection == speed*s.Matrix([s.Rational(1,7),-s.Rational(18,7)])
    shape_mass_rate = chi.dot(K*advection)/2-speed*chi.dot(K*chi)/4
    physical_moving_wall = speed # last slope=-2, so speed*slope^2/4
    assert shape_mass_rate == speed*s.Rational(13,7)
    assert shape_mass_rate-physical_moving_wall == speed*s.Rational(6,7)
