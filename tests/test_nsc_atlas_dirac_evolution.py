"""Independent tetrad, density, Green identity and moving-wall controls."""
import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_dirac_evolution import (
    ALPHA, initial_wall_time_data, moving_wall_density_jump, radial_time_coefficients,
)
from recursive_horizons.nsc_atlas_dirac_slice import source_anchor_slice


def test_connection_contracted_from_full_diagonal_tetrad():
    N, a, r = s.symbols('N a r', positive=True)
    Nt, Nl, at, al, rt, rl = s.symbols('Nt Nl at al rt rl', real=True)
    g = s.diag(-N*N, a*a, r*r, r*r)  # At theta=pi/2; angular warp remains.
    derivatives = [s.diag(-2*N*Nt, 2*a*at, 2*r*rt, 2*r*rt),
                   s.diag(-2*N*Nl, 2*a*al, 2*r*rl, 2*r*rl), s.zeros(4), s.zeros(4)]
    inverse = g.inv()
    christoffel = [[[sum(inverse[i, k]*(derivatives[j][k, m]+derivatives[m][k, j]
        -derivatives[k][j, m])/2 for k in range(4)) for m in range(4)]
        for j in range(4)] for i in range(4)]
    tetrad = s.diag(1/N, 1/a, 1/r, 1/r)
    d_tetrad = [s.diag(-Nt/N**2, -at/a**2, -rt/r**2, -rt/r**2),
                s.diag(-Nl/N**2, -al/a**2, -rl/r**2, -rl/r**2), s.zeros(4), s.zeros(4)]
    pauli = [s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)]
    gamma = [s.I*s.diag(1, 1, -1, -1)]
    for sigma in pauli:
        gamma.append(s.I*s.BlockMatrix([[s.zeros(2), sigma], [-sigma, s.zeros(2)]]).as_explicit())
    contracted = s.zeros(4)
    for mu in range(4):
        connection = s.zeros(4)
        for A in range(4):
            for B in range(4):
                omega = sum(g[i, i]*tetrad[i, A]*(d_tetrad[mu][i, B]
                    +sum(christoffel[i][mu][j]*tetrad[j, B] for j in range(4))) for i in range(4))
                connection += omega*gamma[A]*gamma[B]/4
        contracted += sum((tetrad[mu, A]*gamma[A] for A in range(4)), s.zeros(4))*connection
    expected = gamma[0]/N*(at/(2*a)+rt/r)+gamma[1]/a*(Nl/(2*N)+rl/r)
    assert all(s.simplify(v)==0 for v in contracted-expected)


def test_half_density_removes_volume_time_term_and_retains_lapse_gradient():
    t, l = s.symbols('t l', real=True)
    N, a, r = [s.Function(v)(t, l) for v in ('N', 'a', 'r')]
    chi = 1+t*t+l**3
    weight = s.sqrt(a)*r
    psi = chi/weight
    temporal = weight*(s.diff(psi, t)+(s.diff(a, t)/(2*a)+s.diff(r, t)/r)*psi)
    radial = weight*N/a*(s.diff(psi, l)+(s.diff(N, l)/(2*N)+s.diff(r, l)/r)*psi)
    c = N/a
    assert s.simplify(temporal-s.diff(chi, t)) == 0
    assert s.simplify(radial-c*s.diff(chi, l)-s.diff(c, l)*chi/2) == 0
    row = radial_time_coefficients(lapse=2., lapse_l=.4, radial_length=1.3,
                                  radial_length_l=-.2, radius=.45)
    assert row['speed'] == pytest.approx(2/1.3)
    assert row['derivative_companion'] == pytest.approx((.4/1.3+2*.2/1.3**2)/2)


def test_initial_cap_time_generator_matches_optical_cap_after_unitary_map():
    r = s.symbols('r', positive=True)
    A = 1-r*r
    root_wall = s.sqrt(1-s.Rational(45, 100)**2)
    N = s.sqrt(A)/root_wall
    # d_l=sqrt(A)*d_r and chi=A^(-1/4)*u_optical.
    u = 1+r+r**3
    f = A**(-s.Rational(1, 4))
    derivative = N*s.sqrt(A)*s.diff(f*u, r)+s.sqrt(A)*s.diff(N, r)*f*u/2
    assert s.simplify(derivative/f-A/root_wall*s.diff(u, r)) == 0
    assert s.simplify((N/r)*f/f-(s.sqrt(A)/r)/root_wall) == 0


def test_variable_lapse_generator_obeys_green_identity():
    l = s.symbols('l', real=True)
    N = 1+l*l/3
    W = (1+l)/2
    alpha = s.Matrix([[0, -s.I], [s.I, 0]])
    sigma1 = s.Matrix([[0, 1], [1, 0]])
    f, g = s.Matrix([1+s.I*l, l*l]), s.Matrix([l**3, 2-s.I*l])
    def H(v): return -s.I*alpha*(N*v.diff(l)+s.diff(N, l)*v/2)+sigma1*W*v
    difference = (s.conjugate(f).T*H(g)-s.conjugate(H(f)).T*g)[0]
    boundary = -s.I*s.diff((s.conjugate(f).T*alpha*N*g)[0], l)
    assert s.simplify(difference-boundary) == 0


def test_moving_transparency_recovers_covariant_normal_jump():
    k, b = s.symbols('k b', positive=True)
    gamma = s.sqrt(1+b*b)
    speed = b/gamma
    alpha = s.Matrix([[0, -s.I], [s.I, 0]])
    psi = s.Matrix([1+s.I/3, s.Rational(-2, 7)+s.I/5])
    delta_a = -k/gamma
    delta_l = k/2*(gamma*s.eye(2)+b*alpha)*psi
    delta_t = -alpha*(delta_l+delta_a*psi/2)
    assert all(s.simplify(v)==0 for v in delta_t+speed*delta_l)
    psi_l = delta_l+k*gamma*psi
    psi_t = delta_t-k*b*psi
    assert all(s.simplify(v)==0 for v in gamma*psi_t+b*psi_l)
    covariant_normal = b*psi_t+gamma*psi_l+b*delta_a*alpha*psi/2
    assert all(s.simplify(v)==0 for v in covariant_normal-3*k*psi/2)
    # A constant lapse omits the temporal Israel jump and gives k instead
    # of 3*k/2, exposing the lost spacetime curvature term.
    wrong = b*(-k*b*psi)+gamma*(k*gamma*psi)
    assert all(s.simplify(v)==0 for v in wrong-k*psi)
    assert any(s.simplify(v)!=0 for v in wrong-3*k*psi/2)


def test_actual_temporal_israel_and_nonzero_moving_density_jets():
    wall = source_anchor_slice()
    row = initial_wall_time_data(wall)
    assert row['d_log_lapse_jump'] == pytest.approx(-2/wall['gamma'])
    assert row['reconstructed_Ktau_plus'] == pytest.approx(row['Ktau_plus'])
    assert row['Ktau_plus']-row['Ktau_minus'] == pytest.approx(-2)
    trace = np.array([1+.2j, -.3+.7j])
    dl = moving_wall_density_jump(trace, wall)
    dt = -ALPHA @ (dl+row['d_log_lapse_jump']*trace/2)
    assert dt+row['wall_coordinate_speed_over_N']*dl == pytest.approx([0., 0.], abs=3e-15)
    assert row['density_derivative_jump_alpha_eigenvalues'][0] > 0
    assert row['parallel_frame_normal_derivative_jump'] == 3.


@pytest.mark.parametrize('field,value', [('lapse', 0.), ('radius', 0.),
                                      ('radial_length', -1.), ('lapse_l', np.inf)])
def test_generator_rejects_invalid_geometry(field, value):
    args = dict(lapse=1., lapse_l=0., radius=.45, radial_length=1., radial_length_l=0.)
    args[field] = value
    with pytest.raises(ValueError): radial_time_coefficients(**args)
