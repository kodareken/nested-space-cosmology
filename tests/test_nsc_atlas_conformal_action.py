"""Exact action/boundary/Legendre identities and physical branch controls."""
import math

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_conformal_action import (
    junction_matrix, reduced_bulk_density, tension_junction,
    time_kinetic_matrix, velocities_from_momenta,
)


def test_original_four_dimensional_action_and_GHY_reduction():
    r, chi, k = s.symbols('r chi k', positive=True)
    R4, box_r, dr2, dc2, mixed, V0, Kbar, nr = s.symbols('R4 box_r dr2 dc2 mixed V0 Kbar nr')
    X = r*r/(2*k)-chi*chi/12
    # Direct conformal substitutions into the ORIGINAL bulk integrand.
    Rg = (R4-6*box_r/r)/r**2
    dphi2 = dc2/r**2-2*chi*mixed/r**3+chi*chi*dr2/r**4
    original = r**4*((1/(2*k)-chi*chi/(12*r*r))*Rg-dphi2/(2*r*r)-V0)
    reduced = X*R4+3*dr2/k-dc2/2-V0*r**4
    divergence = (-3*r/k+chi*chi/(2*r))*box_r+(-3/k-chi*chi/(2*r*r))*dr2+chi*mixed/r
    assert s.expand(original-reduced-divergence) == 0
    original_GHY = (r*r/k-chi*chi/6)*(Kbar+3*nr/r)
    reduced_GHY = 2*X*Kbar
    assert s.expand(original_GHY-reduced_GHY+(-3*r/k+chi*chi/(2*r))*nr) == 0
    values = {r: .6, chi: .1, k: 2., R4: 5., dr2: -.3, dc2: .7, V0: .2}
    assert reduced_bulk_density(radius=.6, field=.1, kappa=2., radial_curvature=3.,
        grad_radius_squared=-.3, grad_field_squared=.7, vacuum_density=.2) == pytest.approx(float(reduced.subs(values)))


def test_direct_collar_variation_gives_the_three_junction_rows():
    r, chi, k, nr, nc, K = s.symbols('r chi k nr nc K')
    X = r*r/(2*k)-chi*chi/12
    Xn = s.diff(X, r)*nr+s.diff(X, chi)*nc
    # In d nbar^2-e^(2omega)d tau^2, Rbar2=-2(omega''+omega'^2).
    # After its GHY cancellation the first-normal-derivative density is:
    Lcollar = 2*Xn*K+3*nr*nr/k-nc*nc/2+2*X
    assert s.expand(s.diff(Lcollar, nr)-(6*nr+2*r*K)/k) == 0
    assert s.expand(s.diff(Lcollar, nc)+nc+chi*K/3) == 0
    assert s.expand(s.diff(Lcollar, K)-(2*r*nr/k-chi*nc/3)) == 0
    # Opposite outward normals subtract these momenta; the wall is -sigma*r^3.
    A, B, C, b = s.symbols('A B C b')
    matrix = s.Matrix([[r, -k*chi/6, 0], [3, 0, r], [0, 1, chi/3]])
    F = 1-k*chi*chi/(6*r*r)
    assert s.factor(matrix.det()+r*r*F) == 0
    assert matrix*s.Matrix([-b*r*r, 0, 0]) == s.Matrix([-b*r**3, -3*b*r*r, 0])


def test_time_legendre_matrix_is_the_hessian_of_the_same_action():
    a, r, chi, k, N = s.symbols('a r chi k N', positive=True)
    va, vr, vc = s.symbols('va vr vc')
    Lkin = (-2*(r*vr/k-chi*vc/6)*va-3*a*vr*vr/k+a*vc*vc/2)/N
    matrix = s.hessian(Lkin, (va, vr, vc))
    expected = s.Matrix([[0, -2*r/k, chi/3], [-2*r/k, -6*a/k, 0], [chi/3, 0, a]])/N
    assert matrix == expected
    F = 1-k*chi*chi/(6*r*r)
    assert s.factor(matrix.det()+4*a*r*r*F/(k*k*N**3)) == 0
    cfg = dict(radius=.6, field=.07, kappa=8*math.pi, radial_length=1.4, lapse=.9)
    numeric = time_kinetic_matrix(**cfg)
    evaluated = np.array(matrix.subs({r:.6, chi:.07, k:8*math.pi, a:1.4, N:.9}), dtype=float)
    assert numeric == pytest.approx(evaluated)
    velocity = np.array([.13, -.2, .35])
    recovered = velocities_from_momenta(numeric @ velocity, **cfg)
    assert recovered == pytest.approx(velocity)
    assert np.linalg.norm(numeric @ velocity-velocity) > .1


def test_ADM_kinetic_sign_from_direct_shifted_two_dimensional_curvature():
    t, y = s.symbols('t y', real=True)
    u, v = s.symbols('u v', real=True)
    a, N, b = s.exp(u*t), s.Integer(1), v*y
    metric = s.Matrix([[-N*N+a*a*b*b, a*a*b], [a*a*b, a*a]])
    inverse = metric.inv()
    coords = (t, y)
    Gamma = [[[s.simplify(sum(inverse[i, k]*(s.diff(metric[k, j], coords[m])
        +s.diff(metric[k, m], coords[j])-s.diff(metric[j, m], coords[k]))/2
        for k in range(2))) for m in range(2)] for j in range(2)] for i in range(2)]
    Ricci = s.zeros(2)
    for i in range(2):
        for j in range(2):
            Ricci[i, j] = sum(s.diff(Gamma[k][i][j], coords[k])-s.diff(Gamma[k][i][k], coords[j])
                +sum(Gamma[k][k][m]*Gamma[m][i][j]-Gamma[k][j][m]*Gamma[m][i][k]
                    for m in range(2)) for k in range(2))
    R2 = s.simplify(sum(inverse[i, j]*Ricci[i, j] for i in range(2) for j in range(2)))
    assert s.simplify(R2-2*(u-v)**2) == 0
    X = 1+t+y*y
    va, vX = s.diff(a, t)-s.diff(a*b, y), s.diff(X, t)-b*s.diff(X, y)
    boundary = 2*s.diff(X*va, t)-2*s.diff(X*b*va, y)
    assert s.simplify(a*X*R2+2*vX*va-boundary) == 0
    # The independent static-lapse limit retains the spatial derivative term.
    N = 1+y*y
    assert s.expand(N*X*(-2*s.diff(N, y, 2)/N)-2*s.diff(X, y)*s.diff(N, y)
                    +2*s.diff(X*s.diff(N, y), y)) == 0


@pytest.mark.parametrize('r,chi,k,b', [(.4, 0., 1., 2.), (.6, .07, 8*math.pi, 2.), (.8, .1, 2., 0.)])
def test_exact_junction_retains_radius_kink_and_homogeneous_invariant_field(r, chi, k, b):
    result = tension_junction(radius=r, field=chi, kappa=k, wall_k=b)
    assert result['normal_radius_jump'] == pytest.approx(-b*r*r)
    roundoff = 32*np.finfo(float).eps*max(1., abs(b*r*r))
    assert result['normal_field_jump'] == pytest.approx(0., abs=roundoff)
    assert result['extrinsic_curvature_jump'] == pytest.approx(0., abs=roundoff)
    # Transform back, including the derivative of phi_J=chi/r.
    normal_phi = result['normal_field_jump']/r**2-chi*result['normal_radius_jump']/r**3
    assert normal_phi == pytest.approx(b*chi/r)
    assert np.linalg.det(junction_matrix(radius=r, field=chi, kappa=k)) == pytest.approx(result['determinant'])


@pytest.mark.parametrize('r,chi,k', [(0., 0., 1.), (.5, 2., 1.), (.5, .1, -1.), (math.nan, .1, 1.)])
def test_reject_outside_physical_reduction_branch(r, chi, k):
    with pytest.raises(ValueError): junction_matrix(radius=r, field=chi, kappa=k)
