"""Conformal-weight, transformed-junction and independent curvature controls."""
import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_dirac_conformal import (
    conformal_interface_data, radial_product_curvature,
)


def test_four_dimensional_conformal_spin_connection_and_normal_weight():
    n = s.symbols('n', real=True)
    sigma = n/7+n*n/13
    Omega = s.exp(sigma)
    u = 1+n+n**3
    # In a four-dimensional conformal coframe Omega*dx^A the contracted
    # connection is 3/2 gamma^A d_A(log Omega). A direct normal mode
    # tests D_(Omega^2 g)=Omega^-5/2 D_g Omega^3/2.
    psi = Omega**(-s.Rational(3, 2))*u
    transformed = (s.diff(psi, n)+s.Rational(3, 2)*s.diff(sigma, n)*psi)/Omega
    assert s.simplify(transformed-Omega**(-s.Rational(5, 2))*s.diff(u, n)) == 0
    # The Dirac norm on a spacelike three-slice is unchanged.
    assert s.simplify(Omega**3*(Omega**(-s.Rational(3, 2)))**2-1) == 0


def test_direct_normal_warp_tetrad_becomes_regular_at_tension_wall():
    n = s.symbols('n', real=True)
    # Two conformally warped half-metrics with the actual mixed K jump -k.
    # Choose r=rwall*warp, so Omega=rwall/r removes the entire first jet.
    rw = s.Rational(9, 20)
    k = s.Integer(2)
    slopes = (s.Rational(1, 7), s.Rational(1, 7)-k)
    spin_derivatives, transformed_derivatives = [], []
    for slope in slopes:
        warp = s.exp(slope*n)
        r = rw*warp
        psi = warp**(-s.Rational(3, 2))
        transformed = (r/rw)**s.Rational(3, 2)*psi
        spin_derivatives.append(s.diff(psi, n).subs(n, 0))
        transformed_derivatives.append(s.diff(transformed, n).subs(n, 0))
        assert s.simplify((rw/r)**2*warp**2-1) == 0
    assert spin_derivatives[1]-spin_derivatives[0] == 3*k/2
    assert transformed_derivatives == [0, 0]
    for radius in (.45, .5, .55):
        data = conformal_interface_data(radius)
        assert data['conformal_trace_K_jump'] == pytest.approx(0., abs=2e-15)
        assert data['conformal_normal_spinor_derivative_jump'] == pytest.approx(0., abs=2e-15)
        assert data['conformal_angular_K_jump'] == 0.
        assert data['conformal_temporal_K_jump'] == pytest.approx(0., abs=2e-15)


def test_curvature_from_direct_two_dimensional_conformal_metric():
    r = s.symbols('r', positive=True)
    H, M = s.Rational(3, 7), s.Rational(2, 5)
    for f, expected in ((1-H*H*r*r, -2), (1-2*M/r, -2+12*M/r)):
        # Direct two-dimensional Ricci contraction of diag(-f/r2,1/(f*r2)).
        g = s.diag(-f/r**2, 1/(f*r**2))
        inverse = g.inv()
        d = [s.zeros(2), g.diff(r)]
        Gamma = [[[s.simplify(sum(inverse[a, c]*(d[b][c, j]+d[j][c, b]-d[c][b, j])/2
                     for c in range(2))) for j in range(2)] for b in range(2)] for a in range(2)]
        Ricci = s.zeros(2)
        def deriv(expr, index): return s.diff(expr, r) if index == 1 else 0
        for i in range(2):
            for j in range(2):
                Ricci[i, j] = sum(deriv(Gamma[a][i][j], a)-deriv(Gamma[a][i][a], j)
                    +sum(Gamma[a][a][b]*Gamma[b][i][j]-Gamma[a][j][b]*Gamma[b][i][a]
                        for b in range(2)) for a in range(2))
        curvature = s.simplify(sum(inverse[i, j]*Ricci[i, j] for i in range(2) for j in range(2)))
        assert s.simplify(curvature-expected) == 0
        value = .5
        row = radial_product_curvature(radius=value, f=float(f.subs(r, value)),
            f_r=float(s.diff(f, r).subs(r, value)), f_rr=float(s.diff(f, r, 2).subs(r, value)))
        assert row['radial_scalar_curvature'] == pytest.approx(float(s.sympify(expected).subs(r, value)))


def test_curvature_is_regular_at_horizon_and_covariant_in_sphere_scale():
    for r in (.45, .8, 1.3):
        row = radial_product_curvature(radius=r, f=1-.8/r, f_r=.8/r**2, f_rr=-1.6/r**3)
        scaled = radial_product_curvature(radius=r, f=1-.8/r, f_r=.8/r**2, f_rr=-1.6/r**3, sphere_radius=7.)
        assert row['radial_scalar_curvature'] == pytest.approx(-2+4.8/r)
        assert scaled['radial_scalar_curvature'] == pytest.approx(row['radial_scalar_curvature']/49)
        assert scaled['four_dimensional_scalar_curvature'] == pytest.approx(row['four_dimensional_scalar_curvature']/49)
    # Minkowski and de Sitter share this radial product curvature; the physical
    # metric/clock and global endpoints are still different.
    assert radial_product_curvature(radius=.5, f=1., f_r=0., f_rr=0.)['radial_scalar_curvature'] == -2.


@pytest.mark.parametrize('r,L', [(0., 1.), (-1., 1.), (np.inf, 1.), (1., 0.)])
def test_reject_singular_conformal_map(r, L):
    with pytest.raises(ValueError):
        radial_product_curvature(radius=r, sphere_radius=L, f=1., f_r=0., f_rr=0.)
