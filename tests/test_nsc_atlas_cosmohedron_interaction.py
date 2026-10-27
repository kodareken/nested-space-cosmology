"""Independent action, harmonic, actual-mode and wavefunction controls."""
import math

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_cosmohedron_interaction import (
    body_with_child, first_axial_mode, interaction, wall_vertex,
)
from recursive_horizons.nsc_atlas_gravitational_response import Budget, _advance


@pytest.fixture(scope="module")
def mode():
    return first_axial_mode(body_with_child(), rtol=2e-12)


def test_surface_metric_variation_and_common_clock_normalization():
    N, R, sine = sp.symbols("N R sine", positive=True)
    e, h, X, pt, pa = sp.symbols("e h X pt pa", real=True)
    # sine=sin(theta)>0 is the ordinary sphere chart 0<theta<pi.
    background = sp.diag(-N*N, R*R, R*R*sine**2)
    dg = sp.zeros(3)
    dg[0, 2] = dg[2, 0] = h*X*sine**2
    g = background+e*dg
    deriv = sp.Matrix([pt, 0, pa])
    L = -sp.sqrt(-g.det())*(deriv.T*g.inv()*deriv)[0]/2
    assert sp.simplify(sp.diff(L, e).subs(e, 0)+h*X*sine*pt*pa/N) == 0
    Qdot, Q, angular_eigenvalue = sp.symbols("Qdot Q angular_eigenvalue")
    canonical = R*R/(2*N)*(sp.sqrt(N)*Qdot/R)**2-N*angular_eigenvalue/2*(sp.sqrt(N)*Q/R)**2
    assert sp.simplify(canonical-(Qdot**2-N*N*angular_eigenvalue*Q*Q/R**2)/2) == 0


def test_real_harmonic_selection_antisymmetry_and_gravity_normalization():
    t, p = sp.symbols("t p", real=True)
    N20 = sp.sqrt(5/(16*sp.pi))
    N1, N2 = sp.sqrt(3/(4*sp.pi)), sp.sqrt(15/(4*sp.pi))
    mu = sp.symbols("mu", real=True)
    # Use mu=cos(theta) so the exact sphere integrals are polynomial.
    C12 = sp.simplify(6*N20*N1*N2*sp.integrate(mu**2*(1-mu**2), (mu, -1, 1))*sp.pi)
    C21 = -C12
    assert C12 == 3/(2*sp.sqrt(sp.pi)) and C21 == -C12
    assert N1**2*sp.pi*sp.integrate(1-mu**2, (mu, -1, 1)) == 1
    assert N2**2*sp.pi*sp.integrate(mu**2*(1-mu**2), (mu, -1, 1)) == 1
    # X_AB = D_(A X_B), with unit-sphere harmonic conventions.
    Xcov = sp.sin(t)**2*6*N20*sp.cos(t)
    Xtp = sp.diff(Xcov, t)/2-sp.cot(t)*Xcov
    assert sp.trigsimp(Xtp+3*N20*sp.sin(t)**3) == 0
    norm = 2*sp.pi*18*N20**2*sp.integrate((1-mu**2)**2, (mu, -1, 1))
    assert norm == 12  # (ell-1)ell(ell+1)(ell+2)/2
    G = sp.symbols("G", positive=True)
    assert norm/(32*sp.pi*G) == 24/(64*sp.pi*G)


def test_actual_mode_is_the_first_dirichlet_mode_with_independent_unit_norm(mode):
    assert 0 < mode.omega < mode.rayleigh_upper
    assert abs(mode.end[0]-math.pi) < 2e-11
    for _, lo, hi, sol in mode.parts:
        phases = sol.sol(np.linspace(lo, hi, 35))[0]
        assert np.all(phases > 0) and np.all(phases < math.pi+2e-11)
    assert mode.norm_quadrature(48) == pytest.approx(1., abs=3e-10)
    assert mode.norm_quadrature(96) == pytest.approx(1., abs=3e-10)


@pytest.mark.parametrize("index", [0, 1])
def test_canonical_mode_and_metric_vertex_match_both_sides_of_actual_wall(mode, index):
    shell = (mode.body.inner, mode.body.outer)[index]
    fm, pm = mode.wall_pair(index, "minus")
    fp, pp = mode.wall_pair(index, "plus")
    assert fp == pytest.approx(fm, abs=2e-14)
    assert pp-pm == pytest.approx(shell.jump_strength*fm, abs=2e-14)
    minus_h = shell.inside.speed(shell.radius)*fm+shell.radius*pm
    plus_h = shell.outside.speed(shell.radius)*fp+shell.radius*pp
    assert minus_h == pytest.approx(plus_h, abs=3e-14)
    assert wall_vertex(mode, index)["current_coupling"] != 0


def test_phase_mode_against_direct_regge_wheeler_propagation(mode):
    b = mode.body
    state = mode.wall_pair(0)
    budget = Budget()
    state = _advance(b.parent, b.inner.radius, b.outer.radius, state,
                     mode.omega, 2, budget, rtol=2e-12)
    state[1] += b.outer.jump_strength*state[0]
    state = _advance(b.exterior, b.outer.radius, b.observer, state,
                     mode.omega, 2, budget, rtol=2e-12)
    expected = mode.normalized_pair(mode.end)
    assert state == pytest.approx(expected, rel=2e-10, abs=2e-11)


def test_actual_momentum_vertices_have_the_cosmohedron_exchange_kernel():
    q1, q2, q3, q4, z = q = sp.symbols("q1 q2 q3 q4 z")
    a, b, c, d, y = w = sp.symbols("a b c d y", positive=True)
    L, R = sp.symbols("L R", real=True)
    EL, ER, Et = a+b+y, c+d+y, a+b+c+d
    S0 = sum(wi*qi*qi/2 for wi, qi in zip(w, q))
    SL = sp.I*L*(a-b)*q1*q2*z/EL
    SR = sp.I*R*(c-d)*q3*q4*z/ER
    SLR = L*R*(a-b)*(c-d)*q1*q2*q3*q4/(Et*EL*ER)
    S = S0+SL+SR+SLR
    residual = sum((wi**2*qi**2+sp.diff(S, qi, 2)-sp.diff(S, qi)**2)/2
                   for wi, qi in zip(w, q))-sum(w)/2
    residual += sp.I*z*(L*(q2*sp.diff(S, q1)-q1*sp.diff(S, q2))
                        +R*(q4*sp.diff(S, q3)-q3*sp.diff(S, q4)))
    mixed = sp.diff(residual, L, R).subs({L: 0, R: 0})
    assert sp.simplify(mixed) == 0
    assert SL.subs(b, a) == 0  # degenerate scalar pair cannot create from this vacuum


def test_fixed_G_scaling_and_inner_geometry_control(mode):
    base = interaction(mode)
    other_mode = first_axial_mode(body_with_child(.25, 2.), rtol=2e-12)
    other = interaction(other_mode)
    assert 2*other_mode.omega == pytest.approx(mode.omega, rel=2e-10)
    assert other["coefficient_shape"] == pytest.approx(base["coefficient_shape"], rel=3e-9)
    assert other["correlation_shape"] == pytest.approx(base["correlation_shape"], rel=3e-9)
    assert 16*other["log_wavefunction_exchange_coefficient"] == pytest.approx(
        base["log_wavefunction_exchange_coefficient"], rel=3e-9)
    control_mode = first_axial_mode(body_with_child(.30), rtol=2e-12)
    control = interaction(control_mode)
    assert control_mode.body.outer == mode.body.outer
    assert abs(control_mode.omega/mode.omega-1) < 2e-7
    assert control["coefficient_shape"]/base["coefficient_shape"] > 1.6
