"""Exact checks for the unevolved conservative scalar-response formulation.

These checks use coordinate changes, constraints and manufactured fields.
They do not solve the prepared physical scalar or introduce another metric
gauge. The analytic admission permits 30 process-CPU seconds per check.
"""
import time

import pytest
import sympy as sp


@pytest.fixture(autouse=True)
def analytic_cpu_budget():
    started = time.process_time()
    yield
    assert time.process_time() - started < 30, "analytic check exceeded its CPU admission"


def test_background_pullback_of_the_areal_operator_needs_no_metric_reconstruction():
    U, V, r = sp.symbols("U V r", real=True)
    ru = sp.symbols("r_U", nonzero=True)
    radius = sp.Function("radius")(U, V)
    c = sp.Function("c")
    P = sp.Function("P")(r, V)
    composed = c(radius, V)
    # The background ingoing characteristics obey r_V=-P/2. Their
    # U derivative fixes r_UV, rather than supplying an independent datum.
    mixed = sp.diff(composed, U, V).subs({
        sp.diff(radius, U, V): -ru * sp.diff(P, r) / 2,
        sp.diff(radius, U): ru,
        sp.diff(radius, V): -P / 2,
    }, simultaneous=True).subs(radius, r).doit()
    K = ru * sp.diff(P, r) / (2 * r)
    areal = (-2 * sp.diff(c(r, V), r, V)
             + sp.diff(P * sp.diff(c(r, V), r), r)
             - sp.diff(P, r) * c(r, V) / r)
    assert sp.simplify(areal + 2 * (mixed + K * c(r, V)) / ru) == 0
    # Reversing the mixed-coordinate sign fails even on a simple field.
    assert sp.simplify(areal - 2 * (mixed + K * c(r, V)) / ru) != 0


def test_mapped_metric_forcing_is_a_flux_plus_a_first_jet_density():
    r, V = sp.symbols("r V", real=True)
    ru, Z, curvature = sp.symbols("r_U Z R1", nonzero=True)
    c = sp.Function("c")(r, V)
    P1 = sp.Function("P1")(r, V)
    C = 2 * Z * ru
    areal_forcing = (-sp.diff(P1 * sp.diff(c, r), r)
                     + (sp.diff(P1, r) / r + Z * curvature / 6) * c)
    flux = P1 * sp.diff(c, r) / 2
    density = (ru * sp.diff(P1, r) / (2 * r) + C * curvature / 24) * c
    # d_U=r_U*d_r at fixed V; the flux can be integrated without c_rr.
    assert sp.simplify(-ru * areal_forcing / 2
                       - ru * sp.diff(flux, r) + density) == 0


def test_radial_constraints_and_null_jets_give_the_simple_density():
    r, m, Z, ru = sp.symbols("r m Z r_U", positive=True)
    pU, pV, c = sp.symbols("p_U p_V c", real=True)
    mu, zeta = (sp.Function(name)(r) for name in ("mu", "zeta"))
    A = 1 - 2 * m / r
    C = 2 * Z * ru
    rv = -Z * A / 2
    pr = pU / ru
    pV_areal = pV - rv * pU / ru
    curvature = sp.factor(A * pr**2 - 2 * pV_areal * pr / Z)
    assert sp.factor(curvature + 4 * pU * pV / C) == 0
    P1 = Z * (A * zeta - 2 * mu / r)
    P1r = sp.diff(P1, r).subs({
        sp.diff(mu, r): r**2 * A * pr**2 / 4,
        sp.diff(zeta, r): r * pr**2 / 2,
    })
    density = (ru * P1r / (2 * r) + C * curvature / 24) * c
    expected = (C * (m * zeta + mu) / (2 * r**3) - pU * pV / 6) * c
    assert sp.factor(density - expected) == 0
    without_curvature = ru * P1r * c / (2 * r)
    assert sp.factor(without_curvature - expected) == pU * pV * c / 6
    # c_U=0 on the wall does not imply p_U=0: p=c/r still varies radially.
    p = sp.symbols("p", nonzero=True)
    assert sp.factor((without_curvature - expected).subs(pU, -p * ru / r)) != 0


def test_source_pullback_derivative_uses_only_known_transverse_transport():
    V = sp.symbols("V", real=True)
    H, chiU, ru, chi, K, Pr = (sp.Function(name)(V) for name in
                               ("H", "chi_U", "r_U", "chi", "K", "P_r"))
    source = -H * chiU / ru
    derived = sp.diff(source, V).subs({
        sp.diff(chiU, V): -K * chi,
        sp.diff(ru, V): -Pr * ru / 2,
    })
    expected = (-(sp.diff(H, V) + H * Pr / 2) * chiU / ru
                + H * K * chi / ru)
    assert sp.factor(derived - expected) == 0
    wrong_fixed_areal_derivative = -sp.diff(H, V) * chiU / ru
    assert sp.factor(derived - wrong_fixed_areal_derivative) != 0


def test_time_independent_birth_lift_on_a_nonuniform_wall():
    U = sp.symbols("U", positive=True)
    V = sp.symbols("V", nonnegative=True)
    t = sp.symbols("t", nonnegative=True)
    wall = V + V**2 / 2
    lower = sp.sqrt(1 + 2 * U) - 1
    b = V**2
    # B(U_W(V))=integral U_W'(v)b(v)dv; B is a function of U only.
    primitive = sp.integrate((1 + t) * t**2, (t, 0, t))
    lift = primitive.subs(t, lower)
    assert sp.simplify(wall.subs(V, lower) - U) == 0
    assert sp.simplify(sp.diff(lift, U) - b.subs(V, lower)) == 0
    assert sp.diff(lift, V) == 0
    wall_lift = primitive.subs(t, V)
    assert sp.expand(sp.diff(wall_lift, V) - sp.diff(wall, V) * b) == 0
    assert lift.subs(U, 0) == 0
    assert sp.diff(lift, U).subs(U, 0) == 0


def test_lifted_volterra_rate_and_wall_extraction_on_a_manufactured_solution():
    U, V = sp.symbols("U V", positive=True)
    wall = 1 + V
    birth_lift = sp.Piecewise((0, U <= 1), ((U - 1)**3 / 3, True))
    source = V**2
    z = V * (U**2 / 2 - (1 + V) * U)
    y = z + source + birth_lift
    # Nonzero K, source flux and wall data, with a closed polynomial forcing.
    K = 1
    flux = V * (U + 1)**2
    density = (2 * V - 1) * U + 1 + 4 * V - y
    assert sp.simplify(sp.diff(y, U, V) + K * y
                       - sp.diff(flux, U) + density) == 0
    integrand = sp.simplify(K * y + density)
    rate = flux - flux.subs(U, 0) - sp.integrate(integrand, (U, 0, U))
    assert sp.expand(rate - sp.diff(z, V)) == 0
    assert z.subs(U, 0) == 0
    assert sp.diff(z, U).subs(U, wall) == 0
    assert sp.simplify(sp.diff(y, U).subs(U, wall)) == V**2
    assert y.subs({U: sp.Rational(1, 2), V: 0}) == 0
    wall_y = sp.simplify(y.subs(U, wall))
    reconstructed = rate.subs(U, wall) + sp.diff(source, V) + V**2
    assert sp.expand(sp.diff(wall_y, V) - reconstructed) == 0
    # The extended-grid node must carry the wall's total, not partial, rate.
    missing_wall_transport = rate.subs(U, wall) + sp.diff(source, V)
    assert sp.expand(sp.diff(wall_y, V) - missing_wall_transport) == V**2
    assert flux.subs(U, 0) == V  # Dropping the source flux also changes the solve.


def test_cap_optical_time_derivative_follows_from_radial_constraints_and_flux():
    T, r, h = sp.symbols("T r h", real=True)
    p, nu, delta = (sp.Function(name)(T, r) for name in ("p", "nu", "delta"))
    A = 1 - h**2 * r**2
    pT, pr = sp.diff(p, T), sp.diff(p, r)
    mass_rate = r**2 * (pT**2 / A + A * pr**2) / 4 + h**2 * r**2 * p**2 / 2
    lapse_rate = r * (pr**2 + pT**2 / A**2) / 2
    D = delta - 2 * nu / (r * A)
    Dr = sp.diff(D, r).subs({sp.diff(nu, r): mass_rate,
                             sp.diff(delta, r): lapse_rate})
    optical = 2 * nu * (1 - 3 * h**2 * r**2) / (r**2 * A**2) - h**2 * r * p**2 / A
    assert sp.factor(Dr - optical) == 0
    mixed = sp.diff(optical, T).subs(sp.diff(nu, T), r**2 * A * pT * pr / 2)
    expected = (1 - 3 * h**2 * r**2) * pT * pr / A - 2 * h**2 * r * p * pT / A
    assert sp.factor(mixed - expected) == 0


def test_independent_wave_equation_eliminates_second_jets_from_D_time_derivative():
    T, r, h = sp.symbols("T r h", real=True)
    p = sp.Function("p")(T, r)
    A = 1 - h**2 * r**2
    pT, pr = sp.diff(p, T), sp.diff(p, r)
    lapse_rate = r * (pr**2 + pT**2 / A**2) / 2
    mass_flux = r**2 * A * pT * pr / 2
    # Differentiate the lapse constraint and subtract d_r[2 nu_T/(r A)].
    mixed = sp.diff(lapse_rate, T) - sp.diff(2 * mass_flux / (r * A), r)
    pTT = A**2 * sp.diff(p, r, 2) + (2 * A**2 / r + A * sp.diff(A, r)) * pr - 2 * h**2 * A * p
    expected = (1 - 3 * h**2 * r**2) * pT * pr / A - 2 * h**2 * r * p * pT / A
    assert sp.factor(mixed.subs(sp.diff(p, T, 2), pTT) - expected) == 0
    wrong = mixed.subs(sp.diff(p, T, 2), pTT + 2 * h**2 * A * p)
    assert sp.factor(wrong - expected) == sp.factor(2 * h**2 * r * p * pT / A)
    assert sp.factor(wrong - expected) != 0


def test_regular_manufactured_field_fixes_the_optical_time_integration_constant():
    T, r, q = sp.symbols("T r q", real=True)
    # h=0 and f(q)=q^3 give the regular solution p=6T^2+2r^2.
    p = 6 * T**2 + 2 * r**2
    pT, pr = sp.diff(p, T), sp.diff(p, r)
    assert sp.diff(p, T, 2) - sp.diff(r**2 * pr, r) / r**2 == 0
    p0 = p.subs(r, 0)
    mass_rate = r**2 * (pT**2 + pr**2) / 4
    lapse_rate = r * (pr**2 + pT**2) / 2
    nu = sp.integrate(mass_rate.subs(r, q), (q, 0, r))
    delta = -p0**2 / 12 + sp.integrate(lapse_rate.subs(r, q), (q, 0, r))
    D = sp.expand(delta - 2 * nu / r)
    integral = sp.integrate((pT * pr).subs(r, q), (q, 0, r))
    proposed = -p0 * sp.diff(p0, T) / 6 + integral
    assert sp.expand(sp.diff(D, T) - proposed) == 0
    assert proposed == -12 * T**3 + 24 * T * r**2
    assert sp.limit(proposed, r, 0) == -p0 * sp.diff(p0, T) / 6
    assert sp.expand(sp.diff(D, T) - integral) == -12 * T**3
