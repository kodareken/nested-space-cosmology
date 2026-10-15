"""Independent exact checks for the unevolved next-order scalar problem.

These tests start from four-dimensional divergence operators, conformal
covariance, or manufactured smooth fields.  They do not evolve the physical
preparation and do not depend on a proposed scalar-response implementation.
The admission in docs/prior-art-and-open-claim.md permits 30 process-CPU
seconds per analytic/manufactured check.
"""
import time

import pytest
import sympy as sp


@pytest.fixture(autouse=True)
def analytic_cpu_budget():
    started = time.process_time()
    yield
    assert time.process_time() - started < 30, "analytic check exceeded its CPU admission"


def test_four_dimensional_conformal_covariance_and_reduced_field_invariance():
    t, r = sp.symbols("t r", positive=True)
    omega, phi = (sp.Function(name)(t, r) for name in ("omega", "phi"))

    def flat_box(f):
        return -sp.diff(f, t, 2) + sp.diff(f, r, 2) + 2 * sp.diff(f, r) / r

    def contraction(f, g):
        return -sp.diff(f, t) * sp.diff(g, t) + sp.diff(f, r) * sp.diff(g, r)

    # The exact metric omega^2 diag(-1,1,r^2,r^2 sin(theta)^2).
    transformed = phi / omega
    curved_box = (flat_box(transformed)
                  + 2 * contraction(omega, transformed) / omega) / omega**2
    curvature = -6 * flat_box(omega) / omega**3
    assert sp.simplify(curved_box - curvature * transformed / 6
                       - flat_box(phi) / omega**3) == 0
    assert sp.cancel((omega * r) * transformed - r * phi) == 0
    # Propagating the untransformed scalar is not the same equation.
    wrong_box = (flat_box(phi) + 2 * contraction(omega, phi) / omega) / omega**2
    assert sp.simplify(wrong_box - curvature * phi / 6
                       - flat_box(phi) / omega**3) != 0


def test_canonical_einstein_scalar_differs_from_conformal_scalar_at_cubic_order():
    phi, kappa = sp.symbols("phi kappa", real=True)
    conformal_factor = 1 - kappa * phi**2 / 6
    kinetic = (1 / conformal_factor
               + sp.Rational(3, 2) * sp.diff(conformal_factor, phi)**2
               / (kappa * conformal_factor**2))
    assert sp.factor(kinetic - 1 / conformal_factor**2) == 0
    canonical = sp.integrate(sp.series(1 / conformal_factor, phi, 0, 4).removeO(), phi)
    conformal = sp.series(phi / sp.sqrt(conformal_factor), phi, 0, 5).removeO()
    assert canonical == phi + kappa * phi**3 / 18
    assert conformal == phi + kappa * phi**3 / 12
    assert sp.expand(conformal - canonical) == kappa * phi**3 / 36


def test_einstein_curvature_perturbation_follows_from_the_fixed_jordan_trace():
    alpha, h, p, gradient_squared, box_p = sp.symbols("alpha h p gradient_squared box_p")
    omega = 1 - alpha * p**2 / 12
    # R_E=omega^-2[R_J-6 omega^-1 box_J(omega)], R_J=12h^2.
    # At first order box_J is the reference operator; the product rule supplies
    # box(p^2)=2p box(p)+2(gradient p)^2 without a metric field equation.
    box_omega = -alpha * (2 * p * box_p + 2 * gradient_squared) / 12
    curvature = (12 * h**2 - 6 * box_omega / omega) / omega**2
    coefficient = sp.diff(curvature, alpha).subs(alpha, 0)
    on_shell = coefficient.subs(box_p, 2 * h**2 * p)
    assert sp.expand(on_shell - gradient_squared - 4 * h**2 * p**2) == 0
    assert on_shell.subs(h, 0) == gradient_squared


def test_exact_exterior_reduction_from_four_dimensional_divergence():
    V, r = sp.symbols("V r", positive=True)
    A, Z, c, curvature = (sp.Function(name)(V, r) for name in ("A", "Z", "c", "R"))
    phi = c / r
    # Inverse of [[-Z^2 A,-Z],[-Z,0]] and sqrt(-g)=Z r^2 sin(theta).
    inverse = sp.Matrix([[0, -1 / Z], [-1 / Z, A]])
    metric = sp.Matrix([[-Z**2 * A, -Z], [-Z, 0]])
    assert sp.simplify(metric * inverse) == sp.eye(2)
    box4 = (sp.diff(-r**2 * sp.diff(phi, r), V)
            + sp.diff(-r**2 * sp.diff(phi, V) + Z * r**2 * A * sp.diff(phi, r), r)) / (Z * r**2)
    reduced = -2 * sp.diff(c, V, r) + sp.diff(Z * A * sp.diff(c, r), r)
    reduced -= (sp.diff(Z * A, r) / r + Z * curvature / 6) * c
    assert sp.simplify(Z * r * (box4 - curvature * phi / 6) - reduced) == 0


def test_exact_cap_reduction_from_four_dimensional_divergence():
    T, r = sp.symbols("T r", positive=True)
    A, E, c, curvature = (sp.Function(name)(T, r) for name in ("A", "E", "c", "R"))
    phi = c / r
    # E=exp(delta), inverse h^{TT}=-1/(E^2 A), h^{rr}=A.
    box4 = (sp.diff(-r**2 * sp.diff(phi, T) / (E * A), T)
            + sp.diff(E * r**2 * A * sp.diff(phi, r), r)) / (E * r**2)
    reduced = -sp.diff(sp.diff(c, T) / (E * A), T) + sp.diff(E * A * sp.diff(c, r), r)
    reduced -= (sp.diff(E * A, r) / r + E * curvature / 6) * c
    assert sp.simplify(E * r * (box4 - curvature * phi / 6) - reduced) == 0


def test_exterior_linearization_constraints_and_fixed_source_pullback():
    alpha, r, V, m = sp.symbols("alpha r V m", real=True)
    Z0 = sp.Function("Z0")(V)
    mu, zeta, c0, c1 = (sp.Function(name)(V, r) for name in ("mu", "zeta", "c0", "c1"))
    R1, pr = sp.symbols("R1 pr")
    A0 = 1 - 2 * m / r
    Z = Z0 * (1 + alpha * zeta)
    P = Z * (A0 - 2 * alpha * mu / r)
    c = c0 + alpha * c1
    exact = -2 * sp.diff(c, V, r) + sp.diff(P * sp.diff(c, r), r)
    exact -= (sp.diff(P, r) / r + Z * alpha * R1 / 6) * c
    coefficient = sp.diff(exact, alpha).subs(alpha, 0)
    P0, P1 = Z0 * A0, Z0 * (A0 * zeta - 2 * mu / r)
    L1 = -2 * sp.diff(c1, V, r) + sp.diff(P0 * sp.diff(c1, r), r) - sp.diff(P0, r) * c1 / r
    forcing = -sp.diff(P1 * sp.diff(c0, r), r) + (sp.diff(P1, r) / r + Z0 * R1 / 6) * c0
    assert sp.expand(coefficient - L1 + forcing) == 0
    constrained = sp.diff(P1, r).subs({sp.diff(mu, r): r**2 * A0 * pr**2 / 4,
                                      sp.diff(zeta, r): r * pr**2 / 2})
    assert sp.factor(constrained - Z0 * (sp.diff(A0, r) * zeta + 2 * mu / r**2)) == 0
    # Exact reduced data at the moved source fixes a nonzero background-curve trace.
    H = sp.symbols("H")
    pulled_back = sp.diff(c.subs(r, r + alpha * H), alpha).subs(alpha, 0)
    assert sp.simplify(pulled_back - c1 - H * sp.diff(c0, r)) == 0


def test_cap_linearization_and_radial_constraints_determine_the_potential():
    alpha, T, r, h = sp.symbols("alpha T r h", real=True)
    delta, nu, c0, c1 = (sp.Function(name)(T, r) for name in ("delta", "nu", "c0", "c1"))
    p, pT, pr = sp.symbols("p pT pr")
    A0 = 1 - h**2 * r**2
    A = A0 - 2 * alpha * nu / r
    E = 1 + alpha * delta
    B = E * A
    R1 = -pT**2 / A0 + A0 * pr**2 + 4 * h**2 * p**2
    c = c0 + alpha * c1
    exact = -sp.diff(sp.diff(c, T) / B, T) + sp.diff(B * sp.diff(c, r), r)
    exact -= (sp.diff(B, r) / r + E * (12 * h**2 + alpha * R1) / 6) * c
    coefficient = sp.diff(exact, alpha).subs(alpha, 0)
    B1 = A0 * delta - 2 * nu / r
    W = sp.diff(B1, r) / r + 2 * h**2 * delta + R1 / 6
    L1 = -sp.diff(c1, T, 2) / A0 + sp.diff(A0 * sp.diff(c1, r), r)
    forcing = -sp.diff(B1 * sp.diff(c0, T) / A0**2, T) - sp.diff(B1 * sp.diff(c0, r), r) + W * c0
    assert sp.cancel(coefficient - L1 + forcing) == 0
    mass_density = r**2 * (pT**2 / A0 + A0 * pr**2) / 4 + h**2 * r**2 * p**2 / 2
    lapse_rate = r * (pr**2 + pT**2 / A0**2) / 2
    constrained = W.subs({sp.diff(nu, r): mass_density, sp.diff(delta, r): lapse_rate})
    expected = 2 * nu / r**3 - h**2 * p**2 / 3 - pT**2 / (6 * A0) + A0 * pr**2 / 6
    assert sp.factor(constrained - expected) == 0
    # Removing the cap vacuum-induced quadratic potential changes the forcing.
    wrong = W.subs({sp.diff(nu, r): mass_density - h**2 * r**2 * p**2 / 2,
                    sp.diff(delta, r): lapse_rate})
    assert sp.factor(wrong - expected) == h**2 * p**2


def test_central_null_coordinate_pullback_cancels_both_metric_fluxes():
    u, v = sp.symbols("u v", real=True)
    f = sp.Function("f")
    c = f(v) - f(u)
    D, lm, lp = (sp.Function(name)(u, v) for name in ("D", "lambda_minus", "lambda_plus"))
    pullback = lm * sp.diff(c, u) + lp * sp.diff(c, v)
    mixed = sp.diff(pullback, u, v).subs({
        sp.diff(lm, u, v): sp.diff(D, u) / 2,
        sp.diff(lp, u, v): sp.diff(D, v) / 2,
        sp.diff(lm, v): D / 2,
        sp.diff(lp, u): D / 2,
    })
    flux = (sp.diff(D * sp.diff(c, u), u) + sp.diff(D * sp.diff(c, v), v)) / 2
    assert sp.expand(mixed - flux) == 0
    K = sp.symbols("K")
    # -d_T^2+d_xi^2=-4d_u d_v fixes the sign of the surviving potential.
    static_mixed = flux - K * c
    assert sp.expand(static_mixed - mixed + K * c) == 0
    assert sp.expand(static_mixed - mixed) != 0


def test_null_label_integrals_are_null_and_equal_the_clock_at_the_centre():
    T, xi, y = sp.symbols("T xi y", real=True)
    # Nonseparable manufactured D tests both moving paths and endpoint terms.
    D = 1 + T + 2 * xi + T * xi + xi**2
    lm = sp.integrate(D.subs({T: T - xi + y, xi: y}, simultaneous=True), (y, 0, xi))
    lp = -sp.integrate(D.subs({T: T + xi - y, xi: y}, simultaneous=True), (y, 0, xi))
    assert sp.expand((sp.diff(lm, T) + sp.diff(lm, xi)) / 2 - D / 2) == 0
    assert sp.expand((sp.diff(lp, T) - sp.diff(lp, xi)) / 2 - D / 2) == 0
    assert lm.subs(xi, 0) == lp.subs(xi, 0) == 0
    # Expanded eikonal equation for Q=T+sigma*xi+alpha*lambda.
    for sigma, shift in ((-1, lm), (1, lp)):
        eikonal_coefficient = 2 * (D - sp.diff(shift, T) + sigma * sp.diff(shift, xi))
        assert sp.expand(eikonal_coefficient) == 0


def test_regular_centre_potential_and_jordan_waveform_frame_cancellation():
    r, h, p0, pT0, p2 = sp.symbols("r h p0 pT0 p2", real=True)
    A0 = 1 - h**2 * r**2
    p, pr = p0 + p2 * r**2, 2 * p2 * r
    density = r**2 * (pT0**2 / A0 + A0 * pr**2) / 4 + h**2 * r**2 * p**2 / 2
    nu = sp.integrate(sp.series(density, r, 0, 6).removeO(), r)
    W = 2 * nu / r**3 - h**2 * p**2 / 3 - pT0**2 / (6 * A0) + A0 * pr**2 / 6
    assert sp.expand(nu).coeff(r, 3) == pT0**2 / 12 + h**2 * p0**2 / 6
    assert sp.limit(W, r, 0) == 0
    assert sp.limit(sp.diff(W, r), r, 0) == 0
    # At fixed T, lambda_minus,xi=D0 and lambda_plus,xi=-D0.
    # c0,u=-p0/2, c0,v=p0/2, while dr/dxi=1 at the centre.
    D0, central_null_derivative = -p0**2 / 12, sp.symbols("central_null_derivative")
    pullback_derivative = D0 * (-p0 / 2) - D0 * (p0 / 2)
    static_derivative = central_null_derivative + pullback_derivative
    jordan_coefficient = static_derivative - p0**3 / 12
    assert sp.expand(jordan_coefficient - central_null_derivative) == 0


def test_static_regular_cap_solution_has_a_nonzero_scattering_potential():
    r = sp.symbols("r", positive=True)
    A0 = 1 - r**2
    p = 2 * sp.atanh(r) / r
    # This is a smooth leading scalar solution, not an arbitrary metric source.
    wave = sp.diff(r**2 * A0 * sp.diff(p, r), r) / r**2 - 2 * p
    assert sp.simplify(wave) == 0
    assert sp.limit(p, r, 0) == 2
    density = r**2 * A0 * sp.diff(p, r)**2 / 4 + r**2 * p**2 / 2
    # Regular centre selects the zero integration constant.  Terms through r^7
    # determine K through r^4 without numerical quadrature or a PDE evolution.
    nu = sp.integrate(sp.series(density, r, 0, 8).removeO(), r)
    K = A0 * (2 * nu / r**3 - p**2 / 3 + A0 * sp.diff(p, r)**2 / 6) / 4
    assert sp.series(K, r, 0, 6).removeO() == 4 * r**2 / 135 + 4 * r**4 / 189
    # Consequently the mixed correction equation has a nonzero source even
    # for a static leading field; a homogeneous no-return correction is false.
    source = -K * r * p
    assert sp.series(source, r, 0, 4).removeO() == -8 * r**3 / 135


def test_no_return_outgoing_response_is_fixed_by_bulk_forcing_not_incoming_trace():
    u, v, z, w, u0, amplitude = sp.symbols("u v z w u0 amplitude", real=True)
    # A local no-return rectangle u0<=u<qa=0<=v.  The profile f(v)=v^4
    # joins a C3 quiet profile for v<=0; f(u)=0 throughout this rectangle.
    K = u**2 + u * v + 2 * v**2
    leading = v**4
    source = -K * leading
    retarded = sp.integrate(sp.integrate(source.subs({u: z, v: w}), (w, 0, v)), (z, u0, u))
    incoming = amplitude * v**4 + v**5
    correction = incoming + retarded
    assert sp.expand(sp.diff(correction, u, v) - source) == 0
    assert correction.subs(v, 0) == sp.diff(correction, u).subs(v, 0) == 0
    known_outgoing = sp.integrate(source.subs(v, w), (w, 0, v))
    assert sp.expand(sp.diff(correction, u) - known_outgoing) == 0
    assert sp.diff(correction, amplitude) != 0
    assert sp.diff(correction, u, amplitude) == 0
    assert known_outgoing != 0


def test_tangent_and_normal_transmission_select_equal_weighted_null_traces():
    Us, Vs, qs_minus, qs_plus, cu, cv, cm, cp = sp.symbols("Us Vs qs_minus qs_plus cu cv cm cp")
    tangent_jump = Us * cu + Vs * cv - qs_minus * cm - qs_plus * cp
    normal_jump = -Us * cu + Vs * cv + qs_minus * cm - qs_plus * cp
    assert sp.expand((tangent_jump - normal_jump) / 2 - Us * cu + qs_minus * cm) == 0
    assert sp.expand((tangent_jump + normal_jump) / 2 - Vs * cv + qs_plus * cp) == 0
    # In the exterior chart r_V=-ZA/2 and n^r=beta=-v-ZA V_s.
    # Consequently U_s r_U=(v-beta)/2, fixing the areal boundary factor.
    speed, beta, ZA = sp.symbols("speed beta ZA", nonzero=True)
    rV, wall_Vs = -ZA / 2, -(speed + beta) / ZA
    assert sp.expand(speed - wall_Vs * rV - (speed - beta) / 2) == 0


def test_fixed_clock_pullback_is_applied_once_when_starting_from_emission_label():
    label, alpha = sp.symbols("label alpha", real=True)
    f, shift, correction = (sp.Function(name) for name in ("f", "shift", "correction"))
    # If T=label+alpha*shift(label), invert at fixed T before reading the field.
    inverse_label = label - alpha * shift(label)
    at_fixed_clock = f(inverse_label) + alpha * correction(inverse_label)
    coefficient = sp.diff(at_fixed_clock, alpha).subs(alpha, 0)
    assert sp.simplify(coefficient - correction(label) + shift(label) * sp.diff(f(label), label)) == 0
    # A correction already expressed at fixed T is this coefficient itself.
    double_shifted = coefficient - shift(label) * sp.diff(f(label), label)
    assert sp.simplify(double_shifted - coefficient) == -shift(label) * sp.diff(f(label), label)
