"""First-order spherical feedback action in the owned (+---) chart.

The Weyl auxiliary is required: C_W = 0 is rejected and no pure-Einstein
chart is provided. The six canonical right-hand sides are functional
derivatives of the first-order density. They are not a Cauchy evolution:
initial data and discrete constraint preservation are not owned.

The 4D Riemann tensor is not recomputed here. Schwarzschild and de Sitter
use the already accepted 2D identity sqrt|h| R_h = -2 D_t + 2 (beta D + L_x/Q)_x
and C^2 = (R_h - 2)^2 / (3 r^4).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import isfinite, pi as math_pi

import sympy as sp

CAUCHY_EVOLUTION_OWNED = False
PURE_EINSTEIN_CHART = False
SOURCE_COUPLING_OWNED = False
# Naive Euler current without the mixed 16 v s_x / 16 v s_t terms is known wrong.
EULER_NAIVE_CURRENT_ALLOWED = False
# Periodic x and finite time caps. Dirichlet data on the caps are Q, r, chi.
BOUNDARY_ASSUMPTIONS = {
    "x_domain": "periodic",
    "time_domain": "finite_caps",
    "dirichlet_fields": ("Q", "r", "chi"),
    "spatial_corners": False,
}

_PI = sp.pi


def _uses_symbol(*values):
    return any(isinstance(value, sp.Expr) for value in values)


def _pi(*values):
    return _PI if _uses_symbol(*values) else math_pi


def require_weyl_chart(C_W, A=None):
    """Reject the singular auxiliary chart. Do not fall back to pure Einstein."""
    if C_W == 0:
        raise ValueError(
            "C_W=0 is rejected; this chart keeps the Weyl auxiliary and does not implement pure Einstein"
        )
    if A is not None and A == 0:
        raise ValueError("A=0 makes Z vanish; inverse velocities are outside this chart")


def alpha_of(C_W):
    require_weyl_chart(C_W)
    return -4 * _pi(C_W) * C_W / 3


def feedback_F(r, chi, A, C_W):
    """F = -4 pi A r^2 + 2 alpha chi, alpha = -4 pi C_W / 3."""
    alpha = alpha_of(C_W)
    return -4 * _pi(r, chi, A, C_W) * A * r**2 + 2 * alpha * chi


def feedback_Z(A):
    """Z = -24 pi A. Constant on this chart."""
    return -24 * _pi(A) * A


def feedback_V(r, chi, A, C_W, C_F, flux):
    """V = 8 pi A r^2 - alpha (4 chi + chi^2) - 2 pi C_F flux^2."""
    alpha = alpha_of(C_W)
    return (
        8 * _pi(r, chi, A, C_W, C_F, flux) * A * r**2
        - alpha * (4 * chi + chi**2)
        - 2 * _pi(r, chi, A, C_W, C_F, flux) * C_F * flux**2
    )


def partial_F(r, A, C_W):
    """(F_r, F_chi). F_chi = 2 alpha is constant and nonzero on this chart."""
    alpha = alpha_of(C_W)
    return -8 * _pi(r, A, C_W) * A * r, 2 * alpha


def shell_D(Q, Q_t, Q_x, beta, beta_x, L):
    """D = (Q_t - (beta Q)_x) / L, with (beta Q)_x = beta_x Q + beta Q_x."""
    if L == 0:
        raise ValueError("lapse factor L must be nonzero")
    return (Q_t - beta_x * Q - beta * Q_x) / L


def convective(time_derivative, space_derivative, beta):
    return time_derivative - beta * space_derivative


def first_order_density(
    L, L_x, Q, Q_t, Q_x, r, r_t, r_x, chi, chi_t, chi_x, beta, beta_x, A, C_W, C_F, flux,
):
    """Bulk density after removing div(-2 F D, 2 F (beta D + L_x/Q))."""
    F_r, F_chi = partial_F(r, A, C_W)
    Z = feedback_Z(A)
    V = feedback_V(r, chi, A, C_W, C_F, flux)
    D = shell_D(Q, Q_t, Q_x, beta, beta_x, L)
    u = convective(r_t, r_x, beta)
    v = convective(chi_t, chi_x, beta)
    F_x = F_r * r_x + F_chi * chi_x
    return (
        2 * D * (F_r * u + F_chi * v)
        - 2 * F_x * L_x / Q
        + Z * (Q / L * u**2 - L / Q * r_x**2)
        + L * Q * V
    )


def momenta(
    L, Q, Q_t, Q_x, r, r_t, r_x, chi, chi_t, chi_x, beta, beta_x, A, C_W,
):
    """p_Q, p_r, p_chi from the first-order density. C_F does not enter."""
    F_r, F_chi = partial_F(r, A, C_W)
    Z = feedback_Z(A)
    D = shell_D(Q, Q_t, Q_x, beta, beta_x, L)
    u = convective(r_t, r_x, beta)
    v = convective(chi_t, chi_x, beta)
    p_Q = 2 * (F_r * u + F_chi * v) / L
    p_r = 2 * F_r * D + 2 * Z * Q * u / L
    p_chi = 2 * F_chi * D
    return p_Q, p_r, p_chi


def momentum_pi(p_r, p_chi, r, A, C_W):
    """Pi = p_r - (F_r / F_chi) p_chi."""
    F_r, F_chi = partial_F(r, A, C_W)
    return p_r - (F_r / F_chi) * p_chi


def inverse_velocities(
    L, Q, Q_x, r, r_x, chi_x, beta, beta_x, p_Q, p_r, p_chi, A, C_W,
):
    """Solve D, u, v and the three time derivatives. Requires C_W != 0 and A != 0."""
    require_weyl_chart(C_W, A)
    if L == 0 or Q == 0:
        raise ValueError("L and Q must be nonzero")
    F_r, F_chi = partial_F(r, A, C_W)
    Z = feedback_Z(A)
    Pi = momentum_pi(p_r, p_chi, r, A, C_W)
    D = p_chi / (2 * F_chi)
    u = Pi * L / (2 * Z * Q)
    v = (L * p_Q / 2 - F_r * u) / F_chi
    Q_t = L * D + beta_x * Q + beta * Q_x
    r_t = u + beta * r_x
    chi_t = v + beta * chi_x
    return {"D": D, "u": u, "v": v, "Pi": Pi, "Q_t": Q_t, "r_t": r_t, "chi_t": chi_t}


def _dx_F_over_Q(r, r_x, r_xx, chi_x, chi_xx, Q, Q_x, A, C_W):
    F_r, F_chi = partial_F(r, A, C_W)
    # F_rr = -8 pi A, F_chi is constant.
    pi = _pi(r, A, C_W)
    F_rr = -8 * pi * A
    F_x = F_r * r_x + F_chi * chi_x
    dF_x = F_rr * r_x**2 + F_r * r_xx + F_chi * chi_xx
    return (dF_x * Q - F_x * Q_x) / Q**2, F_x


def lapse_constraint(
    p_Q, p_r, p_chi, Q, r, chi, r_x, r_xx, chi_x, chi_xx, Q_x, A, C_W, C_F, flux,
):
    """C = p_Q p_chi/(2 F_chi) + Pi^2/(4 Z Q) + Z r_x^2/Q - Q V - 2 dx(F_x/Q)."""
    require_weyl_chart(C_W, A)
    F_r, F_chi = partial_F(r, A, C_W)
    Z = feedback_Z(A)
    V = feedback_V(r, chi, A, C_W, C_F, flux)
    Pi = p_r - (F_r / F_chi) * p_chi
    dx_F_over_Q, _ = _dx_F_over_Q(r, r_x, r_xx, chi_x, chi_xx, Q, Q_x, A, C_W)
    return (
        p_Q * p_chi / (2 * F_chi)
        + Pi**2 / (4 * Z * Q)
        + Z * r_x**2 / Q
        - Q * V
        - 2 * dx_F_over_Q
    )


def shift_constraint(p_Q_x, p_r, p_chi, Q, r_x, chi_x):
    """Dcon = p_r r_x + p_chi chi_x - Q p_Q_x. Holds for nonconstant beta."""
    return p_r * r_x + p_chi * chi_x - Q * p_Q_x


def velocity_rhs(L, Q, Q_x, r, r_x, chi_x, beta, beta_x, p_Q, p_r, p_chi, A, C_W):
    """Three kinematic Hamilton equations, before any gauge fix."""
    solved = inverse_velocities(
        L, Q, Q_x, r, r_x, chi_x, beta, beta_x, p_Q, p_r, p_chi, A, C_W,
    )
    return solved["Q_t"], solved["r_t"], solved["chi_t"]


def ibp_current(F, D, beta, L_x, Q):
    """Removed current (-2 F D, 2 F (beta D + L_x/Q))."""
    return -2 * F * D, 2 * F * (beta * D + L_x / Q)


def einstein_conformal_current(A, N, q, r, beta, r_t, r_x):
    """Boundary between -4 pi A N q r^2 R and the conformal Einstein density.

    Locked from the prior Christoffel reduction. Not recomputed from a 4D tensor.
    """
    pi = _pi(A, N, q, r, beta, r_t, r_x)
    u = r_t - beta * r_x
    Jt = 24 * pi * A * q * r / N * u
    Jx = (
        -24 * pi * A * beta * q * r / N * r_t
        - 24 * pi * A * (N**2 - beta**2 * q**2) * r / (N * q) * r_x
    )
    return Jt, Jx


def box_r_flux(L, Q, r, beta, R_t, R_x):
    """Accepted Box-R flux. The counterterm cancels -C_box times this divergence.

    Jt = 4 pi Q r^2 / L * (R_t - beta R_x)
    Jx = 4 pi r^2 * (-beta Q/L * (R_t - beta R_x) - L/Q * R_x)
    """
    pi = _pi(L, Q, r, beta, R_t, R_x)
    convective_R = R_t - beta * R_x
    Jt = 4 * pi * Q * r**2 / L * convective_R
    Jx = 4 * pi * r**2 * (-beta * Q / L * convective_R - L / Q * R_x)
    return Jt, Jx


def euler_candidate_current(N, q, r, beta, *, derivatives):
    """Confirmed Euler current, including the mixed terms. Myers charge Q+ = -4 pi jt.

    k = (q_t - (beta q)_x)/N, v = (r_t - beta r_x)/N, s = r_x/q, P = N_x/q,
    A0 = 1 + v^2 - s^2,
    jt = 8 k A0 - 16 v s_x, jx = -8 (P + beta k) A0 + 16 v s_t.
    """
    k, v, s, P, A0, s_t, s_x = derivatives
    jt = 8 * k * A0 - 16 * v * s_x
    jx = -8 * (P + beta * k) * A0 + 16 * v * s_t
    return jt, jx


def euler_naive_current(k, A0, P, beta):
    """The current with the mixed terms dropped. Verified wrong; do not use."""
    return 8 * k * A0, -8 * (P + beta * k) * A0


def areal_curvature_scalar(N, q, r, beta, t, x):
    """R_h of h = L^2 dt^2 - Q^2 (dx+beta dt)^2, L=N/r, Q=q/r, from the 2D identity."""
    L = N / r
    Q = q / r
    D = (sp.diff(Q, t) - sp.diff(beta * Q, x)) / L
    divergence = -2 * sp.diff(D, t) + 2 * sp.diff(beta * D + sp.diff(L, x) / Q, x)
    return sp.simplify(divergence / (L * Q))


def weyl_square_from_rh(R_h, r):
    """Accepted product identity. Not a new 4D contraction."""
    return (R_h - 2) ** 2 / (3 * r**4)


def homogeneous_einstein_density(N, q, r, t, A=1):
    """Owned GHY-completed homogeneous Einstein density per unit axial coordinate."""
    q_t = sp.diff(q, t)
    r_t = sp.diff(r, t)
    return 8 * _PI * A * (N * q - 2 * r * r_t * q_t / N - q * r_t**2 / N)


class JetBundle:
    """Independent jet symbols and total derivatives up to a fixed order."""

    def __init__(self, fields, max_order=3):
        self.max_order = max_order
        self.table = {}
        for field in fields:
            jets = {}
            for nt in range(max_order + 1):
                for nx in range(max_order + 1 - nt):
                    suffix = "t" * nt + "x" * nx
                    jets[(nt, nx)] = sp.symbols(f"{field}{suffix}")
            self.table[field] = jets

    def __call__(self, field, nt=0, nx=0):
        try:
            return self.table[field][(nt, nx)]
        except KeyError as exc:
            raise RuntimeError(f"jet {field} t^{nt} x^{nx} exceeds order {self.max_order}") from exc

    def total(self, expr, which):
        if expr == 0:
            return sp.Integer(0)
        free = expr.free_symbols
        acc = sp.Integer(0)
        for jets in self.table.values():
            for (nt, nx), sym in jets.items():
                if sym not in free:
                    continue
                nxt = (nt + 1, nx) if which == "t" else (nt, nx + 1)
                higher = jets.get(nxt)
                if higher is None:
                    raise RuntimeError(f"total {which} derivative exceeds order at {sym}")
                acc += sp.diff(expr, sym) * higher
        return acc


def _total_over(expr, which, bundles):
    acc = sp.Integer(0)
    for bundle in bundles:
        acc += bundle.total(expr, which)
    return acc


def _numerator(expr):
    num, _den = sp.fraction(sp.together(expr))
    return sp.expand(num)


def _is_zero(expr):
    return _numerator(expr) == 0


def _geometry_jets(max_order=3):
    return JetBundle(("Q", "r", "chi", "L", "beta"), max_order=max_order)


def _lagrangian_on(jets, A, C_W, C_F, flux):
    return first_order_density(
        jets("L"), jets("L", 0, 1),
        jets("Q"), jets("Q", 1, 0), jets("Q", 0, 1),
        jets("r"), jets("r", 1, 0), jets("r", 0, 1),
        jets("chi"), jets("chi", 1, 0), jets("chi", 0, 1),
        jets("beta"), jets("beta", 0, 1),
        A, C_W, C_F, flux,
    )


def _forward_momenta_on(jets, A, C_W):
    return momenta(
        jets("L"), jets("Q"), jets("Q", 1, 0), jets("Q", 0, 1),
        jets("r"), jets("r", 1, 0), jets("r", 0, 1),
        jets("chi"), jets("chi", 1, 0), jets("chi", 0, 1),
        jets("beta"), jets("beta", 0, 1),
        A, C_W,
    )


@lru_cache(maxsize=1)
def legendre_identity():
    """p·velocity - density - (L C + beta Dcon) - dx(beta Q p_Q + 2 L F_x/Q) = 0."""
    jets = _geometry_jets(3)
    A, C_W, C_F, flux = sp.symbols("A C_W C_F Phi", nonzero=True)
    density = _lagrangian_on(jets, A, C_W, C_F, flux)
    p_Q, p_r, p_chi = _forward_momenta_on(jets, A, C_W)
    p_Q_x = jets.total(p_Q, "x")
    F_r, F_chi = partial_F(jets("r"), A, C_W)
    F_x = F_r * jets("r", 0, 1) + F_chi * jets("chi", 0, 1)
    dx_F_over_Q = jets.total(F_x / jets("Q"), "x")
    Pi = p_r - (F_r / F_chi) * p_chi
    Z = feedback_Z(A)
    V = feedback_V(jets("r"), jets("chi"), A, C_W, C_F, flux)
    C = (
        p_Q * p_chi / (2 * F_chi)
        + Pi**2 / (4 * Z * jets("Q"))
        + Z * jets("r", 0, 1) ** 2 / jets("Q")
        - jets("Q") * V
        - 2 * dx_F_over_Q
    )
    Dcon = shift_constraint(p_Q_x, p_r, p_chi, jets("Q"), jets("r", 0, 1), jets("chi", 0, 1))
    H = jets("L") * C + jets("beta") * Dcon
    phase = p_Q * jets("Q", 1, 0) + p_r * jets("r", 1, 0) + p_chi * jets("chi", 1, 0)
    boundary = jets.total(jets("beta") * jets("Q") * p_Q + 2 * jets("L") * F_x / jets("Q"), "x")
    residual = sp.together(phase - density - H - boundary)
    return {"residual_is_zero": _is_zero(residual), "residual": residual}


@lru_cache(maxsize=1)
def constraint_sign_identities():
    """Lapse EL equals -C and shift EL equals -Dcon, including nonconstant beta."""
    jets = _geometry_jets(3)
    A, C_W, C_F, flux = sp.symbols("A C_W C_F Phi", nonzero=True)
    density = _lagrangian_on(jets, A, C_W, C_F, flux)
    # First-order density: derivatives of L and beta stop at first order.
    el_L = sp.diff(density, jets("L")) - jets.total(sp.diff(density, jets("L", 0, 1)), "x")
    el_beta = sp.diff(density, jets("beta")) - jets.total(sp.diff(density, jets("beta", 0, 1)), "x")
    p_Q, p_r, p_chi = _forward_momenta_on(jets, A, C_W)
    F_r, F_chi = partial_F(jets("r"), A, C_W)
    F_x = F_r * jets("r", 0, 1) + F_chi * jets("chi", 0, 1)
    dx_F_over_Q = jets.total(F_x / jets("Q"), "x")
    Pi = p_r - (F_r / F_chi) * p_chi
    Z = feedback_Z(A)
    V = feedback_V(jets("r"), jets("chi"), A, C_W, C_F, flux)
    C = (
        p_Q * p_chi / (2 * F_chi)
        + Pi**2 / (4 * Z * jets("Q"))
        + Z * jets("r", 0, 1) ** 2 / jets("Q")
        - jets("Q") * V
        - 2 * dx_F_over_Q
    )
    Dcon = shift_constraint(
        jets.total(p_Q, "x"), p_r, p_chi, jets("Q"), jets("r", 0, 1), jets("chi", 0, 1),
    )
    return {
        "lapse_el_plus_C_is_zero": _is_zero(el_L + C),
        "shift_el_plus_Dcon_is_zero": _is_zero(el_beta + Dcon),
        "lapse_sign": "EL_L = -C",
        "shift_sign": "EL_beta = -Dcon",
    }


@lru_cache(maxsize=1)
def hamiltonian_right_hand_sides():
    """Six bulk Hamilton equations derived from H = L C + beta Dcon.

    Velocity equations are the functional derivatives with respect to the momenta.
    Momentum equations are checked against the Lagrangian Euler operators.
    """
    metric = _geometry_jets(3)
    phase = JetBundle(("pQ", "pr", "pchi"), max_order=2)
    A, C_W, C_F, flux = sp.symbols("A C_W C_F Phi", nonzero=True)
    F_r, F_chi = partial_F(metric("r"), A, C_W)
    Z = feedback_Z(A)
    V = feedback_V(metric("r"), metric("chi"), A, C_W, C_F, flux)
    F_x = F_r * metric("r", 0, 1) + F_chi * metric("chi", 0, 1)
    dx_F_over_Q = metric.total(F_x / metric("Q"), "x")
    Pi = phase("pr") - (F_r / F_chi) * phase("pchi")
    C = (
        phase("pQ") * phase("pchi") / (2 * F_chi)
        + Pi**2 / (4 * Z * metric("Q"))
        + Z * metric("r", 0, 1) ** 2 / metric("Q")
        - metric("Q") * V
        - 2 * dx_F_over_Q
    )
    Dcon = (
        phase("pr") * metric("r", 0, 1)
        + phase("pchi") * metric("chi", 0, 1)
        - metric("Q") * phase("pQ", 0, 1)
    )
    H = metric("L") * C + metric("beta") * Dcon

    bundles = (metric, phase)

    def delta(bundle, density, field, spatial_order):
        # Total derivatives hit both geometry and momentum jets.
        acc = sp.diff(density, bundle(field))
        if spatial_order >= 1 and (0, 1) in bundle.table[field]:
            acc -= _total_over(sp.diff(density, bundle(field, 0, 1)), "x", bundles)
        if spatial_order >= 2 and (0, 2) in bundle.table[field]:
            d2 = sp.diff(density, bundle(field, 0, 2))
            acc += _total_over(_total_over(d2, "x", bundles), "x", bundles)
        return acc

    # Momenta have no second x-derivative in H. Q, r, chi: C contains dx(F_x/Q).
    Q_t = delta(phase, H, "pQ", 1)
    r_t = delta(phase, H, "pr", 0)
    chi_t = delta(phase, H, "pchi", 0)
    pQ_t = -delta(metric, H, "Q", 1)
    pr_t = -delta(metric, H, "r", 2)
    pchi_t = -delta(metric, H, "chi", 2)

    solved = inverse_velocities(
        metric("L"), metric("Q"), metric("Q", 0, 1), metric("r"), metric("r", 0, 1),
        metric("chi", 0, 1), metric("beta"), metric("beta", 0, 1),
        phase("pQ"), phase("pr"), phase("pchi"), A, C_W,
    )
    velocity_match = (
        _is_zero(Q_t - solved["Q_t"])
        and _is_zero(r_t - solved["r_t"])
        and _is_zero(chi_t - solved["chi_t"])
    )

    density = _lagrangian_on(metric, A, C_W, C_F, flux)
    p_Q, p_r, p_chi = _forward_momenta_on(metric, A, C_W)
    el = {}
    for field in ("Q", "r", "chi"):
        el[field] = (
            sp.diff(density, metric(field))
            - metric.total(sp.diff(density, metric(field, 1, 0)), "t")
            - metric.total(sp.diff(density, metric(field, 0, 1)), "x")
        )
    replacement = {}
    for name, expr in (("pQ", p_Q), ("pr", p_r), ("pchi", p_chi)):
        replacement[phase(name)] = expr
        replacement[phase(name, 0, 1)] = metric.total(expr, "x")
        replacement[phase(name, 0, 2)] = metric.total(metric.total(expr, "x"), "x")
    # p_t + δH/δfield should cancel the Lagrangian EL on the Legendre image.
    momentum_residuals = {
        "Q": metric.total(p_Q, "t") - pQ_t.xreplace(replacement) + el["Q"],
        "r": metric.total(p_r, "t") - pr_t.xreplace(replacement) + el["r"],
        "chi": metric.total(p_chi, "t") - pchi_t.xreplace(replacement) + el["chi"],
    }
    momentum_match = {name: _is_zero(expr) for name, expr in momentum_residuals.items()}
    return {
        "Q_t": Q_t,
        "r_t": r_t,
        "chi_t": chi_t,
        "pQ_t": pQ_t,
        "pr_t": pr_t,
        "pchi_t": pchi_t,
        "velocity_rhs_matches_inverse": velocity_match,
        "momentum_rhs_matches_euler": momentum_match,
        "momentum_residuals_are_zero": all(momentum_match.values()),
    }


def auxiliary_square_identity():
    """chi = R_h - 2 recovers alpha (R_h - 2)^2 and is the unique critical point."""
    alpha, chi, R_h = sp.symbols("alpha chi R_h")
    packed = alpha * (2 * chi * (R_h - 2) - chi**2)
    restored = sp.simplify(packed.subs(chi, R_h - 2) - alpha * (R_h - 2) ** 2)
    stationarity = sp.diff(packed, chi)
    solved = sp.solve(stationarity, chi)
    return {
        "restored_is_zero": restored == 0,
        "critical_point": solved,
        "critical_point_is_curvature": solved == [R_h - 2],
    }


def homogeneous_ghy_match():
    """Einstein sector of the first-order density against the owned homogeneous GHY density."""
    t = sp.symbols("t")
    N, q, r = (sp.Function(name, positive=True)(t) for name in ("N", "q", "r"))
    A, C_W = sp.symbols("A C_W", nonzero=True)
    L = N / r
    Q = q / r
    density = first_order_density(
        L, 0,
        Q, sp.diff(Q, t), 0,
        r, sp.diff(r, t), 0,
        0, 0, 0,
        0, 0,
        A, C_W, 0, 0,
    )
    residual = sp.simplify(density - homogeneous_einstein_density(N, q, r, t, A))
    return {"residual_is_zero": residual == 0, "residual": residual}


def box_flux_identity():
    """div J equals the angular integral of d(sqrt-g grad R) built from the 2-metric inverse."""
    t, x = sp.symbols("t x")
    N, q, r, beta, R = (sp.Function(name)(t, x) for name in ("N", "q", "r", "beta", "R"))
    L = N / r
    Q = q / r
    Jt, Jx = box_r_flux(L, Q, r, beta, sp.diff(R, t), sp.diff(R, x))
    divergence = sp.diff(Jt, t) + sp.diff(Jx, x)
    # g^{ab} = r^{-2} h^{ab}, angular integral of sqrt-g g^{a b} R_b is J.
    h_inv_tt = 1 / L**2
    h_inv_tx = -beta / L**2
    h_inv_xx = -1 / Q**2 + beta**2 / L**2
    measure = 4 * _PI * N * q  # equals 4 pi r^2 * L * Q * r^0? N q = L Q r^2; times r^0 from r^2 * r^{-2}
    # sqrt-g angular = 4 pi N q r^2, times g^{ab}=r^{-2} h^{ab}, so factor 4 pi N q.
    factor = 4 * _PI * N * q
    expected_t = factor * (h_inv_tt * sp.diff(R, t) + h_inv_tx * sp.diff(R, x))
    expected_x = factor * (h_inv_tx * sp.diff(R, t) + h_inv_xx * sp.diff(R, x))
    expected = sp.diff(expected_t, t) + sp.diff(expected_x, x)
    residual = sp.simplify(divergence - expected)
    # Homogeneous owned endpoint: 4 pi q r^2 R_t / N.
    th = sp.symbols("tau")
    Nh, qh, rh, Rh = (sp.Function(name, positive=True)(th) for name in ("N", "q", "r", "R"))
    Jt_h, Jx_h = box_r_flux(Nh / rh, qh / rh, rh, 0, sp.diff(Rh, th), 0)
    owned = 4 * _PI * qh * rh**2 * sp.diff(Rh, th) / Nh
    C_box = sp.symbols("C_box")
    cap_density = -C_box * Jt + C_box * Jt
    return {
        "divergence_matches_reduced_box": residual == 0,
        "homogeneous_endpoint_matches_owner": sp.simplify(Jt_h - owned) == 0 and sp.simplify(Jx_h) == 0,
        "cap_cancellation_is_zero": sp.simplify(cap_density) == 0,
    }


def _euler_sections(jets):
    """Kinematic sections used by the candidate Euler density. n(v) = (v_t - beta v_x)/N."""
    N, q, beta = jets("N"), jets("q"), jets("beta")
    k = (jets("q", 1, 0) - jets.total(beta * q, "x")) / N
    v = (jets("r", 1, 0) - beta * jets("r", 0, 1)) / N
    s = jets("r", 0, 1) / q
    P = jets("N", 0, 1) / q
    A0 = 1 + v**2 - s**2
    normal_v = (jets.total(v, "t") - beta * jets.total(v, "x")) / N
    A_sec = normal_v - (P / N) * s
    B_sec = jets.total(s, "x") / q - (k / q) * v
    C_sec = jets.total(v, "x") / q - (k / q) * s
    R2 = (-2 * jets.total(k, "t") + 2 * jets.total(beta * k + P, "x")) / (N * q)
    density = N * q * (-4 * R2 * A0 + 16 * (C_sec**2 - A_sec * B_sec))
    s_t, s_x = jets.total(s, "t"), jets.total(s, "x")
    jt, jx = euler_candidate_current(N, q, jets("r"), beta, derivatives=(k, v, s, P, A0, s_t, s_x))
    naive_t, naive_x = euler_naive_current(k, A0, P, beta)
    return {
        "density": density,
        "jt": jt,
        "jx": jx,
        "naive_jt": naive_t,
        "naive_jx": naive_x,
        "divergence": jets.total(jt, "t") + jets.total(jx, "x"),
        "naive_divergence": jets.total(naive_t, "t") + jets.total(naive_x, "x"),
    }


@lru_cache(maxsize=1)
def euler_candidate_status():
    """Symbolic residual of the candidate current. Pending unless the numerator vanishes.

    The section formula is E_4 r^2 sqrt|gamma| = sqrt|gamma| [-4 R_2 A0 + 16 (C^2 - A B)]
    with sqrt|gamma| = N q. A = n(v) - (P/N) s, B = s_x/q - (k/q) v, C = v_x/q - (k/q) s.
    """
    jets = JetBundle(("N", "q", "r", "beta"), max_order=3)
    parts = _euler_sections(jets)
    residual = sp.together(parts["divergence"] - parts["density"])
    closed = _is_zero(residual)
    naive_same = _is_zero(parts["divergence"] - parts["naive_divergence"])
    return {
        "status": "checked" if closed else "pending",
        "residual_is_zero": closed,
        "naive_current_is_distinct": not naive_same,
        "formula": "jt=8 k A0-16 v s_x; jx=-8 (P+beta k) A0+16 v s_t",
        "rejected_naive_formula": "jt=8 k A0; jx=-8 (P+beta k) A0",
    }


def _symbol_named(expr, name):
    for symbol in expr.free_symbols:
        if symbol.name == name:
            return symbol
    raise KeyError(name)


@lru_cache(maxsize=1)
def explicit_hamilton_rhs():
    """Factored six Hamilton equations, derived from H = L C + beta Dcon."""
    raw = hamiltonian_right_hand_sides()
    keys = ("Q_t", "r_t", "chi_t", "pQ_t", "pr_t", "pchi_t")
    displayed = {key: sp.factor(sp.simplify(raw[key])) for key in keys}
    names = {}
    for expr in displayed.values():
        for symbol in expr.free_symbols:
            names[symbol.name] = symbol
    A, C_W, L, Q = names["A"], names["C_W"], names["L"], names["Q"]
    beta, pchi, pr, pQ = names["beta"], names["pchi"], names["pr"], names["pQ"]
    r, chi = names["r"], names["chi"]
    F_r, F_chi = partial_F(r, A, C_W)
    Z = feedback_Z(A)
    Pi = pr - (F_r / F_chi) * pchi
    Q_x, beta_x = names["Qx"], names["betax"]
    r_x, chi_x = names["rx"], names["chix"]
    expected_Q = L * pchi / (2 * F_chi) + beta_x * Q + beta * Q_x
    expected_r = L * Pi / (2 * Z * Q) + beta * r_x
    expected_chi = L * pQ / (2 * F_chi) - (F_r / F_chi) * (L * Pi / (2 * Z * Q)) + beta * chi_x
    Lx, Lxx, Qx = names["Lx"], names["Lxx"], names["Qx"]
    V_chi = sp.diff(feedback_V(r, chi, A, C_W, 0, 0), chi)
    dx_Lx_over_Q = (Lxx * Q - Lx * Qx) / Q**2
    expected_pchi = (
        L * Q * V_chi
        + 2 * F_chi * dx_Lx_over_Q
        + names["betax"] * pchi
        + beta * names["pchix"]
    )
    matches = {
        "Q_t": _is_zero(displayed["Q_t"] - expected_Q),
        "r_t": _is_zero(displayed["r_t"] - expected_r),
        "chi_t": _is_zero(displayed["chi_t"] - expected_chi),
        "pchi_t": _is_zero(displayed["pchi_t"] - expected_pchi),
        "pQ_t": _is_zero(displayed["pQ_t"] - raw["pQ_t"]),
        "pr_t": _is_zero(displayed["pr_t"] - raw["pr_t"]),
    }
    return {
        "displayed": {key: displayed[key] for key in keys},
        "matches_derived_equations": matches,
        "all_match": all(matches.values()),
    }


def boundary_completion_ledger():
    """Cap counterterms on periodic x. No spatial corners and no C_W=0 chart.

    Dirichlet data on the finite time caps are Q, r and chi. The geometric
    counterterms below cancel the corresponding bulk divergences. Einstein is
    already the GHY completion, so the only extra time flux is the Weyl piece
    4 alpha chi D.
    """
    C_E, C_box, jt, Jt, F, D = sp.symbols("C_E C_box jt Jt_R F D")
    euler_cap = -4 * _PI * C_E * jt + 4 * _PI * C_E * jt
    box_cap = -C_box * Jt + C_box * Jt
    removed_time = -2 * F * D
    restored_time = 2 * F * D
    Q_plus = -4 * _PI * jt
    p_Q, p_r, p_chi = sp.symbols("p_Q p_r p_chi")
    d_Q, d_r, d_chi = sp.symbols("delta_Q delta_r delta_chi")
    cap_form = p_Q * d_Q + p_r * d_r + p_chi * d_chi
    dirichlet_cap = cap_form.subs({d_Q: 0, d_r: 0, d_chi: 0})
    angle = sp.symbols("x")
    periodic_density = sp.sin(angle) + sp.cos(2 * angle)
    periodic_flux = sp.integrate(sp.diff(periodic_density, angle), (angle, 0, 2 * _PI))

    A, radius, C_W, chi, Dsym = sp.symbols("A r C_W chi D", nonzero=True)
    alpha = alpha_of(C_W)
    F_eh = -4 * _PI * A * radius**2
    F_full = F_eh + 2 * alpha * chi
    weyl_extra = sp.simplify(2 * F_full * Dsym - 2 * F_eh * Dsym - 4 * alpha * chi * Dsym)
    rebuild = sp.simplify(2 * F_eh * Dsym + 4 * alpha * chi * Dsym - 2 * F_full * Dsym)

    time = sp.symbols("t")
    lapse, radial, areal = (sp.Function(name, positive=True)(time) for name in ("N", "q", "r"))
    mass = sp.symbols("A", nonzero=True)
    shell_l = lapse / areal
    shell_q = radial / areal
    shell_d = sp.diff(shell_q, time) / shell_l
    radial_velocity = sp.diff(areal, time)
    einstein_F = -4 * _PI * mass * areal**2
    raw_Jt = 24 * _PI * mass * radial * areal / lapse * radial_velocity
    raw_counter = -raw_Jt + 2 * einstein_F * shell_d
    ghy_primitive = -8 * _PI * mass * (
        areal**2 * sp.diff(radial, time) / lapse + 2 * radial * areal * radial_velocity / lapse
    )
    return {
        "assumptions": dict(BOUNDARY_ASSUMPTIONS),
        "euler_cap_cancellation_is_zero": sp.simplify(euler_cap) == 0,
        "myers_Q_plus_is_minus_4pi_jt": sp.simplify(Q_plus + 4 * _PI * jt) == 0,
        "box_cap_cancellation_is_zero": sp.simplify(box_cap) == 0,
        "removed_F_Rh_time_flux": removed_time,
        "restored_F_Rh_time_flux_is_zero_sum": sp.simplify(removed_time + restored_time) == 0,
        "weyl_extra_time_flux_is_4_alpha_chi_D": weyl_extra == 0,
        "ghy_plus_weyl_rebuilds_removed_flux": rebuild == 0,
        "raw_eh_counterterm_matches_ghy_primitive": sp.simplify(raw_counter - ghy_primitive) == 0,
        "dirichlet_cap_form_is_zero": dirichlet_cap == 0,
        "periodic_spatial_flux_is_zero": sp.simplify(periodic_flux) == 0,
        "spatial_corners": False,
        "pure_einstein_chart": False,
    }


def _slice_derivative(value, coordinate):
    if isinstance(value, sp.Expr) and coordinate in value.free_symbols:
        return sp.diff(value, coordinate)
    return sp.Integer(0)


def conditional_initial_slice(
    rho, j, p_Q, *, sign, r, Q, chi, p_chi, A, C_W, C_F, flux, coordinate,
):
    """Conditional slice with constant r, Q, chi, p_chi and periodic mean-zero j.

    K = p_Q p_chi/(2 F_chi) - Q V + rho, Pi = sign sqrt(-4 Z Q K),
    p_r = (F_r/F_chi) p_chi + Pi, and p_Q_x = j/Q.
    Then C_geom + rho = 0 and D_geom + j = 0. This is not a physical source
    covariance, not an imposed future pulse, and not the separate source-coupling branch.
    """
    require_weyl_chart(C_W, A)
    if sign not in (1, -1, sp.Integer(1), sp.Integer(-1)):
        raise ValueError("Pi sign must be +1 or -1")
    F_r, f = partial_F(r, A, C_W)
    if f == 0:
        raise ValueError("F_chi vanishes; C_W=0 is rejected")
    Z = feedback_Z(A)
    potential = feedback_V(r, chi, A, C_W, C_F, flux)
    K = p_Q * p_chi / (2 * f) - Q * potential + rho
    discriminant = sp.simplify(-4 * Z * Q * K)
    Pi = sign * sp.sqrt(discriminant)
    p_r = (F_r / f) * p_chi + Pi
    r_x = _slice_derivative(r, coordinate)
    chi_x = _slice_derivative(chi, coordinate)
    Q_x = _slice_derivative(Q, coordinate)
    r_xx = _slice_derivative(r_x, coordinate)
    chi_xx = _slice_derivative(chi_x, coordinate)
    p_Q_x = _slice_derivative(p_Q, coordinate)
    geometric_C = lapse_constraint(
        p_Q, p_r, p_chi, Q, r, chi, r_x, r_xx, chi_x, chi_xx, Q_x, A, C_W, C_F, flux,
    )
    geometric_D = shift_constraint(p_Q_x, p_r, p_chi, Q, r_x, chi_x)
    return {
        "K": sp.simplify(K),
        "Pi": Pi,
        "p_r": p_r,
        "p_Q_x_minus_j_over_Q": sp.simplify(p_Q_x - j / Q),
        "lapse_residual": sp.simplify(geometric_C + rho),
        "shift_residual": sp.simplify(geometric_D + j),
        "physical_source_covariance": False,
        "imposed_future_pulses": False,
        "source_coupling_owned": False,
    }


def normalized_conditional_examples():
    """Two conditional slices at the reviewed normalization. Both Pi signs."""
    coordinate = sp.symbols("x")
    coefficients = {
        "A": 1 / (4 * _PI),
        "C_W": 3 / (4 * _PI),
        "C_F": 1 / _PI,
        "flux": sp.Integer(1),
    }
    rho = 1 + sp.cos(coordinate) / 4
    families = (
        {"j": sp.Integer(0), "p_Q": sp.Integer(2), "K_expected": 2 + sp.cos(coordinate) / 4,
         "minimum": sp.Rational(7, 4)},
        {"j": sp.sin(coordinate) / 5, "p_Q": 2 - sp.cos(coordinate) / 5,
         "K_expected": 2 + 3 * sp.cos(coordinate) / 20, "minimum": sp.Rational(37, 20)},
    )
    rows = []
    for family in families:
        for sign in (1, -1):
            built = conditional_initial_slice(
                rho, family["j"], family["p_Q"], sign=sign, r=sp.Integer(1), Q=sp.Integer(1),
                chi=sp.Integer(0), p_chi=sp.Integer(-2), coordinate=coordinate, **coefficients,
            )
            gap = sp.simplify(built["K"] - family["minimum"])
            rows.append({
                **built,
                "j": family["j"],
                "coordinate": coordinate,
                "K_matches_review": sp.simplify(built["K"] - family["K_expected"]) == 0,
                "minimum": family["minimum"],
                "K_minus_minimum": gap,
                "pi_square": sp.simplify(sp.together(built["Pi"]**2)),
            })
    return rows


def continuum_matter_hamiltonian():
    """Direct conformal continuum operator on u = r sqrt(q) psi. Not the discrete sqrt(N) matrix."""
    return {
        "frame": "u=r*sqrt(q)*psi, covariance C held fixed",
        "operator": "sigma2/2*{L/Q,P}+sigma1*L*kappa-1/2*{beta,P}",
        "momentum": "P=-i*dx",
        "old_matrix_sqrtN_ordering_fixed": False,
        "reason": "existing sqrt(N) ordering fails the finite-grid conformal Ward identity; a separate owner has the source operator",
        "isotropic_copy_weight": "4*kappa multiplies one angular block and is not recomputed here",
        "radius_force_at_fixed_L_Q": "cancels in this continuum chart only",
    }


def continuum_radius_force_residual(L, Q, r, kappa):
    """Coefficient cancellation of δH/δr at fixed L, Q. Continuum symbols only."""
    N = L * r
    q = Q * r
    delta_r = sp.symbols("delta_r")
    delta_N = L * delta_r
    delta_q = Q * delta_r
    sigma2 = delta_N / q - N * delta_q / q**2
    sigma1 = kappa * delta_N / r - N * kappa * delta_r / r**2
    return {
        "sigma2_coefficient": sp.simplify(sigma2),
        "sigma1_coefficient": sp.simplify(sigma1),
        "both_zero": sp.simplify(sigma2) == 0 and sp.simplify(sigma1) == 0,
    }


def curvature_limits():
    """Schwarzschild, de Sitter and flat limits from the 2D formula only."""
    x, M, ell = sp.symbols("x M ell", positive=True)
    t = sp.symbols("t")

    def pack(N, q, radius):
        beta = sp.Integer(0)
        R_h = areal_curvature_scalar(N, q, radius, beta, t, x)
        return R_h, sp.simplify(weyl_square_from_rh(R_h, radius))

    f_s = 1 - 2 * M / x
    R_s, C_s = pack(sp.sqrt(f_s), 1 / sp.sqrt(f_s), x)
    f_d = 1 - x**2 / ell**2
    R_d, C_d = pack(sp.sqrt(f_d), 1 / sp.sqrt(f_d), x)
    R_f, C_f = pack(sp.Integer(1), sp.Integer(1), x)
    return {
        "schwarzschild_C2_is_48M2_over_r6": sp.simplify(C_s - 48 * M**2 / x**6) == 0,
        "de_sitter_Rh_is_2": sp.simplify(R_d - 2) == 0,
        "de_sitter_C2_is_0": sp.simplify(C_d) == 0,
        "flat_C2_is_0": sp.simplify(C_f) == 0,
        "schwarzschild_C2": C_s,
        "de_sitter_Rh": R_d,
    }


def bulk_status():
    """Checked bulk identities and the boundary flags still open."""
    legendre = legendre_identity()
    signs = constraint_sign_identities()
    hamilton = hamiltonian_right_hand_sides()
    auxiliary = auxiliary_square_identity()
    ghy = homogeneous_ghy_match()
    box = box_flux_identity()
    euler = euler_candidate_status()
    limits = curvature_limits()
    matter = continuum_matter_hamiltonian()
    boundary = boundary_completion_ledger()
    initial = normalized_conditional_examples()
    return {
        "cauchy_evolution_owned": CAUCHY_EVOLUTION_OWNED,
        "pure_einstein_chart": PURE_EINSTEIN_CHART,
        "coefficients": "F=-4*pi*A*r**2+2*alpha*chi; Z=-24*pi*A; V=8*pi*A*r**2-alpha*(4*chi+chi**2)-2*pi*C_F*flux**2; alpha=-4*pi*C_W/3",
        "legendre_identity": legendre["residual_is_zero"],
        "constraint_signs": signs["lapse_el_plus_C_is_zero"] and signs["shift_el_plus_Dcon_is_zero"],
        "constraint_sign_convention": signs["lapse_sign"] + "; " + signs["shift_sign"],
        "hamilton_velocity_rhs": hamilton["velocity_rhs_matches_inverse"],
        "hamilton_momentum_rhs": hamilton["momentum_residuals_are_zero"],
        "auxiliary_elimination": auxiliary["restored_is_zero"] and auxiliary["critical_point_is_curvature"],
        "homogeneous_ghy": ghy["residual_is_zero"],
        "box_r_divergence": box["divergence_matches_reduced_box"],
        "box_r_homogeneous_endpoint": box["homogeneous_endpoint_matches_owner"],
        "box_cap_cancellation": box["cap_cancellation_is_zero"],
        "euler_cap_cancellation": boundary["euler_cap_cancellation_is_zero"],
        "raw_eh_matches_ghy": boundary["raw_eh_counterterm_matches_ghy_primitive"],
        "weyl_extra_time_flux": boundary["weyl_extra_time_flux_is_4_alpha_chi_D"],
        "conditional_initial_residuals": all(
            row["lapse_residual"] == 0 and row["shift_residual"] == 0 and row["K_matches_review"]
            for row in initial
        ),
        "euler_current": euler["status"],
        "euler_naive_rejected": euler["naive_current_is_distinct"] and not EULER_NAIVE_CURRENT_ALLOWED,
        "schwarzschild": limits["schwarzschild_C2_is_48M2_over_r6"],
        "de_sitter": limits["de_sitter_Rh_is_2"] and limits["de_sitter_C2_is_0"],
        "flat": limits["flat_C2_is_0"],
        "matter_operator": matter["operator"],
        "old_matrix_source_fixed": matter["old_matrix_sqrtN_ordering_fixed"],
        "source_coupling_owned": SOURCE_COUPLING_OWNED,
        "spatial_corners": boundary["spatial_corners"],
        "open_boundary_flags": (
            "x is periodic and the time interval has finite caps; Dirichlet data there are Q, r, chi.",
            "Spatial corner terms are absent on this domain and are not an extra requirement.",
            "Euler, Box-R and the removed F R_h time flux cancel by the cap counterterms; Einstein stays the GHY completion plus the Weyl flux 4 alpha chi D.",
            "Conditional initial slices are not a physical source covariance and are not evolved.",
            "The continuum matter operator is not coupled; the old discrete sqrt(N) source is not fixed.",
            "No discrete constraint preservation, so the six geometric equations are not a Cauchy evolution.",
        ),
    }


@dataclass(frozen=True)
class FeedbackChart:
    """Numeric coefficient chart. C_W = 0 and A = 0 are rejected."""

    A: float
    C_W: float
    C_F: float
    flux: float

    def __post_init__(self):
        for name in ("A", "C_W", "C_F", "flux"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
                raise ValueError("finite real chart coefficients required")
        if self.C_W == 0:
            raise ValueError(
                "C_W=0 is rejected; this chart keeps the Weyl auxiliary and does not implement pure Einstein"
            )
        if self.A == 0:
            raise ValueError("A=0 makes Z vanish; inverse velocities are outside this chart")
