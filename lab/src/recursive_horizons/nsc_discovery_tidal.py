"""Observer tidal field of the realized conformal metric.

The line element is the owned (+---) chart

    ds^2 = N^2 (dt^2 - dx^2) - r^2 dOmega^2,
    N = r Q,  L = Q,  beta = 0.

No new field or geometry equation is integrated. Time jets are the
directional derivative of the existing projected Hamilton rate along that
same rate. Spatial jets are the owned Nyquist-zero Fourier derivative.
Chi is never inserted into a curvature component.

Riemann convention, matched to ``nsc_finite_terms`` and to a direct
contraction of this metric::

    R^a_bcd = d_c Gamma^a_db - d_d Gamma^a_cb
              + Gamma^a_ce Gamma^e_db - Gamma^a_de Gamma^e_cb
    Ricci_bd = R^a_bad
    R = g^{bd} Ricci_bd

The static sections of that reference reduce to the conformal formulas
below when time derivatives vanish. The candidate scalar

    R4 = (R_h - 2) / r^2 - 6 (r_tt - r_xx) / (r^3 Q^2)

is that contraction, not a separate convention.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np
from threadpoolctl import threadpool_limits

from . import nsc_discovery_episode as episode
from . import nsc_nested_parent_child as nested
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_episode_assessment as metric
from . import nsc_spherical_feedback_action as action
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-TIDAL-v1"
SINGULARITY_DECLARED = False
THEORY_DECLARED_DEAD = False
CHI_SUBSTITUTED = False
NEW_ODE = False
PRODUCTION_RECORD_WRITTEN_BY_LIBRARY = False
PRIMARY_CASE = "nf256_coupled_dt0.0005"
COARSE_CASE = "nf128_coupled_dt0.0005"
TIMESTEP_CASE = "nf256_coupled_dt0.001"
FROZEN_CASE = "nf256_frozen_geometry_dt0.0005"
STATION_TIMES = (0.3, 1.0, 3.0)
WORLD_LINE = 2.0
CLOCK_LOCATIONS = (1.0, 2.0, 3.0)
FROZEN_PRODUCING_COMMITS = {"episode": "1c9e770", "confirmation": "528efc7"}
_PHYSICS_NAMES = (
    "nsc_nested_parent_child.py", "nsc_spherical_coupling.py",
    "nsc_spherical_galerkin_coupling.py", "nsc_spherical_feedback_action.py",
    "nsc_spherical_episode_assessment.py", "nsc_conformal_adm_source.py",
    "nsc_covariant_operator.py", "nsc_nested_qualities.py", "nsc_adm_source.py",
)

EPISODE_DIR = episode.LAB / "results" / "development" / "nsc-discovery-episode-v1"
CONFIRMATION_DIR = episode.LAB / "results" / "development" / "nsc-discovery-confirmation-v1"

RIEMANN_CONVENTION = (
    "R^a_bcd = d_c Gamma^a_db - d_d Gamma^a_cb "
    "+ Gamma^a_ce Gamma^e_db - Gamma^a_de Gamma^e_cb; "
    "Ricci_bd = R^a_bad; signature (+---)"
)
FORMULA_RH = "R_h = 2/Q^2 [(Q_xx/Q - Q_x^2/Q^2) - (Q_tt/Q - Q_t^2/Q^2)]"
FORMULA_R4 = "R4 = (R_h - 2)/r^2 - 6 (r_tt - r_xx)/(r^3 Q^2)"
FORMULA_R0101 = (
    "R_0101 = [(r_tt - r_xx)/r - (r_t^2 - r_x^2)/r^2 "
    "+ (Q_tt/Q - Q_t^2/Q^2) - (Q_xx/Q - Q_x^2/Q^2)] /(Q^2 r^2)"
)
FORMULA_R0202 = (
    "R_0202 = R_0303 = (r r_tt - r_t^2 - r_x^2)/(Q^2 r^4) "
    "- (Q_t r_t + Q_x r_x)/(Q^3 r^3); R_0102 = 0"
)
FORMULA_OWNED_W = "C^2 = (R_h - 2)^2 / (3 r^4)"
FORMULA_WEYL_ALGEBRA = "K - 2 Ricci^2 + R^2/3"
FORMULA_ACCELERATION = (
    "Q_t = P[Q p_chi / (2 F_chi)], "
    "Q_tt = P[(Q_t p_chi + Q p_chi_t)/(2 F_chi)], "
    "r_t = P[pi/(2 Z)], pi = p_r - F_r p_chi/F_chi, F_r = -8 pi A r, "
    "r_tt = P[d/dt(pi)/(2 Z)] along the same projected rate"
)

_PRODUCER_MODULES = (
    episode,
    nested,
    coupling,
    metric,
    action,
    galerkin,
)


def sha256_file(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def producer_hashes(extra_paths=()):
    """Hashes of the modules that define the rate, the derivative, and this consumer."""
    paths = [Path(__file__).resolve()]
    for module in _PRODUCER_MODULES:
        paths.append(Path(module.__file__).resolve())
    paths.extend(Path(__file__).parent / name for name in _PHYSICS_NAMES)
    paths.extend(Path(path).resolve() for path in extra_paths)
    return {str(path): sha256_file(path) for path in paths}


def _positive(values, name):
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.isfinite(array).all() or np.min(array) <= 0.0:
        raise ValueError(name + " must be finite and positive")
    return array


def _as_float(values):
    array = np.asarray(values, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError("curvature jets must be finite")
    return array


def conformal_curvature_jets(Q, r, Qt, Qx, Qtt, Qxx, rt, rx, rtt, rxx, *, Qtx=None, rtx=None, Nx=None):
    """Pointwise 4D curvature of N = r Q, beta = 0.

    Derivatives are supplied by the caller. This function does not difference
    a trajectory and does not read chi. ``Nx`` is the owned spatial derivative
    of the metric factor N when the caller has one; otherwise the product rule
    in x is used and the gap is zero by construction.
    """
    Q = _positive(Q, "Q")
    r = _positive(r, "r")
    Qt, Qx = _as_float(Qt), _as_float(Qx)
    Qtt, Qxx = _as_float(Qtt), _as_float(Qxx)
    rt, rx = _as_float(rt), _as_float(rx)
    rtt, rxx = _as_float(rtt), _as_float(rxx)
    shape = Q.shape
    for array, name in (
        (r, "r"), (Qt, "Q_t"), (Qx, "Q_x"), (Qtt, "Q_tt"), (Qxx, "Q_xx"),
        (rt, "r_t"), (rx, "r_x"), (rtt, "r_tt"), (rxx, "r_xx"),
    ):
        if array.shape != shape:
            raise ValueError(name + " does not match Q")

    bracket_t = Qtt / Q - Qt ** 2 / Q ** 2
    bracket_x = Qxx / Q - Qx ** 2 / Q ** 2
    raw_qtt = 2.0 * Qtt / Q ** 3
    raw_qt2 = 2.0 * Qt ** 2 / Q ** 4
    raw_qxx = 2.0 * Qxx / Q ** 3
    raw_qx2 = 2.0 * Qx ** 2 / Q ** 4
    term_t = raw_qtt - raw_qt2
    term_x = raw_qxx - raw_qx2
    R_h = term_x - term_t

    shell = (R_h - 2.0) / r ** 2
    acceleration = -6.0 * (rtt - rxx) / (r ** 3 * Q ** 2)
    R4 = shell + acceleration

    r_accel = (rtt - rxx) / r - (rt ** 2 - rx ** 2) / r ** 2
    q_wave = bracket_t - bracket_x
    R0101 = (r_accel + q_wave) / (Q ** 2 * r ** 2)
    angular_principal = (r * rtt - rt ** 2 - rx ** 2) / (Q ** 2 * r ** 4)
    angular_shift = (Qt * rt + Qx * rx) / (Q ** 3 * r ** 3)
    R0202 = angular_principal - angular_shift
    R0303 = R0202
    ricci_normal = -(R0101 + R0202 + R0303)

    N = r * Q
    Nt = rt * Q + r * Qt
    Ntt = rtt * Q + 2.0 * rt * Qt + r * Qtt
    Nx_product = rx * Q + r * Qx
    if Nx is None:
        Nx_metric = Nx_product
        product_gap = np.zeros(shape)
    else:
        Nx_metric = _as_float(Nx)
        if Nx_metric.shape != shape:
            raise ValueError("N_x does not match Q")
        product_gap = Nx_metric - Nx_product
    gamma_scalar = np.full(shape, 2.0)
    return {
        "Q": Q,
        "r": r,
        "N": N,
        "R_h": R_h,
        "R4": R4,
        "shell": shell,
        "acceleration": acceleration,
        "R_0101": R0101,
        "R_0202": R0202,
        "R_0303": R0303,
        "R_0102": np.zeros(shape),
        "ricci_normal_00": ricci_normal,
        "owned_W": (R_h - 2.0) ** 2 / (3.0 * r ** 4),
        "term_t": term_t,
        "term_x": term_x,
        "raw_qtt": raw_qtt,
        "raw_qt2": raw_qt2,
        "raw_qxx": raw_qxx,
        "raw_qx2": raw_qx2,
        "angular_principal": angular_principal,
        "angular_shift": angular_shift,
        "gamma_scalar": gamma_scalar,
        "N_t": Nt,
        "N_tt": Ntt,
        "N_x_product": Nx_product,
        "product_gap": product_gap,
        "bracket_t": bracket_t,
        "bracket_x": bracket_x,
    }


def base_geometry(N, r, Nt, Nx, Ntt, Ntx, Nxx, rt, rx, rtt, rtx, rxx):
    """Warped-product pieces of g2 = N^2 diag(1, -1) with fibre radius r.

    ``gamma_scalar`` of the unit sphere is 2. Its contribution to the 4D
    scalar is ``-2 (1 + |grad r|^2) / r^2``. The covariant hessian is
    ``∇_a ∇_b r`` on the Lorentzian base. These are the owned
    ``spherical_invariants`` ingredients specialized to beta = 0 and N = q.
    """
    N = _positive(N, "N")
    r = _positive(r, "r")
    Nt, Nx = _as_float(Nt), _as_float(Nx)
    Ntt, Ntx, Nxx = _as_float(Ntt), _as_float(Ntx), _as_float(Nxx)
    rt, rx = _as_float(rt), _as_float(rx)
    rtt, rtx, rxx = _as_float(rtt), _as_float(rtx), _as_float(rxx)
    bracket_t = Ntt / N - Nt ** 2 / N ** 2
    bracket_x = Nxx / N - Nx ** 2 / N ** 2
    R2 = 2.0 / N ** 2 * (bracket_x - bracket_t)
    htt = rtt - (Nt / N) * rt - (Nx / N) * rx
    hxx = rxx - (Nt / N) * rt - (Nx / N) * rx
    htx = rtx - (Nx / N) * rt - (Nt / N) * rx
    box = (htt - hxx) / N ** 2
    grad2 = (rt ** 2 - rx ** 2) / N ** 2
    hess2 = (htt ** 2 + hxx ** 2 - 2.0 * htx ** 2) / N ** 4
    gamma_piece = -2.0 * (1.0 + grad2) / r ** 2
    scalar = R2 - 4.0 * box / r + gamma_piece
    ricci2 = (
        R2 ** 2 / 2.0
        - 2.0 * R2 * box / r
        + 4.0 * hess2 / r ** 2
        + 2.0 * (1.0 + grad2 + r * box) ** 2 / r ** 4
    )
    kretschmann = R2 ** 2 + 8.0 * hess2 / r ** 2 + 4.0 * (1.0 + grad2) ** 2 / r ** 4
    assembled = kretschmann - 2.0 * ricci2 + scalar ** 2 / 3.0
    A = R2 / 2.0
    B00, B11, B01 = htt / (r * N ** 2), hxx / (r * N ** 2), htx / (r * N ** 2)
    S = (1.0 + grad2) / r ** 2
    return {
        "A": A, "B00": B00, "B11": B11, "B01": B01, "S": S,
        "normal_tidal_radial": -A,
        "normal_tidal_angular": B00,
        "normal_worldline_acceleration_radial": Nx / N ** 2,
        "angular_tidal_condition_numerator": (
            np.abs(rtt) + np.abs(Nt * rt / N) + np.abs(Nx * rx / N)
        ) / (r * N ** 2),
        "R4_condition_numerator": np.abs(2*A) + np.abs(4*(B00-B11)) + np.abs(2*S),
        "factorized_W": 4.0 / 3.0 * (A + B00 - B11 - S) ** 2,
        "R2": R2,
        "box": box,
        "grad2": grad2,
        "hess_tt": htt,
        "hess_xx": hxx,
        "hess_tx": htx,
        "hess2": hess2,
        "gamma_scalar": np.full(N.shape, 2.0),
        "gamma_piece": gamma_piece,
        "R4": scalar,
        "Ricci2": ricci2,
        "K": kretschmann,
        "assembled_W": assembled,
    }


def owner_static_sections(N, q, r, Nx, Nxx, qx, rx, rxx):
    """Owned static (+---) sections. Time derivatives are absent.

    A = D(D(N))/N, B = D(N) D(r)/(N r), C = D(D(r))/r, Dsec = (1 - D(r)^2)/r^2
    with D = q^{-1} d_x. Returns R, Ricci^2 and Kretschmann.
    """
    N = _positive(N, "N")
    q = _positive(q, "q")
    r = _positive(r, "r")
    Nx, Nxx = _as_float(Nx), _as_float(Nxx)
    qx, rx, rxx = _as_float(qx), _as_float(rx), _as_float(rxx)
    A = (Nxx - Nx * qx / q) / (N * q ** 2)
    B = Nx * rx / (N * q ** 2 * r)
    C = (rxx - rx * qx / q) / (q ** 2 * r)
    Dsec = (1.0 - rx ** 2 / q ** 2) / r ** 2
    scalar = 2.0 * (A + 2.0 * B + 2.0 * C - Dsec)
    ricci2 = (A + 2.0 * B) ** 2 + (A + 2.0 * C) ** 2 + 2.0 * (B + C - Dsec) ** 2
    kretschmann = 4.0 * (A ** 2 + 2.0 * B ** 2 + 2.0 * C ** 2 + Dsec ** 2)
    return {
        "A": A, "B": B, "C": C, "D": Dsec,
        "R": scalar, "Ricci2": ricci2, "K": kretschmann,
        "C2": 4.0 * (A - B - C - Dsec) ** 2 / 3.0,
    }


def _zeros_like(Q):
    return np.zeros(np.asarray(Q).shape, dtype=float)


def curvature_from_local_jets(Q, r, Qt, Qx, Qtt, Qtx, Qxx, rt, rx, rtt, rtx, rxx):
    """Full conformal package from local jets. No periodic derivative is applied."""
    Q = _positive(Q, "Q")
    r = _positive(r, "r")
    Qt, Qx = _as_float(Qt), _as_float(Qx)
    Qtt, Qtx, Qxx = _as_float(Qtt), _as_float(Qtx), _as_float(Qxx)
    rt, rx = _as_float(rt), _as_float(rx)
    rtt, rtx, rxx = _as_float(rtt), _as_float(rtx), _as_float(rxx)
    tides = conformal_curvature_jets(Q, r, Qt, Qx, Qtt, Qxx, rt, rx, rtt, rxx)
    N = tides["N"]
    Nt = tides["N_t"]
    Ntt = tides["N_tt"]
    Nx = rx * Q + r * Qx
    Nxx = rxx * Q + 2.0 * rx * Qx + r * Qxx
    Ntx = rtx * Q + rt * Qx + rx * Qt + r * Qtx
    base = base_geometry(N, r, Nt, Nx, Ntt, Ntx, Nxx, rt, rx, rtt, rtx, rxx)
    static = owner_static_sections(N, N, r, Nx, Nxx, Nx, rx, rxx)
    return {"tides": tides, "base": base, "static_sections_at_frozen_time": static}


def analytic_accelerations(pair, state, control_mode):
    """Second time jets of the controller's actual projected rate.

    Frozen geometry passes the controller rate, whose geometry components are
    zero, so ``Q_t = r_t = Q_tt = r_tt = 0``. The field column rate is left
    intact. Coupled evolution uses the projected bilinear for ``Q`` and the
    same projection of ``d/dt(pi/(2 Z))`` for ``r``. ``F_chi`` and ``Z`` are
    constants on this chart; ``F_r = -8 pi A r`` is differentiated.
    """
    nodal, rate, bundle, mode = episode.actual_control_rate(pair, state, control_mode)
    fine = bundle["fine_state"]
    grid = pair.grid
    if grid.gauge != "conformal":
        raise ValueError("tidal jets use the conformal chart L=Q, beta=0")
    area = float(grid.fine.A)
    _force_r, f_chi = action.partial_F(fine.r, area, grid.fine.C_W)
    f_chi = float(f_chi)
    if f_chi == 0.0:
        raise ValueError("F_chi vanished; the Weyl chart is singular")
    Z = float(action.feedback_Z(area))
    if Z == 0.0:
        raise ValueError("Z vanished")
    q_dot = galerkin.prolong_geometry(grid, rate.Q)
    r_dot = galerkin.prolong_geometry(grid, rate.r)
    p_r_dot = galerkin.prolong_geometry(grid, rate.p_r)
    p_chi_dot = galerkin.prolong_geometry(grid, rate.p_chi)
    directional_q = (q_dot * fine.p_chi + fine.Q * p_chi_dot) / (2.0 * f_chi)
    q_ddot = galerkin.prolong_geometry(grid, galerkin.pull_geometry(grid, directional_q))
    force_r = -8.0 * math.pi * area * fine.r
    force_along = -8.0 * math.pi * area * r_dot
    pi_state = fine.p_r - force_r * fine.p_chi / f_chi
    unprojected_r = pi_state / (2.0 * Z)
    pi_along = p_r_dot - (force_along * fine.p_chi + force_r * p_chi_dot) / f_chi
    directional_r = pi_along / (2.0 * Z)
    r_ddot = galerkin.prolong_geometry(grid, galerkin.pull_geometry(grid, directional_r))
    if mode == "frozen_geometry":
        geometry_rate = max(
            float(np.max(np.abs(q_dot))),
            float(np.max(np.abs(r_dot))),
            float(np.max(np.abs(q_ddot))),
            float(np.max(np.abs(r_ddot))),
            float(np.max(np.abs(p_r_dot))),
            float(np.max(np.abs(p_chi_dot))),
        )
        if geometry_rate != 0.0:
            raise ValueError("frozen controller returned a nonzero geometry jet")
    return {
        "mode": mode,
        "nodal": nodal,
        "rate": rate,
        "bundle": bundle,
        "fine": fine,
        "Q_t": q_dot,
        "Q_tt": q_ddot,
        "Q_tt_before_projection": directional_q,
        "r_t": r_dot,
        "r_tt": r_ddot,
        "r_tt_before_projection": directional_r,
        "unprojected_r_formula": unprojected_r,
        "unprojected_r_owner": np.asarray(bundle["unprojected_r"], dtype=float),
        "field_rate_max": float(np.max(np.abs(rate.phi0))),
        "F_chi": f_chi,
        "Z": Z,
        "A": area,
    }


def _summary(values, spacing):
    array = np.asarray(values, dtype=float)
    return {
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "max_abs": float(np.max(np.abs(array))),
        "l2": float(np.sqrt(float(spacing) * np.sum(array ** 2))),
    }


def _condition(numerator, residual):
    num = float(np.max(np.abs(np.asarray(numerator, dtype=float))))
    res = float(np.max(np.abs(np.asarray(residual, dtype=float))))
    return {
        "numerator_max_abs": num,
        "residual_max_abs": res,
        "condition": num / max(res, 1e-300),
    }


def _nearest(coordinate, target):
    position = np.asarray(coordinate, dtype=float)
    index = int(np.argmin(np.abs(position - float(target))))
    return index, float(position[index]), float(abs(position[index] - float(target)))


def profiles_on_grid(pair, jets):
    """Fourier spatial jets of the realized fine metric, then the curvature package."""
    fine = jets["fine"]
    grid = pair.grid
    period = float(grid.length)
    spacing = float(grid.dx_q)
    Q = np.asarray(fine.Q, dtype=float)
    r = np.asarray(fine.r, dtype=float)
    derivative = lambda values: metric.spectral_dx(values, period)
    Q_t = np.asarray(jets["Q_t"], dtype=float)
    r_t = np.asarray(jets["r_t"], dtype=float)
    Q_tt = np.asarray(jets["Q_tt"], dtype=float)
    r_tt = np.asarray(jets["r_tt"], dtype=float)
    Q_x = derivative(Q)
    r_x = derivative(r)
    Q_xx = derivative(Q_x)
    r_xx = derivative(r_x)
    Q_tx = derivative(Q_t)
    r_tx = derivative(r_t)
    N = r * Q
    N_x = derivative(N)
    N_xx = derivative(N_x)
    N_t = r_t * Q + r * Q_t
    N_tx = derivative(N_t)
    N_tt = r_tt * Q + 2.0 * r_t * Q_t + r * Q_tt
    tides = conformal_curvature_jets(
        Q, r, Q_t, Q_x, Q_tt, Q_xx, r_t, r_x, r_tt, r_xx, Nx=N_x,
    )
    base = base_geometry(N, r, N_t, N_x, N_tt, N_tx, N_xx, r_t, r_x, r_tt, r_tx, r_xx)
    zero = np.zeros(grid.nq, dtype=float)
    owned_direct = metric.direct_rh_grid(
        Q, Q, zero, Q_t, Q_tt, period, L_dot=Q_t, beta_dot=zero,
    )["R_h"]
    dense_Q_x = grid.derivative @ Q
    chi = np.asarray(fine.chi, dtype=float)
    proxy = metric.weyl_proxy_from_chi(chi, r)
    shell_gap = tides["R_h"] - chi - 2.0
    return {
        "spacing": spacing,
        "period": period,
        "x": np.asarray(grid.xi_q, dtype=float),
        "tides": tides,
        "base": base,
        "direct_R_h_gap": owned_direct - tides["R_h"],
        "dense_versus_fourier_Q_x": dense_Q_x - Q_x,
        "chi": chi,
        "chi_proxy": proxy,
        "shell_gap": shell_gap,
        "algebra_gap": base["assembled_W"] - tides["owned_W"],
        "scalar_route_gap": base["R4"] - tides["R4"],
        "Q_t_max": float(np.max(np.abs(Q_t))),
        "Q_tt_max": float(np.max(np.abs(Q_tt))),
        "r_t_max": float(np.max(np.abs(r_t))),
        "r_tt_max": float(np.max(np.abs(r_tt))),
        "r_formula_gap": float(np.max(np.abs(
            jets["unprojected_r_formula"] - jets["unprojected_r_owner"]
        ))),
        "projection_gap_Q": float(np.max(np.abs(Q_tt - jets["Q_tt_before_projection"]))),
        "projection_gap_r": float(np.max(np.abs(r_tt - jets["r_tt_before_projection"]))),
    }


def _at(values, index):
    return float(np.asarray(values, dtype=float)[index])


def station_observables(pair, jets, clocks):
    """Scalars, the x=2 worldline, and the angular-tide peak. Chi stays auxiliary."""
    profile = profiles_on_grid(pair, jets)
    tides, base = profile["tides"], profile["base"]
    spacing = profile["spacing"]
    x = profile["x"]
    world_index, world_x, world_miss = _nearest(x, WORLD_LINE)
    peak_index = int(np.argmax(np.abs(tides["R_0202"])))
    names = {
        "R_h": tides["R_h"],
        "R4": tides["R4"],
        "owned_W": tides["owned_W"],
        "R_0101": tides["R_0101"],
        "R_0202": tides["R_0202"],
        "ricci_normal_00": tides["ricci_normal_00"],
        "R2": base["R2"],
        "Ricci2": base["Ricci2"],
        "K": base["K"],
        "grad2": base["grad2"],
        "hess2": base["hess2"],
        "box": base["box"],
        "gamma_piece": base["gamma_piece"],
        "assembled_W": base["assembled_W"],
        "factorized_W": base["factorized_W"],
        "normal_worldline_acceleration_radial": base["normal_worldline_acceleration_radial"],
    }
    clock_index = CLOCK_LOCATIONS.index(WORLD_LINE)
    clock_values = np.asarray(clocks, dtype=float).reshape(-1)
    if clock_values.size < 3:
        raise ValueError("station clocks must carry the worldlines x=1, 2, 3")

    def point(index):
        return {
            "x": float(x[index]),
            **{name: _at(values, index) for name, values in names.items()},
            "Q": _at(tides["Q"], index),
            "r": _at(tides["r"], index),
            "N": _at(tides["N"], index),
            "angular_tidal_condition": _at(base["angular_tidal_condition_numerator"], index)
                / max(abs(_at(tides["R_0202"], index)), 1e-300),
            "R4_condition": _at(base["R4_condition_numerator"], index)
                / max(abs(_at(base["R4"], index)), 1e-300),
        }

    conditioning = {
        "temporal_raw": _condition(tides["raw_qtt"], tides["term_t"]),
        "temporal_versus_spatial": _condition(tides["term_t"], tides["R_h"]),
        "full_Rh": _condition(tides["raw_qtt"], tides["R_h"]),
        "R4_shell_versus_acceleration": _condition(tides["shell"], tides["R4"]),
        "angular_tide": _condition(
            base["angular_tidal_condition_numerator"], tides["R_0202"],
        ),
        "weyl_assembly": _condition(base["K"], base["assembled_W"]),
        "direct_Rh_gap_max": float(np.max(np.abs(profile["direct_R_h_gap"]))),
        "fourier_versus_dense_Qx_max": float(np.max(np.abs(profile["dense_versus_fourier_Q_x"]))),
        "product_rule_Nx_gap_max": float(np.max(np.abs(tides["product_gap"]))),
        "scalar_route_gap_max": float(np.max(np.abs(profile["scalar_route_gap"]))),
        "algebra_gap_max": float(np.max(np.abs(profile["algebra_gap"]))),
        "tidal_route_gap_max": float(max(
            np.max(np.abs(base["normal_tidal_radial"] - tides["R_0101"])),
            np.max(np.abs(base["normal_tidal_angular"] - tides["R_0202"])),
        )),
        "r_formula_gap_max": profile["r_formula_gap"],
        "projection_gap_Q_max": profile["projection_gap_Q"],
        "projection_gap_r_max": profile["projection_gap_r"],
    }
    finite = all(np.isfinite(np.asarray(values)).all() for values in names.values())
    chart_positive = bool(np.min(tides["Q"]) > 0.0 and np.min(tides["r"]) > 0.0)
    return {
        "summaries": {name: _summary(values, spacing) for name, values in names.items()},
        "worldline_x2": point(world_index),
        "worldline_miss": world_miss,
        "peak_R0202": point(peak_index),
        "Q_min": float(np.min(tides["Q"])),
        "Q_max": float(np.max(tides["Q"])),
        "r_min": float(np.min(tides["r"])),
        "r_max": float(np.max(tides["r"])),
        "Q_t_max": profile["Q_t_max"],
        "Q_tt_max": profile["Q_tt_max"],
        "r_t_max": profile["r_t_max"],
        "r_tt_max": profile["r_tt_max"],
        "field_rate_max": jets["field_rate_max"],
        "control_mode": jets["mode"],
        "proper_time_x2": float(clock_values[clock_index]),
        "proper_times": [float(value) for value in clock_values[:3]],
        "auxiliary_shell_gap_max": float(np.max(np.abs(profile["shell_gap"]))),
        "auxiliary_proxy_gap_max": float(np.max(np.abs(tides["owned_W"] - profile["chi_proxy"]))),
        "chi_substituted_for_curvature": False,
        "chart_positive": chart_positive,
        "finite": finite,
        "conditioning": conditioning,
        "gamma_scalar": 2.0,
        "_curvature_profiles": {name: values for name, values in names.items()},
    }


def _hash_tree(directory):
    directory = Path(directory)
    hashed = {}
    if not directory.is_dir():
        return hashed
    for path in sorted(directory.iterdir()):
        if path.is_file() and path.suffix in {".json", ".npz", ".jsonl"}:
            hashed[str(path.resolve())] = sha256_file(path)
    return hashed


def _stations(directory, case_id):
    rows = episode.read_station_states(directory, case_id)
    selected = []
    for row in rows:
        time_value = float(row["coordinate_time"])
        if any(abs(time_value - target) <= 1e-8 for target in STATION_TIMES):
            selected.append(row)
    return selected


def measure_case(directory, case_id):
    """Read-only curvature of the saved stations of one case. No step is taken."""
    directory = Path(directory)
    stations = _stations(directory, case_id)
    if not stations:
        raise FileNotFoundError("no T=0.3, 1, 3 station for " + case_id)
    measured = []
    declared_pins = None
    for row in stations:
        record, arrays = episode.load_checkpoint(directory, case_id, int(row["ordinal"]))
        if declared_pins is None:
            declared_pins = dict(record.get("source_pins") or {})
        pair = episode.pair_from_arrays(arrays, record)
        state = episode.state_from_arrays(arrays, record["momentum_representation"])
        jets = analytic_accelerations(pair, state, record.get("control_mode") or record.get("geometry"))
        observables = station_observables(pair, jets, arrays["normal_clocks"])
        observables.update({
            "case_id": case_id,
            "ordinal": int(record["ordinal"]) if "ordinal" in record else int(row["ordinal"]),
            "coordinate_time": float(record["coordinate_time"]),
            "nf": int(record["nf"]),
            "step_cap": float(record["step_cap"]),
            "snapshot_kind": row.get("snapshot_kind"),
            "arrays_sha256": record.get("arrays_sha256"),
            "geometry_evolved_here": False,
        })
        measured.append(observables)
        del pair, state, jets, arrays
    measured.sort(key=lambda item: item["coordinate_time"])
    return {"case_id": case_id, "stations": measured, "declared_source_pins": declared_pins}


def _relative(later, earlier):
    later = float(later)
    earlier = float(earlier)
    scale = max(abs(earlier), 1e-300)
    return {
        "earlier": earlier,
        "later": later,
        "difference": later - earlier,
        "ratio": None if earlier == 0.0 else later / earlier,
        "relative_difference": (later - earlier) / scale,
    }


def _station_map(case):
    return {float(row["coordinate_time"]): row for row in case["stations"]}


def _growth_series(case):
    rows = _station_map(case)
    times = [time for time in STATION_TIMES if time in rows or any(abs(time - key) < 1e-8 for key in rows)]
    # Keys are the stored times. Match targets.
    keyed = {}
    for target in STATION_TIMES:
        for key, row in rows.items():
            if abs(key - target) <= 1e-8:
                keyed[target] = row
    if set(keyed) != set(STATION_TIMES):
        raise ValueError("growth series needs stations 0.3, 1 and 3")
    channels = ("R_0202", "R_0101", "ricci_normal_00", "K", "Ricci2", "owned_W", "R4", "grad2")

    def series(place):
        out = {}
        for channel in channels:
            values = {str(time): keyed[time][place][channel] for time in STATION_TIMES}
            out[channel] = {
                "values": values,
                "T0.3_to_1": _relative(values["1.0"], values["0.3"]),
                "T1_to_3": _relative(values["3.0"], values["1.0"]),
                "T0.3_to_3": _relative(values["3.0"], values["0.3"]),
            }
        clocks = {str(time): keyed[time]["proper_time_x2"] for time in STATION_TIMES}
        locations = {str(time): keyed[time][place]["x"] for time in STATION_TIMES}
        return {"channels": out, "proper_time_x2": clocks, "x": locations}

    clock = {}
    for place in ("worldline_x2", "peak_R0202"):
        block = series(place)["channels"]["R_0202"]
        clock[place] = {}
        for label, left, right in (("T0.3_to_1", 0.3, 1.0), ("T1_to_3", 1.0, 3.0)):
            dt = right - left
            dtau = keyed[right]["proper_time_x2"] - keyed[left]["proper_time_x2"]
            dvalue = keyed[right][place]["R_0202"] - keyed[left][place]["R_0202"]
            clock[place][label] = {
                "coordinate_dt": dt,
                "proper_dtau_x2": dtau,
                "proper_time_source": "stored trapezoid of r Q on the fixed worldline x=2",
                "delta_R0202": dvalue,
                "per_coordinate_time": dvalue / dt,
                "per_proper_time_x2": None if dtau == 0.0 else dvalue / dtau,
            }
    return {
        "worldline_x2": series("worldline_x2"),
        "peak_R0202": series("peak_R0202"),
        "clock_comparison": clock,
    }


def _pair_gap(left, right, channel):
    a = left["summaries"][channel]["max_abs"]
    b = right["summaries"][channel]["max_abs"]
    return _relative(b, a)


def _at_time(case, time_value):
    for row in case["stations"]:
        if abs(row["coordinate_time"] - time_value) <= 1e-8:
            return row
    raise KeyError((case["case_id"], time_value))


def resolution_table(coarser, finer):
    """nf refinement of max-abs and of the fixed worldline x=2. Later minus earlier."""
    channels = ("R_0202", "R_0101", "K", "Ricci2", "owned_W", "R4", "R_h", "grad2")
    table = {}
    def shape_gap(left, right, channel):
        a = left["_curvature_profiles"][channel]
        b = right["_curvature_profiles"][channel]
        if b.size % a.size:
            raise ValueError("refinement profiles do not share the saved periodic nodes")
        sampled = b[::b.size // a.size]
        delta = sampled - a
        return {
            "protocol": "fine profile restricted to the shared coarse periodic nodes",
            "max_relative": float(np.max(np.abs(delta)) / max(np.max(np.abs(sampled)), 1e-300)),
            "l2_relative": float(np.linalg.norm(delta) / max(np.linalg.norm(sampled), 1e-300)),
        }
    for time_value in STATION_TIMES:
        left = _at_time(coarser, time_value)
        right = _at_time(finer, time_value)
        table[str(time_value)] = {
            channel: {
                "max_abs": _pair_gap(left, right, channel),
                "worldline_x2": _relative(right["worldline_x2"][channel], left["worldline_x2"][channel]),
                "peak_R0202": _relative(right["peak_R0202"][channel], left["peak_R0202"][channel]),
                "shape": shape_gap(left, right, channel),
            }
            for channel in channels
        }
    return table


def compare_cases(cases):
    """Spatial, timestep, and clock comparisons. Peaks are each case's own peak."""
    by_id = {case["case_id"]: case for case in cases}
    primary = by_id[PRIMARY_CASE]
    coarse = by_id[COARSE_CASE]
    stepped = by_id[TIMESTEP_CASE]
    frozen = by_id[FROZEN_CASE]
    channels = ("R_0202", "R_0101", "K", "Ricci2", "owned_W", "R4", "R_h", "grad2")
    space = resolution_table(coarse, primary)
    timestep = {}
    for time_value in STATION_TIMES:
        timestep[str(time_value)] = {
            channel: _pair_gap(_at_time(primary, time_value), _at_time(stepped, time_value), channel)
            for channel in channels
        }
    frozen_rows = frozen["stations"]
    frozen_tide = [row["summaries"]["R_0202"]["max_abs"] for row in frozen_rows]
    frozen_field = [row["field_rate_max"] for row in frozen_rows]
    return {
        "primary_case": PRIMARY_CASE,
        "space_nf256_versus_nf128_dt0.0005": space,
        "timestep_nf256_dt0.001_versus_dt0.0005": timestep,
        "growth": _growth_series(primary),
        "coarse_growth": _growth_series(coarse),
        "frozen": {
            "case_id": FROZEN_CASE,
            "Q_t_max": [row["Q_t_max"] for row in frozen_rows],
            "r_t_max": [row["r_t_max"] for row in frozen_rows],
            "Q_tt_max": [row["Q_tt_max"] for row in frozen_rows],
            "r_tt_max": [row["r_tt_max"] for row in frozen_rows],
            "R_0202_max_abs": frozen_tide,
            "field_rate_max": frozen_field,
            "geometry_jets_zero": all(
                row["Q_t_max"] == 0.0 and row["r_t_max"] == 0.0
                and row["Q_tt_max"] == 0.0 and row["r_tt_max"] == 0.0
                for row in frozen_rows
            ),
            "tide_constant": max(frozen_tide) - min(frozen_tide) == 0.0,
            "field_still_flows": all(value > 0.0 for value in frozen_field) and max(frozen_field) != min(frozen_field),
        },
    }


def measure_directory(directory, case_ids):
    directory = Path(directory)
    if not directory.is_dir():
        return {"present": False, "directory": str(directory), "cases": []}
    before = _hash_tree(directory)
    with threadpool_limits(limits=1):
        cases = [measure_case(directory, case_id) for case_id in case_ids]
    after = _hash_tree(directory)
    return {
        "present": True,
        "directory": str(directory.resolve()),
        "cases": cases,
        "input_hashes_before": before,
        "input_hashes_after": after,
        "input_hashes_unchanged": before == after,
    }


def frozen_binding(directory, campaign):
    """Authenticate the saved run envelope and its physics dependencies.

    The runner is a current reader; its current bytes are bound separately.
    Physics dependencies must match the immutable producing envelope.
    """
    path = Path(directory) / "run-binding-envelope.json"
    binding = json.loads(path.read_text(encoding="utf-8"))
    commit = binding.get("producing_commit", "")
    if not commit.startswith(FROZEN_PRODUCING_COMMITS[campaign]):
        raise ValueError("unexpected producing commit for " + campaign)
    imported = binding.get("import_closure_sha256") or binding.get("imports") or {}
    physics = {}
    for name in _PHYSICS_NAMES:
        relative = "lab/src/recursive_horizons/" + name
        expected = imported.get(relative)
        target = episode.REPO / relative
        if expected is None or sha256_file(target) != expected:
            raise ValueError("frozen physics binding changed: " + relative)
        physics[str(target.resolve())] = expected
    return {
        "envelope": str(path.resolve()), "envelope_sha256": sha256_file(path),
        "producing_commit": commit, "physics_hashes": physics,
        "physics_matches_saved_producer": True,
    }


def exact_suite():
    """Local jets of three geometries. No periodic seam is differenced."""
    results = {}

    # Spherical Minkowski: r = x, Q = 1/x, N = 1, at x = 2.
    x = 2.0
    Q = 1.0 / x
    Qx = -1.0 / x ** 2
    Qxx = 2.0 / x ** 3
    minkowski = curvature_from_local_jets(Q, x, 0.0, Qx, 0.0, 0.0, Qxx, 0.0, 1.0, 0.0, 0.0, 0.0)
    results["minkowski"] = _exact_pack(minkowski, expected_flat=True)

    # Milne: r = e^t sinh x, Q = 1/sinh x, at t = 0.3, x = 1.2.
    t, x = 0.3, 1.2
    sh, ch = math.sinh(x), math.cosh(x)
    et = math.exp(t)
    r = et * sh
    Q = 1.0 / sh
    Qx = -ch / sh ** 2
    Qxx = Q * (ch ** 2 / sh ** 2 + 1.0 / sh ** 2)
    milne = curvature_from_local_jets(
        Q, r, 0.0, Qx, 0.0, 0.0, Qxx, r, et * ch, r, et * ch, r,
    )
    results["milne"] = _exact_pack(milne, expected_flat=True)

    # Bertotti–Robinson product: r constant, Q = Q0 sech(Q0 t), R_h = 2.
    Q0, t, radius = 0.4, 0.35, 3.0
    th = math.tanh(Q0 * t)
    Q = Q0 / math.cosh(Q0 * t)
    Qt = -Q0 * Q * th
    Qtt = Q0 ** 2 * Q * (2.0 * th ** 2 - 1.0)
    br = curvature_from_local_jets(Q, radius, Qt, 0.0, Qtt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    expected = {
        "R_h": 2.0,
        "R4": 0.0,
        "owned_W": 0.0,
        "Ricci2": 4.0 / radius ** 4,
        "K": 8.0 / radius ** 4,
        "R_0101": -1.0 / radius ** 2,
        "R_0202": 0.0,
        "assembled_W": 0.0,
    }
    results["bertotti_robinson"] = _exact_pack(br, expected=expected)

    # R x S^2, static, r = 2, N = q = 1. Tides vanish; scalar invariants do not.
    product = curvature_from_local_jets(0.5, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    results["product_time_sphere"] = _exact_pack(product, expected={
        "R_h": 0.0,
        "R4": -0.5,
        "R_0101": 0.0,
        "R_0202": 0.0,
        "Ricci2": 0.125,
        "K": 0.25,
        "owned_W": 4.0 / (3.0 * 16.0),
    })
    return results


def _scalar_of(package, name):
    tides, base = package["tides"], package["base"]
    if name in tides:
        return float(np.asarray(tides[name]).reshape(-1)[0])
    return float(np.asarray(base[name]).reshape(-1)[0])


def _exact_pack(package, *, expected_flat=False, expected=None):
    names = ("R_h", "R4", "owned_W", "R_0101", "R_0202", "Ricci2", "K", "grad2", "assembled_W")
    values = {name: _scalar_of(package, name) for name in names}
    values["algebra_gap"] = values["assembled_W"] - values["owned_W"]
    values["static_section_R"] = float(np.asarray(package["static_sections_at_frozen_time"]["R"]).reshape(-1)[0])
    if expected_flat:
        # R_h = 2 and grad2 = -1 are the flat spherical chart. The 4D
        # curvature tensors, the Weyl scalar and the observer tides vanish.
        curved = ("R4", "owned_W", "R_0101", "R_0202", "Ricci2", "K", "assembled_W")
        values["flat_curvature_max_abs"] = max(abs(values[name]) for name in curved)
    if expected is not None:
        values["expected_gap"] = {
            name: values[name] - target for name, target in expected.items()
        }
    return values


def static_conformal_section_gap():
    """A static conformal jet: dynamic formulas versus the owned sections."""
    r, Q = 2.2, 1.5
    rx, rxx = 0.3, -0.4
    Qx, Qxx = 0.1, -0.2
    package = curvature_from_local_jets(Q, r, 0.0, Qx, 0.0, 0.0, Qxx, 0.0, rx, 0.0, 0.0, rxx)
    N = r * Q
    Nx = rx * Q + r * Qx
    Nxx = rxx * Q + 2.0 * rx * Qx + r * Qxx
    sections = owner_static_sections(N, N, r, Nx, Nxx, Nx, rx, rxx)
    return {
        "R": abs(_scalar_of(package, "R4") - float(sections["R"])),
        "Ricci2": abs(_scalar_of(package, "Ricci2") - float(sections["Ricci2"])),
        "K": abs(_scalar_of(package, "K") - float(sections["K"])),
    }


def directional_radius_chain_gap():
    """Chain rule of r_t = pi/(2 Z) against a one-step probe, not an evolution."""
    area = 0.37
    radius = 1.8
    p_r = 0.2
    p_chi = -0.15
    radius_rate = 0.05
    p_r_rate = -0.07
    p_chi_rate = 0.03
    f_chi = float(action.partial_F(radius, area, -1.0)[1])
    Z = float(action.feedback_Z(area))
    force = -8.0 * math.pi * area * radius
    pi_value = p_r - force * p_chi / f_chi
    r_dot = pi_value / (2.0 * Z)
    force_along = -8.0 * math.pi * area * radius_rate
    analytic = (p_r_rate - (force_along * p_chi + force * p_chi_rate) / f_chi) / (2.0 * Z)
    step = 1e-6

    def velocity(delta):
        r = radius + delta * radius_rate
        pr = p_r + delta * p_r_rate
        pc = p_chi + delta * p_chi_rate
        fr = -8.0 * math.pi * area * r
        return (pr - fr * pc / f_chi) / (2.0 * Z)

    probe = (velocity(step) - velocity(-step)) / (2.0 * step)
    return {
        "r_dot": r_dot,
        "analytic": analytic,
        "probe": probe,
        "gap": abs(analytic - probe),
    }


def measurement_report(episode_dir=None, confirmation_dir=None, extra_paths=(), *, include_nf512=True):
    """Read saved stations and bind every file hash that was opened.

    The library does not write a record and does not step the state. The
    saved nf=512 confirmation can be omitted explicitly by the caller.
    """
    import time
    started = time.process_time()
    episode_dir = EPISODE_DIR if episode_dir is None else Path(episode_dir)
    confirmation_dir = CONFIRMATION_DIR if confirmation_dir is None else Path(confirmation_dir)
    producers_before = producer_hashes(extra_paths)
    bindings = {"episode": frozen_binding(episode_dir, "episode")}
    first = measure_directory(episode_dir, episode.baseline_case_ids())
    confirmation_ids = tuple(spec["case_id"] for spec in episode.confirmation_specs())
    if include_nf512:
        bindings["confirmation"] = frozen_binding(confirmation_dir, "confirmation")
        confirmation = measure_directory(confirmation_dir, confirmation_ids)
        confirmation["status"] = "measured"
    else:
        confirmation = {
            "present": False,
            "status": "not_selected",
            "directory": str(Path(confirmation_dir)),
            "expected_cases": list(confirmation_ids),
            "cases": [],
            "input_hashes_unchanged": True,
        }
    producers_after = producer_hashes(extra_paths)
    comparisons = compare_cases(first["cases"]) if first["present"] else {}
    finer = None
    if confirmation.get("present"):
        by_finer = {case["case_id"]: case for case in confirmation.get("cases", [])}
        candidate = by_finer.get("nf512_coupled_dt0.0005")
        if candidate is not None and len(candidate["stations"]) == 3:
            primary = {case["case_id"]: case for case in first["cases"]}[PRIMARY_CASE]
            finer = resolution_table(primary, candidate)
    if comparisons:
        comparisons["space_nf512_versus_nf256_dt0.0005"] = finer
    all_cases = first.get("cases", []) + confirmation.get("cases", [])
    for case in all_cases:
        for row in case["stations"]:
            row.pop("_curvature_profiles", None)
    finite = bool(all_cases) and all(
        row["finite"] and row["chart_positive"]
        for case in all_cases
        for row in case["stations"]
    )
    cpu = time.process_time() - started
    return {
        "schema": SCHEMA,
        "frozen_producer_bindings": bindings,
        "observer": {
            "frame": "e0=N^-1 partial_t; e1=N^-1 partial_x; e2=r^-1 partial_theta; e3=(r sin(theta))^-1 partial_phi",
            "normal_covariant_tidal_tensor": "diag(-A, B00, B00); components R_0i0j in signature +---",
            "normal_worldline_acceleration": "a^hat1=N_x/N^2; fixed x normals are generally accelerated",
            "freefall_interpretation": "The tidal tensor gives instantaneous geodesic deviation for freely falling particles momentarily comoving with the normal. Fixed-x worldlines require acceleration; their relative acceleration also includes the acceleration field.",
        },
        "riemann_convention": RIEMANN_CONVENTION,
        "formulas": {
            "R_h": FORMULA_RH,
            "R4": FORMULA_R4,
            "R_0101": FORMULA_R0101,
            "R_0202": FORMULA_R0202,
            "owned_W": FORMULA_OWNED_W,
            "weyl_algebra": FORMULA_WEYL_ALGEBRA,
            "acceleration": FORMULA_ACCELERATION,
            "gamma_scalar": "unit sphere scalar curvature 2",
            "candidate_checked_against_owner": True,
        },
        "chi_substituted_for_curvature": CHI_SUBSTITUTED,
        "new_ode": NEW_ODE,
        "geometry_evolved": False,
        "production_record_written": PRODUCTION_RECORD_WRITTEN_BY_LIBRARY,
        "singularity_declared": SINGULARITY_DECLARED,
        "theory_declared_dead": THEORY_DECLARED_DEAD,
        "finite_window": {
            "stations": list(STATION_TIMES),
            "chart_positive": finite,
            "values_finite": finite,
            "statement": (
                "The saved window reaches coordinate time 3 with a positive chart "
                "and finite tidal components. This finite window is the measured domain."
            ),
        },
        "exact": exact_suite(),
        "static_section_gap": static_conformal_section_gap(),
        "radius_chain": directional_radius_chain_gap(),
        "episode_v1": {
            "directory": first.get("directory"),
            "input_hashes_unchanged": first.get("input_hashes_unchanged"),
            "declared_source_pins": {
                case["case_id"]: case["declared_source_pins"] for case in first.get("cases", [])
            },
            "cases": [
                {
                    "case_id": case["case_id"],
                    "stations": [
                        {key: row[key] for key in row if key != "declared_source_pins"}
                        for row in case["stations"]
                    ],
                }
                for case in first.get("cases", [])
            ],
        },
        "input_hashes": {
            "before": first.get("input_hashes_before", {}),
            "after": first.get("input_hashes_after", {}),
            "unchanged": first.get("input_hashes_unchanged"),
        },
        "confirmation_nf512": confirmation,
        "comparisons": comparisons,
        "producer_hashes": producers_after,
        "producer_hashes_unchanged": producers_before == producers_after,
        "saved_state_hashes_unchanged": (
            bool(first.get("input_hashes_unchanged"))
            and bool(confirmation.get("input_hashes_unchanged", True))
            and producers_before == producers_after
        ),
        "cpu_seconds": cpu,
    }


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value
