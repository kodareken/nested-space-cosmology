"""The kappa = 0 lowest-Landau charged column on the dyonic slice.

The [charged source](nsc-rn-charged-source.md) note declares the lowest
multiplet: j = (|q_mag| - 1)/2, multiplicity |q_mag|, and Dirac angular
eigenvalue kappa = 0. The neutral owner `nsc_rn_source` owns the kappa = 1
shell only and refuses kappa = 0, so this module is a separate owner for the
kappa = 0 sector. It does not change the neutral producer.

Domain. One fixed dyonic slice of the fixed-charge branch, on the vacuum
chart, with a declared magnetic flux and a declared black-hole Gauss charge.
No time step is taken and no metric is evolved.

The mode. On the characteristic basis the kappa = 0 system decouples, because
the angular coupling alpha = kappa N / r vanishes:

    a+ d_plus  = ( i E_eff - a+'/2 ) plus
    a- d_minus = ( i E_eff - a-'/2 ) minus

The future-horizon regular root forces plus(r+) = 0, and the plus equation is
linear and homogeneous, so **plus vanishes identically** and the column is
carried entirely by the throughgoing minus characteristic. That is a property
of kappa = 0, not an approximation.

The electric window enters as a local shift of the frequency,

    E_eff(r) = omega - s Phi(r),      Phi(r) = Q_E / (16 pi C_F r),

which is the note's closed form on the vacuum chart N = 1. A coupled lapse is
refused rather than silently given a manufactured tail: the closed form is
exact only at unit lapse, and the note derives it there.
"""
from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np
from scipy.integrate import solve_ivp

from . import nsc_rn_charged_window as window
from . import nsc_rn_pg as pg
from . import nsc_rn_reference as reference
from . import nsc_rn_source as source

API_VERSION = "nsc-rn-charged-lll-v1"
LLL_KAPPA = 0
#: Same relative resolution as the neutral production route: 876 nodes across
#: the declared parent cap 32 r_m,tot of the charged solution itself.
DEFAULT_POINTS = 876
CAP_MULTIPLE = 32.0
HORIZON_TOL_FRACTION = 0.05


def dyonic_reference(*, A: float, C_F: float, q_mag: float, Q_E: float,
                     mass_over_total_radius: float) -> Dict[str, Any]:
    """The source-free dyonic slice, with P_tot^2 in the RN function.

    ``RNReference`` computes its geometry from ``magnetic_r2``. For a dyonic
    slice that slot holds the total charge square, so the reference is built
    with P_tot^2 and the magnetic and electric pieces are recorded beside it.
    The horizons are checked against the window algebra rather than assumed.
    """
    squares = window.dyonic_charge_squares(A, C_F, q_mag, Q_E)
    horizons = window.dyonic_horizons(mass_over_total_radius, squares["total_charge_radius"])
    model = reference.RNReference(
        A=float(A),
        C_F=float(C_F),
        flux=float(q_mag),
        mass=float(horizons["mass"]),
        magnetic_r2=float(squares["P_total_squared"]),
        G_N=reference.newton_constant(A),
        V4=reference.AUDITED_V4,
    )
    got = model.horizons
    if got["r_plus"] is None or got["r_minus"] is None:
        raise RuntimeError("the dyonic slice has no horizon pair")
    if abs(got["r_plus"] - horizons["r_plus"]) > 1e-9 * max(1.0, horizons["r_plus"]):
        raise RuntimeError("reference horizons disagree with the window algebra")
    if abs(got["r_minus"] - horizons["r_minus"]) > 1e-9 * max(1.0, horizons["r_minus"]):
        raise RuntimeError("inner reference horizon disagrees with the window algebra")
    return {
        "model": model,
        "squares": squares,
        "horizons": horizons,
        "Phi_H": window.horizon_potential(Q_E, C_F, horizons["r_plus"]),
        "r_m_tot": float(squares["total_charge_radius"]),
        "extremality_ratio": float(mass_over_total_radius),
    }


def electric_potential(radius, *, Q_E: float, C_F: float) -> np.ndarray:
    """Phi(r) = Q_E/(16 pi C_F r), the note's unit-lapse closed form."""
    if C_F <= 0.0:
        raise ValueError("C_F must be positive")
    r = np.asarray(radius, dtype=float)
    if np.any(r <= 0.0) or not np.isfinite(r).all():
        raise ValueError("positive finite radii required")
    return float(Q_E) / (16.0 * math.pi * float(C_F) * r)


def charged_chart(model, radius) -> Dict[str, np.ndarray]:
    """Vacuum chart with the unit-lapse precondition made explicit."""
    chart = source.vacuum_chart(model, np.asarray(radius, dtype=float))
    if not np.allclose(np.asarray(chart["lapse"], dtype=float), 1.0, rtol=0.0, atol=1e-12):
        raise RuntimeError(
            "the charged window closed form is derived at unit lapse; a coupled "
            "lapse needs its own potential integral and is not assumed here"
        )
    return chart


def charged_grid(*, model, r_m_tot: float, points: int = DEFAULT_POINTS,
                 cap_multiple: float = CAP_MULTIPLE) -> Dict[str, Any]:
    """Uniform SBP grid from one cell above r_minus to the declared cap."""
    points = int(points)
    if points < 8:
        raise ValueError("at least eight nodes are required")
    r_plus = float(model.horizons["r_plus"])
    r_minus = float(model.horizons["r_minus"])
    outer = float(cap_multiple) * float(r_m_tot)
    if not r_minus < r_plus < outer:
        raise ValueError("the cap must lie outside the outer horizon")
    spacing = (outer - r_minus) / float(points - 1)
    inner = r_minus + HORIZON_TOL_FRACTION * spacing
    radius = inner + spacing * np.arange(points, dtype=float)
    radius[-1] = outer
    weights, derivative, q_matrix = pg.sbp_operator(points, spacing, layout="csr")
    if abs(float(radius[0]) - inner) > 1e-12 or abs(float(radius[-1]) - outer) > 1e-12:
        raise RuntimeError("grid ends drifted from the declared inner edge and cap")
    return {
        "radius": radius,
        "weights": weights,
        "derivative": derivative,
        "q_matrix": q_matrix,
        "spacing": float(spacing),
        "points": points,
        "inner_edge": float(inner),
        "outer": float(outer),
        "cap_multiple": float(cap_multiple),
        "r_plus": r_plus,
        "r_minus": r_minus,
    }


def _chart_at(model, radius: float) -> Dict[str, float]:
    chart = charged_chart(model, np.asarray([float(radius)], dtype=float))
    return {name: (float(value[0]) if isinstance(value, np.ndarray) else value)
            for name, value in chart.items() if name != "reference"}


def _minus_rhs(model, radius: float, minus: complex, *, omega: float, sign: int,
               Q_E: float, C_F: float) -> complex:
    point = _chart_at(model, radius)
    speed = float(point["speed_minus"])
    slope = float(point["speed_minus_derivative"])
    if abs(speed) < 1e-14:
        raise RuntimeError("throughgoing characteristic vanished on this node")
    phi = float(electric_potential(radius, Q_E=Q_E, C_F=C_F))
    effective = float(omega) - float(sign) * phi
    return (1j * effective * minus - 0.5 * slope * minus) / speed


def horizon_jet(*, model, omega: float, sign: int, Q_E: float, C_F: float) -> Dict[str, Any]:
    """Future-horizon regular initial data for the kappa = 0 column.

    At r+ the outgoing speed vanishes. The regular root of the plus row is
    plus(r+) = 0, so the whole column sits on the minus characteristic. The
    minus value is a free normalization and its slope follows from the minus
    row at the horizon.
    """
    r_plus = float(model.horizons["r_plus"])
    point = _chart_at(model, r_plus)
    if abs(float(point["speed_plus"])) > 1e-8:
        raise RuntimeError("the outer root is not a zero of the outgoing speed")
    speed_minus = float(point["speed_minus"])
    if abs(speed_minus) < 1e-14:
        raise ValueError("throughgoing characteristic vanished at the horizon")
    phi_H = float(electric_potential(r_plus, Q_E=Q_E, C_F=C_F))
    effective = float(omega) - float(sign) * phi_H
    minus = 1.0 + 0.0j
    minus_derivative = (
        1j * effective * minus - 0.5 * float(point["speed_minus_derivative"]) * minus
    ) / speed_minus
    return {
        "radius": r_plus,
        "plus": 0.0 + 0.0j,
        "minus": minus,
        "minus_derivative": minus_derivative,
        "effective_frequency_at_horizon": effective,
        "Phi_H": phi_H,
        "speed_minus": speed_minus,
        "window_member": bool(effective < 0.0),
    }


def _continuation_gap(r_plus: float, r_minus: float, radius: np.ndarray) -> float:
    gap = 1.0e-4 * r_plus
    caps = [gap, 0.25 * (r_plus - r_minus)]
    if float(radius[0]) < r_plus:
        caps.append(0.25 * (r_plus - float(radius[0])))
    if float(radius[-1]) > r_plus:
        caps.append(0.25 * (float(radius[-1]) - r_plus))
    gap = float(min(caps))
    if not np.isfinite(gap) or gap <= 0.0:
        raise ValueError("no room to step off the future horizon")
    return gap


def _taylor(jet: Dict[str, Any], radius: float) -> complex:
    delta = float(radius) - float(jet["radius"])
    return jet["minus"] + delta * jet["minus_derivative"]


def _integrate(model, jet, radius, indices, start, end, column, *, omega, sign,
               Q_E, C_F, rtol, atol) -> int:
    if indices.size == 0:
        return 0
    initial = _taylor(jet, start)

    def rhs(r, state):
        value = state[0] + 1j * state[1]
        rate = _minus_rhs(model, r, value, omega=omega, sign=sign, Q_E=Q_E, C_F=C_F)
        return [float(np.real(rate)), float(np.imag(rate))]

    solution = solve_ivp(
        rhs, (float(start), float(end)), [initial.real, initial.imag],
        method="DOP853", rtol=rtol, atol=atol, dense_output=True,
    )
    if not solution.success or solution.sol is None:
        raise RuntimeError(solution.message)
    samples = solution.sol(radius[indices])
    column[1, indices] = samples[0] + 1j * samples[1]
    return int(solution.nfev)


def charged_lll_mode(*, model, radius, omega: float, sign: int, Q_E: float, C_F: float,
                     rtol: float = 1.0e-7, atol: float = 1.0e-9) -> Dict[str, Any]:
    """The kappa = 0 charged column on the caller grid. No resampling."""
    grid = np.asarray(radius, dtype=float)
    if grid.ndim != 1 or grid.size < 8 or np.any(np.diff(grid) <= 0.0):
        raise ValueError("an increasing grid of at least eight nodes is required")
    r_plus = float(model.horizons["r_plus"])
    r_minus = float(model.horizons["r_minus"])
    if grid[0] <= r_minus:
        raise ValueError("the grid must stay strictly above the inner horizon")
    jet = horizon_jet(model=model, omega=omega, sign=sign, Q_E=Q_E, C_F=C_F)
    gap = _continuation_gap(r_plus, r_minus, grid)
    column = np.zeros((2, grid.size), dtype=complex)
    near = np.abs(grid - r_plus) <= gap
    for index in np.flatnonzero(near):
        column[1, index] = _taylor(jet, float(grid[index]))
    outward = np.flatnonzero(grid > r_plus + gap)
    inward = np.flatnonzero(grid < r_plus - gap)
    evaluations = _integrate(model, jet, grid, outward, r_plus + gap, float(grid[-1]),
                             column, omega=omega, sign=sign, Q_E=Q_E, C_F=C_F,
                             rtol=rtol, atol=atol)
    evaluations += _integrate(model, jet, grid, inward, r_plus - gap, float(grid[0]),
                              column, omega=omega, sign=sign, Q_E=Q_E, C_F=C_F,
                              rtol=rtol, atol=atol)
    if int(near.sum() + outward.size + inward.size) != grid.size:
        raise RuntimeError("radial nodes were dropped during horizon continuation")
    if not np.isfinite(column).all():
        raise RuntimeError("the charged column is not finite on this grid")
    return {
        "api_version": API_VERSION,
        "kappa": LLL_KAPPA,
        "mode": column,
        "omega": float(omega),
        "fermion_sign": int(sign),
        "Q_E": float(Q_E),
        "C_F": float(C_F),
        "horizon_jet": jet,
        "continuation_gap": float(gap),
        "evaluations": int(evaluations),
        "plus_identically_zero": bool(not np.any(column[0])),
        "nodes_dropped": 0,
        "node_count": int(grid.size),
        "resampling": False,
    }


def charged_mode_form(mode, *, chart, weights, derivative, radius, Q_E: float,
                      C_F: float, sign: int, omega: float) -> Dict[str, Any]:
    """Kappa = 0 Killing form: the angular row is identically zero.

    Mirrors ``nsc_rn_source.mode_quadratic_form`` at kappa = 0, where the
    angular bilinear carries the factor kappa and therefore drops. A test pins
    the kinetic rows against the neutral owner's own helper.
    """
    column = np.asarray(mode, dtype=complex)
    if column.shape != (2, np.asarray(radius).shape[0]):
        raise ValueError("the charged column must have shape (2, points)")
    lapse = np.asarray(chart["lapse"], dtype=float)
    shift = np.asarray(chart["shift"], dtype=float)
    speed_plus = (lapse - shift)
    speed_minus = -(lapse + shift)
    weight_row = np.asarray(weights, dtype=float)
    radial = np.asarray(radius, dtype=float)
    flux_plus = _hermitian_flux(column[0], column[0], derivative)
    flux_minus = _hermitian_flux(column[1], column[1], derivative)
    density = np.real(speed_plus * flux_plus + speed_minus * flux_minus)
    effective = np.asarray(omega, dtype=float) - float(sign) * electric_potential(
        radial, Q_E=Q_E, C_F=C_F
    )
    norm = float(np.sqrt(np.real(np.sum(weight_row * np.sum(np.abs(column) ** 2, axis=0)))))
    return {
        "density": density,
        "energy": float(np.sum(weight_row * density)),
        "probability_current": speed_plus * np.abs(column[0]) ** 2 + speed_minus * np.abs(column[1]) ** 2,
        "weighted_norm": norm,
        "effective_frequency": effective,
        "effective_frequency_min": float(np.min(effective)),
        "window_entered": bool(np.any(effective < 0.0)),
        "angular_row_dropped_by_kappa": True,
    }


def _hermitian_flux(left, right, derivative) -> np.ndarray:
    """The neutral owner's split flux, mirrored at this owner's boundary."""
    a = np.asarray(left, dtype=complex)
    b = np.asarray(right, dtype=complex)
    return (np.conj(a) * (derivative @ b) - np.conj(derivative @ a) * b) / (2j)


def normalized_column(mode, weights) -> Dict[str, Any]:
    """Weighted-Gram normalization and the centre phase alignment."""
    column = np.asarray(mode, dtype=complex)
    weight_row = np.asarray(weights, dtype=float)
    norm = float(np.sqrt(np.real(np.sum(weight_row * np.sum(np.abs(column) ** 2, axis=0)))))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("the charged column has vanishing weighted norm")
    scaled = column / norm
    through = scaled[1]
    index = int(np.argmax(np.abs(through)))
    value = through[index]
    if abs(value) > 0.0:
        scaled = scaled * np.conj(value / abs(value))
    return {"column": scaled, "raw_norm": norm, "phase_index": index}


def infall_winding(*, model, omega: float, sign: int, Q_E: float, C_F: float,
                   r_outer: float, points: int = 1200) -> Dict[str, Any]:
    """Phase accumulated by the infalling kappa = 0 wavefront, quadratic in
    omega.

    The throughgoing minus characteristic obeys dr/dt = -(1 + sqrt(f)) on
    the unit-lapse chart, so along it the phase is

        W(omega) = int [omega - s Phi(r)] dt
                 = omega T_total - s Q_E/(16 pi C_F) I_1

    with dt/dr = -1/(1 + sqrt(f)). Both quadratures are elementary in the
    integration variable r, so W is linear in omega and the exact zero is

        omega_star = s (Q_E/(16 pi C_F)) I_1 / T_total = s Phi_mean.

    This is the phase-averaged window potential: below omega_star the whole
    infall winds opposite to the exterior Killing phase.
    """
    if r_outer <= float(model.horizons["r_plus"]):
        raise ValueError("the outer radius must lie outside the horizon")
    radius = np.linspace(float(model.horizons["r_plus"]), float(r_outer), int(points))
    with np.errstate(invalid="ignore"):
        weight = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, 1.0 - 2.0 * float(model.mass) / radius
                                                  + float(model.magnetic_r2) / radius ** 2)))
    # Numerically stable replacement where sqrt(f) vanishes at r_+:
    weight = np.asarray(weight, dtype=float)
    f_guess = 1.0 - 2.0 * float(model.mass) / radius + float(model.magnetic_r2) / radius ** 2
    if np.any(f_guess < -1e-12):
        raise ValueError("the quadrature radius range left the horizon chart")
    dt_dr = -weight
    phi = float(Q_E) / (16.0 * math.pi * float(C_F) * radius)
    t_total = float(-np.trapezoid(dt_dr, radius))
    i_one = float(-np.trapezoid(dt_dr * phi, radius))
    phi_mean = i_one / t_total
    omega_star = float(sign) * phi_mean
    effective = float(omega) - float(sign) * phi
    in_window = effective < 0.0
    negative_part = float(-np.trapezoid((dt_dr * effective)[in_window], radius[in_window])) \
        if np.any(in_window) else 0.0
    phase_total = float(omega) * t_total - float(sign) * i_one
    return {
        "T_total": t_total,
        "I_one": i_one,
        "Phi_mean": phi_mean,
        "omega_star": omega_star,
        "sign": int(sign),
        "total_phase_at_omega": phase_total,
        "negative_window_phase": negative_part,
        "negative_window_fraction": float(negative_part / phase_total) if phase_total else None,
        "net_winding": "negative" if phase_total < 0.0 else "positive",
        "horizon_potential": float(electric_potential(float(model.horizons["r_plus"]), Q_E=Q_E, C_F=C_F)),
        "r_outer": float(r_outer),
        "points": int(points),
    }


def charged_twisted_column(*, model, radius, omega: float, sign: int, Q_E: float,
                           C_F: float, kappa: float, rtol: float = 1.0e-7,
                           atol: float = 1.0e-9) -> Dict[str, Any]:
    """The coupled nonzero-kappa charged column on the caller grid.

    The characteristic system with angular coupling alpha = kappa N/r, N = 1:

        a_+ d(chi_+)/dr = ( i E_eff - a_+'/2 ) chi_+ - alpha chi_-
        a_- d(chi_-)/dr = ( i E_eff - a_-'/2 ) chi_- + alpha chi_+

    At r_+ the outgoing speed vanishes. The plus row then constrains the
    regular root at the horizon:

        0 = ( i E_eff - a_+'/2 ) chi_+(r_+) - alpha chi_-(r_+)

    For kappa = 0 that fixes chi_+ = 0 (the previous owner). For nonzero
    kappa both rows are populated; the plus characteristic is the channel
    through which a horizon-mode can radiate outward, which is the first
    reflection amplitude available to this branch. The horizon derivative
    of chi_+ is obtained by differentiating the plus row (L'Hopital).
    """
    grid = np.asarray(radius, dtype=float)
    if grid.ndim != 1 or grid.size < 8 or np.any(np.diff(grid) <= 0.0):
        raise ValueError("an increasing grid of at least eight nodes is required")
    r_plus = float(model.horizons["r_plus"])
    r_minus = float(model.horizons["r_minus"])
    if grid[0] <= r_minus:
        raise ValueError("the grid must stay strictly above the inner horizon")
    kappa = float(kappa)
    if not np.isfinite(kappa) or abs(kappa) <= 0.0:
        raise ValueError("the twisted column needs a strictly nonzero angular scale")

    point_h = _chart_at(model, r_plus)
    a_plus_prime = float(point_h["speed_plus_derivative"])
    a_minus = float(point_h["speed_minus"])
    a_minus_prime = float(point_h["speed_minus_derivative"])
    phi_h = float(electric_potential(r_plus, Q_E=Q_E, C_F=C_F))
    e_h = float(omega) - float(sign) * phi_h
    alpha_h = kappa / r_plus

    # Taylor data at the horizon: chi_- free (set to 1), chi_+ fixed by the
    # plus row, both derivatives by L'Hopital / the minus row.
    plus_h = alpha_h / (1j * e_h - 0.5 * a_plus_prime)
    minus_h = 1.0 + 0.0j
    minus_p = ((1j * e_h - 0.5 * a_minus_prime) * minus_h + alpha_h * plus_h) / a_minus
    # Derivatives of the coefficient functions at the horizon.
    e_prime = float(sign) * float(Q_E) / (16.0 * math.pi * float(C_F) * r_plus ** 2)
    g = lambda r: 1.0 - 2.0 * float(model.mass) / r + float(model.magnetic_r2) / r ** 2
    gp = float(model.mass) / r_plus ** 2 - float(model.magnetic_r2) / r_plus ** 3  # g'/2
    gpp = -float(model.mass) / r_plus ** 3 + float(model.magnetic_r2) / r_plus ** 4  # g''/4
    # a_+' = -beta', beta' = gp (with g=1 at the horizon)
    a_pp = -gpp
    alpha_prime = -kappa / r_plus ** 2
    plus_p = (
        (1j * e_prime - 0.5 * a_pp) * plus_h - alpha_prime * minus_h - alpha_h * minus_p
    ) / (a_plus_prime - (1j * e_h - 0.5 * a_plus_prime))

    pivot = {
        "plus": plus_h, "minus": minus_h,
        "plus_derivative": plus_p, "minus_derivative": minus_p,
        "effective_frequency_at_horizon": e_h,
        "alpha": alpha_h, "a_plus_prime": a_plus_prime, "a_minus": a_minus,
    }

    def rhs(r, state):
        pt = _chart_at(model, float(r))
        a_plus = float(pt["speed_plus"])
        a_minus = float(pt["speed_minus"])
        a_plus_p = float(pt["speed_plus_derivative"])
        a_minus_p = float(pt["speed_minus_derivative"])
        phi = float(electric_potential(float(r), Q_E=Q_E, C_F=C_F))
        e = float(omega) - float(sign) * phi
        alpha = float(kappa) / float(r)
        plus = state[0] + 1j * state[1]
        minus = state[2] + 1j * state[3]
        plus_rate = ((1j * e - 0.5 * a_plus_p) * plus - alpha * minus) / a_plus
        minus_rate = ((1j * e - 0.5 * a_minus_p) * minus + alpha * plus) / a_minus
        return [plus_rate.real, plus_rate.imag, minus_rate.real, minus_rate.imag]

    gap = _continuation_gap(r_plus, r_minus, grid)
    column = np.zeros((2, grid.size), dtype=complex)

    def taylor(r):
        delta = float(r) - r_plus
        return (
            plus_h + delta * plus_p,
            minus_h + delta * minus_p,
        )

    near = np.abs(grid - r_plus) <= gap
    for index in np.flatnonzero(near):
        plus_v, minus_v = taylor(float(grid[index]))
        column[0, index] = plus_v
        column[1, index] = minus_v
    outward = np.flatnonzero(grid > r_plus + gap)
    initial = [plus_h.real, plus_h.imag, minus_h.real, minus_h.imag]
    evaluations = 0
    for indices, start, end in ((outward, r_plus + gap, float(grid[-1])),
                                (np.flatnonzero(grid < r_plus - gap), r_plus - gap, float(grid[0]))):
        if indices.size == 0:
            continue
        init_plus, init_minus = taylor(start)
        init = [init_plus.real, init_plus.imag, init_minus.real, init_minus.imag]
        solution = solve_ivp(
            rhs, (start, end), init, method="DOP853", rtol=rtol, atol=atol,
            dense_output=True,
        )
        if not solution.success or solution.sol is None:
            raise RuntimeError(solution.message)
        samples = solution.sol(grid[indices])
        column[0, indices] = samples[0] + 1j * samples[1]
        column[1, indices] = samples[2] + 1j * samples[3]
        evaluations += int(solution.nfev)
    if not np.isfinite(column).all():
        raise RuntimeError("the twisted column is not finite on this grid")
    return {
        "api_version": API_VERSION,
        "kappa": kappa,
        "mode": column,
        "omega": float(omega),
        "fermion_sign": int(sign),
        "Q_E": float(Q_E),
        "C_F": float(C_F),
        "horizon_pivot": pivot,
        "continuation_gap": float(gap),
        "evaluations": int(evaluations),
        "nodes": int(grid.size),
        "plus_populated": bool(np.any(np.abs(column[0]) > 0.0)),
        "nodes_dropped": 0,
        "resampling": False,
    }


def reflection_ratio(column, chart, *, speed_plus=None, speed_minus=None) -> Dict[str, Any]:
    """Outgoing/ingoing characteristic flux ratio at the outermost node.

    The plus speed is positive outside the horizon, so |a_+| |chi_+|^2 is
    the outgoing flux carried by the reflected part and |a_-| |chi_-|^2
    the ingoing flux. For kappa = 0 the numerator vanishes identically.
    """
    array = np.asarray(column, dtype=complex)
    if array.shape[1] < 2:
        raise ValueError("need at least two radial nodes")
    sp = np.asarray(chart["speed_plus"], dtype=float)
    sm = np.asarray(chart["speed_minus"], dtype=float)
    outgoing = float(sp[-1]) * float(np.abs(array[0, -1]) ** 2)
    ingoing = float(-sm[-1]) * float(np.abs(array[1, -1]) ** 2)
    ratio = outgoing / ingoing if ingoing > 0.0 else None
    return {
        "outgoing_flux": outgoing,
        "ingoing_flux": ingoing,
        "reflection_ratio": ratio,
        "radius_outer": float(array.shape[1]),
        "note": "ratios taken at the outermost node; kappa=0 makes the numerator vanish",
    }
