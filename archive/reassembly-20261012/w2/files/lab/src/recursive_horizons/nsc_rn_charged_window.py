"""The dyonic electric window and its parent-cap budget.

One finite algebraic gate, taken before any mode solver or march is written
for the charged branch. The [charged source](nsc-rn-charged-source.md) note
declares the window and the cap; this module evaluates them on the audited
coefficients and answers a single yes/no question:

    does the electric window B_< contain a frequency that fits the parent cap?

If it does not, the charged campaign cannot be run as declared, and the next
step is a changed declaration rather than a mode solver.

Domain. The source-free Reissner-Nordstrom chart of the fixed-charge branch
with a declared magnetic flux q_mag and a declared black-hole Gauss charge
Q_E. Neither is fitted. No time step is taken and no field is propagated.

The algebra is the note's own, restated here so the driver cannot drift from
it:

    P_mag^2 = C_F q_mag^2 / (4 A)
    P_elec^2 = Q_E^2 / (256 pi^2 C_F A)
    P_tot^2 = P_mag^2 + P_elec^2
    r_m,tot = sqrt(P_tot^2)
    M       = eta r_m,tot,               eta > 1 subextremal
    r_+     = r_m,tot (eta + sqrt(eta^2 - 1))
    Phi_H   = Q_E / (16 pi C_F r_+)
    room    = 32 r_m,tot - r_+
    omega_min = 4 / room

A packet centred at omega carries the budget width sigma_omega = 1/omega and
must satisfy 4 sigma_omega <= room, so omega >= omega_min. The window is
B_< = {omega : 0 < omega < s Phi_H}. The gate is

    omega_min < s Phi_H,

with every quantity in one self-consistent unit system. Equality is not a
window: it is a single endpoint and is reported as closed.
"""
from __future__ import annotations

import math
from typing import Any, Dict

API_VERSION = "nsc-rn-charged-window-v1"

#: The note's parent cap: the outer boundary sits at this multiple of the
#: solution's own total-charge radius, not at 32 magnetic radii.
PARENT_CAP_MULTIPLE = 32.0
#: The declared budget width relation 4 sigma_omega <= room with
#: sigma_omega = 1/omega.
BUDGET_NUMERATOR = 4.0


def dyonic_charge_squares(A: float, C_F: float, q_mag: float, Q_E: float) -> Dict[str, float]:
    """Magnetic, electric and total charge squares of the declared slice."""
    if not all(math.isfinite(v) for v in (A, C_F, q_mag, Q_E)):
        raise ValueError("A, C_F, q_mag and Q_E must be finite")
    if A <= 0.0 or C_F <= 0.0:
        raise ValueError("A and C_F must be positive")
    p_mag2 = C_F * q_mag * q_mag / (4.0 * A)
    p_elec2 = Q_E * Q_E / (256.0 * math.pi ** 2 * C_F * A)
    return {
        "P_mag_squared": p_mag2,
        "P_electric_squared": p_elec2,
        "P_total_squared": p_mag2 + p_elec2,
        "magnetic_radius": math.sqrt(p_mag2),
        "total_charge_radius": math.sqrt(p_mag2 + p_elec2),
    }


def dyonic_horizons(mass_over_total_radius: float, total_charge_radius: float) -> Dict[str, float]:
    """Outer root and discriminant of f = 1 - 2M/r + P_tot^2/r^2.

    ``mass_over_total_radius`` is the extremality ratio M/sqrt(P_tot^2). The
    subextremal branch requires it strictly above one.
    """
    eta = float(mass_over_total_radius)
    r_m = float(total_charge_radius)
    if not math.isfinite(eta) or not math.isfinite(r_m) or r_m <= 0.0:
        raise ValueError("finite positive total-charge radius required")
    if eta <= 1.0:
        raise ValueError("the dyonic arm requires a subextremal ratio M/sqrt(P_tot^2) > 1")
    mass = eta * r_m
    root = math.sqrt(eta * eta - 1.0)
    return {
        "mass": mass,
        "extremality_ratio": eta,
        "r_plus": r_m * (eta + root),
        "r_minus": r_m * (eta - root),
        "discriminant_over_r_m_squared": eta * eta - 1.0,
    }


def horizon_potential(Q_E: float, C_F: float, r_plus: float) -> float:
    """Phi_H = Q_E/(16 pi C_F r_+), the note's derived horizon potential."""
    if not math.isfinite(Q_E) or not math.isfinite(r_plus) or r_plus <= 0.0:
        raise ValueError("finite positive r_plus required")
    if C_F <= 0.0:
        raise ValueError("C_F must be positive")
    return Q_E / (16.0 * math.pi * C_F * r_plus)


def cap_room(total_charge_radius: float, r_plus: float) -> float:
    """32 r_m,tot - r_+, the radial room the declared cap leaves a packet."""
    room = PARENT_CAP_MULTIPLE * float(total_charge_radius) - float(r_plus)
    if not math.isfinite(room):
        raise ValueError("the cap room is not finite")
    return room


def window_gate(*, A: float, C_F: float, q_mag: float, Q_E: float,
                mass_over_total_radius: float, sign: int) -> Dict[str, Any]:
    """Decide whether B_< holds a frequency the parent cap can carry.

    ``sign`` is the fermion charge sign s = +-1. Only s Phi_H > 0 opens the
    electric window; the other sign keeps every omega > 0 in B_>.
    """
    if sign not in (-1, 1):
        raise ValueError("the fermion charge sign is +1 or -1")
    squares = dyonic_charge_squares(A, C_F, q_mag, Q_E)
    r_m_tot = squares["total_charge_radius"]
    horizons = dyonic_horizons(mass_over_total_radius, r_m_tot)
    phi_H = horizon_potential(Q_E, C_F, horizons["r_plus"])
    room = cap_room(r_m_tot, horizons["r_plus"])
    shifted = sign * phi_H
    if room > 0.0:
        omega_min = BUDGET_NUMERATOR / room
    else:
        omega_min = math.inf
    window_upper = shifted
    admitted = bool(window_upper > 0.0 and omega_min < window_upper)
    if window_upper <= 0.0:
        reason = "the shifted potential is not positive for this fermion sign"
    elif room <= 0.0:
        reason = "the outer horizon leaves no room inside the declared cap"
    elif not admitted:
        reason = "every frequency in B_< is narrower than the declared budget width"
    else:
        reason = "B_< holds a frequency the declared cap can carry"
    return {
        "api_version": API_VERSION,
        "charge_squares": squares,
        "extremality_ratio": float(mass_over_total_radius),
        "mass": horizons["mass"],
        "r_plus": horizons["r_plus"],
        "r_minus": horizons["r_minus"],
        "Phi_H": phi_H,
        "fermion_sign": int(sign),
        "shifted_potential": shifted,
        "window_open": bool(window_upper > 0.0),
        "window_upper": window_upper,
        "cap_room": room,
        "omega_min": omega_min,
        "admitted": admitted,
        "reason": reason,
        "omega_central": (0.5 * (omega_min + window_upper)) if admitted else None,
        "band_width": (window_upper - omega_min) if admitted else None,
        "neutral_frequency": (1.0 / r_m_tot) if r_m_tot > 0.0 else None,
        "neutral_frequency_band": (
            "B_<" if 0.0 < (1.0 / r_m_tot) < window_upper else "B_>"
        ) if r_m_tot > 0.0 and window_upper > 0.0 else "no window",
    }


def pure_magnetic_control(*, A: float, C_F: float, q_mag: float,
                          mass_over_total_radius: float, sign: int) -> Dict[str, Any]:
    """The same slice with Q_E = 0. Phi_H vanishes, so no window opens."""
    return window_gate(
        A=A, C_F=C_F, q_mag=q_mag, Q_E=0.0,
        mass_over_total_radius=mass_over_total_radius, sign=sign,
    )
