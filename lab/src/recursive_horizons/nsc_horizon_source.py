"""State-defined source of the magnetic Dirac zero-mode sector.

Known 2D conformal stress and Einstein spherical reduction are applied to
the declared unwrapped geometry. Units: hbar=c_light=kB=1. The central
charge is a mode count, not the speed of light. This is the free massless
LLL sector; magnetic flux, full UV matching and global state are not fixed.
"""
from math import pi


def state_constants(surface_gravity, state, central_charge=1, circle_constant=pi):
    scale = central_charge * surface_gravity ** 2 / (48 * circle_constant)
    if state == "Boulware":
        return 0 * scale, 0 * scale
    if state == "Unruh":
        return scale, 0 * scale
    if state == "Hartle-Hawking":
        return scale, scale
    raise ValueError("unknown stationary conformal-state control")


def conformal_stress(A, A_prime, A_second, t_u, t_v,
                     central_charge=1, circle_constant=pi):
    """Covariant null components for ds2=A du dv, signature +−."""
    vacuum = central_charge * (2 * A * A_second - A_prime ** 2) / (192 * circle_constant)
    mixed = central_charge * A * A_second / (96 * circle_constant)
    return vacuum + t_u, mixed, vacuum + t_v


def pg_source(A, beta, radius, T_uu, T_uv, T_vv, circle_constant=pi):
    """Project to the existing PG chart away from A=0.

    The 4D result is the spherical projection of this effective LLL action,
    not the complete renormalized 4D stress. Horizon limits are taken
    analytically in the application record.
    """
    area = 4 * circle_constant * radius ** 2
    du_rho, dv_rho = -(1 + beta) / A, (1 - beta) / A
    T_tt = T_uu + 2 * T_uv + T_vv
    T_tr = du_rho * T_uu + (du_rho + dv_rho) * T_uv + dv_rho * T_vv
    T_rr = du_rho ** 2 * T_uu + 2 * du_rho * dv_rho * T_uv + dv_rho ** 2 * T_vv
    current = -beta * T_tt - A * T_tr
    return {
        "T_tau_tau_2D": T_tt, "T_tau_rho_2D": T_tr, "T_rho_rho_2D": T_rr,
        "parent_Killing_power": current,
        "fourD_null_plus": 4 * T_vv / ((1 + beta) ** 2 * area),
        "fourD_null_minus": 4 * T_uu / ((1 - beta) ** 2 * area),
    }


def static_reduced_equations(A, A_prime, radius, radius_prime,
                             xi, potential_u, potential_du, state_kappa):
    """Static Einstein + potential + conformal-Dirac source equations.

    xi=G*c_eff/(12*pi); U(r)=8*pi*G*r²*rho_pot(r).
    U and dU/dr must come from matching the SAME functional. They are
    inputs to this interface, not fitted outputs or an inserted fluid.
    The state has t_u=t_v=c_eff*state_kappa²/(48*pi).
    Domain: A!=0, r>0, r²!=xi. Global horizon/state conditions remain open.
    """
    B = A_prime ** 2 - 4 * state_kappa ** 2
    A_second = (-radius * potential_du - 2 * A_prime * radius * radius_prime
                - xi * B / (2 * A)) / (radius ** 2 - xi)
    radius_second = -xi * (2 * A * A_second - B) / (4 * A ** 2 * radius)
    constraint = (1 - A * radius_prime ** 2 - A_prime * radius * radius_prime
                  - potential_u - xi * B / (4 * A))
    return A_second, radius_second, constraint


def reconstructed_neck(radius, A_prime, xi, potential_u, potential_du, state_kappa):
    """Required local neck data at r'=0; no potential/coupling is selected."""
    deficit = 1 - potential_u
    A = xi * (A_prime ** 2 - 4 * state_kappa ** 2) / (4 * deficit)
    A_second, radius_second, constraint = static_reduced_equations(
        A, A_prime, radius, 0, xi, potential_u, potential_du, state_kappa)
    return {"A": A, "A_second": A_second, "radius_second": radius_second,
            "constraint": constraint}
