"""Fermion pair production by shear in a Kasner collapse (the stretching interior of a black hole).

Model. Near the centre of a Schwarzschild hole the interior is Kasner,
    ds^2 = dtau^2 - (tau/tau0)^(-2/3) dx^2 - (tau/tau0)^(4/3) (dy^2 + dz^2)   (c = hbar = 1),
with tau the proper time left to the crunch (tau decreasing = collapse). One direction stretches,
two shrink. A massless Dirac field in a diagonal Bianchi-I metric obeys, after rescaling by the
volume factor, i d(chi)/dt = [ sum_i alpha_i k_i / a_i(tau) ] chi: exactly, with no spin-connection
remainder. The sphere is treated as locally flat (valid for modes with angular momentum >> 1).

Reduction. Each Weyl block is H = sigma . q(tau), q = (k_x rho^(1/3), k_perp rho^(-2/3)), rho = tau/tau0.
Rotating to the frame of q-hat gives an exact two-level problem. With u = rho / r, r = k_perp/k_x and
lambda = (k_x tau0) r^(4/3), the problem depends on lambda alone:
    dc/d ln u = i [ lambda u^(1/3) sqrt(1+u^2) sigma_z - (1/2) u/(1+u^2) sigma_y ] c .
F(lambda) = probability that a negative-energy state ends in the positive-energy state as u -> 0.
Dirac field = two Weyl blocks, so the pairs per mode are 2 F.
Density of states: d^3k = 4 pi lambda^2 r^(-3) dlambda dr / tau0^3 (both signs of k_x, 2D transverse k; Jacobian
checked numerically). Physical volume = comoving volume * rho. A mode counts as produced when u_f = rho_f / r < u_c.
Pair density:  n_phys * tau^3 = C = (1/(2 pi^2)) u_c^2 Lambda_2,  Lambda_2 = integral of lambda^2 F dlambda,
u_c ~ 1 marking where a mode counts as produced (scanned 0.5 .. 2).
(An earlier version of this file and its note had 1/(4 pi^2): a factor of 2 low. The Dirac field has 2F pairs per k,
confirmed by the 4x4 trace above, and the 2 F was not carried through. Corrected 7 Oct 2026.)
Energy density of the produced pairs (each pair carries 2|q|, |q| tau = lambda u^(1/3) sqrt(1+u^2)):
    eps * tau^4 = (2/pi^2) Lambda_3 J(u_c),  Lambda_3 = integral lambda^3 F dlambda (cut at lambda <= 5),
    J(u_c) = integral_0^u_c u^(4/3) sqrt(1+u^2) du.
"""
import json
import math

import numpy as np
from scipy.integrate import solve_ivp, quad

SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)


def _rhs(s, c, lam):
    u = math.exp(s)
    a = lam * u ** (1 / 3) * math.sqrt(1 + u * u)
    b = 0.5 * u / (1 + u * u)
    return 1j * (a * (SZ @ c) - b * (SY @ c))


def evolve(lam, u_i, u_f, rtol=1e-10, atol=1e-12):
    sol = solve_ivp(_rhs, [math.log(u_i), math.log(u_f)], np.array([0, 1], dtype=complex),
                    args=(lam,), method="DOP853", rtol=rtol, atol=atol)
    return sol.y[:, -1]


def F(lam):
    """Asymptotic positive-energy probability starting adiabatically from u_i -> infinity."""
    u_i = max(30.0, (3e3 / lam) ** (3 / 7))
    u_f = min(1e-3, (lam * 1e-3) ** 1.5)
    c = evolve(lam, u_i, u_f)
    return float(abs(c[0]) ** 2), abs(c[0]) ** 2 + abs(c[1]) ** 2 - 1.0


def direct_4x4(kx, r, rho_f, steps=400000):
    """Full 4x4 massless Dirac evolution in the original frame, from rho = 1 down to rho_f."""
    s0 = np.eye(2, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = SY
    sz = SZ
    Z = np.zeros((2, 2), dtype=complex)
    alpha = lambda s: np.block([[Z, s], [s, Z]])
    ax, ay, az = alpha(sx), alpha(sy), alpha(sz)
    kp = kx * r

    def H(tau):
        return ax * kx * tau ** (1 / 3) + ay * kp * tau ** (-2 / 3)      # phi = 0 (y-direction)

    def rhs(t, y):
        U = y.reshape(4, 4)
        tau = 1.0 - t
        return (-1j * H(tau) @ U).ravel()

    def proj(tau, sign):
        h = H(tau)
        w = np.linalg.eigvalsh(h)[-1]
        return 0.5 * (np.eye(4) + sign * h / w)

    sol = solve_ivp(rhs, [0, 1 - rho_f], np.eye(4, dtype=complex).ravel(), method="DOP853",
                    rtol=1e-10, atol=1e-12)
    U = sol.y[:, -1].reshape(4, 4)
    Pm0, Pp1 = proj(1.0, -1), proj(rho_f, +1)
    return float(np.real(np.trace(Pp1 @ U @ Pm0 @ U.conj().T)))


def two_level_same_start(lam, u_i, u_f):
    c = evolve(lam, u_i, u_f)
    return float(abs(c[0]) ** 2)


def main():
    out = {}
    # 1. validate the reduction: full 4x4 against 2 x two-level, same finite start
    kx, r, rho_f = 2.0, 0.3, 0.01
    lam = kx * r ** (4 / 3)
    direct = direct_4x4(kx, r, rho_f)
    u_i, u_f = 1.0 / r, rho_f / r
    # in the reduced problem the start is the adiabatic negative-energy state at u_i
    reduced = 2 * two_level_same_start(lam, u_i, u_f)
    out["reduction_check"] = {"kx_tau0": kx, "r": r, "rho_f": rho_f, "lambda": lam,
                              "direct_4x4_pairs": direct, "two_times_two_level": reduced,
                              "abs_diff": abs(direct - reduced)}
    # 2. F(lambda) and Lambda_2
    lams = np.geomspace(1e-3, 60.0, 90)
    Fs, norm_err = [], []
    for L in lams:
        f, e = F(L)
        Fs.append(f)
        norm_err.append(abs(e))
    Fs = np.array(Fs)
    out["F_small_lambda"] = float(Fs[0])
    out["max_norm_error"] = float(max(norm_err))
    out["lambda_grid"] = lams.tolist()
    out["F"] = Fs.tolist()
    lam2 = float(np.trapezoid(lams**3 * Fs, np.log(lams)))        # int lambda^2 F dlambda = int lambda^3 F dln(lambda)
    lam3 = float(np.trapezoid(lams**4 * Fs, np.log(lams)))
    out["Lambda2"] = lam2
    out["Lambda3"] = lam3
    out["F_at_lambda_max"] = float(Fs[-1])
    C = {uc: lam2 * uc**2 / (2 * math.pi**2) for uc in (0.5, 1.0, 2.0)}
    out["C_pairs_density_coeff"] = C
    m5 = lams <= 5.0
    lam3_5 = float(np.trapezoid((lams**4 * Fs)[m5], np.log(lams[m5])))
    Jf = lambda uc: quad(lambda u: u ** (4 / 3) * math.sqrt(1 + u * u), 0, uc)[0]
    out["Lambda3_cut5"] = lam3_5
    out["eps_tau4_coeff"] = {uc: (2 / math.pi**2) * lam3_5 * Jf(uc) for uc in (0.5, 1.0, 2.0)}
    out["mean_pair_energy_times_tau"] = {uc: out["eps_tau4_coeff"][uc] / C[uc] for uc in (0.5, 1.0, 2.0)}
    json.dump(out, open("kasner_shear_production.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    o = main()
    print(json.dumps({k: v for k, v in o.items() if k not in ("lambda_grid", "F")}, indent=1))
