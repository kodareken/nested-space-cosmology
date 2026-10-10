"""Self-consistency check of the Einstein-Cartan black-hole-universe collapse in
Poplawski, arXiv:1110.5019 (Eqs. 1-4), using the paper's own assumptions.

Assumption taken from the paper: fermion number grows as N ~ a^-3 during the
collapse (stiff matter, eps ~ a^-6, n ~ a^-6). The question asked here is
whether the dynamics that assumption implies would allow the production that
assumption needs. Nothing is fitted. SI units; dot = d/d(ct).

Step-by-step substitution of Eqs. 1-4 is used, not the paper's closed forms
(5)-(9): Eqs. 5-6 agree with direct substitution only if m_Pl in them is the
non-reduced Planck mass, and the constant in Eq. 8 does not reproduce.
Scalings are therefore trusted; prefactors are not.
"""
import math

import numpy as np

G = 6.67430e-11
c = 2.99792458e8
hbar = 1.054571817e-34
Msun = 1.98847e30
mn = 1.67492750e-27
kappa = 8 * math.pi * G / c**4
omega_n = mn * c**2 / hbar            # neutron Compton frequency [1/s]


def setup(m_bh_sun, f=1.0):
    """Initial data at the turnaround a0 = r_g and the bounce scale a* (Eqs. 1-4)."""
    M = m_bh_sun * Msun
    a0 = 2 * G * M / c**2
    eps0 = 3 * M * c**2 / (4 * math.pi * a0**3)
    n0 = 9 * M / (4 * math.pi * a0**3 * mn)
    alpha = f * kappa * (hbar * c) ** 2 / 32
    astar = (alpha * n0**2 / eps0) ** (1 / 6) * a0

    def adot2(a):                       # Eq. (4) with k = 1
        return (kappa / 3) * (eps0 * a0**6 / a**4 - alpha * n0**2 * a0**12 / a**10) - 1.0

    return dict(M=M, a0=a0, eps0=eps0, n0=n0, alpha=alpha, astar=astar, adot2=adot2)


def profile(m_bh_sun, npts=800001):
    s = setup(m_bh_sun)
    a = np.geomspace(s["astar"] * (1 + 1e-9), s["a0"] * (1 - 1e-9), npts)
    ad2 = np.array([s["adot2"](x) for x in a])
    ok = ad2 > 0
    a, ad = a[ok], np.sqrt(ad2[ok])
    return s, a, ad, c * ad / a


def check(m_bh_sun=10):
    s, a, ad, H = profile(m_bh_sun)
    i = int(np.argmax(H))
    n_need = (s["M"] / mn) * (s["a0"] / s["astar"]) ** 3
    n_max = (2 / 3) * ad[i] ** 3        # comoving modes with k_phys <~ H; Pauli: occupation <= 1
    n_b = s["n0"] * (s["a0"] / s["astar"]) ** 6
    eps_b = s["eps0"] * (s["a0"] / s["astar"]) ** 6
    kF = (3 * math.pi**2 * n_b) ** (1 / 3)
    eps_fermi = 0.75 * hbar * c * kF * n_b
    return {
        "M_BH_solar": m_bh_sun,
        "omega_n_per_s": omega_n,
        "H_max_per_s": float(H[i]),
        "H_max_over_omega_n": float(H[i] / omega_n),
        "a_at_H_max_over_astar": float(a[i] / s["astar"]),
        "exp_minus_pi_omega_over_H": math.exp(-math.pi * omega_n / H[i]),
        "N_initial": s["M"] / mn,
        "N_needed_at_bounce": n_need,
        "N_max_excitable_modes": n_max,
        "gap_orders_of_magnitude": math.log10(n_need / n_max),
        "n_bounce_per_m3": n_b,
        "energy_per_fermion_J": eps_b / n_b,
        "hbar_c_kF_J": hbar * c * kF,
        "eps_bounce_J_per_m3": eps_b,
        "eps_min_fermi_J_per_m3": eps_fermi,
        "fermi_over_eps_orders": math.log10(eps_fermi / eps_b),
    }


if __name__ == "__main__":
    import json
    out = {m: check(m) for m in (3, 10, 30)}
    json.dump(out, open("eck_collapse_consistency_check.json", "w"), indent=1)
    for m, r in out.items():
        print(m, {k: float(f"{v:.4g}") for k, v in r.items()})
