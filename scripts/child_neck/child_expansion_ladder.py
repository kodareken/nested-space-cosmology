"""After the torsion bounce: how far does the child expand, and what does 'its mass' mean at each epoch?

Closed child (k = +1), N fermions of mass m_n in a degenerate gas, adiabatic expansion at fixed N:
    adot^2 = (kappa/3) (eps(n) - alpha n^2) a^2 - 1,   n = N / ((4 pi / 3) a^3)   (dot = d/d(ct)).
Bounce at eps = alpha n^2 (a_b, shifted by ~1e-38 relative by the curvature term, neglected).
Turnaround where (kappa/3) eps a^2 = 1; there eps V = a c^4 / (2 G), i.e. the child's size equals its own
Schwarzschild radius: the same condition that defines the parent's turnaround a0 = r_g.
The Fermi-gas energy density uses the exact relativistic expression for x = k_F / (m c) > 0.05 and the series
1 + 3x^2/10 - 3x^4/56 below it (the exact form loses all precision there).
"""
import json
import math

from scipy.optimize import brentq

import child_energy_floor as cef

G, c, hbar, Msun, mn = cef.G, cef.c, cef.hbar, cef.Msun, cef.mn
kappa, alpha, lP, EP, SUN = cef.kappa, cef.alpha, cef.lP, cef.EP, cef.SUN


def eps(n):
    kF = (3 * math.pi**2 * n) ** (1 / 3)
    x = kF * cef.lam_n
    if x < 0.05:
        return n * mn * c**2 * (1 + 0.3 * x**2 - 3 * x**4 / 56)
    return cef.eps_fermi(n)


def ladder(N):
    nb = cef.bounce_density()
    ab = (3 * N / (4 * math.pi * nb)) ** (1 / 3)
    n_of = lambda a: N / ((4 / 3) * math.pi * a**3)
    per = lambda a: eps(n_of(a)) / n_of(a)                    # energy per fermion [J]
    F = lambda la: (kappa / 3) * (eps(n_of(math.exp(la))) - alpha * n_of(math.exp(la)) ** 2) * math.exp(la) ** 2 - 1
    la_max = brentq(F, math.log(ab * 30), math.log(ab * 1e45), xtol=1e-12)
    amax = math.exp(la_max)
    out = {"N": N, "a_b_m": ab, "a_max_m": amax, "a_max_over_a_b": amax / ab,
           "M_turn_solar": amax * c**2 / (2 * G) / Msun, "M_turn_over_N_mn": amax * c**2 / (2 * G) / (N * mn),
           "per_fermion_at_bounce_J": alpha * nb, "per_fermion_at_turn_J": per(amax)}
    for label, target in (("paper_1p67J", 1.672), ("two_rest_energies", 2 * mn * c**2)):
        la = brentq(lambda l: math.log(per(math.exp(l))) - math.log(target), math.log(ab * 1.01), la_max)
        out[f"a_over_a_b_at_{label}"] = math.exp(la) / ab
        out[f"a_m_at_{label}"] = math.exp(la)
    # accelerating phase: q = -(1/6)(kappa)(eps~ + 3 p~) sign; for the degenerate/radiation limit eps~ = eps - alpha n^2, p~ = eps/3 - alpha n^2
    # a'' > 0 for alpha n^2 > eps/2, i.e. only for a / a_b < sqrt(2) (ultra-relativistic limit)
    return out


def main():
    f = json.load(open("child_energy_floor.json"))
    cases = {"baseline_10_solar": f["baseline"]["10"]["N"], "shear_u_c_1_L_1": f["shear_stage2"]["1.0"]["fermions_br"],
             "N_need": f["N_need"], "shear_u_c_1_L_150": 150 * f["shear_stage2"]["1.0"]["fermions_br"]}
    res = {k: ladder(v) for k, v in cases.items()}
    # analytic ultra-relativistic turnaround for the N_need child, for comparison
    nb = cef.bounce_density()
    eb = alpha * nb**2
    a_b = (3 * f["N_need"] / (4 * math.pi * nb)) ** (1 / 3)
    res["N_need"]["a_max_massless_analytic_m"] = math.sqrt(kappa * eb / 3) * a_b**2
    json.dump(res, open("child_expansion_ladder.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    r = main()
    for k, v in r.items():
        print(k, {a: float(f"{b:.4g}") for a, b in v.items()})
