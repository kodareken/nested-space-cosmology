"""Child energy at a torsion bounce versus the number of fermions it holds.

1. No-production baseline: N = M_BH / m_n neutrons, degenerate free neutron gas (full relativistic EoS),
   spin-torsion contact term eps_eff = eps - alpha n^2, alpha = kappa (hbar c)^2 / 32, bounce at eps = alpha n^2.
   Turnaround radius a0 solved from (kappa/3) eps(n0) a0^2 = 1 with n0 = N / ((4 pi / 3) a0^3).
2. Floor: any degenerate-fermion bounce at eps = alpha n^2 has energy per fermion alpha n_b = 3.99 E_Planck, so
   E_child >= 3.99 E_P * N whatever made the fermions.
3. Shear-produced pairs (kasner_stopping_scales.json): fermion number = 2 x pairs; per-fermion energy at the stage-2 bounce.
4. Thread length: at the singularity r* -> 0, so the Killing range equals the advanced-time range Delta v of the infall.
SI units. Nothing is fitted.
"""
import json
import math

import numpy as np
from scipy.optimize import brentq

G, c, hbar, Msun, mn = 6.67430e-11, 2.99792458e8, 1.054571817e-34, 1.98847e30, 1.67492750e-27
kappa = 8 * math.pi * G / c**4
alpha = kappa * (hbar * c) ** 2 / 32
lP = math.sqrt(hbar * G / c**3)
EP = math.sqrt(hbar * c**5 / G)
epsP = hbar * c / lP**4
lam_n = hbar / (mn * c)
SUN = Msun * c**2


def eps_fermi(n):
    kF = (3 * math.pi**2 * n) ** (1 / 3)
    x = kF * lam_n
    return mn**4 * c**5 / (8 * math.pi**2 * hbar**3) * (x * (2 * x * x + 1) * math.sqrt(1 + x * x) - math.asinh(x))


def bounce_density():
    return math.exp(brentq(lambda ln: math.log(eps_fermi(math.exp(ln))) - math.log(alpha) - 2 * ln,
                           math.log(1e60), math.log(1e120)))


def baseline(M_sun, nb):
    N0 = M_sun * Msun / mn
    rg = 2 * G * M_sun * Msun / c**2
    f = lambda la: (kappa / 3) * eps_fermi(N0 / ((4 * math.pi / 3) * math.exp(la) ** 3)) * math.exp(la) ** 2 - 1
    a0 = math.exp(brentq(f, math.log(rg * 1e-3), math.log(rg * 1e3)))
    ab = (N0 / nb / (4 * math.pi / 3)) ** (1 / 3)
    a = np.geomspace(a0 * (1 - 1e-9), ab * (1 + 1e-9), 100000)
    n = N0 / ((4 * math.pi / 3) * a**3)
    ad2 = np.array([(kappa / 3) * (eps_fermi(ni) - alpha * ni**2) * ai**2 - 1 for ni, ai in zip(n, a)])
    Eb = alpha * nb**2 * (4 * math.pi / 3) * ab**3
    return dict(N=N0, a0_over_rg=a0 / rg, a_b_m=ab, adot2_min=float(ad2.min()), E_b_solar=Eb / SUN,
                child_over_parent=Eb / (M_sun * SUN))


def main():
    nb = bounce_density()
    out = {"n_b_per_m3": nb, "n_b_per_lP3": nb * lP**3, "eps_b_over_epsP": alpha * nb**2 / epsP,
           "E_per_fermion_over_EP": alpha * nb / EP, "baseline": {str(m): baseline(m, nb) for m in (3, 10, 30)}}
    sc = json.load(open("kasner_stopping_scales.json"))
    eck = json.load(open("eck_collapse_consistency_check.json"))["10"]
    N_need = eck["N_needed_at_bounce"]
    out["N_need"] = N_need
    out["N_need_energy_paper_stage1_solar"] = N_need * mn / Msun
    out["floor_energy_at_N_need_solar"] = alpha * nb * N_need / SUN
    out["paper_M_b_solar"] = 1.7e31
    out["paper_energy_per_fermion_J"] = 1.7e31 * SUN / N_need
    t2M = 2 * G * 10 * Msun / c**3
    out["light_crossing_2GM_over_c3_s"] = t2M
    shear = {}
    for uc in ("0.5", "1.0", "2.0"):
        r = sc[f"species1.0_uc{uc}"]
        Nf = 2 * r["N_at_x_br"]
        shear[uc] = dict(pairs_br=r["N_at_x_br"], fermions_br=Nf, fermions_over_N_need=Nf / N_need,
                         thread_length_needed=N_need / Nf, delta_v_needed_ms=1e3 * t2M * N_need / Nf,
                         E_b_solar_L1=r["E_b_solar"], E_per_fermion_over_EP=r["E_b_solar"] * SUN / Nf / EP,
                         fermions_per_second_of_infall=Nf / t2M)
    out["shear_stage2"] = shear
    json.dump(out, open("child_energy_floor.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    o = main()
    print(json.dumps({k: v for k, v in o.items() if k != "baseline"}, indent=1, default=float))
