"""Net matter fraction eta from the Hehl-Datta level splitting (Poplawski, PRD 83, 084033, arXiv:1101.4012) at the torsion bounce.

Taken as the working assumption (the paper's Eqs. 14-17): fermion levels E = |p| + delta, antifermion levels E = |p| - delta,
with delta = (3/8) kappa n_aligned = 12 alpha n p_pol (alpha = kappa/32 in hbar = c = 1, p_pol = net polarisation fraction of the
background, n = fermion number density). Equilibrium asymmetry from Fermi-Dirac occupation (the paper's Eq. 16).
At the bounce alpha n_f equals the energy per fermion (eps = alpha n_f^2), and for a thermalised massless Fermi gas the mean
energy is 3.151 T, so delta / T = 12 x 3.151 p_pol = 37.8 p_pol, independent of the stage.
Parent angular momentum: J_ext = G M^2 / c, J / hbar = (M / m_P)^2 = R^2 / 4 with R = r_g / l_P; child fermions N_f = c1 L R^2,
so the polarisation that the parent's rotation (spin parameter a*) can supply is p_max = 2 J a* / (N_f hbar) = a* / (2 c1 L).
Hidden antimatter as dark matter (the paper's 4-5 GeV, taken as 4.5 GeV) adds dust mass (m_n + m_DM) per net unit.
Writes torsion_splitting_asymmetry.json. Run after child_energy_floor.py, kasner_stopping_scales.py, child_matter_fraction.py.
"""
import json
import math

import numpy as np
from scipy.integrate import quad

import child_energy_floor as cef
import child_lifetime_entropy as cle
import child_matter_fraction as cmf

G, c, hbar, Msun, mn, EP, lP = cef.G, cef.c, cef.hbar, cef.Msun, cef.mn, cef.EP, cef.lP
MEAN_E_OVER_T = 7 * math.pi**4 / (180 * 1.2020569031595942)       # <E>/T, massless Fermi-Dirac


def eta_equilibrium(x):
    """(n_antifermion - n_fermion) / (n_antifermion + n_fermion), levels p +- delta, x = delta / T, p >= 0, x < 1."""
    nX = quad(lambda p: p * p / (math.exp(p + x) + 1), 0, 60)[0]
    nXb = quad(lambda p: p * p / (math.exp(p - x) + 1), 0, 60)[0]
    return (nXb - nX) / (nXb + nX)


def main():
    out = {"mean_E_over_T": MEAN_E_OVER_T}
    xs = [0.01, 0.05, 0.2, 0.5, 0.9]
    out["eta_vs_delta_over_T"] = {str(x): eta_equilibrium(x) for x in xs}
    slope = eta_equilibrium(1e-3) / 1e-3
    out["eta_slope_per_delta_over_T"] = slope
    k = 12 * MEAN_E_OVER_T                                          # delta / T per unit polarisation
    out["delta_over_T_per_polarisation"] = k
    out["eta_per_polarisation"] = slope * k
    out["polarisation_needed"] = {str(e): e / (slope * k) for e in (6e-10, 1.56e-8, 9e-8, 1e-3)}
    sc = json.load(open("kasner_stopping_scales.json"))
    R = 2 * G * 10 * Msun / c**2 / lP
    out["c1"] = {uc: 2 * sc[f"species1.0_uc{uc}"]["N_at_x_br"] / R**2 for uc in ("0.5", "1.0", "2.0")}
    out["p_max_over_a_star_L1"] = {uc: 1 / (2 * v) for uc, v in out["c1"].items()}
    out["J_ext_over_hbar_10Msun"] = (10 * Msun / (EP / c**2)) ** 2
    p6 = out["polarisation_needed"]["6e-10"]
    out["spin_AM_needed_fraction_of_extremal_for_eta_6e-10_u_c_1"] = p6 * 2 * out["c1"]["1.0"] / 1.0   # J_spin / J_ext = (p N / 2) / (R^2/4)
    # hidden antimatter as dark matter: dust mass per net unit
    mDM = 4.5 * 1.602176634e-10 / c**2                              # kg
    dust = (mn + mDM) / mn
    out["dust_mass_factor_with_hidden_antimatter"] = dust
    out["eta_crossover_with_hidden_antimatter"] = cmf.turn(json.load(open("child_energy_floor.json"))["N_need"], 0)[0] / (2 * G * json.load(open("child_energy_floor.json"))["N_need"] * mn / c**2) / dust
    N = json.load(open("child_energy_floor.json"))["N_need"]
    out["turnaround_and_lifetime_N_need"] = {
        str(e): {"turnaround_solar": cmf.turn(N, e * dust)[1], "lifetime_yr": cle.life(N, e * dust)[0] / cle.YR}
        for e in (6e-10, 1.56e-8, 1e-6, 1e-3)}
    json.dump(out, open("torsion_splitting_asymmetry.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
