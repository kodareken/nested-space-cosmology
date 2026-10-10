"""Child lifetime versus net matter fraction, entanglement entropy of the produced pairs, and the baryon-number bound.

1. Lifetime (child proper time, bounce to recollapse), closed child, N fermions, net matter fraction eta:
       (da/d(ct))^2 = r_rad^2 / a^2 + r_m / a - 1,   r_m = 2 G eta N m_n / c^2,   r_rad^2 = (2G/c^4) N (3.99 E_P) a_b
   T = 2 t_half, t_half = (1/c) integral_0^{sqrt(a_max)} 2 (a_max - s^2) / sqrt(a_max - s^2 + a_2) ds, a_2 = r_rad^2 / a_max.
   Checks: eta = 0 gives 2 a_max / c; eta large gives pi a_max / c.
2. Entanglement entropy of the pairs if one partner is traced out: h(F) = -F ln F - (1-F) ln(1-F) per mode block, same density
   of states as the pair count, so S / N_pairs = Lambda_2[h] / Lambda_2[F] (cut at lambda <= 5; the same to 0.5% for cuts 2..10).
3. Baryon-number bound: child net matter eta N_f <= parent's baryon number M / m_n, with N_f = c1 (r_g / l_P)^2,
   so eta <= M* / M, M* = m_Pl^2 / (4 c1 m_n).
Run after child_energy_floor.py and kasner_stopping_scales.py. Writes child_lifetime_entropy.json and child_lifetime_entropy.png.
"""
import json
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

import child_energy_floor as cef

G, c, Msun, mn, EP, lP = cef.G, cef.c, cef.Msun, cef.mn, cef.EP, cef.lP
YR = 3.15576e7
nb = cef.bounce_density()
fb = cef.alpha * nb


def life(N, eta):
    ab = (3 * N / (4 * math.pi * nb)) ** (1 / 3)
    k = 2 * G / c**4
    rm = k * eta * N * mn * c**2
    rr2 = k * N * fb * ab
    amax = 0.5 * (rm + math.sqrt(rm * rm + 4 * rr2))
    a2 = rr2 / amax
    f = lambda s: 2 * (amax - s * s) / math.sqrt(amax - s * s + a2)
    th = quad(f, 0, math.sqrt(amax), limit=200, epsabs=0, epsrel=1e-10)[0]
    return 2 * th / c, amax


def main():
    res = json.load(open("kasner_shear_production.json"))
    lam, F = np.array(res["lambda_grid"]), np.array(res["F"])
    h = np.where((F > 0) & (F < 1), -F * np.log(np.clip(F, 1e-300, 1)) - (1 - F) * np.log(np.clip(1 - F, 1e-300, 1)), 0.0)
    cut = lambda f, lm: np.trapezoid((lam**3 * f)[lam <= lm], np.log(lam[lam <= lm]))
    out = {"S_over_Npairs": {str(lm): cut(h, lm) / cut(F, lm) for lm in (2, 3, 5, 10)}}
    sc = json.load(open("kasner_stopping_scales.json"))
    cen = sc["species1.0_uc1.0"]
    out["entropy_over_S_BH_at_x_br"] = {uc: out["S_over_Npairs"]["5"] * sc[f"species1.0_uc{uc}"]["N_at_x_br"] / sc[f"species1.0_uc{uc}"]["S_BH"]
                                        for uc in ("0.5", "1.0", "2.0")}
    out["infall_span_for_S_equals_S_BH_in_2M"] = {uc: 1 / v for uc, v in out["entropy_over_S_BH_at_x_br"].items()}
    N_need = json.load(open("child_energy_floor.json"))["N_need"]
    out["lifetime_yr_N_need"] = {str(e): life(N_need, e)[0] / YR for e in (0, 9e-21, 6e-10, 9e-8, 1e-6, 1e-3, 1.0)}
    c1 = 2 * cen["N_at_x_br"] / (2 * G * 10 * Msun / c**2 / lP) ** 2
    Mstar = (EP / c**2) ** 2 / (4 * c1 * mn)
    out.update({"c1": c1, "M_star_kg": Mstar,
                "eta_max": {str(m): Mstar / (m * Msun) for m in (3, 10, 30)},
                "parent_mass_for_eta_6e-10_kg": Mstar / 6e-10, "parent_mass_for_eta_9e-8_kg": Mstar / 9e-8,
                "torsion_energy_per_fermion_over_EP": 2 * fb / EP})
    json.dump(out, open("child_lifetime_entropy.json", "w"), indent=1)

    etas = np.logspace(-22, 0, 160)
    T = np.array([life(N_need, e)[0] / YR for e in etas])
    try:
        apply_figure_style()
    except NameError:
        pass
    fig, ax = plt.subplots(figsize=(6.6, 4.6), layout="constrained")
    ax.plot(etas, T, color="#1f5fa8", lw=2.2)
    for x, col, lab in ((out["eta_max"]["10"], "#c0392b", "most a 10 solar-mass parent can supply"),
                        (6e-10, "0.35", "baryon-to-photon ratio (order only)"),
                        (9e-8, "#2e7d32", "matter = radiation")):
        ax.axvline(x, color=col, lw=1.1, ls=":")
        ax.text(x * 0.7, 1.6e2, lab, color=col, rotation=90, ha="right", va="bottom", fontsize=8.3)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-22, 3); ax.set_ylim(30, 1e10)
    ax.set_xlabel(r"net matter fraction $\eta$ (surviving matter per fermion)")
    ax.set_ylabel("child lifetime, bounce to recollapse (years)")
    ax.set_title("How long the child lasts depends on its net matter", loc="left", fontsize=10)
    fig.savefig("child_lifetime_entropy.png", dpi=300)
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
