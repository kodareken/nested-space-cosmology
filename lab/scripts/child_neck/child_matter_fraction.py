"""Child turnaround mass versus the fraction of fermions that survive as massive matter, and a two-generation map.

Closed child with N fermions released at the torsion bounce (a_b from n_b = 5.08 l_P^-3, E = 3.99 E_P per fermion).
Radiation energy falls as 1/a. A fraction eta of the fermions is net matter of mass m_n and does not redshift.
Turnaround where the child's size equals its own Schwarzschild radius, a = 2 G E / c^4:
    a^2 - r_m a - r_rad^2 = 0,   r_m = 2 G eta N m_n / c^2,   r_rad^2 = (2G/c^4) (N 3.99 E_P) a_b.
Generation map: N_f = c1 (r_g / l_P)^2 with c1 = 0.3975 (shear production, u_c = 1, thread length 2M).
Writes child_matter_fraction.json and child_expansion_ladder.png. Run after child_energy_floor.py.
"""
import json
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import child_energy_floor as cef
import child_expansion_ladder as cx

G, c, Msun, mn, EP, SUN, lP = cef.G, cef.c, cef.Msun, cef.mn, cef.EP, cef.SUN, cef.lP
nb = cef.bounce_density()
fb = cef.alpha * nb


def turn(N, eta):
    ab = (3 * N / (4 * math.pi * nb)) ** (1 / 3)
    k = 2 * G / c**4
    rm = k * eta * N * mn * c**2
    rr2 = k * (N * fb) * ab
    a = 0.5 * (rm + math.sqrt(rm * rm + 4 * rr2))
    return a, a * c**2 / (2 * G) / Msun, ab


def main():
    f = json.load(open("child_energy_floor.json"))
    N = f["N_need"]
    out = {"eta_table": {}, "generations": {}}
    for eta in (0, 1e-10, 6e-10, 1e-8, 1e-7, 1e-6, 1e-3, 1):
        a, M, ab = turn(N, eta)
        out["eta_table"][str(eta)] = {"a_max_m": a, "M_turn_solar": M}
    a0, _, ab = turn(N, 0)
    out["eta_crossover_matter_equals_radiation"] = a0 / (2 * G * N * mn / c**2)
    c1 = f["shear_stage2"]["1.0"]["fermions_br"] / (2 * G * 10 * Msun / c**2 / lP) ** 2
    out["c1"] = c1
    for label, eta in (("matter", 1.0), ("radiation", 0.0)):
        M, seq = 10.0, [10.0]
        for _ in range(3):
            Nf = c1 * (2 * G * M * Msun / c**2 / lP) ** 2
            M = turn(Nf, eta)[1]
            seq.append(M)
        out["generations"][label] = seq
    K = c1 * mn * (2 * G / c**2 / lP) ** 2
    out["matter_map_fixed_point_kg"] = 1 / K
    out["m_Pl_squared_over_m_n_kg"] = (EP / c**2) ** 2 / mn
    json.dump(out, open("child_matter_fraction.json", "w"), indent=1)

    ab = (3 * N / (4 * math.pi * nb)) ** (1 / 3)
    a = np.geomspace(ab, 3e27, 600)
    n_of = lambda x: N / ((4 / 3) * math.pi * x**3)
    M1 = np.array([N * cx.eps(n_of(x)) / n_of(x) for x in a]) / SUN
    Mrad = N * fb * (ab / a) / SUN
    mix = lambda eta: Mrad + eta * N * mn * c**2 / SUN
    line = a * c**2 / (2 * G) / Msun
    blue, red, green, grey = "#1f5fa8", "#c0392b", "#2e7d32", "0.4"
    try:
        apply_figure_style()
    except NameError:
        pass
    fig, ax = plt.subplots(figsize=(6.8, 4.9), layout="constrained")
    ax.plot(a, line, color=grey, lw=1.3, ls="--")
    ax.plot(a, M1, color=green, lw=2.2)
    ax.plot(a, mix(1e-7), color=blue, lw=1.8)
    ax.plot(a, mix(6e-10), color=red, lw=2.0)
    for eta, col in ((1.0, green), (1e-7, blue), (6e-10, red)):
        at, Mt, _ = turn(N, eta)
        ax.plot([at], [Mt], "o", color=col, ms=6, mfc="white", mew=1.6)
    a_pap = ab * 4666e6
    ax.plot([a_pap], [N * 1.672 / SUN], "D", color="k", ms=6)
    ax.plot([ab], [N * fb / SUN], "o", color="k", ms=5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(ab * 0.8, 3e27); ax.set_ylim(1e12, 3e43)
    ax.set_xlabel("size of the child, a (m)"); ax.set_ylabel("mass-energy of the child (solar masses)")
    ax.text(2.5e-9, 3e26, "dot: the bounce,\n3.99 Planck energies\nper fermion", color="0.2", va="center", ha="left", fontsize=8.5)
    ax.text(a_pap * 4, N * 1.672 / SUN * 1500, "paper's recorded child mass\n(size 3 m, 1.67 J per fermion)", color="k", ha="left", va="center")
    h = [Line2D([0], [0], color=green, lw=2.2), Line2D([0], [0], color=blue, lw=1.8),
         Line2D([0], [0], color=red, lw=2.0), Line2D([0], [0], color=grey, lw=1.3, ls="--")]
    ax.legend(h, ["every fermion survives as massive matter: turnaround 1.5e21",
                  "1e-7 of the fermions survive: turnaround 2.3e14",
                  "6e-10 or fewer survive: turnaround 1.4e14",
                  "child size = its own Schwarzschild radius"], loc="upper right", frameon=False, fontsize=8.2, handlelength=2.2)
    ax.set_title("What the child weighs depends on when it is read and what it holds", loc="left", fontsize=10)
    fig.savefig("child_expansion_ladder.png", dpi=300)
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
