"""What sets the child's thread length, and what the generation map then does.

1. Entropy cap on L (thread length in units of 2GM/c^3). Thermal upper bound on the child's entropy: s/n = (4/3) <E>/T per fermion
   = 4.20 for a massless Fermi-Dirac gas, so S_th = 4.20 N_f. Require S_th <= S_BH(parent) = pi R^2, R = r_g / l_P.
   N_f = c1 L R^2 (shear production at x_br, c1 = 2 N_br / R^2 per u_c).  The paper's own N_need = sqrt(8/27) R^2 gives S_th/S_BH = 4.20 x 0.1733 = 0.73 for any u_c.
2. Closure of a closed (S^3) radiation child: isotropic collapse until cT* = x a_max, then Kasner, V proportional to tau.
   a(eta) = a_max sin eta, c tau = a_max (1 + cos eta), so a_* / a_max = sqrt(x (2 - x)), V_* = 2 pi^2 a_*^3, and
   L_c = V / (6 pi r_g,c^2 c tau) = (pi/3) [x (2 - x)]^(3/2) / x with r_g,c = a_max.
3. Generation map M_{n+1} = turnaround mass of a closed child holding N_f(M_n, L) fermions (child_matter_fraction.turn), for all
   fermions surviving as matter (eta = 1) and for radiation (eta = 0); N_f / N_need and r_g(M_{n+1}) / r_g(M_n).
4. Parent mass whose child has the Hubble mass M_H = c^3 / (2 G H0), H0 = 70 km/s/Mpc, and the L that keeps the second generation
   inside the observable radius 4.4e26 m.
Run after restoring the child_* scripts and JSON files. Writes nesting_thread_length.json and nesting_thread_length.png.
"""
import json
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

import child_energy_floor as cef
import child_matter_fraction as cmf

G, c, Msun, mn, lP = cef.G, cef.c, cef.Msun, cef.mn, cef.lP
S_PER_N = (4 / 3) * 7 * math.pi**4 / (180 * 1.2020569031595942)       # 4.20 for massless Fermi-Dirac
H0 = 70e3 / 3.0856775814913673e22
M_H = c**3 / (2 * G * H0) / Msun
R_OBS = 4.4e26
SQ = math.sqrt(8 / 27)


def c1_of(uc):
    sc = json.load(open("kasner_stopping_scales.json"))[f"species1.0_uc{uc}"]
    R = 2 * G * 10 * Msun / c**2 / lP
    return 2 * sc["N_at_x_br"] / R**2


def child_mass(M_sun, L, uc, eta):
    R = 2 * G * M_sun * Msun / c**2 / lP
    N = c1_of(uc) * L * R**2
    return cmf.turn(N, eta)[1], N, R


def solve_parent_for(M_target, L, uc, eta):
    f = lambda lm: math.log(child_mass(math.exp(lm), L, uc, eta)[0]) - math.log(M_target)
    return math.exp(brentq(f, math.log(1e-14), math.log(1e12)))


def main():
    out = {"S_per_fermion": S_PER_N, "M_H_solar": M_H, "H0_per_s": H0}
    out["c1"] = {uc: c1_of(uc) for uc in ("0.5", "1.0", "2.0")}
    out["L_max_thermal_entropy"] = {uc: math.pi / (S_PER_N * v) for uc, v in out["c1"].items()}
    ent = json.load(open("child_lifetime_entropy.json"))["infall_span_for_S_equals_S_BH_in_2M"]
    out["L_max_entanglement_entropy"] = ent
    out["L_needed_for_N_need"] = {uc: SQ / v for uc, v in out["c1"].items()}
    out["S_th_over_S_BH_at_N_need"] = S_PER_N * SQ / math.pi
    out["L_closure_radiation_child"] = {str(x): (math.pi / 3) * (x * (2 - x)) ** 1.5 / x for x in (0.01, 0.1, 0.3, 0.5, 1.0)}
    # generation maps
    gens = {}
    for label, eta in (("matter", 1.0), ("radiation", 0.0)):
        for L in (0.3, 1.0, 1.37):
            M, seq, ratio_N, omega = 10.0, [10.0], [], []
            for _ in range(2):
                Mn, N, R = child_mass(M, L, "1.0", eta)
                ratio_N.append(N / (SQ * R**2))
                omega.append(Mn / M)
                M = Mn
                seq.append(M)
            gens[f"{label}_L{L}"] = {"M_solar": seq, "N_over_N_need": ratio_N, "r_g_ratio_next_over_previous": omega,
                                     "r_g_last_m": 2 * G * seq[-1] * Msun / c**2}
    out["generations_u_c_1"] = gens
    out["parent_mass_for_Hubble_mass_child"] = {
        f"{lab}_u_c{uc}_L{L}": solve_parent_for(M_H, L, uc, eta)
        for lab, eta in (("matter", 1.0), ("matter_with_hidden_antimatter_eta_6e-10", 5.79 * 6e-10), ("radiation", 0.0))
        for uc in ("0.5", "1.0", "2.0") for L in (0.3, 1.0)}
    out["lifetime_Hubble_mass_matter_child_Gyr"] = math.pi / H0 / 3.15576e7 / 1e9
    # L needed to keep generation 2 inside the observable radius
    for label, eta in (("radiation", 0.0), ("matter", 1.0)):
        def r2(lL):
            L = math.exp(lL)
            M1 = child_mass(10.0, L, "1.0", eta)[0]
            M2 = child_mass(M1, L, "1.0", eta)[0]
            return math.log(2 * G * M2 * Msun / c**2) - math.log(R_OBS)
        out[f"L_for_generation2_inside_observable_{label}"] = math.exp(brentq(r2, math.log(1e-120), math.log(1.0)))
    json.dump(out, open("nesting_thread_length.json", "w"), indent=1)
    return out


def figure():
    Mp = np.logspace(-38, 12, 260)
    curve = lambda L, eta, uc="1.0": np.array([child_mass(m, L, uc, eta)[0] for m in Mp])
    fig, ax = plt.subplots(figsize=(6.8, 4.9), layout="constrained")
    blue, red, green, grey = "#1f5fa8", "#c0392b", "#2e7d32", "0.4"
    ax.fill_between(Mp, curve(0.3, 1.0), curve(1.4, 1.0), color=green, alpha=0.18, lw=0)
    ax.fill_between(Mp, curve(0.3, 0.0), curve(1.4, 0.0), color=blue, alpha=0.18, lw=0)
    ax.plot(Mp, curve(1.0, 1.0), color=green, lw=2.2)
    ax.plot(Mp, curve(1.0, 0.0), color=blue, lw=2.2)
    ax.plot(Mp, Mp, color=grey, lw=1.3, ls="--")
    ax.axhline(M_H, color=red, lw=1.2, ls=":")
    ax.plot([10], [child_mass(10, 1.0, "1.0", 1.0)[0]], "o", color=green, mfc="white", mew=1.6, ms=6)
    ax.plot([10], [child_mass(10, 1.0, "1.0", 0.0)[0]], "o", color=blue, mfc="white", mew=1.6, ms=6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-38, 1e12); ax.set_ylim(1e-38, 1e40)
    ax.set_xlabel("parent mass (solar masses)"); ax.set_ylabel("child mass at turnaround (solar masses)")
    from matplotlib.lines import Line2D
    h = [Line2D([0], [0], color=green, lw=2.2), Line2D([0], [0], color=blue, lw=2.2),
         Line2D([0], [0], color=grey, lw=1.3, ls="--"), Line2D([0], [0], color=red, lw=1.2, ls=":")]
    ax.legend(h, ["every fermion survives as matter (band: thread length 0.3 to 1.4)", "radiation only (same band)",
                  "child = parent", "Hubble mass, 4.4e22"], loc="upper left", frameon=False, fontsize=8.2)
    ax.set_title("From the Planck mass upward every generation outweighs its parent", loc="left", fontsize=10)
    fig.savefig("nesting_thread_length.png", dpi=300)
    rend = fig.canvas.get_renderer(); lb = ax.get_legend().get_window_extent(rend)
    pts = sum(int(((p[:, 0] >= lb.x0) & (p[:, 0] <= lb.x1) & (p[:, 1] >= lb.y0) & (p[:, 1] <= lb.y1)).sum())
              for p in (ax.transData.transform(np.column_stack([Mp, y])) for y in (curve(1.0, 1.0), curve(1.0, 0.0), Mp)))
    return pts


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
    print("curve points inside legend:", figure())
