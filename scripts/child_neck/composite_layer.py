"""Composite neck layer: Dirac vacuum layer + torsion matter, matched to the thread path (sigma = -s, w = P/sigma = -1/2).   Geometric units (metres).

Thread path: |E| = sqrt(2R(M_p - m_-)), s = |E| / (4 pi R^2) (surface energy, magnitude), required surface pressure P = +s/2.
Vacuum layer (m species x one Dirac field), thickness d_v:  sigma_v = -x,  P_v = +x,  x = m (7 pi^2/2880) l_P^2 / d_v^3.
Torsion matter, thickness d_m, density n:  sigma_m = eps_eff(n) d_m,  P_m = p_eff(n) d_m  (w_m = p_eff/eps_eff).
Totals:  sigma_m = x - s,  P_m = s/2 - x   =>   w_m = (x - s/2)/(s - x).
  Branch A (x > s): sigma_m > 0, n < n_b, w_m < -1 (always reachable below n_b).     Branch B (3s/4 < x < s): sigma_m < 0, n > n_b, w_m > 1.
Thickness d_m = |sigma_m| / |eps_eff(n)| (geometric, G eps/c^4 = (eps/eps_b) EPSB / l_P^2).
Question: is there an x for which both thicknesses are at least one Planck length?  Reports the best min(d_v, d_m)/l_P over x, and the number of species m needed to reach 1.
Writes composite_layer.json and composite_layer.png.  Needs torsion_layer_eos.py, israel_matching.py, junction_neck.py and the restored child_* files.
"""
import json
import math

import matplotlib.pyplot as plt
import numpy as np

import child_energy_floor as cef
import israel_matching as im
import torsion_layer_eos as tl

lP, nb, eb, EPSB = cef.lP, tl.nb, tl.eb, tl.EPSB
C_DIRAC = 7 * math.pi**2 / 2880


def tables():
    nu_A = np.geomspace(1e-10, 0.99, 3000)                   # n = n_b (1 - nu)
    eA = np.array([tl.eff(nb * (1 - v)) for v in nu_A])
    nu_B = np.geomspace(1e-10, 1e3, 3000)                    # n = n_b (1 + nu)
    eB = np.array([tl.eff(nb * (1 + v)) for v in nu_B])
    # w_m and |eps_eff|/eps_b, monotone in nu
    wA, EA = eA[:, 1] / eA[:, 0], eA[:, 0] / eb
    wB, EB = eB[:, 1] / eB[:, 0], eB[:, 0] / eb
    return dict(nu_A=nu_A, wA=wA, EA=EA, nu_B=nu_B, wB=wB, EB=EB)


def best_thickness(s, m, T):
    """max over x of min(d_v, d_m)/l_P for the two branches; returns (value, branch, x/s)."""
    best = (0.0, None, None)
    # branch A: x = s (1 + delta)
    delta = np.geomspace(1e-9, 1e60, 1500)
    x = s * (1 + delta)
    wm = -(x - s / 2) / (x - s)                                # < -1
    iw = np.argsort(T["wA"])                                    # wA decreasing with nu_A decreasing: sort ascending for interp
    nu = np.exp(np.interp(wm, T["wA"][iw], np.log(T["nu_A"])[iw]))
    Eeff = np.exp(np.interp(np.log(nu), np.log(T["nu_A"]), np.log(T["EA"]))) * EPSB / lP**2
    dm = (x - s) / Eeff
    dv = (m * C_DIRAC * lP**2 / x) ** (1 / 3)
    v = np.minimum(dv, dm) / lP
    i = int(np.argmax(v))
    best = (float(v[i]), "A", float(x[i] / s))
    # branch B: x = s (1 - eps), eps in (0, 1/4)
    eps = np.geomspace(1e-9, 0.25, 600)
    x = s * (1 - eps)
    wm = (x - s / 2) / (s - x)                                  # > 1
    iw = np.argsort(T["wB"])
    nu = np.exp(np.interp(wm, T["wB"][iw], np.log(T["nu_B"])[iw]))
    Eeff = np.exp(np.interp(np.log(nu), np.log(T["nu_B"]), np.log(np.abs(T["EB"])))) * EPSB / lP**2
    dm = (s - x) / Eeff
    dv = (m * C_DIRAC * lP**2 / x) ** (1 / 3)
    v = np.minimum(dv, dm) / lP
    i = int(np.argmax(v))
    if v[i] > best[0]:
        best = (float(v[i]), "B", float(x[i] / s))
    return best


def species_needed(s, T):
    lo, hi = 0.0, 14.0
    if best_thickness(s, 10**hi, T)[0] < 1:
        return None
    if best_thickness(s, 1.0, T)[0] >= 1:
        return 1.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if best_thickness(s, 10**mid, T)[0] >= 1:
            hi = mid
        else:
            lo = mid
    return 10**hi


def main():
    rg = 2 * cef.G * 10 * cef.Msun / cef.c**2
    o = im.main()
    th = o["_th"]
    Rg, Eneed, _, _ = im.volume_curves(th, rg)
    s_of = lambda R: np.interp(R, Rg, Eneed) / (4 * np.pi * R**2)
    T = tables()
    Rs = np.geomspace(1e-30, rg * 0.99, 90)
    out = {"vacuum_coefficient": C_DIRAC, "eps_eff_max_below_nb_over_eps_b": float(T["EA"].max())}
    curves = {}
    for m in (1.0, 22.5, 1e3):
        curves[m] = [best_thickness(s_of(R), m, T) for R in Rs]
    need = [species_needed(s_of(R), T) for R in Rs]
    out["R_grid_m"] = list(map(float, Rs))
    out["best_min_thickness_over_lP"] = {f"m={m:g}": [c_[0] for c_ in v] for m, v in curves.items()}
    out["best_branch_m1"] = [c_[1] for c_ in curves[1.0]]
    out["species_needed"] = [None if n_ is None else float(n_) for n_ in need]
    pick = {}
    for R in (1e-26, 1e-24, 1e-22, 1e-18, 1.0, 1e3, rg * 0.99):
        v = best_thickness(s_of(R), 1.0, T)
        d = {"s_per_m": float(s_of(R)), "best_min_thickness_over_lP_m1": v[0], "branch": v[1], "x_over_s": v[2],
             "species_needed": None if species_needed(s_of(R), T) is None else float(species_needed(s_of(R), T))}
        mneed = d["species_needed"]
        if mneed is not None:
            vv = best_thickness(s_of(R), mneed * 1.0001, T)
            d["cancellation_x_over_s_at_species_needed"] = vv[2]
            # integrated normal stress of the composite: vacuum layer p_perp = 3 rho_v gives -3x, torsion layer (isotropic) gives P_m = s/2 - x  => -4x + s/2, against the shell's own scale s
            d["integrated_normal_stress_over_s"] = float((-4 * vv[2] + 0.5))
        pick[f"R={R:g}"] = d
    out["selected"] = pick
    valid = [(R, n_) for R, n_ in zip(Rs, need) if n_ is not None]
    out["common_thickness_plateau_over_lP"] = {f"m={m:g}": float((m * C_DIRAC / (T["EA"].max() * EPSB)) ** 0.25) for m in (1.0, 22.5, 123.0, 1e3)}
    out["species_floor_for_plateau_equal_lP"] = float(T["EA"].max() * EPSB / C_DIRAC)
    out["species_needed_min_over_R"] = float(min(n_ for _, n_ in valid)) if valid else None
    out["R_at_min_species_m"] = float(min(valid, key=lambda t: t[1])[0]) if valid else None
    json.dump(out, open("composite_layer.json", "w"), indent=1, default=float)
    out["_curves"], out["_need"], out["_Rs"] = curves, need, Rs
    return out


def figure(out, fname="composite_layer.png"):
    Rs, curves, need = out["_Rs"], out["_curves"], out["_need"]
    blue, red, grey, green = "#1f5fa8", "#c0392b", "0.4", "#2e8b57"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0))
    for m, col, lab in ((1.0, blue, "1 field"), (22.5, green, "22.5 fields"), (1e3, red, "1000 fields")):
        a1.loglog(Rs, [max(c_[0], 1e-60) for c_ in curves[m]], color=col, lw=2.5)
    a1.axhline(1.0, color=grey, ls=":", lw=1.2)
    a1.set_xlim(1e-30, 6e4)
    a1.set_ylim(1e-6, 10)
    a1.text(2e4, 1.9, "1000 fields", color=red, fontsize=8, ha="right", va="bottom")
    a1.text(2e4, 0.7, "22.5 fields", color=green, fontsize=8, ha="right", va="bottom")
    a1.text(2e4, 0.32, "1 field", color=blue, fontsize=8, ha="right", va="bottom")
    a1.text(1e-29, 1.25, "dotted: one\nPlanck length", color=grey, fontsize=8, ha="left", va="bottom")
    a1.set_xlabel("neck radius (m)")
    a1.set_ylabel("thickness of the thinner layer (Planck lengths)")
    a1.set_title("a   Best thickness of the thinner of the two layers", loc="left", fontsize=9)
    a2.loglog(Rs, [np.nan if n_ is None else n_ for n_ in need], color=blue, lw=2.5)
    a2.axhline(22.5, color=grey, ls=":", lw=1.2)
    a2.axhline(125.06, color=green, ls="--", lw=1.2)
    a2.set_xlim(1e-30, 6e4)
    a2.set_ylim(5, 1e12)
    a2.text(1e-29, 29, "dotted: 22.5 fields (the paper's count)", color=grey, fontsize=8, ha="left", va="bottom")
    a2.text(2e4, 170, "dashed: 125 fields, the floor above 1e-22 m", color=green, fontsize=8, ha="right", va="bottom")
    a2.text(1e-12, 1e6, "blue: fields needed", color=blue, fontsize=8, ha="left", va="center")
    a2.set_xlabel("neck radius (m)")
    a2.set_ylabel("Dirac fields needed for both layers to reach one Planck length")
    a2.set_title("b   Number of fields needed", loc="left", fontsize=9)
    fig.tight_layout()
    fig.savefig(fname, dpi=220)
    return fig


if __name__ == "__main__":
    o = main()
    print(json.dumps({k: v for k, v in o.items() if not k.startswith("_") and k not in ("R_grid_m", "best_min_thickness_over_lP", "best_branch_m1", "species_needed")}, indent=1, default=float))
