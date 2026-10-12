"""Can torsion matter above the bounce density supply the neck's negative-energy layer?   (israel-matching-note.md, section 4)

Torsion fermion fluid: eps_eff = eps_F(n) - alpha n^2,  p_eff = p_F(n) - alpha n^2  (contact term p = n d(eps)/dn - eps = -alpha n^2).
Layer needed by the junction (israel_matching.py): negative surface energy sigma < 0 with surface pressure P; thread path w = P/sigma = -1/2 for R << r_g, mass-matched path w from +4 to -1/2.
A perfect-fluid layer has T_parallel = p_eff and surface quantities sigma = eps_eff d, P = p_eff d, so w = p_eff/eps_eff.
Tests: (1) region of the (eps_eff, p_eff) plane the fluid reaches; (2) w(n) above n_b; (3) layer thickness d = |E|/(4 pi R^2 |eps_eff|) against the Planck length and R;
(4) dynamical time of the layer, 1/sqrt((8 pi/3) G eps_b/c^2), against the proper time the thread path needs the layer to last; (5) mass-matched-path window where w >= 1.
Writes torsion_layer_eos.json and torsion_layer_eos.png.  Needs israel_matching.py, junction_neck.py and the restored child_* files.
"""
import json
import math

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.optimize import brentq

import child_energy_floor as cef
import israel_matching as im

G, c, hbar, lP, mn, alpha = cef.G, cef.c, cef.hbar, cef.lP, cef.mn, cef.alpha
nb = cef.bounce_density()
eb = cef.eps_fermi(nb)
EPSB = im.EPSB                                              # eps_b in Planck densities
TB = lP / (c * math.sqrt(8 * math.pi / 3 * EPSB))           # torsion dynamical time at the bounce density


def p_fermi(n):
    x = (3 * math.pi**2 * n) ** (1 / 3) * cef.lam_n
    return mn**4 * c**5 / (24 * math.pi**2 * hbar**3) * (x * (2 * x * x - 3) * math.sqrt(1 + x * x) + 3 * math.asinh(x))


def eff(n):
    return cef.eps_fermi(n) - alpha * n**2, p_fermi(n) - alpha * n**2


def w_of_f(f):
    e, p = eff(f * nb)
    return p / e


def main():
    out = {"n_b_per_m3": nb, "eps_b_J_per_m3": eb, "eps_b_over_eps_P": EPSB, "torsion_dynamical_time_s": TB}
    # (0) pressure formula check: p = n d eps/dn - eps
    h = 1e-6
    n0 = 3 * nb
    dnum = (cef.eps_fermi(n0 * (1 + h)) - cef.eps_fermi(n0 * (1 - h))) / (2 * n0 * h)
    out["pressure_formula_rel_error"] = abs(n0 * dnum - cef.eps_fermi(n0) - p_fermi(n0)) / p_fermi(n0)
    # (1) the (eps_eff, p_eff) values the fluid reaches
    fs = np.concatenate([np.linspace(0.05, 0.999, 300), np.geomspace(1.0001, 1e3, 600)])
    E_, P_ = np.array([eff(f * nb) for f in fs]).T / eb
    neg = E_ < 0
    out["fluid_eps_eff_negative_means"] = dict(max_p_eff_over_eps_b=float(P_[neg].max()), min_w=float((P_[neg] / E_[neg]).min()), w_at_2nb=w_of_f(2.0), w_at_5nb=w_of_f(5.0), w_at_1p01nb=w_of_f(1.01),
                                              w_at_100nb=w_of_f(100.0), eps_eff_over_eb_at_2nb=eff(2 * nb)[0] / eb, p_eff_over_eb_at_2nb=eff(2 * nb)[1] / eb)
    out["NEC_eps_plus_p_over_eb"] = {f"{f}nb": (eff(f * nb)[0] + eff(f * nb)[1]) / eb for f in (1.0, 1.5, 2.0, 5.0)}
    # w = -1/2 would need n/n_b = ... with eps_eff > 0 (wrong sign of sigma)
    out["w_minus_half_at_n_over_nb"] = brentq(lambda f: w_of_f(f) + 0.5, 0.2, 0.9)
    # (3) thickness of the layer on the thread path
    th = im.thread_neck(2 * G * 10 * cef.Msun / c**2, math.sqrt(im.jn.child_params(json.load(open("child_energy_floor.json"))["N_need"], 0.0)[1]))
    rg = 2 * G * 10 * cef.Msun / c**2
    Rg, Eneed, _, _ = im.volume_curves(th, rg)
    rows = []
    for f in (1.5, 2.0, 5.0):
        epsg = abs(eff(f * nb)[0]) / eb * EPSB / lP**2       # |eps_eff| in geometric units (1/m^2)
        d = Eneed / (4 * np.pi * Rg**2 * epsg)
        Rthin = float(np.exp(np.interp(0.0, np.log(d / Rg)[::-1], np.log(Rg)[::-1])))
        Rplanck = float(np.exp(np.interp(0.0, np.log(d / lP)[::-1], np.log(Rg)[::-1])))
        rows.append(dict(n_over_nb=f, thickness_at_1e_minus_23_m=float(np.interp(1e-23, Rg, d)), thickness_at_1m=float(np.interp(1.0, Rg, d)), thickness_at_r_g=float(d[-1]),
                         R_where_thickness_equals_R_m=Rthin, R_where_thickness_equals_lP_m=Rplanck))
    out["thickness"] = rows
    # thick filling: density excess needed so that the whole neck volume carries E
    k = -(eff(nb * 1.0001)[0] / eb) / 1e-4
    out["eps_eff_slope_per_unit_delta"] = k
    out["thick_fill_delta"] = {f"R={R:g}": float(np.interp(R, Rg, Eneed) / ((4 * math.pi / 3) * R**3 * EPSB / lP**2 * k)) for R in (1e-20, 1.0, rg)}
    # (4) lifetime needed: parent proper time from the singularity, tau(R) on the thread
    tau = th["tau"] / c
    out["lifetime_ratio"] = {f"R={R:g}": float(np.interp(R, th["R"], tau) / TB) for R in (1e-22, 1e-18, 1.0, 1e3, rg * 0.99)}
    out["thread_total_proper_time_s"] = float(tau[-1])
    out["lifetime_ratio_total"] = float(tau[-1] / TB)
    R_tb = float(np.interp(TB, tau, th["R"]))
    out["R_reached_in_one_dynamical_time_m"] = R_tb
    # (5) mass-matched path: window where w >= 1
    N = json.load(open("child_energy_floor.json"))["N_need"]
    F_open = im.build(True)
    ab, rr2, rm = im.jn.child_params(N, 0.0)
    aex = im.jn.exit_point(10.0, N, 0.0)["a_exit_m"]
    args = lambda av: (mp.mpf(av), mp.mpf(rr2), mp.mpf(rm), mp.mpf(ab), mp.mpf(rg))
    facs = np.geomspace(0.1935, 0.45, 160)
    W, RR, Es = [], [], []
    for fac in facs:
        sg, P, Rj = F_open["sigma"](*args(fac * aex)), F_open["P"](*args(fac * aex)), F_open["Rj"](*args(fac * aex))
        W.append(float(P / sg)); RR.append(float(Rj / rg)); Es.append(float(-4 * mp.pi * Rj**2 * sg))
    W, RR, Es = np.array(W), np.array(RR), np.array(Es)
    ia = int(np.argmax(W < 1.0))
    fac_w1 = float(np.exp(np.interp(1.0, W[::-1], np.log(facs)[::-1])))
    t_of_a = lambda f: f**2 * rg / 2                         # radiation child: t = a^2/(2 r_rad), exit at t = r_g/2
    win = dict(a_over_aexit_start=float(facs[0]), a_over_aexit_w_equals_1=fac_w1, R_over_rg_w_equals_1=float(np.interp(1.0, W[::-1], RR[::-1])),
               window_duration_s=float((t_of_a(fac_w1) - t_of_a(3 ** -1.5)) / c), duration_over_torsion_time=float((t_of_a(fac_w1) - t_of_a(3 ** -1.5)) / c / TB))
    f_w = [2.0, 4.0]
    win["n_over_nb_for_w"] = {f"w={w}": brentq(lambda f: w_of_f(f) - w, 1.0000001, 1e6) for w in f_w}
    out["mass_matched_window_w_ge_1"] = win
    json.dump(out, open("torsion_layer_eos.json", "w"), indent=1, default=float)
    out["_th"], out["_ms"] = th, (RR, W)
    out["_curve"] = (fs, E_, P_)
    return out


def figure(out, fname="torsion_layer_eos.png"):
    th = out["_th"]
    rg = 2 * G * 10 * cef.Msun / c**2
    blue, red, grey = "#1f5fa8", "#c0392b", "0.4"
    RR, W = out["_ms"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0))
    sel = th["R"] < 0.8 * rg
    a1.semilogx(th["R"][sel] / rg, (th["P"] / th["sig"])[sel], color=blue, lw=2.5)
    a1.semilogx(RR, W, color=red, lw=2.5)
    a1.axhspan(1.0, 6.0, color="0.88", zorder=0)
    a1.set_xlim(1e-8, 8)
    a1.set_ylim(-1.0, 6.0)
    a1.set_xlabel("neck radius (units of the parent's horizon radius)")
    a1.set_ylabel("pressure / energy density of the layer")
    a1.set_title("a   Layer the neck needs, and what torsion matter can give", loc="left", fontsize=9)
    a1.text(2e-8, 3.5, "grey: all that torsion matter above the bounce\ndensity can give (ratio of at least 1; both are negative)", color=grey, fontsize=8, ha="left", va="center")
    a1.text(3e-8, -0.2, "blue: neck follows the parent's thread (-1/2)", color=blue, fontsize=8, ha="left", va="center")
    a1.text(0.3, 2.6, "red: mass-\nmatched\nneck, from\nr_g/9", color=red, fontsize=8, ha="left", va="center")
    tau = th["tau"] / c
    ok = th["R"] > 1e-24
    a2.loglog(th["R"][ok], tau[ok] / TB, color=blue, lw=2.5)
    a2.axhline(1.0, color=grey, ls=":", lw=1.2)
    a2.set_ylim(1e-4, 1e42)
    a2.set_xlabel("neck radius (m)")
    a2.set_ylabel("time the layer must last / its own dynamical time")
    a2.set_title("b   Does the layer last long enough?", loc="left", fontsize=9)
    a2.text(0.98, 0.2, "dotted line: dynamical time of matter\nat the bounce density (4e-45 s)", transform=a2.transAxes, color=grey, fontsize=8, ha="right", va="center")
    fig.tight_layout()
    fig.savefig(fname, dpi=220)
    return fig


if __name__ == "__main__":
    o = main()
    print(json.dumps({k: v for k, v in o.items() if not k.startswith("_")}, indent=1, default=float))
