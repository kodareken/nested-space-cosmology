"""A closed child glued to a Schwarzschild exterior of mass M_p through a neck (standard GR junction, Misner-Sharp mass).

Closed FRW, areal radius R = a sin(chi). The Misner-Sharp mass of the sphere at (t, chi) is
    2 G m / c^2 = a sin^3(chi) (1 + adot^2),   adot = da / d(ct),
which is purely geometric, so it holds with the torsion-modified Friedmann equation too. Requiring m = M_p (the exterior mass seen
by the parent's outside) fixes the neck: sin^3(chi) = r_gp / (a (1 + adot^2)), R_j^3 = r_gp a^2 / (1 + adot^2).
The neck is hidden (inside the parent's horizon, R_j < r_gp) iff  a < r_gp sqrt(1 + adot^2).
Child dynamics after the bounce:  adot^2 = r_rad^2/a^2 + r_m/a - 1 - (r_rad^2/a_b^2)(a_b/a)^4  (the last term is the torsion contact term;
it cancels the first at a_b), with r_rad^2 = (2G/c^4) N (3.99 E_P) a_b and r_m = 2 G eta_eff N m_n / c^2 as in child_matter_fraction.
CORRECTION (7 Oct 2026, israel_matching.py): the "neck at the bounce" R_b = a_b (M_p / M_b)^(1/3) below used the bare density; with the geometric mass (torsion included) the child has
almost no effective mass at a_b, and the mass-matched neck is a SPACELIKE surface for t < t_exit/27 (a < 0.192 a_exit, R < r_g/9). The neck_at_bounce_* entries are superseded.
Writes junction_neck.json and junction_neck.png. Needs the restored child_* scripts and JSON files.
"""
import json
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import child_energy_floor as cef
import child_matter_fraction as cmf

G, c, Msun, mn, lP = cef.G, cef.c, cef.Msun, cef.mn, cef.lP
nb = cef.bounce_density()
fb = cef.alpha * nb
epsb = cef.alpha * nb**2
DUST = 5.79


def child_params(N, eta):
    ab = (3 * N / (4 * math.pi * nb)) ** (1 / 3)
    k = 2 * G / c**4
    return ab, k * N * fb * ab, k * eta * N * mn * c**2        # a_b, r_rad^2, r_m


def adot2(a, ab, rr2, rm):
    return rr2 / a**2 + rm / a - 1 - (rr2 / ab**2) * (ab / a) ** 4


def exit_point(Mp_sun, N, eta):
    ab, rr2, rm = child_params(N, eta)
    rgp = 2 * G * Mp_sun * Msun / c**2
    g = lambda la: math.exp(la) - rgp * math.sqrt(max(1 + adot2(math.exp(la), ab, rr2, rm), 0))
    amax = 0.5 * (rm + math.sqrt(rm * rm + 4 * rr2))
    lo, hi = math.log(ab * 1e3), math.log(amax * (1 - 1e-12))
    if g(lo) > 0:
        return None
    if g(hi) < 0:
        return {"hidden_to_turnaround": True, "a_max_m": amax, "r_gp_m": rgp}
    aex = math.exp(brentq(g, lo, hi, xtol=1e-14))
    f = lambda la: math.exp(la) / math.sqrt(adot2(math.exp(la), ab, rr2, rm)) / c
    t = quad(f, math.log(ab * (1 + 1e-6)), math.log(aex), limit=400, epsrel=1e-10)[0]
    return {"hidden_to_turnaround": False, "a_exit_m": aex, "a_exit_over_a_max": aex / amax, "t_exit_s": t,
            "t_exit_over_GMp_c3": t / (G * Mp_sun * Msun / c**3), "a_max_m": amax, "r_gp_m": rgp}


def main():
    Nn = json.load(open("child_energy_floor.json"))["N_need"]
    Mp = 10.0
    rgp = 2 * G * Mp * Msun / c**2
    out = {"r_gp_m": rgp, "GMp_over_c3_s": G * Mp * Msun / c**3}
    out["SUPERSEDED_neck_at_bounce_note"] = "bare-density value; see israel_matching.json (timelike only from t_exit/27, R = r_g/9)"
    out["neck_at_bounce_m"] = (3 * Mp * Msun * c**2 / (4 * math.pi * epsb)) ** (1 / 3)
    out["neck_at_bounce_over_lP"] = out["neck_at_bounce_m"] / lP
    out["neck_at_bounce_for_N"] = {}
    for label, N in (("N_f_L1", json.load(open("child_energy_floor.json"))["shear_stage2"]["1.0"]["fermions_br"]), ("N_need", Nn)):
        ab, rr2, rm = child_params(N, 0.0)
        Mb = N * fb / c**2
        out["neck_at_bounce_for_N"][label] = ab * (Mp * Msun / Mb) ** (1 / 3)
    out["parent_sphere_radius_coefficient"] = 12 ** (2 / 3) / 4          # r = coeff r_g^(1/3) (c tau)^(2/3), exact Schwarzschild interior
    out["neck_radius_coefficient"] = 2 ** (2 / 3)                         # R_j = coeff r_gp^(1/3) (c t)^(2/3), radiation child
    out["exit"] = {}
    for label, eta in (("radiation", 0.0), ("eta_6e-10_with_hidden_antimatter", DUST * 6e-10), ("eta_1e-6", DUST * 1e-6), ("all_matter", 1.0)):
        out["exit"][label] = exit_point(Mp, Nn, eta)
    # largest N for which a child stays hidden through turnaround (turnaround mass <= M_p)
    for label, eta in (("radiation", 0.0), ("all_matter", 1.0)):
        f = lambda lN: math.log(cmf.turn(math.exp(lN), eta)[1]) - math.log(Mp)
        out[f"N_max_hidden_{label}"] = math.exp(brentq(f, math.log(1e40), math.log(1e70)))
    out["N0_parent_baryons"] = Mp * Msun / mn
    json.dump(out, open("junction_neck.json", "w"), indent=1)
    return out


def figure(out):
    rg = out["r_gp_m"]
    blue, red, grey = "#1f5fa8", "#c0392b", "0.4"
    xbr = json.load(open("kasner_stopping_scales.json"))["species1.0_uc1.0"]["x_br"] * lP
    tau = np.geomspace(xbr, rg / 3, 200)
    tb = (rg / 2) / 27                      # first moment the mass-matched neck is a timelike surface (R = r_g/9)
    tex = rg / 2
    tc = np.geomspace(tb, tex, 200)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 4.3), layout="constrained")
    a2.sharey(None) if False else None
    a1.plot(tau, out["parent_sphere_radius_coefficient"] * rg ** (1 / 3) * tau ** (2 / 3), color=blue, lw=2.2)
    a2.plot(tc, out["neck_radius_coefficient"] * rg ** (1 / 3) * tc ** (2 / 3), color=red, lw=2.2)
    a2.plot([tex], [rg], "o", color=red, mfc="white", mew=1.6, ms=6)
    a1.plot([xbr], [out["parent_sphere_radius_coefficient"] * rg ** (1 / 3) * xbr ** (2 / 3)], "o", color=blue, ms=5)
    a2.plot([tb], [rg / 9], "o", color=red, ms=5)
    for ax in (a1, a2):
        ax.axhline(rg, color=grey, lw=1.0, ls=":")
        ax.set_xscale("log"); ax.set_yscale("log")
    a1.set_xlim(rg / 2, xbr * 0.15); a2.set_xlim(tb * 0.5, rg * 1.2); a2.set_ylim(rg / 20, rg * 2)
    a1.set_ylim(1e-24, 1e6)
    a1.set_xlabel("distance light travels before the crunch (m)")
    a2.set_xlabel("distance light travels after the bounce (m)")
    a2.set_ylabel("neck radius (m)")
    a1.set_ylabel("radius of the sphere (m)")
    a1.text(0.74, 0.74, "dotted line: horizon radius\nof a 10 solar-mass parent\n(29.5 km)", transform=a1.transAxes, color=grey, va="center", ha="center", fontsize=8)
    a1.text(0.3, 0.2, "dot: where the stretch stops,\n0.26 Planck lengths\nbefore the crunch", transform=a1.transAxes, color=blue, va="center", ha="center", fontsize=8)
    a2.text(0.04, 0.58, "filled dot: the neck is a timelike\nsurface only from here on\n(R = r_g/9, 1.8 microseconds);\nearlier it is spacelike", transform=a2.transAxes, color=red, va="center", ha="left", fontsize=8)
    a2.text(0.68, 0.16, "open dot: the neck reaches\nthe horizon radius\nafter 49 microseconds", transform=a2.transAxes, color=red, va="center", ha="center", fontsize=8)
    a1.set_title("Parent: sphere shrinks", loc="left", fontsize=9)
    a2.set_title("Child: neck grows (mass-matched path)", loc="left", fontsize=9)
    fig.savefig("junction_neck.png", dpi=300)


if __name__ == "__main__":
    o = main()
    print(json.dumps(o, indent=1))
    figure(o)
