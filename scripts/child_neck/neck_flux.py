"""What the child gives to the neck layer: energy and quanta swept up as the neck moves through the child's radiation.   Geometric units for lengths (metres).

The neck moves through the child fluid (comoving, radiation) with speed v = |a chi'| relative to comoving observers.  The layer absorbs
    dE   = A (eps + p) v gamma dt   (eps + p = (4/3) eps; energy in metres, G = c = 1),
    dN_f = A n_f v dt               (quanta; A = 4 pi R^2, n_f the child's number density of quanta),
    dB   = eta dN_f                 (net baryon number carried, eta = net matter per quantum).
Child: closed radiation FRW, a^2 = 2 r_rad t - t^2.  Quanta: one Dirac field, g = 3.5, <E> = 3.151 kT, eps = (pi^2/30) g (kT)^4/(hbar c)^3.
Paths: thread (beta_+ = 0), mass-matched, intrinsic P = -sigma/2 layer launched on the thread, and the P = -sigma layer (casimir_layer_check.py).
Writes neck_flux.json and neck_flux.png.  Needs israel_matching.py, casimir_layer_check.py, junction_neck.py and the restored child_* files.
"""
import json
import math

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

import casimir_layer_check as cl
import child_energy_floor as cef
import israel_matching as im

G, c, hbar, eV, Msun, mn = cef.G, cef.c, cef.hbar, 1.602176634e-19, cef.Msun, cef.mn
GEFF = 3.5


def quanta_density(eps_geom):
    """child's number density of quanta (1/m^3) and kT (J) from eps in geometric units (1/m^2)."""
    eps = eps_geom * c**4 / G
    kT = (eps * 30 * (hbar * c) ** 3 / (math.pi**2 * GEFF)) ** 0.25
    return eps / (3.151 * kT), kT


def cumint(f, x):
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])


def radiation_child(t, rr):
    a2 = 2 * rr * t - t**2
    a = np.sqrt(a2)
    ad = (rr - t) / a
    eps_g = 3 * (1 + ad**2) / (8 * np.pi * a2)
    return a, ad, eps_g


def fluxes(R, t, gam, v, rr, rg):
    a, ad, eps_g = radiation_child(t, rr)
    nf = np.array([quanta_density(e)[0] for e in eps_g])
    A = 4 * np.pi * R**2
    dE = A * (4 / 3) * eps_g * v * gam
    dN = A * nf * v
    return cumint(dE, t), cumint(dN, t)


def summarize(label, R, Eabs, Nabs, rg, N0, Elayer=None):
    Mp = rg / 2
    out = dict(label=label, R_end_over_rg=float(R[-1] / rg), E_abs_over_Mp=float(Eabs[-1] / Mp), N_abs=float(Nabs[-1]), N_abs_over_N0=float(Nabs[-1] / N0),
               net_baryon_over_N0={f"eta={eta:g}": float(eta * Nabs[-1] / N0) for eta in (6e-10, 9e-8, 1.0)})
    half = np.argmax(Nabs >= 0.5 * Nabs[-1])
    out["R_where_half_of_quanta_absorbed_over_rg"] = float(R[half] / rg)
    out["R_where_half_of_energy_absorbed_over_rg"] = float(R[np.argmax(Eabs >= 0.5 * Eabs[-1])] / rg)
    if Elayer is not None:
        out["layer_energy_end_over_Mp"] = float(Elayer / Mp)
    out["cumulative"] = {f"R/rg={f:g}": dict(E_abs_over_Mp=float(np.interp(f * rg, R, Eabs) / Mp), N_abs_over_N0=float(np.interp(f * rg, R, Nabs) / N0)) for f in (1e-9, 1e-6, 1e-3, 0.1, 0.5)}
    return out


def main():
    N = json.load(open("child_energy_floor.json"))["N_need"]
    o = im.main()
    rg, rr = o["r_gp_m"], o["child_r_rad_m"]
    Mp = rg / 2
    N0 = 10 * Msun / mn                                                  # 10 solar masses of neutrons
    res = {}
    arr = {}
    # thread
    th = o["_th"]
    R, t = th["R"], th["t"]
    Rdot = 1 / np.tan(th["eta"] / 2)
    a, ad, _ = radiation_child(t, rr)
    F = 1 - R**2 * (1 + ad**2) / a**2
    gam = (Rdot**2 + 1 - R**2 / a**2) / (R * Rdot * ad / a + np.sqrt((1 - R**2 / a**2) * (F + Rdot**2)))
    v = np.abs(th["axd"])
    Ea, Na = fluxes(R, t, gam, v, rr, rg)
    res["thread"] = summarize("thread", R, Ea, Na, rg, N0, Elayer=th["E"][-1])
    arr["thread"] = (R, Ea, Na)
    res["thread"]["T_MeV_at_R_equals_rg"] = float(quanta_density(radiation_child(np.array([t[-1]]), rr)[2][0])[1] / eV / 1e6)
    # mass-matched path
    F_open = im.build(True)
    ab, rr2, rm = im.jn.child_params(N, 0.0)
    aex = im.jn.exit_point(10.0, N, 0.0)["a_exit_m"]
    facs = np.geomspace(0.1935, 1.0, 400)
    args = lambda av: (mp.mpf(av), mp.mpf(rr2), mp.mpf(rm), mp.mpf(ab), mp.mpf(rg))
    Rj = np.array([float(F_open["Rj"](*args(f * aex))) for f in facs])
    gm = np.array([float(F_open["gam"](*args(f * aex))) for f in facs])
    vv = np.array([abs(float(mp.re(F_open["aXd"](*args(f * aex))))) for f in facs])
    aa = facs * aex
    tt = rr - np.sqrt(rr**2 - aa**2)                                   # radiation child: a^2 = 2 r t - t^2
    Ea, Na = fluxes(Rj, tt, gm, vv, rr, rg)
    res["mass_matched"] = summarize("mass-matched (from the timelike onset, R = r_g/9)", Rj, Ea, Na, rg, N0)
    arr["mass_matched"] = (Rj, Ea, Na)
    # layers with intrinsic ratios on the thread: -1/2 and -1
    fn = cl.build_symbols()

    def run(minus_half):
        scale = lambda rho: math.sqrt(rg / rho)

        def rhs(x, y):
            rho = math.exp(x)
            b, t_ = y[0] * scale(rho), math.exp(y[1])
            v_, td, Tun = fn["v"](rho, b, t_, rr, Mp), fn["tdot"](rho, b, t_, rr, Mp), fn["Tun"](rho, b, t_, rr, Mp)
            sg = fn["sigma"](rho, b, t_, rr, Mp)
            sr, sb, st = fn["s_rho"](rho, b, t_, rr, Mp), fn["s_b"](rho, b, t_, rr, Mp), fn["s_t"](rho, b, t_, rr, Mp)
            extra = (sg * v_ / rho) if minus_half else 0.0
            bd = (-Tun - extra - sr * v_ - st * td) / sb
            return [rho * bd / v_ / scale(rho) + 0.5 * y[0], rho * td / (v_ * t_)]

        rho0 = 1e-25
        end = 0.9 if minus_half else 0.99
        sol = solve_ivp(rhs, [math.log(rho0), math.log(rg * end)], [0.0, math.log(rho0 / 2)], method="LSODA", rtol=1e-9, atol=1e-12, dense_output=True, max_step=0.5)
        xs = np.linspace(sol.t[0], sol.t[-1], 600)
        rho = np.exp(xs)
        b = sol.sol(xs)[0] * np.array([scale(r_) for r_ in rho])
        t_ = np.exp(sol.sol(xs)[1])
        g_ = np.array([fn["tdot"](r_, b_, tt_, rr, Mp) for r_, b_, tt_ in zip(rho, b, t_)])
        v_ = np.array([abs(fn["axd"](r_, b_, tt_, rr, Mp)) for r_, b_, tt_ in zip(rho, b, t_)])
        E_ = np.array([4 * np.pi * r_**2 * fn["sigma"](r_, b_, tt_, rr, Mp) for r_, b_, tt_ in zip(rho, b, t_)])
        return rho, t_, g_, v_, E_

    for key, mh in (("layer_minus_half", True), ("layer_minus_one", False)):
        rho, t_, g_, v_, E_ = run(mh)
        Ea, Na = fluxes(rho, t_, g_, v_, rr, rg)
        res[key] = summarize("intrinsic P = -sigma/2 layer" if mh else "intrinsic P = -sigma layer", rho, Ea, Na, rg, N0, Elayer=E_[-1])
        res[key]["child_time_end_over_GMp"] = float(t_[-1] / Mp)
        arr[key] = (rho, Ea, Na)
    res["N0_parent_baryons"] = float(N0)
    res["quanta_energy_note"] = "one Dirac field g=3.5, <E>=3.151 kT"
    json.dump(res, open("neck_flux.json", "w"), indent=1, default=float)
    res["_arr"], res["_rg"], res["_N0"] = arr, rg, N0
    return res


def figure(res, fname="neck_flux.png"):
    arr, rg, N0 = res["_arr"], res["_rg"], res["_N0"]
    Mp = rg / 2
    blue, red, green, grey = "#1f5fa8", "#c0392b", "#2e8b57", "0.4"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0))
    for key, col in (("thread", blue), ("mass_matched", red), ("layer_minus_half", green)):
        R, Ea, Na = arr[key]
        k = (Na / N0 > 1e-4) if key == "mass_matched" else (R > 0)
        a1.loglog(R[k] / rg, Na[k] / N0, color=col, lw=2.5)
        a2.loglog(R[k] / rg, Ea[k] / Mp, color=col, lw=2.5)
    for ax in (a1, a2):
        ax.axhline(1.0, color=grey, ls=":", lw=1.2)
        ax.set_xlabel("neck radius (units of the parent's horizon radius)")
        ax.set_xlim(1e-9, 1.5)
    a1.set_ylim(1e-14, 10)
    a2.set_ylim(1e-5, 30)
    a1.set_ylabel("quanta swept up by the neck (units of the parent's baryon number)")
    a2.set_ylabel("energy swept up by the neck (units of the parent's mass)")
    a1.set_title("a   Quanta", loc="left", fontsize=9)
    a2.set_title("b   Energy", loc="left", fontsize=9)
    a1.text(2e-9, 1.5, "dotted: the parent's baryon number", color=grey, fontsize=8, ha="left", va="bottom")
    a1.text(2e-9, 0.2, "blue: neck follows the thread", color=blue, fontsize=8, ha="left", va="center")
    a1.text(2e-9, 0.05, "green: layer with pressure = -energy/2", color=green, fontsize=8, ha="left", va="center")
    a1.text(2e-9, 0.012, "red: mass-matched neck, counted from R = r_g/9", color=red, fontsize=8, ha="left", va="center")
    a2.text(2e-9, 1.5, "dotted: the parent's mass", color=grey, fontsize=8, ha="left", va="bottom")
    fig.tight_layout()
    fig.savefig(fname, dpi=220)
    return fig


if __name__ == "__main__":
    r = main()
    print(json.dumps({k: v for k, v in r.items() if not k.startswith("_")}, indent=1, default=float))
