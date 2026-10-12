"""Dirac vacuum (Casimir-type) layer at the neck, and the path-independent relation between the layer's pressure/energy ratio and its energy growth.   G = c = 1 (metres).

(1) Casimir coefficient of a massless Dirac field between plates (bag boundary conditions): E/A = -(7 pi^2 / 2880) hbar c / d^3, from zeta regularisation of the half-integer spectrum, 4 modes x (-1/2) omega.
    Boost invariance along the plates gives p_parallel = -rho, i.e. surface pressure P = -sigma (w = -1) for any planar vacuum layer.
(2) Exact relation along any outgoing neck (child fluid at rest, exterior vacuum), from the Codazzi identity  d sigma/d tau + (sigma + P) Adot/A = -(eps + p) a chi' gamma^2 :
        sigma (q + 2 w) Rdot / R = -T_un,   q = d ln|E| / d ln R,  E = 4 pi R^2 sigma,  T_un = (eps + p) a chi' gamma^2 <= 0 for a chi' <= 0   ==>   w <= -q/2.
(3) Dirac layer thickness on the thread path that carries the required |E| (sqrt(2R(M_p - m_-))): d = [ (7 pi^2/2880) l_P^2 4 pi R^2 / |E| ]^(1/3).
(4) A w = -1 layer launched on the thread (b = beta_+ = 0, R0 = 1e-25 m): integrate sigma_dot = -T_un for the radiation child (state R, b, t).  Shows E growth, Rdot = 2 pi |sigma| R and E tau ~ R^2.
Writes casimir_layer_check.json and casimir_layer_check.png.  Needs israel_matching.py, junction_neck.py and the restored child_* files.
"""
import json
import math

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

import child_energy_floor as cef
import israel_matching as im

G, c, hbar, lP = cef.G, cef.c, cef.hbar, cef.lP
C_DIRAC = 7 * math.pi**2 / 2880


def casimir_coefficients():
    """Zeta regularisation: modes k_n = (n + s) pi/d, E/A = K zeta_H(-3, s) with K = -pi^2/12 for one +1/2 omega oscillator (scalar Dirichlet check: -pi^2/1440)."""
    K = -mp.pi**2 / 12
    scalar = K * mp.zeta(-3)
    half = K * mp.zeta(-3, 0.5)                                # half-integer modes, boson sign
    dirac = -4 * half                                          # 4 fermionic oscillators (2 spins, particle and antiparticle): opposite sign
    return {"scalar_dirichlet": float(scalar), "scalar_dirichlet_expected": -math.pi**2 / 1440, "dirac_bag": float(dirac), "dirac_bag_expected": -C_DIRAC}


def thread_layer_sizing(th, rg):
    Rg, Eneed, _, _ = im.volume_curves(th, rg)
    d = ((C_DIRAC * lP**2 * 4 * np.pi * Rg**2 / Eneed)) ** (1 / 3)
    rho_dot = np.sqrt(np.maximum(rg / Rg - 1, 0))
    adiab = d * rho_dot / (2 * Rg)                              # (d/c) d ln d / d tau with d ~ sqrt(R)
    out = {}
    for R in (1e-22, 1e-18, 1e-10, 1.0, 1e3, rg * 0.99):
        out[f"R={R:g}"] = dict(d_m=float(np.interp(R, Rg, d)), d_over_R=float(np.interp(R, Rg, d) / R), adiabatic_parameter=float(np.interp(R, Rg, adiab)))
    out["R_where_d_equals_lP_m"] = float(np.exp(np.interp(0.0, np.log(d / lP), np.log(Rg))))
    out["R_where_adiabatic_param_equals_1_m"] = float(Rg[np.argmax(adiab < 1)])
    out["d_coefficient_m_per_sqrtm"] = float(d[0] / math.sqrt(Rg[0]))
    return out


def build_symbols():
    rho, b, t = sp.symbols("rho b t", positive=True)
    r, M = sp.symbols("r_rad M_p", positive=True)
    a = sp.sqrt(2 * r * t - t**2)
    ad = (r - t) / a
    mm = rho**3 * (1 + ad**2) / (2 * a**2)
    B = sp.sqrt(b**2 + 2 * (M - mm) / rho)
    sigma = -(b + B) / (4 * sp.pi * rho)
    v = sp.sqrt(b**2 + 2 * M / rho - 1)
    ct = sp.sqrt(1 - rho**2 / a**2)
    tdot = (v**2 + ct**2) / (rho * v * ad / a + ct * B)
    axd = -(v / tdot - rho * ad / a) / ct
    Tun = (1 + ad**2) / (2 * sp.pi * a**2) * axd * tdot**2
    args = (rho, b, t, r, M)
    ex = dict(sigma=sigma, B=B, v=v, tdot=tdot, axd=axd, Tun=Tun, s_rho=sp.diff(sigma, rho), s_b=sp.diff(sigma, b), s_t=sp.diff(sigma, t), mm=mm)
    return {k: sp.lambdify(args, e, "math") for k, e in ex.items()}


def w_minus_one_run(fn, rg, rr, rho0=1e-25):
    Mp = rg / 2
    scale = lambda rho: math.sqrt(rg / rho)

    def rhs(x, y):
        rho = math.exp(x)
        b, t = y[0] * scale(rho), math.exp(y[1])
        v, td, Tun = fn["v"](rho, b, t, rr, Mp), fn["tdot"](rho, b, t, rr, Mp), fn["Tun"](rho, b, t, rr, Mp)
        sr, sb, st = fn["s_rho"](rho, b, t, rr, Mp), fn["s_b"](rho, b, t, rr, Mp), fn["s_t"](rho, b, t, rr, Mp)
        bd = (-Tun - sr * v - st * td) / sb                      # sigma_dot = -T_un   (w = -1)
        return [rho * bd / v / scale(rho) + 0.5 * y[0], rho * td / (v * t)]

    sol = solve_ivp(rhs, [math.log(rho0), math.log(rg * 0.99)], [0.0, math.log(rho0 / 2)], method="LSODA", rtol=1e-9, atol=1e-12, dense_output=True, max_step=0.5)
    xs = np.linspace(sol.t[0], sol.t[-1], 400)
    rho = np.exp(xs)
    b = sol.sol(xs)[0] * np.array([scale(r_) for r_ in rho])
    t = np.exp(sol.sol(xs)[1])
    sg = np.array([fn["sigma"](r_, b_, t_, rr, Mp) for r_, b_, t_ in zip(rho, b, t)])
    v = np.array([fn["v"](r_, b_, t_, rr, Mp) for r_, b_, t_ in zip(rho, b, t)])
    B = np.array([fn["B"](r_, b_, t_, rr, Mp) for r_, b_, t_ in zip(rho, b, t)])
    E = 4 * np.pi * rho**2 * sg
    tau = np.concatenate([[rho[0] / v[0]], rho[0] / v[0] + np.cumsum(0.5 * (1 / v[1:] + 1 / v[:-1]) * np.diff(rho))])
    return dict(rho=rho, b=b, t=t, sigma=sg, v=v, B=B, E=E, tau=tau, status=int(sol.status))


def main():
    N = json.load(open("child_energy_floor.json"))["N_need"]
    o = im.main()
    th, rg, rr = o["_th"], o["r_gp_m"], o["child_r_rad_m"]
    Mp = rg / 2
    out = {"dirac_casimir_coefficient": C_DIRAC, "dirac_casimir_coefficient_check": casimir_coefficients()}
    out["dirac_layer_on_thread_path"] = thread_layer_sizing(th, rg)
    # identity check on the thread path: sigma (q + 2 w) Rdot/R = -T_un
    R, E, P, sig = th["R"], th["E"], th["P"], th["sig"]
    q = np.gradient(np.log(np.abs(E)), np.log(R))
    w = P / sig
    Rdot = 1 / np.tan(th["eta"] / 2)
    t_ = th["t"]
    a2 = 2 * rr * t_ - t_**2
    aa = np.sqrt(a2)
    ada = (rr - t_) / aa
    F = 1 - R**2 * (1 + ada**2) / a2
    tdot = (Rdot**2 + 1 - R**2 / a2) / (R * Rdot * ada / aa + np.sqrt((1 - R**2 / a2) * (F + Rdot**2)))
    Tun = (1 + ada**2) / (2 * np.pi * a2) * th["axd"] * tdot**2
    lhs = sig * (q + 2 * w) * Rdot / R
    sel = (R > 1e-20 * rg) & (R < 0.5 * rg)
    ok = (R > 1e-20 * rg) & (R < 0.5 * rg)
    out["thread_identity_rel_dev_median"] = float(np.median(np.abs(lhs[ok] + Tun[ok]) / np.abs(Tun[ok])))
    out["thread_T_un_negative_everywhere"] = bool(np.all(Tun[ok] < 0))
    out["thread_path_q_and_w"] = dict(q_median=float(np.median(q[sel])), w_median=float(np.median(w[sel])), w_over_minus_q_half=float(np.median(w[sel] / (-q[sel] / 2))))
    # (4) w = -1 layer launched on the thread
    fn = build_symbols()
    run = w_minus_one_run(fn, rg, rr)
    rho, E_run, sg, v, tau = run["rho"], run["E"], run["sigma"], run["v"], run["tau"]
    qrun = np.gradient(np.log(np.abs(E_run)), np.log(rho))
    out["w_minus_one_run"] = dict(status=run["status"], R_end_m=float(rho[-1]), E_end_over_Mp=float(E_run[-1] / Mp), q_late=float(np.median(qrun[200:380])),
                                  b_over_B_end=float(run["b"][-1] / run["B"][-1]), v_end=float(v[-1]), tau_total_s=float(tau[-1] / c), tau_total_over_torsion_time=float(tau[-1] / c / (lP / (c * math.sqrt(8 * math.pi / 3 * im.EPSB)))),
                                  Rdot_over_2pi_sigma_R_end=float(v[-1] / (2 * np.pi * abs(sg[-1]) * rho[-1])), E_tau_over_R2_end=float(abs(E_run[-1]) * tau[-1] / rho[-1] ** 2),
                                  sigma_start=float(sg[0]), sigma_end=float(sg[-1]))
    sP = C_DIRAC / lP                                          # |sigma| of a Dirac layer one Planck length thick
    iP = int(np.argmax(np.abs(sg) <= sP))
    dth = lambda s_: (C_DIRAC * lP**2 / abs(s_)) ** (1 / 3)
    out["w_minus_one_run"].update(dict(R_where_vacuum_layer_thickness_equals_lP_m=float(rho[iP]), vacuum_layer_thickness_start_over_lP=float(dth(sg[0]) / lP),
                                       vacuum_layer_thickness_end_over_lP=float(dth(sg[-1]) / lP), child_time_at_end_m=float(run["t"][-1]), R_over_child_time_at_end=float(rho[-1] / run["t"][-1])))
    out["thread_E_tau_over_R2"] = {f"R/rg={f:g}": float(abs(np.interp(f * rg, R, E)) * np.interp(f * rg, R, th["tau"]) / (f * rg) ** 2) for f in (1e-12, 1e-6, 1e-2, 0.3)}
    json.dump(out, open("casimir_layer_check.json", "w"), indent=1, default=float)
    out["_run"], out["_th"] = run, th
    out["_qw"] = (q, w)
    return out


def figure(out, fname="casimir_layer_check.png"):
    th, run = out["_th"], out["_run"]
    rg = 2 * G * 10 * cef.Msun / c**2
    Mp = rg / 2
    blue, red, grey, green = "#1f5fa8", "#c0392b", "0.4", "#2e8b57"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0))
    Rg, Eneed, _, _ = im.volume_curves(th, rg)
    a1.loglog(Rg / rg, Eneed / Mp, color=blue, lw=2.5)
    a1.loglog(run["rho"] / rg, -run["E"] / Mp, color=red, lw=2.5)
    a1.axhline(1.0, color=grey, ls=":", lw=1.2)
    a1.set_xlim(1e-30, 3)
    a1.set_ylim(1e-18, 1e44)
    a1.set_xlabel("neck radius (units of the parent's horizon radius)")
    a1.set_ylabel("negative energy at the neck (units of the parent's mass)")
    a1.set_title("a   Thread path against a layer with pressure = -energy", loc="left", fontsize=9)
    a1.text(1e-17, 2.5, "dotted: the parent's mass", color=grey, fontsize=8, ha="left", va="bottom")
    a1.text(2.0, 3e-14, "blue: neck follows the parent's thread", color=blue, fontsize=8, ha="right", va="center")
    a1.text(1e-27, 1e30, "red: layer with pressure = -energy,\nlaunched on the thread; the neck\nexpands at near light speed", color=red, fontsize=8, ha="left", va="center")
    qq = np.linspace(-6, 4, 100)
    a2.fill_between(qq, -q_half(qq), -8, color="0.88")
    a2.plot(qq, -qq / 2, color=grey, lw=1.2)
    a2.axhline(-1.0, color=green, lw=2.0, xmin=(0.5 + 6) / 10, xmax=(2.0 + 6) / 10)
    a2.plot([0.5], [-0.5], "o", color=blue, ms=7)
    qt = np.linspace(-6, -2, 50)
    a2.fill_between(qt, 1.0, -qt / 2, color=red, alpha=0.35)
    a2.set_xlim(-6, 4)
    a2.set_ylim(-3, 5)
    a2.set_xlabel("how fast the neck's negative energy grows: d ln|E| / d ln R")
    a2.set_ylabel("pressure / energy density of the layer")
    a2.set_title("b   What each layer type allows", loc="left", fontsize=9)
    a2.text(-3.0, -2.2, "grey: allowed region,\nratio below -(growth)/2", color=grey, fontsize=8, ha="center", va="center")
    a2.text(-3.6, 4.2, "red: torsion matter (ratio of at least 1)\nfits only where the energy\nfalls faster than R^-2", color=red, fontsize=8, ha="center", va="center")
    a2.text(0.8, -2.0, "green: vacuum layer,\nratio -1, growth\n1/2 to 2", color=green, fontsize=8, ha="center", va="center")
    a2.text(0.9, -0.25, "blue dot: thread path", color=blue, fontsize=8, ha="left", va="center")
    fig.tight_layout()
    fig.savefig(fname, dpi=220)
    return fig


def q_half(q):
    return q / 2


if __name__ == "__main__":
    o = main()
    print(json.dumps({k: v for k, v in o.items() if not k.startswith("_")}, indent=1, default=float))
