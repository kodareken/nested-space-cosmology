"""Normal-stress balance of the composite neck layer, and the layer whose own ratio is -1/2.   Geometric units (metres).

Static balance: T^zz must vanish at every point of the layer (conservation of the planar static stress tensor, exterior vacuum). A constituent with normal pressure p_n is therefore
embedded in a vacuum-like fill (rho_b, p_b = -rho_b, isotropic) with p_n + p_b = 0.
  Dirac vacuum slab (plates, m fields): rho_f = -x/d, p_perp = 3 rho_f, p_par = -rho_f.  Fill rho_b = -3x/d  =>  surface energy -4x, tangential pressure +4x, ratio -1.
  Isotropic torsion fluid (eps_eff, p_eff): fill rho_b = p_eff  =>  surface energy (eps_eff + p_eff) d_m, tangential pressure 0 (negative-density dust where eps_eff + p_eff < 0, i.e. n > 0.544 n_b).
Any balanced mixture has w = -(fraction of the surface energy in the vacuum-like layer).  Thread path: w = -1/2  =>  x = s/8 and (eps_eff + p_eff) d_m = -s/2, s = |E|/(4 pi R^2).
Also: (a) adiabaticity of the normal stress (d/c against the neck proper time tau = (2/3) R^2/|E|);
      (b) a layer with intrinsic P = -sigma/2 integrated from thread initial data stays on the thread (the thread path is its exact solution);
      (c) sizes of an isotropic planar network of negative-tension lines (w = -1/2, T^zz = 0).
Writes balanced_layers.json and balanced_layers.png.  Needs casimir_layer_check.py, torsion_layer_eos.py, israel_matching.py, junction_neck.py and the restored child_* files.
"""
import json
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import casimir_layer_check as cl
import child_energy_floor as cef
import israel_matching as im
import torsion_layer_eos as tl

lP, c, nb, eb, EPSB = cef.lP, cef.c, tl.nb, tl.eb, tl.EPSB
C_DIRAC = 7 * math.pi**2 / 2880


def h_of(f):                                                   # (eps_eff + p_eff)/eps_b at n = f n_b
    e, p = tl.eff(f * nb)
    return (e + p) / eb


def main():
    rg = 2 * cef.G * 10 * cef.Msun / c**2
    o = im.main()
    th, rr = o["_th"], o["child_r_rad_m"]
    Rg, Eneed, _, _ = im.volume_curves(th, rg)
    s = Eneed / (4 * np.pi * Rg**2)
    out = {}
    f0 = brentq(h_of, 0.3, 0.9)
    hp = (h_of(f0 * (1 + 1e-6)) - h_of(f0 * (1 - 1e-6))) / (2e-6 * f0)
    out["torsion_zero_of_eps_plus_p"] = dict(n_over_nb=f0, dh_df=hp, h_at_nb=h_of(1.0), h_at_2nb=h_of(2.0), analytic_f0=(2 / 3) ** 1.5)
    # balanced composite sizing, thread path
    dv = {}
    for m in (1.0, 22.5):
        d = (8 * m * C_DIRAC * lP**2 / s) ** (1 / 3)
        dv[f"m={m:g}"] = dict(R_where_dv_equals_lP_m=float(np.exp(np.interp(0.0, np.log(d / lP), np.log(Rg)))), dv_over_lP_at={f"R={R:g}": float(np.interp(R, Rg, d) / lP) for R in (1e-22, 1e-18, 1.0, 1e3, rg * 0.99)})
    out["balanced_vacuum_layer"] = dv
    # torsion layer: d_m = (s/2)/|h| (EPSB/lP^2 units) ; Planck-thick requires |h| <= s lP/(2 EPSB)
    dh = s * lP / (2 * EPSB)
    out["balanced_torsion_layer"] = {f"R={R:g}": dict(s_per_m=float(np.interp(R, Rg, s)), delta_f_needed_for_dm_equals_lP=float(np.interp(R, Rg, dh) / abs(hp)),
                                                      dm_over_lP_if_n_equals_2nb=float(np.interp(R, Rg, s) / (2 * abs(h_of(2.0)) * EPSB / lP**2) / lP),
                                                      dm_over_lP_if_n_equals_1p1nb=float(np.interp(R, Rg, s) / (2 * abs(h_of(1.1)) * EPSB / lP**2) / lP))
                                     for R in (1e-22, 1e-18, 1.0, 1e3, rg * 0.99)}
    # adiabaticity: relaxation time d/c of the normal stress against the neck proper time tau = (2/3) R^2/|E|
    tau_neck = (2 / 3) * Rg**2 / Eneed
    ratio = lP / tau_neck                                      # d = l_P
    out["adiabatic"] = dict(R_where_lP_over_tau_equals_1_m=float(np.exp(np.interp(0.0, np.log(ratio)[::-1], np.log(Rg)[::-1]))),
                            ratio_at={f"R={R:g}": float(np.interp(R, Rg, ratio)) for R in (1e-22, 1e-18, 1.0, 1e3)})
    # (b) intrinsic P = -sigma/2 layer launched on the thread
    fn = cl.build_symbols()
    Mp = rg / 2
    scale = lambda rho: math.sqrt(rg / rho)

    def rhs(x, y):
        rho = math.exp(x)
        b, t = y[0] * scale(rho), math.exp(y[1])
        v, td, Tun = fn["v"](rho, b, t, rr, Mp), fn["tdot"](rho, b, t, rr, Mp), fn["Tun"](rho, b, t, rr, Mp)
        sg = fn["sigma"](rho, b, t, rr, Mp)
        sr, sb, st = fn["s_rho"](rho, b, t, rr, Mp), fn["s_b"](rho, b, t, rr, Mp), fn["s_t"](rho, b, t, rr, Mp)
        bd = (-Tun - sg * v / rho - sr * v - st * td) / sb         # sigma_dot + (sigma + P) Adot/A = -T_un with P = -sigma/2
        return [rho * bd / v / scale(rho) + 0.5 * y[0], rho * td / (v * t)]

    rho0 = 1e-25
    sol = solve_ivp(rhs, [math.log(rho0), math.log(rg * 0.9)], [0.0, math.log(rho0 / 2)], method="LSODA", rtol=1e-10, atol=1e-13, dense_output=True, max_step=0.5)
    xs = np.linspace(sol.t[0], sol.t[-1], 300)
    rho = np.exp(xs)
    bb = sol.sol(xs)[0] * np.array([scale(r_) for r_ in rho])
    tt = np.exp(sol.sol(xs)[1])
    B = np.array([fn["B"](r_, b_, t_, rr, Mp) for r_, b_, t_ in zip(rho, bb, tt)])
    E = np.array([4 * np.pi * r_**2 * fn["sigma"](r_, b_, t_, rr, Mp) for r_, b_, t_ in zip(rho, bb, tt)])
    Eth = np.interp(rho, th["R"], th["E"])
    out["intrinsic_minus_half_layer_run"] = dict(status=int(sol.status), max_abs_b_over_B=float(np.max(np.abs(bb / B))), max_rel_dev_E_from_thread=float(np.max(np.abs(E / Eth - 1))),
                                                 R_end_over_rg=float(rho[-1] / rg), E_end_over_Mp=float(E[-1] / Mp), child_time_end_over_GMp=float(tt[-1] / Mp))
    # (c) isotropic planar network of negative-tension lines: sigma = -g lambda, P = +g lambda/2, g = G mu / c^4
    out["line_network"] = {f"R={R:g}": dict(min_tension_g_for_spacing_lP=float(np.interp(R, Rg, s) * lP), spacing_m_for_g1=float(1 / np.interp(R, Rg, s)), total_line_length_m_for_g1=float(np.interp(R, Rg, Eneed)))
                           for R in (1e-22, 1e-18, 1e-10, 1.0, 1e3, rg * 0.99)}
    # energy exchange between the two layers if x = s/8 follows s(R): vacuum layer must gain 0.75 s Rdot/R
    out["energy_exchange_per_s_Rdot_over_R"] = dict(vacuum_layer_gain=0.75, dust_like_layer_net=-0.25, total_T_un=-0.5)
    json.dump(out, open("balanced_layers.json", "w"), indent=1, default=float)
    out["_Rg"], out["_s"], out["_run"] = Rg, s, (rho, E, Eth, bb, B)
    return out


def figure(out, fname="balanced_layers.png"):
    Rg, s = out["_Rg"], out["_s"]
    blue, red, grey, green = "#1f5fa8", "#c0392b", "0.4", "#2e8b57"
    rg = Rg[-1]
    hp = out["torsion_zero_of_eps_plus_p"]["dh_df"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0))
    for m, col in ((1.0, blue), (22.5, green)):
        a1.loglog(Rg, (8 * m * C_DIRAC * lP**2 / s) ** (1 / 3) / lP, color=col, lw=2.5)
    a1.axhline(1.0, color=grey, ls=":", lw=1.2)
    a1.set_xlim(1e-30, 6e4)
    a1.set_ylim(1e-6, 1e14)
    a1.set_xlabel("neck radius (m)")
    a1.set_ylabel("thickness of the vacuum layer (Planck lengths)")
    a1.set_title("a   Vacuum layer once the normal stress is balanced", loc="left", fontsize=9)
    a1.text(1e-12, 0.1, "dotted: one Planck length", color=grey, fontsize=8, ha="left", va="top")
    a1.text(1e-6, 3e2, "blue: 1 field", color=blue, fontsize=8, ha="left", va="center")
    a1.text(1e-26, 1e8, "green: 22.5 fields", color=green, fontsize=8, ha="left", va="center")
    a2.loglog(Rg, s * lP / (2 * EPSB) / abs(hp), color=red, lw=2.5)
    a2.set_xlim(1e-30, 6e4)
    a2.set_xlabel("neck radius (m)")
    a2.set_ylabel("how close n must be to 0.544 n_b (fractional)")
    a2.set_title("b   Torsion layer to stay one Planck length thick", loc="left", fontsize=9)
    fig.tight_layout()
    fig.savefig(fname, dpi=220)
    return fig


if __name__ == "__main__":
    o = main()
    print(json.dumps({k: v for k, v in o.items() if not k.startswith("_")}, indent=1, default=float))
