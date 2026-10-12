"""Full Israel (thin-shell) matching of a closed FRW child to a Schwarzschild exterior of mass M_p.   G = c = 1 inside formulas (lengths in metres).

Child: ds^2 = -dt^2 + a(t)^2 [dchi^2 + sin^2(chi) dOmega^2], torsion-modified Friedmann  adot^2 = r_rad^2/a^2 + r_m/a - 1 - (r_rad^2/a_b^2)(a_b/a)^4.
Neck (junction surface): chi = chi(t), areal radius R = a sin(chi), proper time tau.  Exterior: Schwarzschild, f = 1 - 2 M_p / R.
Both sides have K^theta_theta = beta / R with beta = n . grad R, and beta^2 = 1 - 2 m / R + Rdot^2 (m = Misner-Sharp mass, child side m = a sin^3(chi)(1 + adot^2)/2).
Shell content (normal from child side to exterior side, [X] = X_+ - X_-):
    sigma = -[beta] / (4 pi R),     P = ([K^tau_tau] + [K^theta_theta]) / (8 pi),
    K^tau_tau(child) = (a chi'' + 2 adot chi' - a^2 adot chi'^3) / (1 - a^2 chi'^2)^(3/2),   K^tau_tau(exterior) = (Rddot + M_p/R^2) / beta_+ .
Conservation (Codazzi) identity used as the check:  d sigma/d tau + (sigma + P) Adot / A = [T_{mu nu} n^mu u^nu] = - (eps + p) a chi' gamma^2 (child fluid, exterior vacuum).
Topology: 'cap' = the child is the ball chi < chi_j (beta_- > 0).  'open' = the child is the big complement (beta_- < 0), so the exterior sees a neck.
Writes israel_matching.json and israel_matching.png.  Needs the restored child_* scripts, junction_neck.py and their JSON files.
"""
import json
import math

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

import child_energy_floor as cef
import junction_neck as jn

G, c, Msun, lP = cef.G, cef.c, cef.Msun, cef.lP
mp.mp.dps = 40
EPSB = jn.epsb / (1.054571817e-34 * c / lP**4)          # bounce energy density in Planck densities (20.2)


# ----------------------------------------------------------------------------------------- symbolic checks
def symbolic_checks():
    t = sp.symbols("t", real=True)
    a = sp.Function("a")(t)
    x = sp.Function("chi")(t)
    ad, xd = sp.diff(a, t), sp.diff(x, t)
    gam = 1 / sp.sqrt(1 - a**2 * xd**2)
    n_lo = [-a * xd * gam, a * gam]
    n_up = [a * xd * gam, gam / a]
    u = [gam, gam * xd]
    acc_t = gam * sp.diff(u[0], t) + a * ad * u[1] ** 2
    acc_x = gam * sp.diff(u[1], t) + 2 * (ad / a) * u[0] * u[1]
    Ktt = sp.simplify(n_lo[0] * acc_t + n_lo[1] * acc_x)
    R = a * sp.sin(x)
    beta = n_up[0] * ad * sp.sin(x) + n_up[1] * a * sp.cos(x)
    Rdot = gam * sp.diff(R, t)
    ident = sp.simplify((beta**2 - Rdot**2) - (sp.cos(x) ** 2 - ad**2 * sp.sin(x) ** 2))
    return {"beta2_minus_Rdot2_minus_(1-2m/R)": str(ident), "K_tau_tau_child": str(Ktt)}


def conservation_check():
    """Radiation FRW a = sqrt(1-(t-1)^2), arbitrary timelike neck chi = 1.1 + 0.2 t, arbitrary M = 0.37: Codazzi identity to ~1e-15."""
    t = sp.symbols("t", real=True)
    a = sp.sqrt(1 - (t - 1) ** 2)
    x = sp.Rational(11, 10) + sp.Rational(1, 5) * t
    M = sp.Rational(37, 100)
    ad, xd = sp.diff(a, t), sp.diff(x, t)
    gam = 1 / sp.sqrt(1 - a**2 * xd**2)
    R = a * sp.sin(x)
    D = lambda f: gam * sp.diff(f, t)
    Rd = D(R)
    Rdd = D(Rd)
    beta_m = gam * (a * ad * xd * sp.sin(x) + sp.cos(x))
    Ktt_m = (a * sp.diff(x, t, 2) + 2 * ad * xd - a**2 * ad * xd**3) / (1 - a**2 * xd**2) ** sp.Rational(3, 2)
    beta_p = sp.sqrt(1 - 2 * M / R + Rd**2)
    Ktt_p = (Rdd + M / R**2) / beta_p
    sigma = -(beta_p - beta_m) / (4 * sp.pi * R)
    P = ((Ktt_p - Ktt_m) + (beta_p - beta_m) / R) / (8 * sp.pi)
    A = 4 * sp.pi * R**2
    eps = 3 / (8 * sp.pi * a**4)
    lhs = D(sigma) + (sigma + P) * D(A) / A
    rhs = -(eps + eps / 3) * a * xd * gam**2
    errs = []
    for tv in (sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(9, 10), sp.Rational(13, 10)):
        l, r = sp.N(lhs.subs(t, tv), 25), sp.N(rhs.subs(t, tv), 25)
        errs.append(float(abs(l - r) / abs(r)))
    return {"max_rel_error_codazzi": max(errs)}


def os_dust_check():
    """Oppenheimer-Snyder cap: closed dust a_m, comoving neck chi0, M = a_m sin^3(chi0)/2  ->  sigma = P = 0 exactly."""
    out = []
    for chi0 in (0.3, 0.9, 1.3):
        for eta in (0.7, 1.9, 3.0):
            am = 2.0
            a = am / 2 * (1 - math.cos(eta))
            ad2 = am / a - 1
            ad = math.sqrt(ad2)
            add = -am / (2 * a * a)
            R = a * math.sin(chi0)
            M = am * math.sin(chi0) ** 3 / 2
            beta_m = math.cos(chi0)
            beta_p = math.sqrt(1 - 2 * M / R + (ad * math.sin(chi0)) ** 2)
            Ktt_p = (add * math.sin(chi0) + M / R**2) / beta_p
            out.append((beta_p - beta_m, Ktt_p))
    return {"max_abs_beta_jump": max(abs(o[0]) for o in out), "max_abs_K_tau_tau_exterior": max(abs(o[1]) for o in out)}


# ----------------------------------------------------------------------------------- mass-matched neck path
A_, RR2, RM, AB, RG = sp.symbols("A rr2 rm ab rg", positive=True)
gA = RR2 / A_**2 + RM / A_ - 1 - RR2 * AB**2 / A_**4
RjA = (RG * A_**2 / (1 + gA)) ** sp.Rational(1, 3)


def build(open_topology=True):
    s = RjA / A_
    chi = sp.pi - sp.asin(s) if open_topology else sp.asin(s)
    ad = sp.sqrt(gA)
    D = lambda f: ad * sp.diff(f, A_)
    xd = D(chi)
    xdd = D(xd)
    gam = 1 / sp.sqrt(1 - A_**2 * xd**2)
    Rd = gam * D(RjA)
    Rdd = gam * D(Rd)
    beta_m = gam * (A_ * ad * xd * sp.sin(chi) + sp.cos(chi))
    Ktt_m = (A_ * xdd + 2 * ad * xd - A_**2 * ad * xd**3) / (1 - A_**2 * xd**2) ** sp.Rational(3, 2)
    Mp = RG / 2
    beta_p = -beta_m if open_topology else beta_m
    Ktt_p = (Rdd + Mp / RjA**2) / beta_p
    sigma = -(beta_p - beta_m) / (4 * sp.pi * RjA)
    P = ((Ktt_p - Ktt_m) + (beta_p - beta_m) / RjA) / (8 * sp.pi)
    ex = dict(sigma=sigma, P=P, beta_m=beta_m, beta_p=beta_p, gam=gam, aXd=A_ * xd, Rj=RjA)
    return {k: sp.lambdify((A_, RR2, RM, AB, RG), v, "mpmath") for k, v in ex.items()}


def ms_onset(F_open, N, eta, Mp_sun=10.0):
    ab, rr2, rm = jn.child_params(N, eta)
    rg = 2 * G * Mp_sun * Msun / c**2
    ar = lambda av: (mp.mpf(av), mp.mpf(rr2), mp.mpf(rm), mp.mpf(ab), mp.mpf(rg))
    ex = jn.exit_point(Mp_sun, N, eta)
    aex = ex["a_exit_m"]

    def h(la):
        try:
            v = F_open["aXd"](*ar(math.exp(la)))
            return float(abs(mp.re(v))) - 1.0
        except Exception:
            return 1.0

    prev = None
    at = None
    for la in np.linspace(math.log(aex), math.log(aex) - 12, 600):
        if h(la) > 0:
            at = math.exp(brentq(h, la, prev, xtol=1e-14))
            break
        prev = la
    g = lambda a_: rr2 / a_**2 + rm / a_ - 1
    tt = lambda a_: quad(lambda la: math.exp(la) / math.sqrt(g(math.exp(la))), math.log(ab * (1 + 1e-6)), math.log(a_), limit=400, epsrel=1e-10)[0]
    Rj = lambda a_: (rg * a_**2 / (1 + g(a_))) ** (1 / 3)
    return dict(a_exit_m=aex, a_timelike_m=at, a_ratio=at / aex, t_exit_over_t_timelike=(ex["t_exit_s"] * c) / tt(at), t_ratio=tt(at) / (ex["t_exit_s"] * c), R_timelike_over_rg=Rj(at) / rg)


def ms_table(F_open, F_cap, N, eta, facs, Mp_sun=10.0):
    ab, rr2, rm = jn.child_params(N, eta)
    rg = 2 * G * Mp_sun * Msun / c**2
    aex = jn.exit_point(Mp_sun, N, eta)["a_exit_m"]
    ar = lambda av: (mp.mpf(av), mp.mpf(rr2), mp.mpf(rm), mp.mpf(ab), mp.mpf(rg))
    rows = []
    for fac in facs:
        vo = {k: F_open[k](*ar(fac * aex)) for k in ("aXd", "Rj", "beta_p", "sigma", "P")}
        vc = {k: F_cap[k](*ar(fac * aex)) for k in ("sigma", "P")}
        A4 = 4 * mp.pi * vo["Rj"] ** 2
        Mp = rg / 2
        rows.append(dict(a_over_a_exit=fac, R_over_rg=float(vo["Rj"] / rg), abs_a_chi_dot=float(abs(vo["aXd"])), beta=float(vo["beta_p"]),
                         open_E_over_Mp=float(A4 * vo["sigma"] / Mp), open_PA_over_Mp=float(A4 * vo["P"] / Mp), open_w=float(vo["P"] / vo["sigma"]),
                         cap_sigma=float(vc["sigma"]), cap_PA_over_Mp=float(A4 * vc["P"] / Mp)))
    return rows


def early_neck_correction(N, Mp_sun=10.0):
    """Earlier note: neck at the bounce 3.57e-23 m used the bare energy density.  With the geometric mass: R_j^3 = a_b^3 (M_p/M_b) x^6/(x^2-1), x = a/a_b."""
    ab, rr2, rm = jn.child_params(N, 0.0)
    Mb = N * jn.fb / c**2
    eps_m = Mp_sun * Msun / Mb
    xs = math.sqrt(1.5)
    Rmin = ab * (eps_m * xs**6 / (xs**2 - 1)) ** (1 / 3)
    out = dict(a_b_m=ab, M_p_over_M_b=eps_m, x_of_min=xs, R_min_m=Rmin, R_min_over_old=Rmin / (ab * eps_m ** (1 / 3)), old_value_m=ab * eps_m ** (1 / 3))
    F_open = build(True)
    rg = 2 * G * Mp_sun * Msun / c**2
    pts = {}
    for u in (1e-39, 1e-30, 1e-20, 1e-10, 1e-3, 1.0):
        av = mp.mpf(ab) * mp.sqrt(1 + mp.mpf(u))
        v = F_open["aXd"](av, mp.mpf(rr2), mp.mpf(rm), mp.mpf(ab), mp.mpf(rg))
        pts[f"x2m1={u:g}"] = float(abs(mp.re(v)))
    out["abs_a_chi_dot_near_bounce"] = pts
    return out


# ------------------------------------------------------------------------------- neck that follows the parent's thread
def thread_neck(rg, rr, eta0=3.7e-15):
    """Exterior side follows the radial geodesic T = const of the white-hole region (beta_+ = 0):  R = (r_g/2)(1-cos eta), tau = (r_g/2)(eta - sin eta), Rdot = cot(eta/2).
    Child is radiation, a^2 = 2 r_rad t - t^2.  Child time t(eta) from the timelike condition of the neck (quadratic for dt/dtau)."""
    def rho_of(eta):
        return 0.5 * rg * (1 - math.cos(eta)) if eta > 1e-3 else 0.25 * rg * eta**2 * (1 - eta**2 / 12)

    def rhs(s, y):
        eta, t_ = math.exp(s), math.exp(y[0])
        rho, rhod = rho_of(eta), 1 / math.tan(eta / 2)
        a2 = 2 * rr * t_ - t_**2
        a_ = math.sqrt(a2)
        ad_ = (rr - t_) / a_
        X = rho * rhod * ad_ / a_
        F = 1 - rho**2 * (1 + ad_**2) / a2
        Dd = (1 - rho**2 / a2) * (F + rhod**2)
        tdot = (rhod**2 + 1 - rho**2 / a2) / (X + math.sqrt(max(Dd, 0)))
        return [eta * rho * tdot / t_]

    rho0 = 0.125 * rg * eta0**2
    sol = solve_ivp(rhs, [math.log(eta0), math.log(math.pi * (1 - 1e-9))], [math.log(rho0)], method="LSODA", rtol=1e-11, atol=1e-13, dense_output=True)
    ss = np.linspace(sol.t[0], sol.t[-1], 4000)
    eta = np.exp(ss)
    tt = np.exp(sol.sol(ss)[0])
    rho = np.where(eta > 1e-3, 0.5 * rg * (1 - np.cos(eta)), 0.25 * rg * eta**2 * (1 - eta**2 / 12))
    rhod = 1 / np.tan(eta / 2)
    a2 = 2 * rr * tt - tt**2
    aa = np.sqrt(a2)
    ada = (rr - tt) / aa
    F = 1 - rho**2 * (1 + ada**2) / a2
    mminus = rho * (1 - F) / 2
    E = -rho * np.sqrt(F + rhod**2)
    sig = E / (4 * np.pi * rho**2)
    sigdot = sig * (np.gradient(np.log(np.abs(sig)), ss) / eta) / rho
    X = rho * rhod * ada / aa
    tdot = (rhod**2 + 1 - rho**2 / a2) / (X + np.sqrt((1 - rho**2 / a2) * (F + rhod**2)))
    axd = -(rhod / tdot - rho * ada / aa) / np.sqrt(1 - rho**2 / a2)
    Tun = (1 + ada**2) / (2 * np.pi * a2) * axd * tdot**2
    P = (-Tun - sigdot) / (2 * rhod / rho) - sig
    return dict(eta=eta, t=tt, R=rho, E=E, P=P, sig=sig, mminus=mminus, axd=axd, tau=np.where(eta > 1e-2, 0.5 * rg * (eta - np.sin(eta)), rg / 12 * eta**3 * (1 - eta**2 / 20)))


def volume_curves(th, rg):
    """Negative energy needed on the thread path (sqrt(R r_g) below the integrated range, where m_- << M_p) versus the energy the neck volume holds at the bounce density
    (geometric units: G eps / c^4 = EPSB / lP^2 per m^2)."""
    Rg = np.geomspace(1e-32, rg, 600)
    Eneed = np.where(Rg < th["R"][0], np.sqrt(Rg * rg), np.interp(Rg, th["R"], -th["E"]))
    Evol = (4 * np.pi / 3) * Rg**3 * (EPSB / lP**2)
    return Rg, Eneed, Evol, Evol / Eneed


# --------------------------------------------------------------------------------------------------- main + figure
def main():
    Nn = json.load(open("child_energy_floor.json"))["N_need"]
    Mp_sun = 10.0
    rg = 2 * G * Mp_sun * Msun / c**2
    Mp = rg / 2
    ab, rr2, rm = jn.child_params(Nn, 0.0)
    rr = math.sqrt(rr2)
    out = {"r_gp_m": rg, "GMp_over_c3_s": G * Mp_sun * Msun / c**3, "N_need": Nn, "child_r_rad_m": rr}
    out["symbolic"] = symbolic_checks()
    out["codazzi_identity"] = conservation_check()
    out["oppenheimer_snyder_cap"] = os_dust_check()
    F_open, F_cap = build(True), build(False)
    out["timelike_onset"] = {lab: ms_onset(F_open, Nn, eta) for lab, eta in (("radiation", 0.0), ("eta_5.79e-6", 5.79e-6), ("all_matter", 1.0))}
    out["timelike_onset"]["analytic_radiation"] = dict(a_ratio=3 ** -1.5, t_ratio=3 ** -3, R_over_rg=3 ** -2)
    out["ms_path_radiation"] = ms_table(F_open, F_cap, Nn, 0.0, (0.2, 0.25, 0.3, 0.5, 0.75, 1.0, 1.5, 3.0, 10.0))
    out["early_neck_correction"] = early_neck_correction(Nn)
    th = thread_neck(rg, rr)
    sel = (500, 1000, 2000, 3000, 3500, 3900)
    out["thread_path"] = [dict(R_over_rg=float(th["R"][i] / rg), R_m=float(th["R"][i]), E_over_Mp=float(th["E"][i] / Mp), PA_over_Mp=float(4 * np.pi * th["R"][i] ** 2 * th["P"][i] / Mp),
                               w=float(th["P"][i] / th["sig"][i]), mminus_over_Mp=float(th["mminus"][i] / Mp), child_time_over_GMp=float(th["t"][i] / Mp), abs_a_chi_dot=float(abs(th["axd"][i]))) for i in sel]
    out["thread_end"] = dict(R_over_rg=float(th["R"][-1] / rg), E_over_Mp=float(th["E"][-1] / Mp), mminus_over_Mp=float(th["mminus"][-1] / Mp),
                             child_time_over_GMp=float(th["t"][-1] / Mp), parent_tau_over_GMp=float(th["tau"][-1] / Mp), Rdot_at_end=float(1 / np.tan(th["eta"][-1] / 2)))
    out["thread_E_formula_check"] = dict(max_rel_dev_from_sqrt_2R_Mp_minus_m=float(np.max(np.abs(th["E"] + np.sqrt(2 * th["R"] * (Mp - th["mminus"]))) / np.abs(th["E"]))))
    # volume check: energy held by the neck volume at the bounce density versus the negative energy needed (thread path, geometric units: G eps/c^4 = 20.2/lP^2 per m^2)
    Rg, Eneed, Evol, ratio = volume_curves(th, rg)
    Rstar_an = ((3 / (4 * math.pi)) * math.sqrt(rg) * lP**2 / EPSB) ** 0.4          # E_need = sqrt(R r_g) (m_- << M_p)
    Rstar = float(np.exp(np.interp(0.0, np.log(ratio), np.log(Rg))))
    out["volume_check"] = dict(R_star_m=Rstar, R_star_analytic_m=Rstar_an, R_star_over_lP=Rstar / lP, eps_b_over_eps_P=EPSB,
                               layer_thickness_at_r_g_m=float(Eneed[-1] * lP**2 / (EPSB * 4 * math.pi * rg**2)),
                               layer_thickness_at_1m_m=float(math.sqrt(1.0 * rg) * lP**2 / (EPSB * 4 * math.pi)),
                               thread_start_R_m=float(th["R"][0]))
    out["energy_unit_J_per_m"] = c**4 / G
    json.dump(out, open("israel_matching.json", "w"), indent=1, default=float)
    out["_th"] = th
    out["_ms_dense"] = ms_table(F_open, F_cap, Nn, 0.0, list(np.geomspace(0.1935, 12, 90)))
    return out


def figure(out, fname="israel_matching.png"):
    th = out["_th"]
    rg = out["r_gp_m"]
    Mp = rg / 2
    blue, red, grey = "#1f5fa8", "#c0392b", "0.4"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.2))
    Rg0, Eneed0, _, _ = volume_curves(th, rg)
    a1.loglog(Rg0 / rg, Eneed0 / Mp, color=blue, lw=2.5)
    rows = out["_ms_dense"]
    a1.loglog([r["R_over_rg"] for r in rows], [-r["open_E_over_Mp"] for r in rows], color=red, lw=2.5)
    a1.axhline(1.0, color=grey, ls=":", lw=1.2)
    a1.set_xlim(1e-3, 40)
    a1.set_ylim(1e-2, 400)
    a1.set_xlabel("neck radius (units of the parent's horizon radius)")
    a1.set_ylabel("negative energy at the neck (units of the parent's mass)")
    a1.set_title("a   Negative energy needed at the neck", loc="left", fontsize=9)
    a1.text(1.3e-3, 1.45, "dotted line: the parent's own mass", color=grey, fontsize=8, ha="left", va="bottom")
    a1.text(4e-3, 60, "red: neck holds exactly the parent's mass;\nit exists only from R = r_g/9 on", color=red, fontsize=8, ha="left", va="center")
    a1.text(36, 0.05, "blue: neck follows the parent's thread,\nthe smallest demand of any path;\nit ends at R = r_g", color=blue, fontsize=8, ha="right", va="center")
    Rg, Eneed, Evol, ratio = volume_curves(th, rg)
    a2.loglog(Rg, ratio, color=blue, lw=2.5)
    a2.axhline(1.0, color=grey, ls=":", lw=1.2)
    Rs = out["volume_check"]["R_star_m"]
    a2.plot([Rs], [1.0], "o", color=blue, ms=7)
    a2.set_xlabel("neck radius (m)")
    a2.set_ylabel("energy in the neck volume at bounce density / energy needed")
    a2.set_title("b   Does the neck volume hold enough?", loc="left", fontsize=9)
    a2.text(2.5e-26, 1e-3, "dot: %.0e m (%.0e Planck lengths),\nabove it the volume holds enough" % (Rs, Rs / lP), color=blue, fontsize=8, ha="left", va="center")
    fig.tight_layout()
    fig.savefig(fname, dpi=220)
    return fig


if __name__ == "__main__":
    o = main()
    print(json.dumps({k: v for k, v in o.items() if k != "_th"}, indent=1, default=float)[:6000])
