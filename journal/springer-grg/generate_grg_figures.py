#!/usr/bin/env python3
"""Generate 600-DPI and vector PDF figures for the Springer GRG manuscript."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_offset_generator import (
    nested_zero_crossing_map,
    nonlocal_to_local_atlas_readout,
    unified_least_resistance_cycle,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = Path(__file__).resolve().parent / "figures"


def save_both(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIG_DIR / f"{stem}.png"
    pdf_path = FIG_DIR / f"{stem}.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path.relative_to(REPO_ROOT)} (600 DPI) and {pdf_path.relative_to(REPO_ROOT)}")


def build_fig1_bounce_geometry() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.5))
    fig.patch.set_facecolor("#ffffff")

    # Panel (a): Signed radial extension vs. unsigned cusp at Sigma_0
    ax = axes[0]
    phase = np.linspace(0.0, 2.0 * np.pi, 500)
    u = phase / np.pi
    r_signed = np.cos(phase)
    r_unsigned = np.abs(r_signed)
    det_f = r_signed ** 3

    ax.plot(u, r_signed, color="#1f4e79", lw=2.3, label=r"Signed coordinate $R(\tau)/R_0 = \cos(\omega\tau)$")
    ax.plot(u, r_unsigned, color="#b83b2a", lw=1.9, ls="--", label=r"Unsigned radius $r(\tau)/R_0 = |R(\tau)|/R_0$")
    ax.plot(u, det_f, color="#1b7866", lw=1.8, ls=":", label=r"Orientation Jacobian $\det F = (R/R_0)^3$")
    ax.axhline(0.0, color="#444444", lw=0.9)
    ax.scatter([0.5, 1.5], [0.0, 0.0], color="#b83b2a", s=42, zorder=5)
    ax.annotate(
        r"Throat $\Sigma_0$: $\Delta\dot{R}=0$" + "\n" + r"(unsigned cusp $\Delta\dot{r}=+2R_0\omega$)",
        xy=(0.5, 0.0),
        xytext=(0.06, -0.84),
        arrowprops=dict(arrowstyle="->", color="#222222", lw=1.0),
        fontsize=8.2,
        bbox=dict(boxstyle="round,pad=0.22", fc="#ffffff", ec="#cccccc", alpha=0.92),
    )
    ax.annotate(
        r"Two crossings: $\mathcal{P}^2 = (-I_3)^2 = +I_3$",
        xy=(1.5, 0.0),
        xytext=(1.02, 0.80),
        arrowprops=dict(arrowstyle="->", color="#1b7866", lw=1.0),
        fontsize=8.0,
        color="#1b7866",
        bbox=dict(boxstyle="round,pad=0.20", fc="#ffffff", ec="none", alpha=0.90),
    )
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(-1.22, 1.46)
    ax.set_xlabel(r"Normalized phase $\omega\tau / \pi$", fontsize=9.5)
    ax.set_ylabel("Normalized radius / Jacobian", fontsize=9.5)
    ax.set_title(r"(a) Orientation-reversing extension across $\Sigma_0$", fontsize=10.2, fontweight="bold")
    ax.grid(True, alpha=0.26)
    ax.legend(fontsize=7.7, loc="upper left", framealpha=0.92)

    # Panel (b): Regular two-sheeted bounce metric & orthonormal Einstein curvature
    ax = axes[1]
    r_grid = np.linspace(-4.0, 4.0, 500)
    r_core = 1.0
    denom2 = (r_grid ** 2 + r_core ** 2) ** 2
    m_rad = 1.0 / (1.0 + (r_grid ** 2) / 3.0)
    g_tt = -(r_grid ** 4 + (4.0 / 3.0) * (r_grid ** 2) * (r_core ** 2) + r_core ** 2) / denom2
    g_rr = ((1.0 / 3.0) * (r_grid ** 4) - r_core ** 2) / denom2
    g_thth = ((1.0 / 3.0) * (r_grid ** 4) + (2.0 / 3.0) * (r_grid ** 2) * (r_core ** 2) + r_core ** 2) / denom2

    ax.plot(r_grid, m_rad, color="#6c3483", lw=1.9, ls="-.", label=r"Hybrid radial factor $M(R) = (1+R^2/3)^{-1}$")
    ax.plot(r_grid, g_tt, color="#1f4e79", lw=2.2, label=r"$R_{\mathrm{core}}^2 G_{\hat{\tau}\hat{\tau}}(R)$")
    ax.plot(r_grid, g_rr, color="#b83b2a", lw=2.1, ls="--", label=r"$R_{\mathrm{core}}^2 G_{\hat{R}\hat{R}}(R)$")
    ax.plot(r_grid, g_thth, color="#1b7866", lw=2.1, label=r"$R_{\mathrm{core}}^2 G_{\hat{\theta}\hat{\theta}}(R)$")
    ax.axhline(0.0, color="#444444", lw=0.9)
    ax.axvline(0.0, color="#777777", lw=0.9, ls=":")
    ax.annotate(
        r"Regular throat $R=0$:" + "\n" + r"$\mathcal{R}(0)=R_{\mathrm{core}}>0$, $\mathcal{K}(0)=12/R_{\mathrm{core}}^4$",
        xy=(0.0, -1.0),
        xytext=(-3.75, -1.32),
        arrowprops=dict(arrowstyle="->", color="#222222", lw=0.95),
        fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.22", fc="#ffffff", ec="#cccccc", alpha=0.92),
    )
    ax.set_xlim(-4.0, 4.0)
    ax.set_ylim(-1.48, 1.48)
    ax.set_xlabel(r"Signed radial coordinate $R / R_{\mathrm{core}}$", fontsize=9.5)
    ax.set_ylabel("Dimensionless metric / Einstein components", fontsize=9.5)
    ax.set_title(r"(b) Analytic curvature invariants across $R \in (-\infty,\infty)$", fontsize=10.2, fontweight="bold")
    ax.grid(True, alpha=0.26)
    ax.legend(fontsize=7.6, loc="upper right", framealpha=0.92)

    plt.tight_layout()
    save_both(fig, "fig1-bounce-and-geometry")


def build_fig2_variational_budget() -> None:
    atlas = nonlocal_to_local_atlas_readout()
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))
    fig.patch.set_facecolor("#ffffff")

    # Panel (a): 5-variable variational minimum
    ax = axes[0]
    chis = np.linspace(0.0, 3.6, 300)
    rows = [unified_least_resistance_cycle(float(c)) for c in chis]
    v_plus = [r["plus_infall"]["parcel_volume"] for r in rows]
    f_minus = [r["derived_cooker_takeover"]["minus_volume_fraction"] for r in rows]
    theta_star = [(((math.sqrt(16.0 * math.exp(2.0 * float(c)) + 17.0) - 3.0) / (4.0 * math.exp(float(c)))) ** 2) for c in chis]
    aspect_scaled = [0.05 * r["plus_infall"]["filament_aspect_ratio"] for r in rows]
    chi_spag = 0.5 * math.log(2.0)

    ax.plot(chis, v_plus, color="#1f4e79", lw=2.2, label=r"Parent infall volume $V_+^\star(\chi)$")
    ax.plot(chis, theta_star, color="#6c3483", lw=2.0, ls="--", label=r"Transferred winding fraction $\theta^\star(\chi)$")
    ax.plot(chis, f_minus, color="#1b7866", lw=2.2, label=r"Interior sphere share $V_-^\star/(V_+^\star+V_-^\star)$")
    ax.plot(chis, aspect_scaled, color="#c86419", lw=1.9, ls="-.", label=r"Infall aspect ratio $0.05\times(\lambda_3^\star/e^{-\chi})$")
    ax.axvline(chi_spag, color="#8e2949", lw=1.1, ls=":")
    ax.annotate(
        f"Prolate filament onset\n$\\chi_{{\\mathrm{{spag}}}}=\\frac{{1}}{{2}}\\ln 2={chi_spag:.4f}$ ($\\theta^\\star=1/2$)",
        xy=(chi_spag, 0.5),
        xytext=(0.56, 0.24),
        arrowprops=dict(arrowstyle="->", color="#8e2949", lw=1.0),
        fontsize=8.0,
        color="#8e2949",
        bbox=dict(boxstyle="round,pad=0.22", fc="#ffffff", ec="#e0c8d0", alpha=0.94),
    )
    ax.set_xlim(0.0, 3.6)
    ax.set_ylim(0.0, 1.12)
    ax.set_xlabel(r"Logarithmic compression depth $\chi$", fontsize=9.5)
    ax.set_ylabel("Equilibrium fraction / scaled ratio", fontsize=9.5)
    ax.set_title(r"(a) Five-variable variational minimum of $E = Q_{\mathrm{stretch}} + PV$", fontsize=10.1, fontweight="bold")
    ax.grid(True, alpha=0.26)
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.92)

    # Panel (b): Cosmological density shares along 9 D V B = 4 (V + B/2)^3
    ax = axes[1]
    chis_c = np.linspace(0.2, 4.2, 260)
    c_rows = [nested_zero_crossing_map(1, float(c))["nested_chain"][0] for c in chis_c]
    om_de = [100.0 * r["dark_energy_like_fraction"] for r in c_rows]
    om_dm = [100.0 * r["dark_matter_like_fraction"] for r in c_rows]
    om_vis = [100.0 * r["visible_filament_like_fraction"] for r in c_rows]
    om_dark = [100.0 * r["dark_sector_combined_fraction"] for r in c_rows]
    fp = atlas["two_valve_fixed_point"]
    chi_star = fp["chi_star"]

    ax.plot(chis_c, om_dark, color="#122033", lw=1.9, ls="--", label=r"Combined dark sector $\Omega_{\mathrm{DE}}+\Omega_{\mathrm{DM}}$")
    ax.plot(chis_c, om_de, color="#1f4e79", lw=2.3, label=r"Open-horizon vacuum share $\Omega_{\mathrm{DE}}$ ($w_{\mathrm{DE}}=-1$)")
    ax.plot(chis_c, om_dm, color="#1b7866", lw=2.2, label=r"Boundary/infall dust share $\Omega_{\mathrm{DM}}$ ($w_{\mathrm{DM}}=0$)")
    ax.plot(chis_c, om_vis, color="#c86419", lw=2.2, label=r"Baryon/filament share $\Omega_{\mathrm{vis}}$")
    ax.axvspan(3.0, chi_star, color="#8e2949", alpha=0.15, label=r"3-axis freeze-out $\chi \in [3, 3.067]$")
    ax.axvline(3.0, color="#8e2949", lw=1.1, ls=":")
    ax.scatter(
        [chi_star] * 3,
        [100.0 * fp["dark_energy_fraction"], 100.0 * fp["dark_matter_fraction"], 100.0 * fp["visible_filament_fraction"]],
        color="#8e2949",
        s=32,
        zorder=5,
    )
    ax.annotate(
        f"$\\Omega_b = e^{{-3}} = 4.98\\%$\n($\\Omega_{{\\mathrm{{DE}}}}=71.14\\%$, $\\Omega_{{\\mathrm{{DM}}}}=23.88\\%$)",
        xy=(chi_star, 71.14),
        xytext=(1.92, 31.5),
        arrowprops=dict(arrowstyle="->", color="#8e2949", lw=1.0),
        fontsize=7.8,
        color="#8e2949",
        bbox=dict(boxstyle="round,pad=0.22", fc="#ffffff", ec="#e0c8d0", alpha=0.94),
    )
    ax.set_xlim(0.2, 4.2)
    ax.set_ylim(0.0, 108.0)
    ax.set_xlabel(r"Logarithmic compression depth $\chi$", fontsize=9.5)
    ax.set_ylabel("Energy budget fraction (%)", fontsize=9.5)
    ax.set_title(r"(b) Exact budget surface $9DVB = 4(V + B/2)^3$", fontsize=10.1, fontweight="bold")
    ax.grid(True, alpha=0.26)
    ax.legend(fontsize=7.2, loc="upper left", bbox_to_anchor=(0.01, 0.72), framealpha=0.92)

    plt.tight_layout()
    save_both(fig, "fig2-variational-budget-surface")


def build_fig3_neck_cascade_and_observables() -> None:
    daughter_data = json.loads((REPO_ROOT / "results/summaries/nsc-inherited-daughter-v1.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.2))
    fig.patch.set_facecolor("#ffffff")

    # Panel (a): Whole-exterior characteristic cubic p^3 - 3p - 3 = 0
    ax = axes[0, 0]
    p_vals = np.linspace(-2.1, 2.6, 400)
    poly = p_vals ** 3 - 3.0 * p_vals - 3.0
    phi_gold = 0.5 * (1.0 + math.sqrt(5.0))
    p_plus = (phi_gold ** (2.0 / 3.0)) + (phi_gold ** (-2.0 / 3.0))

    ax.plot(p_vals, poly, color="#1f4e79", lw=2.3, label=r"$P(p) = p^3 - 3p - 3$")
    ax.axhline(0.0, color="#444444", lw=1.0)
    ax.axvline(p_plus, color="#b83b2a", lw=1.2, ls="--")
    ax.scatter([p_plus], [0.0], color="#b83b2a", s=50, zorder=5)
    ax.annotate(
        f"Unique unstable root\n$p_+ = \\varphi^{{2/3}} + \\varphi^{{-2/3}} = {p_plus:.6f}$",
        xy=(p_plus, 0.0),
        xytext=(-1.85, 4.2),
        arrowprops=dict(arrowstyle="->", color="#b83b2a", lw=1.0),
        fontsize=8.2,
        color="#b83b2a",
        bbox=dict(boxstyle="round,pad=0.24", fc="#ffffff", ec="#e6b8b0", alpha=0.94),
    )
    ax.set_xlim(-2.1, 2.6)
    ax.set_ylim(-6.5, 8.5)
    ax.set_xlabel(r"Dimensionless Laplace growth exponent $p$", fontsize=9.3)
    ax.set_ylabel(r"Boundary release residual $P(p)$", fontsize=9.3)
    ax.set_title(r"(a) Whole-exterior quadrupole neck onset $p^3 - 3p - 3 = 0$", fontsize=10.0, fontweight="bold")
    ax.grid(True, alpha=0.26)
    ax.legend(fontsize=8.0, loc="lower right")

    # Panel (b): Sequential 3-axis harmonic cascade (1 -> 2 -> 4 -> 8 lobes)
    ax = axes[0, 1]
    four_path = daughter_data["four_center_path"]
    b_sep = [row["daughter_half_separation"] for row in four_path]
    e_rel = [row["relaxed_energy_over_original"] for row in four_path]
    avail = [row["available_energy_fraction"] for row in four_path]

    ax.plot(b_sep, e_rel, "o-", color="#1f4e79", lw=2.2, markersize=5.5, label=r"Relaxed energy $E_{\mathrm{rel}}(b)/E_0$")
    ax.plot(b_sep, avail, "s--", color="#1b7866", lw=2.0, markersize=5.0, label=r"Released cascade fraction $1 - E_{\mathrm{rel}}/E_0$")
    ax.axvline(0.6274604, color="#8e2949", lw=1.2, ls=":", label=r"Transverse neck pinch $b/R = 0.6275$")
    ax.annotate(
        r"Laplace trace $\sum_k \partial_{\zeta_k} v = 0$:" + "\n"
        + r"$\mathbf{e}_z\to\mathbf{e}_x$ ratio $-1/2$" + "\n"
        + r"$\mathbf{e}_x\to\mathbf{e}_y$ ratio $-2$",
        xy=(0.6275, 0.7446),
        xytext=(0.04, 0.44),
        arrowprops=dict(arrowstyle="->", color="#8e2949", lw=1.0),
        fontsize=8.0,
        color="#122033",
        bbox=dict(boxstyle="round,pad=0.24", fc="#ffffff", ec="#cccccc", alpha=0.94),
    )
    ax.set_xlim(-0.03, 1.05)
    ax.set_ylim(0.10, 0.92)
    ax.set_xlabel(r"Transverse daughter half-separation $b / R$", fontsize=9.3)
    ax.set_ylabel("Normalized energy fraction", fontsize=9.3)
    ax.set_title(r"(b) Three-axis orthogonal cascade ($\mathbf{e}_z \to \mathbf{e}_x \to \mathbf{e}_y$)", fontsize=10.0, fontweight="bold")
    ax.grid(True, alpha=0.26)
    ax.legend(fontsize=7.6, loc="center right")

    # Panel (c): Invariant drift slope d Omega_DE / d Omega_vis along cubic surface
    ax = axes[1, 0]
    chis_d = np.linspace(1.8, 4.4, 220)
    d_rows = [nested_zero_crossing_map(1, float(c))["nested_chain"][0] for c in chis_d]
    v_arr = np.array([100.0 * r["visible_filament_like_fraction"] for r in d_rows])
    de_arr = np.array([100.0 * r["dark_energy_like_fraction"] for r in d_rows])

    row_chi3 = nested_zero_crossing_map(1, 3.0)["nested_chain"][0]
    v3 = 100.0 * row_chi3["visible_filament_like_fraction"]
    de3 = 100.0 * row_chi3["dark_energy_like_fraction"]
    v_tan = np.linspace(1.5, 11.5, 80)
    de_tan = de3 - 0.805654 * (v_tan - v3)

    ax.plot(v_arr, de_arr, color="#1f4e79", lw=2.4, label=r"Exact cubic locus $\Omega_{\mathrm{DE}}(\Omega_{\mathrm{vis}})$")
    ax.plot(v_tan, de_tan, color="#b83b2a", lw=1.7, ls="--", label=r"Tangent at $\chi=3$: $d\Omega_{\mathrm{DE}}/d\Omega_{\mathrm{vis}} = -0.806$")
    ax.scatter([v3, 100.0 * math.exp(-3.0)], [de3, 71.1441], color="#8e2949", s=42, zorder=5)
    ax.errorbar([4.93], [68.47], xerr=[0.06], yerr=[0.73], fmt="D", color="#222222", capsize=3.5, label=r"Planck 2018 base-$\Lambda$CDM")
    ax.set_xlim(1.0, 14.0)
    ax.set_ylim(63.5, 74.8)
    ax.set_xlabel(r"Visible / baryon fraction $\Omega_{\mathrm{vis}}$ (%)", fontsize=9.3)
    ax.set_ylabel(r"Vacuum stretch fraction $\Omega_{\mathrm{DE}}$ (%)", fontsize=9.3)
    ax.set_title(r"(c) Zero-parameter consistency locus & drift slope", fontsize=10.0, fontweight="bold")
    ax.grid(True, alpha=0.26)
    ax.legend(fontsize=7.6, loc="lower left")

    # Panel (d): H_0 (calibrated by omega_b) and S_8 weak-lensing comparison
    ax = axes[1, 1]
    labels = [
        "Planck 2018\nCMB ($\\Lambda$CDM)",
        "NSC ($\\chi_{\\mathrm{freeze}}=3$,\n$3.11\\%$ slip)",
        "DES-Y3 + KiDS\nJoint Lensing",
        "KiDS-1000\nCosmic Shear",
    ]
    s8_vals = [0.832, 0.776, 0.790, 0.759]
    s8_errs = [0.013, 0.010, 0.016, 0.022]
    y_pos = np.arange(len(labels))
    colors = ["#5d6d7e", "#1b7866", "#2874a6", "#2e86c1"]

    ax.axvspan(0.776 - 0.010, 0.776 + 0.010, color="#1b7866", alpha=0.16)
    for yp, val, err, col in zip(y_pos, s8_vals, s8_errs, colors):
        ax.errorbar([val], [yp], xerr=[err], fmt="o", color=col, markersize=7.5, lw=2.1, capsize=4.0)
        ax.text(val, yp - 0.25, f"{val:.3f} $\\pm$ {err:.3f}", ha="center", va="top", fontsize=7.8, fontweight="bold", color=col)

    h0_note = (
        r"Single empirical anchor $\omega_b = 0.02237 \pm 0.00015$:" + "\n"
        + r"$\Omega_b = e^{-3} \Rightarrow H_0 = 100\sqrt{\omega_b e^3} = 67.03 \pm 0.22\ \mathrm{km\,s^{-1}\,Mpc^{-1}}$" + "\n"
        + r"(Planck 2018: $67.36 \pm 0.54$; ACT DR6: $67.9 \pm 1.1$)"
    )
    ax.text(
        0.723,
        4.42,
        h0_note,
        fontsize=7.2,
        va="bottom",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.28", fc="#f8f9fa", ec="#1b7866", lw=1.0),
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8.2)
    ax.invert_yaxis()
    ax.set_xlim(0.720, 0.865)
    ax.set_ylim(4.65, -0.55)
    ax.set_xlabel(r"Clustering amplitude $S_8 \equiv \sigma_8 \sqrt{\Omega_m / 0.3}$", fontsize=9.3)
    ax.set_title(r"(d) Sub-horizon $S_8$ suppression & $\omega_b$-calibrated $H_0$", fontsize=10.0, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.26)

    plt.tight_layout()
    save_both(fig, "fig3-neck-cubic-and-cascade")


def main() -> None:
    build_fig1_bounce_geometry()
    build_fig2_variational_budget()
    build_fig3_neck_cascade_and_observables()


if __name__ == "__main__":
    main()
