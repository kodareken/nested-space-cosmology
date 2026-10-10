"""Figure for the Kasner shear-production result, from kasner_shear_production.json (corrected C) and
kasner_stopping_scales.json. (a) probability a mode is excited vs the adiabaticity parameter;
(b) pairs created vs distance light travels before the crunch, 10 solar-mass hole, one Dirac field."""
import json
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

G, c, hbar, Msun = 6.67430e-11, 2.99792458e8, 1.054571817e-34, 1.98847e30
lP = math.sqrt(hbar * G / c**3)
rg = 2 * G * 10 * Msun / c**2
res = json.load(open("kasner_shear_production.json"))
sc = json.load(open("kasner_stopping_scales.json"))["species1.0_uc1.0"]
eck = json.load(open("eck_collapse_consistency_check.json"))["10"]
lam, F = np.array(res["lambda_grid"]), np.array(res["F"])
Cs = {k: v for k, v in res["C_pairs_density_coeff"].items()}
N = lambda x, C: 6 * math.pi * C * rg**2 / x**2

try:
    apply_figure_style()
except NameError:
    pass
blue, red = "#1f5fa8", "#c0392b"
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3), layout="constrained")
ax = axs[0]
ax.plot(lam, np.maximum(F, 1e-12), color=blue, lw=2.2)
ax.axhline(0.5, color="0.55", lw=0.9, ls="--")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-3, 1e1); ax.set_ylim(1e-6, 2)
ax.set_xlabel(r"adiabaticity parameter $\lambda = k_x \tau_0\, r^{4/3}$")
ax.set_ylabel("probability a mode is excited")
ax.text(1.4e-3, 0.75, "sudden limit 1/2", color="0.35", va="bottom")
ax.text(1.0, 2.5e-2, "stays empty", color="0.35", ha="left")
ax.set_title("Only slowly varying modes\nare produced", loc="left")
ax.text(-0.14, 1.05, "a", transform=ax.transAxes, fontweight="bold", fontsize=12)

ax = axs[1]
xs = np.geomspace(1e-7, 1e-38, 300)
x_br = sc["x_br"] * lP
lo, mid, hi = (N(xs, Cs[k]) for k in ("0.5", "1.0", "2.0"))
ax.fill_between(xs, lo, hi, color=blue, alpha=0.2, lw=0)
ok = xs >= x_br
ax.plot(xs[ok], mid[ok], color=blue, lw=2.2)
ax.plot(xs[~ok], mid[~ok], color=blue, lw=2.2, ls=(0, (3, 2)))
ax.axhline(eck["N_needed_at_bounce"], color=red, lw=1.5)
ax.axhline(eck["N_max_excitable_modes"], color="0.5", lw=1.0, ls="--")
ax.axvline(lP, color="0.5", lw=1.0, ls=":")
ax.axvline(x_br, color="0.5", lw=1.0, ls=":")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-7, 1e-38); ax.set_ylim(1e10, 1e82)
ax.set_xlabel("distance light travels before the crunch (m); later is to the right")
ax.set_ylabel("pairs created (one Dirac field)")
ax.text(1e-7 * 0.6, eck["N_needed_at_bounce"] * 0.35, "needed in the paper: 1.8e78", color=red, va="top", ha="left")
ax.text(1e-34, eck["N_max_excitable_modes"] / 2.2, "isotropic bound: 2.8e39", color="0.35", ha="right", va="top")
ax.text(2e-14, 3e20, "shear production,\n10 solar-mass hole", color=blue, ha="left", va="center")
ax.text(lP * 1.25, 3e11, "Planck length", rotation=90, color="0.35", ha="right", va="bottom")
ax.text(x_br * 0.8, 3e11, "shear energy used up", rotation=90, color="0.35", ha="left", va="bottom")
ax.set_title("Count reaches the needed size only\nnear where the shear energy runs out", loc="left")
ax.text(-0.14, 1.05, "b", transform=ax.transAxes, fontweight="bold", fontsize=12)
fig.savefig("kasner_shear_production.png", dpi=300)

rend = fig.canvas.get_renderer()
axb = axs[1]
for t in axb.texts:
    bb = t.get_window_extent(rend)
    pts = [axb.transData.transform(np.column_stack([xs, y])) for y in (lo, mid, hi)]
    n = sum(int(((p[:, 0] >= bb.x0) & (p[:, 0] <= bb.x1) & (p[:, 1] >= bb.y0) & (p[:, 1] <= bb.y1)).sum()) for p in pts)
    inside = bb.x0 >= axb.get_window_extent(rend).x0 and bb.x1 <= axb.get_window_extent(rend).x1
    print(repr(t.get_text()[:22]), "curve points inside:", n, "inside axes:", inside)
txs = [t for t in axb.texts]
print("overlaps:", [(a.get_text()[:12], b.get_text()[:12]) for i, a in enumerate(txs) for b in txs[i + 1:]
                    if a.get_window_extent(rend).overlaps(b.get_window_extent(rend))])
