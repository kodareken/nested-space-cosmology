"""Illustrate the analytic initial data and retained transport coefficients."""
from pathlib import Path
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_vortex_transfer import KAPPA, exact_theta_rate, jets


ROOT = Path(__file__).resolve().parents[1]


def main():
    record = json.loads((ROOT/"results/summaries/nsc-atlas-vortex-transfer-v1.json").read_text(encoding="utf-8"))
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(11, 7.3))
    grid = fig.add_gridspec(2, 2, hspace=.48, wspace=.30)
    blue, orange = "#176d8a", "#b46132"
    theta = np.linspace(0, 2*math.pi, 385)
    radius = np.linspace(.5, 1.5, 121)
    rr, tt = np.meshgrid(radius, theta, indexing="ij")
    for column, mode in enumerate((4, 16)):
        ax = fig.add_subplot(grid[0, column], projection="polar")
        W = jets(rr, np.zeros_like(rr), mode)[0]
        velocity = np.real(W[1]*np.exp(1j*(mode*tt+KAPPA*mode*np.log(rr))))
        mesh = ax.pcolormesh(tt, rr, velocity, shading="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        ax.set_ylim(0, 1.5)
        ax.set_yticks([.5, 1., 1.5])
        ax.set_yticklabels([".5", "1", "1.5"], fontsize=8)
        ax.set_xticklabels([])
        ax.grid(alpha=.25)
        ax.set_title(f"Prepared angular mode m = {mode}\nZero mean velocity", pad=12)
    color = fig.colorbar(mesh, ax=fig.axes[:2], shrink=.7, pad=.08, aspect=22)
    color.set_label("Initial tangential velocity / U")

    ax = fig.add_subplot(grid[1, 0])
    r = np.linspace(.5, 1.5, 301)
    ax.plot(r, exact_theta_rate(r), color=blue, lw=2,
            label="Exact mean rate for every m and viscosity")
    targets = np.array([.75, 1.25])
    rates = exact_theta_rate(targets)
    ax.scatter(targets, rates, color=orange, zorder=3)
    ax.axhline(0, color="#888888", lw=.7)
    ax.set(title="Fine correlations drive opposite rotating bands",
           xlabel="Radius / L, at the midplane", ylabel="Initial mean acceleration / (U squared / L)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.grid(alpha=.15)

    ax = fig.add_subplot(grid[1, 1])
    m = [row["mode"] for row in record["mode_results"]]
    ratio = [1/row["energy_over_initial_loss"] for row in record["mode_results"]]
    ax.loglog(m, ratio, "o-", color=orange, lw=2)
    ax.set_xticks(m, labels=[str(j) for j in m])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set(title="Finer motion pays a larger viscous cost", xlabel="Angular mode m",
           ylabel="Initial fractional energy loss per L/U")
    ax.text(.05, .93, "Energy stays finite; loss grows as m squared", transform=ax.transAxes,
            va="top", fontsize=9)
    ax.grid(alpha=.15, which="both")
    output = ROOT/"docs/figures/nsc-atlas-vortex-transfer.png"
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
