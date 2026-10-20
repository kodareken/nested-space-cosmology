"""Draw an original comparison plate from the parent-observable equations."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_atlas_parent_observables import (
    DustParent, critical_tension_wall, de_sitter_bulk_slice,
    shadow_lensing, tension_potential,
)

ROOT = Path(__file__).resolve().parents[1]


def main():
    ink, blue, gold, red = "#172e43", "#247e8d", "#bb8132", "#b44f4c"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.labelcolor": ink, "text.color": ink,
                         "axes.edgecolor": "#b6bdc3", "xtick.color": ink,
                         "ytick.color": ink, "axes.titlesize": 12,
                         "legend.fontsize": 9.5})
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.4))
    fig.patch.set_facecolor("#faf9f5")
    fig.subplots_adjust(left=.085, right=.975, bottom=.16, top=.87, hspace=.56, wspace=.29)
    fig.text(.085, .955, "SURROUNDINGS, WALLS AND WHAT A CLOCK CAN SEE",
             fontsize=16, weight="bold", color=ink)
    fig.text(.085, .919, "Four calculated views | each panel retains its own physical assumptions",
             fontsize=11, color="#526673")
    p = DustParent()
    rec = p.arrival()
    t = np.linspace(0., rec["reception_time"], 180)
    ray = solve_ivp(lambda tau, y: [(1-.8/y[0])/p.state(tau)["f"]],
                    (0., t[-1]), [2.], t_eval=t, method="DOP853",
                    rtol=2e-12, atol=2e-13)
    a = ax[0, 0]
    a.plot(t, [p.state(s)["radius"] for s in t], color=blue, lw=2.4, label="Comoving dust join")
    a.plot(t, ray.y[0], color=gold, lw=2.4, label="Outgoing light")
    a.axhline(2., color=ink, lw=1, ls=":", label="Static emitter")
    a.scatter([t[-1]], [rec["reception_radius"]], color=red, zorder=5)
    a.set(xlabel=r"Parent proper time $\tau/\ell$", ylabel=r"Areal radius $R/\ell$",
          title="A  |  Expanding parent, retained 4D exterior")
    a.text(.04, .76, r"$1+z=$"+f'{rec["one_plus_redshift"]:.6f}',
           transform=a.transAxes, color=red, fontsize=12)
    a.legend(loc="lower right", frameon=False)

    a = ax[0, 1]
    f = np.linspace(0., 1., 100)
    a.plot(f, [shadow_lensing(s)["lensing_over_GR_same_dynamical_mass"] for s in f],
           color=blue, lw=2.4)
    a.axhline(1., color=ink, lw=1, ls=":")
    a.scatter([1.], [.75], color=red, zorder=5)
    a.text(.05, .12, "Pure shadow source: 0.75\nZero mode; unstabilized radion",
           transform=a.transAxes, fontsize=10, color=ink)
    a.set(xlabel="Shadow fraction of dynamical potential", ylabel="Lensing / GR at equal dynamical mass",
          ylim=(.72, 1.035), title="B  |  Separate 5D two-brane response")

    a = ax[1, 0]
    x = np.linspace(.55, 2.1, 250)
    crit = critical_tension_wall(.25)
    a.plot(x, tension_potential(x, mu=.8, u=.25)["V"], color=blue, lw=2.4,
           label=r"Retained $kM=0.8$")
    a.plot(x, tension_potential(x, mu=crit["mu"], u=.25)["V"], color=gold, lw=2.4,
           label=f'Critical kM = {crit["mu"]:.6f}')
    a.axhline(0., color=ink, lw=1, ls=":")
    a.scatter([crit["x"]], [0.], color=red, zorder=5)
    a.annotate("Unstable maximum", xy=(crit["x"], 0.), xytext=(1.38, .3),
               fontsize=10, arrowprops={"arrowstyle": "-", "color": red}, color=red)
    a.set(xlabel=r"Wall radius $x=kR$", ylabel=r"$V(x)$, with $\dot R^2+V=0$",
          ylim=(-2., .55), title="C  |  Original constant-tension wall")
    a.legend(loc="lower center", frameon=False)

    a = ax[1, 1]
    y = np.linspace(0., 1.48, 150)
    a.plot(y, np.exp(-y), color=gold, lw=2.4, label=r"$s=1$: flat induced geometry")
    warp = de_sitter_bulk_slice(y)
    a.plot(y, warp["warp"], color=blue, lw=2.4, label=r"$s=1.1$: de Sitter induced geometry")
    a.text(.43, .64, r"$\Lambda_4L_5^2=3(s^2-1)$"+"\n"+r"$=0.63$ at supplied $s=1.1$",
           transform=a.transAxes, fontsize=11, color=ink)
    a.set(xlabel=r"Bulk normal distance $y/L_5$", ylabel=r"Warp factor $A(y)$",
          title="D  |  Separate 5D constant-curvature bulk")
    a.legend(loc="lower left", frameon=False, fontsize=9)
    for a in ax.flat:
        a.set_facecolor("#ffffff")
        a.spines[["top", "right"]].set_visible(False)
        a.grid(alpha=.12)
    fig.text(.085, .048,
             "A: Einstein--Straus matching   B: Garriga--Tanaka   C: vacuum-bubble junctions\n"
             "D: Gauss--Israel projection. Curves are recomputed here; the four panels are not one combined solution.",
             fontsize=9.5, color="#526673")
    output = ROOT/"docs/figures/nsc-atlas-parent-observables.png"
    fig.savefig(output, dpi=190, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
