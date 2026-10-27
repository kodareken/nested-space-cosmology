"""Plot retained mode/vertex evidence; never solve a new mode in the renderer."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def main():
    d = json.loads((ROOT/"results/summaries/nsc-atlas-cosmohedron-interaction-v1.json").read_text(encoding="utf-8"))
    a, b = [f["levels"][-1] for f in d["families"]]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.0))
    blue, orange, gray = "#176d8a", "#b46132", "#737a82"
    ax = axes[0, 0]
    ax.set(xlim=(0, 1), ylim=(0, 1), title="One actual gravitational mode connects two currents")
    for x, text, color in [(.27, "Inner wall\nR = .25 M", blue), (.73, "Outer wall\nR = 2.5 M", orange)]:
        ax.scatter([x], [.5], s=70, color=color, zorder=3)
        ax.text(x, .23, text, ha="center", fontsize=9, color=color)
    ax.plot([.27, .73], [.5, .5], color=gray, lw=2)
    ax.text(.5, .6, "M omega = .392804", ha="center", fontsize=9)
    for x, outer in [(.27, .04), (.73, .96)]:
        for y in (.78, .22):
            ax.plot([outer, x], [y, .5], color=gray, lw=1)
    ax.text(.1, .87, "Two scalar\nharmonics", ha="center", fontsize=8)
    ax.text(.9, .87, "Two scalar\nharmonics", ha="center", fontsize=8)
    ax.text(.5, .03, "Graph kernel: 1 / (Et EL ER)", ha="center", fontsize=10)
    ax.axis("off")

    ax = axes[0, 1]
    x = np.arange(2)
    for data, offset, color, label in [(a, -.17, blue, "Child radius .25 M"),
                                        (b, .17, orange, "Child radius .30 M")]:
        ax.bar(x+offset, [v["current_coupling"] for v in data["vertices"]], width=.3,
               color=color, label=label)
    ax.set_yscale("log")
    ax.set(title="Vertices are calculated from the field action", xticks=x,
           xticklabels=["Inner current", "Outer current"], ylabel="c M to the power 5/2 / sqrt(G)")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.grid(alpha=.15, axis="y")

    ax = axes[1, 0]
    ax.bar([0, 1], [-a["coefficient_shape"]*1e7, -b["coefficient_shape"]*1e7],
           color=[blue, orange], width=.5)
    ax.set(title="A hidden-interior change alters the selected exchange",
           xticks=[0, 1], xticklabels=[".25 M child", ".30 M child"],
           ylabel="Magnitude of Gamma M^4/G, in 1e-7 units", ylim=(0, 8.))
    ax.text(.5, 7.1, "64.84% coefficient contrast", ha="center", fontsize=10)
    ax.grid(alpha=.15, axis="y")

    ax = axes[1, 1]
    scale = np.linspace(1, 2, 60)
    ax.loglog(scale, scale**-4, color=blue, label="Selected coefficient: M^-4")
    ax.loglog(scale, scale**-2, color=orange, label="Normalized selected moment: M^-2")
    c = d["mass_rescaling_control"]
    ax.scatter([1, 2], [1, c["log_wavefunction_exchange_coefficient"]/a["log_wavefunction_exchange_coefficient"]], color=blue)
    ax.scatter([1, 2], [1, c["selected_exchange_correlation"]/a["selected_exchange_correlation"]], color=orange)
    ax.set(title="Geometric shape inherits; quantum strength runs", xlabel="Mass-length ratio at fixed G",
           ylabel="Value relative to the M = 1 case", xticks=[1, 2], xticklabels=["1", "2"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.grid(alpha=.15, which="both")
    fig.tight_layout(h_pad=2.7, w_pad=2.3)
    output = ROOT/"docs/figures/nsc-atlas-cosmohedron-interaction.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
