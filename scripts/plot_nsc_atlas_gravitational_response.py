"""Plot the declared background, response and causal front; no new wave solve."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from recursive_horizons.nsc_atlas_gravitational_response import Region


ROOT = Path(__file__).resolve().parents[1]


def main():
    data = json.loads((ROOT/"results/summaries/nsc-atlas-gravitational-response-v1.json").read_text())
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.1))
    blue, orange, gray = "#176d8a", "#b46132", "#737a82"
    ax = axes[0, 0]
    for name, bounds, color in [("core", (.025, .25), blue),
                                 ("parent", (.25, 2.5), orange),
                                 ("exterior", (2.5, 10.), gray)]:
        g = Region(**data["regions"][name])
        r = np.geomspace(*bounds, 100)
        ax.semilogx(r, g.clock*np.sqrt(g.f(r)), color=color, lw=2, label=name.capitalize())
    for R in (.25, 2.5):
        ax.axvline(R, color="#999999", lw=.8, ls=":")
    ax.set(title="One geometry, one matched clock", xlabel="Areal radius / M",
           ylabel="Proper clock rate relative to infinity", xlim=(.025, 10.))
    ax.set_xticks([.025, .1, .25, 1., 2.5, 10.], labels=[".025", ".1", ".25", "1", "2.5", "10"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(frameon=False, loc="upper left", fontsize=9)

    ax = axes[0, 1]
    w = np.array([row["omega"] for row in data["curve"]])
    ax.plot(w, np.ones(len(w)), color=blue, lw=2, label="Regular child: lossless axial sector")
    ax.plot(w, [1-row["absorbed"] for row in data["curve"]], ".-", color=orange,
            label="Same parent, ingoing inner horizon", ms=4)
    for target in data["targets"]:
        value = target["comparison"]["horizon"]["final"]["reflected_power"]
        ax.plot(target["omega"], value, "o", color=orange, ms=7)
    ax.set(title="Same exterior mass, different returned energy", xlabel="M times global angular frequency",
           ylabel="Returned / incident canonical power", xlim=(.15, 2.45), ylim=(-.04, 1.08))
    ax.legend(frameon=False, fontsize=8, loc="lower left")

    ax = axes[1, 0]
    x = np.arange(2)
    contrast = [t["complex_response_difference"] for t in data["targets"]]
    errors = [sum(t["comparison"][k]["complex_reflection_error_indicator"] for k in ("cap", "horizon"))
              for t in data["targets"]]
    ax.semilogy(x, contrast, "o", color=blue, ms=8, label="Complex response difference")
    ax.semilogy(x, errors, "_", color=gray, ms=18, mew=2, label="Combined numerical indicator")
    for i in range(2):
        ax.plot([i, i], [errors[i], contrast[i]], color="#d8e1e4", zorder=-1)
    ax.set(title="The interior contrast is numerically resolved", ylabel="Dimensionless amplitude",
           xticks=x, xticklabels=["M omega = .4", "M omega = 2"], xlim=(-.6, 1.6))
    ax.legend(frameon=False, fontsize=9, loc="upper left")

    ax = axes[1, 1]
    first = data["causal_times"]["first_possible_difference_proper_over_M"]
    center = data["causal_times"]["center_null_roundtrip_proper_over_M"]
    xi, xc = first/2, center/2
    ax.plot([0, xi, 0], [0, xi, first], color=orange, lw=2,
            label="First changed interface")
    ax.plot([0, xc, 0], [0, xc, center], color=blue, lw=1.5, ls="--",
            label="Null round trip to center")
    ax.axvline(0, color=gray, lw=.8)
    ax.axvline(xi, color="#999999", lw=.8, ls=":")
    ax.scatter([0, 0], [first, center], color=[orange, blue], s=30)
    ax.annotate("34.31868 M", (0, first), (2.5, 29), color=orange,
                arrowprops={"arrowstyle": "-", "color": orange}, fontsize=9)
    ax.set(title="A finite causal delay before any difference", xlabel="Inward optical distance / M",
           ylabel="Receiving-clock time / M", xlim=(-.7, 20.5), ylim=(-1, 40))
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    for ax in axes.flat:
        ax.grid(alpha=.15)
    fig.tight_layout(h_pad=2.3, w_pad=2.2)
    output = ROOT/"docs/figures/nsc-atlas-gravitational-response.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
