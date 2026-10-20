"""Plot intrinsic initial profiles of the two checked compact completions."""
from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_compact_parent import CompactParent


ROOT = Path(__file__).resolve().parents[1]


def side_profile(model, nested):
    specs = [("cap", 0., .8 if nested else 1.5, 1)]
    if nested:
        specs += [
            ("host_left", .8, .1, 1),
            ("host_bridge", -math.sqrt(.1-model.rb), 0., 1),
            ("host_bridge", 0., math.sqrt(.5-model.rb), 1),
            ("host_transition", .5, .9, 1), ("host_sheet", .9, 1.5, 1)]
    specs += [("outer_taper", 1.5, 1., 0), ("outer_bridge", -math.sqrt(.2), 0., 0)]
    pieces, walls, offset = [], [], 0.
    for name, low, high, H in specs:
        points = np.linspace(low, high, 41)
        if name == "cap":
            distance, radius = points-low, points
        else:
            radius = model.initial_rays(name, points)[0][:, 0]
            distance = np.array([model.quad(
                lambda x: model.initial_rays(name, x)[1], low, p) for p in points])
        pieces.append((offset+distance, radius, H))
        offset += float(distance[-1])
        if name == "host_sheet" or name == "cap":
            walls.append(offset)
    return pieces, walls, offset


def main():
    model = CompactParent(64)
    nested = side_profile(model.nested, True)
    regular = side_profile(model.regular, False)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.5), constrained_layout=True)
    for ax, right, label in zip(axes, (regular, nested),
                               ("Nested + regular: three physical walls",
                                "Nested + nested: four physical walls")):
        left_pieces, left_walls, Ll = nested
        right_pieces, right_walls, Lr = right
        for pieces, mirror in ((left_pieces, False), (right_pieces, True)):
            for x, R, H in pieces:
                if mirror:
                    x, R = Ll+Lr-x[::-1], R[::-1]
                color = "#f0c96b" if H else "#bccbd9"
                ax.fill_between(x, 0, R, color=color, alpha=.48, linewidth=0)
                ax.plot(x, R, color="#254b70", linewidth=1.8)
        for wall in left_walls+[Ll+Lr-x for x in right_walls]:
            ax.axvline(wall, color="#a44c26", linewidth=.9, linestyle="--")
        ax.scatter([0., Ll+Lr], [0., 0.], color="#254b70", s=28, zorder=5)
        ax.scatter([Ll], [.8], color="#254b70", s=24, zorder=5)
        ax.annotate("Shared parent throat", (Ll, .8), (Ll, 1.23),
                    ha="center", arrowprops={"arrowstyle": "-", "color": "#526777"},
                    color="#254b70", fontsize=9)
        ax.set_title(label, loc="left", color="#183c54", fontsize=11, fontweight="bold")
        ax.set(xlim=(-.08, Ll+Lr+.08), ylim=(-.03, 1.7), ylabel=r"$R/\ell$",
               xlabel=r"Initial proper-distance label $\lambda/\ell$")
        ax.grid(alpha=.15)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Complete compact parents: two regular centres, no exterior cutoff",
                 fontsize=14, color="#183c54")
    fig.supxlabel("Intrinsic radial profiles, not Euclidean embeddings or causal diagrams. "
                  "Gold: supplied H = 1/ell; blue-grey: Schwarzschild connector.",
                  fontsize=9, color="#526777")
    target = ROOT/"docs/figures/nsc-atlas-compact-parent.png"
    fig.savefig(target, dpi=160)
    plt.close(fig)
    print(f"Wrote {target}; {model.budget.count} charged geometry/plot queries")


if __name__ == "__main__":
    main()
