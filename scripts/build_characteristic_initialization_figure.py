"""Plot stored null-constraint coefficients; no scientific integration occurs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RECORD = "results/development/nsc-atlas-characteristic-initialization-v1/initialization.json"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def render(output: Path) -> dict:
    raw = (ROOT / RECORD).read_bytes()
    record = json.loads(raw)
    calculation = record["calculation"]
    data, geometry = calculation["samples"], calculation["geometry"]
    x = data["x"]
    if (len(x) != 129 or len(data["h"]) != len(x)
            or len(data["h_prime"]) != len(x)
            or data["h"][-1] >= 0 or data["h_prime"][-1] >= 0):
        raise ValueError("Stored data do not support the characteristic figure.")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.unicode_minus": False,
                         "svg.hashsalt": "nsc-initial-characteristic-2026-10-13"})
    fig = plt.figure(figsize=(7, 5.0), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[.8, 1.2])
    source = fig.add_subplot(grid[0, :])
    radius = fig.add_subplot(grid[1, 0])
    slope = fig.add_subplot(grid[1, 1])
    axes = (source, radius, slope)
    for ax in axes:
        ax.axvspan(*geometry["support_x"], color="#eef2f5", zorder=0)
        ax.axhline(0, color="#9badb8", lw=.8, zorder=1)
        ax.set_xlim(x[0], x[-1])
        ax.grid(axis="y", color="#dbe3e8", lw=.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=10)
    source.plot(x, data["p"], color="#2879a0", lw=2)
    source.set(ylabel=r"$p=b/x$", title="One initial null surface: prescribed field and forced geometry")
    source.tick_params(labelbottom=False)
    radius.plot(x, data["h"], color="#21816c", lw=2)
    slope.plot(x, data["h_prime"], color="#865a98", lw=2)
    radius.set(title="Radius correction", ylabel=r"$h$", xlabel=r"Anchored affine parameter $x=R_0/\ell$")
    slope.set(title="Radial-slope correction", ylabel=r"$h'$", xlabel=r"Anchored affine parameter $x=R_0/\ell$")
    for ax, key, color in ((radius, "h", "#21816c"), (slope, "h_prime", "#865a98")):
        ax.plot(x[-1], data[key][-1], "o", color=color, ms=4, clip_on=False)
    fig.supxlabel(r"$R/\ell=x+\alpha h$,  $R_x/\ell=1+\alpha h'$;  $\alpha=\kappa_4\varepsilon^2/\ell^2$", fontsize=11)
    output.mkdir(parents=True, exist_ok=True)
    files = {}
    for suffix in ("png", "svg"):
        path = output / ("characteristic-initialization." + suffix)
        metadata = {"Date": None} if suffix == "svg" else {"Software": "NSC stored characteristic data"}
        fig.savefig(path, dpi=220, metadata=metadata)
        files[path.name] = sha(path.read_bytes())
    plt.close(fig)
    manifest = {
        "schema": "nsc.characteristic-initialization-figure.v1",
        "producer": "scripts/build_characteristic_initialization_figure.py",
        "producer_sha256": sha(Path(__file__).read_bytes()),
        "inputs": [{"path": RECORD, "sha256": sha(raw)}],
        "outputs": files,
        "semantics": {
            "source": "Original C3 characteristic preparation; shaded interval is its support",
            "corrections": "Quadratic coefficients at equal corner-anchored affine parameter",
            "no_new_scientific_calculation": True,
            "finite_amplitude_selected": False,
            "mass_or_wall_response": False,
            "full_coupled_evolution": False,
            "zero_trace_is_not_vacuum_region": True,
        },
    }
    (output / "characteristic-initialization.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "paper/figures")
    args = parser.parse_args()
    result = render(args.output_dir.resolve())
    print(json.dumps({"outputs": list(result["outputs"]), "new_scientific_calculations": 0}))
