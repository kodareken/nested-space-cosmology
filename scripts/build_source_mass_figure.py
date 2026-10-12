"""Render the recorded geometric mass coefficients without new integration."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RECORD = "results/development/nsc-atlas-source-mass-v1/mass.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def render(output: Path) -> dict:
    from record_nsc_atlas_source_mass import check
    check(ROOT / Path(RECORD).parent)
    raw = (ROOT / RECORD).read_bytes()
    calculation = json.loads(raw)["calculation"]
    data, geometry = calculation["samples"], calculation["geometry"]
    x = data["x"]
    if (len(x) != 129 or any(len(data[k]) != len(x)
            for k in ("p", "mu_jordan", "mu_einstein"))
            or not data["mu_jordan"][-1] < 0):
        raise ValueError("Stored mass samples do not support the figure.")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.unicode_minus": False,
                         "svg.hashsalt": "nsc-source-mass-2026-10-13"})
    fig, axes = plt.subplots(2, 1, figsize=(7, 4.7), sharex=True,
                             height_ratios=[.7, 1.3], layout="constrained")
    for ax in axes:
        ax.axvspan(*geometry["support_x"], color="#edf1f4", zorder=0)
        ax.axhline(0, color="#9badb8", lw=.8, zorder=1)
        ax.set_xlim(x[0], x[-1])
        ax.grid(axis="y", color="#dbe3e8", lw=.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=10)
    axes[0].plot(x, data["p"], color="#2879a0", lw=2)
    axes[0].set(ylabel=r"$p=b/x$", title="A finite preparation leaves a local geometric mass correction")
    axes[1].plot(x, data["mu_jordan"], color="#21816c", lw=2.2,
                 label="Jordan geometric mass")
    axes[1].plot(x, data["mu_einstein"], color="#865a98", lw=1.8, ls="--",
                 label="Einstein-frame comparison")
    axes[1].plot(x[-1], data["mu_jordan"][-1], "o", color="#21816c",
                 ms=4, clip_on=False)
    axes[1].legend(frameon=False, fontsize=10, loc="upper right")
    axes[1].set(ylabel=r"Quadratic coefficient $\mu$",
                xlabel=r"Anchored affine parameter $x=R_0/\ell$")
    fig.supxlabel(r"$\mathfrak{m}_{J,E}/\ell=0.4+\alpha\mu_{J,E}$,  $\alpha=\kappa_4\varepsilon^2/\ell^2$"
                  "\nInitial null surface only; no finite amplitude selected", fontsize=10)
    output.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for suffix in ("png", "svg"):
        path = output / ("source-mass." + suffix)
        metadata = {"Date": None} if suffix == "svg" else {"Software": "NSC recorded source mass"}
        fig.savefig(path, dpi=220, metadata=metadata)
        outputs[path.name] = sha(path.read_bytes())
    plt.close(fig)
    manifest = {
        "schema": "nsc.source-mass-figure.v1",
        "producer": "scripts/build_source_mass_figure.py",
        "producer_sha256": sha(Path(__file__).read_bytes()),
        "inputs": [{"path": RECORD, "sha256": sha(raw)}],
        "outputs": outputs,
        "semantics": {
            "shading": "compact support of the original characteristic trace",
            "comparison": "geometric mass in two conformally related metrics; equal after compact input",
            "negative_coefficient_is_not_negative_total_mass": True,
            "coordinate_invariant_but_frame_dependent": True,
            "no_new_scientific_calculation": True,
            "finite_amplitude_selected": False,
            "global_energy_or_evolved_wall": False,
            "zero_trace_is_not_vacuum_region": True,
        },
    }
    (output / "source-mass.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "paper/figures")
    result = render(parser.parse_args().output_dir.resolve())
    print(json.dumps({"outputs": list(result["outputs"]), "new_scientific_calculations": 0}))
