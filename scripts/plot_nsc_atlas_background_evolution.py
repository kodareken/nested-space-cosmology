#!/usr/bin/env python3
"""Render one saved-array comparison; no source, metric, or time evolution."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
import scipy

ROOT = Path(__file__).resolve().parents[1]
PRODUCER = "scripts/plot_nsc_atlas_background_evolution.py"
DEFAULT_OUTPUT = "docs/figures/nsc-atlas-background-evolution.png"
INPUTS = (
    ("nsc-atlas-background-evolution-v1", "background.json"),
    ("nsc-atlas-compatible-evolution-v1", "compatible.json"),
)
SCHEMA = "nsc-atlas-background-evolution-figure-v1"


def identity(path):
    raw = path.read_bytes()
    return dict(sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw))


class Allowance:
    """Count array operands/query points, not files or rendering pixels."""
    def __init__(self):
        self.started = time.process_time()
        self.points = 0

    def check(self, points=0):
        self.points += int(points)
        if self.points > 100000 or time.process_time()-self.started > 30:
            raise RuntimeError("saved-array figure allowance exhausted")


def retained_record(directory, filename, allowance):
    folder = ROOT/"results/development"/directory
    path, manifest_path = folder/filename, folder/"manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    file_identity = identity(path)
    if (manifest.get("schema") != directory+"-manifest"
            or manifest.get("outputs") != {filename: file_identity}):
        raise ValueError("retained figure input hash mismatch")
    for relative, expected in manifest["source_bindings"].items():
        source = (ROOT/relative).resolve()
        if not source.is_relative_to(ROOT) or identity(source)["sha256"] != expected:
            raise ValueError("retained figure source binding mismatch: "+relative)
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != directory:
        raise ValueError("retained figure schema mismatch")
    allowance.check()
    return data, {str(p.relative_to(ROOT)).replace("\\", "/"): identity(p)
                  for p in (path, manifest_path)}


def selected(data, amplitude, cells, steps):
    matches = [c for c in data["cases"] if
               (c["amplitude"], c["cells"], c["steps"]) == (amplitude, cells, steps)]
    if len(matches) != 1:
        raise ValueError("figure case selection is not unique")
    return matches[0]


def wall_arrays(case, allowance):
    rows = case["snapshots"]
    out = {key: np.asarray([s["physical"]["wall"][key] for s in rows])
           for key in ("r", "T", "c")}
    out["x"] = np.asarray([s["x"] for s in rows])
    if any(not np.all(np.isfinite(a)) for a in out.values()) or np.any(np.diff(out["x"]) <= 0):
        raise ValueError("non-finite or unordered saved wall history")
    allowance.check(4*len(rows))
    return out


def render(output):
    allowance = Allowance()
    output = Path(output).resolve()
    manifest_path = output.with_suffix(".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError("figure outputs are creation-only")
    if output.suffix.lower() != ".png" or not output.is_relative_to(ROOT):
        raise ValueError("figure must be a PNG inside the repository")
    producer = identity(ROOT/PRODUCER)
    new, new_inputs = retained_record(*INPUTS[0], allowance)
    old, old_inputs = retained_record(*INPUTS[1], allowance)
    parent = new["parent_input"]
    parent_path = "results/development/"+INPUTS[1][0]+"/"
    if (parent["sha256"] != old_inputs[parent_path+"compatible.json"]["sha256"]
            or parent["manifest_sha256"] != old_inputs[parent_path+"manifest.json"]["sha256"]):
        raise ValueError("background figure input has a different compatible parent")

    grids = ((16, 256), (32, 256), (32, 512))
    histories = []
    for cells, steps in grids:
        zero = wall_arrays(selected(new, "zero", cells, steps), allowance)
        finite = wall_arrays(selected(new, "finite", cells, steps), allowance)
        if not np.array_equal(zero["x"], finite["x"]):
            raise ValueError("paired wall coordinates differ")
        histories.append(dict(x=zero["x"], r=finite["r"]-zero["r"], c=finite["c"]-zero["c"]))
        allowance.check(4*len(zero["x"]))
    if any(not np.array_equal(h["x"], histories[0]["x"]) for h in histories[1:]):
        raise ValueError("comparison histories do not share x")
    spatial = np.abs(histories[1]["c"]-histories[0]["c"])
    temporal = np.abs(histories[2]["c"]-histories[1]["c"])
    allowance.check(4*len(spatial))

    new_case = selected(new, "zero", 32, 256)
    old_case = selected(old, "zero", 32, 256)
    wall = wall_arrays(new_case, allowance)
    q = wall["T"]+np.arctanh(wall["r"])
    if np.any(np.diff(q) <= 0) or np.max(wall["T"]) >= q[0]:
        raise ValueError("saved incoming oracle has no certified quiet reflected branch")
    y = np.linspace(0., 1., 33)
    old_final, new_final = old_case["snapshots"][-1], new_case["snapshots"][-1]
    if any(old_final["physical"]["wall"][k] != new_final["physical"]["wall"][k] for k in ("T", "r")):
        raise ValueError("zero final cap coordinates differ")
    old_chi = np.r_[0., np.cumsum(old_final["state"]["cap"]["Iq"])]
    new_chi = np.r_[0., np.cumsum(new_final["physical"]["cap"]["Iq"])]
    interpolant = PchipInterpolator(q, wall["c"], extrapolate=False)

    def incoming(y_values):
        labels = wall["T"][-1]+np.arctanh(wall["r"][-1]*y_values)
        if np.any(labels > q[-1]):
            raise ValueError("incoming diagnostic would extrapolate the retained wall trace")
        result = np.zeros_like(labels)
        active = labels >= q[0]
        result[active] = interpolant(labels[active])
        allowance.check(5*len(labels))
        return result

    y_smooth = np.linspace(0., 1., 513)
    oracle_smooth, oracle_nodes = incoming(y_smooth), incoming(y)
    allowance.check(6*len(y)+3*len(q))
    colors = ("#b36b24", "#167d8d", "#273d75")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 11, "axes.labelsize": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.5))
    fig.subplots_adjust(left=.085, right=.98, top=.865, bottom=.13, wspace=.27, hspace=.39)
    fig.suptitle("Nested-space evolution: saved finite-window comparison", fontsize=16,
                 fontweight="bold", y=.965, color="#17283a")
    fig.text(.5, .92, "Dimensionless coordinates and normalized fields  •  one common prescribed source",
             ha="center", fontsize=11, color="#4b5968")
    ax = axes[0, 0]
    ax.plot(y_smooth, oracle_smooth, color="#222c36", linewidth=2, label="Incoming characteristic (new trace; PCHIP)")
    ax.plot(y, old_chi, "s--", color=colors[0], markersize=3.5, linewidth=1.4, label="Previous compatible cap")
    ax.plot(y, new_chi, "o", color=colors[1], markersize=3.5, label="Characteristic-background cap")
    ax.set(title=r"A  Final zero-amplitude cap: 32 cells, 256 steps", xlabel=r"Cap position $y=r/R_W$", ylabel=r"Reduced field $\chi$ ($\alpha=0$)", xlim=(0, 1))
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    labels = ("16 cells / 256 steps", "32 cells / 256 steps", "32 cells / 512 steps")
    styles = ("--", "-", ":")
    for j, h in enumerate(histories):
        axes[0, 1].plot(h["x"], 1e4*h["r"], styles[j], color=colors[j], linewidth=1.8, label=labels[j])
        axes[1, 0].plot(h["x"], 1e5*h["c"], styles[j], color=colors[j], linewidth=1.8, label=labels[j])
    axes[0, 1].set(title="B  Paired wall-radius response", xlabel=r"Fixed incoming source coordinate $x$", ylabel=r"$\Delta r_W\;[10^{-4}]$")
    axes[0, 1].legend(fontsize=8.5, frameon=False, loc="lower left")
    axes[1, 0].set(title="C  Paired wall-field response", xlabel=r"Fixed incoming source coordinate $x$", ylabel=r"$\Delta c_W\;[10^{-5}]$")
    axes[1, 0].legend(fontsize=8.5, frameon=False, loc="upper left")
    axes[1, 1].plot(histories[0]["x"], 1e6*spatial, color=colors[1], linewidth=1.8,
                    label="Spatial: 32 − 16 cells (256 steps)")
    axes[1, 1].plot(histories[0]["x"], 1e6*temporal, color=colors[2], linewidth=1.8,
                    label="Temporal: 512 − 256 steps (32 cells)")
    axes[1, 1].set(title="D  Absolute differences of paired wall fields", xlabel=r"Fixed incoming source coordinate $x$", ylabel=r"$|\Delta c_W\ \mathrm{difference}|\;[10^{-6}]$")
    axes[1, 1].legend(fontsize=8, frameon=False, loc="upper left")
    for ax in axes.flat:
        ax.grid(alpha=.18, linewidth=.7)
        ax.tick_params(labelsize=9)
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.set_xlim(histories[0]["x"][[0, -1]])
    alpha = selected(new, "finite", 32, 256)["alpha"]
    fig.text(.085, .065,
             rf"$\Delta=$ finite $\alpha$ − zero $\alpha$ at matching $x$; finite $\alpha={alpha:.5g}$. "
             "Cap curves share the final zero-case radius and central time.", fontsize=9, color="#46525e")
    fig.text(.085, .035, "Incoming curve interpolates the new retained wall trace. These are finite-grid diagnostics, not error bounds or a central-arrival result.",
             fontsize=9, color="#46525e")
    allowance.check(20*(513+6*129+2*33))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, facecolor="white", metadata={"Software": PRODUCER})
    plt.close(fig)
    allowance.check()
    if identity(ROOT/PRODUCER) != producer:
        raise RuntimeError("figure producer changed during rendering")
    metadata = dict(schema=SCHEMA, producer={PRODUCER: producer}, inputs={**new_inputs, **old_inputs},
        outputs={output.relative_to(ROOT).as_posix(): identity(output)},
        runtime=dict(numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
        execution=dict(cpu_seconds=time.process_time()-allowance.started, arithmetic_points=allowance.points,
            cpu_limit_seconds=30, point_limit=100000, renders=1, source_integrations=0, time_evolutions=0,
            recursive_scientific_authentications=0, point_count_scope="array operands and interpolation query points; excludes file bytes and rendering pixels"),
        comparison=dict(alpha=alpha, grids=[dict(cells=n, steps=k) for n,k in grids],
            cap_case="zero / 32 cells / 256 steps", incoming="PCHIP of new zero32/256 retained wall c versus T+atanh(R); quiet reflected branch",
            max_abs_new_cap_minus_incoming=float(np.max(abs(new_chi-oracle_nodes))),
            max_abs_old_cap_minus_new_incoming=float(np.max(abs(old_chi-oracle_nodes))),
            max_abs_paired_c_spatial=float(np.max(spatial)), max_abs_paired_c_time=float(np.max(temporal)),
            scope="saved finite-window diagnostic; no error-bound or central-arrival claim"))
    manifest_path.write_text(json.dumps(metadata, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(render(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
