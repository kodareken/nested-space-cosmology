"""Vector figures from pinned finite-regeneration records; no simulation."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import io
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "paper/finite-regeneration/figures"
PRIMARY = "nf512_dt_0.0005"
SOURCE = "lab/results/development/"
FILES = {
    "episode": SOURCE + "nsc-spherical-conformal-episode-v2.json",
    "curvature": SOURCE + "nsc-spherical-conformal-curvature-v1.json",
    "local": SOURCE + "nsc-conformal-local-response-v2.json",
    "local_npz": SOURCE + "nsc-conformal-local-response-v2.npz",
    "memory": SOURCE + "nsc-conformal-memory-control-v1.json",
    "memory_npz": SOURCE + "nsc-conformal-memory-control-v1.npz",
}
NAMES = ("regional-transfer.pdf", "geometry-curvature.pdf", "normal-ledger.pdf", "local-response.pdf")
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#333333")


def git_read(commit, path):
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)


def digest(data):
    return sha256(data).hexdigest()


def cumulative(times, values):
    return np.r_[0., np.cumsum(.5 * np.diff(times) * (values[:-1] + values[1:]))]


def axes():
    fig, panels = plt.subplots(1, 2, figsize=(7.2, 2.65), layout="constrained")
    for ax in panels:
        ax.grid(True, color="#d1d5db", linewidth=.4, alpha=.65)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("Coordinate time $T$")
        ax.tick_params(labelsize=8)
    return fig, panels


def local_window(ax):
    ax.axvspan(.2, .3, color="#9ca3af", alpha=.13, linewidth=0, zorder=-1)


def generate(science_commit):
    # Every plotted value is read from this commit, never from a mutable live file.
    raw = {name: git_read(science_commit, path) for name, path in FILES.items()}
    episode = json.loads(raw["episode"])
    curvature = json.loads(raw["curvature"])
    local = json.loads(raw["local"])
    memory = json.loads(raw["memory"])
    if not episode["all_target_reached"] or not local["all_numerical_movements_within_one_percent"]:
        raise ValueError("completed numerical sources are required")
    if curvature["auxiliary_substituted"] is not False or not memory["error_within_one_percent"]:
        raise ValueError("metric curvature and independent memory reference required")
    if local["payload_sha256"] != digest(raw["local_npz"]) or memory["payload_sha256"] != digest(raw["memory_npz"]):
        raise ValueError("local source payload binding differs")
    with np.load(io.BytesIO(raw["local_npz"]), allow_pickle=False) as stored:
        series = {name: np.array(stored[name]) for name in stored.files}
    with np.load(io.BytesIO(raw["memory_npz"]), allow_pickle=False) as stored:
        memory_reference = np.array(stored["occupation_reference_substeps_4"])
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
        "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 7.2,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
        "mathtext.fontset": "dejavusans", "axes.linewidth": .6})
    FOLDER.mkdir(parents=True, exist_ok=True)
    outputs = {}

    def save(fig, name, title, dependencies, transformations):
        epoch = datetime(1970, 1, 1, tzinfo=timezone.utc)
        path = FOLDER / name
        fig.savefig(path, metadata={"Title": title, "Author": "Douglas Ek",
            "Creator": f"Matplotlib {matplotlib.__version__}", "CreationDate": epoch,
            "ModDate": epoch, "Subject": "Pinned finite spherical realization; numerical indicators"})
        plt.close(fig)
        outputs[name] = {"sha256": digest(path.read_bytes()), "bytes": path.stat().st_size,
            "sources": {FILES[key]: digest(raw[key]) for key in dependencies},
            "transformations": transformations, "vector_pdf": True,
            "metadata_epoch": 0, "embedded_searchable_font_requested": True}

    physical = episode["physical_series"][PRIMARY]
    times = np.array([row["time"] for row in physical])
    fig, (left, right) = axes()
    for color, x, window, channel in zip(COLORS, (0, 2, 4), (0, 0, 1),
            ("normal_left_outward_flux", "normal_right_outward_flux", "normal_right_outward_flux")):
        values = np.abs([row["windows"][window][channel] for row in physical])
        left.plot(times, values, color=color, lw=1.5, label=f"$x={x}$")
    left.set(title="(a) Resolved regional surfaces", ylabel="Normal energy flux magnitude")
    left.legend(loc="upper left", frameon=False, ncol=3)
    right.plot(times, [row["probability_leader_share"] for row in physical], color=COLORS[0], lw=1.5, label="Probability share")
    right.plot(times, [row["shell_leader_share"] for row in physical], color=COLORS[1], lw=1.5, ls="--", label="Normal-energy share")
    right.set(title="(b) Region 0 stays the leader", ylabel="Leader's regional share")
    right.legend(loc="best", frameon=False)
    local_window(left); local_window(right)
    save(fig, NAMES[0], "Regional transfer and maintained localization", ("episode",),
         "Absolute physical surface fluxes from physical_series, whose producer already divides flux_nodal by dx; stored region-0 shares. Grey interval is local-response window.")

    fig, (left, right) = axes()
    trajectory = episode["series"][PRIMARY]
    trajectory_times = np.array([row["time"] for row in trajectory])
    lower = np.array([row["proper_min"] for row in trajectory])
    upper = np.array([row["proper_max"] for row in trajectory])
    left.fill_between(trajectory_times, lower, upper, color=COLORS[0], alpha=.18)
    left.plot(trajectory_times, lower, color=COLORS[0], lw=1.2, label="Minimum / maximum")
    left.plot(trajectory_times, upper, color=COLORS[0], lw=1.2)
    left.axhline(0, color="#6b7280", lw=.6)
    left.set(title="(a) Generated normal motion", ylabel="Normal radial velocity")
    left.legend(frameon=False, loc="best")
    metric = curvature["series"][PRIMARY]
    right.plot([row["time"] for row in metric], [row["weyl_C2_max"] for row in metric], color=COLORS[1], lw=1.6)
    right.set(title="(b) Actual metric curvature", ylabel=r"Maximum $C_{abcd}C^{abcd}$")
    local_window(left); local_window(right)
    save(fig, NAMES[1], "Normal motion and metric Weyl invariant", ("episode", "curvature"),
         "Recorded normal-velocity extrema and metric Weyl invariant maxima from actual projected ODE jets; no auxiliary-field substitution.")

    fig, (left, right) = axes()
    energy = np.array([row["windows"][0]["normal_shell"] for row in physical])
    boundary = np.array([row["windows"][0]["normal_boundary_flux"] for row in physical])
    pressure = np.array([row["windows"][0]["proper_pressure_work"] for row in physical])
    lapse = np.array([row["windows"][0]["momentum_lapse_work"] for row in physical])
    left.plot(times, energy - energy[0], color=COLORS[4], lw=1.6, label=r"$\Delta E_0$")
    left.plot(times, cumulative(times, boundary + pressure + lapse), color=COLORS[0], lw=1.25, ls="--", label="All balance terms")
    left.plot(times, cumulative(times, boundary), color=COLORS[1], lw=1.15, label="Boundary contribution")
    left.plot(times, cumulative(times, pressure + lapse), color=COLORS[2], lw=1.15, label="Pressure + lapse")
    left.set(title="(a) Normal-energy ledger", ylabel="Energy change / contribution")
    left.legend(frameon=False, loc="best", fontsize=6.8)
    norm = np.array([row["constraint_norm"] for row in trajectory])
    forcing = np.array([row["forcing_norm"] for row in trajectory])
    budget = norm[0] + cumulative(trajectory_times, forcing)
    right.plot(trajectory_times, norm, color=COLORS[4], lw=1.5, label=r"$R(T)$")
    right.plot(trajectory_times, budget, color=COLORS[0], lw=1.2, ls="--", label=r"$R(T_0)+\int \|f\|_2\,dT$")
    right.set(title="(b) Sampled constraint forcing", ylabel="Weighted constraint norm")
    right.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    right.legend(frameon=False, loc="best", fontsize=7)
    local_window(left); local_window(right)
    save(fig, NAMES[2], "Normal-energy and constraint-forcing ledgers", ("episode",),
         "Independent frame-trapezoid cumulative normal balance terms; all numerical-step samples for R and forcing. Sampled integrals are indicators, not continuum/time enclosures.")

    fig, (left, right) = axes()
    t = series["times"]
    cross = series["covariance_retained"] - series["covariance_without_cross"]
    left.plot(t, series["occupation_autonomous"][:, 0], color=COLORS[4], lw=1.7, label="Actual mode 0")
    left.plot(t, series["occupation_autonomous"][:, 1], color="#9ca3af", lw=.85, ls=":", label="Actual mode 1")
    left.plot(t, series["occupation_memory_off"][:, 0], color=COLORS[0], lw=1.25, label="Memory off")
    left.plot(t, series["occupation_drive_off"][:, 0], color=COLORS[1], lw=1.25, label="Initial drive off")
    left.plot(t, np.diagonal(series["covariance_without_cross"], axis1=1, axis2=2).real[:, 0], color=COLORS[2], lw=1.25, label="Initial cross removed")
    left.set(title="(a) Same fixed observer", ylabel="Occupation per mode")
    left.legend(frameon=False, loc="best", fontsize=6.7)
    errors = {
        "Reduced / full": np.max(np.abs(series["occupation_retained"] - series["occupation_full"]), axis=1),
        "Drive control": np.max(np.abs(series["occupation_drive_off"] - series["occupation_drive_off_full"]), axis=1),
        "Cross control": np.max(np.abs(np.diagonal(cross - series["cross_covariance_full"], axis1=1, axis2=2)), axis=1),
        "Memory control": np.max(np.abs(series["occupation_memory_off"] - memory_reference), axis=1),
    }
    for color, (label, error) in zip(COLORS, errors.items()):
        positive = (error > 0) & (t > t[0])
        right.semilogy(t[positive], error[positive], color=color, lw=1.25, label=label)
    right.set(title="(b) Each path has its own reference", ylabel="Maximum occupation error")
    right.legend(frameon=False, loc="best", fontsize=6.7)
    save(fig, NAMES[3], "Fixed-observer response and independent control errors", ("local", "local_npz", "memory", "memory_npz"),
         "Left: actual autonomous occupations of modes 0/1, with mode-0 causal omission curves. Right: absolute max errors over both retained modes, using independent full drive/cross and triangular memory references. Initial roundoff sample and exact zeros omitted from logarithmic panel.")
    manifest = {"schema": "NSC-FINITE-REGENERATION-FIGURES-v1", "science_commit": science_commit,
        "matplotlib": matplotlib.__version__, "numpy": np.__version__, "font": "DejaVu Sans",
        "figure_inches": [7.2, 2.65], "producer_sha256": digest(Path(__file__).read_bytes()), "figures": outputs}
    path = ROOT / "paper/finite-regeneration/figure-manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--science-commit", required=True)
    args = parser.parse_args()
    result = generate(args.science_commit)
    print(json.dumps({name: row["sha256"] for name, row in result["figures"].items()}, indent=2))
