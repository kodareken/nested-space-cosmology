#!/usr/bin/env python3
"""Plot one saved NSC-P001 three-grid endpoint comparison.

Reads a completed comparison and its three receipts. It does not evolve the
system, rebuild a history, or overwrite an existing image.

    python scripts/run_python.py scripts/plot_nsc_atlas_fixed_endpoint.py --comparison results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json --output docs/figures/nsc-atlas-fixed-endpoint-comparison.png
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from analyze_nsc_atlas_fixed_endpoint import compare
from check_prediction_register import (
    ENVELOPE_STATUS, UNCONVERGED_STATUS, assert_comparison_matches,
    endpoint_receipts, finite_amplitude_block, load_comparison, local_file,
)
from recursive_horizons.nsc_atlas_endpoint_comparison import (
    ENDPOINT_REFERENCE_RADIUS, ENDPOINT_SOURCE_X, FIDUCIAL_LENGTH_METRES,
)
from register_nsc_atlas_fixed_endpoint import endpoint_rows

ROOT = Path(__file__).resolve().parents[1]
PRODUCER = "scripts/plot_nsc_atlas_fixed_endpoint.py"
MESHES = ((2, 32), (4, 64), (8, 128))
RESIDUALS = (
    ("mass_trace", "A", "Cap-trace mass residual"),
    ("mass_matched", "B", "Matched-wall mass residual"),
    ("wall_flux", "C", "Conformal wall-flux residual"),
)
STATUS_TEXT = {
    UNCONVERGED_STATUS: "unconverged numerical estimate",
    ENVELOPE_STATUS: "finite-amplitude numerical estimate with a resolution envelope",
}
INCREASE = "#9b2c2c"
UNCHANGED = "#8a5a12"
INK = "#17283a"
MUTED = "#46525e"


def yes_no(flag):
    if type(flag) is not bool:
        raise ValueError("contraction flag must be a boolean")
    return "yes" if flag else "no"


def as_real(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(label + " must be a finite real")
    return float(value)


def format_saved(value):
    value = as_real(value, "saved number")
    if value == 0:
        return "0"
    return json.dumps(value, allow_nan=False)


def compact_real(value):
    """Readable metres. Binary noise past eight digits is not a new digit of radius."""
    value = as_real(value, "saved number")
    if value == 0:
        return "0"
    return format(value, ".8g")


def distinct_labels(values):
    """Short labels that still separate unequal saved numbers."""
    full = [format_saved(value) for value in values]
    for digits in range(4, 17):
        short = ["0" if value == 0 else format(value, "." + str(digits) + "g") for value in values]
        if all((short[i] == short[j]) == (values[i] == values[j])
               for i in range(len(values)) for j in range(i + 1, len(values))):
            return short
    return full


def residual_scale(values):
    """Use a log axis only when every residual is strictly positive.

    An exact zero has no logarithm. It stays on a linear axis.
    """
    if any(value == 0 for value in values):
        return "linear"
    if min(values) <= 0:
        return "linear"
    if max(values) / min(values) >= 50:
        return "log"
    return "linear"


def hidden_increase(values, scale):
    """Return the successive increase that a shared axis would flatten."""
    span = max(values) - min(values)
    for index, (left, right) in enumerate(zip(values, values[1:])):
        if right <= left:
            continue
        if scale == "log":
            hidden = right / left < 1.25
        else:
            hidden = span == 0 or (right - left) < 0.12 * span
        if hidden:
            return index
    return None


def linear_limits(values):
    low, high = min(values), max(values)
    if low == high:
        if low == 0:
            return -0.2, 0.55
        pad = abs(low) * 0.3
        return low - pad, high + pad
    pad = 0.28 * (high - low)
    bottom = -0.08 * (high - low) if low == 0 else low - pad
    return bottom, high + pad


def log_limits(values):
    low, high = min(values), max(values)
    if low <= 0 or high <= 0:
        raise RuntimeError("log axis requested for a nonpositive residual")
    return low / (10 ** 0.4), high * (10 ** 0.4)


def prepare_output(path):
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    if candidate.suffix.lower() != ".png":
        raise ValueError("output must be a PNG path")
    if not candidate.parent.is_dir():
        raise ValueError("output parent directory must already exist")
    resolved = candidate.resolve()
    if not resolved.is_relative_to(ROOT.resolve()):
        raise ValueError("output must stay inside the repository")
    if candidate.exists() or resolved.exists():
        raise FileExistsError("output already exists")
    return resolved


def comparison_binding(path):
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    resolved = candidate.resolve()
    if not resolved.is_file() or not resolved.is_relative_to(ROOT.resolve()):
        raise ValueError("comparison must be an existing repository file")
    relative = resolved.relative_to(ROOT.resolve()).as_posix()
    raw = local_file(ROOT, relative).read_bytes()
    return dict(path=relative, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def require_unchanged(binding):
    raw = local_file(ROOT, binding["path"]).read_bytes()
    if len(raw) != binding["bytes"] or hashlib.sha256(raw).hexdigest() != binding["sha256"]:
        raise RuntimeError("saved input changed during rendering: " + binding["path"])


def mesh_labels(records):
    found = [(record["geometry"]["cap_cells"], record["geometry"]["exterior_cells"]) for record in records]
    if found != list(MESHES):
        raise ValueError("unexpected mesh family")
    return [f"{cap} cap\n{exterior} exterior" for cap, exterior in found], [f"{cap}/{exterior}" for cap, exterior in found]


def authenticate(comparison_path):
    binding = comparison_binding(comparison_path)
    saved, found = load_comparison(ROOT, binding)
    rows = endpoint_rows(found)
    records = endpoint_receipts(found, rows)
    recomputed = compare(records)
    assert_comparison_matches(saved, recomputed)
    # The same scale, event, and envelope claims required before register binding.
    finite_amplitude_block(saved, binding, rows)
    if saved["prediction_status"] not in STATUS_TEXT:
        raise ValueError("unknown comparison status")
    if saved["physical_evolutions"] != 0:
        raise ValueError("endpoint figure cannot record a new evolution")
    return binding, saved, rows, records


def draw_values(ax, values, labels, scale, color):
    xs = (0, 1, 2)
    hidden = hidden_increase(values, scale)
    for index, (left, right) in enumerate(zip(values, values[1:])):
        if right > left:
            style, segment, width = "--", INCREASE, 2.3
        elif right == left:
            style, segment, width = "--", UNCHANGED, 1.8
        else:
            style, segment, width = "-", color, 1.8
        ax.plot(xs[index:index + 2], (left, right), style, color=segment, lw=width, zorder=2)
        if right > left and index != hidden:
            y = math.sqrt(left * right) if scale == "log" else (left + right) / 2
            ax.annotate("increases", ((xs[index] + xs[index + 1]) / 2, y),
                        textcoords="offset points", xytext=(0, 8), ha="center", va="bottom",
                        fontsize=8, fontweight="bold", color=INCREASE, clip_on=False, zorder=5)
        elif right == left:
            ax.annotate("unchanged", ((xs[index] + xs[index + 1]) / 2, left),
                        textcoords="offset points", xytext=(0, 8), ha="center", va="bottom",
                        fontsize=8, color=UNCHANGED, clip_on=False, zorder=5)
    ax.plot(xs, values, "o", color=color, ms=7, mew=0.7, mec="white", zorder=3)
    top = max(values)
    for x, value, label in zip(xs, values, labels):
        above = value < top or top == 0
        ax.annotate(label, (x, value), textcoords="offset points",
                    xytext=(0, 8 if above else -11), ha="center",
                    va="bottom" if above else "top", fontsize=7.5, color=INK, zorder=4)
    return hidden


def add_increase_inset(ax, index, values, short_meshes, labels):
    left, right = values[index], values[index + 1]
    anchor = (0.52, 0.34, 0.45, 0.5) if index else (0.04, 0.34, 0.45, 0.5)
    inset = ax.inset_axes(anchor)
    inset.plot((0, 1), (left, right), "--", color=INCREASE, lw=1.8)
    inset.plot((0, 1), (left, right), "o", color=INCREASE, ms=5)
    span = right - left
    pad = span * 0.55 if span else 1
    bottom = -0.12 * (right if right else 1) if left == 0 else left - pad
    inset.set_ylim(bottom, right + pad)
    inset.set_xlim(-0.4, 1.4)
    inset.set_xticks((0, 1), (short_meshes[index], short_meshes[index + 1]))
    inset.tick_params(labelsize=6.5, colors=INCREASE)
    inset.set_title("increases", fontsize=8, color=INCREASE, fontweight="bold", pad=2)
    for x, value, label in ((0, left, labels[index]), (1, right, labels[index + 1])):
        inset.annotate(label, (x, value), textcoords="offset points", xytext=(0, 4),
                       ha="center", fontsize=6.5, color=INCREASE)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 3))
    inset.yaxis.set_major_formatter(formatter)
    inset.grid(alpha=.15, linewidth=.5)
    for spine in inset.spines.values():
        spine.set_color(INCREASE)
        spine.set_linewidth(.7)
    inset.set_facecolor("white")


def draw_residual(ax, letter, title, row, ticks, short_meshes, color):
    values = [as_real(value, title) for value in row["values"]]
    if len(values) != 3 or any(value < 0 for value in values):
        raise ValueError(title + " must contain three nonnegative maxima")
    strictly = row["strictly_decreasing"]
    scale = residual_scale(values)
    if scale == "log":
        ax.set_yscale("log")
        ax.set_ylim(*log_limits(values))
    else:
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 3))
        ax.yaxis.set_major_formatter(formatter)
        ax.set_ylim(*linear_limits(values))
    labels = distinct_labels(values)
    hidden = draw_values(ax, values, labels, scale, color)
    if hidden is not None:
        add_increase_inset(ax, hidden, values, short_meshes, labels)
    ax.set_xlim(-0.45, 2.45)
    ax.set_xticks((0, 1, 2), ticks)
    ax.set_title(f"{letter}  {title}", loc="left", color=INK, pad=8)
    flag = "strictly decreasing: " + yes_no(strictly)
    on_right = values[0] >= values[1] and values[0] >= values[2]
    ax.text(0.98 if on_right else 0.02, 0.96, flag, transform=ax.transAxes,
            ha="right" if on_right else "left", va="top", fontsize=8,
            color=INK if strictly else INCREASE)
    if any(value == 0 for value in values):
        note, note_color = "linear axis (exact zero)", MUTED
    elif scale == "log":
        note, note_color = "logarithmic axis", MUTED
    else:
        note, note_color = None, MUTED
    if note:
        # Keep the scale note off the inset that exposes a flattened increase.
        left_note = hidden == 1
        ax.text(0.02 if left_note else 0.98, 0.03, note, transform=ax.transAxes,
                ha="left" if left_note else "right", va="bottom", fontsize=7.5, color=note_color)
    ax.set_ylabel("maximum |residual|")
    ax.yaxis.offsetText.set_fontsize(8)


def draw_radius(ax, observable, ticks, contracts, paired):
    grid = [as_real(value, "Jordan radius") for value in observable["grid_values_metres"]]
    if len(grid) != 3 or any(value <= 0 for value in grid):
        raise ValueError("three positive Jordan radii required")
    fine = as_real(observable["value_metres"], "fine Jordan radius")
    half = as_real(observable["resolution_half_width_metres"], "resolution half-width")
    low, high = (as_real(value, "resolution range") for value in observable["resolution_range_metres"])
    if half < 0 or low > high:
        raise ValueError("resolution envelope is not ordered")
    quiet = as_real(observable["quiet_reference_metres"], "quiet reference")
    view_low = min(low, min(grid))
    view_high = max(high, max(grid))
    span = view_high - view_low
    # A zero spread still needs a readable metre axis. The pad is display room,
    # not an enlargement of the saved envelope.
    pad = 1.0 if span == 0 else 0.62 * span
    y0, y1 = view_low - pad, view_high + pad
    ax.axhspan(low, high, color="#167d8d", alpha=0.18, zorder=0)
    ax.axhline(fine, color="#273d75", lw=1.15, zorder=1)
    quiet_inside = y0 <= quiet <= y1
    if quiet_inside:
        ax.axhline(quiet, color="#8a7560", ls=":", lw=1.2, zorder=1)
    ax.plot((0, 1, 2), grid, "o", color="#273d75", ms=7, mew=0.7, mec="white", zorder=3)
    for x, value, label in zip((0, 1, 2), grid, distinct_labels(grid)):
        above = value < max(grid) or span == 0
        ax.annotate(label + " m", (x, value), textcoords="offset points",
                    xytext=(0, 8 if above else -11), ha="center",
                    va="bottom" if above else "top", fontsize=7.5, color=INK, zorder=4)
    ax.set_xlim(-0.45, 2.45)
    ax.set_ylim(y0, y1)
    ax.set_xticks((0, 1, 2), ticks)
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_title("D  Jordan wall areal radius", loc="left", color=INK, pad=8)
    ax.set_ylabel(r"$R_J$ (m)")
    if len(paired) != 2:
        raise ValueError("two radius pair differences required")
    saved_left, saved_right = (as_real(value, "radius pair difference") for value in paired)
    # Convert the authenticated dimensionless differences directly. Subtracting
    # already rounded metre values adds cancellation and can change an ulp-sized
    # contraction decision without changing the underlying saved comparison.
    ell = as_real(observable["scale"]["ell_metres"], "physical length scale")
    gaps_m = [ell * saved_left, ell * saved_right]
    if saved_right > saved_left:
        movement = "successive $|\\Delta R_J|$ increases"
    elif saved_right < saved_left:
        movement = "successive $|\\Delta R_J|$ decreases"
    else:
        movement = "successive $|\\Delta R_J|$ unchanged"
    gap_labels = distinct_labels(gaps_m)
    lines = ["pair-gap contracts: " + yes_no(contracts),
             movement + ": " + gap_labels[0] + " m → " + gap_labels[1] + " m",
             "half-width " + compact_real(half) + " m; horizontal line is the fine-grid value"]
    if span == 0:
        lines.append("all three grid values coincide")
    if quiet_inside:
        lines.append("dotted line: quiet reference " + compact_real(quiet) + " m")
    else:
        lines.append("quiet reference " + compact_real(quiet) + " m is outside this window")
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes, ha="left", va="top",
            fontsize=7.4, color=INK if contracts else INCREASE, zorder=5,
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#e1e6ec", lw=0.6, alpha=0.94))


def render(comparison_path, output_path):
    output = prepare_output(output_path)
    binding, saved, rows, records = authenticate(comparison_path)
    ticks, short_meshes = mesh_labels(records)
    observable = saved["physical_observable"]
    status = saved["prediction_status"]
    mass_flux = saved["requested_mass_and_flux_contraction"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
                         "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(11.8, 8.6))
    try:
        grid = fig.add_gridspec(2, 2, left=.075, right=.985, top=.80, bottom=.22,
                                hspace=.48, wspace=.28)
        colors = ("#b36b24", "#167d8d", "#273d75")
        axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
                fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
        for ax, (name, letter, title), color in zip(axes, RESIDUALS, colors):
            draw_residual(ax, letter, title, saved["residual_maxima"][name], ticks, short_meshes, color)
        draw_radius(axes[3], observable, ticks, saved["endpoint_radius_pair_difference_contracts"],
                    saved["endpoint_radius_pair_differences"])
        for ax in axes:
            ax.tick_params(labelsize=8)
        for ax in axes:
            ax.grid(alpha=.18, linewidth=.7)
            ax.set_axisbelow(True)
        fig.suptitle("NSC-P001 fixed incoming-label endpoint", fontsize=16, fontweight="bold",
                     y=.975, color=INK)
        fig.text(.5, .925, r"Declared fiducial $\ell$ = 10 km  ·  fixed incoming-label endpoint"
                          r"  ·  resolution envelope, not a total error",
                 ha="center", fontsize=10.5, color="#4b5968")
        status_color = INCREASE if status == UNCONVERGED_STATUS else MUTED
        fig.text(.5, .885, "Saved status: " + STATUS_TEXT[status]
                          + ".  Requested mass and flux contraction: " + yes_no(mass_flux) + ".",
                 ha="center", fontsize=10, color=status_color)
        alpha = as_real(saved["shared_alpha"], "shared alpha")
        fig.text(.075, .128,
                 rf"Declared fiducial $\ell$ = {FIDUCIAL_LENGTH_METRES / 1000:g} km. "
                 rf"Fixed incoming-label endpoint $x$ = {format_saved(ENDPOINT_SOURCE_X)}, "
                 rf"$V$ = 1, reference radius {format_saved(ENDPOINT_REFERENCE_RADIUS)} $\ell$. "
                 rf"$\alpha$ = {format_saved(alpha)}.",
                 ha="left", fontsize=8.4, color=MUTED)
        fig.text(.075, .086,
                 "The radius band covers the three saved endpoints around the fine-grid value. "
                 "It is an observed resolution envelope, not a total error, confidence interval, or observational test.",
                 ha="left", fontsize=8.4, color=MUTED)
        fig.text(.075, .044, saved["numerical_method_note"], ha="left", fontsize=8.4, color=MUTED,
                 parse_math=False)
        stream = io.BytesIO()
        fig.savefig(stream, dpi=180, facecolor="white", metadata={"Software": PRODUCER})
    finally:
        plt.close(fig)
    image = stream.getvalue()
    if not image.startswith(b"\x89PNG\r\n\x1a\n"):
        raise RuntimeError("renderer did not produce a PNG")
    require_unchanged(binding)
    for row in rows:
        require_unchanged(row)
    if output.exists():
        raise FileExistsError("output already exists")
    with output.open("xb") as target:
        target.write(image)
    return dict(output=output.relative_to(ROOT.resolve()).as_posix(), bytes=len(image),
                sha256=hashlib.sha256(image).hexdigest(), prediction_status=status,
                requested_mass_and_flux_contraction=mass_flux,
                endpoint_radius_pair_difference_contracts=saved["endpoint_radius_pair_difference_contracts"],
                residual_strictly_decreasing={name: saved["residual_maxima"][name]["strictly_decreasing"]
                                              for name, _letter, _title in RESIDUALS},
                physical_evolutions=0, time_evolutions=0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(render(args.comparison, args.output), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
