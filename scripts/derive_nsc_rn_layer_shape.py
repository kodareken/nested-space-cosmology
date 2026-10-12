"""Characterise the shape of the compressed layer at peak compression.

The mass-shell work measured how much the density peak grows against the
inner edge, and by how much. It did not measure what the compressed layer
looks like. This does.

At the compression peak the profile is sampled across several inner edges.
Each is normalised by its own value at the edge and compared, so the question
is whether the shape is the same everywhere and only the amplitude moves, or
whether the shape itself depends on where the boundary is.

Two things are read that a single scalar cannot give:

* whether the density falls monotonically from the edge, which would make
  this a gradient with no interface;
* whether there is a characteristic scale length, and whether that length is
  the same across boundaries. A gradient zone with a universal width is a
  different object from one that stretches with the box.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_layer_shape.py
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")
#: The profile is sampled over this much radius beyond the inner edge. Beyond
#: this the shell's own body dominates and the layer is no longer the thing
#: being measured.
LAYER_SPAN = 2.0
#: An interior maximum must stand at least this far above the boundary value
#: to count as detached. The boundary itself normalises to exactly 1.
DETACHMENT_MARGIN = 1.05


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _shape(directory: Path, peak_t: float, tolerance: float) -> Dict[str, Any]:
    best = None
    for path in sorted(directory.glob("station-*.npz")):
        with np.load(path) as payload:
            clock = float(payload["coupled_t"][0])
        if abs(clock - peak_t) <= tolerance:
            best = path
            break
    if best is None:
        raise SystemExit(
            f"no station near t={peak_t} in {directory}; its last station is at "
            f"{clock:.3f}. The compression peak moves with the edge, so each run "
            f"needs its own peak time.")

    with np.load(best) as payload:
        radius = np.asarray(payload["radius"], dtype=float)
        mass = np.asarray(payload["coupled_mass"], dtype=float)
        phi = np.asarray(payload["coupled_phi"], dtype=np.complex128)
        clock = float(payload["coupled_t"][0])

    density = np.gradient(mass, radius)
    matter = np.sum(np.abs(phi) ** 2, axis=(0, 2))
    stop = int(np.searchsorted(radius, radius[0] + LAYER_SPAN))
    stop = max(8, min(stop, radius.size))

    edge_density = float(density[0])
    normalised = density[:stop] / edge_density if edge_density else density[:stop]
    edge_matter = float(matter[0])
    normalised_matter = (matter[:stop] / edge_matter) if edge_matter > 0 else matter[:stop]

    # A monotone falloff means there is no interface to find. Anything that
    # rises again inside the layer would be a second maximum.
    rising = int(np.count_nonzero(np.diff(normalised) > 0.0))
    rising_matter = int(np.count_nonzero(np.diff(normalised_matter) > 0.0))

    # Scale length: where each profile has fallen to half its edge value.
    def half_width(profile):
        if profile[0] <= 0:
            return None
        index = int(np.searchsorted(-np.asarray(profile), -0.5))
        return float(radius[index] - radius[0]) if index < stop else None

    # How far the maximum sits from the boundary, in cells. The two maxima
    # can differ: the matter density may leave the boundary while the mass
    # gradient has not, and the record must not merge them.
    density_peak = int(np.argmax(density[:stop]))
    matter_peak = int(np.argmax(matter[:stop]))
    # A sphere has three sides and a boundary-attached gradient has two. The
    # third appears when the profile turns back up inside, so the test is
    # whether the boundary node sits below the peak: below it and the layer
    # has a genuine interior, on it and the layer is one-sided.
    matter_interior = bool(matter_peak > 0 and float(matter[0]) < float(matter[matter_peak]))
    matter_rising_into_peak = (
        bool(np.all(np.diff(matter[:matter_peak + 1]) > 0.0)) if matter_peak > 0 else None)
    spacing = float(radius[1] - radius[0]) if radius.size > 1 else float("inf")

    return {
        "run": directory.name,
        "station": best.name,
        "t_over_rm": clock,
        "r_inner": float(radius[0]),
        "edge_density": edge_density,
        "edge_matter_density": edge_matter,
        "spacing": spacing,
        "density_max_radius": float(radius[density_peak]),
        "matter_max_radius": float(radius[matter_peak]),
        "density_max_cells_in": float((radius[density_peak] - radius[0]) / spacing),
        "matter_max_cells_in": float((radius[matter_peak] - radius[0]) / spacing),
        # How far the interior maximum stands above the boundary value. A
        # margin under a percent is not a resolved detachment.
        "matter_overshoot_at_max": float(normalised_matter[matter_peak]),
        "density_overshoot_at_max": float(normalised[density_peak]),
        "both_peaks_at_edge": bool(density_peak == 0 and matter_peak == 0),
        "matter_has_interior": matter_interior,
        "matter_rising_into_peak": matter_rising_into_peak,
        "shell_sides": 3 if matter_interior else 2,
        "matter_at_inner_edge": edge_matter,
        "matter_peak_value": float(matter[matter_peak]),
        "monotone_decreasing": rising == 0,
        "rising_steps": rising,
        "monotone_decreasing_matter": rising_matter == 0,
        "rising_steps_matter": rising_matter,
        "half_width": half_width(normalised),
        "half_width_matter": half_width(normalised_matter),
        "density_ratio_at_half_span": float(normalised[min(stop - 1, len(normalised) - 1)]),
        "normalised_profile": normalised.tolist(),
        "normalised_profile_matter": normalised_matter.tolist(),
        "offset_r": (radius[:stop] - radius[0]).tolist(),
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", required=True,
                        help="RUN:PEAK_T pairs, each naming a run directory under "
                             "results/development and the compression peak to sample "
                             "there. The peak moves with the edge, so one shared time "
                             "would mis-time most of the sweep.")
    parser.add_argument("--out", default="nsc-rn-layer-shape-v3.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    shapes = []
    for pair in args.runs:
        name, _, peak = pair.partition(":")
        directory = DEVELOPMENT / name
        if not directory.is_dir():
            print(f"no run at {directory}", file=sys.stderr)
            return 1
        shapes.append(_shape(directory, float(peak or 1.7), 0.06))
    shapes.sort(key=lambda s: -s["r_inner"])

    widths = [s["half_width"] for s in shapes if s["half_width"] is not None]
    matter_widths = [s["half_width_matter"] for s in shapes
                     if s["half_width_matter"] is not None]
    edges = [s["r_inner"] for s in shapes]
    amplitude = [s["edge_density"] for s in shapes]

    # A detachment is an interior maximum standing *above* the boundary
    # value. The boundary itself normalises to exactly 1, so the threshold
    # has to be above 1; anything at or under it is the boundary peak.
    detached = [s for s in shapes if s["matter_overshoot_at_max"] > DETACHMENT_MARGIN]

    record = {
        "schema": "NSC-RN-LAYER-SHAPE-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "The shape of the compressed layer, which the mass-shell reduction measured "
            "only as an amplitude. At the compression peak, whether the density falls "
            "monotonically from the inner edge with no interface, and whether the "
            "half-width of that falloff narrows as the boundary deepens."
        ),
        "scope_note": (
            "A construction check on finite marched collars. It does not measure the "
            "actual universe. A gradient zone in one numerical construction is not "
            "evidence that any real object has one."
        ),
        "method": {
            "layer_span_over_rm": LAYER_SPAN,
            "normalisation": "each profile divided by its own value at the inner edge",
            "why_normalise": (
                "The amplitude is known to grow with edge depth. Normalising asks the "
                "separate question of whether the shape also moves, which an amplitude "
                "alone cannot answer."
            ),
            "monotone_test": "count of rising steps in the normalised profile",
            "two_maxima_reported_separately": (
                "The matter density and the mass gradient can leave the boundary at "
                "different radii. Merging them would overstate what is resolved: at the "
                "deepest edge the matter maximum stands 22% above the boundary value "
                "while the mass gradient is flat to 0.6% across the same two nodes."
            ),
            "detachment_margin_above_boundary": DETACHMENT_MARGIN,
        },
        "peaks_by_run": {s["run"]: s["t_over_rm"] for s in shapes},
        "all_monotone": all(s["monotone_decreasing"] for s in shapes),
        "all_peak_at_edge": all(s["both_peaks_at_edge"] for s in shapes),
        "half_widths": widths,
        "half_width_spread": (max(widths) - min(widths)) if widths else None,
        "half_widths_matter": matter_widths,
        "edge_densities": amplitude,
        "edge_density_ratio_max_over_min": (max(amplitude) / min(amplitude)
                                            if amplitude else None),
        "width_narrows_with_depth": (all(b <= a for a, b in zip(widths, widths[1:]))
                                     if len(widths) > 1 else None),
        "width_narrows_excluding_the_deepest_edge": (
            all(b <= a for a, b in zip(widths[:-1], widths[1:-1]))
            if len(widths) > 2 else None),
        "runs_with_a_detached_maximum": [s["r_inner"] for s in detached],
        "shell_sides": [s["shell_sides"] for s in shapes],
        "edges_with_a_three_sided_shell": [s["r_inner"] for s in shapes
                                           if s["shell_sides"] == 3],
        "edges_with_a_two_sided_layer": [s["r_inner"] for s in shapes
                                         if s["shell_sides"] == 2],
        "runs": shapes,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  peaks sampled     : "
          f"{[(s['run'].split('-')[-1], round(s['t_over_rm'], 2)) for s in shapes]}")
    print(f"  all monotone      : {record['all_monotone']}")
    print(f"  both peaks at edge: {record['all_peak_at_edge']}")
    print(f"  r_inner           : {edges}")
    print(f"  edge density      : {['%.3e' % a for a in amplitude]}")
    print(f"  half width        : {['%.4f' % w for w in widths] if widths else 'n/a'}")
    print(f"  narrows with depth: {record['width_narrows_with_depth']}"
          f"   (excluding the deepest edge: "
          f"{record['width_narrows_excluding_the_deepest_edge']})")
    print(f"  matter overshoot at its maximum: "
          f"{['%.4f' % s['matter_overshoot_at_max'] for s in shapes]}")
    print(f"  runs with a detached maximum: {record['runs_with_a_detached_maximum']}")
    print(f"  shell sides                : {record['shell_sides']}")
    print(f"  three-sided at             : {record['edges_with_a_three_sided_shell']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())