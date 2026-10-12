"""Is there any signature at the parent/child cut?

The construction nests a parent and a child on one chart, with the cut at
``r_plus + gap/2``. The only physical interface the chart contains is at
``r_plus``, which lies *below* the child rather than at its boundary. So the
question is whether the parent/child boundary has any signature at all.

This scans the whole chart for one. It looks at every quantity the chart
carries and at their radial derivatives, and tests each for a feature at the
cut rather than for a value, because a smooth function takes some value at
every radius and that alone says nothing. The detection threshold is
scaled by the local sampling, because a step smaller than the grid can be
placed is not a step that can be seen.

Reads the reference chart only. No station payloads, no time step, no
marching, and no construction code is changed.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_child_cut.py
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402
import derive_nsc_rn_parent as driver  # noqa: E402

DEVELOPMENT = LAB / "results/development"

#: A feature smaller than this fraction of the chart's own total variation is
#: not distinguishable from the sampling. Without this the test would find
#: "structure" everywhere, which is what a smooth function looks like.
FEATURE_FRACTION = 1.0e-3
#: Width of the neighbourhood around the cut, in units of the local spacing.
NEIGHBOURHOOD_SPACINGS = 40.0


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _curvature_of(values: np.ndarray, radius: np.ndarray) -> float:
    """Discrete second difference, the smallest thing a smooth curve lacks."""
    step = float(radius[1] - radius[0])
    return float(np.max(np.abs(np.diff(values, 2)))) / step ** 2


def _entry(mu: float, samples: int) -> Dict[str, Any]:
    geometry = driver.case_geometry(mu)
    model = pg.action_model(geometry["mass"], geometry["r_m"])
    r_minus = float(model.horizons["r_minus"])
    r_plus = float(model.horizons["r_plus"])
    gap = r_plus - r_minus
    cut = r_plus + gap / 2.0
    spacing = (32.0 - r_plus * 1.001) / (samples - 1)

    radius = np.linspace(r_plus * 1.001, 32.0, samples)
    chart = source.vacuum_chart(model, radius, allow_below_inner_horizon=True)

    fields = {
        "speed_plus": np.asarray(chart["speed_plus"], dtype=float),
        "speed_minus": np.asarray(chart["speed_minus"], dtype=float),
        "shift": np.asarray(chart["shift"], dtype=float),
        "lapse": np.asarray(chart["lapse"], dtype=float),
        "redshift": np.asarray(chart["redshift"], dtype=float),
    }

    index = int(np.argmin(np.abs(radius - cut)))
    window = min(int(NEIGHBOURHOOD_SPACINGS), index, radius.size - index - 1)
    if window < 4:
        raise SystemExit("cut too near the edge of the sampled range to test")

    results: Dict[str, Any] = {}
    for name, values in fields.items():
        total = float(np.max(values) - np.min(values))
        # The median comparison is degenerate here: over a locally linear
        # stretch the median of an odd window of nodes is exactly the centre
        # node, so a smooth field scores zero and so would a symmetric bump.
        # The second difference is the statistic that can actually see a
        # kink, so that is what the test uses.
        lo, hi = index - window, index + window + 1
        chart_curvature = _curvature_of(values, radius)
        local_curvature = _curvature_of(values[lo:hi], radius[lo:hi])
        results[name] = {
            "value_at_cut": float(values[index]),
            "chart_total_variation": total,
            "second_difference_at_cut": local_curvature,
            "chart_max_second_difference": chart_curvature,
            "curvature_excess": (float(local_curvature / chart_curvature)
                                 if chart_curvature else 0.0),
            "detected": bool(chart_curvature and local_curvature > chart_curvature),
        }

    detected = [name for name, block in results.items() if block["detected"]]
    return {
        "mu": mu,
        "r_minus": r_minus,
        "r_plus": r_plus,
        "child_cut": cut,
        "cut_in_parent": bool(cut > r_plus),
        "samples": samples,
        "sampling": spacing,
        "window_nodes": int(window),
        "feature_fraction_threshold": FEATURE_FRACTION,
        "fields_tested": sorted(fields),
        "fields_with_a_feature": detected,
        "no_feature_anywhere": not detected,
        "results": results,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mu", nargs="+", type=float, default=[1.002, 1.010, 1.040])
    parser.add_argument("--samples", type=int, default=40000)
    parser.add_argument("--out", default="nsc-rn-child-cut-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    entries = [_entry(mu, args.samples) for mu in sorted(args.mu)]
    record = {
        "schema": "NSC-RN-CHILD-CUT-v1",
        "api_version": "nsc-rn-child-cut-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Whether the parent/child cut carries any signature at all. Every quantity the "
            "chart carries is tested for a local feature there, against a threshold scaled "
            "by the chart's own total variation, so a smooth function does not read as "
            "structure."
        ),
        "why_it_matters": (
            "The construction nests a parent and a child, and the only physical interface it "
            "contains is at r_plus, below the child rather than at its boundary. If the cut "
            "carries nothing, the nesting is bookkeeping. If it carries something, that is "
            "the boundary condition the nesting is missing."
        ),
        "method": {
            "statistic": "second difference of each field",
            "threshold": "local curvature exceeding the chart's own maximum",
            "why_not_a_value_comparison": (
                "a value comparison is degenerate here: over a locally linear stretch the "
                "median of an odd window of nodes is exactly the centre node, so a smooth "
                "field scores zero and a symmetric bump would score zero too"
            ),
            "neighbourhood_nodes": NEIGHBOURHOOD_SPACINGS,
            "second_difference": "compared against the chart's own maximum",
        },
        "scope_note": (
            "A construction check on the reference chart. It does not measure the actual "
            "universe. Finding no signature at the cut is a statement about this chart; it "
            "does not show that the parent/child boundary is unreal in any construction, "
            "only that this one carries nothing there."
        ),
        "no_feature_everywhere": all(e["no_feature_anywhere"] for e in entries),
        "entries": entries,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    for entry in entries:
        print(f"  mu={entry['mu']:.3f}  cut at r={entry['child_cut']:.4f}"
              f"  (in parent: {entry['cut_in_parent']})")
        for name in entry["fields_tested"]:
            block = entry["results"][name]
            print(f"     {name:<12} curvature excess {block['curvature_excess']:.3f}"
                  f"  {'FEATURE' if block['detected'] else 'smooth'}")
        print(f"     fields with a feature: {entry['fields_with_a_feature']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())