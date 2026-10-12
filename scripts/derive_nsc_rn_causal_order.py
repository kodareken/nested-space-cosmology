"""Does the causal ordering run one way, nest to nest, and never back?

A one-way flow is a statement about causal structure, not about a profile.
This asks it directly, from the recorded states, using the characteristic
speeds in the evolved metric.

The test has three parts, and all three have to hold for the claim to stand:

* **Sealed.** The outgoing speed is strictly negative inside the innermost
  zero of it, so nothing inside can get out.
* **Stable.** The number of zeros does not change over time. A second zero
  appearing outside would be a new sealed region forming, and losing one
  would be a sealed region opening.
* **One-way in time.** The located zero does not reverse. A horizon that
  moves outward and then moves back is not a one-way membrane.

None of this is a global result and none of it is about a sequence of
nests. The chart in this construction defines one child interval and one
parent, so "nest to nest to nest" can be tested here only as a single
step. Whether the nesting iterates is a property of the construction, not
something these states can answer.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_causal_order.py
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
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_rn_observables as obs  # noqa: E402
from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
import derive_nsc_rn_parent as driver  # noqa: E402

DEVELOPMENT = LAB / "results/development"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _station_rows(directory: Path) -> List[Dict[str, Any]]:
    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)),
                      key=lambda p: float(STATION.search(p.name).group(1)))
    rows: List[Dict[str, Any]] = []
    for station in stations:
        with np.load(station) as payload:
            radius = np.asarray(payload["radius"], dtype=float)
            lapse = np.asarray(payload["frozen_lapse"], dtype=float)
            shift = np.asarray(payload["frozen_shift"], dtype=float)
            clock = float(payload["coupled_t"][0])
        located = obs.trapping_horizon(radius, lapse, shift)
        outgoing = np.asarray(located["outgoing_speed"], dtype=float)
        ingoing = np.asarray(located["ingoing_speed"], dtype=float)
        zeros = int(np.count_nonzero(np.diff(np.sign(outgoing)) != 0))
        rows.append({
            "t_over_rm": clock,
            "r_inner": float(radius[0]),
            "spacing": float(radius[1] - radius[0]),
            "trapping_horizon": located["trapping_horizon"],
            "zeros": zeros,
            "outgoing_at_inner_edge": float(outgoing[0]),
            "ingoing_at_inner_edge": float(ingoing[0]),
            "inner_edge_sealed": bool(outgoing[0] < 0.0),
            "ingoing_non_zero": bool(abs(ingoing[0]) > 0.0),
        })
    return rows


#: A reversal smaller than this fraction of a grid cell is below what the
#: sampled profile can resolve and is not counted against one-way flow. The
#: first version of this test used an absolute tolerance near 1e-9, which
#: flagged a quarter-cell excursion as a reversal.
RESOLUTION_FRACTION = 0.5


def _verdict(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not rows:
        return {"assessable": False}
    radii = [r["trapping_horizon"] for r in rows if r["trapping_horizon"] is not None]
    spacing = rows[0].get("spacing", 0.0)
    threshold = RESOLUTION_FRACTION * spacing if spacing else 0.0

    # Largest single step that moves the located zero inward, and the total
    # excursion, both in cells, so the reader can judge resolution directly.
    inward = [a - b for a, b in zip(radii, radii[1:]) if b < a]
    largest_inward = max(inward) if inward else 0.0
    excursion = max(radii) - min(radii) if radii else 0.0

    resolved_reversal = largest_inward > threshold
    return {
        "assessable": True,
        "grid_spacing": spacing,
        "reversal_threshold": threshold,
        "sealed_everywhere": all(r["inner_edge_sealed"] for r in rows),
        "zero_count_constant": len({r["zeros"] for r in rows}) == 1,
        "zero_count": rows[0]["zeros"],
        "largest_inward_step": largest_inward,
        "largest_inward_step_in_cells": (largest_inward / spacing if spacing else None),
        "total_excursion": excursion,
        "total_excursion_in_cells": (excursion / spacing if spacing else None),
        "excursion_resolved": bool(spacing and excursion > threshold),
        "resolved_reversal": bool(resolved_reversal),
        "horizon_moves_outward_only": not resolved_reversal,
        "horizon_drifts_outward": bool(radii and radii[-1] > radii[0]),
        "ingoing_never_vanishes": all(r["ingoing_non_zero"] for r in rows),
        "one_way": (all(r["inner_edge_sealed"] for r in rows)
                    and len({r["zeros"] for r in rows}) == 1
                    and not resolved_reversal),
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", required=True,
                        help="RUN:MU pairs, each naming a run directory under "
                             "results/development and the mass ratio it was built at")
    parser.add_argument("--out", default="nsc-rn-causal-order-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    runs = []
    for pair in args.runs:
        name, _, mu = pair.partition(":")
        directory = DEVELOPMENT / name
        if not directory.is_dir():
            print(f"no run at {directory}", file=sys.stderr)
            return 1
        rows = _station_rows(directory)
        geometry = driver.case_geometry(float(mu or 1.04))
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        runs.append({
            "run": name,
            "mu": float(mu or 1.04),
            "r_plus": float(model.horizons["r_plus"]),
            "r_minus": float(model.horizons["r_minus"]),
            "stations": len(rows),
            "verdict": _verdict(rows),
            "rows": rows,
        })

    one_way = [r["run"] for r in runs if r["verdict"].get("one_way")]
    record = {
        "schema": "NSC-RN-CAUSAL-ORDER-v1",
        "api_version": "nsc-rn-causal-order-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "question": (
            "Does the causal ordering run one way, nest to nest, and never back? Answered "
            "from the characteristic speeds in the evolved metric rather than from a "
            "profile, because one-way flow is a claim about causal structure."
        ),
        "the_three_tests": {
            "sealed": "the outgoing speed is strictly negative inside the innermost zero",
            "stable": "the number of zeros does not change over time",
            "one_way_in_time": "the located zero does not reverse",
        },
        "what_this_cannot_answer": (
            "Whether the nesting iterates. This construction defines one child interval "
            "and one parent, so nest-to-nest-to-nest can be tested here only as a single "
            "step. A hierarchy is a property of the construction, not of these states, and "
            "no claim about one is made."
        ),
        "scope_note": (
            "A construction check on recorded finite states. It does not measure the "
            "actual universe. A simple zero of the outgoing characteristic on a sampled "
            "slice is the standard numerical proxy for a trapping surface; a sequence of "
            "such slices is not a proof of a global event horizon."
        ),
        "runs_one_way": one_way,
        "runs_total": len(runs),
        "all_one_way": len(one_way) == len(runs),
        "runs": runs,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    for run in runs:
        v = run["verdict"]
        rows = run["rows"]
        radii = [r["trapping_horizon"] for r in rows if r["trapping_horizon"] is not None]
        print(f"  {run['run']:<44} zeros={v['zero_count']}  sealed={v['sealed_everywhere']}"
              f"  one-way={v['one_way']}")
        if radii:
            print(f"      horizon {radii[0]:.4f} -> {radii[-1]:.4f}"
                  f"   (r_plus = {run['r_plus']:.4f})")
            print(f"      excursion {v['total_excursion']:.3e} = "
                  f"{v['total_excursion_in_cells']:.3f} cells"
                  f"   resolved={v['excursion_resolved']}")
            print(f"      largest inward step {v['largest_inward_step']:.3e} = "
                  f"{v['largest_inward_step_in_cells']:.4f} cells"
                  f"   threshold {v['reversal_threshold']:.3e}")
    print(f"  all runs one-way: {record['all_one_way']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())