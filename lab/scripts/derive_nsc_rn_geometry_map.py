"""Map the geometry of the RN radial family across the mass ratio.

The ambition behind this line is a map rather than a discovery: take one
space, describe its structure completely, and see how that structure moves
as the one free parameter changes. This is that map for the family in hand.

Nothing here is quantum and nothing here is exotic. The geometry is a
charged spherical solution in areal coordinates, the quantities are the two
horizons, the two nesting cuts the chart defines, the shift, and the causal
membrane located from the metric. All of it is classical.

The map answers a question the individual reductions cannot: as the mass
ratio grows, do the features keep their order, and where does anything
cross anything else. Ordering is what a map is for. If the child cut ever
crossed the horizon, or the horizon crossed the inner edge, the nesting
picture would have a transition in it, and that would be worth knowing.

Reads the reference chart only. No station payloads, no time step, no
marching.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_geometry_map.py
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


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _entry(mu: float) -> Dict[str, Any]:
    geometry = driver.case_geometry(mu)
    model = pg.action_model(geometry["mass"], geometry["r_m"])
    r_minus = float(model.horizons["r_minus"])
    r_plus = float(model.horizons["r_plus"])
    gap = r_plus - r_minus
    child_start = r_plus + gap / 16.0
    child_end = r_plus + gap / 2.0
    excision = 0.5 * (r_minus + r_plus)

    # Where the shift turns over, and where the outgoing speed vanishes.
    grid = np.linspace(r_minus + 1.0e-4, 40.0, 8000)
    chart = source.vacuum_chart(model, grid)
    outgoing = np.asarray(chart["speed_plus"], dtype=float)
    shift = np.asarray(chart["shift"], dtype=float)
    turn = float(grid[int(np.argmax(shift))])
    zero = None
    for i in range(1, grid.size):
        if outgoing[i - 1] < 0.0 <= outgoing[i]:
            span = outgoing[i] - outgoing[i - 1]
            zero = float(grid[i - 1]
                         + (-outgoing[i - 1]) / span * (grid[i] - grid[i - 1]))
            break

    # The ordering a map exists to show. The located zero is an interpolated
    # crossing on a finite grid, so it is known only to within a cell and can
    # sit below the horizon it coincides with. The tolerance is therefore
    # the local grid spacing, scaled rather than chosen: the observed
    # deficit runs from 3.6e-06 to 4.6e-05 against a spacing near 5e-03.
    tol = float(grid[1] - grid[0])

    return {
        "mu": mu,
        "mass": float(model.mass),
        "charge_squared": float(model.magnetic_r2),
        "r_minus": r_minus,
        "r_plus": r_plus,
        "gap": gap,
        "excision_midpoint": excision,
        "child_start": child_start,
        "child_end": child_end,
        "sampling_spacing": float(grid[1] - grid[0]),
        "shift_turn_radius": turn,
        "shift_turn_above_excision": bool(turn > excision),
        "outgoing_zero": zero,
        "zero_at_r_plus": (None if zero is None else bool(abs(zero - r_plus) < 1e-4)),
        "ordering_tolerance": tol,
        "order_holds": bool(
            r_minus < excision
            and excision < r_plus
            and zero is not None
            and r_plus - tol <= zero <= child_start + tol
            and child_start < child_end),
        "zero_deficit_below_r_plus": (None if zero is None else zero - r_plus),
        "zero_deficit_in_cells": (None if zero is None
                                  else (zero - r_plus) / float(grid[1] - grid[0])),
        "order": ["r_minus", "excision", "r_plus", "outgoing_zero",
                  "child_start", "child_end"],
        "separations": {
            "excision_to_horizon": r_plus - excision,
            "horizon_to_child_start": child_start - (zero if zero else r_plus),
            "child_width": child_end - child_start,
            "parent_width": 32.0 - child_end,
        },
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mu", nargs="+", type=float,
                        default=[1.0005, 1.001, 1.002, 1.005, 1.010, 1.020,
                                 1.030, 1.040, 1.060, 1.080, 1.100])
    parser.add_argument("--out", default="nsc-rn-geometry-map-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    entries = [_entry(mu) for mu in sorted(args.mu)]
    orders_hold = all(e["order_holds"] for e in entries)
    zeros_at_r_plus = all(e["zero_at_r_plus"] for e in entries)
    turns_below = sum(1 for e in entries if not e["shift_turn_above_excision"])
    turn_ratios = [e["shift_turn_radius"] / e["excision_midpoint"] for e in entries]

    record = {
        "schema": "NSC-RN-GEOMETRY-MAP-v1",
        "api_version": "nsc-rn-geometry-map-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_is": (
            "A map of one space and its nesting, as the one free parameter moves. Every "
            "quantity is classical: a charged spherical solution in areal coordinates, its "
            "two horizons, the two cuts the chart defines, the shift, and the causal "
            "membrane located from the metric."
        ),
        "reading_the_map": (
            "The ordering column is the point. If the child cut ever crossed the horizon, or "
            "the horizon crossed the inner edge, the nesting would have a transition in it. "
            "The shift turn moving relative to the excision is the other thing worth "
            "watching, because the default excision sits above it for small mu and below it "
            "for large ones."
        ),
        "scope_note": (
            "A map of one family, not of nested spaces in general. It does not measure the "
            "actual universe, and one chart with one child interval is not a hierarchy. "
            "What the map can say is whether the ordering of the features is stable as the "
            "parameter moves."
        ),
        "summary": {
            "entries": len(entries),
            "mu_min": entries[0]["mu"],
            "mu_max": entries[-1]["mu"],
            "order_holds_everywhere": orders_hold,
            "outgoing_zero_at_r_plus_everywhere": zeros_at_r_plus,
            "entries_with_shift_turn_below_excision": turns_below,
            "shift_turn_over_excision_range": [min(turn_ratios), max(turn_ratios)],
            "the_default_edge_never_samples_the_shift_maximum": bool(
                turns_below == len(entries)),
        },
        "entries": entries,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print()
    print("%7s %9s %9s %9s %9s %9s %9s %6s"
          % ("mu", "r_minus", "excision", "r_plus", "zero", "child0", "child1", "order"))
    for entry in entries:
        print("%7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %6s"
              % (entry["mu"], entry["r_minus"], entry["excision_midpoint"],
                 entry["r_plus"], entry["outgoing_zero"] or float("nan"),
                 entry["child_start"], entry["child_end"],
                 "yes" if entry["order_holds"] else "NO"))
    print()
    print(f"  ordering holds everywhere : {orders_hold}")
    print(f"  zero sits at r_plus       : {zeros_at_r_plus}")
    print(f"  shift turn below excision : {turns_below} of {len(entries)} entries")
    deficits = [e["zero_deficit_in_cells"] for e in entries
                if e["zero_deficit_in_cells"] is not None]
    print(f"  zero deficit below r_plus : {min(deficits):.4f} to {max(deficits):.4f} cells"
          f"  (tolerance is one cell)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())