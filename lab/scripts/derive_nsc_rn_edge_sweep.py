"""Sweep the grid's inner edge inboard and ask whether the trend continues.

Moving the excision from the horizon midpoint to 0.85 changed two numbers in
the same direction: the compression grew from x1.20 to x1.39, and the
relaxation weakened from x0.0145 to x0.69. Both had been truncation
artifacts of a domain that ended in the middle of the accumulation.

This asks whether that continues as the edge approaches the inner horizon at
r_minus = 0.754. If the compression keeps growing and the relaxation keeps
weakening, the effect is a function of how much of the accumulation the
domain contains, and the shallow-edge numbers were lower bounds throughout.

If it turns over instead, the deepening is running into something the chart
stops supporting well before r_minus, which is equally worth knowing.

Runs are named explicitly. A glob here would be dangerous twice over: it
would mix edge positions, and those runs differ only in an argument to one
function call.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_edge_sweep.py
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
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_rn_mass_shell as shell  # noqa: E402

DEVELOPMENT = LAB / "results/development"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


r_outer = 32.0


def _clock(path: Path) -> float:
    return float(STATION.search(path.name).group(1))


def _reduce(directory: Path) -> Dict[str, Any]:
    summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)), key=_clock)
    rows: List[Dict[str, Any]] = []
    r4_edge: List[float] = []
    for station in stations:
        with np.load(station) as payload:
            rows.append(shell.shell_state(payload, base_mass=float(summary["mu"])))
            rows[-1]["t_over_rm"] = _clock(station)
            if "R4" in payload.files:
                field = np.asarray(payload["R4"], dtype=float)
                r4_edge.append(abs(float(field[0])) if np.all(np.isfinite(field)) else float("nan"))
    motion = shell.shell_motion(rows)
    accumulation = shell.accumulation_after_arrival(rows, motion)
    first = r4_edge[0] if r4_edge else None
    last = r4_edge[-1] if r4_edge else None
    requested = summary.get("excision_requested")
    # Runs made before the override existed carry no request, and the edge
    # they used was the default midpoint. That is read, not assumed.
    is_default = bool(summary.get("excision_default", requested is None))
    r_inner = summary["r_inner"]
    r_minus = summary["r_minus"]
    points = summary["points"]
    spacing = (r_outer - r_inner) / (points - 1)
    gap = r_inner - r_minus
    gap_cells = gap / spacing if spacing > 0.0 else float("inf")
    return {
        "run": directory.name,
        "r_inner": r_inner,
        "r_minus": r_minus,
        "r_outer": r_outer,
        "points": points,
        "grid_spacing": spacing,
        "gap_to_inner_horizon": gap,
        "gap_in_cells": gap_cells,
        # A boundary less than one cell from the horizon is not
        # distinguishable from the horizon on this grid. Reporting depth in
        # r_m alone would hide that, because the gaps look meaningfully
        # different there and are not.
        "edge_resolvable_from_horizon": bool(gap_cells >= 1.0),
        "is_default_edge": is_default,
        "excision_requested": requested,
        "mu": summary["mu"],
        "eps": summary["eps"],
        "points": summary["points"],
        "target_over_rm": summary["target_over_rm"],
        "stop_reason": summary["stop_reason"],
        "stations": len(rows),
        "motion": motion,
        "accumulation_after_arrival": accumulation,
        "fraction_below_inner_edge": (
            rows[-1]["gain_below_inner_edge"] / rows[-1]["total_gain"]
            if rows[-1]["total_gain"] else None),
        "R4_edge_first": first,
        "R4_edge_last": last,
        "R4_edge_final_over_initial": (last / first) if first else None,
        "R4_edge_max": max(r4_edge) if r4_edge else None,
        "R4_edge_max_over_initial": (max(r4_edge) / first) if first else None,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", required=True,
                        help="run directory names under results/development")
    parser.add_argument("--out", default="nsc-rn-edge-sweep-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    reduced = []
    for name in args.runs:
        directory = DEVELOPMENT / name
        if not (directory / "summary.json").is_file():
            print(f"no completed run at {directory}", file=sys.stderr)
            return 1
        reduced.append(_reduce(directory))
    # Deepest edge last, so the series reads as the edge moves inboard.
    reduced.sort(key=lambda r: -r["r_inner"])

    for key in ("mu", "eps", "points", "target_over_rm"):
        values = [r[key] for r in reduced]
        if len(set(values)) != 1:
            print(f"runs differ in {key}: {values}. This is not an edge sweep.",
                  file=sys.stderr)
            return 2
    edges = [r["r_inner"] for r in reduced]
    if len(set(edges)) != len(edges):
        print(f"runs do not span distinct edge positions: {edges}", file=sys.stderr)
        return 2

    gains = [r["accumulation_after_arrival"]["amplitude_gain"]
             for r in reduced if r["accumulation_after_arrival"]]
    retained = [r["fraction_below_inner_edge"] for r in reduced
                if r["fraction_below_inner_edge"] is not None]
    relax = [r["R4_edge_final_over_initial"] for r in reduced
             if r["R4_edge_final_over_initial"] is not None]

    monotonic_gain = all(b >= a for a, b in zip(gains, gains[1:]))
    monotonic_relax = all(b >= a for a, b in zip(relax, relax[1:]))
    gaps = [r["gap_in_cells"] for r in reduced]
    resolvable = [r["edge_resolvable_from_horizon"] for r in reduced]
    resolvable_gains = [g for g, ok in zip(gains, resolvable) if ok]

    monotonic_where_resolvable = all(
        b >= a for a, b in zip(resolvable_gains, resolvable_gains[1:]))

    record = {
        "schema": "NSC-RN-EDGE-SWEEP-v1",
        "api_version": shell.API_VERSION,
        "creation_only": True,
        "producing_commit": _git_head(),
        "run_selection": list(args.runs),
        "what_this_measures": (
            "Whether the effect of moving the grid's inner edge inboard continues as the "
            "edge approaches the inner horizon, or turns over because the chart stops "
            "supporting the region well before r_minus."
        ),
        "scope_note": (
            "A construction check on finite marched collars. It does not measure the "
            "actual universe. Every edge in this sweep is still a boundary and none is "
            "the inner horizon, so none of these readings is a statement about regularity "
            "at a horizon."
        ),
        "held_fixed": {key: sorted({r[key] for r in reduced})
                       for key in ("mu", "eps", "points", "target_over_rm")},
        "sweep": {
            "r_inner": edges,
            "gap_to_inner_horizon": [r["gap_to_inner_horizon"] for r in reduced],
            "gap_in_cells": gaps,
            "edge_resolvable_from_horizon": resolvable,
            "deepest_resolvable_edge": max(
                (r["r_inner"] for r in reduced if r["edge_resolvable_from_horizon"]),
                default=None),
            "amplitude_gain": gains,
            "gain_increases_with_depth": monotonic_gain,
            "gain_increases_with_depth_where_resolvable": monotonic_where_resolvable,
            "fraction_below_inner_edge": retained,
            "R4_edge_final_over_initial": relax,
            "relaxation_weakens_with_depth": monotonic_relax,
            "verdict": (
                "effect_continues_to_the_deepest_edge" if monotonic_gain and monotonic_relax
                else "effect_does_not_continue_monotonically"
            ),
            "resolution_limit": (
                "The sweep is resolution-limited before it is depth-limited. Edges below "
                "about one grid cell from the inner horizon are not distinguishable from "
                "the horizon on this grid, so the deepest runs measure the same place. "
                "Depth in r_m alone overstates how far the sweep has actually gone, and "
                "the trend cannot be pushed further without more points near the edge."
            ),
        },
        "runs": reduced,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  r_inner      : {edges}")
    print(f"  gap to r_-   : {['%.4f' % g for g in record['sweep']['gap_to_inner_horizon']]}")
    print(f"  gap in CELLS : {['%.3f' % g for g in gaps]}")
    print(f"  resolvable   : {resolvable}   deepest resolvable edge: "
          f"{record['sweep']['deepest_resolvable_edge']}")
    print(f"  gain         : {['%.4f' % g for g in gains]}  monotone={monotonic_gain}"
          f"  monotone_where_resolvable={monotonic_where_resolvable}")
    print(f"  gain below   : {['%.4f' % g for g in retained]}")
    print(f"  relax (final/initial): {['%.4f' % g for g in relax]}  monotone={monotonic_relax}")
    print(f"  verdict      : {record['sweep']['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())