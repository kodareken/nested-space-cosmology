"""Locate the trapping horizon in the evolved metric and compare it to the chart.

A horizon is not a boundary this programme chose. The grid's inner edge is
where the march stopped computing; a horizon is an outermost simple zero of
the outgoing speed in the evolved metric, and it has to be found rather than
placed. This locates it in recorded states and sets it against the three
radii the chart already defines: the reference outer horizon, the start of the
labelled child interval, and the end of it.

The comparison matters for the nested picture. If the chart's child interval
were a trapped region, its boundaries and the trapping horizon would
coincide. If they do not, the labelled interval and the causal region are
different objects and saying one is the other would be wrong.

Reads only recorded station payloads. No time step is taken and nothing is
marched. The lapse and shift are already in the recorded states.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_trapping_horizon.py
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


def chart_radii(mu: float) -> Dict[str, float]:
    geometry = driver.case_geometry(mu)
    model = pg.action_model(geometry["mass"], geometry["r_m"])
    r_minus = float(model.horizons["r_minus"])
    r_plus = float(model.horizons["r_plus"])
    gap = r_plus - r_minus
    return {
        "mu": mu,
        "r_minus": r_minus,
        "r_plus": r_plus,
        "excision_midpoint": 0.5 * (r_minus + r_plus),
        "child_start": r_plus + gap / 16.0,
        "child_end": r_plus + gap / 2.0,
    }


def _reduce(directory: Path, mu: float) -> Dict[str, Any]:
    radii = chart_radii(mu)
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
        zeros = int(np.count_nonzero(np.diff(np.sign(outgoing)) != 0))
        trapping = located["trapping_horizon"]
        rows.append({
            "t_over_rm": clock,
            "r_inner": float(radius[0]),
            "trapping_horizon": trapping,
            "located": bool(located["located"]),
            "zeros_of_outgoing_speed": zeros,
            "outgoing_speed_at_inner_edge": float(outgoing[0]),
            "inner_edge_inside_trapping_region": bool(outgoing[0] < 0.0),
            "trapping_minus_r_plus": (trapping - radii["r_plus"]
                                      if trapping is not None else None),
            "trapping_minus_child_start": (trapping - radii["child_start"]
                                           if trapping is not None else None),
            "trapping_minus_child_end": (trapping - radii["child_end"]
                                         if trapping is not None else None),
        })

    found = [r["trapping_horizon"] for r in rows if r["trapping_horizon"] is not None]
    return {
        "run": directory.name,
        "chart": radii,
        "stations": len(rows),
        "stations_with_a_trapping_horizon": len(found),
        "located_everywhere": len(found) == len(rows),
        "single_zero_everywhere": all(r["zeros_of_outgoing_speed"] == 1 for r in rows),
        "inner_edge_trapped_everywhere": all(
            r["inner_edge_inside_trapping_region"] for r in rows),
        "trapping_first": found[0] if found else None,
        "trapping_last": found[-1] if found else None,
        "trapping_drift": (found[-1] - found[0]) if found else None,
        "offset_from_r_plus_first": rows[0]["trapping_minus_r_plus"] if rows else None,
        "offset_from_r_plus_last": rows[-1]["trapping_minus_r_plus"] if rows else None,
        "offset_from_child_start_last": (rows[-1]["trapping_minus_child_start"]
                                         if rows else None),
        "offset_from_child_end_last": (rows[-1]["trapping_minus_child_end"]
                                       if rows else None),
        "below_child_start_everywhere": all(
            (r["trapping_horizon"] is not None
             and r["trapping_minus_child_start"] < 0.0) for r in rows),
        "rows": rows,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", required=True,
                        help="RUN:MU pairs, each naming a run directory under "
                             "results/development and the mass ratio it was built at")
    parser.add_argument("--out", default="nsc-rn-trapping-horizon-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    reduced = []
    for pair in args.runs:
        name, _, mu = pair.partition(":")
        directory = DEVELOPMENT / name
        if not directory.is_dir():
            print(f"no run at {directory}", file=sys.stderr)
            return 1
        reduced.append(_reduce(directory, float(mu or 1.04)))

    record = {
        "schema": "NSC-RN-TRAPPING-HORIZON-v1",
        "api_version": "nsc-rn-trapping-horizon-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Where the trapping horizon is in the evolved metric, as an outermost simple "
            "zero of the outgoing speed, against the reference outer horizon and the two "
            "radii that define the chart's labelled child interval."
        ),
        "why_it_matters": (
            "The grid's inner edge is where this programme stopped computing and is not a "
            "horizon. If the inner edge sits inside the trapping region, then the matter "
            "measured compressing against it is matter that cannot escape outward, which "
            "is a different and much stronger statement than piling up against a boundary "
            "that was chosen."
        ),
        "scope_note": (
            "A construction check on recorded finite states. It does not measure the "
            "actual universe. A simple zero of the outgoing characteristic on a sampled "
            "slice is the standard numerical proxy for a trapping surface; it is not a "
            "proof of a global event horizon, and no claim of one is made here."
        ),
        "criterion": {
            "definition": "outermost simple zero of the outgoing speed in the evolved metric",
            "source": "nsc_rn_observables.trapping_horizon",
            "event_horizon_copied_from_trapping": False,
            "inner_edge_is_not_a_horizon": (
                "The inner edge is a grid boundary chosen by the run, carried by no causal "
                "argument. It is reported separately and never called a horizon."
            ),
        },
        "all_located": all(r["located_everywhere"] for r in reduced),
        "all_single_zero": all(r["single_zero_everywhere"] for r in reduced),
        "all_inner_edges_trapped": all(r["inner_edge_trapped_everywhere"] for r in reduced),
        "all_below_child_start": all(r["below_child_start_everywhere"] for r in reduced),
        "child_interval_is_the_trapped_region": any(
            r["offset_from_child_start_last"] is not None
            and r["offset_from_child_start_last"] > 0.0 for r in reduced),
        "runs": reduced,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    for run in reduced:
        chart = run["chart"]
        print(f"  {run['run']}")
        print(f"     r_plus={chart['r_plus']:.4f}  child_start={chart['child_start']:.4f}  "
              f"child_end={chart['child_end']:.4f}  inner_edge={run['rows'][0]['r_inner']:.4f}")
        print(f"     trapping: {run['trapping_first']:.4f} -> {run['trapping_last']:.4f}"
              f"   (drift {run['trapping_drift']:+.5f},"
              f" offset from r_plus {run['offset_from_r_plus_last']:+.5f})")
        print(f"     below child_start everywhere: {run['below_child_start_everywhere']}"
              f"   single zero: {run['single_zero_everywhere']}"
              f"   edge trapped: {run['inner_edge_trapped_everywhere']}")
    print(f"  child interval is the trapped region: "
          f"{record['child_interval_is_the_trapped_region']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())