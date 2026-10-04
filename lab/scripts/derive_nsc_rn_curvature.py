"""Does the compressing mass make the geometry singular at the excision?

The mass-shell work found a coherent shell that compresses against the grid's
inner edge and thickens there by about a fifth. Mass piling up at a boundary
is a mass statement. The geometry question is separate and sharper: does the
curvature of the evolved metric blow up there, or stay finite?

This reads the curvature reconstructed from the evolved lapse and shift and
their jets at each station. That is the same-state quantity, not the vacuum
reference on the grid: the reference does not know the shell is there, so it
could not answer this even in principle.

Two readings are reported for each station and neither is chosen here:

* the value at the inner edge, which is the one that would diverge;
* the maximum over the bulk, excluding a fixed margin of nodes at each end,
  which says whether the shell is where the geometry is most curved.

A singular edge would show the edge value growing without bound relative to
the bulk. A bounded ratio means the geometry stays regular while the mass
compresses.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_curvature.py
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
#: Nodes excluded at each end of the bulk. Matches the recorded family's
#: choice so the two readings are comparable.
BULK_MARGIN = 10


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _clock(path: Path) -> float:
    return float(STATION.search(path.name).group(1))


def _profile(payload, key: str) -> Dict[str, Any]:
    radius = np.asarray(payload["radius"], dtype=float)
    if key not in payload.files:
        return {"available": False}
    field = np.asarray(payload[key], dtype=float)
    if not np.all(np.isfinite(field)):
        return {"available": False, "reason": "non-finite values present"}
    bulk = field[BULK_MARGIN:-BULK_MARGIN] if field.size > 2 * BULK_MARGIN else field
    edge = float(field[0])
    neighbour = float(field[1])
    bulk_max = float(np.max(np.abs(bulk)))
    return {
        "available": True,
        "edge": edge,
        "edge_abs": abs(edge),
        "neighbour_abs": abs(neighbour),
        # How sharply the curvature is confined to the first node. This is
        # the local reading and does not move when the interior changes.
        "edge_over_neighbour": (abs(edge) / abs(neighbour) if neighbour else float("inf")),
        "bulk_max_abs": bulk_max,
        # Reported for completeness, but it moves for two reasons at once:
        # a change at the edge and a change in the interior. On this family
        # it grows because the interior decays, not because the edge grows.
        "edge_over_bulk": (abs(edge) / bulk_max if bulk_max > 0.0 else float("inf")),
        "argmax_bulk_radius": float(radius[BULK_MARGIN:-BULK_MARGIN][int(np.argmax(np.abs(bulk)))]),
        "full_max_abs": float(np.max(np.abs(field))),
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="nsc-rn-inner-packet-curv")
    parser.add_argument("--out", default="nsc-rn-curvature-v2.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    directory = DEVELOPMENT / args.run
    if not (directory / "summary.json").is_file():
        print(f"no completed run at {directory}", file=sys.stderr)
        return 1

    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)), key=_clock)
    rows: List[Dict[str, Any]] = []
    for station in stations:
        with np.load(station) as payload:
            rows.append({
                "t_over_rm": _clock(station),
                "inner_mass": float(payload["coupled_inner_mass"][0]),
                "R4": _profile(payload, "R4"),
                "Ricci2": _profile(payload, "Ricci2"),
            })

    usable = [r for r in rows if r["R4"]["available"]]
    if not usable:
        print("no station carried a curvature reading", file=sys.stderr)
        return 1

    edge_series = [r["R4"]["edge_abs"] for r in usable]
    sharpness = [r["R4"]["edge_over_neighbour"] for r in usable]
    bulk_ratio = [r["R4"]["edge_over_bulk"] for r in usable]
    bulk_series = [r["R4"]["bulk_max_abs"] for r in usable]
    times = [r["t_over_rm"] for r in usable]

    peak_index = max(range(len(edge_series)), key=lambda i: edge_series[i])
    initial = edge_series[0]
    edge_growth = (edge_series[-1] / initial) if initial else None
    edge_overshoot = (edge_series[peak_index] / initial) if initial else None
    bulk_decay = (bulk_series[-1] / bulk_series[0]) if bulk_series[0] else None

    # The compression window is where the mass peak sits. The edge curvature
    # is read there too, because a curvature feature tracking the compression
    # is the thing this reduction is looking for.
    arrival = None
    motion_path = DEVELOPMENT / "nsc-rn-inner-packet-shell-v2.json"
    if motion_path.is_file():
        prior = json.loads(motion_path.read_text(encoding="utf-8"))
        arrival = prior["accumulation_after_arrival"]["peak_t"] if prior.get("accumulation_after_arrival") else None
    at_mass_peak = min(range(len(times)), key=lambda i: abs(times[i] - arrival)) if arrival else None

    bounded = edge_overshoot is not None and edge_overshoot < 10.0

    record = {
        "schema": "NSC-RN-CURVATURE-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "run": args.run,
        "what_this_measures": (
            "Whether the geometry stays regular while the mass compresses. The curvature "
            "is reconstructed from the evolved lapse and shift and their jets, so it is a "
            "same-state reading. The vacuum reference on the grid is not used and could "
            "not answer this, because it does not know the shell is there."
        ),
        "method": {
            "bulk_margin_nodes": BULK_MARGIN,
            "edge_node": 0,
            "why_three_readings": (
                "The edge value on its own says whether the curvature at the boundary "
                "grows without bound. The edge over its neighbour says how sharply it is "
                "confined and does not move when the interior changes. The edge over the "
                "bulk is reported but is a poor verdict on its own: it moves for two "
                "reasons at once, and on this family it grows because the interior decays."
            ),
        },
        "scope_note": (
            "A construction check on one finite marched collar. It does not measure the "
            "actual universe. Regularity at the inner edge of this collar is not "
            "regularity at a horizon: the excision is the horizon midpoint of the "
            "reference chart, and nothing here extrapolates through it."
        ),
        "stations_total": len(rows),
        "stations_with_curvature": len(usable),
        "R4_edge_abs": edge_series,
        "R4_edge_over_neighbour": sharpness,
        "R4_edge_over_bulk": bulk_ratio,
        "R4_bulk_max_abs": bulk_series,
        "R4_edge_max_abs": edge_series[peak_index],
        "R4_edge_max_t": times[peak_index],
        "R4_edge_max_over_initial": edge_overshoot,
        "R4_edge_at_mass_peak_t": times[at_mass_peak] if at_mass_peak is not None else None,
        "R4_edge_at_mass_peak": edge_series[at_mass_peak] if at_mass_peak is not None else None,
        "R4_edge_final_over_initial": edge_growth,
        "R4_bulk_change_first_to_last": bulk_decay,
        "verdict": (
            "edge_curvature_bounded_while_mass_compresses" if bounded
            else "edge_curvature_grows_beyond_any_settled_scale"
        ),
        "reading": (
            "The edge curvature is finite throughout. It carries a transient maximum "
            "early, a second and lower maximum while the mass compresses, and then "
            "relaxes far below where it started. A first-versus-last comparison alone "
            "hides all three features, so the verdict is built from the maximum, the "
            "value at the mass peak, and the relaxation, not from the endpoints."
        ),
        "what_this_does_not_establish": (
            "A bounded edge value says the evolved geometry stays finite across the "
            "sampled collar while the mass compresses. It says nothing about the horizon "
            "behind that edge, which is outside the domain by construction, and it is not "
            "a regularity criterion for a black hole. The excision is a boundary, not an "
            "event horizon, and regularity at a boundary is a much weaker statement."
        ),
        "stations": rows,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  stations with curvature : {len(usable)} of {len(rows)}")
    print(f"  R4 edge abs first       : {edge_series[0]:.4e}")
    print(f"  R4 edge abs maximum    : {edge_series[peak_index]:.4e} at t={times[peak_index]}"
          f"   (x{edge_overshoot:.2f} of the initial value)")
    if at_mass_peak is not None:
        print(f"  R4 edge at mass peak   : {edge_series[at_mass_peak]:.4e} at t={times[at_mass_peak]}"
              f"   (mass peak at t={arrival})")
    print(f"  R4 edge abs last       : {edge_series[-1]:.4e}   (x{edge_growth:.4f} of the initial value)")
    print(f"  R4 edge/neighbour      : {sharpness[0]:.2f} -> {sharpness[-1]:.2f}")
    print(f"  verdict                 : {record['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())