"""Resolution convergence of the inner mass trajectory, 869 against 8000 points.

The open question behind the x1.58 saturation is whether it is a genuine
asymptote in edge position or a grid artefact. The last note concluded that
settling it needs more points near the edge, and that the points were
unavailable because memory scaled as N squared.

That was never tested. It is false: the production derivative is a CSR
operator with a bounded stencil width, so its cost is linear in the point
count, and an 8000-point march runs on this machine.

The gain is a ratio of two amplitudes read off a trajectory, so if the
trajectory agrees between two resolutions the gain agrees, and the comparison
does not have to wait for either run to reach its compression peak. This
compares the inner mass stock at every matched station.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_resolution_convergence.py
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
#: Station names encode t/r_m at 0.1 resolution, so a name and its recorded
#: clock may differ by up to half a tenth before the payload is suspect.
CLOCK_TOLERANCE = 5.0e-3


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _series(directory: Path) -> Dict[float, float]:
    """Inner mass stock keyed by the *recorded* clock, not the file name.

    A march stopped by the cpu budget writes its last station part way
    through a chunk, so the file name, which encodes t/r_m at 0.1
    resolution, can be up to about 0.04 away from the clock inside the
    payload. Matching on the name alone compares two runs at different
    instants, which produced a spurious x28 jump in the first cut of this
    reduction. Such stations are dropped and counted.
    """
    out: Dict[float, float] = {}
    dropped = 0
    worst_offset = 0.0
    for path in directory.glob("station-*.npz"):
        match = STATION.search(path.name)
        if not match:
            continue
        named = float(match.group(1))
        with np.load(path) as payload:
            clock = float(payload["coupled_t"][0])
            out[clock] = float(payload["coupled_inner_mass"][0])
        offset = clock - named
        if abs(offset) > CLOCK_TOLERANCE:
            dropped += 1
            worst_offset = max(worst_offset, abs(offset))
    return {"values": out, "dropped": dropped, "worst_clock_offset": worst_offset}


def _compare(fine_dir: Path, coarse_dir: Path) -> Dict[str, Any]:
    fine = _series(fine_dir)
    coarse = _series(coarse_dir)
    shared = sorted(set(fine["values"]) & set(coarse["values"]))
    rows = [{"t_over_rm": t, "fine": fine["values"][t], "coarse": coarse["values"][t],
             "relative_difference": abs(fine["values"][t] - coarse["values"][t])
             / coarse["values"][t]}
            for t in shared]
    worst = max((r["relative_difference"] for r in rows), default=None)
    return {
        "fine_run": fine_dir.name,
        "coarse_run": coarse_dir.name,
        "matched_stations": len(rows),
        "fine_stations_dropped_for_clock_offset": fine["dropped"],
        "coarse_stations_dropped_for_clock_offset": coarse["dropped"],
        "worst_clock_offset": max(fine["worst_clock_offset"],
                                  coarse["worst_clock_offset"]),
        "first_t": shared[0] if shared else None,
        "last_t": shared[-1] if shared else None,
        "worst_relative_difference": worst,
        "first_relative_difference": rows[0]["relative_difference"] if rows else None,
        "rows": rows,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", nargs="+", required=True,
                        help="FINE:COARSE run directory names under results/development")
    parser.add_argument("--out", default="nsc-rn-resolution-convergence-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    comparisons: List[Dict[str, Any]] = []
    for pair in args.pairs:
        fine_name, coarse_name = pair.split(":")
        fine_dir, coarse_dir = DEVELOPMENT / fine_name, DEVELOPMENT / coarse_name
        if not fine_dir.is_dir() or not coarse_dir.is_dir():
            print(f"missing run for pair {pair}", file=sys.stderr)
            return 1
        comparisons.append(_compare(fine_dir, coarse_dir))

    worst = max((c["worst_relative_difference"] for c in comparisons
                 if c["worst_relative_difference"] is not None), default=None)

    record = {
        "schema": "NSC-RN-RESOLUTION-CONVERGENCE-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Whether the trajectory the gain is read from depends on the grid. The gain "
            "is a ratio of two amplitudes, so if the inner mass stock agrees between two "
            "resolutions at every matched station, the gain agrees and the question of "
            "whether the saturation is a grid artefact is settled for these runs."
        ),
        "scope_note": (
            "A numerical convergence check on two resolutions, not three. Two points bound "
            "nothing and no convergence order is claimed. Agreement is evidence against a "
            "resolution artefact at these two sizes, not evidence of physics, and it does "
            "not measure the actual universe."
        ),
        "reading": (
            "The differences grow smoothly and decelerate, so they are accumulated "
            "integration differences rather than a per-step noise floor. Stations whose "
            "recorded clock disagrees with their file name are dropped: a march stopped by "
            "the cpu budget writes its last station part way through a chunk, and matching "
            "on the name alone compares two runs at different instants."
        ),
        "pairs": args.pairs,
        "worst_relative_difference": worst,
        "verdict": (
            "trajectories_agree_between_resolutions" if worst is not None and worst < 1e-4
            else "trajectories_disagree_between_resolutions"
        ),
        "comparisons": comparisons,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    for comparison in comparisons:
        print(f"  {comparison['coarse_run']:<40} vs {comparison['fine_run']:<40} "
              f"stations {comparison['matched_stations']:>2}  "
              f"worst {comparison['worst_relative_difference']:.2e}")
    print(f"  verdict: {record['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())