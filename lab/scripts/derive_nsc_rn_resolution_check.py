"""Resolution control on the inner-packet compression.

The density peak that compresses sits at the grid's left edge, and the left
edge is exactly where the one-sided SBP closure lives. A growth of about a
fifth measured while pinned there could therefore be partly the closure
accumulating rather than the geometry. This checks it.

The same release radius, the same packet, the same mass ratio and source
strength, marched at two grid resolutions. If the compression gain agrees
between them to a few percent, it is not a resolution artifact. If it does
not, the single-resolution result is a numerical artifact and says so.

This is the first same-parameter, different-resolution replay in this line.
Earlier notes repeatedly recorded that no such replay existed, which left
every near-throat quantity without a convergence control.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_resolution_check.py
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


def _clock(path: Path) -> float:
    return float(STATION.search(path.name).group(1))


def _reduce(directory: Path) -> Dict[str, Any]:
    summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)), key=_clock)
    rows: List[Dict[str, Any]] = []
    for station in stations:
        with np.load(station) as payload:
            state = shell.shell_state(payload, base_mass=float(summary["mu"]))
            state["t_over_rm"] = _clock(station)
        rows.append(state)
    motion = shell.shell_motion(rows)
    accumulation = shell.accumulation_after_arrival(rows, motion)
    return {
        "run": directory.name,
        "points": summary["points"],
        "release_radius_over_rm": summary["packet_center_over_rm"],
        "mu": summary["mu"],
        "eps": summary["eps"],
        "target_over_rm": summary["target_over_rm"],
        "stop_reason": summary["stop_reason"],
        "cpu_seconds": summary["cpu_seconds"],
        "stations": len(rows),
        "motion": motion,
        "accumulation_after_arrival": accumulation,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse", default="nsc-rn-inner-packet-v1")
    parser.add_argument("--fine", default="nsc-rn-inner-packet-c3-fine")
    parser.add_argument("--out", default="nsc-rn-resolution-check-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    coarse = _reduce(DEVELOPMENT / args.coarse)
    fine = _reduce(DEVELOPMENT / args.fine)
    runs = [coarse, fine]

    held = {}
    for key in ("release_radius_over_rm", "mu", "eps", "target_over_rm"):
        values = [run[key] for run in runs]
        held[key] = values
        if len(set(values)) != 1:
            print(f"resolution control is invalid: {key} differs across runs: {values}",
                  file=sys.stderr)
            return 2

    gains = [run["accumulation_after_arrival"]["amplitude_gain"]
             for run in runs if run["accumulation_after_arrival"]]
    arrivals = [run["motion"]["inner_edge_arrival_t"]
                for run in runs if run["motion"].get("inner_edge_arrival_t") is not None]
    rises = [run["accumulation_after_arrival"]["rise_time"]
             for run in runs if run["accumulation_after_arrival"]]

    relative = None
    if len(gains) == 2 and gains[1]:
        relative = abs(gains[0] - gains[1]) / gains[1]

    record = {
        "schema": "NSC-RN-RESOLUTION-CHECK-v1",
        "api_version": shell.API_VERSION,
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Whether the compression gain measured at one resolution survives the same "
            "construction at another. The growing peak sits at the grid's left edge, where "
            "the one-sided SBP closure lives, so it could be closure accumulating rather "
            "than geometry."
        ),
        "scope_note": (
            "A construction check on two finite marched collars. It does not measure the "
            "actual universe. Agreement at two resolutions is evidence against a "
            "resolution artifact, not evidence of physics. Two points is not a "
            "convergence order and none is claimed."
        ),
        "held_fixed": held,
        "points": [run["points"] for run in runs],
        "cpu_seconds": [run["cpu_seconds"] for run in runs],
        "amplitude_gain": gains,
        "amplitude_gain_relative_difference": relative,
        "edge_arrival_t": arrivals,
        "rise_time": rises,
        "verdict": (
            "gain_holds_under_resolution_change" if relative is not None and relative < 0.05
            else "gain_does_not_hold_under_resolution_change"
        ),
        "what_this_does_not_establish": (
            "Two resolutions bracket nothing. If the gain moves by more than a few percent "
            "between them the single-resolution result is not usable; if it does not, the "
            "result is merely not obviously an artifact. A Richardson extrapolation needs "
            "a third resolution and is not attempted here."
        ),
        "runs": runs,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  points          : {record['points']}")
    print(f"  amplitude gain  : {['%.4f' % g for g in gains]}")
    if relative is not None:
        print(f"  relative diff   : {relative:.4f}")
    print(f"  arrival t       : {arrivals}")
    print(f"  rise time       : {['%.3f' % r for r in rises]}")
    print(f"  verdict         : {record['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())