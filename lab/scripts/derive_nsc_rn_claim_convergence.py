"""Resolution check on the claims that were only ever measured coarsely.

The compression gain was reported as x1.58 from a coarse grid and as x2.05
at 8000 points, because a gain read off dM/dr is not converged when the field
it differentiates is. That correction was made for the gain alone. Three
other headline numbers had not been rechecked under the same refinement,
and one of them, the retained fraction, had been sitting on the mass stock
without a resolution check of its own.

This compares the two resolutions at matched times, which is the only valid
comparison: the retained fraction grows with time, so comparing a run that
stopped at t = 2.0 against one that reached t = 6.0 would measure the clock
rather than the resolution.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_claim_convergence.py
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
#: Times at which both resolutions are recorded.
TIMES = (1.4, 1.8, 2.0)


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _at(directory: Path, wanted: float, tolerance: float, base_mass: float):
    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)),
                      key=lambda p: float(STATION.search(p.name).group(1)))
    for station in stations:
        with np.load(station) as payload:
            clock = float(payload["coupled_t"][0])
            if abs(clock - wanted) <= tolerance:
                mass = np.asarray(payload["coupled_mass"], dtype=float)
                total = float(mass[-1] - base_mass)
                if total == 0.0:
                    return None
                return {"t": clock,
                        "retained": float(mass[0] - base_mass) / total,
                        "stations": len(stations)}
    return None


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", nargs="+", required=True,
                        help="COARSE:FINE run names under results/development, for "
                             "the same inner edge")
    parser.add_argument("--mu", type=float, default=1.04)
    parser.add_argument("--out", default="nsc-rn-claim-convergence-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    comparisons: List[Dict[str, Any]] = []
    for pair in args.pairs:
        coarse_name, _, fine_name = pair.partition(":")
        rows = []
        for wanted in TIMES:
            coarse = _at(DEVELOPMENT / coarse_name, wanted, 0.06, args.mu)
            fine = _at(DEVELOPMENT / fine_name, wanted, 0.06, args.mu)
            if coarse is None or fine is None:
                continue
            rows.append({
                "requested_t": wanted,
                "coarse_t": coarse["t"],
                "fine_t": fine["t"],
                "coarse_retained": coarse["retained"],
                "fine_retained": fine["retained"],
                "absolute_difference": abs(coarse["retained"] - fine["retained"]),
                "relative_difference": (abs(coarse["retained"] - fine["retained"])
                                        / coarse["retained"]),
                "matched": abs(coarse["t"] - fine["t"]) < 0.02,
            })
        if rows:
            comparisons.append({
                "coarse": coarse_name, "fine": fine_name, "rows": rows,
                "worst_relative_difference": max(r["relative_difference"] for r in rows),
                "converged": max(r["relative_difference"] for r in rows) < 1e-3,
            })

    record = {
        "schema": "NSC-RN-CLAIM-CONVERGENCE-v1",
        "api_version": "nsc-rn-claim-convergence-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Whether the retained fraction, which had only ever been read off a coarse "
            "grid, survives refinement. It is compared at matched times, because it grows "
            "with time and a mismatched comparison would measure the clock."
        ),
        "scope_note": (
            "A numerical convergence check on one quantity. It does not measure the actual "
            "universe, and it says nothing about the late-time floor or the relaxation "
            "fraction, which need a longer run than is currently recorded at fine "
            "resolution."
        ),
        "worst_relative_difference": max((c["worst_relative_difference"]
                                         for c in comparisons), default=None),
        "converges_at_every_edge": all(c["converged"] for c in comparisons),
        "reading_the_verdict": (
            "The strict threshold is one part in a thousand. A result that misses it is "
            "drifting, not wrong by much, and the magnitude is reported so that a False "
            "here is not read as the thirty percent the gain moved. For comparison the "
            "compression gain moved by that much under the same refinement."
        ),
        "comparisons": comparisons,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    for comparison in comparisons:
        print(f"  {comparison['coarse']} vs {comparison['fine']}")
        for row in comparison["rows"]:
            print(f"     t={row['coarse_t']:.1f}  coarse {row['coarse_retained']:.4f}"
                  f"  fine {row['fine_retained']:.4f}"
                  f"  rel {row['relative_difference']:.2e}"
                  f"  {'matched' if row['matched'] else 'CLOCK MISMATCH'}")
        print(f"     converged: {comparison['converged']}")
    print(f"  worst relative difference : {record['worst_relative_difference']:.2e}")
    print(f"  under one part in a thousand at every edge: "
          f"{record['converges_at_every_edge']}")
    print("  (compare: the compression gain moved by 3e-1 under the same refinement)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())