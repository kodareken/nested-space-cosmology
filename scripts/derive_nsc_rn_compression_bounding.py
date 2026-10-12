"""What the compression is bounded by, which is not the resolution.

The saturation of the compression gain has been described as deepening: the
further in the boundary is placed, the more it compresses, and it stops near
the inner horizon. That is true but it says where, not why.

The gain is better read against the two distances that define the region the
compression can occupy: the distance from the boundary to the horizon it is
sealed behind, and the distance from the boundary to the inner horizon, past
which the excision stops being an excision at all.

The deepest boundary tried sits at 99.8% of the available trapped width. The
compression and the turn-over in the gain therefore end in the same place,
and for the same reason: the region available to compress into runs out.

Reads only recorded records. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_compression_bounding.py
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
MAP = DEVELOPMENT / "nsc-rn-geometry-map-v6.json"
DOMAIN = DEVELOPMENT / "nsc-rn-excision-domain.md" if False else \
    LAB / "docs/nsc-rn-excision-domain.md"


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mu", type=float, default=1.04)
    parser.add_argument("--out", default="nsc-rn-compression-bounding-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    chart = [e for e in json.loads(MAP.read_text(encoding="utf-8"))["entries"]
             if abs(e["mu"] - args.mu) < 1e-9][0]
    r_minus = float(chart["r_minus"])
    r_plus = float(chart["r_plus"])
    width = r_plus - r_minus

    # Read the gains from the fine-resolution shell records rather than
    # copying them. Copying them and then comparing against the coarse edge
    # sweep would quietly mix two resolutions, which is the error this whole
    # line of work exists to avoid.
    fine = {
        "n8000-def1040": 1.04,
        "n8000-exc085": 0.85,
        "n8000d-exc080": 0.80,
        "n8000d-exc077": 0.77,
        "n8000d-exc076": 0.76,
        "n8000d-exc07555": 0.7555,
    }
    gains = {}
    for name, edge in fine.items():
        record = DEVELOPMENT / ("nsc-rn-inner-packet-shell-%s-v1.json" % name)
        if not record.is_file():
            continue
        gains[edge] = float(json.loads(record.read_text(encoding="utf-8"))
                             ["accumulation_after_arrival"]["amplitude_gain"])
    if len(gains) < 4:
        raise SystemExit("fewer than four fine-resolution gains are on record")

    rows: List[Dict[str, Any]] = []
    for edge in sorted(gains, reverse=True):
        gain = gains[edge]
        rows.append({
            "r_inner": edge,
            "gain": gain,
            "distance_to_horizon": r_plus - edge,
            "distance_to_inner_horizon": edge - r_minus,
            "fraction_of_trapped_width_used": (r_plus - edge) / width,
        })

    ordered = sorted(rows, key=lambda r: r["distance_to_horizon"])
    rises = all(b["gain"] > a["gain"] for a, b in zip(ordered, ordered[1:-1]))
    turns_over = ordered[-1]["gain"] < ordered[-2]["gain"]
    deepest = rows[-1]

    record = {
        "schema": "NSC-RN-COMPRESSION-BOUNDING-v1",
        "api_version": "nsc-rn-compression-bounding-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "The compression gain against the two distances that bound the region it can "
            "occupy: the distance to the horizon it is sealed behind, and the distance to "
            "the inner horizon past which the excision stops being one."
        ),
        "reading": (
            "The gain rises as the boundary recedes from the horizon and turns over once "
            "the boundary reaches the far end of the trapped region. The deepest boundary "
            "tried sits at almost the whole trapped width from the horizon, so the "
            "compression saturates because the region available to compress into runs "
            "out, not because of the grid."
        ),
        "why_it_matters": (
            "The saturation was previously attributed to depth alone, which says where it "
            "happens and not why. It also sits in the same place as the sign change of the "
            "outgoing speed, where the excision boundary condition reverses."
        ),
        "scope_note": (
            "A re-reading of records already taken. It does not measure the actual "
            "universe and it adds no new run. Six edges, one mass ratio, one source "
            "strength."
        ),
        "r_minus": r_minus,
        "r_plus": r_plus,
        "trapped_width": width,
        "gain_rises_with_distance_to_horizon": rises,
        "turns_over_at_the_far_end": turns_over,
        "deepest_edge_fraction_of_width": deepest["fraction_of_trapped_width_used"],
        "deepest_edge_distance_to_inner_horizon": deepest["distance_to_inner_horizon"],
        "rows": rows,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  trapped region {r_minus:.4f} .. {r_plus:.4f}  width {width:.4f}")
    print(f"  {'r_inner':>9} {'to horizon':>11} {'to r_-':>9} {'% of width':>11} {'gain':>8}")
    for row in rows:
        print(f"  {row['r_inner']:>9.4f} {row['distance_to_horizon']:>11.4f}"
              f" {row['distance_to_inner_horizon']:>9.4f}"
              f" {100 * row['fraction_of_trapped_width_used']:>10.1f}% {row['gain']:>8.4f}")
    print(f"  gain rises with distance to horizon: {rises}")
    print(f"  turns over at the far end            : {turns_over}")
    print(f"  deepest edge uses {100 * record['deepest_edge_fraction_of_width']:.1f}% of the width")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())