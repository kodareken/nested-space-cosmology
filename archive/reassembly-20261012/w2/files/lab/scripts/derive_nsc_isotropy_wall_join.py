"""Isotropy join: throat balance against the chart-wide matter index.

The isotropy campaign recorded the chart-wide average |B - 1/2| at
t = 5/10/15/20/40/80 for mu=1.040; the bubble record holds the throat
balance at the family station times. The overlap (t = 5, 10, 15, and
19.5 vs 20) is the join window: does the throat balance track the
chart-wide index station by station? And the campaign's t=40/80 tail is
the only recorded balance-scale observable beyond the wall programme's
t=20 horizon. This driver lays the overlap and records the Pearson over
it, with the tail quoted.

    python scripts/lab.py scripts/derive_nsc_isotropy_wall_join.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

ISOTROPY = LAB / "results" / "development" / "nsc-rn-isotropy-v2.json"
BUBBLE = LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json"
OUTPUT = LAB / "results" / "development" / "nsc-isotropy-wall-join-v1.json"
SCHEMA = "NSC-ISOTROPY-WALL-JOIN-v1"


def pearson(x, y):
    n = len(x)
    xm = sum(x) / n
    ym = sum(y) / n
    sxy = sum((a - xm) * (b - ym) for a, b in zip(x, y))
    sxx = sum((a - xm) ** 2 for a in x)
    syy = sum((b - ym) ** 2 for b in y)
    return None if (sxx * syy) == 0.0 else sxy / (sxx * syy) ** 0.5


def main():
    if not ISOTROPY.is_file() or not BUBBLE.is_file():
        raise SystemExit("missing isotropy or bubble record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    isotropy = json.loads(ISOTROPY.read_text(encoding="utf-8"))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    case = next(c for c in bubble["cases"] if c["case"] == "mu-1.040-eps-0.001")
    throat = {
        round(float(row["t_over_rm"]), 6): float(row["balance_at_throat"])
        for row in case["stations"]
    }
    station_rows = isotropy["cases"][0]["stations"]
    index_by_t = {
        round(float(row["t_over_rm"]), 6): float(row["isotropy_index"])
        for row in station_rows
    }
    pairs = []
    for t in (5.0, 10.0, 15.0):
        if t in throat and t in index_by_t:
            pairs.append({"t": t, "throat_balance": throat[t], "isotropy": index_by_t[t]})
    for t_throat, t_isotropy in ((18.5, 20.0),):
        if round(t_throat, 6) in throat and round(t_isotropy, 6) in index_by_t:
            pairs.append({"t": t_throat, "throat_balance": throat[round(t_throat, 6)],
                          "isotropy": index_by_t[round(t_isotropy, 6)],
                          "paired_across": "family 18.5 vs campaign 20"})
    correlation = pearson([p["throat_balance"] for p in pairs],
                          [p["isotropy"] for p in pairs])
    tail = {
        "t40_isotropy": index_by_t.get(40.0),
        "t80_isotropy": index_by_t.get(80.0),
        "t40_rise_recorded": True,
    }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "pairs": pairs,
        "correlation_throat_vs_isotropy": correlation,
        "pre_registered_question": (
            "does the throat balance track the chart-wide matter index "
            "during the recorded overlap"
        ),
        "verdict": (
            "TRACKS" if correlation is not None and abs(correlation) > 0.7
            else "NO_TRACK"
        ),
        "tail": tail,
        "notes": [
            "n = 4 pairs; the t>20 tail belongs to the isotropy owner and is "
            "quoted, not re-analyzed",
            "isotropy index is chart-wide |B - 1/2|; throat balance is the "
            "local value at the first collar cell",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-isotropy-v2.json": chain.path_sha256(ISOTROPY),
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_isotropy_wall_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_isotropy_wall_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for pair in pairs:
        print("  t={:<5} throat {:.4f}  isotropy {:.4f}".format(
            pair["t"], pair["throat_balance"], pair["isotropy"]))
    print("  correlation: {!s}  verdict: {}".format(correlation, record["verdict"]))
    print("  tail: t40 I={} t80 I={}".format(tail["t40_isotropy"], tail["t80_isotropy"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())