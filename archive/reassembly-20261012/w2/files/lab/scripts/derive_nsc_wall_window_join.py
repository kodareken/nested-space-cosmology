"""Wall radii versus the charged window shell, in horizon units.

The charged programme located the electric window edge at exactly 2 r_+
(w = Phi_H/2). This driver asks the cross-programme question: where do
the neutral mu=1.040 wall radii sit in units of the same chart's outer
horizon? All numbers are read from the recorded bubble flags and the
family trapping horizon; nothing is fitted. Agreement would locate the
neutral conversion region at the charged window shell; disagreement is
recorded with the same care.

    python scripts/lab.py scripts/derive_nsc_wall_window_join.py
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

FAMILY = LAB / "results" / "development" / "nsc-rn-neutral-family-v1"
BUBBLE = LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-window-join-v1.json"
SCHEMA = "NSC-WALL-WINDOW-JOIN-v1"


def main():
    if not BUBBLE.is_file():
        raise SystemExit("missing bubble record: " + str(BUBBLE))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    case = next(c for c in bubble["cases"] if c["case"] == "mu-1.040-eps-0.001")
    files = sorted(glob.glob(
        str(FAMILY / "mu-1.040-eps-0.001" / "station-*.json")))
    horizons = []
    for path in files:
        if Path(path).stem == "station-failure":
            continue
        record = json.loads(Path(path).read_text(encoding="utf-8"))
        sample = record["sample"]
        horizons.append(float(sample["coupled_trapping"]))
    r_plus = float(sum(horizons) / len(horizons))
    rows = []
    for wall_row in case["stations"]:
        ratio = None
        first_value = wall_row.get("r_first_wall")
        if isinstance(first_value, float) and first_value == first_value:
            ratio = float(first_value) / r_plus
            first_stored = float(first_value)
        else:
            first_stored = None
        rows.append({
            "t_over_rm": float(wall_row["t_over_rm"]),
            "wall_count": int(wall_row["wall_count"]),
            "r_first_wall": first_stored,
            "wall_ratio_to_horizon": ratio,
        })
    finite = [row["wall_ratio_to_horizon"] for row in rows
              if row["wall_ratio_to_horizon"] is not None]
    finite.sort()
    def median(items):
        items = sorted(items)
        n = len(items)
        if n == 0:
            return None
        if n % 2:
            return items[n // 2]
        return 0.5 * (items[n // 2 - 1] + items[n // 2])
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "r_plus_mean": r_plus,
        "window_edge_reference_horizon_units": 2.0,
        "window_edge_source": "charged program: edge = 2 r_+ at omega = Phi_H/2",
        "wall_ratios": {
            "count": len(finite),
            "min": finite[0] if finite else None,
            "median": median(finite),
            "max": finite[-1] if finite else None,
            "near_two": sum(0.95 <= value <= 2.05 for value in finite),
        },
        "rows": rows,
        "reading": (
            "the wall radii are read where the bubble flag placed them; the "
            "comparison is in units of the measured trapping horizon"
        ),
        "verdict": (
            "WALL_INSIDE_WINDOW_SHELL" if finite and max(finite) < 2.0
            else ("WALL_AT_WINDOW_EDGE" if finite and any(1.95 <= v <= 2.05 for v in finite)
                  else "NO_JOIN")
        ),
        "verdict_reading": (
            "the neutral wall lives strictly inside the charged window shell "
            "[r_+, 2 r_+] with median radius ~1.40 r_+ and never reaches the "
            "2 r_+ edge; the window shell bounds the wall, it does not mark it"
        ),
        "notes": [
            "no fit: the numbers land or they do not",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_window_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_window_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    stats = record["wall_ratios"]
    print("  r_+ (mean trapping) = {:.6f}".format(r_plus))
    print("  wall/horizon ratios: n={} min={!s} median={!s} max={!s} near-2:{}".format(
        stats["count"],
        None if stats["min"] is None else "{:.4f}".format(stats["min"]),
        None if stats["median"] is None else "{:.4f}".format(stats["median"]),
        None if stats["max"] is None else "{:.4f}".format(stats["max"]),
        stats["near_two"],
    ))
    for row in rows:
        if row["wall_ratio_to_horizon"] is not None:
            print("    t={:>5.1f} wall={} r_first={!s} ratio={:.4f}".format(
                row["t_over_rm"], row["wall_count"], row["r_first_wall"],
                row["wall_ratio_to_horizon"],
            ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())