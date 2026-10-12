"""Epsilon-crossing test of the wall clock pair.

The round-26/30 power laws (rise ~ gap^-1/2, fall ~ gap^-1/4) were
measured on the eps=0.001 trio. The family also holds the eps=0.010
trio with independent wall-epoch records. This driver evaluates the SAME
pre-registered products on the 0.010 trio and records (a) whether the
two laws recur, (b) whether the clock constants scale with the packet
strength epsilon. Nothing is fitted.

    python scripts/lab.py scripts/derive_nsc_wall_clock_epscross.py
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
from recursive_horizons import nsc_rn_source as source  # noqa: E402

FAMILY = LAB / "results" / "development" / "nsc-rn-neutral-family-v1"
BUBBLE = LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-clock-epscross-v1.json"
SCHEMA = "NSC-WALL-CLOCK-EPSCROSS-v1"
CASES = (
    "mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001",
    "mu-1.002-eps-0.010", "mu-1.010-eps-0.010", "mu-1.040-eps-0.010",
)


def stats_pair(names, entries):
    gaps = {name: entries[name]["gap"] for name in names}
    half = [abs(entries[name]["rise_rate"]) * gaps[name] ** 0.5 for name in names]
    quarter = [abs(entries[name]["fall_rate"]) * gaps[name] ** 0.25 for name in names]
    raw_rise = [abs(entries[name]["rise_rate"]) for name in names]
    raw_fall = [abs(entries[name]["fall_rate"]) for name in names]
    def spread(values):
        mean = sum(values) / len(values)
        return (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5 / mean
    return {
        "raw_rise_spread": spread(raw_rise),
        "raw_fall_spread": spread(raw_fall),
        "rise_half_spread": spread(half),
        "fall_quarter_spread": spread(quarter),
    }


def main():
    if not BUBBLE.is_file():
        raise SystemExit("missing bubble record: " + str(BUBBLE))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    entries = {}
    for name in CASES:
        wall_case = next(c for c in bubble["cases"] if c["case"] == name)
        files = sorted(glob.glob(str(FAMILY / name / "station-*.json")))
        horizons = []
        for path in files:
            if Path(path).stem == "station-failure":
                continue
            record = json.loads(Path(path).read_text(encoding="utf-8"))
            sample = record["sample"]
            horizons.append(float(sample["coupled_trapping"]))
        r_plus = sum(horizons) / len(horizons)
        epochs = []
        for wall_row in wall_case["stations"]:
            first = wall_row.get("r_first_wall")
            if not (isinstance(first, float) and first == first):
                continue
            t = round(float(wall_row["t_over_rm"]), 6)
            epochs.append((t, float(first) / r_plus))
        epochs.sort()
        rising = []
        falling = []
        for (t0, r0), (t1, r1) in zip(epochs, epochs[1:]):
            if t1 - t0 <= 0.0:
                continue
            velocity = (r1 - r0) / (t1 - t0)
            (rising if velocity > 0.0 else falling).append(velocity)
        mean = lambda items: (sum(items) / len(items)) if items else None
        masses = {"mu-1.002": 1.002, "mu-1.010": 1.010, "mu-1.040": 1.040}
        mu_key = name.split("-eps")[0]
        mass = masses[mu_key]
        model = source.action_reference(mass)
        gap = float(model.horizons["r_plus"]) - float(model.horizons["r_minus"])
        entries[name] = {
            "epochs": len(epochs),
            "gap": gap,
            "rise_rate": mean(rising),
            "fall_rate": mean(falling),
            "rise_epochs": len(rising),
            "fall_epochs": len(falling),
        }
    triple_0001 = CASES[:3]
    triple_0010 = CASES[3:]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "eps_0001": stats_pair(triple_0001, entries),
        "eps_0010": stats_pair(triple_0010, entries),
        "cases": entries,
        "verdict": (
            "EPS_CROSS_HOLDS"
            if stats_pair(triple_0010, entries)["fall_quarter_spread"] < 0.08
            and stats_pair(triple_0010, entries)["rise_half_spread"] < 0.08
            else "EPS_CROSS_BROKEN"
        ),
        "reading": (
            "the same pre-registered power pair is evaluated on the strong "
            "packet trio; gapped spreads below 8% mean the laws recur"
        ),
        "notes": [
            "epsilon multiplies the packet energy target tenfold; the wall "
            "flags and station files are the recorded ones",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-gradient-bubble-v5.json":
                chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_clock_epscross.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_clock_epscross.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name in CASES:
        entry = entries[name]
        print("  {:<20} epochs {:>2} rise {: .4f} fall {: .4f}".format(
            name, entry["epochs"],
            0.0 if entry["rise_rate"] is None else entry["rise_rate"],
            0.0 if entry["fall_rate"] is None else entry["fall_rate"],
        ))
    for label, stat in (("eps=0.001", record["eps_0001"]), ("eps=0.010", record["eps_0010"])):
        print("  {:<10} raw rise {:.1%} -> half {:.1%}   raw fall {:.1%} -> quarter {:.1%}".format(
            label, stat["raw_rise_spread"], stat["rise_half_spread"],
            stat["raw_fall_spread"], stat["fall_quarter_spread"],
        ))
    print("  verdict:", record["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())