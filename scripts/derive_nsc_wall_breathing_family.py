"""Family breathing: rise/fall asymmetry across all three mu cases.

Round-23 measured the mu=1.040 wall breathing (rise +0.131, fall -0.242
r_+/t). This driver repeats the measurement for mu=1.002 and mu=1.010
and re-derives mu=1.040 inside the same file, asking whether the
rise/fall asymmetry and the radius-inflow slaving recur, or whether they
are another property that locks only at mu=1.040.

    python scripts/lab.py scripts/derive_nsc_wall_breathing_family.py
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
OUTPUT = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
SCHEMA = "NSC-WALL-BREATHING-FAMILY-v1"
CASES = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")


def pearson(x, y):
    n = len(x)
    if n < 3:
        return None, n
    xm = sum(x) / n
    ym = sum(y) / n
    sxy = sum((a - xm) * (b - ym) for a, b in zip(x, y))
    sxx = sum((a - xm) ** 2 for a in x)
    syy = sum((b - ym) ** 2 for b in y)
    if sxx <= 0.0 or syy <= 0.0:
        return None, n
    return sxy / (sxx * syy) ** 0.5, n


def main():
    if not BUBBLE.is_file():
        raise SystemExit("missing bubble record: " + str(BUBBLE))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    cases = {}
    for name in CASES:
        wall_case = next(c for c in bubble["cases"] if c["case"] == name)
        files = sorted(glob.glob(str(FAMILY / name / "station-*.json")))
        horizons = {}
        for path in files:
            if Path(path).stem == "station-failure":
                continue
            record = json.loads(Path(path).read_text(encoding="utf-8"))
            sample = record["sample"]
            horizons[round(float(sample["t_over_rm"]), 6)] = float(sample["coupled_trapping"])
        r_plus = sum(horizons.values()) / len(horizons) if horizons else None
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
            dt = t1 - t0
            if dt <= 0.0:
                continue
            velocity = (r1 - r0) / dt
            (rising if velocity > 0.0 else falling).append(velocity)
        mean = lambda items: (sum(items) / len(items)) if items else None
        cases[name] = {
            "epochs": len(epochs),
            "radial_range": [min(r for _, r in epochs), max(r for _, r in epochs)]
            if epochs else None,
            "rise_epochs": len(rising),
            "fall_epochs": len(falling),
            "mean_rise_rate": mean(rising),
            "mean_fall_rate": mean(falling),
            "asymmetry_fall_over_rise": (
                None if (mean(rising) and mean(falling)) is None
                else abs(mean(falling)) / abs(mean(rising))
                if mean(rising) != 0.0 else None
            ),
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "cases": cases,
        "notes": [
            "same definition as nsc-wall-breathing-v1, applied per case",
            "wall epochs come from the recorded bubble v5 flags",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_breathing_family.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_breathing_family.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in cases.items():
        print(
            "  {:<20} epochs {:>3} range {!s} rise {: .4f} (n={}) fall {: .4f} (n={}) asymmetry {!s}".format(
                name, entry["epochs"],
                "none" if entry["radial_range"] is None
                else "[{:.3f},{:.3f}]".format(*entry["radial_range"]),
                0.0 if entry["mean_rise_rate"] is None else entry["mean_rise_rate"],
                entry["rise_epochs"],
                0.0 if entry["mean_fall_rate"] is None else entry["mean_fall_rate"],
                entry["fall_epochs"],
                "none" if entry["asymmetry_fall_over_rise"] is None
                else "{:.2f}".format(entry["asymmetry_fall_over_rise"]),
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())