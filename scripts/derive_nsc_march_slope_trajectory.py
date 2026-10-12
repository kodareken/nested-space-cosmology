"""The march throat slope's per-station trajectory and the sign flip's time.

Round-66 joined the march's last-station throat slope to the regime
boundary (+36.6, +48.3, -23.3). This driver recomputes the SAME slope
at every recorded march station (t = 5, 10, 15, 20) per case, so the
mu=1.040 sign crossing gets a time window and the light cases'
temporal shape is visible.

    python scripts/lab.py scripts/derive_nsc_march_slope_trajectory.py
"""
from __future__ import annotations

import glob
import json
import math
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402
from recursive_horizons import nsc_rn_inward as inward  # noqa: E402

MARCH = LAB / "results" / "development" / "nsc-rn-inward-march-v1"
OUTPUT = LAB / "results" / "development" / "nsc-march-slope-trajectory-v1.json"
SCHEMA = "NSC-MARCH-SLOPE-TRAJECTORY-v1"
CASES = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")


def main():
    if not MARCH.is_dir():
        raise SystemExit("missing march campaign directory")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    series = {}
    for name in CASES:
        files = sorted(
            glob.glob(str(MARCH / name / "station-*.npz")),
            key=lambda path: float(path.rsplit("station-", 1)[1][:-4]),
        )
        rows = []
        mass = float(name.split("-eps")[0].split("mu-")[1])
        r_minus = mass - math.sqrt(mass * mass - 1.0)
        for path in files:
            t = float(path.rsplit("station-", 1)[1][:-4])
            payload = np.load(path)
            decision = inward.inner_region_decision(
                np.asarray(payload["radius"], dtype=float),
                np.asarray(payload["radial_pressure"], dtype=float),
                r_minus,
            )
            rows.append({
                "t_over_rm": t,
                "throat_slope_inner": decision["throat_slope_inner"],
                "inner_region_slope_min": decision["inner_region_slope_min"],
            })
        series[name] = rows
    crossing = []
    for left, right in zip(series["mu-1.040-eps-0.001"], series["mu-1.040-eps-0.001"][1:]):
        if left["throat_slope_inner"] * right["throat_slope_inner"] < 0.0:
            crossing.append([left["t_over_rm"], right["t_over_rm"]])
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "series": series,
        "mu_1040_sign_crossing_windows": crossing,
        "reading": (
            "the march's own throat slope, recomputed per recorded station; "
            "the mu=1.040 sign flip is timed by its crossing window"
        ),
        "notes": [
            "same decision helper the march assessment itself used",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_march_slope_trajectory.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_march_slope_trajectory.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, rows in series.items():
        print("  {:<20} ".format(name) + "  ".join(
            "t{}:{:+.1f}".format(row["t_over_rm"], row["throat_slope_inner"])
            for row in rows))
    print("  mu=1.040 sign crossing windows:", crossing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())