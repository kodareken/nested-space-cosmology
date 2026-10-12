"""The wall's post-20 morphology, read from the relaxation payloads.

Round-57 quoted the isotropy tail but missed that the relaxation
campaign directory holds the full station payloads at t = 5..80 for
mu=1.040, complete with radial_pressure, angular_pressure and R4 - the
exact ingredients of the wall machinery. This driver runs the
gradient-bubble reducer on every recorded station, delivering the wall
counts, throat balances and wall radii across the whole campaign,
including the t=40 and t=80 epochs that were the wall programme's
blind spot.

    python scripts/lab.py scripts/derive_nsc_relaxation_wall.py
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402
from recursive_horizons import nsc_rn_gradient_bubble as bubble  # noqa: E402

RELAXATION = LAB / "results" / "development" / "nsc-rn-relaxation-v1" / "mu-1.040-eps-0.001"
OUTPUT = LAB / "results" / "development" / "nsc-relaxation-wall-v1.json"
SCHEMA = "NSC-RELAXATION-WALL-v1"


def main():
    if not RELAXATION.is_dir():
        raise SystemExit("missing relaxation campaign directory")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    files = sorted(
        glob.glob(str(RELAXATION / "station-*.npz")),
        key=lambda path: float(path.rsplit("station-", 1)[1][:-4]),
    )
    rows = []
    for path in files:
        t = float(path.rsplit("station-", 1)[1][:-4])
        payload = np.load(path)
        profile = bubble.gradient_bubble_profile(payload)
        first_wall = profile["r_first_wall"]
        if isinstance(first_wall, float) and not np.isfinite(first_wall):
            first_wall = None
        rows.append({
            "t_over_rm": t,
            "balance_at_throat": profile["balance_at_throat"],
            "balance_min": profile["balance_min"],
            "wall_count": profile["wall_count"],
            "r_first_wall": first_wall,
            "throat_slope_inner": profile["throat_slope_inner"],
        })
    beyond = [row for row in rows if row["t_over_rm"] > 20.0]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "beyond_t20": beyond,
        "wall_survives_to_t40": bool(
            any(row["t_over_rm"] == 40.0 and row["wall_count"] > 0 for row in beyond)),
        "wall_survives_to_t80": bool(
            any(row["t_over_rm"] == 80.0 and row["wall_count"] > 0 for row in beyond)),
        "notes": [
            "the payloads are the recorded relaxation campaign bytes; the "
            "wall reducer is the same one the bubble record used",
        ],
        "input_sha256": {},
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_relaxation_wall.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_relaxation_wall.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        first_text = "none" if row["r_first_wall"] is None else "{:.3f}".format(row["r_first_wall"])
        print("  t={:<5} balance {:.4f} walls {} first {} slope {:.2f}".format(
            row["t_over_rm"], row["balance_at_throat"], row["wall_count"],
            first_text, row["throat_slope_inner"]))
    print("  survives to t40: {}  to t80: {}".format(
        record["wall_survives_to_t40"], record["wall_survives_to_t80"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())