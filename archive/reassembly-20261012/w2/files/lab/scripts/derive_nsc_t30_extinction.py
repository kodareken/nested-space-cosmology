"""The light wall's fate beyond t=20: read the t25/t30 march stations.

Round-88 launched the deep-domain march of mu=1.002 to t=30. This
driver reduces its new stations through the wall machinery once they
land and settles the round-77 open question: does the light wall die
between 20 and 30, and which window?

    python scripts/lab.py scripts/derive_nsc_t30_extinction.py
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

CAMPAIGN = LAB / "results" / "development" / "nsc-rn-inward-march-t30-v1" \
    / "mu-1.002-eps-0.001"
OUTPUT = LAB / "results" / "development" / "nsc-t30-extinction-v1.json"
SCHEMA = "NSC-T30-EXTINCTION-v1"


def main():
    if not CAMPAIGN.is_dir():
        raise SystemExit("t30 campaign not landed yet")
    files = sorted(glob.glob(str(CAMPAIGN / "station-*.npz")))
    if not files:
        raise SystemExit("no station payloads yet in the t30 campaign")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    rows = []
    for path in files:
        t = float(path.rsplit("station-", 1)[1][:-4])
        payload = np.load(path)
        profile = bubble.gradient_bubble_profile(payload)
        rows.append({
            "t_over_rm": t,
            "balance_at_throat": profile["balance_at_throat"],
            "wall_count": profile["wall_count"],
            "r_first_wall": (None if not np.isfinite(profile["r_first_wall"])
                             else profile["r_first_wall"]),
        })
    rows.sort(key=lambda row: row["t_over_rm"])
    t25 = next((row for row in rows if abs(row["t_over_rm"] - 25.0) < 1e-6), None)
    t30 = next((row for row in rows if abs(row["t_over_rm"] - 30.0) < 1e-6), None)
    wall_at_25 = t25 is not None and t25["wall_count"] > 0
    wall_at_30 = t30 is not None and t30["wall_count"] > 0
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "wall_at_25": bool(wall_at_25),
        "wall_at_30": bool(wall_at_30),
        "verdict": (
            "WALL_SURVIVES_TO_30" if wall_at_30
            else ("DIES_IN_25_30" if wall_at_25
                  else "DIES_IN_20_25" if t25 is not None
                  else "NOT_YET_FULLY_LANDED")
        ),
        "reading": (
            "the light wall's extinction window on the deep grid, settled "
            "from the freshly marched stations"
        ),
        "notes": [
            "run only after the t30 campaign has landed both stations",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-inward-march-t30-v1/": None,
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_extinction.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_extinction.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  t={:<4} balance {:.4f} walls {} first {!s}".format(
            row["t_over_rm"], row["balance_at_throat"], row["wall_count"],
            None if row["r_first_wall"] is None else "{:.3f}".format(row["r_first_wall"])))
    print("  wall at 25: {}  at 30: {}  verdict: {}".format(
        record["wall_at_25"], record["wall_at_30"], record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())