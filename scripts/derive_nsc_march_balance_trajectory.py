"""Wall flags on the horizon-reaching march stations.

The march campaign (the v2/two-ends domain) never had its wall flags
read: round-73's relaxation payloads carried them, and these do too.
This driver runs the gradient-bubble reducer on all twelve march
stations (three cases, t = 5..20), giving the wall count and throat
balance trajectory on the horizon-reaching grid and the t=10 wall epoch
as seen by the wall machinery itself.

    python scripts/lab.py scripts/derive_nsc_march_balance_trajectory.py
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

MARCH = LAB / "results" / "development" / "nsc-rn-inward-march-v1"
OUTPUT = LAB / "results" / "development" / "nsc-march-balance-trajectory-v1.json"
SCHEMA = "NSC-MARCH-BALANCE-TRAJECTORY-v1"
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
                "wall_count": profile["wall_count"],
                "r_first_wall": first_wall,
            })
        series[name] = rows
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "series": series,
        "reading": (
            "the march grid reaches one cell above r_minus; these are the "
            "first wall flags on that domain"
        ),
        "notes": [
            "same reducer as the bubble record; the march payloads carry "
            "the pressures and were written by the march campaign",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_march_balance_trajectory.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_march_balance_trajectory.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, rows in series.items():
        print("  {:<20} ".format(name) + "  ".join(
            "t{}:B{:.2f} w{}".format(row["t_over_rm"], row["balance_at_throat"],
                                     row["wall_count"]) for row in rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())