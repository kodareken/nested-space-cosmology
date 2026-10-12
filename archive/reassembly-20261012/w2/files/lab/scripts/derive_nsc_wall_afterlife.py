"""The wall's afterlife: the t=40 event and the single-force approach.

Round-73 read the wall's full post-20 morphology from the relaxation
payloads. This driver holds the recorded rows against the isotropy
tail to characterize the afterlife: the death window [10, 15], the
cross-observable t=40 event (isotropy rise, high balance, slope wiggle),
and the single-force relaxation rate from the balance deficit at t=40
and t=80. All arithmetic is over recorded bytes.

    python scripts/lab.py scripts/derive_nsc_wall_afterlife.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

WALL = LAB / "results" / "development" / "nsc-relaxation-wall-v1.json"
ISOTROPY = LAB / "results" / "development" / "nsc-isotropy-wall-join-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-afterlife-v1.json"
SCHEMA = "NSC-WALL-AFTERLIFE-v1"


def main():
    if not WALL.is_file() or not ISOTROPY.is_file():
        raise SystemExit("missing wall or isotropy record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    wall = json.loads(WALL.read_text(encoding="utf-8"))
    isotropy = json.loads(ISOTROPY.read_text(encoding="utf-8"))
    rows = wall["rows"]
    deaths = [(left["t_over_rm"], right["t_over_rm"])
              for left, right in zip(rows, rows[1:])
              if left["wall_count"] > 0 and right["wall_count"] == 0]
    deficits = {row["t_over_rm"]: 1.0 - row["balance_at_throat"] for row in rows}
    t20_wiggle = None
    by_t = {row["t_over_rm"]: row for row in rows}
    if all(t in by_t for t in (15.0, 20.0, 40.0)):
        t20_wiggle = {
            "slope_15": by_t[15.0]["throat_slope_inner"],
            "slope_20": by_t[20.0]["throat_slope_inner"],
            "slope_40": by_t[40.0]["throat_slope_inner"],
            "balance_40": by_t[40.0]["balance_at_throat"],
        }
    deficit_40 = deficits.get(40.0)
    deficit_80 = deficits.get(80.0)
    relaxation_rate = (
        None if deficit_40 is None or deficit_80 is None or deficit_40 <= 0
        else math.log(deficit_40 / deficit_80) / (80.0 - 40.0)
    )
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "death_window": deaths,
        "t40_event": {
            "isotropy_40": isotropy["tail"]["t40_isotropy"],
            "balance_40": deficits.get(40.0),
            "slope_wiggle": t20_wiggle,
        },
        "single_force_deficits": {"t40": deficit_40, "t80": deficit_80},
        "single_force_relaxation_rate_per_t": relaxation_rate,
        "deficit_halving_time": None if relaxation_rate is None
        else math.log(2.0) / relaxation_rate,
        "reading": (
            "the wall dies once, in the window [10, 15]; the t=40 event "
            "is a balance dip/wiggle missing from a pure relaxation, and "
            "the t=80 state sits inside a percent of the single-force limit"
        ),
        "notes": [
            "the relaxation rate is two-point over t=40..80",
        ],
        "input_sha256": {
            "lab/results/development/nsc-relaxation-wall-v1.json": chain.path_sha256(WALL),
            "lab/results/development/nsc-isotropy-wall-join-v1.json":
                chain.path_sha256(ISOTROPY),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_afterlife.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_afterlife.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  death windows:", deaths)
    print("  t40 event: isotropy {} balance {:.4f} slopes {} -> {} -> {}".format(
        record["t40_event"]["isotropy_40"], record["t40_event"]["balance_40"],
        t20_wiggle["slope_15"], t20_wiggle["slope_20"], t20_wiggle["slope_40"]))
    print("  single-force rate {:.4f}/t  halving {:.2f} t  deficit at 80: {:.4f}".format(
        relaxation_rate, record["deficit_halving_time"], deficit_80))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())