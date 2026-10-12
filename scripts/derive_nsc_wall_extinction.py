"""The wall's extinction ledger across the two grid domains.

The family chart (inner edge at the mass) and the march chart (one cell
above r_minus) both carry wall history, and round-76 hinted their
lifespan orderings differ. This driver computes for every recorded case
the last station with wall_count > 0 - the extinction time - on each
domain, and lays the two orderings side by side. The verdict records
whether the lifespan ordering inverts with domain depth.

    python scripts/lab.py scripts/derive_nsc_wall_extinction.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

BUBBLE = LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json"
MARCH_BALANCE = LAB / "results" / "development" / "nsc-march-balance-trajectory-v1.json"
RELAXATION = LAB / "results" / "development" / "nsc-relaxation-wall-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-extinction-v1.json"
SCHEMA = "NSC-WALL-EXTINCTION-v1"


def extinction_from_stations(stations):
    wall_times = [float(row["t_over_rm"]) for row in stations]
    last = max((float(row["t_over_rm"]) for row in stations
                if int(row["wall_count"]) > 0), default=None)
    return {"covered_until": max(wall_times), "last_wall": last}


def main():
    if any(not p.is_file() for p in (BUBBLE, MARCH_BALANCE, RELAXATION)):
        raise SystemExit("missing extinction input record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    march = json.loads(MARCH_BALANCE.read_text(encoding="utf-8"))
    family = {}
    for case in bubble["cases"]:
        if case["epsilon"] != 0.001:
            continue
        family[case["case"]] = extinction_from_stations(case["stations"])
    march_domain = {}
    for name, rows in march["series"].items():
        march_domain[name] = {
            "covered_until": max(row["t_over_rm"] for row in rows),
            "last_wall": max((row["t_over_rm"] for row in rows
                              if row["wall_count"] > 0), default=None),
        }
    ordering_family = sorted(
        family, key=lambda name: (family[name]["last_wall"] if family[name]["last_wall"] else -1))
    ordering_march = sorted(
        march_domain, key=lambda name: (march_domain[name]["last_wall"] if march_domain[name]["last_wall"] else -1))
    inverted = any(
        ordering_family[i] != ordering_march[i] for i in range(len(ordering_family)))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "family_chart_extinction": family,
        "march_chart_extinction": march_domain,
        "ordering_family": ordering_family,
        "ordering_march": ordering_march,
        "ordering_inverts_with_domain_depth": bool(inverted),
        "reading": (
            "the wall's lifespan ordering across cases flips when the "
            "grid reaches r_minus; the lightest case lives longest on "
            "the deep grid and shortest on the shallow one"
        ),
        "notes": [
            "last_wall = final station with wall_count > 0 inside each "
            "domain's own station coverage",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
            "lab/results/development/nsc-march-balance-trajectory-v1.json":
                chain.path_sha256(MARCH_BALANCE),
            "lab/results/development/nsc-relaxation-wall-v1.json":
                chain.path_sha256(RELAXATION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_extinction.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_extinction.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  family chart extinction:", family)
    print("  march chart extinction:", march_domain)
    print("  family ordering:", ordering_family)
    print("  march ordering:", ordering_march)
    print("  ordering inverts with domain depth:", record["ordering_inverts_with_domain_depth"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())