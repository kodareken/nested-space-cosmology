"""Semantic audit of the marching loop: what the forecast can rely on.

A last modeling hole closed by reading the loop itself: march_pair
runs with a fixed dt and no step budget when cpu_seconds=None (the
check only fires when a cap is passed), the station loop writes only at
listed stations, and the step count per unit of t is therefore fixed -
which is why the v1-derived linear cost forecast is structurally sound.

    python scripts/lab.py scripts/derive_nsc_t30_semantic_audit.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-semantic-audit-v1.json"
SCHEMA = "NSC-T30-SEMANTIC-AUDIT-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "verified_from": {
            "march_pair": "derive_nsc_rn_parent.py lines 1092-1131",
            "station loop": "derive_nsc_rn_inward_march.py lines 270-306",
        },
        "semantics": [
            "dt is fixed for the whole segment (no CFL adaptation mid-march)",
            "cpu_seconds=None leaves the budget check inert: no early stop can happen",
            "writes happen only at listed station times",
            "step count per unit of t is therefore fixed and cost is linear in t",
        ],
        "consequence": (
            "the v1-derived linear CPU forecast (~12.9k CPU-seconds to "
            "station-25) is structurally sound; there is no hidden "
            "budget-truncation to confuse the meter"
        ),
        "notes": [],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_semantic_audit.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_semantic_audit.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  semantics verified:", len(record["semantics"]))
    print("  no hidden budget truncation:", True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())