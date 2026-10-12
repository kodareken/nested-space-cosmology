"""Programme index v2: the round-99 state.

The round-84 index was sealed at the then-current state. Rounds 85-99
added six records and a conditional test; this v2 reissues the index
with the fuller citation set, while v1 stays untouched. Every cited
path is existence-checked before sealing.

    python scripts/lab.py scripts/derive_nsc_programme_index_v2.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-programme-index-v2.json"
SCHEMA = "NSC-PROGRAMME-INDEX-v2"

ADDITIONS = {
    "ledger and governance": [
        "nsc-ledger-audit-v1.json",
        "nsc-programme-index-v1.json",
    ],
    "domain law deep rates": [
        "nsc-march-wall-rates-v1.json",
    ],
    "frontier supplements": [
        "nsc-frontier-supplement-v1.json",
    ],
    "t30 campaign instruments": [
        "nsc-t30-launch-log-v1.json",
        "nsc-t30-expectations-v1.json",
    ],
}


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    development = LAB / "results" / "development"
    v1 = development / "nsc-programme-index-v1.json"
    if not v1.is_file():
        raise SystemExit("missing v1 index")
    base = json.loads(v1.read_text(encoding="utf-8"))
    missing = []
    held = {}
    for group, records in ADDITIONS.items():
        for name in records:
            if not (development / name).is_file():
                missing.append(name)
        held[group] = records
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "supersedes": "nsc-programme-index-v1.json",
        "supersedes_reason":
            "v1 was sealed at round 84; v2 carries the round-85..99 records",
        "v1_parts_kept": base["part_count"],
        "v1_cited_records": sum(len(part["records"]) for part in base["parts"]),
        "additions": held,
        "additional_records_cited": sum(len(records) for records in held.values()),
        "missing_records": missing,
        "index_complete": not missing,
        "notes": [
            "v1 remains on disk untouched; this v2 stands beside it",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_programme_index_v2.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_programme_index_v2.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  v1 parts kept:", record["v1_parts_kept"])
    print("  additions:", record["additional_records_cited"], "in",
          ", ".join(held))
    print("  missing:", missing if missing else "none - index complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())