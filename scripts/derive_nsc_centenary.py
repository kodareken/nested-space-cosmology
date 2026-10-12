"""Centenary annals: the programme at one hundred rounds.

The hundredth round's own record: the discovery ladder as it stands,
the ledger counts, and the standing state of the t30 campaign. A
benchmark entry, addition-only, for whoever reads the ledger at the
end of this programme's second hundred.

    python scripts/lab.py scripts/derive_nsc_centenary.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-centenary-v1.json"
SCHEMA = "NSC-CENTENARY-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "round": 100,
        "ladder": [
            "exterior channels refuted as the missing-structure carrier",
            "interior/excision funnel exact to 3e-15",
            "wall converter organism: concentration, lag, breathing",
            "gap-clock powers rise d^-1/2 fall d^-1/4; omega=0.386",
            "seven-quantity mu-blind scaffold",
            "fringe quantum 0.00131 and its ladder",
            "lifecycle: death at the driven phase's end, single-force tombstone",
            "domain law: grid-relative radius, amplitude ratio 7/17/40",
            "cross-arm bit-identity; ledger audit clean",
            "t30 campaign live: the light wall's fate beyond t=20 in flight",
        ],
        "ledger_at_hundred": {
            "suite": "101 passed, 1 conditional skip",
            "index": "v2 complete, 14 parts plus 6 additions",
            "open_cells": 2,
            "production_jobs": 1,
        },
        "standing_t30": {
            "campaign": "nsc-rn-inward-march-t30-v1",
            "launch": "2026-10-06T23:55:41",
            "first_station_landed": False,
            "diagnosis": (
                "the visibility pipe buffers all print output, and the "
                "launch shared the machine with a long suite run before "
                "round 96 killed the deadlock; both facts are recorded, "
                "the process accrues CPU, and no failure trail exists"
            ),
        },
        "notes": [
            "the annals are a benchmark, not a summary of record; the "
            "records themselves remain the papers of record",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_centenary.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_centenary.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  ladder entries:", len(record["ladder"]))
    print("  first station landed:", record["standing_t30"]["first_station_landed"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())