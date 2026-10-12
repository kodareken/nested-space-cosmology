"""Adoption checklist: the landing procedure, written before the landing.

When the station files appear, the adoption must follow one fixed
sequence, so the verdict enters the ledger exactly as pre-registered.
This record is that sequence, written while the march still integrates.

    python scripts/lab.py scripts/derive_nsc_t30_adoption_checklist.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-adoption-checklist-v1.json"
SCHEMA = "NSC-T30-ADOPTION-CHECKLIST-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "checklist": [
            "detect station-25/30 files in the campaign directory",
            "run scripts/derive_nsc_t30_extinction.py once",
            "read the verdict from nsc-t30-extinction-v1.json",
            "adopt the reading per nsc-t30-expectations-v1.json (no new framing)",
            "update the code-map corridor row with the landed verdict",
            "extend the owning-note part 15 with the closing sentence",
            "re-run both suites (conditional test must now pass, not skip)",
        ],
        "guardrails": [
            "never edit the campaign files; they are the campaign's sealed writes",
            "never re-run the extinction reducer if the record exists",
            "never let a verdict enter the synthesis before its record exists",
        ],
        "notes": [
            "written at round 114 while the march integrating continues",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_adoption_checklist.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_adoption_checklist.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  checklist steps:", len(record["checklist"]))
    print("  guardrails:", len(record["guardrails"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())