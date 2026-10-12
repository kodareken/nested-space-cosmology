"""Extinction extension: merge the t30 landing into the extinction ledger.

The round-77 extinction ledger is sealed. When the t30 campaign lands,
this driver (fired after the t30 extinction record exists) issues the
extension: the same table with the light case's extinction updated to
the landing's verdict. It refuses to run before the landing record
exists, so the sealed ledger cannot be pre-empted.

    python scripts/lab.py scripts/derive_nsc_extinction_extension.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-extinction-extension-v1.json"
SCHEMA = "NSC-EXTINCTION-EXTENSION-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    t30 = LAB / "results" / "development" / "nsc-t30-extinction-v1.json"
    if not t30.is_file():
        raise SystemExit("t30 extinction record not landed yet; refusing to pre-empt")
    landing = json.loads(t30.read_text(encoding="utf-8"))
    if landing["verdict"] not in ("WALL_SURVIVES_TO_30", "DIES_IN_25_30", "DIES_IN_20_25"):
        raise SystemExit("unexpected landing verdict: " + str(landing["verdict"]))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "extends": "nsc-wall-extinction-v1.json",
        "landing": {
            "campaign": "nsc-rn-inward-march-t30-v1",
            "verdict": landing["verdict"],
            "rows_t25_t30": [row for row in landing["rows"]
                             if row["t_over_rm"] >= 25.0],
        },
        "light_case_extinction_now": None,
        "reading": (
            "the light wall's extinction, updated by the landed stations; "
            "the sealed round-77 table remains untouched and this record "
            "stands beside it"
        ),
        "notes": [],
        "input_sha256": {
            "lab/results/development/nsc-t30-extinction-v1.json":
                chain.path_sha256(t30),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_extinction_extension.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_extinction_extension.py"
            ),
        },
    }
    if landing["verdict"] == "WALL_SURVIVES_TO_30":
        record["light_case_extinction_now"] = {
            "status": "alive through t=30",
            "next_named_question": "a longer campaign (t=35/40 or beyond)",
        }
    elif landing["verdict"] == "DIES_IN_25_30":
        record["light_case_extinction_now"] = {"status": "dies in (25, 30]"}
    else:
        record["light_case_extinction_now"] = {"status": "dies in (20, 25]"}
        record["notes"].append(
            "the wall dies between the march-v1 coverage (20) and the "
            "landed t=25 station")
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  landing verdict:", landing["verdict"])
    print("  light case now:", record["light_case_extinction_now"]["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())