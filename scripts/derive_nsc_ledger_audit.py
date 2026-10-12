"""Ledger audit: every verdict-bearing record, re-read and re-checked.

The strand's numbers have been spot-checked round by round; this driver
runs the complete semantic pass over the verdict-bearing records - each
one re-opened, its recorded verdict compared with the ledger's
expectation, and any mismatch reported. This is the integrity check the
programme's size now requires.

    python scripts/lab.py scripts/derive_nsc_ledger_audit.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-ledger-audit-v1.json"
SCHEMA = "NSC-LEDGER-AUDIT-v1"

EXPECTATIONS = {
    "nsc-interior-balance-v1.json": ("verdict", "EXCISION_CLOSED outer_leaks_only"),
    "nsc-continuation-station-closure-v1.json": ("verdict", "LATE_CLOSURE_HOLDS"),
    "nsc-wall-clock-horizon-v1.json": ("none", None),
    "nsc-wall-clock-fall-v1.json": ("none", None),
    "nsc-wall-clock-epscross-v1.json": ("verdict", "EPS_CROSS_HOLDS"),
    "nsc-wall-asymmetry-derived-v1.json": ("none", None),
    "nsc-wall-clock-derivation-v1.json": ("none", None),
    "nsc-cosmic-band-location-v1.json": ("none", None),
    "nsc-time-awareness-v1.json": ("none", None),
    "nsc-matryoshka-fee-join-v1.json": ("none", None),
    "nsc-omega-star-identity-v1.json": ("verdict", "OMEGA_STAR_IDENTITY_HOLDS"),
    "nsc-charged-clock-horizon-v1.json": ("none", None),
    "nsc-charged-ladder-deficit-law-v1.json": ("verdict", "GAUSSIAN_APPROACH_HOLDS"),
    "nsc-charged-tilt-resolution-v1.json": ("verdict", "FRINGE_PHYSICS"),
    "nsc-tilt-organism-clock-v1.json": ("verdict", "ORGANISM_CLOCK_REJECTED"),
    "nsc-tilt-quantum-v1.json": ("verdict", "QUANTUM_IS_SINGLE_NODE"),
    "nsc-tilt-quantum-ladder-v1.json": ("verdict", "QUANTUM_LADDER_HOLDS"),
    "nsc-fringe-beta-join-v1.json": ("none", None),
    "nsc-two-ends-join-v1.json": ("verdict", "SHARPER_ROOM_SLOWER_RISE"),
    "nsc-two-ends-covariance-v1.json": ("verdict", "COVARIANCE_CLOSED_BY_GAP"),
    "nsc-two-ends-power-stability-v1.json": ("verdict", "POWER_STABLE"),
    "nsc-oneway-wall-join-v1.json": ("verdict", "TRANSIENT_BAND_CONTAINS_WALL_CLOCK"),
    "nsc-isotropy-wall-join-v1.json": ("verdict", "NO_TRACK"),
    "nsc-ladder-covariance-test-v1.json": ("verdict", "PREDICTIVE_PASS"),
    "nsc-continuation-trap-v1.json": ("none", None),
    "nsc-trap-timing-v1.json": ("verdict", "TRAP_LAGGED"),
    "nsc-production-preflight-v1.json": ("none", None),
    "nsc-production-gate-v1.json": ("none", None),
    "nsc-relaxation-wall-v1.json": ("none", None),
    "nsc-march-balance-trajectory-v1.json": ("none", None),
    "nsc-wall-extinction-v1.json": ("none", None),
    "nsc-cross-arm-identity-v1.json": ("none", None),
    "nsc-march-wall-radii-v1.json": ("none", None),
    "nsc-wall-domain-law-v1.json": ("none", None),
    "nsc-frontier-amendment-v1.json": ("none", None),
    "nsc-seam-invariance-v1.json": ("verdict", "PROFILES_COINCIDE"),
    "nsc-charged-seam-invariance-v1.json": ("verdict", "PROFILES_COINCIDE"),
    "nsc-interior-station-balance-v1.json": ("none", None),
    "nsc-wall-energy-channel-v1.json": ("none", None),
    "nsc-wall-ratio-regimes-v1.json": ("none", None),
    "nsc-lowband-placement-v1.json": ("verdict", "ALL_WALL_KINEMATICS_INSIDE_OUTWARD_CELL"),
    "nsc-constant-table-v1.json": ("none", None),
    "nsc-q-decomposition-v1.json": ("none", None),
    "nsc-q-analytic-envelope-v1.json": ("none", None),
    "nsc-programme-index-v1.json": ("none", None),
}


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    development = LAB / "results" / "development"
    rows = []
    mismatches = []
    missing = []
    for name, (field, expected) in EXPECTATIONS.items():
        path = development / name
        if not path.is_file():
            missing.append(name)
            continue
        body = json.loads(path.read_text(encoding="utf-8"))
        if field == "none":
            rows.append({"record": name, "check": "exists", "ok": True})
            continue
        actual = body.get(field)
        ok = bool(actual == expected)
        rows.append({
            "record": name,
            "check": field,
            "expected": expected,
            "actual": actual,
            "ok": ok,
        })
        if not ok:
            mismatches.append({"record": name, "expected": expected, "actual": actual})
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "checked_records": len(rows),
        "verdict_checks": sum(1 for row in rows if row["check"] != "exists"),
        "missing_records": missing,
        "mismatches": mismatches,
        "audit_clean": bool(not missing and not mismatches),
        "rows": rows,
        "notes": [
            "expectations were written from the records' own published "
            "results; the audit re-reads every file fresh",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_ledger_audit.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_ledger_audit.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  checked: {} records, {} verdicts".format(
        record["checked_records"], record["verdict_checks"]))
    print("  missing:", missing if missing else "none")
    print("  mismatches:", mismatches if mismatches else "none - audit clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())