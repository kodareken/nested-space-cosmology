"""The programme's master index: fourteen parts, their records, one table.

The strand's synthesis lives in the owning note; this driver writes the
machine-readable counterpart: every part, its headline numbers, and its
records, with every record path existence-checked before the index is
sealed. Navigation for a human, provenance for a script.

    python scripts/lab.py scripts/derive_nsc_programme_index.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-programme-index-v1.json"
SCHEMA = "NSC-PROGRAMME-INDEX-v1"

PARTS = [
    ("exterior channels closed",
     "fees 0.0446/0.2517/0.2610 mass 0; LLL 103x; boundary blocked",
     ["nsc-influence-chain-v1.json", "nsc-influence-shallow-power-v1.json",
      "nsc-influence-response-v1.json", "nsc-influence-seam-v1.json"]),
    ("interior funnel exact",
     "dM = d(excision) 3.2e-15; outer leak-only every interval",
     ["nsc-interior-balance-v1.json", "nsc-continuation-station-closure-v1.json"]),
    ("wall converter",
     "16-22x concentration; lag +1; E/M conjugate -0.99; mu=1.040 lock",
     ["nsc-interior-station-balance-v1.json", "nsc-wall-lead-lag-v1.json",
      "nsc-wall-energy-channel-v1.json", "nsc-wall-ratio-regimes-v1.json"]),
    ("gap-clock organism",
     "rise d^-1/2 4.0%, fall d^-1/4 3.3%; eps-cross 3.5/6.8%; omega 0.386",
     ["nsc-wall-clock-horizon-v1.json", "nsc-wall-clock-fall-v1.json",
      "nsc-wall-clock-epscross-v1.json", "nsc-wall-clock-derivation-v1.json"]),
    ("universal scaffold",
     "seven mu-blind quantities incl. redshift 8e-7, I_first 0.24%",
     ["nsc-seam-invariance-v1.json", "nsc-charged-seam-invariance-v1.json",
      "nsc-cone-wall-join-v1.json", "nsc-isotropy-driven-join-v1.json"]),
    ("founding riddle positives",
     "band at seams [25.9, 29.0]; target 0.049 at 28.17; shell weights",
     ["nsc-cosmic-band-location-v1.json", "nsc-cosmic-shell-weight-v1.json",
      "nsc-time-awareness-v1.json"]),
    ("cross-instrument anchors",
     "omega=4 join; omega_star identity 2.2e-16; two clocks",
     ["nsc-matryoshka-fee-join-v1.json", "nsc-omega-star-identity-v1.json",
      "nsc-charged-clock-horizon-v1.json"]),
    ("fringe staircase",
     "quantum 0.00131; ladder {3,3,2,1}; q*N 1.149; beta threshold",
     ["nsc-tilt-quantum-v1.json", "nsc-tilt-quantum-ladder-v1.json",
      "nsc-neutral-fringe-quantum-v1.json", "nsc-fringe-beta-join-v1.json"]),
    ("the joins",
     "two-ends corr -0.93; one-way band placement; cone nesting; isotropy locality",
     ["nsc-two-ends-join-v1.json", "nsc-oneway-wall-join-v1.json",
      "nsc-cone-wall-join-v1.json", "nsc-isotropy-wall-join-v1.json"]),
    ("late harvest",
     "identity predictive 2.06 vs 2.07%; composed 0.948 of clock; trap r=0.835/0.999",
     ["nsc-ladder-covariance-test-v1.json", "nsc-composed-vs-clock-v1.json",
      "nsc-continuation-trap-v1.json", "nsc-trap-timing-v1.json"]),
    ("production door",
     "preflight passed; freeze gate refused (pre-existing dirty producer)",
     ["nsc-production-preflight-v1.json", "nsc-production-gate-v1.json"]),
    ("wall lifecycle",
     "death [10,15]; t=40 event 3 witnesses; halving 7.7 t; single force 0.9975",
     ["nsc-relaxation-wall-v1.json", "nsc-wall-afterlife-v1.json",
      "nsc-march-balance-trajectory-v1.json", "nsc-wall-extinction-v1.json"]),
    ("calibration standard",
     "cross-arm bit-identity; frontier amendment 2 cells left",
     ["nsc-cross-arm-identity-v1.json", "nsc-frontier-amendment-v1.json"]),
    ("domain law",
     "deep-grid radii inside horizons; amplitude ratio 7/17/40; standoff law",
     ["nsc-march-wall-radii-v1.json", "nsc-wall-domain-law-v1.json"]),
]


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    development = LAB / "results" / "development"
    missing = []
    parts = []
    for title, headline, records in PARTS:
        records_ok = []
        for name in records:
            if not (development / name).is_file():
                missing.append(name)
            records_ok.append({"record": name, "exists": (development / name).is_file()})
        parts.append({"title": title, "headline": headline, "records": records_ok})
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "parts": parts,
        "part_count": len(parts),
        "missing_records": missing,
        "index_complete": not missing,
        "notes": [
            "headlines are summaries; the papers of record are the "
            "individual creation-only records and the owning note",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_programme_index.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_programme_index.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  parts:", len(parts), " records cited:", sum(len(p['records']) for p in parts))
    print("  missing:", missing if missing else "none - index complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())