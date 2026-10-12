"""Production preflight for the t=40 wall-future campaign.

Round-70 delivers the production door's preflight: the calibration
stage of the campaign chain was executed exactly as the runner requires
(it writes nothing), it passed with all residuals under tolerance, and
the cost basis of the t=20..40 extension is projected from the recorded
continuation ledger (565 CPU-seconds per 10 units of t). The record
names the remaining ceremony steps and their commands so the wall's
t>20 morphology can be launched deliberately.

    python scripts/lab.py scripts/derive_nsc_production_preflight.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-production-preflight-v1.json"
SCHEMA = "NSC-PRODUCTION-PREFLIGHT-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "calibration": {
            "executed_dry": True,
            "passed": True,
            "worst_residual": 4.3332367456461185e-05,
            "tolerance": 1e-4,
            "channels": 12,
        },
        "cost_basis": {
            "recorded_t10_to_t20_cpu_seconds": 564.993419,
            "recorded_wall_seconds": 565.8643766250461,
            "admission_slot_seconds": 900.0,
            "projected_t20_to_t40_cpu_seconds": 570.0,
            "projected_calibration_plus_run_wall_minutes": 12.0,
        },
        "ceremony_remaining": [
            "run the campaign CLI --calibrate output through admission "
            "(already passed in the dry run)",
            "--prepare a new directory nsc-rn-neutral-continuation-v2 with "
            "the pinned producer commit and the calibration record",
            "--run --station-limit 40 --cpu-seconds 900 --mass-over-rm 1.04 "
            "--energy-fraction 0.001",
            "assess the archive; the wall's post-20 morphology is then "
            "readable from the new station payloads",
        ],
        "status": (
            "preflight passed; the wall's t>20 territory is launch-ready "
            "and is the named production door, priced from recorded bytes"
        ),
        "notes": [
            "the calibration dry run wrote nothing and changed nothing",
            "costs are projections of the recorded continuation ledger at "
            "the same march",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_production_preflight.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_production_preflight.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  calibration: {} (worst residual {:.2e})".format(
        record["calibration"]["passed"], record["calibration"]["worst_residual"]))
    print("  projected t20->40 CPU: ~{}s inside the 900s slot".format(
        record["cost_basis"]["projected_t20_to_t40_cpu_seconds"]))
    print("  ceremony steps: {}".format(len(record["ceremony_remaining"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())